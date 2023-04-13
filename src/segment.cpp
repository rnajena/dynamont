#include <iostream>
#include <fstream> // file io
#include <sstream> // file io
#include <string>
#include <map> // dictionary
#include <tuple>
#include <algorithm> // reverse strings
#include <vector>
#include <cmath> // exp
#include <limits> // for inifinity
#include <math.h> // log1p
// #include <iterator>

using namespace std;

map<char, int> BASE2ID;
int ALPHABET_SIZE;

// TODO move to config file?
const string MODELPATH = "data/template_median69pA.model";
const string TERM_STRING = "€";
const int K = 5; // our model works with this kmer size

// Asserts floating point compatibility at compile time
// necessary for INFINITY usage
static_assert(numeric_limits<float>::is_iec559, "IEEE 754 required");

/**
 * Fill the BASE2ID map with the base and ID pairs.
 */
void fillBASE2ID() {
    ALPHABET_SIZE = 5;
    BASE2ID.insert(pair<char, int>('A', 0));
    BASE2ID.insert(pair<char, int>('C', 1));
    BASE2ID.insert(pair<char, int>('G', 2));
    BASE2ID.insert(pair<char, int>('T', 3));
    BASE2ID.insert(pair<char, int>('N', 4));
    BASE2ID.insert(pair<char, int>('a', 0));
    BASE2ID.insert(pair<char, int>('c', 1));
    BASE2ID.insert(pair<char, int>('g', 2));
    BASE2ID.insert(pair<char, int>('t', 3));
    BASE2ID.insert(pair<char, int>('n', 4));
}

/**
 * Converts a number of a given base (ALPHABET_SIZE) to a decimal number.
 * Works ONLY if ALPHABET_SIZE is smaller or equal to 10!
 * 
 * @param i input number in the given base
*/
int toDeci(int i) {
    int ret = 0, r = 0;
    int m = 1;
    while(i > 0) {
        r = i % 10;
        ret += m*r;
        m *= ALPHABET_SIZE;
        i /= 10;
    }
    return ret;
}

/**
 * Converts a number of a given base (ALPHABET_SIZE) to a decimal number.
 * Works ONLY if ALPHABET_SIZE is smaller or equal to 10!
 * 
 * @param i input number in the given base as an array
*/
int toDeci(int *i) {
    int* c = new int[K];
    copy(i, i+K, c);
    reverse(c, c+K); //array is a kmer, so should have length K
    int ret = 0;
    int m = 1;
    for(int r = 0; r < ALPHABET_SIZE; r++) {
        ret += m*c[r];
        m *= ALPHABET_SIZE;
    }
    return ret;
}

/**
 * Converts the kmer to the integer representation using the BASE2ID map
 *
 * @param s kmer containing nucleotides 
 * @return integer representation of the given kmer
 */
int kmer2int(const string &s) {
    int i = 0;
    for (char const &c:s){
        i *= 10; // move the number to the left
        i+=BASE2ID.at(c);
    }
    return toDeci(i);
}

// https://stackoverflow.com/questions/10847007/using-the-gaussian-probability-density-function-in-c
/**
 * Calculate probabilty value for a given x, mean and standard deviation of a normal distribution
 *
 * @param x value
 * @param m mean
 * @param s standard deviation 
 * @return probabily density at position x for N~(m, s²)
 */
float normal_pdf(const float &x, const float &m, const float &s) {
    static const float INV_SQRT_2PI = 0.3989422804014327;
    const float a = (x - m) / s;

    return INV_SQRT_2PI / s * exp(-0.5f * a * a);
}

/**
 * Return probability density for a given value and a given normal distribution
 *
 * @param x point to calculate probability density
 * @param kmer key for the model kmer:(mean, stdev) map
 * @param model map containing kmers as keys and (mean, stdev) tuples as values
 * @return probability density value for x in the given normal distribution
 */
float scoreKmer(const float &x, const int &kmer, vector<tuple<float, float>>* model) {
    tuple<float, float> kmerModel = (*model)[kmer];
    return normal_pdf(x, get<0>(kmerModel), get<1>(kmerModel));
}

// /**
//  * Slices seq and writes the sliced interval [j,k) into s
//  * 
//  * @param seq array to be sliced
//  * @param s array to write slice into
//  * @param j start of slice (inclusive)
//  * @param k end of slice (exclusive)
// */
// void sliceSeq(const int* seq, int* s, const int &start, const int &end) {
//     int i = 0;
//     for(int j=start; j<end; j++, i++) {
//         s[i] = seq[j];
//     }
// }

/**
 * Calculate addition of a+b in log space as efficiently as possible
 *
 * @param a first value
 * @param b second value
 * @return log(exp(a) + exp(b))
 */
float logPlus(const float &a, const float &b) {
    // safety check
    if(a==b && isinf(a) && isinf(b)) {
        return a;
    } 
    if(a>=b){
        return a + log1p(exp(b-a));
    }
    return b + log1p(exp(a-b));
}

/**
 * Calculate forward matrices using logarithmic values
 *
 * @param sig ONT raw signal with pA values
 * @param seq nucleotide sequence represented by the ONT signal
 * @param MM matrix containing forward-values for segment borders 
 * @param MC matrix containing forward-values for extending segment
 * @param T length of the ONT raw signal + 1
 * @param N length of nucleotide sequence
 * @param model map containing kmers as keys and (mean, stdev) tuples as values
 */
void logF(float* sig, int* seq, float* MM, float* MC, const int &T, const int &N, vector<tuple<float, float>>* model){
    int* tempKmer = new int[K];
    // fill_n(tempKmer, K, -1);
    int kmer = -1;
    float mm, mc;
    // i iterates through the signal, T = len(sig) + 1
    for(int i=0; i<T; i++){
        // j iterates through the read sequence
        for(int j=0; j<N; j++){
            copy(seq + (j-(K/2)+1), seq + (j+(K/2)+2), tempKmer); // add +2 bc of index transformation
            kmer = toDeci(tempKmer); // convert kmer of tokens to integer ID
            mm=-INFINITY;
            if(i>0 && j>0){
              mm = logPlus(mm, MC[(i-1)*N+(j-1)] + log(scoreKmer(sig[i-1], kmer, model)));
            }
            MM[i*N+j] = mm;
            mc=-INFINITY;
            if(i>0){
              mc = logPlus(mc, MC[(i-1)*N+j] + log(scoreKmer(sig[i-1], kmer, model)));
              mc = logPlus(mc, MM[(i-1)*N+j] + log(scoreKmer(sig[i-1], kmer, model)));
            }
            if(i==0 && j==0){
                mc = 0; // initialize with log(1) 
            }
            MC[i*N+j]=mc;
        }
    }
}

/**
 * Calculate backward matrices using logarithmic values
 *
 * @param sig ONT raw signal with pA values
 * @param seq nucleotide sequence represented by the ONT signal
 * @param MM matrix containing backward-values for segment borders
 * @param MC matrix containing backward-values for extending segment
 * @param T length of the ONT raw signal + 1
 * @param N length of nucleotide sequence
 * @param model map containing kmers as keys and (mean, stdev) tuples as values
 */
void logB(float* sig, int* seq, float* MM, float* MC, const int &T, const int &N, vector<tuple<float, float>>* model) {
    int* tempKmer = new int[K];
    int kmer = -1;
    float mm, mc;
    for(int i=T-1; i>=0; i--){
        for(int j=N-1; j>=0; j--){
            copy(seq + (j-(K/2)+1), seq + (j+(K/2)+2), tempKmer);
            kmer = toDeci(tempKmer);
            mm=-INFINITY;
            if(i<T-1 && j<N-1){
                mm = logPlus(mm, MC[(i+1)*N+(j+1)] + log(scoreKmer(sig[i-1], kmer, model)));
            }
            MM[i*N+j] = mm;
            mc=-INFINITY;
            if(i<T-1){
                mc = logPlus(mc, MC[(i+1)*N+j] + log(scoreKmer(sig[i-1], kmer, model)));
                mc = logPlus(mc, MM[(i+1)*N+j] + log(scoreKmer(sig[i-1], kmer, model)));
            }
            if(i==T-1 && j==N-1){
                mc = 0; // initialize with log(1) 
            }
            MC[i*N+j]=mc;
        }
    }
}

/**
 * Calculate the logarithmic probability matrix
 *
 * @param forMM matrix containing forward-values for segment borders
 * @param forMC matrix containing forward-values for extending segment
 * @param backMM matrix containing backward-values for extending segment
 * @param backMC matrix containing backward-values for extending segment
 * @param T length of the ONT raw signal + 1
 * @param N length of nucleotide sequence
 * @return matrix containing logarithmic probabilities for segment borders
 */
float* logP(float* forMM, float* forMC, float* backMM, float* backMC, const int &T, const int &N) {
    float forZ = forMC[T*N -2]; // -2 because we exclude the first and last 2 nucleotides for segmentation for now
    // float backZ = backMC[2]; // 2 because we exclude the first and last 2 nucleotides for segmentation for now
    // TODO? check if forZ ~= backZ
    float* LP = new float[T*N];
    fill_n(LP, T*N, -INFINITY);
    // LP = {-1.0f};
    for(int i=0; i<T; i++){
        for(int j=0; j<N; j++){
            int x = i*N+j;
            LP[x] = forMM[x] + backMM[x] - forZ;
        }
    }
    // TODO? check LP: no value > log(1) (=0) allowed
    return LP;
}

/**
 * Read the normal distribution parameters from a given TSV file
 *
 * @param file path to the TSV file containing the parameters
 * @param model kmer model to fill
 */
void readKmerModel(const string &file, vector<tuple<float, float>>* model) {
    ifstream inputFile(file);
    string line, kmer, tmp;
    float mean, stdev;
    getline(inputFile, line);

    while(getline(inputFile, line)) { // read line
        stringstream buffer(line); // parse line to stringstream for getline
        getline(buffer, kmer, '\t');
        reverse(kmer.begin(), kmer.end()); // 3-5 -> 5-3 orientation
        getline(buffer, tmp, '\t');
        mean = atof(tmp.c_str());
        getline(buffer, tmp, '\t');
        stdev = atof(tmp.c_str());
        (*model)[kmer2int(kmer)]=make_tuple(mean, stdev);
    }
    inputFile.close();
}

int main(int argc, char* argv[]) {
    cout<<"START\n";
    // read data from stdin
    string input;
    string signal;
    string read;
    fillBASE2ID();
    int truish = 1;

    while(truish) {
        truish = 0;
        // 107.2,108.0,108.9,111.2,105.7,104.3,107.1,105.7 CAAAAA
        // read input, signal and read whitespace separated in single line
        getline(cin, signal, ' ');
        getline(cin, read);

        std::cout << "SSS " << signal << std::endl;
        std::cout << "RRR " << read << std::endl;

        // break loop if termination character ...
        if (signal.find(TERM_STRING) != string::npos) {
            break;
        // ... or signal or read is missing
        } else if (signal.empty()) {
            cout<<"Signal missing!\n";
            break;
        } else if (read.empty()) {
            cout<<"Read missing!\n";
            break;
        }
        
        // process signal: convert string to float array
        int T = count(signal.begin(), signal.end(), ',')+2; // len(sig) + 1
        float* sig = new float[T - 1];
        fill_n(sig, T-1, -INFINITY);
        string value;
        stringstream ss(signal);
        int i = 0;
        while(getline(ss, value, ',')) {
            sig[i++] = stof(value);
        }

        // process read: convert string to int array
        // we exclude the first and last 2 nucleotides for segmentation for now
        int N = read.size() - (K-2); // ignore first and last nucleotides that cannot form a Kmer
        int* seq = new int[read.size()];
        fill_n(seq, read.size(), -1);
        i = 0;
        for (const char &c: read) { //string::size_type i = 0; i < read.size(); i++) {
            try {
                seq[i++] = BASE2ID.at(c);
            } catch (const std::out_of_range) {
                // TODO handle unknown nucleotide in read
                cout<<"Unknown nucleotide: "<<c<<endl;
                break;
            }
        }

        // initialize matrices
        float* forMM = new float[T*N];
        fill_n(forMM, T*N, -INFINITY);
        float* forMC = new float[T*N];
        fill_n(forMC, T*N, -INFINITY);
        float* backMM = new float[T*N];
        fill_n(backMM, T*N, -INFINITY);
        float* backMC = new float[T*N];
        fill_n(backMC, T*N, -INFINITY);
        
        vector<tuple<float, float>> model(pow(ALPHABET_SIZE, K), make_tuple(-INFINITY, -INFINITY));
        readKmerModel(MODELPATH, &model);

        // calculate segmentation probabilities, fill matrices
        logF(sig, seq, forMM, forMC, T, N, &model);

        cout<<"forMM\n";
        for(int i=0; i<T; i++){
            for(int j=0; j<N; j++){
                cout<<forMM[i*N+j]<<", ";
            }
            cout<<endl;
        }

        cout<<"forMC\n";
        for(int i=0; i<T; i++){
            for(int j=0; j<N; j++){
                cout<<forMC[i*N+j]<<", ";
            }
            cout<<endl;
        }

        //logB(sig, seq, backMM, backMC, T, N, &model);
        //
        //cout<<"backMM\n";
        //for(int i=0; i<T; i++){
        //    for(int j=0; j<N; j++){
        //        cout<<backMM[i*N+j]<<", ";
        //    }
        //    cout<<endl;
        //}

        //cout<<"backMC\n";
        //for(int i=0; i<T; i++){
        //    for(int j=0; j<N; j++){
        //        // cout<<backMC[i*N+j]<<", ";
        //    }
        //    cout<<endl;
        //}

        // float* LP = logP(forMM, forMC, backMM, backMC, T, N);

        // temporary testcode
        // cout<<"LP: "<<*LP<<endl;
        // cout<<"LP\n";
        // for(int i=0; i<T; i++){
        //     for(int j=0; j<N; j++){
        //         int x = i*N+j;
        //         cout<<LP[x]<<", ";
        //     }
        //     cout<<endl;
        // }


        // TODO: calculate where to put segment the signal

        // TODO: write solution to stdout

        // TODO: remove temporary data

    }

    cout<<"End\n";

    return 0;
}
