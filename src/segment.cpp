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
    int ret = 0;
    int m = 1;
    for(int r = K - 1; r >= 0; r--) {
        ret += m*i[r];
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
            if (j>0) {
                copy(seq + (j+1-(K/2)), seq + (j+2+(K/2)), tempKmer);
                kmer = toDeci(tempKmer); // convert kmer of tokens to integer ID
            }
            mm=-INFINITY;
            mc=-INFINITY;
            if (i>0 && j>0){
                mm = logPlus(mm, MC[(i-1)*N+(j-1)] + log(scoreKmer(sig[i-1], kmer, model)));
            }
            MM[i*N+j] = mm;
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
    delete[] tempKmer;
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
            // cerr<<"IT: "<<i<<", "<<j<<endl;
            if (j>0) {
                copy(seq + (j+1-(K/2)), seq + (j+2+(K/2)), tempKmer);
                kmer = toDeci(tempKmer); // convert kmer of tokens to integer ID
            }
            mm=-INFINITY;
            mc=-INFINITY;
            if(i==T-1 && j==N-1){
                mc = 0; // initialize with log(1) 
            }
            if(i<T-1){
                mc = logPlus(mc, MC[(i+1)*N+j] + log(scoreKmer(sig[i], kmer, model)));
                mm = logPlus(mm, MC[(i+1)*N+j] + log(scoreKmer(sig[i], kmer, model)));
                // cerr<<"log(score) "<<log(scoreKmer(sig[i], kmer, model))<<endl;
                // cerr<<"sig[i]: "<<sig[i]<<", kmer: "<<kmer<<endl;
            }
            if(i<T-1 && j<N-1){
                copy(seq + (j+2-(K/2)), seq + (j+3+(K/2)), tempKmer);
                kmer = toDeci(tempKmer); // convert kmer of tokens to integer ID
                mc = logPlus(mc, MM[(i+1)*N+(j+1)] + log(scoreKmer(sig[i], kmer, model)));
            }
            MC[i*N+j]=mc;
            MM[i*N+j] = mm;
        }
    }
    delete[] tempKmer;
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
    float forZ = forMC[T*N-1]; // -2 because we exclude the first and last 2 nucleotides for segmentation for now
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
    // cout<<"Z: "<<forZ<<endl;
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
    // cout<<"START\n";
    // read data from stdin
    string input;
    string signal;
    string read;
    fillBASE2ID();
    int truish = 1;

    while(truish) {
        truish = 0;
        // echo 107.2,108.0,108.9,111.2,105.7,104.3,107.1,105.7 CAAAAA| src\segment.exe
        // read input, signal and read whitespace separated in single line
        getline(cin, signal, ' ');
        getline(cin, read);

        // std::cout << "SSS " << signal << std::endl;
        // std::cout << "RRR " << read << std::endl;

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

        // cout<<"forMM\n";
        // for(int i=0; i<T; i++){
        //     for(int j=0; j<N; j++){
        //         cout<<forMM[i*N+j]<<", ";
        //     }
        //     cout<<endl;
        // }

        // cout<<"forMC\n";
        // for(int i=0; i<T; i++){
        //     for(int j=0; j<N; j++){
        //         cout<<forMC[i*N+j]<<", ";
        //     }
        //     cout<<endl;
        // }

        logB(sig, seq, backMM, backMC, T, N, &model);
        
        // cout<<"backMM\n";
        // for(int i=0; i<T; i++){
        //    for(int j=0; j<N; j++){
        //        cout<<backMM[i*N+j]<<", ";
        //    }
        //    cout<<endl;
        // }

        // cout<<"backMC\n";
        // for(int i=0; i<T; i++){
        //    for(int j=0; j<N; j++){
        //        cout<<backMC[i*N+j]<<", ";
        //    }
        //    cout<<endl;
        // }

        float* LP = logP(forMM, forMC, backMM, backMC, T, N); // log probs
        
        // temporary testcode
        // cout<<"LP\n";
        // int x = 0;
        // for(int i=0; i<T; i++){
        //     for(int j=0; j<N; j++){
        //         x = i*N+j;
        //         cout<<LP[x]<<", ";
        //     }
        //     cout<<endl;
        // }

        // calculate where to put segment the signal
        // calculate max or sum of rows, sum better if uncertain
        float* LSP = new float[T]; // log segment probs
        int idx = 0;
        float sum = -INFINITY;
        for(int i=0;i<T;i++) {
            sum = -INFINITY;
            for(int j=0; j<N; j++) {
                idx = i*N+j;
                sum = logPlus(sum, LP[idx]);
            }
            LSP[i] = sum;
            cout<<sum<<",";
        }
        cout<<endl;

        // TODO: remove data
        delete[] LP;
        LP = 0;
        delete[] LSP;
        LSP = 0;
        delete[] sig;
        sig = 0;
        delete[] seq;
        seq = 0;
        delete[] forMM;
        forMM = 0;
        delete[] forMC;
        forMC = 0;
        delete[] backMM;
        backMM = 0;
        delete[] backMC;
        backMC = 0;
    }

    // cout<<"End\n";

    return 0;
}
