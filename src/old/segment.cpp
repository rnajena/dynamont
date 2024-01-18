// author: Jannes Spangenberg
// e-mail: jannes.spangenberg@uni-jena.de
// github: https://github.com/JannesSP
// website: https://jannessp.github.io

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
#include <assert.h>
// #include <iomanip>
// #include <iterator>

using namespace std;

map<char, int> BASE2ID;
int ALPHABET_SIZE;
double EPSILON = pow(10, -5);

// TODO move to config file?
// const string MODELPATH = "/home/yi98suv/projects/ont_segmentation/data/template_median69pA.model";
const string MODELPATH = "/home/yi98suv/projects/ont_segmentation/data/template_median69pA_extended.model";
const string TERM_STRING = "$";
const int K = 5; // our model works with this kmer size

// Asserts doubleing point compatibility at compile time
// necessary for INFINITY usage
static_assert(numeric_limits<double>::is_iec559, "IEEE 754 required");

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
 * Converts a number of base ALPHABET_SIZE to a decimal number.
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
 * Converts a number of base ALPHABET_SIZE to a decimal number.
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
        assert (BASE2ID.at(c)>=0);
        i *= 10; // move the number to the left
        i+=BASE2ID.at(c);
    }
    int ret = toDeci(i);
    return ret;
}

/**
 * Convert the read sequence to a kmer sequence which is represented by integers.
 * 
 * @param seq read sequence
 * @param N length of the read sequence, number of nucleotides
 * @return kmer sequence in integer representation
*/
int* seq2kmer(int* seq, const int &N) {
    int* kmer_seq = new int[N];
    int* tempKmer = new int[K];
    for (int i=0; i<N; i++){ // extend loop to ad 2 Ns at start and end of read        
        copy(seq + i, seq + i+K, tempKmer);
        kmer_seq[i] = toDeci(tempKmer);
    }
    return kmer_seq;
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
double normal_pdf(const double &x, const double &m, const double &s) {
    static const double INV_SQRT_2PI = 0.3989422804014327;
    const double a = (x - m) / s;
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
double scoreKmer(const double &x, const int &kmer, vector<tuple<double, double>>* model) {
    tuple<double, double> kmerModel = (*model)[kmer];
    return normal_pdf(x, get<0>(kmerModel), get<1>(kmerModel));
}

/**
 * Calculate addition of a+b in log space as efficiently as possible
 *
 * @param a first value
 * @param b second value
 * @return log(exp(a) + exp(b))
 */
double logPlus(const double &a, const double &b) {
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
 * @param N length of nucleotide sequence + 1
 * @param model map containing kmers as keys and (mean, stdev) tuples as values
 */
void logF(double* sig, int* kmer_seq, double* MM, double* MC, const int &T, const int &N, vector<tuple<double, double>>* model){
    double mm, mc;
    for(int i=0; i<T; i++){
        for(int j=0; j<N; j++){
            mm=-INFINITY;
            mc=-INFINITY;
            if (i>0 && j>0){
                mm = logPlus(mm, MC[(i-1)*N+(j-1)] + log(scoreKmer(sig[i-1], kmer_seq[j-1], model)));
                mc = logPlus(mc, MC[(i-1)*N+j] + log(scoreKmer(sig[i-1], kmer_seq[j-1], model)));
                mc = logPlus(mc, MM[(i-1)*N+j] + log(scoreKmer(sig[i-1], kmer_seq[j-1], model)));
                // mm = logPlus(mm, MM[(i-1)*N+(j-1)] + log(scoreKmer(sig[i-1], kmer_seq[j-1], model)));
            }
            MM[i*N+j] = mm;
            // if(i>0){
            // }
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
 * @param N length of nucleotide sequence + 1
 * @param model map containing kmers as keys and (mean, stdev) tuples as values
 */
void logB(double* sig, int* kmer_seq, double* MM, double* MC, const int &T, const int &N, vector<tuple<double, double>>* model) {
    double mm, mc;
    for(int i=T-1; i>=0; i--){
        for(int j=N-1; j>=0; j--){
            mm=-INFINITY;
            mc=-INFINITY;
            if(i==T-1 && j==N-1){
                mc = 0; // initialize with log(1)
            }
            if(i<T-1 && j>0){
                mc = logPlus(mc, MC[(i+1)*N+j] + log(scoreKmer(sig[i], kmer_seq[j-1], model)));
                mm = logPlus(mm, MC[(i+1)*N+j] + log(scoreKmer(sig[i], kmer_seq[j-1], model)));
            }
            if(i<T-1 && j<N-1){
                mc = logPlus(mc, MM[(i+1)*N+(j+1)] + log(scoreKmer(sig[i], kmer_seq[j], model)));
                // mm = logPlus(mm, MM[(i+1)*N+(j+1)] + log(scoreKmer(sig[i], kmer_seq[j], model)));
            }
            MC[i*N+j] = mc;
            MM[i*N+j] = mm;
        }
    }
}

/**
 * Calculate the logarithmic probability matrix
 *
 * @param FOR matrix containing forward-values for segment borders
 * @param BACK matrix containing backward-values for extending segment
 * @param T length of the ONT raw signal + 1
 * @param N length of nucleotide sequence + 1
 * @return matrix containing logarithmic probabilities for segment borders
 */
double* logP(double* FOR, double* BACK, const double &Z, const int &T, const int &N) {
    // cerr<<"Z: "<<Z;
    double* LP = new double[T*N];
    fill_n(LP, T*N, -INFINITY);
    for(int i=0; i<T; i++){
        for(int j=0; j<N; j++){
            int x = i*N+j;
            LP[x] = FOR[x] + BACK[x] - Z;
        }
    }
    // cout<<"Z: "<<forZ<<endl;
    // TODO? check LP: no value > log(1) (=0) allowed
    return LP;
}

/**
 * Calculate the maximum a posteriori path through LP
 *
 */
int* getBorders(double* LPM, double* LPC, const int &T, const int &N){
    
    double* MM = new double[T*N];
    fill_n(MM, T*N, -INFINITY);

    double mm;
    for(int i=0; i<T; i++){
        for(int j=0; j<N; j++){
            mm=-INFINITY;
            if(i>0 && j>0){
                mm = MM[(i-1)*N+(j-1)] + LPM[i*N+j];
            }
            if(i>0){
                mm = max(mm, MM[(i-1)*N+j] + LPC[i*N+j]);
            }
            if(i==0 && j==0){
                mm = 0; // initialize with log(1)
            }
            MM[i*N+j] = mm;
        }
    }

    int* borderMap = new int[N];
    fill_n(borderMap, N, 0);
    int i = T-1;
    int j = N-1;
    while(j>0 && i>0){
        if(i>0 && j>0 && (MM[(i-1)*N+(j-1)] + LPM[i*N+j] == MM[i*N+j])){
            borderMap[j]=i;
            i--;
            j--;
        } else {
            i--;
        }
    }

    delete[] MM;
    return borderMap;
}

/**
 * Read the normal distribution parameters from a given TSV file
 *
 * @param file path to the TSV file containing the parameters
 * @param model kmer model to fill
 */
void readKmerModel(const string &file, vector<tuple<double, double>>* model) {
    ifstream inputFile(file);
    string line, kmer, tmp;
    double mean, stdev;
    getline(inputFile, line);
    while(getline(inputFile, line)) { // read line
        stringstream buffer(line); // parse line to stringstream for getline
        getline(buffer, kmer, '\t');
        // reverse(kmer.begin(), kmer.end()); // 3-5 -> 5-3 orientation
        getline(buffer, tmp, '\t');
        mean = atof(tmp.c_str());
        getline(buffer, tmp, '\t');
        stdev = atof(tmp.c_str());
        (*model)[kmer2int(kmer)]=make_tuple(mean, stdev);
    }
    inputFile.close();
}

/**
 * Read signal and read from stdin until the TERM_STRING is seen
*/
int main(int argc, char* argv[]) {
    // cout<<"START\n";
    // read data from stdin
    string input;
    string signal;
    string read;
    fillBASE2ID();
    int truish = 1;
    vector<tuple<double, double>> model(pow(ALPHABET_SIZE, K), make_tuple(-INFINITY, -INFINITY));
    readKmerModel(MODELPATH, &model);

    while(truish) {
        // truish = 0;
        // echo 107,107,107.2,108.0,108.9,111.2,105.7,104.3,107.1,105.7,105,105 CAAAAA| src\segment.exe
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
        
        // process signal: convert string to double array
        int T = count(signal.begin(), signal.end(), ',')+2; // len(sig) + 1
        // cerr<<"T: "<<T<<endl;
        double* sig = new double[T - 1];
        fill_n(sig, T-1, -INFINITY);
        string value;
        stringstream ss(signal);
        int i = 0;
        while(getline(ss, value, ',')) {
            sig[i++] = stof(value);
        }
        // process read: convert string to int array
        int N = read.size() + 1; // operate on base transitions
        // int seq_size = read.size() + (2*K+2); // add NNNNNAAAAA in front and NN to end
        int seq_size = read.size() + (K-1); // add 2 nucleotides to front and back of read
        int* seq = new int[seq_size];
        // fill_n(seq, seq_size, 0); // default: fill with A
        fill_n(seq, seq_size, 0); // default: fill with A
        // add NNNNN to front
        // for (int i = 0; i<K; i++) {
        //     seq[i] = 4;
        // }
        // i = 2*K; // start at position K because of As at front
        i = floor(K/2);
        for (const char &c: read) {
            try {
                seq[i] = BASE2ID.at(c);
                i++;
            } catch (const std::out_of_range) {
                // TODO handle unknown nucleotide in read
                cerr<<"Unknown nucleotide: "<<c<<endl;
                exit(10);
            }
        }
        // add NN to end of sequence
        seq[i] = 4;
        seq[i+1] = 4;

        int* kmer_seq = seq2kmer(seq, N-1);

        // initialize matrices
        double* forMM = new double[T*N];
        fill_n(forMM, T*N, -INFINITY);
        double* forMC = new double[T*N];
        fill_n(forMC, T*N, -INFINITY);
        double* backMM = new double[T*N];
        fill_n(backMM, T*N, -INFINITY);
        double* backMC = new double[T*N];
        fill_n(backMC, T*N, -INFINITY);

        // calculate segmentation probabilities, fill matrices
        logF(sig, kmer_seq, forMM, forMC, T, N, &model);
        logB(sig, kmer_seq, backMM, backMC, T, N, &model);
        
        // Check if Z value matches, must match between matrices
        if (fabs(forMC[N*T-1] - backMC[0])>EPSILON) {
            cerr.precision(20);
            cerr<<"Z values between matrices do not match!\n";
            cerr<<"forMC: "<<forMC[N*T-1]<<", backMC: "<<backMC[0]<<"\n";
            cerr<<backMC[0]-forMC[N*T-1]<<"\n";
            cerr.flush();
            exit(11);
        }

        double* LPM = logP(forMM, backMM, backMC[0], T, N); // log probs
        double* LPC = logP(forMC, backMC, backMC[0], T, N); // log probs
        int* borderMap = getBorders(LPM, LPC, T, N);

        for(int j = 0; j<N; j++){
            cout<<borderMap[j]<<",";
        }
        cout<<endl;
        // cout.flush();

        // calculate sum of segment probabilities
        double* LSP = new double[T]; // log segment probs
        fill_n(LSP, T, -INFINITY);
        int idx;
        double sum;
        for(int i=0;i<T;i++) {
            sum = -INFINITY;
            for(int j=0; j<N; j++) {
                idx = i*N+j;
                sum = logPlus(sum, LPM[idx]);
            }
            LSP[i] = sum;
            cout<<sum<<",";
        }
        cout<<endl;
        cout.flush();

// ============================================

        // cerr.precision(4);
        // cerr<<"LP\n";
        // for(int i=0; i<T; i++){
        //    for(int j=0; j<N; j++){
        //        cerr<<LP[i*N+j]<<", ";
        //    }
        //    cerr<<endl;
        // }
        // cerr.precision(4);
        // cerr<<"forMM\n";
        // for(int i=0; i<T; i++){
        //    for(int j=0; j<N; j++){
        //        cerr<<forMM[i*N+j]<<", ";
        //    }
        //    cerr<<endl;
        // }
        // cerr<<"forMC\n";
        // for(int i=0; i<T; i++){
        //    for(int j=0; j<N; j++){
        //        cerr<<forMC[i*N+j]<<", ";
        //    }
        //    cerr<<endl;
        // }
        // cerr<<"backMM\n";
        // for(int i=0; i<T; i++){
        //    for(int j=0; j<N; j++){
        //        cerr<<backMM[i*N+j]<<", ";
        //    }
        //    cerr<<endl;
        // }
        // cerr<<"backMC\n";
        // for(int i=0; i<T; i++){
        //    for(int j=0; j<N; j++){
        //        cerr<<backMC[i*N+j]<<", ";
        //    }
        //    cerr<<endl;
        // }

        // Clean up
        delete[] LPM, LPC,borderMap, sig, seq, kmer_seq, forMM, forMC, backMM, backMC;
    }
    return 0;
}
