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
#include "argparse.hpp"
// #include <iomanip>
// #include <iterator>

using namespace std;

void funcM(int i, int j, double* P, double* M, double* E, double* I, double* S, double* LPP, double* LPM, double* LPE, double* LPI, double* LPS, list<string>* segString, const int &T, const int &N);
void funcE(int i, int j, double* P, double* M, double* E, double* I, double* S, double* LPP, double* LPM, double* LPE, double* LPI, double* LPS, list<string>* segString, const int &T, const int &N);
void funcI(int i, int j, double* P, double* M, double* E, double* I, double* S, double* LPP, double* LPM, double* LPE, double* LPI, double* LPS, list<string>* segString, const int &T, const int &N);
void funcS(int i, int j, double* P, double* M, double* E, double* I, double* S, double* LPP, double* LPM, double* LPE, double* LPI, double* LPS, list<string>* segString, const int &T, const int &N);

map<char, int> BASE2ID;
int ALPHABET_SIZE;
double EPSILON = pow(10, -8);
double a, b, g; // parameters for DP
// a: stay in insertion, 1-a: is going into insertion
// b: suffix
// g: deletion 

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
    for (int n=0; n<N; n++){ // extend loop to ad 2 Ns at start and end of read        
        copy(seq + n, seq + n+K, tempKmer);
        kmer_seq[n] = toDeci(tempKmer);
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
inline double scoreKmer(const double &x, const int &kmer, vector<tuple<double, double>>* model) {
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
double logPlus(const double &x, const double &y) {
    // safety check
    if(x==y && isinf(x) && isinf(y)) {
        return x;
    } 
    if(x>=y){
        return x + log1p(exp(y-x));
    }
    return y + log1p(exp(x-y));
}

// inline double m1() {
//     return log(1);
// }

// inline double m3() {
//     return log(1-0.025);
// }

inline double i1() { // insertion
    return log(a); // TODO find good value here
}

// inline double i2() { // extend insertion
//     return log(a); // TODO find good value here
// }

/**
 * Calculte prefix score, RNA read should start with polyA
*/
inline double prefix(const double &signal_dp) {
    // AAAAA mean ~= 108.9, stdev ~= 2.7
    // check 3 stdevs
    if (signal_dp < 100.8 && 117.0 < signal_dp){
        return 0; // TODO parameterise
    }
    return -INFINITY;
}

inline double s1() { // suffix
    return log(b); // TODO find good value here
}

inline double delete(const double &sig, const int &kmer, const int &suckmer, vector<tuple<double, double>>* model) { // deletion
    double mean = get<0>((*model)[kmer]);
    double stdev = get<1>((*model)[kmer]);
    double smean = get<0>((*model)[suckmer]);
    double sstdev = get<1>((*model)[suckmer]);
    // current kmer && successor kmer
    if ((sig<mean-3*stdev || sig>mean+3*stdev) && (sig<smean-3*sstdev || sig>smean+3*sstdev)) {
        return log(g);
    }
    return -INFINITY;
}

inline double error(const double &signal_dp) {
    if (signal_dp <= 50.0 || signal_dp >= 150.0){
        return 0;
    }
    return -INFINITY;
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
void logF(double* sig, int* kmer_seq, double* P, double* M, double* I, double* E, double* S, const int &T, const int &N, vector<tuple<double, double>>* model){
    double mat, ins, ext, pre, suf;
    for(int t=0; t<T; t++){
        for(int n=0; n<N; n++){
            pre=-INFINITY;
            mat=-INFINITY;
            ins=-INFINITY;
            ext=-INFINITY;
            suf=-INFINITY;
            if (t>0 && n>0){
                mat=logPlus(mat, I[(t-1)*N+(n-1)] + log(scoreKmer(sig[t-1], kmer_seq[n-1], model)));
                mat=logPlus(mat, P[(t-1)*N+(n-1)] + log(scoreKmer(sig[t-1], kmer_seq[n-1], model)));
                mat=logPlus(mat, E[(t-1)*N+(n-1)] + log(scoreKmer(sig[t-1], kmer_seq[n-1], model)));

                ext=logPlus(ext, M[(t-1)*N+n] + log(scoreKmer(sig[t-1], kmer_seq[n-1], model)));//E1 first extend
                ext=logPlus(ext, E[(t-1)*N+n] + log(scoreKmer(sig[t-1], kmer_seq[n-1], model)));//E2 extend further
                ext=logPlus(ext, E[(t-1)*N+n] + delete(sig[t-1], kmer_seq[n-2], kmer_seq[n-1], model));//E3 delete/kill
                ext=logPlus(ext, E[(t-1)*N+n] + error(sig[t-1]));//E4 error
            }
            if (t>0) {
                pre=logPlus(pre, P[(t-1)*N+n] + prefix(sig[t-1]));
                suf=logPlus(suf, E[(t-1)*N+n] + s1());//S1
                suf=logPlus(suf, S[(t-1)*N+n] + s1());//S2
            }
            if (n>0){
                ins=logPlus(ins, E[t*N+(n-1)] + i1());//I1 open insertion
                ins=logPlus(ins, I[t*N+(n-1)] + i1());//I2 extend insertion
            }
            if(t==0 && n==0){
                pre = 0; // initialize with log(1)
            }
            P[t*N+n] = pre;
            M[t*N+n] = mat;
            I[t*N+n] = ins;
            E[t*N+n] = ext;
            S[t*N+n] = suf;
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
void logB(double* sig, int* kmer_seq, double* P, double* M, double* I, double* E, double* S, const int &T, const int &N, vector<tuple<double, double>>* model) {
    double mat, ins, ext, pre, suf;
    for(int t=T-1; t>=0; t--){
        for(int n=N-1; n>=0; n--){
            pre=-INFINITY;
            mat=-INFINITY;
            ins=-INFINITY;
            ext=-INFINITY;
            suf=-INFINITY;
            if(t==T-1 && n==N-1){
                suf = 0; // initialize with log(1)
            }
            if(t<T-1 && n<N-1){
                pre=logPlus(pre, M[(t+1)*N+(n+1)] + log(scoreKmer(sig[t], kmer_seq[n], model)));
                ins=logPlus(ins, M[(t+1)*N+(n+1)] + log(scoreKmer(sig[t], kmer_seq[n], model)));
                ext=logPlus(ext, M[(t+1)*N+(n+1)] + log(scoreKmer(sig[t], kmer_seq[n], model)));
            }
            if (t<T-1 && n>0) {
                mat=logPlus(mat, E[(t+1)*N+n] + log(scoreKmer(sig[t], kmer_seq[n-1], model)));//E1 first extend
                ext=logPlus(ext, E[(t+1)*N+n] + log(scoreKmer(sig[t], kmer_seq[n-1], model)));//E2 extend further
                ext=logPlus(ext, E[(t+1)*N+n] + delete(sig[t], kmer_seq[n-2], kmer_seq[n-1], model));//E3 delete/kill
                ext=logPlus(ext, E[(t+1)*N+n] + error(sig[t]));//E4 error
            }
            if(n<N-1){
                ext=logPlus(ext, I[t*N+(n+1)] + i1());//I1 open insertion
                ins=logPlus(ins, I[t*N+(n+1)] + i1());//I2 extend insertion
            }
            if(t<T-1) {
                suf=logPlus(suf, S[(t+1)*N+n] + s1());//S1
                ext=logPlus(ext, S[(t+1)*N+n] + s1());//S2
                pre=logPlus(pre, P[(t+1)*N+n] + prefix(sig[t]));
            }
            P[t*N+n] = pre;
            M[t*N+n] = mat;
            I[t*N+n] = ins;
            E[t*N+n] = ext;
            S[t*N+n] = suf;
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
    for(int t=0; t<T; t++){
        for(int n=0; n<N; n++){
            int x = t*N+n;
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
list<string> getBorders(double* LPP, double* LPM, double* LPE, double*LPI, double* LPS, const int &T, const int &N){
    double* P = new double[T*N];
    double* M = new double[T*N];
    double* E = new double[T*N];
    double* I = new double[T*N];
    double* S = new double[T*N];
    fill_n(P, T*N, -INFINITY);
    fill_n(M, T*N, -INFINITY);
    fill_n(E, T*N, -INFINITY);
    fill_n(I, T*N, -INFINITY);
    fill_n(S, T*N, -INFINITY);

    double pre, mat, ext, ins, suf;
    for(int t=0; t<T; t++){
        for(int n=0; n<N; n++){
            pre=-INFINITY;
            mat=-INFINITY;
            ins=-INFINITY;
            ext=-INFINITY;
            suf=-INFINITY;
            if (t>0 && n>0){
                mat=max(mat, I[(t-1)*N+(n-1)] + LPM[t*N+n]); // M1
                mat=max(mat, P[(t-1)*N+(n-1)] + LPM[t*N+n]); // M2
                mat=max(mat, E[(t-1)*N+(n-1)] + LPM[t*N+n]); // M3

                ext=max(ext, M[(t-1)*N+n] + LPE[t*N+n]); // E1
                ext=max(ext, E[(t-1)*N+n] + LPE[t*N+n]); // E2, E3, E4 // TODO vielleicht spaeter mit eigenem nicht terminal U (undesired) abessen oder nochmal ruber iterieren und markieren
            }
            if (t>0) {
                pre=max(pre, P[(t-1)*N+n] + LPP[t*N+n]); // P
                suf=max(suf, S[(t-1)*N+n] + LPS[t*N+n]); // S1
                suf=max(suf, E[(t-1)*N+n] + LPS[t*N+n]); // S2
            }
            if (n>0){
                ins=max(ins, E[t*N+(n-1)] + LPI[t*N+n]); // I1
                ins=max(ins, I[t*N+(n-1)] + LPI[t*N+n]); // I2
            }
            if(t==0 && n==0){
                pre=LPP[0]; // initialize
            }
            P[t*N+n]=pre;
            M[t*N+n]=mat;
            I[t*N+n]=ins;
            E[t*N+n]=ext;
            S[t*N+n]=suf;
        }
    }
    
    list<string> segString;
    funcS(T-1, N-1, P, M, E, I, S, LPP, LPM, LPE, LPI, LPS, &segString, T, N);
    delete[] P, M, E, I, S;
    return segString;
}

// double Srec(int i, int j){
//     return 
// }

// void funcP(int i, int j, double* P, double* M, double* E, double* I, double* S, double* LPP, double* LPM, double* LPE, double* LPI, double* LPS, list<string>* segString, const int &T, const int &N){
//     double score = P[i*N+j];
//     if (i>0 && score == P[(i-1)*N+j] + LPP[i*N+j]){

//     }
// }

void funcM(int t, int n, double* P, double* M, double* E, double* I, double* S, double* LPP, double* LPM, double* LPE, double* LPI, double* LPS, list<string>* segString, const int &T, const int &N){
    double score = M[t*N+n];
    if (t>0 && n>0 && score == I[(t-1)*N+(n-1)] + LPM[t*N+n]){
        segString->push_front("M"+to_string(n-1)+","+to_string(t-1));
        return funcI(t-1, n-1, P, M, E, I, S, LPP, LPM, LPE, LPI, LPS, segString, T, N);
    }
    if (t>0 && n>0 && score == P[(t-1)*N+(n-1)] + LPM[t*N+n]){
        segString->push_front("M"+to_string(n-1)+","+to_string(t-1));
        return;
        // TODO brauche ich ja nicht mehr checken, weil bleibt in P
        // funcP(i-1, j-1, P, M, E, I, S, LPP, LPM, LPE, LPI, LPS, segString, T, N);
    }
    if (t>0 && n>0 && score == E[(t-1)*N+(n-1)] + LPM[t*N+n]){
        segString->push_front("M"+to_string(n-1)+","+to_string(t-1));
        return funcE(t-1, n-1, P, M, E, I, S, LPP, LPM, LPE, LPI, LPS, segString, T, N);
    }
}

void funcE(int t, int n, double* P, double* M, double* E, double* I, double* S, double* LPP, double* LPM, double* LPE, double* LPI, double* LPS, list<string>* segString, const int &T, const int &N){
    double score = E[t*N+n];
    if (t>0 && score == M[(t-1)*N+n] + LPE[t*N+n]){
        return funcM(t-1, n, P, M, E, I, S, LPP, LPM, LPE, LPI, LPS, segString, T, N);
    }
    if (t>0 && score == E[(t-1)*N+n] + LPE[t*N+n]){
        return funcE(t-1, n, P, M, E, I, S, LPP, LPM, LPE, LPI, LPS, segString, T, N);
    }
}

void funcI(int t, int n, double* P, double* M, double* E, double* I, double* S, double* LPP, double* LPM, double* LPE, double* LPI, double* LPS, list<string>* segString, const int &T, const int &N){
    double score = I[t*N+n];
    if (n>0 && score == I[t*N+(n-1)] + LPI[t*N+n]){
        segString->push_front("I"+to_string(n-1)+","+to_string(t));
        return funcI(t, n-1, P, M, E, I, S, LPP, LPM, LPE, LPI, LPS, segString, T, N);
    }
    if (n>0 && score == E[t*N+(n-1)] + LPI[t*N+n]){
        segString->push_front("I"+to_string(n-1)+","+to_string(t));
        return funcE(t, n-1, P, M, E, I, S, LPP, LPM, LPE, LPI, LPS, segString, T, N);
    }
}

void funcS(int t, int n, double* P, double* M, double* E, double* I, double* S, double* LPP, double* LPM, double* LPE, double* LPI, double* LPS, list<string>* segString, const int &T, const int &N){
    double score = S[t*N+n];
    if (n>0 && score == S[(t-1)*N+n] + LPS[t*N+n]){
        return funcS(t-1, n, P, M, E, I, S, LPP, LPM, LPE, LPI, LPS, segString, T, N);
    }
    if (n>0 && score == E[(t-1)*N+n] + LPS[t*N+n]){
        segString->push_front("S"+to_string(n)+","+to_string(t-1));
        return funcE(t-1, n, P, M, E, I, S, LPP, LPM, LPE, LPI, LPS, segString, T, N);
    }
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
    // Argparser
    argparse::ArgumentParser program("program_name", "0.0.1");
    // 0.00001 is just some starting default value
    program.add_argument("-a", "--alpha").help("Score for insertion").default_value(0.00001); // insertion
    program.add_argument("-b", "--beta").help("Score for suffix").default_value(0.00001); // suffix
    program.add_argument("-g", "--gamma").help("Score for deletion").default_value(0.00001); // deletion
    
    try {
        program.parse_args(argc, argv);
    }
    catch (const std::runtime_error& err) {
        cerr << err.what() << std::endl;
        cerr << program;
        return 1;
    }
    
    a = program.get<double>("alpha");
    b = program.get<double>("beta");
    g = program.get<double>("gamma");

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
        int seq_size = read.size() + (K-1); 
        int* seq = new int[seq_size];
        fill_n(seq, seq_size, 0); // default: fill with A add 2 As to 3' of read
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

        // initialize for matrices
        double* forM = new double[T*N];
        fill_n(forM, T*N, -INFINITY);
        double* forI = new double[T*N];
        fill_n(forI, T*N, -INFINITY);
        double* forP = new double[T*N];
        fill_n(forP, T*N, -INFINITY);
        double* forS = new double[T*N];
        fill_n(forS, T*N, -INFINITY);
        double* forE = new double[T*N];
        fill_n(forE, T*N, -INFINITY);
        // initialize back matrices
        double* backM = new double[T*N];
        fill_n(backM, T*N, -INFINITY);
        double* backI = new double[T*N];
        fill_n(backI, T*N, -INFINITY);
        double* backP = new double[T*N];
        fill_n(backP, T*N, -INFINITY);
        double* backS = new double[T*N];
        fill_n(backS, T*N, -INFINITY);
        double* backE = new double[T*N];
        fill_n(backE, T*N, -INFINITY);

        // calculate segmentation probabilities, fill matrices
        logF(sig, kmer_seq, forP, forM, forI, forE, forS, T, N, &model);
        logB(sig, kmer_seq, backP, backM, backI, backE, backS, T, N, &model);
        
        // cerr.precision(4);
        // cerr<<"backP\n";
        // for(int t=0; t<T; t++){
        //    for(int n=0; n<N; n++){
        //        cerr<<backP[t*N+n]<<", ";
        //    }
        //    cerr<<endl;
        // }
        // cerr<<"backM\n";
        // for(int t=0; t<T; t++){
        //    for(int n=0; n<N; n++){
        //        cerr<<backM[t*N+n]<<", ";
        //    }
        //    cerr<<endl;
        // }cerr<<"backI\n";
        // for(int t=0; t<T; t++){
        //    for(int n=0; n<N; n++){
        //        cerr<<backI[t*N+n]<<", ";
        //    }
        //    cerr<<endl;
        // }cerr<<"backE\n";
        // for(int t=0; t<T; t++){
        //    for(int n=0; n<N; n++){
        //        cerr<<backE[t*N+n]<<", ";
        //    }
        //    cerr<<endl;
        // }cerr<<"backS\n";
        // for(int t=0; t<T; t++){
        //    for(int n=0; n<N; n++){
        //        cerr<<backS[t*N+n]<<", ";
        //    }
        //    cerr<<endl;
        // }

        // cerr<<logPlus(logPlus(logPlus(forM[T*N-1], forI[T*N-1]), forP[T*N-1]), forE[T*N-1])<<endl;
        // cerr<<"forS: "<<forS[T*N-1]<<", forP: "<<forP[T*N-1]<<", forE: "<<forE[T*N-1]<<", forI: "<<forI[T*N-1]<<", forM: "<<forM[T*N-1]<<endl;
        // cerr<<"backS: "<<backS[0]<<", backP: "<<backP[0]<<", backE: "<<backE[0]<<", backI: "<<backI[0]<<", backM: "<<backM[0]<<endl;

        // Check if Z value matches, must match between matrices
        if (fabs(forS[N*T-1] - backP[0])>EPSILON) {
            cerr.precision(20);
            cerr<<"Z values between matrices do not match!\n";
            cerr<<"forS: "<<forS[N*T-1]<<", backP: "<<backP[0]<<"\n";
            cerr<<backP[0]-forS[N*T-1]<<"\n";
            cerr<<fabs(forS[N*T-1] - backP[0])<<" < "<<EPSILON<<endl;
            cerr.flush();
            exit(11);
        }

        double* LPP = logP(forP, backP, backP[0], T, N); // log probs
        double* LPM = logP(forM, backM, backP[0], T, N); // log probs
        double* LPE = logP(forE, backE, backP[0], T, N); // log probs
        double* LPI = logP(forI, backI, backP[0], T, N); // log probs
        double* LPS = logP(forS, backS, backP[0], T, N); // log probs
        list<string> segString = getBorders(LPP, LPM, LPE, LPI, LPS, T, N);

        // cout<<"Segmentation"<<endl;
        for (auto const& seg : segString) {
            cout<<seg;
        }
        cout<<endl;
        cout.flush();

        // calculate sum of segment probabilities
        double* LSP = new double[T]; // log segment probs
        fill_n(LSP, T, -INFINITY);
        int idx;
        double sum;
        for(int t=0;t<T;t++) {
            sum = -INFINITY;
            for(int n=0; n<N; n++) {
                idx = t*N+n;
                sum = logPlus(sum, LPM[idx]);
            }
            LSP[i] = sum;
            cout<<sum<<",";
        }
        cout<<endl;
        cout.flush();

        // Clean up
        delete[] forP, forM, forI, forE, backS, backP, backM, backI, backE, backS, sig, seq, kmer_seq, LPP, LPM, LPE, LPI, LPS, LSP;
    }
    return 0;
}
