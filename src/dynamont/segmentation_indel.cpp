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

void funcM(int i, int j, double* M, double* E, double* I, double* D, double* LPM, double* LPE, double* LPI, double* LPD, list<string>* segString, const int &N);
void funcE(int i, int j, double* M, double* E, double* I, double* D, double* LPM, double* LPE, double* LPI, double* LPD, list<string>* segString, const int &N);
void funcI(int i, int j, double* M, double* E, double* I, double* D, double* LPM, double* LPE, double* LPI, double* LPD, list<string>* segString, const int &N);
void funcD(int i, int j, double* M, double* E, double* I, double* D, double* LPM, double* LPE, double* LPI, double* LPD, list<string>* segString, const int &N);

map<char, int> BASE2ID;
string modelpath;
int ALPHABET_SIZE;
int pore;
double EPSILON = pow(10, -5);
bool atrain, calcZ, train;
double m2, m3, m4, e1, e2, e3, i1, i2, d1, d2; // transition parameters
int K; // our model works with this kmer size

// TODO move to config file?
// const string MODELPATH = "/home/yi98suv/projects/ont_segmentation/data/template_median69pA.model";
// const string MODELPATH = "/home/yi98suv/projects/dynamont/data/template_median69pA_extended.model";
// const string MODELPATH = "/home/yi98suv/projects/dynamont/data/template_median69pA_reduced.model";
const string MODELPATH = "/home/yi98suv/projects/dynamont/data/template_median69pA_reduced_halfed.model";
const string TERM_STRING = "$";

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
    BASE2ID.insert(pair<char, int>('U', 3));
    BASE2ID.insert(pair<char, int>('N', 4));
    BASE2ID.insert(pair<char, int>('a', 0));
    BASE2ID.insert(pair<char, int>('c', 1));
    BASE2ID.insert(pair<char, int>('g', 2));
    BASE2ID.insert(pair<char, int>('t', 3));
    BASE2ID.insert(pair<char, int>('u', 3));
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
 * @param signal point to calculate probability density
 * @param kmer key for the model kmer:(mean, stdev) map
 * @param model map containing kmers as keys and (mean, stdev) tuples as values
 * @return probability density value for x in the given normal distribution
 */
inline double scoreKmer(const double &signal, const int &kmer, vector<tuple<double, double>>* model) {
    tuple<double, double> kmerModel = (*model)[kmer];
    return 10*(log(normal_pdf(signal, get<0>(kmerModel), get<1>(kmerModel))) + 2);

    // norm signal with kmer model
    // double sig = (signal - get<0>(kmerModel)) / get<1>(kmerModel);
    // return 2*(log(normal_pdf(sig, 0.0, 1.0)) + 6);
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
    if (x==y && isinf(x) && isinf(y)) {
        return x;
    } 
    if (x>=y){
        return x + log1p(exp(y-x));
    }
    return y + log1p(exp(x-y));
}

/**
 * Calculte prefix score, RNA read should start with polyA
*/
// inline double prefix(const double &signal, vector<tuple<double, double>>* model) {
//     int kmer = 0;
//     tuple<double, double> kmerModel = (*model)[kmer];
//     return normal_pdf(signal, get<0>(kmerModel), get<1>(kmerModel));
// }

// inline double indel(const double &sig, const int &kmer, const int &suckmer, vector<tuple<double, double>>* model) {
inline double indel(const double &sig, const int &kmer, const int &suckmer, vector<tuple<double, double>>* model) {
    // double mean = get<0>((*model)[kmer]);
    // double stdev = get<1>((*model)[kmer]);
    // double smean = get<0>((*model)[suckmer]);
    // double sstdev = get<1>((*model)[suckmer]);
    // current kmer && successing kmer
    // if ((sig<mean-5*stdev || sig>mean+5*stdev) && (sig<smean-5*sstdev || sig>smean+5*sstdev)) {
        // return log(normal_pdf(mean, mean, stdev)); // TODO? parameterize
    return -1000;
    // }
    // return log(0);
}

// inline double deletion(const double &sig, const int &kmer, const int &suckmer, vector<tuple<double, double>>* model) {
//     double mean = get<0>((*model)[kmer]);
//     double stdev = get<1>((*model)[kmer]);
//     double smean = get<0>((*model)[suckmer]);
//     double sstdev = get<1>((*model)[suckmer]);
//     // current kmer && successing kmer
//     if ((sig<mean-3*stdev || sig>mean+3*stdev) && (sig>smean-3*sstdev || sig<smean+3*sstdev)) {
//         return log(normal_pdf(mean, mean, stdev)); // TODO? parameterize
//     }
//     return -1000;
// }

// inline double indel(const double &signal, const int &kmer, const int &suckmer, vector<tuple<double, double>>* model) {
//     double mean = get<0>((*model)[kmer]);
//     double stdev = get<1>((*model)[kmer]);
//     double smean = get<0>((*model)[suckmer]);
//     double sstdev = get<1>((*model)[suckmer]);
//     // current kmer && successing kmer
//     return (normal_pdf(mean, mean, stdev) - normal_pdf(signal, mean, stdev)) * (normal_pdf(smean, smean, sstdev) - normal_pdf(signal, smean, sstdev));
//     }

inline double error(const double &signal_dp) {
    if (signal_dp <= 50.0 || signal_dp >= 150.0){
        return log(1);
    }
    return -1000;
}

/**
 * Calculate forward matrices using logarithmic values
 *
 * @param sig ONT raw signal with pA values
 * @param seq nucleotide sequence represented by the ONT signal
 * @param T length of the ONT raw signal + 1
 * @param N length of nucleotide sequence + 1
 * @param model map containing kmers as keys and (mean, stdev) tuples as values
 */
void logF(double* sig, int* kmer_seq, double* M, double* I, double* D, double* E, const int &T, const int &N, vector<tuple<double, double>>* model){
    double mat, ins, del, ext;
    for(int t=0; t<T; t++){
        for(int n=0; n<N; n++){
            mat=-INFINITY;
            ins=-INFINITY;
            del=-INFINITY;
            ext=-INFINITY;
            if (t>0 && n>0){
                mat=logPlus(mat, E[(t-1)*N+(n-1)] + scoreKmer(sig[t-1], kmer_seq[n-1], model) + m2); // m2
                mat=logPlus(mat, I[(t-1)*N+(n-1)] + scoreKmer(sig[t-1], kmer_seq[n-1], model) + m3); // m3
                mat=logPlus(mat, D[(t-1)*N+(n-1)] + scoreKmer(sig[t-1], kmer_seq[n-1], model) + m4); // m4

                ext=logPlus(ext, M[(t-1)*N+n] + scoreKmer(sig[t-1], kmer_seq[n-1], model) + e1); // e1 first extend
                ext=logPlus(ext, E[(t-1)*N+n] + scoreKmer(sig[t-1], kmer_seq[n-1], model) + e2); // e2 extend further
                ext=logPlus(ext, E[(t-1)*N+n] + error(sig[t-1]) + e3); // e3 error
            }
            if (n<N-1 && t>0 && n>0) { // kmer_seq has size N-1: [0 1 ... N-3 N-2]
                del=logPlus(del, E[(t-1)*N+n] + indel(sig[t-1], kmer_seq[n-1], kmer_seq[n], model) + d1); // d1 open deletion
                del=logPlus(del, D[(t-1)*N+n] + indel(sig[t-1], kmer_seq[n-1], kmer_seq[n], model) + d2); // d2 extend deletion
                ins=logPlus(ins, E[t*N+(n-1)] + indel(sig[t-1], kmer_seq[n-1], kmer_seq[n], model) + i1); // i1 open insertion
                ins=logPlus(ins, I[t*N+(n-1)] + indel(sig[t-1], kmer_seq[n-1], kmer_seq[n], model) + i2); // i2 extend insertion
            }
            if (t==0 && n==1){
                mat = 0; // initialize with log(1)
            }
            M[t*N+n] = mat;
            I[t*N+n] = ins;
            D[t*N+n] = del;
            E[t*N+n] = ext;
        }
    }
}

// ext=logPlus(ext, E[(t-1)*N+n] + deletion(sig[t-1], kmer_seq[n-2], kmer_seq[n-1], model) + e3);//E3 delete/deletion

/**
 * Calculate backward matrices using logarithmic values
 *
 * @param sig ONT raw signal with pA values
 * @param seq nucleotide sequence represented by the ONT signal
 * @param T length of the ONT raw signal + 1
 * @param N length of nucleotide sequence + 1
 * @param model map containing kmers as keys and (mean, stdev) tuples as values
 */
void logB(double* sig, int* kmer_seq, double* M, double* I, double* D, double* E, const int &T, const int &N, vector<tuple<double, double>>* model) {
    double mat, ins, del, ext;
    for(int t=T-1; t>=0; t--){
        for(int n=N-1; n>=0; n--){
            mat=-INFINITY;
            ins=-INFINITY;
            del=-INFINITY;
            ext=-INFINITY;
            if (t==T-1 && n==N-1){
                ext = 0; // initialize with log(1)
            }
            if (t<T-1 && n<N-1) {
                ext=logPlus(ext, M[(t+1)*N+(n+1)] + scoreKmer(sig[t], kmer_seq[n], model) + m2); // m2
                ins=logPlus(ins, M[(t+1)*N+(n+1)] + scoreKmer(sig[t], kmer_seq[n], model) + m3); // m3
                del=logPlus(del, M[(t+1)*N+(n+1)] + scoreKmer(sig[t], kmer_seq[n], model) + m4); // m4
            }
            if (t<T-1 && n>0) {
                mat=logPlus(mat, E[(t+1)*N+n] + scoreKmer(sig[t], kmer_seq[n-1], model) + e1); // e1 first extend
                ext=logPlus(ext, E[(t+1)*N+n] + scoreKmer(sig[t], kmer_seq[n-1], model) + e2); // e2 extend further
                ext=logPlus(ext, E[(t+1)*N+n] + error(sig[t]) + e3); // e3 error
            }
            if (t<T-1 && n>0 && n<N-1) {
                ext=logPlus(ext, D[(t+1)*N+n] + indel(sig[t], kmer_seq[n-1], kmer_seq[n], model) + d1); // d1 deletion
                del=logPlus(del, D[(t+1)*N+n] + indel(sig[t], kmer_seq[n-1], kmer_seq[n], model) + d2); // d2 deletion
            }
            if (t>0 && n<N-2){
                ext=logPlus(ext, I[t*N+(n+1)] + indel(sig[t-1], kmer_seq[n], kmer_seq[n+1], model) + i1); // i1 open insertion
                ins=logPlus(ins, I[t*N+(n+1)] + indel(sig[t-1], kmer_seq[n], kmer_seq[n+1], model) + i2); // i2 extend insertion
            }
            M[t*N+n] = mat;
            I[t*N+n] = ins;
            D[t*N+n] = del;
            E[t*N+n] = ext;
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
list<string> getBorders(double* LPM, double* LPE, double* LPI, double* LPD, const int &T, const int &N){
    double* M = new double[T*N];
    double* E = new double[T*N];
    double* I = new double[T*N];
    double* D = new double[T*N];
    fill_n(M, T*N, -INFINITY);
    fill_n(E, T*N, -INFINITY);
    fill_n(I, T*N, -INFINITY);
    fill_n(D, T*N, -INFINITY);
    double mat, ext, ins, del;
    for(int t=0; t<T; t++){
        for(int n=0; n<N; n++){
            mat=-INFINITY;
            ins=-INFINITY;
            del=-INFINITY;
            ext=-INFINITY;
            if (t>0 && n>0){
                mat=max(mat, I[(t-1)*N+(n-1)] + LPM[t*N+n]); // m1
                mat=max(mat, E[(t-1)*N+(n-1)] + LPM[t*N+n]); // m3
                mat=max(mat, D[(t-1)*N+(n-1)] + LPM[t*N+n]); // m4
                ext=max(ext, M[(t-1)*N+n] + LPE[t*N+n]); // e1
                ext=max(ext, E[(t-1)*N+n] + LPE[t*N+n]); // e2, e3
            }
            if (n<N-1 && t>0 && n>0) {
                del=max(del, E[(t-1)*N+n] + LPD[t*N+n]); // d1 open deletion
                del=max(del, D[(t-1)*N+n] + LPD[t*N+n]); // d2 extend deletion
                ins=max(ins, E[t*N+(n-1)] + LPI[t*N+n]); // i1
                ins=max(ins, I[t*N+(n-1)] + LPI[t*N+n]); // i2
            }
            if (t==0 && n==1){
                mat=LPM[1]; // initialize
            }
            M[t*N+n]=mat;
            I[t*N+n]=ins;
            D[t*N+n]=del;
            E[t*N+n]=ext;
        }
    }
    
    list<string> segString;
    funcE(T-1, N-1, M, E, I, D, LPM, LPE, LPI, LPD, &segString, N);
    delete[] M, E, I, D;
    return segString;
}

void funcM(int t, int n, double* M, double* E, double* I, double* D, double* LPM, double* LPE, double* LPI, double* LPD, list<string>* segString, const int &N){
    double score = M[t*N+n];
    if (t==0 && n==1){ // Start value
        segString->push_front("M"+to_string(n-1)+","+to_string(0)); // n-1 because N is 1 larger than the sequences
        return;
    }
    if (t>0 && n>0 && score == I[(t-1)*N+(n-1)] + LPM[t*N+n]){
        segString->push_front("M"+to_string(n-1)+","+to_string(t-1)); // n-1, t-1 because T and N are 1 larger than the sequences
        return funcI(t-1, n-1, M, E, I, D, LPM, LPE, LPI, LPD, segString, N);
    }
    if (t>0 && n>0 && score == E[(t-1)*N+(n-1)] + LPM[t*N+n]){
        segString->push_front("M"+to_string(n-1)+","+to_string(t-1));
        return funcE(t-1, n-1, M, E, I, D, LPM, LPE, LPI, LPD, segString, N);
    }
    if (t>0 && n>0 && score == D[(t-1)*N+(n-1)] + LPM[t*N+n]){
        segString->push_front("M"+to_string(n-1)+","+to_string(t-1));
        return funcD(t-1, n-1, M, E, I, D, LPM, LPE, LPI, LPD, segString, N);
    }
}

void funcE(int t, int n, double* M, double* E, double* I, double* D, double* LPM, double* LPE, double* LPI, double* LPD, list<string>* segString, const int &N){
    double score = E[t*N+n];
    if (t>0 && score == M[(t-1)*N+n] + LPE[t*N+n]){
        return funcM(t-1, n, M, E, I, D, LPM, LPE, LPI, LPD, segString, N);
    }
    if (t>0 && score == E[(t-1)*N+n] + LPE[t*N+n]){
        return funcE(t-1, n, M, E, I, D, LPM, LPE, LPI, LPD, segString, N);
    }
}

void funcI(int t, int n, double* M, double* E, double* I, double* D, double* LPM, double* LPE, double* LPI, double* LPD, list<string>* segString, const int &N){
    double score = I[t*N+n];
    if (n>0 && score == I[t*N+(n-1)] + LPI[t*N+n]){
        // segString->push_front("I"+to_string(n-1)+","+to_string(t));
        return funcI(t, n-1, M, E, I, D, LPM, LPE, LPI, LPD, segString, N);
    }
    if (n>0 && score == E[t*N+(n-1)] + LPI[t*N+n]){
        segString->push_front("I"+to_string(n-1)+","+to_string(t));
        return funcE(t, n-1, M, E, I, D, LPM, LPE, LPI, LPD, segString, N);
    }
}

void funcD(int t, int n, double* M, double* E, double* I, double* D, double* LPM, double* LPE, double* LPI, double* LPD, list<string>* segString, const int &N){
    double score = D[t*N+n];
    if (t>0 && score == D[(t-1)*N+n] + LPD[t*N+n]){
        // segString->push_front("D"+to_string(n)+","+to_string(t-1));
        return funcD(t-1, n, M, E, I, D, LPM, LPE, LPI, LPD, segString, N);
    }
    if (t>0 && score == E[(t-1)*N+n] + LPD[t*N+n]){
        segString->push_front("D"+to_string(n)+","+to_string(t-1));
        return funcE(t-1, n, M, E, I, D, LPM, LPE, LPI, LPD, segString, N);
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
        getline(buffer, tmp, '\t'); // level_mean
        mean = atof(tmp.c_str());
        getline(buffer, tmp, '\t'); // level_stdv
        stdev = atof(tmp.c_str());
        (*model)[kmer2int(kmer)]=make_tuple(mean, stdev);
    }
    inputFile.close();
}

/**
 * Train parameter
*/
tuple<double, double, double, double, double, double, double, double, double, double> trainBaumWelch(double* sig, int* kmer_seq, double* forM, double* forE, double* forI, double* forD, double* backM, double* backE, double* backI, double* backD, double Z, const int &T, const int &N, vector<tuple<double, double>>* model) {
    // double m2, m3, m4, e1, e2, e3, i1, i2, d1, d2; // parameters for DP
    double newM2 = -INFINITY;
    double newM3 = -INFINITY;
    double newM4 = -INFINITY;
    double newE1 = -INFINITY;
    double newE2 = -INFINITY;
    double newE3 = -INFINITY;
    double newD1 = -INFINITY;
    double newD2 = -INFINITY;
    double newI1 = -INFINITY;
    double newI2 = -INFINITY;
    // double A = log(T*N); // number of transitions
    double countM = 0;

    for(int t=0; t<T; t++){
        for(int n=0; n<N; n++){
            if (t<T-1 && n<N-1) {
                // m1:                 forward(i)     a    e(i+1)                                                backward(i+1)
                newM2 = logPlus(newM2, (forE[t*N+n] + m2 + scoreKmer(sig[t], kmer_seq[n], model)            + backM[(t+1)*N+(n+1)]));
                newM3 = logPlus(newM3, (forI[t*N+n] + m3 + scoreKmer(sig[t], kmer_seq[n], model)            + backM[(t+1)*N+(n+1)]));
                newM4 = logPlus(newM4, (forD[t*N+n] + m4 + scoreKmer(sig[t], kmer_seq[n], model)            + backM[(t+1)*N+(n+1)]));
            }
            if (t<T-1 && n>0) {
                newE1 = logPlus(newE1, (forM[t*N+n] + e1 + scoreKmer(sig[t], kmer_seq[n-1], model)          + backE[(t+1)*N+n]));
                newE2 = logPlus(newE2, (forE[t*N+n] + e2 + scoreKmer(sig[t], kmer_seq[n-1], model)          + backE[(t+1)*N+n]));
                if (n<N-1) {
                    newD1 = logPlus(newD1, (forE[t*N+n] + d1 + indel(sig[t], kmer_seq[n-1], kmer_seq[n], model)  + backD[(t+1)*N+n]));
                    newD2 = logPlus(newD2, (forD[t*N+n] + d2 + indel(sig[t], kmer_seq[n-1], kmer_seq[n], model)  + backD[(t+1)*N+n]));
                }
            }
            if (t<T-1) {
                newE3 = logPlus(newE3, (forE[t*N+n] + e3 + error(sig[t])                                         + backE[(t+1)*N+n]));
            }
            if (n<N-2 && t>0){
                newI1 = logPlus(newI1, (forE[t*N+n] + i1 + indel(sig[t-1], kmer_seq[n], kmer_seq[n+1], model)    + backI[t*N+(n+1)]));
                newI2 = logPlus(newI2, (forI[t*N+n] + i2 + indel(sig[t-1], kmer_seq[n], kmer_seq[n+1], model)    + backI[t*N+(n+1)]));
            }  
        }
    }
    // divide by Z
    // average over the number of transitions
    newM2 = newM2 - Z;
    newM3 = newM3 - Z;
    newM4 = newM4 - Z;
    newE1 = newE1 - Z;
    newE2 = newE2 - Z;
    newE3 = newE3 - Z;
    newD1 = newD1 - Z;
    newD2 = newD2 - Z;
    newI1 = newI1 - Z;
    newI2 = newI2 - Z;
    double Am = newE1;
    double Ae = logPlus(newD1, logPlus(newE2, logPlus(newE3, logPlus(newI1, newM2))));
    double Ad = logPlus(newM4, newD2);
    double Ai = logPlus(newM3, newI2);
    newM2 = newM2 - Ae;
    newM3 = newM3 - Ai;
    newM4 = newM4 - Ad;
    newE1 = newE1 - Am;
    newE2 = newE2 - Ae;
    newE3 = newE3 - Ae;
    newD1 = newD1 - Ae;
    newD2 = newD2 - Ad;
    newI1 = newI1 - Ae;
    newI2 = newI2 - Ai;
    return {newM2, newM3, newM4, newE1, newE2, newE3, newD1, newD2, newI1, newI2};
}

void printTrainedTransitionParams(double* sig, int* kmer_seq, double* forM, double* forE, double* forI, double* forD, double* backM, double* backE, double* backI, double* backD, const int &T, const int &N, vector<tuple<double, double>>* model) {
    auto [newM2, newM3, newM4, newE1, newE2, newE3, newD1, newD2, newI1, newI2] = trainBaumWelch(sig, kmer_seq, forM, forE, forI, forD, backM, backE, backI, backD, backM[1], T, N, model);

            double newParams[16] = {newM2, newM3, newM4, newE1, newE2, newE3, newD1, newD2, newI1, newI2};
            double sumLogs = -INFINITY;
            for (int i = 0; i < 16; i++) {
                sumLogs = logPlus(sumLogs, newParams[i]);
            }

            cout<<"m2:"<<exp(newM2)<<";m3:"<<exp(newM3)<<";m4:"<<exp(newM4)<<";e1:"<<exp(newE1)<<";e2:"<<exp(newE2)<<";e3:"<<exp(newE3)<<";d1:"<<exp(newD1)<<";d2:"<<exp(newD2)<<";i1:"<<exp(newI1)<<";i2:"<<exp(newI2)<<endl;
            cout<<"Z:"<<backM[1]<<endl;
            cout.flush();
}

/**
 * Read signal and read from stdin until the TERM_STRING is seen
*/
int main(int argc, char* argv[]) {
    // Argparser
    argparse::ArgumentParser program("program_name", "0.0.2");
    // 0.00001 is just some starting default value
    // double m2, m3, m4, e1, e2, e3, i1, i2, d, s1, s2; // parameters for DP
    program.add_argument("-e1", "--extendscore1").help("Transition probability for extend rule 1").default_value(    1.00).scan<'g', double>(); // e1
    program.add_argument("-m2", "--matchscore2").help("Transition probability for match rule 2").default_value(      0.16).scan<'g', double>(); // m2
    program.add_argument("-d1", "--deletescore1").help("Transition probability for deletion rule 1").default_value( 0.16).scan<'g', double>(); // d1
    program.add_argument("-e2", "--extendscore2").help("Transition probability for extend rule 2").default_value(    0.16).scan<'g', double>(); // e2
    program.add_argument("-e3", "--extendscore3").help("Transition probability for extend rule 3").default_value(    0.16).scan<'g', double>(); // e3
    program.add_argument("-i1", "--insertscore1").help("Transition probability for insertion rule 1").default_value( 0.16).scan<'g', double>(); // i1
    program.add_argument("-m3", "--matchscore3").help("Transition probability for match rule 3").default_value(      0.50).scan<'g', double>(); // m3
    program.add_argument("-i2", "--insertscore2").help("Transition probability for insertion rule 2").default_value( 0.50).scan<'g', double>(); // i2
    program.add_argument("-m4", "--matchscore4").help("Transition probability for match rule 4").default_value(      0.50).scan<'g', double>(); // m4
    program.add_argument("-d2", "--deletescore2").help("Transition probability for deletion rule 2").default_value(  0.50).scan<'g', double>(); // d2
    program.add_argument("-at", "--atrain").help("Switch algorithm to transition parameter training mode").default_value(false).implicit_value(true);
    program.add_argument("-t", "--train").help("Switch algorithm to transition and emission parameter training mode").default_value(false).implicit_value(true);
    program.add_argument("-z", "--calcZ").help("Switch algorithm to only calculate Z").default_value(false).implicit_value(true);
    program.add_argument("-m", "--model").help("Path to kmer model table").default_value(MODELPATH);
    program.add_argument("-r", "--pore").help("Pore generation used to sequence the data").default_value(9).choices(9, 10);
    
    try {
        program.parse_args(argc, argv);
    }
    catch (const std::runtime_error& err) {
        cerr << err.what() << std::endl;
        cerr << program;
        return 1;
    }
    
    pore = program.get<int>("pore");
    modelpath = program.get<string>("model");
    atrain = program.get<bool>("atrain");
    train = program.get<bool>("train");
    calcZ = program.get<bool>("calcZ");
    m2 = log(program.get<double>("matchscore2"));
    m3 = log(program.get<double>("matchscore3"));
    m4 = log(program.get<double>("matchscore4"));
    e1 = log(program.get<double>("extendscore1"));
    e2 = log(program.get<double>("extendscore2"));
    e3 = log(program.get<double>("extendscore3"));
    i1 = log(program.get<double>("insertscore1"));
    i2 = log(program.get<double>("insertscore2"));
    d1 = log(program.get<double>("deletescore1"));
    d2 = log(program.get<double>("deletescore2"));

    if (pore == 9) {
        K = 5;
    } else if (pore == 10) {
        K = 9;
    }
    // cerr<<p<<", "<<e1<<endl;
    // cerr.flush();
    
    // cout<<"START\n";
    // read data from stdin

    fillBASE2ID();
    vector<tuple<double, double>> model(pow(ALPHABET_SIZE, K), make_tuple(-INFINITY, -INFINITY));
    readKmerModel(modelpath, &model);
    string input;
    string signal;
    string read;
    int truish = 1;

    while(truish) {
        // truish = 0;
        // echo 107,107,107.2,108.0,108.9,111.2,105.7,104.3,107.1,105.7,105,105 CAAAAA| src\segment.exe
        // read input, signal and read whitespace separated in single line
        getline(cin, signal, ' ');
        getline(cin, read);

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
        double* forD = new double[T*N];
        fill_n(forD, T*N, -INFINITY);
        double* forE = new double[T*N];
        fill_n(forE, T*N, -INFINITY);
        // initialize back matrices
        double* backM = new double[T*N];
        fill_n(backM, T*N, -INFINITY);
        double* backI = new double[T*N];
        fill_n(backI, T*N, -INFINITY);
        double* backD = new double[T*N];
        fill_n(backD, T*N, -INFINITY);
        double* backE = new double[T*N];
        fill_n(backE, T*N, -INFINITY);

        // calculate segmentation probabilities, fill matrices
        logF(sig, kmer_seq, forM, forI, forD, forE, T, N, &model);
        logB(sig, kmer_seq, backM, backI, backD, backE, T, N, &model);

        // Check if Z value matches, must match between matrices
        cerr.precision(20);
        if (fabs(forE[N*T-1] - backM[1])>EPSILON) {
            cerr<<"Z values between matrices do not match! forE: "<<forE[N*T-1]<<", backP: "<<backM[1]<<"\n";
            cerr.flush();
            exit(11);
        }
        if (calcZ){
            cout<<backM[1]<<"\n";

        // train only transitions
        } else if (atrain) {
            printTrainedTransitionParams(sig, kmer_seq, forM, forE, forI, forD, backM, backE, backI, backD, T, N, &model);

        } else {

            // train both Transitions and Emissions
            if (train) {
                printTrainedTransitionParams(sig, kmer_seq, forM, forE, forI, forD, backM, backE, backI, backD, T, N, &model);
            }
            
            // return segmentations for Emission training
            double* LPM = logP(forM, backM, backM[1], T, N); // log probs
            double* LPE = logP(forE, backE, backM[1], T, N); // log probs
            double* LPI = logP(forI, backI, backM[1], T, N); // log probs
            double* LPD = logP(forD, backD, backM[1], T, N); // log probs
            list<string> segString = getBorders(LPM, LPE, LPI, LPD, T, N);

            // cout<<"Segmentation"<<endl;
            for (auto const& seg : segString) {
                cout<<seg;
            }
            cout<<endl;
            cout.flush();

            // calculate sum of segment probabilities
            // double* LSP = new double[T]; // log segment probs
            // fill_n(LSP, T, -INFINITY);
            // int idx;
            // double sum;
            // for(int t=0;t<T;t++) {
            //     sum = -INFINITY;
            //     for(int n=0; n<N; n++) {
            //         idx = t*N+n;
            //         sum = logPlus(sum, LPM[idx]);
            //     }
            //     LSP[i] = sum;
            //     cout<<sum<<",";
            // }
            // cout<<endl;
            // cout.flush();

            // Clean up
            delete[] LPM, LPE, LPI, LPD; //, LSP;
        }
        delete[] forM, forI, forD, forE, backM, backI, backD, backE, sig, seq, kmer_seq;
    }
    return 0;
}