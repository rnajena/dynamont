// author: Jannes Spangenberg
// e-mail: jannes.spangenberg@uni-jena.de
// github: https://github.com/JannesSP
// website: https://jannessp.github.io

#include <iostream>
#include <iomanip>
#include <fstream> // file io
#include <sstream> // file io
#include <string>
#include <map> // dictionary
#include <tuple>
#include <bits/stdc++.h> // reverse strings
#include <vector>
#include <cmath> // exp
// #include <limits> // for inifinity
#include <assert.h>
#include <stdlib.h>
#include "argparse.hpp"
#include <algorithm>

using namespace std;

void funcM(int t, int n, double* M, double* E, double* D, double* I, double* LPM, double* LPE, double* LPD, double* LPI, list<string>* segString, const int &N);
void funcE(int t, int n, double* M, double* E, double* D, double* I, double* LPM, double* LPE, double* LPD, double* LPI, list<string>* segString, const int &N);
void funcD(int t, int n, double* M, double* E, double* D, double* I, double* LPM, double* LPE, double* LPD, double* LPI, list<string>* segString, const int &N);
void funcI(int t, int n, double* M, double* E, double* D, double* I, double* LPM, double* LPE, double* LPD, double* LPI, list<string>* segString, const int &N);

map<char, int> BASE2ID;
map<char, int> ID2BASE;
string modelpath;
int ALPHABET_SIZE;
double EPSILON = pow(10, -2);
bool train, calcZ; // atrain
double m1, e1, e2, e3, d1, d2, m2, i1, i2, m3; // transition parameters
int K; // our model works with this kmer size
int C;
int numKmers;

const string MODELPATH = "/home/yi98suv/projects/dynamont/data/template_median69pA_extended.model";
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
    ID2BASE.insert(pair<char, int>('0', 'A'));
    ID2BASE.insert(pair<char, int>('1', 'C'));
    ID2BASE.insert(pair<char, int>('2', 'G'));
    ID2BASE.insert(pair<char, int>('3', 'T'));
    ID2BASE.insert(pair<char, int>('4', 'N'));
}

/**
 *
 * C++ version 0.4 std::string style "itoa":
 * Contributions from Stuart Lowe, Ray-Yuan Sheu,
 * Rodrigo de Salvo Braz, Luc Gallant, John Maloney
 * and Brian Hunt
 * 
 * Converts a decimal to number to a number of base ALPHABET_SIZE.
 * TODO Works for base between 2 and 16 (included)
 * 
 * @param value input number in decimal to convert to base
*/
string itoa(int value) {
    string buf;
    int base = K;

    // check that the base if valid
    if (base < 2 || base > 16) return to_string(value);

    enum { kMaxDigits = 35 };
    buf.reserve( kMaxDigits ); // Pre-allocate enough space.
    int quotient = value;

    // Translating number to string with base:
    do {
        buf += ID2BASE.at("0123456789abcdef"[ abs( quotient % base ) ]);
        quotient /= base;
    } while ( quotient );

    // Append the negative sign
    if ( value < 0) buf += '-';

    while ((int) buf.length() < base) {
        buf += ID2BASE.at('0');
    }

    reverse( buf.begin(), buf.end() );
    return buf;
}

/**
 * Converts a number of base ALPHABET_SIZE to a decimal number.
 * TODO Works ONLY if ALPHABET_SIZE is smaller or equal to 10!
 * 
 * @param value input number in the given base
*/
int toDeci(int value) {
    int ret = 0, r = 0;
    int m = 1;
    while(value > 0) {
        r = value % 10;
        ret += m*r;
        m *= ALPHABET_SIZE;
        value /= 10;
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
 * Converts the kmers of the model file to the integer representation using the BASE2ID map
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
    delete[] tempKmer;
    return kmer_seq;
}

// https://ethz.ch/content/dam/ethz/special-interest/mavt/dynamic-systems-n-control/idsc-dam/Lectures/Stochastic-Systems/Statistical_Methods.pdf
/**
 * Calculate log pdf for a given x, mean and standard deviation
 * 
 * @param x value
 * @param m mean
 * @param s standard deviation 
 * @return probabily density at position x for N~(m, s²)
*/
double log_normal_pdf(const double &x, const double &m, const double &s) {
    if(s==0) {
        return -INFINITY;
    }
    return -0.5*(log(2*M_PI*s*s)+((x - m)*(x - m)/(s*s)));
}

/**
 * Return log probability density for a given value and a given normal distribution
 *
 * @param signal point to calculate probability density
 * @param kmer key for the model kmer:(mean, stdev) map
 * @param model map containing kmers as keys and (mean, stdev) tuples as values
 * @return log probability density value for x in the given normal distribution
 */
inline double scoreKmer(const double &signal, const int &kmer, vector<tuple<double, double>>* model) {
    tuple<double, double> kmerModel = (*model)[kmer];
    // return 4.*(log_normal_pdf(signal, get<0>(kmerModel), 1) + 4.);
    // return 4*(log_normal_pdf(signal, get<0>(kmerModel), get<1>(kmerModel)) + 4);
    return 2*(log_normal_pdf(signal, get<0>(kmerModel), get<1>(kmerModel)) + 6);
    // return log_normal_pdf(signal, get<0>(kmerModel), get<1>(kmerModel));

    // norm signal with kmer model
    // double sig = (signal - get<0>(kmerModel)) / get<1>(kmerModel);
    // return 2*(log_normal_pdf(sig, 0.0, 1.0) + 6);
}

// https://en.wikipedia.org/wiki/Log_probability
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

inline double error(const double &signal_dp) {
    if (signal_dp <= 40.0 || signal_dp >= 160.0){
        return 0;
    }
    return -INFINITY;
}

inline double indel(const double &signal_dp, const int &kmer, const int &sucKmer, vector<tuple<double, double>>* model) {
    tuple<double, double> kmerModel = (*model)[kmer];
    tuple<double, double> sucKmerModel = (*model)[sucKmer];
    bool aboveCurDist = (signal_dp >= get<0>(kmerModel) + 3*get<1>(kmerModel));
    bool belowCurDist = (signal_dp <= get<0>(kmerModel) - 3*get<1>(kmerModel));
    bool aboveNextDist = (signal_dp >= get<0>(sucKmerModel) + 3*get<1>(sucKmerModel));
    bool belowNextDist = (signal_dp <= get<0>(sucKmerModel) - 3*get<1>(sucKmerModel));
    
    if ((aboveCurDist || belowCurDist) && (aboveNextDist || belowNextDist)) {
        return scoreKmer(get<0>(kmerModel), kmer, model);
    }
    return -INFINITY;
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
void logF(double* sig, int* kmer_seq, double* M, double* E, double* D, double* I, const int &T, const int &N, vector<tuple<double, double>>* model){
    double mat, ext, del, ins, tmp;
    for(int t=0; t<T; t++){
        for(int n=0; n<N; n++){
            mat=-INFINITY;
            ext=-INFINITY;
            del=-INFINITY;
            ins=-INFINITY;
            if(t==0 && n==0) {
                ext = 0;
            }

            if (t>0 && n>0) {
                if (t-C>0){
                    tmp=E[(t-C-1)*N+(n-1)] + scoreKmer(sig[t-1], kmer_seq[n-1], model) + m1;
                    for(int l=1; l<=C; l++){
                        tmp+=scoreKmer(sig[t-l-1], kmer_seq[n-1], model);
                    }
                    mat=logPlus(mat, tmp);
                }

                mat=logPlus(mat, D[(t-1)*N+(n-1)] + scoreKmer(sig[t-1], kmer_seq[n-1], model) + m2);
                mat=logPlus(mat, I[(t-1)*N+(n-1)] + scoreKmer(sig[t-1], kmer_seq[n-1], model) + m3);
                ext=logPlus(ext, M[(t-1)*N+n] + scoreKmer(sig[t-1], kmer_seq[n-1], model) + e1); // e1 first extend
                ext=logPlus(ext, E[(t-1)*N+n] + scoreKmer(sig[t-1], kmer_seq[n-1], model) + e2); // e2 extend further
                ext=logPlus(ext, E[(t-1)*N+n] + error(sig[t-1]) + e3); // e3 error

                if (n<N-1) {
                    del=logPlus(del, E[(t-1)*N+n] + indel(sig[t-1], kmer_seq[n-1], kmer_seq[n], model) + d1);
                    del=logPlus(del, D[(t-1)*N+n] + indel(sig[t-1], kmer_seq[n-1], kmer_seq[n], model) + d2);
                    ins=logPlus(ins, E[t*N+(n-1)] + indel(sig[t-1], kmer_seq[n-1], kmer_seq[n], model) + i1);
                    ins=logPlus(ins, I[t*N+(n-1)] + indel(sig[t-1], kmer_seq[n-1], kmer_seq[n], model) + i2);
                }
            }

            M[t*N+n]=mat;
            E[t*N+n]=ext;
            D[t*N+n]=del;
            I[t*N+n]=ins;
        }
    }
}

/**
 * Calculate backward matrices using logarithmic values
 *
 * @param sig ONT raw signal with pA values
 * @param seq nucleotide sequence represented by the ONT signal
 * @param T length of the ONT raw signal + 1
 * @param N length of nucleotide sequence + 1
 * @param model map containing kmers as keys and (mean, stdev) tuples as values
 */
void logB(double* sig, int* kmer_seq, double* M, double* E, double* D, double* I, const int &T, const int &N, vector<tuple<double, double>>* model) {
    double mat, ext, del, ins, tmp;
    for(int t=T-1; t>=0; t--){
        for(int n=N-1; n>=0; n--){
            mat=-INFINITY;
            ext=-INFINITY;
            del=-INFINITY;
            ins=-INFINITY;
            if(t==T-1 && n==N-1) {
                ext = 0;
            }

            // m with minimum length C
            if (t+C+1<T && n+1<N) {
                tmp=M[(t+1+C)*N+(n+1)] + scoreKmer(sig[t], kmer_seq[n], model) + m1;
                for (int l=1; l<=C; l++){
                    tmp+=scoreKmer(sig[t+l], kmer_seq[n], model);
                }
                ext=logPlus(ext, tmp);
            }

            if (t+1<T && n+1<N) {
                ins=logPlus(ins, M[(t+1)*N+(n+1)] + scoreKmer(sig[t], kmer_seq[n], model) + m2);
                del=logPlus(del, M[(t+1)*N+(n+1)] + scoreKmer(sig[t], kmer_seq[n], model) + m3);
            }

            if (t+1<T && n>0) {
                mat=logPlus(mat, E[(t+1)*N+n] + scoreKmer(sig[t], kmer_seq[n-1], model) + e1); // e1 first extend
                ext=logPlus(ext, E[(t+1)*N+n] + scoreKmer(sig[t], kmer_seq[n-1], model) + e2); // e2 extend further
                ext=logPlus(ext, E[(t+1)*N+n] + error(sig[t]) + e3); // e3 error
            }

            if (t+1<T && n>0 && n+1<N) {
                ext=logPlus(ext, D[(t+1)*N+n] + indel(sig[t], kmer_seq[n-1], kmer_seq[n], model) + d1);
                del=logPlus(del, D[(t+1)*N+n] + indel(sig[t], kmer_seq[n-1], kmer_seq[n], model) + d2);
            }

            if (n+1<N-1 && t>0) {
                ext=logPlus(ext, I[t*N+(n+1)] + indel(sig[t-1], kmer_seq[n], kmer_seq[n+1], model) + i1);
                ins=logPlus(ins, I[t*N+(n+1)] + indel(sig[t-1], kmer_seq[n], kmer_seq[n+1], model) + i2);
            }
            
            M[t*N+n] = mat;
            E[t*N+n] = ext;
            D[t*N+n] = del;
            I[t*N+n] = ins;
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
    double* LP = new double[T*N];
    fill_n(LP, T*N, -INFINITY);
    for(int t=0; t<T; t++){
        for(int n=0; n<N; n++){
            int x = t*N+n;
            LP[x] = FOR[x] + BACK[x] - Z;
        }
    }
    return LP;
}

/**
 * Calculate the maximum a posteriori path (MAP) through LP
 *
 */
list<string> getBorders(double* LPM, double* LPE, double* LPD, double* LPI, const int &T, const int &N){
    double* M = new double[T*N];
    double* E = new double[T*N];
    double* D = new double[T*N];
    double* I = new double[T*N];
    fill_n(M, T*N, -INFINITY);
    fill_n(E, T*N, -INFINITY);
    fill_n(D, T*N, -INFINITY);
    fill_n(I, T*N, -INFINITY);
    double mat, ext, del, ins;
    for(int t=0; t<T; t++){
        for(int n=0; n<N; n++){
            mat=-INFINITY;
            ext=-INFINITY;
            del=-INFINITY;
            ins=-INFINITY;
            if(t==0 && n==0) {
                ext = 0;
            }
            if (t-C>0 && n>0){
                mat=max(mat, E[(t-C-1)*N+(n-1)] + LPM[t*N+n]); // m1
            }

            if (t>0 && n>0) {
                ext=max(ext, M[(t-1)*N+n] + LPE[t*N+n]); // e1
                ext=max(ext, E[(t-1)*N+n] + LPE[t*N+n]); // e2, e3
                mat=max(mat, D[(t-1)*N+(n-1)] + LPM[t*N+n]); // m2
                mat=max(mat, I[(t-1)*N+(n-1)] + LPM[t*N+n]); // m3

                if (n<N-1) {
                    del=max(del, E[(t-1)*N+n] + LPD[t*N+n]); // d1
                    del=max(del, D[(t-1)*N+n] + LPD[t*N+n]); // d2
                    ins=max(ins, E[t*N+(n-1)] + LPI[t*N+n]); // i1
                    ins=max(ins, I[t*N+(n-1)] + LPI[t*N+n]); // i2
                }
            }

            M[t*N+n]=mat;
            E[t*N+n]=ext;
            D[t*N+n]=del;
            I[t*N+n]=ins;
        }
    }
    list<string> segString;
    funcE(T-1, N-1, M, E, D, I, LPM, LPE, LPD, LPI, &segString, N);
    delete[] M;
    delete[] E;
    delete[] D;
    delete[] I;
    return segString;
}

void funcM(int t, int n, double* M, double* E, double* D, double* I, double* LPM, double* LPE, double* LPD, double* LPI, list<string>* segString, const int &N){
    double score = M[t*N+n];
    // start value
    if (t<=1 && n<=1){
        segString->push_front("M"+to_string(0)+","+to_string(0)); // n-1 because N is 1 larger than the sequences
        return;
    }
    if (score == E[(t-C-1)*N+(n-1)] + LPM[t*N+n]){
        segString->push_front("M"+to_string(n-1)+","+to_string(t-C-1));
        return funcE(t-C-1, n-1, M, E, D, I, LPM, LPE, LPD, LPI, segString, N);
    }
    if (score == D[(t-1)*N+(n-1)] + LPM[t*N+n]) {
        segString->push_front("M"+to_string(n-1)+","+to_string(t-1));
        return funcD(t-1, n-1, M, E, D, I, LPM, LPE, LPD, LPI, segString, N);
    }
    if (score == I[(t-1)*N+(n-1)] + LPM[t*N+n]) {
        segString->push_front("M"+to_string(n-1)+","+to_string(t-1));
        return funcI(t-1, n-1, M, E, D, I, LPM, LPE, LPD, LPI, segString, N);
    }
}

void funcE(int t, int n, double* M, double* E, double* D, double* I, double* LPM, double* LPE, double* LPD, double* LPI, list<string>* segString, const int &N){
    double score = E[t*N+n];
    if (score == M[(t-1)*N+n] + LPE[t*N+n]){
        return funcM(t-1, n, M, E, D, I, LPM, LPE, LPD, LPI, segString, N);
    }
    if (score == E[(t-1)*N+n] + LPE[t*N+n]){
        return funcE(t-1, n, M, E, D, I, LPM, LPE, LPD, LPI, segString, N);
    }
}

void funcI(int t, int n, double* M, double* E, double* D, double* I, double* LPM, double* LPE, double* LPD, double* LPI, list<string>* segString, const int &N){
    double score = I[t*N+n];
    if (score == I[t*N+(n-1)] + LPI[t*N+n]) {
        return funcI(t, n-1, M, E, D, I, LPM, LPE, LPD, LPI, segString, N);
    }
    if (score == E[t*N+(n-1)] + LPI[t*N+n]) {
        segString->push_front("I"+to_string(n-1)+","+to_string(t));
        return funcE(t, n-1, M, E, D, I, LPM, LPE, LPD, LPI, segString, N);
    }
}

void funcD(int t, int n, double* M, double* E, double* D, double* I, double* LPM, double* LPE, double* LPD, double* LPI, list<string>* segString, const int &N){
    double score = D[t*N+n];
    if (score == D[(t-1)*N+n] + LPD[t*N+n]) {
        return funcI(t-1, n, M, E, D, I, LPM, LPE, LPD, LPI, segString, N);
    }
    if (score == E[(t-1)*N+n] + LPD[t*N+n]) {
        segString->push_front("D"+to_string(n)+","+to_string(t-1));
        return funcE(t-1, n, M, E, D, I, LPM, LPE, LPD, LPI, segString, N);
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
 * Train transition parameter with baum welch algorithm
*/
tuple<double, double, double, double, double, double, double, double, double, double> trainTransition(double* sig, int* kmer_seq, double* forM, double* forE, double* forD, double* forI, double* backM, double* backE, double* backD, double* backI, const int &T, const int &N, vector<tuple<double, double>>* model) {
    // Transition parameters
    double newM1 = -INFINITY;
    double newM2 = -INFINITY;
    double newM3 = -INFINITY;
    double newE1 = -INFINITY;
    double newE2 = -INFINITY;
    double newE3 = -INFINITY;
    double newI1 = -INFINITY;
    double newI2 = -INFINITY;
    double newD1 = -INFINITY;
    double newD2 = -INFINITY;
    double tempM = -INFINITY;

    for(int t=0; t<T; t++){
        for(int n=0; n<N; n++){
            if (t+C+1<T && n+1<N) {
                // m1:  forward(i)    a    e(i+1)                                  backward(i+1)
                tempM = forE[t*N+n] + m1 + scoreKmer(sig[t], kmer_seq[n], model) + backM[(t+C+1)*N+(n+1)];
                for(int l=1; l<=C; l++){
                    tempM+=scoreKmer(sig[t+l], kmer_seq[n], model);
                }
                newM1 = logPlus(newM1, tempM);
            }

            if (t+1<T && n+1<N) {
                newM2 = logPlus(newM2, forD[t*N+n] + m2 + scoreKmer(sig[t], kmer_seq[n], model) + backM[(t+1)*N+(n+1)]);
                newM3 = logPlus(newM3, forI[t*N+n] + m3 + scoreKmer(sig[t], kmer_seq[n], model) + backM[(t+1)*N+(n+1)]);
            }

            if (t+1<T && n>0) {
                newE1 = logPlus(newE1, forM[t*N+n] + e1 + scoreKmer(sig[t], kmer_seq[n-1], model) + backE[(t+1)*N+n]);
                newE2 = logPlus(newE2, forE[t*N+n] + e2 + scoreKmer(sig[t], kmer_seq[n-1], model) + backE[(t+1)*N+n]);
                newE3 = logPlus(newE3, forE[t*N+n] + e3 + error(sig[t])                           + backE[(t+1)*N+n]);
            }

            if (t+1<T && n>0 && n+1<N) {
                newD1 = logPlus(newD1, forE[t*N+n] + d1 + indel(sig[t], kmer_seq[n-1], kmer_seq[n], model) + backD[(t+1)*N+n]);
                newD2 = logPlus(newD2, forD[t*N+n] + d2 + indel(sig[t], kmer_seq[n-1], kmer_seq[n], model) + backD[(t+1)*N+n]);
            }

            if (n+1<N-1 && t>0) {
                newI1 = logPlus(newI1, forE[t*N+n] + i1 + indel(sig[t-1], kmer_seq[n], kmer_seq[n+1], model) + backI[t*N+(n+1)]);
                newI2 = logPlus(newI2, forI[t*N+n] + i2 + indel(sig[t-1], kmer_seq[n], kmer_seq[n+1], model) + backI[t*N+(n+1)]);
            }
        }
    }
    // average over the number of transitions
    double Am = newE1;
    newE1 -= Am;
    double Ae = logPlus(newI1, logPlus(newD1, logPlus(newE2, logPlus(newE3, newM1))));
    newM1 -= Ae;
    newE2 -= Ae;
    newE3 -= Ae;
    newI1 -= Ae;
    newD1 -= Ae;
    double Ad = logPlus(newD2, newM2);
    newD2 -= Ad;
    newM2 -= Ad;
    double Ai = logPlus(newI2, newM3);
    newI2 -= Ai;
    newM3 -= Ai;

    return tuple<double, double, double, double, double, double, double, double, double, double>({exp(newM1), exp(newE1), exp(newE2), exp(newE3), exp(newI1), exp(newD1), exp(newD2), exp(newM2), exp(newI2), exp(newM3)});
}

/**
 * Train emission parameter with baum welch algorithm
*/
tuple<double*, double*> trainEmission(double* sig, int* kmer_seq, double* forM, double* forE, double* backM, double* backE, const int &T, const int &N, vector<tuple<double, double>>* model) {
    // Emission
    // https://courses.grainger.illinois.edu/ece417/fa2021/lectures/lec15.pdf
    // https://f.hubspotusercontent40.net/hubfs/8111846/Unicon_October2020/pdf/bilmes-em-algorithm.pdf
    // gamma_t(i) is the probability of being in state i at time t
    // gamma for state M - expected number of transitions of M at given time (T) for all latent states (kmers)
    double* G = new double[T*N];
    fill_n(G, T*N, -INFINITY);

    for(int t=0; t<T; t++){
        // calibrate with the sum of transitions
        double s = -INFINITY;
        for(int n=0; n<N; n++){
            G[t*N+n] = logPlus(forM[t*N+n] + backM[t*N+n], forE[t*N+n] + backE[t*N+n]);
            s = logPlus(s, G[t*N+n]);
        }
        for(int n=0; n<N; n++){
            if (!isinf(s)) {
                G[t*N+n] -= s;
            }
        }
    }
    // decimal space
    double* kmers = new double[N];
    fill_n(kmers, N, 0);
    double* d = new double[N];
    fill_n(d, N, 0);
    // normal space
    double* means = new double[numKmers];
    fill_n(means, numKmers, 0.0);
    int* counts = new int[(int) numKmers];
    fill_n(counts, numKmers, 0);
    for (int n=0; n<N; n++) {
        if (n>0) {
            counts[kmer_seq[n-1]]++;
        }
        for (int t=0; t<T; t++) {
            if (n>0 && t>0) {
                kmers[n] += exp(G[t*N+n]) * sig[t-1];
            }
            d[n] += exp(G[t*N+n]);
        }
        kmers[n] = kmers[n] / d[n];
    }

    for (int n=1; n<N; n++) {
        means[kmer_seq[n-1]] += kmers[n] / counts[kmer_seq[n-1]];
    }

    // Emission (stdev of kmers)
    fill_n(kmers, N, 0);
    double* stdevs = new double[numKmers];
    fill_n(stdevs, numKmers, 0.0);
    for (int n=0; n<N; n++) {
        for (int t=0; t<T; t++) {
            if (n>0 && t>0) {
                kmers[n] += exp(G[t*N+n]) * pow(sig[t-1] - means[kmer_seq[n-1]], 2.);
            }
        }
        kmers[n] = kmers[n] / d[n];
    }

    for (int n=1; n<N; n++) {
        // transform vars to stdevs
        stdevs[kmer_seq[n-1]] += kmers[n] / counts[kmer_seq[n-1]];
        // stdevs[kmer_seq[n-1]] += kmers[n] / max(counts[kmer_seq[n-1]] - 1, 1);
        stdevs[kmer_seq[n-1]] = sqrt(stdevs[kmer_seq[n-1]]);
    }

    delete[] G;
    delete[] kmers;
    delete[] counts;
    delete[] d;
    return tuple<double*, double*>({means, stdevs});
}

void printTrainedTransitionParams(double* sig, int* kmer_seq, double* forM, double* forE, double* forD, double* forI, double* backM, double* backE, double* backD, double* backI, const int &T, const int &N, vector<tuple<double, double>>* model) {

    // exp(newM1), exp(newE1), exp(newE2), exp(newE3), exp(newI1), exp(newD1), exp(newD2), exp(newM2), exp(newI2), exp(newM3)
    auto [newM1, newE1, newE2, newE3, newI1, newD1, newD2, newM2, newI2, newM3] = trainTransition(sig, kmer_seq, forM, forE, forD, forI, backM, backE, backD, backI, T, N, model);

    // pseudocount
    if (newE3 == 0.0) {
        newE3 =  0.00001;
        newE2 -= 0.0000025;
        newM1 -= 0.0000025;
        newD1 -= 0.0000025;
        newI1 -= 0.0000025;
    }

    cout<<"m1:"<<newM1<<";e1:"<<newE1<<";e2:"<<newE2<<";e3:"<<newE3<<";i1:"<<newI1<<";d1:"<<newD1<<";d2:"<<newD2<<";m2:"<<newM2<<";i2:"<<newI2<<";m3:"<<newM3<<endl;

    auto [newMeans, newStdevs] = trainEmission(sig, kmer_seq, forM, forE, backM, backE, T, N, model);
    for (int i=0; i<numKmers; i++){
        if (newMeans[i]!=0.0){
            cout<<itoa(i)<<":"<<newMeans[i]<<","<<newStdevs[i]<<";";
        }
    }
    cout<<endl;
    cout<<"Z:"<<forE[T*N-1]<<endl;
    cout.flush();
}

/**
 * Read signal and read from stdin until the TERM_STRING is seen
*/
int main(int argc, char* argv[]) {
    // Argparser
    argparse::ArgumentParser program("dynamont basic", "0.1");
    // parameters for DP
    program.add_argument("-e1", "--extendscore1").help("Transition probability for extend rule 1").default_value(1.00).scan<'g', double>(); // e1
    program.add_argument("-m1", "--matchscore1").help("Transition probability for match rule 2").default_value(.035).scan<'g', double>(); // m
    program.add_argument("-e2", "--extendscore2").help("Transition probability for extend rule 2").default_value(.964).scan<'g', double>(); // e2
    program.add_argument("-e3", "--extendscore3").help("Transition probability for extend rule 3").default_value(.001).scan<'g', double>(); // e3
    // program.add_argument("-at", "--atrain").help("Switch algorithm to transition parameter training mode").default_value(false).implicit_value(true);
    program.add_argument("-t", "--train").help("Switch algorithm to transition and emission parameter training mode").default_value(false).implicit_value(true);
    program.add_argument("-z", "--calcZ").help("Switch algorithm to only calculate Z").default_value(false).implicit_value(true);
    program.add_argument("-m", "--model").help("Path to kmer model table").default_value(MODELPATH);
    program.add_argument("-r", "--pore").help("Pore generation used to sequence the data").default_value(9).choices(9, 10);
    program.add_argument("-c", "--minSegLen").help("MinSegLen + 1 is the minimal segment length").default_value(0).store_into(C);

    try {
        program.parse_args(argc, argv);
    }
    catch (const std::runtime_error& err) {
        cerr << err.what() << std::endl;
        cerr << program;
        return 1;
    }
    
    int pore = program.get<int>("pore");
    modelpath = program.get<string>("model");
    // atrain = program.get<bool>("atrain");
    train = program.get<bool>("train");
    calcZ = program.get<bool>("calcZ");
    m1 = log(program.get<double>("matchscore1"));
    e1 = log(program.get<double>("extendscore1"));
    e2 = log(program.get<double>("extendscore2"));
    e3 = log(program.get<double>("extendscore3"));
    
    if (pore == 9) {
        K = 5;
    } else if (pore == 10) {
        K = 9;
    }
    fillBASE2ID();
    numKmers = pow(ALPHABET_SIZE, K);
    vector<tuple<double, double>> model(numKmers, make_tuple(-INFINITY, -INFINITY));
    readKmerModel(modelpath, &model);
    string signal;
    string read;
    int truish = 1;

    while(truish) {
        // truish = 0;
        // echo 107,107,107.2,108.0,108.9,111.2,105.7,104.3,107.1,105.7,105,105 CAAAAA| src\segment.exe
        // read input, signal and read whitespace separated in single line
        getline(cin, signal);
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
        double* sig = new double[T-1];
        fill_n(sig, T-1, -INFINITY);
        string value;
        stringstream ss(signal);
        int i = 0;
        while(getline(ss, value, ',')) {
            sig[i++] = stod(value);
        }
        // process read: convert string to int array
        int N = read.size() + 1; // operate on base transitions
        int seq_size = read.size() + (K-1); 
        int* seq = new int[seq_size];
        fill_n(seq, seq_size, 0); // default: fill with A add 2 As to 3' of read
        i = floor(K/2);
        for (const char &c: read) {
            seq[i] = BASE2ID.at(c);
            i++;
        }
        // add NN to end of sequence
        seq[i] = 4;
        seq[i+1] = 4;

        int* kmer_seq = seq2kmer(seq, N-1);

        // initialize for matrices
        double* forM = new double[T*N];
        fill_n(forM, T*N, -INFINITY);
        double* forE = new double[T*N];
        fill_n(forE, T*N, -INFINITY);
        double* forD = new double[T*N];
        fill_n(forD, T*N, -INFINITY);
        double* forI = new double[T*N];
        fill_n(forI, T*N, -INFINITY);
        // calculate segmentation probabilities, fill forward matrices
        logF(sig, kmer_seq, forM, forE, forD, forI, T, N, &model);

        // initialize back matrices
        double* backM = new double[T*N];
        fill_n(backM, T*N, -INFINITY);
        double* backE = new double[T*N];
        fill_n(backE, T*N, -INFINITY);
        double* backD = new double[T*N];
        fill_n(backD, T*N, -INFINITY);
        double* backI = new double[T*N];
        fill_n(backI, T*N, -INFINITY);
        // calculate segmentation probabilities, fill backward matrices
        logB(sig, kmer_seq, backM, backE, backD, backI, T, N, &model);

        // Numeric error is scaled by input size, Z in forward and backward should match by some numeric error EPSILON
        if ((abs(forE[T*N-1] - backE[0])/(T*N))>EPSILON || isinf(forE[T*N-1]) || isinf(backE[0]) || isnan(forE[T*N-1]) || isnan(backE[0])) {
            cerr << fixed << showpoint;
            cerr << setprecision(20);
            cerr<<"Z values between matrices do not match! forE[T*N-1]: "<<forE[T*N-1]<<", backE[0]: "<<backE[0]<<", "<<endl;
            cerr<<abs(forE[T*N-1] - backE[0])/(T*N)<<" > "<<EPSILON<<endl;
            cerr.flush();
            exit(11);
        }

        if (calcZ){
            cout<<forE[T*N-1]<<"\n";
        
        } else {
            
            // train both Transitions and Emissions
            if (train) {
                printTrainedTransitionParams(sig, kmer_seq, forM, forE, forD, forI, backM, backE, backD, backI, T, N, &model);
            }

            double* LPM = logP(forM, backM, forE[T*N-1], T, N); // log probs
            double* LPE = logP(forE, backE, forE[T*N-1], T, N); // log probs
            double* LPD = logP(forD, backD, forE[T*N-1], T, N); // log probs
            double* LPI = logP(forI, backI, forE[T*N-1], T, N); // log probs
            list<string> segString = getBorders(LPM, LPE, LPD, LPI, T, N);

            for (auto const& seg : segString) {
                cout<<seg;
            }
            cout<<endl;
            cout.flush();

            // // calculate sum of segment probabilities
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
            delete[] LPM;
            delete[] LPE;
            delete[] LPD;
            delete[] LPI;
        }
        delete[] forM;
        delete[] forE;
        delete[] forD;
        delete[] forI;
        delete[] backM;
        delete[] backE;
        delete[] backD;
        delete[] backI;
        delete[] sig;
        delete[] seq;
        delete[] kmer_seq;
    }
    return 0;
}