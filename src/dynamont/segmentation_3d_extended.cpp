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

// void funcA(int t, int n, int k, double* A, double* P, double* S, double* E, double* I, double* LPA, double* LPP, double* LPS, double* LPE, double* LPI, list<string>* segString, const int &N, const int &T);
// void funcP(int t, int n, int k, double* A, double* P, double* S, double* E, double* I, double* LPA, double* LPP, double* LPS, double* LPE, double* LPI, list<string>* segString, const int &N, const int &T);
// void funcS(int t, int n, int k, double* A, double* P, double* S, double* E, double* I, double* LPA, double* LPP, double* LPS, double* LPE, double* LPI, list<string>* segString, const int &N, const int &T);
// void funcE(int t, int n, int k, double* A, double* P, double* S, double* E, double* I, double* LPA, double* LPP, double* LPS, double* LPE, double* LPI, list<string>* segString, const int &N, const int &T);
// void funcI(int t, int n, int k, double* A, double* P, double* S, double* E, double* I, double* LPA, double* LPP, double* LPS, double* LPE, double* LPI, list<string>* segString, const int &N, const int &T);

map<char, int> BASE2ID;
map<char, int> ID2BASE;
string modelpath;
int ALPHABET_SIZE, kmerSize, C, numKmers; // size of acceptable base characters, kmer size, min segment length, number of kmers
double EPSILON = pow(10, -2);
bool train, calcZ; // atrain
double a1, a2, p1, p2, p3, s1, s2, s3, e1, e2, e3, e4, i1, i2; // transition parameters

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
    ID2BASE.insert(pair<char, char>('0', 'A'));
    ID2BASE.insert(pair<char, char>('1', 'C'));
    ID2BASE.insert(pair<char, char>('2', 'G'));
    ID2BASE.insert(pair<char, char>('3', 'T'));
    ID2BASE.insert(pair<char, char>('4', 'N'));
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
    int base = kmerSize;

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
 * Works ONLY if ALPHABET_SIZE is smaller or equal to 10!
 * 
 * @param i input number in the given base as an array
*/
int toDeci(int *i) {
    int ret = 0;
    int m = 1;
    for(int r = kmerSize - 1; r >= 0; r--) {
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
    int ret = 0;
    for (char const &c:s){
        // assert (BASE2ID.at(c)>=0); // check if nucleotide is known
        ret*=kmerSize; // move the number in base to the left
        ret+=BASE2ID.at(c);
    }
    return ret;
}

/**
 * Calculates the integer representation of the successing kmer given the current kmer and the upcoming nucleotide
 * k_i+1 = (k_i mod base^(kmerSize-1)) * base + value(nextNt, base)
 * 
 * @param currentKmer
 * @param nextNt
 * @return successing Kmer as integer representation in the current base
 */
int successingKmer(int currentKmer, int nextNt) {
    return (currentKmer % int(pow(ALPHABET_SIZE, kmerSize-1))) * ALPHABET_SIZE + nextNt;
}

/**
 * Calculates the integer representation of the precessor kmer given the current kmer and the precessing nucleotide
 * k_i-1 = int(k_i/base) + value(priorNt, base) * base^(kmerSize-1)
 * 
 * @param currentKmer
 * @param priorNt
 * @return precessing Kmer as integer representation in the current base
 */
int precessingKmer(int currentKmer, int priorNt) {
    return currentKmer/ALPHABET_SIZE + priorNt * int(pow(ALPHABET_SIZE, kmerSize-1));
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
    int* tempKmer = new int[kmerSize];
    for (int n=0; n<N; n++){ // extend loop to ad 2 Ns at start and end of read
        copy(seq + n, seq + n+kmerSize, tempKmer);
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
 * Calculates the Hamming-Distance between two given kmers in their integer base representation
 */
int distanceSequenceKmer(const int &kmer_N, const int &kmer_K) { 
    int acc = 0;
    div_t dv_N{}; dv_N.quot=kmer_N;
    div_t dv_K{}; dv_K.quot=kmer_K;
    for(int i=0; i<kmerSize; i++){
        dv_N = div(dv_N.quot, kmerSize);
        dv_K = div(dv_K.quot, kmerSize);
        acc += (dv_N.rem != dv_K.rem);
    }
    return acc;
}

/**
 * Return log probability density for a given value and a given normal distribution.
 * affineScale = 0.05 -> https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-024-10440-w#Fig1
 * 
 * @param signal point to calculate probability density
 * @param kmer key for the model kmer:(mean, stdev) map
 * @param model map containing kmers as keys and (mean, stdev) tuples as values
 * @return log probability density value for x in the given normal distribution
 */
inline double score(const double &signal_T, const int &kmer_N, const int &kmer_K, const double &affineScale, vector<tuple<double, double>>* model) {
    tuple<double, double> kmerModel_N = (*model)[kmer_N];
    tuple<double, double> kmerModel_K = (*model)[kmer_K];

    double scoreNT = 2*(log_normal_pdf(signal_T, get<0>(kmerModel_N), get<1>(kmerModel_N)) + 6);
    double scoreKT = 2*(log_normal_pdf(signal_T, get<0>(kmerModel_K), get<1>(kmerModel_K)) + 6);
    double scoreNK = -distanceSequenceKmer(kmer_N, kmer_K) * affineScale; // log(exp(-HD * affineCost)) = -HD * affineCost

    return scoreNT + scoreKT + scoreNK;
}

/**
 * Return log probability density for a given value and a given normal distribution
 *
 * @param signal point to calculate probability density
 * @param kmer key for the model kmer:(mean, stdev) map
 * @param model map containing kmers as keys and (mean, stdev) tuples as values
 * @return log probability density value for x in the given normal distribution
 */
inline double score(const double &signal_T, const int &kmer_N, const int &kmer_K, vector<tuple<double, double>>* model) {
    tuple<double, double> kmerModel_N = (*model)[kmer_N];
    tuple<double, double> kmerModel_K = (*model)[kmer_K];

    double scoreNT = 2*(log_normal_pdf(signal_T, get<0>(kmerModel_N), get<1>(kmerModel_N)) + 6);
    double scoreKT = 2*(log_normal_pdf(signal_T, get<0>(kmerModel_K), get<1>(kmerModel_K)) + 6);
    double scoreNK = -distanceSequenceKmer(kmer_N, kmer_K); // log(exp(-HD)) = -HD

    return scoreNT + scoreKT + scoreNK;
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

/**
 * Calculate forward matrices using logarithmic values
 *
 * @param sig ONT raw signal with pA values
 * @param seq nucleotide sequence represented by the ONT signal
 * @param T length of the ONT raw signal + 1
 * @param N length of nucleotide sequence + 1
 * @param model map containing kmers as keys and (mean, stdev) tuples as values
 */
void logF(double* sig, int* kmer_seq, double* A, double* P, double* S, double* E, double* I, const int &T, const int &N, const int &K, vector<tuple<double, double>>* model){
    double a, p, s, e, i;
    int preKmer;
    for(int t=0; t<T; t++){
        for(int n=0; n<N; n++){
            for(int k=0; k<K; k++){
                a=-INFINITY;
                p=-INFINITY;
                s=-INFINITY;
                e=-INFINITY;
                i=-INFINITY;
                
                // init, log(1) for any k
                if(t==0 && n==0){
                    e=0;
                }

                if(t>0 && n>0){
                    for (int prevK=0; prevK<ALPHABET_SIZE; prevK++) {
                        preKmer = precessingKmer(k, prevK);
                        a=logPlus(a, E[(t-1)*(N*K)+(n-1)*K+preKmer] + a1 + score(sig[t-1], kmer_seq[n-1], k, 0.05, model));
                        a=logPlus(a, I[(t-1)*(N*K)+(n-1)*K+preKmer] + a2 + score(sig[t-1], kmer_seq[n-1], k, 0.05, model));

                        p=logPlus(p, S[(t-1)*(N*K)+n*K+preKmer] + p1 + score(sig[t-1], kmer_seq[n-1], k, 0.05, model));
                        p=logPlus(p, E[(t-1)*(N*K)+n*K+preKmer] + p2 + score(sig[t-1], kmer_seq[n-1], k, 0.05, model));
                        p=logPlus(p, I[(t-1)*(N*K)+n*K+preKmer] + p3 + score(sig[t-1], kmer_seq[n-1], k, 0.05, model));
                    }
                    s=logPlus(s, P[(t-1)*(N*K)+(n-1)*K+k] + s1 + score(sig[t-1], kmer_seq[n-1], k, 0.05, model));
                    s=logPlus(s, E[(t-1)*(N*K)+(n-1)*K+k] + s2 + score(sig[t-1], kmer_seq[n-1], k, 0.05, model));
                    s=logPlus(s, I[(t-1)*(N*K)+(n-1)*K+k] + s3 + score(sig[t-1], kmer_seq[n-1], k, 0.05, model));

                    e=logPlus(e, A[(t-1)*(N*K)+n*K+k] + e1 + score(sig[t-1], kmer_seq[n-1], k, model));
                    e=logPlus(e, P[(t-1)*(N*K)+n*K+k] + e2 + score(sig[t-1], kmer_seq[n-1], k, model));
                    e=logPlus(e, S[(t-1)*(N*K)+n*K+k] + e3 + score(sig[t-1], kmer_seq[n-1], k, model));
                    e=logPlus(e, E[(t-1)*(N*K)+n*K+k] + e4 + score(sig[t-1], kmer_seq[n-1], k, model));

                    i=logPlus(i, E[t*(N*K)+(n-1)*K+k] + i1 + score(sig[t-1], kmer_seq[n-1], k, 0.05, model));
                    i=logPlus(i, I[t*(N*K)+(n-1)*K+k] + i2 + score(sig[t-1], kmer_seq[n-1], k, 0.05, model));
                }
                A[t*(N*K)+n*K+k]=a;
                P[t*(N*K)+n*K+k]=p;
                S[t*(N*K)+n*K+k]=s;
                E[t*(N*K)+n*K+k]=e;
                I[t*(N*K)+n*K+k]=i;
            }
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
void logB(double* sig, int* kmer_seq, double* A, double* P, double* S, double* E, double* I, const int &T, const int &N, const int &K, vector<tuple<double, double>>* model){
    double a, p, s, e, i;
    int sucKmer;
    for(int t=T-1; t>=0; t--){
        for(int n=N-1; n>=0; n--){
            for(int k=K-1; k>=0; k--){
                a=-INFINITY;
                p=-INFINITY;
                s=-INFINITY;
                e=-INFINITY;
                i=-INFINITY;

                // init, log(1) for any k
                if (t==T-1 && n==N-1) {
                    e=0;
                }

                if (t<T-1) {
                    if (n>0) {
                        a=logPlus(a, E[(t+1)*(N*K)+n*K+k] + e1 + score(sig[t], kmer_seq[n-1], k, model));
                        p=logPlus(p, E[(t+1)*(N*K)+n*K+k] + e2 + score(sig[t], kmer_seq[n-1], k, model));
                        s=logPlus(s, E[(t+1)*(N*K)+n*K+k] + e3 + score(sig[t], kmer_seq[n-1], k, model));
                        e=logPlus(e, E[(t+1)*(N*K)+n*K+k] + e4 + score(sig[t], kmer_seq[n-1], k, model));
                        if (k<K-1) {
                            for (int nextK=0; nextK<ALPHABET_SIZE; nextK++){
                                sucKmer = successingKmer(k, nextK);
                                s=logPlus(s, P[(t+1)*(N*K)+n*K+sucKmer] + p1 + score(sig[t], kmer_seq[n-1], sucKmer, 0.05, model));
                                e=logPlus(e, P[(t+1)*(N*K)+n*K+sucKmer] + p2 + score(sig[t], kmer_seq[n-1], sucKmer, 0.05, model));
                                i=logPlus(i, P[(t+1)*(N*K)+n*K+sucKmer] + p3 + score(sig[t], kmer_seq[n-1], sucKmer, 0.05, model));
                            }
                        }
                    }

                    if (n<N-1) {
                        p=logPlus(p, S[(t+1)*(N*K)+(n+1)*K+k] + s1 + score(sig[t], kmer_seq[n], k, 0.05, model));
                        e=logPlus(e, S[(t+1)*(N*K)+(n+1)*K+k] + s2 + score(sig[t], kmer_seq[n], k, 0.05, model));
                        i=logPlus(i, S[(t+1)*(N*K)+(n+1)*K+k] + s3 + score(sig[t], kmer_seq[n], k, 0.05, model));
                        
                        for (int nextK=0; nextK<ALPHABET_SIZE; nextK++){
                            sucKmer = successingKmer(k, nextK);
                            e=logPlus(e, A[(t+1)*(N*K)+(n+1)*K+sucKmer] + a1 + score(sig[t], kmer_seq[n], sucKmer, 0.05, model));
                            i=logPlus(i, A[(t+1)*(N*K)+(n+1)*K+sucKmer] + a2 + score(sig[t], kmer_seq[n], sucKmer, 0.05, model));
                        }
                        
                    }
                }

                if (t>0 && n<N-1) {
                    e=logPlus(e, I[t*(N*K)+(n+1)*K+k] + i1 + score(sig[t-1], kmer_seq[n], k, 0.05, model));
                    i=logPlus(i, I[t*(N*K)+(n+1)*K+k] + i2 + score(sig[t-1], kmer_seq[n], k, 0.05, model));                    
                }

                A[t*(N*K)+n*K+k]=a;
                P[t*(N*K)+n*K+k]=p;
                S[t*(N*K)+n*K+k]=s;
                E[t*(N*K)+n*K+k]=e;
                I[t*(N*K)+n*K+k]=i;
            }
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
double* logP(double* FOR, double* BACK, const double &Z, const int &T, const int &N, const int &K) {
    double* LP = new double[T*N*K];
    fill_n(LP, T*N*K, -INFINITY);
    for(int t=0; t<T; t++){
        for(int n=0; n<N; n++){
            for(int k=0; k<K; k++){
                int x = t*(N*K)+n*K+k;
                LP[x] = FOR[x] + BACK[x] - Z;
            }
        }
    }
    return LP;
}

/**
 * Calculate the maximum a posteriori path through LP
 *
 */
// list<string> getBorders(double* LPA, double* LPP, double* LPS, double* LPE, double* LPI, const int &T, const int &N, const int &K){
//     double* A = new double[T*N*K];
//     double* P = new double[T*N*K];
//     double* S = new double[T*N*K];
//     double* E = new double[T*N*K];
//     double* I = new double[T*N*K];
//     fill_n(A, T*N*K, -INFINITY);
//     fill_n(P, T*N*K, -INFINITY);
//     fill_n(S, T*N*K, -INFINITY);
//     fill_n(E, T*N*K, -INFINITY);
//     fill_n(I, T*N*K, -INFINITY);
//     double a, p, s, e, i;
//     int idx;
//     for(int t=0; t<T; t++){
//         for(int n=0; n<N; n++){
//             for(int k=0; k<K; k++){
//                 a=-INFINITY;
//                 p=-INFINITY;
//                 s=-INFINITY;
//                 e=-INFINITY;
//                 i=-INFINITY;
                
//                 idx = t*(N*K)+n*K+k;

//                 if(t==0 && n==0 && k==0){
//                     e=0;
//                 }

//                 if(t>0 && n>0 && k>0){
//                     for (int prevK=0; prevK<ALPHABET_SIZE; prevK++) {
//                         a=max(a, E[(t-1)*(N*K)+(n-1)*K+precessingKmer(k, prevK)] + LPA[idx]);
//                         a=max(a, I[(t-1)*(N*K)+(n-1)*K+precessingKmer(k, prevK)] + LPA[idx]);

//                         p=max(p, S[(t-1)*(N*K)+n*K+precessingKmer(k, prevK)] + LPP[idx]);
//                         p=max(p, E[(t-1)*(N*K)+n*K+precessingKmer(k, prevK)] + LPP[idx]);
//                         p=max(p, I[(t-1)*(N*K)+n*K+precessingKmer(k, prevK)] + LPP[idx]);
//                     }
//                 }

//                 if (t>0 && n>0) {
//                     s=max(s, P[(t-1)*(N*K)+(n-1)*K+k] + LPS[idx]);
//                     s=max(s, E[(t-1)*(N*K)+(n-1)*K+k] + LPS[idx]);
//                     s=max(s, I[(t-1)*(N*K)+(n-1)*K+k] + LPS[idx]);

//                     e=max(e, A[(t-1)*(N*K)+n*K+k] + LPE[idx]);
//                     e=max(e, P[(t-1)*(N*K)+n*K+k] + LPE[idx]);
//                     e=max(e, S[(t-1)*(N*K)+n*K+k] + LPE[idx]);
//                     e=max(e, E[(t-1)*(N*K)+n*K+k] + LPE[idx]);

//                     i=max(i, E[t*(N*K)+(n-1)*K+k] + LPI[idx]);
//                     i=max(i, I[t*(N*K)+(n-1)*K+k] + LPI[idx]);
//                 }
//                 A[idx]=a;
//                 P[idx]=p;
//                 S[idx]=s;
//                 E[idx]=e;
//                 I[idx]=i;
//             }
//         }
//     }
//     list<string> segString;
//     funcE(T-1, N-1, K-1, A, P, S, E, I, LPA, LPP, LPS, LPE, LPI, &segString, N, T);
//     delete[] A;
//     delete[] P;
//     delete[] S;
//     delete[] E;
//     delete[] I;
//     return segString;
// }

// void funcA(int t, int n, int k, double* A, double* P, double* S, double* E, double* I, double* LPA, double* LPP, double* LPS, double* LPE, double* LPI, list<string>* segString, const int &N, const int &T){
//     double score = A[t*(N*K)+n*K+k];
//     if (t<=1 && n<=1){ // Start value
//         segString->push_front("M"+to_string(0)+","+to_string(0)); // n-1 because N is 1 larger than the sequences
//         return;
//     }
//     // TODO fix to get correct k here
//     if (t>0 && n>0 && score == E[(t-1)*N+(n-1)] + LPA[t*(N*K)+n*K+k]){
//         segString->push_front("M"+to_string(n-1)+","+to_string(t-1));
//         return funcE(t-1, n-1, k, A, P, S, E, I, LPA, LPP, LPS, LPE, LPI, segString, N, T);
//     }
//     if (t>0 && n>0 && score == I[(t-1)*N+(n-1)] + LPA[t*(N*K)+n*K+k]){
//         segString->push_front("M"+to_string(n-1)+","+to_string(t-1));
//         return funcI(t-1, n-1, k, A, P, S, E, I, LPA, LPP, LPS, LPE, LPI, segString, N, T);
//     }
// }

// void funcE(int t, int n, int k, double* A, double* P, double* S, double* E, double* I, double* LPA, double* LPP, double* LPS, double* LPE, double* LPI, list<string>* segString, const int &N, const int &T){
    // TODO

    // double score = E[t*N+n];
    // if (t>0 && n>0 && score == M[(t-1)*N+n] + LPE[t*N+n]){
    //     return funcM(t-1, n, M, E, LPM, LPE, segString, N);
    // }
    // if (t>0 && n>0 && score == E[(t-1)*N+n] + LPE[t*N+n]){
    //     return funcE(t-1, n, M, E, LPM, LPE, segString, N);
    // }
// }

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
// tuple<double, double, double, double> trainTransition(double* sig, int* kmer_seq, double* forM, double* forE, double* backM, double* backE, const int &T, const int &N, vector<tuple<double, double>>* model) {
//     // Transition parameters
//     double newM1 = -INFINITY;
//     double newE1 = -INFINITY;
//     double newE2 = -INFINITY;
//     double newE3 = -INFINITY;
//     double tempM = -INFINITY;

//     for(int t=0; t<T; t++){
//         for(int n=0; n<N; n++){
//             if (n+1<N && t+C+1<T) {
//                 // m1:  forward(i)        a    e(i+1)                                  backward(i+1)
//                 tempM = forE[t*N+n] + m1 + scoreKmer(sig[t], kmer_seq[n], model) + backM[(t+C+1)*N+(n+1)];
//                 for(int l=1; l<=C; l++){
//                     tempM+=scoreKmer(sig[t+l], kmer_seq[n], model);
//                 }
//                 newM1 = logPlus(newM1, tempM);
//             }

//             if (t+1<T && n>0) {
//                 newE1 = logPlus(newE1, forM[t*N+n] + e1 + scoreKmer(sig[t], kmer_seq[n-1], model) + backE[(t+1)*N+n]);
//                 newE2 = logPlus(newE2, forE[t*N+n] + e2 + scoreKmer(sig[t], kmer_seq[n-1], model) + backE[(t+1)*N+n]);
//                 newE3 = logPlus(newE3, forE[t*N+n] + e3 + error(sig[t])                           + backE[(t+1)*N+n]);
//             }
//         }
//     }
//     // average over the number of transitions
//     double Am = newE1;
//     newE1 = newE1 - Am;
//     double Ae = logPlus(newE2, logPlus(newE3, newM1));
//     newM1 = newM1 - Ae;
//     newE2 = newE2 - Ae;
//     newE3 = newE3 - Ae;

//     return tuple<double, double, double, double>({exp(newM1), exp(newE1), exp(newE2), exp(newE3)});
// }

/**
 * Train emission parameter with baum welch algorithm
*/
// tuple<double*, double*> trainEmission(double* sig, int* kmer_seq, double* forA, double* forP, double* forS, double* forE, double* forI, double* backA, double* backP, double* backS, double* backE, double* backI, const int &T, const int &N, const int &K, vector<tuple<double, double>>* model) {
//     // Emission
//     // https://courses.grainger.illinois.edu/ece417/fa2021/lectures/lec15.pdf
//     // https://f.hubspotusercontent40.net/hubfs/8111846/Unicon_October2020/pdf/bilmes-em-algorithm.pdf
//     // gamma_t(i) is the probability of being in state i at time t
//     // gamma for state M - expected number of transitions of M at given time (T) for all latent states (kmers)
//     double* G = new double[T*N*K];
//     fill_n(G, T*N*K, -INFINITY);

//     // TODO
//     for(int t=0; t<T; t++){
//         // calibrate with the sum of transitions
//         double s = -INFINITY;
//         for(int n=0; n<N; n++){
//             G[t*N+n] = logPlus(forM[t*N+n] + backM[t*N+n], forE[t*N+n] + backE[t*N+n]);
//             s = logPlus(s, G[t*N+n]);
//         }
//         for(int n=0; n<N; n++){
//             if (!isinf(s)) {
//                 G[t*N+n] -= s;
//             }
//         }
//     }
//     // decimal space
//     double* kmers = new double[N];
//     fill_n(kmers, N, 0);
//     double* d = new double[N];
//     fill_n(d, N, 0);
//     // normal space
//     double* means = new double[numKmers];
//     fill_n(means, numKmers, 0.0);
//     int* counts = new int[(int) numKmers];
//     fill_n(counts, numKmers, 0);
//     for (int n=0; n<N; n++) {
//         if (n>0) {
//             counts[kmer_seq[n-1]]++;
//         }
//         for (int t=0; t<T; t++) {
//             if (n>0 && t>0) {
//                 kmers[n] += exp(G[t*N+n]) * sig[t-1];
//             }
//             d[n] += exp(G[t*N+n]);
//         }
//         kmers[n] = kmers[n] / d[n];
//     }

//     for (int n=1; n<N; n++) {
//         means[kmer_seq[n-1]] += kmers[n] / counts[kmer_seq[n-1]];
//     }

//     // Emission (stdev of kmers)
//     fill_n(kmers, N, 0);
//     double* stdevs = new double[numKmers];
//     fill_n(stdevs, numKmers, 0.0);
//     for (int n=0; n<N; n++) {
//         for (int t=0; t<T; t++) {
//             if (n>0 && t>0) {
//                 kmers[n] += exp(G[t*N+n]) * pow(sig[t-1] - means[kmer_seq[n-1]], 2.);
//             }
//         }
//         kmers[n] = kmers[n] / d[n];
//     }

//     for (int n=1; n<N; n++) {
//         // transform vars to stdevs
//         stdevs[kmer_seq[n-1]] += kmers[n] / counts[kmer_seq[n-1]];
//         // stdevs[kmer_seq[n-1]] += kmers[n] / max(counts[kmer_seq[n-1]] - 1, 1);
//         stdevs[kmer_seq[n-1]] = sqrt(stdevs[kmer_seq[n-1]]);
//     }

//     delete[] G;
//     delete[] kmers;
//     delete[] counts;
//     delete[] d;
//     return tuple<double*, double*>({means, stdevs});
// }

// void printTrainedTransitionParams(double* sig, int* kmer_seq, double* forA, double* forP, double* forS, double* forE, double* forI, double* backA, double* backP, double* backS, double* backE, double* backI, const int &T, const int &N, vector<tuple<double, double>>* model) {

//     auto [a1, a2, p1, p2, p3, s1, s2, s3, e1, e2, e3, e4, i1, i2] = trainTransition(sig, kmer_seq, forA, forP, forS, forE, forI, backA, backP, backS, backE, backI, T, N, model);

//     // pseudocount
//     // if (newE3 == 0.0) {
//     //     newE3 = 0.00001;
//     //     newE2 -= 0.000005;
//     //     newM -= 0.000005;
//     // }
//     cout<<"a1:"<<a1<<";a2:"<<a2<<";p1:"<<p1<<";p2:"<<p2<<";p2:"<<p3<<";s1:"<<s1<<";s2:"<<s2<<";s3:"<<s3<<";e1:"<<e1<<";e2:"<<e2<<";e3:"<<e3<<";e4:"<<e4<<";i1:"<<i1<<";i2:"<<i2<<endl;

//     auto [newMeans, newStdevs] = trainEmission(sig, kmer_seq, forA, forP, forS, forE, forI, backA, backP, backS, backE, backI, T, N, model);
//     for (int i=0; i<numKmers; i++){
//         if (newMeans[i]!=0.0){
//             cout<<itoa(i)<<":"<<newMeans[i]<<","<<newStdevs[i]<<";";
//         }
//     }
//     cout<<endl;
//     cout<<"Z:"<<forE[T*N*K-1]<<endl;
//     cout.flush();
// }

/**
 * Read signal and read from stdin until the TERM_STRING is seen
*/
int main(int argc, char* argv[]) {
    // Argparser
    argparse::ArgumentParser program("dynamont 3d extended", "0.1");
    // parameters for DP

    // double a1, a2, p1, p2, p3, s1, s2, s3, e1, e2, e3, e4, i1, i2; // transition parameters
    program.add_argument("-e1", "--extendscore1").help("Transition parameter").default_value(1.00).scan<'g', double>(); // e1

    program.add_argument("-e2", "--extendscore2").help("Transition parameter").default_value(0.50).scan<'g', double>(); // e2
    program.add_argument("-s1", "--sequencescore1").help("Transition parameter").default_value(0.50).scan<'g', double>(); // s1

    program.add_argument("-e3", "--extendscore3").help("Transition parameter").default_value(0.50).scan<'g', double>(); // e3
    program.add_argument("-p1", "--polishscore1").help("Transition parameter").default_value(0.50).scan<'g', double>(); // p1

    program.add_argument("-a1", "--alignscore1").help("Transition parameter").default_value(0.20).scan<'g', double>(); // a1
    program.add_argument("-p2", "--polishscore2").help("Transition parameter").default_value(0.20).scan<'g', double>(); // p2
    program.add_argument("-e4", "--extendscore4").help("Transition parameter").default_value(0.20).scan<'g', double>(); // e4
    program.add_argument("-s2", "--sequencescore2").help("Transition parameter").default_value(0.20).scan<'g', double>(); // s2
    program.add_argument("-i1", "--insertionscore1").help("Transition parameter").default_value(0.20).scan<'g', double>(); // i1

    program.add_argument("-a2", "--alignscore2").help("Transition parameter").default_value(0.25).scan<'g', double>(); // a2
    program.add_argument("-p3", "--polishscore3").help("Transition parameter").default_value(0.25).scan<'g', double>(); // p3
    program.add_argument("-s3", "--sequencescore3").help("Transition parameter").default_value(0.25).scan<'g', double>(); // s3
    program.add_argument("-i2", "--insertionscore2").help("Transition parameter").default_value(0.25).scan<'g', double>(); // i2

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
    train = program.get<bool>("train");
    calcZ = program.get<bool>("calcZ");

    a1 = log(program.get<double>("alignscore1"));
    a2 = log(program.get<double>("alignscore2"));
    p1 = log(program.get<double>("polishscore1"));
    p2 = log(program.get<double>("polishscore2"));
    p3 = log(program.get<double>("polishscore3"));
    s1 = log(program.get<double>("sequencescore1"));
    s2 = log(program.get<double>("sequencescore2"));
    s3 = log(program.get<double>("sequencescore3"));
    e1 = log(program.get<double>("extendscore1"));
    e2 = log(program.get<double>("extendscore2"));
    e3 = log(program.get<double>("extendscore3"));
    e4 = log(program.get<double>("extendscore4"));
    i1 = log(program.get<double>("insertionscore1"));
    i2 = log(program.get<double>("insertionscore2"));
    
    if (pore == 9) {
        kmerSize = 5;
    } else if (pore == 10) {
        kmerSize = 9;
    }
    fillBASE2ID();
    numKmers = int(pow(ALPHABET_SIZE, kmerSize));
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
            return 0;
        // ... or signal or read is missing
        } else if (signal.empty()) {
            cout<<"Signal missing!\n";
            cout.flush();
            return 1;
        } else if (read.empty()) {
            cout<<"Read missing!\n";
            cout.flush();
            return 2;
        }
        
        // process signal T: convert string to double array
        int T = count(signal.begin(), signal.end(), ',')+2; // len(sig) + 1
        double* sig = new double[T-1];
        fill_n(sig, T-1, -INFINITY);
        string value;
        stringstream ss(signal);
        int i = 0;
        while(getline(ss, value, ',')) {
            sig[i++] = stod(value);
        }
        // process read N: convert string to int array
        int N = read.size() + 1; // operate on base transitions
        int seq_size = read.size() + (kmerSize-1); 
        int* seq = new int[seq_size];
        fill_n(seq, seq_size, 0); // default: fill with A add 2 As to 3' of read
        i = floor(kmerSize/2);
        for (const char &c: read) {
            seq[i] = BASE2ID.at(c);
            i++;
        }
        // add NN to end of sequence
        seq[i] = 4;
        seq[i+1] = 4;

        // polishing dimension K
        int K = pow(4, 5); // currently acceptable A, C, G, T, N

        int* kmer_seq = seq2kmer(seq, N-1);

        // initialize forward matrices
        double* forA = new double[T*N*K];
        fill_n(forA, T*N*K, -INFINITY);
        double* forP = new double[T*N*K];
        fill_n(forP, T*N*K, -INFINITY);
        double* forS = new double[T*N*K];
        fill_n(forS, T*N*K, -INFINITY);
        double* forE = new double[T*N*K];
        fill_n(forE, T*N*K, -INFINITY);
        double* forI = new double[T*N*K];
        fill_n(forI, T*N*K, -INFINITY);
        // calculate segmentation probabilities, fill forward matrices
        logF(sig, kmer_seq, forA, forP, forS, forE, forI, T, N, K, &model);

        // initialize backward matrices
        double* backA = new double[T*N*K];
        fill_n(backA, T*N*K, -INFINITY);
        double* backP = new double[T*N*K];
        fill_n(backP, T*N*K, -INFINITY);
        double* backS = new double[T*N*K];
        fill_n(backS, T*N*K, -INFINITY);
        double* backE = new double[T*N*K];
        fill_n(backE, T*N*K, -INFINITY);
        double* backI = new double[T*N*K];
        fill_n(backI, T*N*K, -INFINITY);
        // calculate segmentation probabilities, fill backward matrices
        logB(sig, kmer_seq, backA, backP, backS, backE, backI, T, N, K, &model);

        double Zf = -INFINITY;
        double Zb = -INFINITY;

        // init, log(1) for any k
        for(int k=0; k<K; k++){
            Zf = logPlus(Zf, forE[(T-1)*(N*K)+(N-1)*K+k]);
            Zb = logPlus(Zb, backE[k]);
        }

        double MaxFe = -INFINITY;
        double MaxFa = -INFINITY;
        double MaxBe = -INFINITY;
        double MaxBa = -INFINITY;
        int idx;
        for(int t=T-1; t>=0; t--){
            for(int n=N-1; n>=0; n--){
                for(int k=K-1; k>=0; k--){
                    idx = t*(N*K)+n*K+k;
                    MaxFe=max(MaxFe, forE[idx]);
                    MaxBe=max(MaxBe, backE[idx]);
                    MaxFa=max(MaxFa, forA[idx]);
                    MaxBa=max(MaxBa, backA[idx]);
                }
            }
        }
        cout<<"MaxFe: "<<MaxFe<<endl;
        cout<<"MaxBe: "<<MaxBe<<endl;
        cout<<"MaxFa: "<<MaxFa<<endl;
        cout<<"MaxBa: "<<MaxBa<<endl;


        // Numeric error is scaled by input size, Z in forward and backward should match by some numeric error EPSILON
        if ((isinf(Zf) || isinf(Zb) || isnan(Zf) || isnan(Zb) || abs(Zf-Zb)/(T*N*K)>EPSILON)) {
            cerr << fixed << showpoint;
            cerr << setprecision(20);
            cerr<<"Z values between matrices do not match! forZ: "<<Zf<<", backZ: "<<Zb<<", "<<abs(Zf-Zb)/(T*N*K)<<" > "<<EPSILON<<endl;
            cerr.flush();
            exit(11);
        }

        cerr<<"forZ: "<<Zf<<", backZ: "<<Zb<<", "<<abs(Zf-Zb)/(T*N*K)<<" > "<<EPSILON<<endl;


        // // Numeric error is scaled by input size, Z in forward and backward should match by some numeric error EPSILON
        // if ((isinf(forE[T*N*K-1]) || isinf(backE[0]) || isnan(forE[T*N*K-1]) || isnan(backE[0] || abs(forE[T*N*K-1] - backE[0])/(T*N*K))>EPSILON)) {
        //     cerr << fixed << showpoint;
        //     cerr << setprecision(20);
        //     cerr<<"Z values between matrices do not match! forE[T*N*K-1]: "<<forE[T*N*K-1]<<", backE[0]: "<<backE[0]<<", "<<abs(forE[T*N*K-1] - backE[0])/(T*N)<<" > "<<EPSILON<<endl;
        //     cerr.flush();
        //     exit(11);
        // }

        // if (calcZ){
        //     cout<<forE[T*N*K-1]<<"\n";
        //     cout.flush();
        // } else {
        //     // train both Transitions and Emissions
        //     if (train) {
        //         printTrainedTransitionParams(sig, kmer_seq, forA, forP, forS, forE, forI, backA, backP, backS, backE, backI, T, N, &model);
        //     } else {
        //         double* LPA = logP(forA, backA, forE[T*N*K-1], T, N); // log probs
        //         double* LPP = logP(forP, backP, forE[T*N*K-1], T, N); // log probs
        //         double* LPS = logP(forS, backS, forE[T*N*K-1], T, N); // log probs
        //         double* LPE = logP(forE, backE, forE[T*N*K-1], T, N); // log probs
        //         double* LPI = logP(forI, backI, forE[T*N*K-1], T, N); // log probs
        //         list<string> segString = getBorders(LPA, LPP, LPS, LPE, LPI, T, N);

        //         for (auto const& seg : segString) {
        //             cout<<seg;
        //         }
        //         cout<<endl;
        //         cout.flush();

        //         // Clean up
        //         delete[] LPA;
        //         delete[] LPP;
        //         delete[] LPS;
        //         delete[] LPE;
        //         delete[] LPI;
        //     }
        // }

        delete[] forA;
        delete[] forP;
        delete[] forS;
        delete[] forE;
        delete[] forI;
        delete[] backA;
        delete[] backP;
        delete[] backS;
        delete[] backE;
        delete[] backI;
        delete[] sig;
        delete[] seq;
        delete[] kmer_seq;
    }
    return 0;
}