// author: Jannes Spangenberg
// e-mail: jannes.spangenberg@uni-jena.de
// github: https://github.com/JannesSP
// website: https://jannessp.github.io

#include <iostream>
#include <iomanip>
#include <fstream> // file io
#include <sstream> // file io
#include <string>
#include <algorithm>
#include <vector>
#include <array>
#include <unordered_map>
#include <tuple>
#include <bits/stdc++.h> // reverse strings
#include <cmath> // exp
// #include <limits> // for inifinity
#include <assert.h>
#include <stdlib.h>
#include "argparse.hpp"

// #include <chrono> // clock

using namespace std;

// https://stackoverflow.com/questions/72807569/set-default-value-of-unordered-map-if-key-doesnt-exist/72807851#72807851
// workaround to change default double value in map from 0 to -INFINITY
class dproxy {
    double value_;
public:
    dproxy(double value = -INFINITY)
    : value_{value} {}
    operator double () { return value_; }
    operator double const () const { return value_; }
};

unordered_map<char, int> BASE2ID, ID2BASE;
// unordered_map<char, int> ID2BASE;
string modelpath;
const int ALPHABET_SIZE = 5;
const int NUMMAT = 5;
const double EPSILON = pow(10, -2); // chose by eye just to distinguish real errors from numeric errors
const double TRAIN_THRESHOLD = pow(10, -6); // plotted one example as histogram and chose by eye
int kmerSize, C; // size of acceptable base characters, kmer size, min segment length, number of kmers
bool train, calcZ; // atrain
double a1, a2, p1, p2, p3, s1, s2, s3, e1, e2, e3, e4, i1, i2; // transition parameters
int T, N, K;
long S, NK;

// init model for unnormalised signals r9
// const string MODELPATH = "/home/yi98suv/projects/dynamont/data/template_median69pA_extended.model";
// init model for normalised signals r9
const string MODELPATH = "/home/yi98suv/projects/dynamont/data/norm_models/rna_r9.4_180mv_70bps_extended_stdev1.model";
// const string MODELPATH = "/home/yi98suv/projects/dynamont/data/norm_models/rna_r9.4_180mv_70bps_extended_stdev0_5.model";
// const string MODELPATH = "/home/yi98suv/projects/dynamont/data/norm_models/rna_r9.4_180mv_70bps_extended_stdev0_25.model";
const string TERM_STRING = "$";

// Asserts doubleing point compatibility at compile time
// necessary for INFINITY usage
static_assert(numeric_limits<double>::is_iec559, "IEEE 754 required");

void funcA(const int t, const int n, const int k, unordered_map<int, array<dproxy, NUMMAT>> &APSEI, unordered_map<int, array<dproxy, NUMMAT>> &logAPSEI, list<string>* segString, const int &N, const int &K);
void funcP(const int t, const int n, const int k, unordered_map<int, array<dproxy, NUMMAT>> &APSEI, unordered_map<int, array<dproxy, NUMMAT>> &logAPSEI, list<string>* segString, const int &N, const int &K);
void funcS(const int t, const int n, const int k, unordered_map<int, array<dproxy, NUMMAT>> &APSEI, unordered_map<int, array<dproxy, NUMMAT>> &logAPSEI, list<string>* segString, const int &N, const int &K);
void funcE(const int t, const int n, const int k, unordered_map<int, array<dproxy, NUMMAT>> &APSEI, unordered_map<int, array<dproxy, NUMMAT>> &logAPSEI, list<string>* segString, const int &N, const int &K);
void funcI(const int t, const int n, const int k, unordered_map<int, array<dproxy, NUMMAT>> &APSEI, unordered_map<int, array<dproxy, NUMMAT>> &logAPSEI, list<string>* segString, const int &N, const int &K);

/**
 * Fill the BASE2ID map with the base and ID pairs.
 */
void fillBASE2ID() {
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
 * C++ version 0.4 std::string style "itoa":
 * Contributions from Stuart Lowe, Ray-Yuan Sheu,
 * Rodrigo de Salvo Braz, Luc Gallant, John Maloney
 * and Brian Hunt
 * 
 * Converts a decimal to number to a number of base ALPHABET_SIZE.
 * TODO Works for base between 2 and 16 (included)
 * 
 * Returns kmer in reversed direction!
 * 
 * @param value input number in decimal to convert to base
*/
string itoa(const int &value) {
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

    // skip this so kmer is in 5' - 3' direction for output
    // reverse( buf.begin(), buf.end() );
    return buf;
}

/**
 * Converts a number of base ALPHABET_SIZE to a decimal number.
 * Works ONLY if ALPHABET_SIZE is smaller or equal to 10!
 * 
 * @param i input number in the given base as an array
*/
int toDeci(const int *i) {
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
inline int successingKmer(const int &currentKmer, const int &nextNt) {
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
inline int precessingKmer(const int &currentKmer, const int &priorNt) {
    return currentKmer/ALPHABET_SIZE + priorNt * int(pow(ALPHABET_SIZE, kmerSize-1));
}

/**
 * Convert the read sequence to a kmer sequence which is represented by integers.
 * 
 * @param seq read sequence
 * @param N length of the read sequence, number of nucleotides
 * @return kmer sequence in integer representation
*/
int* seq2kmer(const int *seq, const int &N) {
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
inline double log_normal_pdf(const double &x, const double &m, const double &s) {
    if(s==0) {
        return -INFINITY;
    }
    return -0.5*(log(2*M_PI*s*s)+((x - m)*(x - m)/(s*s)));
}

/**
 * Calculates the Hamming-Distance between two given kmers in their integer base representation
 */
inline int distanceSequenceKmer(const int &kmer_N, const int &kmer_K) { 
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
 * Return shifted & scaled log probability density for a given value and a given normal distribution
 */
inline double scoreSignalKmer(const double &signal_T, const double &mean, const double &stdev) {
    // return 2*(log_normal_pdf(signal_T, mean, stdev) + 6);
    return log_normal_pdf(signal_T, mean, stdev);
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
inline double score(const double &signal_T, const int &kmer_N, const int &kmer_K, const double &affineScale, const vector<tuple<double, double>> &model) {
    tuple<double, double> kmerModel_N = model[kmer_N];
    tuple<double, double> kmerModel_K = model[kmer_K];

    double scoreNT = scoreSignalKmer(signal_T, get<0>(kmerModel_N), get<1>(kmerModel_N));
    double scoreKT = scoreSignalKmer(signal_T, get<0>(kmerModel_K), get<1>(kmerModel_K));
    double scoreNK = -distanceSequenceKmer(kmer_N, kmer_K) * affineScale; // log(exp(-HD * affineCost)) = -HD * affineCost

    return scoreNT + scoreKT + scoreNK;
}

/**
 * Return combined log probability density for a given value and a given normal distribution
 *
 * @param signal point to calculate probability density
 * @param kmer key for the model kmer:(mean, stdev) map
 * @param model map containing kmers as keys and (mean, stdev) tuples as values
 * @return log probability density value for x in the given normal distribution
 */
inline double score(const double &signal_T, const int &kmer_N, const int &kmer_K, const vector<tuple<double, double>> &model) {
    tuple<double, double> kmerModel_N = model[kmer_N];
    tuple<double, double> kmerModel_K = model[kmer_K];

    double scoreNT = scoreSignalKmer(signal_T, get<0>(kmerModel_N), get<1>(kmerModel_N));
    double scoreKT = scoreSignalKmer(signal_T, get<0>(kmerModel_K), get<1>(kmerModel_K));
    double scoreNK = -distanceSequenceKmer(kmer_N, kmer_K); // log(exp(-HD)) = -HD

    return scoreNT + scoreKT + scoreNK;
}

/**
 * Return score for erroneous signal
 *
 * @param signal point to calculate probability density
 * @return log probability density value for signal_T
 */
// inline double errorScore(const double &signal_T) {
//     // TODO make variable for user input
//     if (signal_T > 4 || signal_T < -4) { // normalised
//         return 0;
//     }
//     return -INFINITY;
// }

// https://en.wikipedia.org/wiki/Log_probability
/**
 * Calculate addition of a+b in log space as efficiently as possible
 *
 * @param a first value
 * @param b second value
 * @return log(exp(a) + exp(b))
 */
inline double logPlus(const double &x, const double &y) {
    if (isinf(x) && isinf(y)) {
        return x;
    }
    if (x>=y){
        return x + log1p(exp(y-x));
    }
    return y + log1p(exp(x-y));
}

/**
 * Preprocess on signal and nucleotides to select plausible indices
 * Speeds up the DP
 */
// void preProcTN(const double *sig, const int *kmer_seq, vector<int> &allowedKeys, const int &T, const int &N, const vector<tuple<double, double>> &model) {
//     int e = 0;
//     // const double THRESHOLD = log(pow(10, -30)); // unnormalised
//     // const double THRESHOLD = log(0.1); // normalised
//     const double THRESHOLD = log(0); // normalised
//     allowedKeys.push_back(0);
//     tuple<double, double> kmerModel;
//     for(int n=1; n<N; n++) { // n=1
//         auto& [mean, stdev] = model[kmer_seq[n-1]];
//         for(int t=1; t<T; t++){ // t=1 dense preprocessing
//         // for(int t=max(0, n*35-5000); t<min(T-1, n*35+5000); t++){ // banded preprocessing
//             // store allowed indices
//             if(scoreSignalKmer(sig[t-1], mean, stdev) > THRESHOLD) {
//                 allowedKeys.push_back(t*N+n);
//             }
//             if(errorScore(sig[t-1]) > -INFINITY) {
//                 allowedKeys.push_back(t*N+n);
//                 e++;
//             }
//         }
//     }
//     cerr<<"sensor error rate: "<<e/float(T*N)<<endl;
//     sort(allowedKeys.begin(), allowedKeys.end());
// }

/**
 * Preprocess on signal and polishing kmers to select plausible indices
 * Speeds up the DP
 */
void preProcTK(const double *sig, const int *kmer_seq, vector<int> &allowedKeys, const int &T, const int &K, const vector<tuple<double, double>> &model) {
    // int e = 0;
    // TODO make variable for user
    const double THRESHOLD = log(0.15); // normalised
    // const double THRESHOLD = log(0); // normalised
    // init first column with all ks possible
    for(int k=0; k<K; k++){
        allowedKeys.push_back(k);
        allowedKeys.push_back((T-1)*K+k);
    }
    tuple<double, double> kmerModel;
    for(int k=0; k<K; k++) {
        auto& [mean, stdev] = model[k];
        for(int t=1; t<T-1; t++){
            if(scoreSignalKmer(sig[t-1], mean, stdev) > THRESHOLD) { [[likely]]
                allowedKeys.push_back(t*K+k);
                for (int prevK=0; prevK<ALPHABET_SIZE; prevK++) {
                    allowedKeys.push_back((t-1)*K+precessingKmer(k, prevK));
                    allowedKeys.push_back((t+1)*K+successingKmer(k, prevK));
                }
            }
                // add a range to fix weird transitions
            // } else if(errorScore(sig[t-1]) > -INFINITY) {
            //     allowedKeys.push_back(t*K+k);
            //     e++;
            // }
        }
    }
    // cerr<<"sensor error rate: "<<e/float(T*K)<<endl;
    sort(allowedKeys.begin(), allowedKeys.end());
    allowedKeys.erase(unique(allowedKeys.begin(), allowedKeys.end()), allowedKeys.end());
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
void logF(const double *sig, const int *kmer_seq, unordered_map<int, array<dproxy, NUMMAT>> &forAPSEI, vector<int> &allowedKeys, const int &T, const int &N, const int &K, const vector<tuple<double, double>> &model){
    double a, p, s, e, i, extendScore, openScore;
    int preKmer;
    int t, k; // , n
    long tNK, nK;
    // int idx, D1;
    //decleration of vector iterator
    for(vector<int>::iterator iter = allowedKeys.begin(); iter<allowedKeys.end(); iter++){
        // t = *iter / N;
        // n = *iter % N;
        // for(int k=0; k<K; k++){
        t = *iter / K;
        k = *iter % K;
        // for(int t=max(0, n*35-5000); t<min(T-1, n*35+5000); t++){ // banded preprocessing
        for(int n=0; n<N; n++){
            // // further sparsify
            // if (extendScore < log(pow(10, -5))) {
            //     continue;
            // }
            a=-INFINITY;
            p=-INFINITY;
            s=-INFINITY;
            e=-INFINITY;
            i=-INFINITY;
            if(t==0 && n==0) { [[unlikely]]
                e=0;
            }
            tNK = (t-1)*NK;
            nK = n*K;
            if(t>0 && n>0) { [[likely]]
                // precalc expensive values
                extendScore = score(sig[t-1], kmer_seq[n-1], k, model);
                openScore = score(sig[t-1], kmer_seq[n-1], k, 0.05, model);

                // only additions here
                for (int prevK=0; prevK<ALPHABET_SIZE; prevK++) {
                    preKmer = precessingKmer(k, prevK);
                    a=logPlus(a, forAPSEI[tNK+nK-K+preKmer][3] + a1 + openScore);
                    a=logPlus(a, forAPSEI[tNK+nK-K+preKmer][4] + a2 + openScore);
                
                    p=logPlus(p, forAPSEI[tNK+nK+preKmer][2] + p1 + openScore);
                    p=logPlus(p, forAPSEI[tNK+nK+preKmer][3] + p2 + openScore);
                    p=logPlus(p, forAPSEI[tNK+nK+preKmer][4] + p3 + openScore);
                }
                s=logPlus(s, forAPSEI[tNK+nK-K+k][1] + s1 + openScore);
                s=logPlus(s, forAPSEI[tNK+nK-K+k][3] + s2 + openScore);
                s=logPlus(s, forAPSEI[tNK+nK-K+k][4] + s3 + openScore);

                e=logPlus(e, forAPSEI[tNK+nK+k][0] + e1 + extendScore);
                e=logPlus(e, forAPSEI[tNK+nK+k][1] + e2 + extendScore);
                e=logPlus(e, forAPSEI[tNK+nK+k][2] + e3 + extendScore);
                e=logPlus(e, forAPSEI[tNK+nK+k][3] + e4 + extendScore);
                
                i=logPlus(i, forAPSEI[tNK-NK+nK-K+k][3] + i1 + openScore);
                i=logPlus(i, forAPSEI[tNK-NK+nK-K+k][4] + i2 + openScore);
            }
            forAPSEI[tNK+NK+nK+k] = {a, p, s, e, i};
        }
    }
}

// /**
//  * Calculate backward matrices using logarithmic values
//  *
//  * @param sig ONT raw signal with pA values
//  * @param seq nucleotide sequence represented by the ONT signal
//  * @param T length of the ONT raw signal + 1
//  * @param N length of nucleotide sequence + 1
//  * @param model map containing kmers as keys and (mean, stdev) tuples as values
//  */
void logB(const double *sig, const int *kmer_seq, unordered_map<int, array<dproxy, NUMMAT>> &backAPSEI, vector<int> &allowedKeys, const int &T, const int &N, const int &K, const vector<tuple<double, double>> &model){
    double a, p, s, e, i, pv, sc;
    int sucKmer;
    int t, k;
    long tNKnK;

    for(vector<int>::reverse_iterator iter = allowedKeys.rbegin(); iter<allowedKeys.rend(); ++iter){
        // t = *iter / N;
        // n = *iter % N;
        // for(int k=K-1; k>=0; k--){
        t = *iter / K;
        k = *iter % K;
        for(int n=N-1; n>=0; n--){
            a=-INFINITY;
            p=-INFINITY;
            s=-INFINITY;
            e=-INFINITY;
            i=-INFINITY;
            if(t==T-1 && n==N-1) { [[unlikely]]
                e=0;
            }
            // precalc expensive value    
            tNKnK = (t+1)*(NK)+n*K;
            if (t<T-1) { [[likely]]
                if (n>0) { [[likely]]
                    // precalc expensive value    
                    pv=backAPSEI[tNKnK+k][3];
                    sc=score(sig[t], kmer_seq[n-1], k, model);

                    // only addition here
                    a=logPlus(a, pv + e1 + sc);
                    p=logPlus(p, pv + e2 + sc);
                    s=logPlus(s, pv + e3 + sc);
                    e=logPlus(e, pv + e4 + sc);
                    for (int nextK=0; nextK<ALPHABET_SIZE; nextK++){
                        sucKmer = successingKmer(k, nextK);
                        pv=backAPSEI[tNKnK+sucKmer][1];
                        sc=score(sig[t], kmer_seq[n-1], sucKmer, 0.05, model);
                        s=logPlus(s, pv + p1 + sc);
                        e=logPlus(e, pv + p2 + sc);
                        i=logPlus(i, pv + p3 + sc);
                    }
                }
                if (n<N-1) { [[likely]]
                    pv=backAPSEI[tNKnK+K+k][2]; //[(t+1)*(NK)+(n+1)*K+k][2];
                    sc=score(sig[t], kmer_seq[n], k, 0.05, model);
                    p=logPlus(p, pv + s1 + sc);
                    e=logPlus(e, pv + s2 + sc);
                    i=logPlus(i, pv + s3 + sc);
                    for (int nextK=0; nextK<ALPHABET_SIZE; nextK++){
                        sucKmer = successingKmer(k, nextK);
                        pv=backAPSEI[tNKnK+K+sucKmer][0];
                        sc=score(sig[t], kmer_seq[n], sucKmer, 0.05, model);
                        e=logPlus(e, pv + a1 + sc);
                        i=logPlus(i, pv + a2 + sc);
                    }
                }
            }
            if (t>0 && n<N-1) { [[likely]]
                pv=backAPSEI[tNKnK-NK+K+k][4]; //[t*(NK)+(n+1)*K+k][4];
                sc=score(sig[t-1], kmer_seq[n], k, 0.05, model);
                e=logPlus(e, pv + i1 + sc);
                i=logPlus(i, pv + i2 + sc);
            }
            // backAPSEI[t*(NK)+n*K+k] = {a, p, s, e, i};
            backAPSEI[tNKnK-NK+k] = {a, p, s, e, i};
        }
    }
}

/**
 * Calculate the logarithmic probability matrix
 *
 * @return matrix containing logarithmic probabilities for segment borders
 */
void logP(unordered_map<int, array<dproxy, NUMMAT>> &logAPSEI, unordered_map<int, array<dproxy, NUMMAT>> &forAPSEI, unordered_map<int, array<dproxy, NUMMAT>> &backAPSEI, const double &Z, vector<int> &allowedKeys, const int &N, const int &K) {
    int t = 0, k;
    long idx;
    for(vector<int>::iterator iter = allowedKeys.begin(); iter<allowedKeys.end(); iter++){
        t = *iter / K;
        k = *iter % K;
        for(int n=0; n<N; n++){
            idx = t*(NK)+n*K+k;
            for(int mat=0; mat<NUMMAT; mat++) {
                logAPSEI[idx][mat] = forAPSEI[idx][mat] + backAPSEI[idx][mat] - Z;
            }
        }
    }

    // // just testing
    // for(int a=0; a<t; a++) {
    //     double s = -INFINITY;
    //     for(int n=0; n<N; n++) {
    //         for (int k=0; k<K; k++) {
    //             idx = a*(NK)+n*K+k; 
    //             s = logPlus(s, logPlus(logAPSEI[idx][0], logPlus(logAPSEI[idx][1], logPlus(logAPSEI[idx][2], logPlus(logAPSEI[idx][3], logAPSEI[idx][4])))));
    //         }
    //     }
    //     cout<<"this should be 1: "<<exp(s)<<endl;
    // }
}

/**
 * Calculate the maximum a posteriori path through LP
 *
 */
list<string> getBorders(unordered_map<int, array<dproxy, NUMMAT>> &logAPSEI, vector<int> &allowedKeys, const int &T, const int &N, const int &K){
    unordered_map<int, array<dproxy, NUMMAT>> APSEI;
    array<dproxy, NUMMAT> curIdx;
    double a, p, s, e, i;
    int t, k;
    long tNK, nK;
    for(vector<int>::iterator iter = allowedKeys.begin(); iter<allowedKeys.end(); iter++){
        t = *iter / K;
        k = *iter % K;
        for(int n=0; n<N; n++){
            a=-INFINITY;
            p=-INFINITY;
            s=-INFINITY;
            e=-INFINITY;
            i=-INFINITY;
            if(t==0 && n==0){ [[unlikely]]
                e=0;
            }
            // precalc expensive multiplications
            tNK = (t-1)*(NK);
            nK = n*K;

            // mostly addition left
            if(t>0 && n>0) { [[likely]]
                curIdx = logAPSEI[tNK+NK+nK+k];
                for (int prevK=0; prevK<ALPHABET_SIZE; prevK++) {
                    a=max(a, APSEI[tNK+nK-K+precessingKmer(k, prevK)][3] + curIdx[0]);
                    a=max(a, APSEI[tNK+nK-K+precessingKmer(k, prevK)][4] + curIdx[0]);

                    p=max(p, APSEI[tNK+nK+precessingKmer(k, prevK)][2] + curIdx[1]);
                    p=max(p, APSEI[tNK+nK+precessingKmer(k, prevK)][3] + curIdx[1]);
                    p=max(p, APSEI[tNK+nK+precessingKmer(k, prevK)][4] + curIdx[1]);
                }
                s=max(s, APSEI[tNK+nK-K+k][1] + curIdx[2]);
                s=max(s, APSEI[tNK+nK-K+k][3] + curIdx[2]);
                s=max(s, APSEI[tNK+nK-K+k][4] + curIdx[2]);

                e=max(e, APSEI[tNK+nK+k][0] + curIdx[3]);
                e=max(e, APSEI[tNK+nK+k][1] + curIdx[3]);
                e=max(e, APSEI[tNK+nK+k][2] + curIdx[3]);
                e=max(e, APSEI[tNK+nK+k][3] + curIdx[3]);

                i=max(i, APSEI[tNK+NK+nK-K+k][3] + curIdx[4]);
                i=max(i, APSEI[tNK+NK+nK-K+k][4] + curIdx[4]);
            }
            APSEI[tNK+NK+nK+k] = {a, p, s, e, i};
        }
    }
    list<string> segString;
    funcE(T-1, N-1, K-1, APSEI, logAPSEI, &segString, N, K);
    return segString;
}

void funcA(const int t, const int n, const int k, unordered_map<int, array<dproxy, NUMMAT>> &APSEI, unordered_map<int, array<dproxy, NUMMAT>> &logAPSEI, list<string>* segString, const int &N, const int &K){
    double score = APSEI[t*(NK)+n*K+k][0];
    if (t<=1 && n<=1){ // Start value
        segString->push_front("M"+to_string(0)+","+to_string(0)+","+itoa(k)+";"); // n-1 because N is 1 larger than the sequences
        return;
    }
    for (int prevK=0; prevK<ALPHABET_SIZE; prevK++) {
        if (t>0 && n>0 && score == APSEI[(t-1)*(NK)+(n-1)*K+precessingKmer(k, prevK)][3] + logAPSEI[t*(NK)+n*K+k][0]){
            segString->push_front("M"+to_string(n-1)+","+to_string(t-1)+","+itoa(k)+";");
            return funcE(t-1, n-1, precessingKmer(k, prevK), APSEI, logAPSEI, segString, N, K);
        }
        if (t>0 && n>0 && score == APSEI[(t-1)*(NK)+(n-1)*K+precessingKmer(k, prevK)][4] + logAPSEI[t*(NK)+n*K+k][0]){
            segString->push_front("M"+to_string(n-1)+","+to_string(t-1)+","+itoa(k)+";");
            return funcI(t-1, n-1, precessingKmer(k, prevK), APSEI, logAPSEI, segString, N, K);
        }
    }
}

void funcE(const int t, const int n, const int k, unordered_map<int, array<dproxy, NUMMAT>> &APSEI, unordered_map<int, array<dproxy, NUMMAT>> &logAPSEI, list<string>* segString, const int &N, const int &K){
    double score = APSEI[t*(NK)+n*K+k][3];
    if (t>0 && n>0 && score == APSEI[(t-1)*(NK)+n*K+k][3] + logAPSEI[t*(NK)+n*K+k][3]){
        return funcE(t-1, n, k, APSEI, logAPSEI, segString, N, K);
    }
    if (t>0 && n>0 && score == APSEI[(t-1)*(NK)+n*K+k][0] + logAPSEI[t*(NK)+n*K+k][3]){
        return funcA(t-1, n, k, APSEI, logAPSEI, segString, N, K);
    }
    if (t>0 && n>0 && score == APSEI[(t-1)*(NK)+n*K+k][2] + logAPSEI[t*(NK)+n*K+k][3]){
        return funcS(t-1, n, k, APSEI, logAPSEI, segString, N, K);
    }
    if (t>0 && n>0 && score == APSEI[(t-1)*(NK)+n*K+k][4] + logAPSEI[t*(NK)+n*K+k][3]){
        return funcI(t-1, n, k, APSEI, logAPSEI, segString, N, K);
    }
}

void funcP(const int t, const int n, const int k, unordered_map<int, array<dproxy, NUMMAT>> &APSEI, unordered_map<int, array<dproxy, NUMMAT>> &logAPSEI, list<string>* segString, const int &N, const int &K){
    double score = APSEI[t*(NK)+n*K+k][1];
    for (int prevK=0; prevK<ALPHABET_SIZE; prevK++) {
        if (t>0 && n>0 && score == APSEI[(t-1)*(NK)+n*K+precessingKmer(k, prevK)][2] + logAPSEI[t*(NK)+n*K+k][1]){
            segString->push_front("P"+to_string(n)+","+to_string(t-1)+","+itoa(k)+";");
            return funcS(t-1, n, precessingKmer(k, prevK), APSEI, logAPSEI, segString, N, K);
        }
        if (t>0 && n>0 && score == APSEI[(t-1)*(NK)+n*K+precessingKmer(k, prevK)][3] + logAPSEI[t*(NK)+n*K+k][1]){
            segString->push_front("P"+to_string(n)+","+to_string(t-1)+","+itoa(k)+";");
            return funcE(t-1, n, precessingKmer(k, prevK), APSEI, logAPSEI, segString, N, K);
        }
        if (t>0 && n>0 && score == APSEI[(t-1)*(NK)+n*K+precessingKmer(k, prevK)][4] + logAPSEI[t*(NK)+n*K+k][1]){
            segString->push_front("P"+to_string(n)+","+to_string(t-1)+","+itoa(k)+";");
            return funcI(t-1, n, precessingKmer(k, prevK), APSEI, logAPSEI, segString, N, K);
        }
    }
}

void funcS(const int t, const int n, const int k, unordered_map<int, array<dproxy, NUMMAT>> &APSEI, unordered_map<int, array<dproxy, NUMMAT>> &logAPSEI, list<string>* segString, const int &N, const int &K){
    double score = APSEI[t*(NK)+n*K+k][2];
    if (t>0 && n>0 && score == APSEI[(t-1)*(NK)+(n-1)*K+k][3] + logAPSEI[t*(NK)+n*K+k][2]){
        segString->push_front("S"+to_string(n-1)+","+to_string(t-1)+","+itoa(k)+";");
        return funcE(t-1, n-1, k, APSEI, logAPSEI, segString, N, K);
    }
    if (t>0 && n>0 && score == APSEI[(t-1)*(NK)+(n-1)*K+k][1] + logAPSEI[t*(NK)+n*K+k][2]){
        segString->push_front("S"+to_string(n-1)+","+to_string(t-1)+","+itoa(k)+";");
        return funcP(t-1, n-1, k, APSEI, logAPSEI, segString, N, K);
    }
    if (t>0 && n>0 && score == APSEI[(t-1)*(NK)+(n-1)*K+k][4] + logAPSEI[t*(NK)+n*K+k][2]){
        segString->push_front("S"+to_string(n-1)+","+to_string(t-1)+","+itoa(k)+";");
        return funcI(t-1, n-1, k, APSEI, logAPSEI, segString, N, K);
    }
}

void funcI(const int t, const int n, const int k, unordered_map<int, array<dproxy, NUMMAT>> &APSEI, unordered_map<int, array<dproxy, NUMMAT>> &logAPSEI, list<string>* segString, const int &N, const int &K){
    double score = APSEI[t*(NK)+n*K+k][4];
    if (t>0 && n>0 && score == APSEI[t*(NK)+(n-1)*K+k][4] + logAPSEI[t*(NK)+n*K+k][4]){
        segString->push_front("I"+to_string(n-1)+","+to_string(t)+","+itoa(k)+";");
        return funcI(t, n-1, k, APSEI, logAPSEI, segString, N, K);
    }
    if (t>0 && n>0 && score == APSEI[t*(NK)+(n-1)*K+k][3] + logAPSEI[t*(NK)+n*K+k][4]){
        segString->push_front("I"+to_string(n-1)+","+to_string(t)+","+itoa(k)+";");
        return funcE(t, n-1, k, APSEI, logAPSEI, segString, N, K);
    }
}

/**
 * Read the normal distribution parameters from a given TSV file
 *
 * @param file path to the TSV file containing the parameters
 * @param model kmer model to fill
 */
void readKmerModel(const string &file, vector<tuple<double, double>>& model) {
    ifstream inputFile(file);
    string line, kmer, tmp;
    double mean, stdev;
    getline(inputFile, line);
    while(getline(inputFile, line)) { // read line
        stringstream buffer(line); // parse line to stringstream for getline
        getline(buffer, kmer, '\t');
        // legacy models are stored from 3' - 5'
        // all other (basically new) models are stored in 5' - 3'
        reverse(kmer.begin(), kmer.end()); // 5-3 -> 3-5 orientation
        getline(buffer, tmp, '\t'); // level_mean
        mean = atof(tmp.c_str());
        getline(buffer, tmp, '\t'); // level_stdv
        stdev = atof(tmp.c_str());
        model[kmer2int(kmer)]=make_tuple(mean, stdev);
    }
    inputFile.close();
}

/**
 * Train transition parameter with baum welch algorithm
*/
tuple<double, double, double, double, double, double, double, double, double, double, double, double, double, double> trainTransition(const double *sig, const int *kmer_seq, unordered_map<int, array<dproxy, NUMMAT>> &forAPSEI, unordered_map<int, array<dproxy, NUMMAT>> &backAPSEI, vector<int> &allowedKeys, const int &T, const int &N, const int &K, vector<tuple<double, double>> &model) {
    // Transition parameters
    double newa1=-INFINITY, newa2=-INFINITY, newp1=-INFINITY, newp2=-INFINITY, newp3=-INFINITY, news1=-INFINITY, news2=-INFINITY, news3=-INFINITY, newe1=-INFINITY, newe2=-INFINITY, newe3=-INFINITY, newe4=-INFINITY, newi1=-INFINITY, newi2=-INFINITY; // transition parameters
    int sucKmer;
    double sc;
    int t, k;
    long tNK, nK;

    for(vector<int>::iterator iter = allowedKeys.begin(); iter<allowedKeys.end(); iter++){
        t = *iter / K;
        k = *iter % K;
        for(int n=0; n<N; n++){
            // precalc expensive multiplications
            tNK = t*NK;
            nK = n*K;

            // mostly addition from here
            if (t<T-1) {
                if (n>0) {
                    sc = score(sig[t], kmer_seq[n-1], k, model);
                    newe1 = logPlus(newe1, forAPSEI[tNK+nK+k][0] + e1 + sc                 + backAPSEI[tNK+NK+nK+k][3]);
                    newe2 = logPlus(newe2, forAPSEI[tNK+nK+k][1] + e2 + sc                 + backAPSEI[tNK+NK+nK+k][3]);
                    newe3 = logPlus(newe3, forAPSEI[tNK+nK+k][2] + e3 + sc                 + backAPSEI[tNK+NK+nK+k][3]);
                    newe4 = logPlus(newe4, forAPSEI[tNK+nK+k][3] + e4 + sc                 + backAPSEI[tNK+NK+nK+k][3]);

                    for (int nextK=0; nextK<ALPHABET_SIZE; nextK++){
                        sucKmer = successingKmer(k, nextK);
                        sc = score(sig[t], kmer_seq[n-1], sucKmer, 0.05, model);
                        newp1 = logPlus(newp1, forAPSEI[tNK+nK+k][2] + p1 + sc + backAPSEI[tNK+NK+nK+sucKmer][1]);
                        newp2 = logPlus(newp2, forAPSEI[tNK+nK+k][3] + p2 + sc + backAPSEI[tNK+NK+nK+sucKmer][1]);
                        newp3 = logPlus(newp3, forAPSEI[tNK+nK+k][4] + p3 + sc + backAPSEI[tNK+NK+nK+sucKmer][1]);
                    }
                }

                if (n<N-1) {
                    sc = score(sig[t], kmer_seq[n], k, 0.05, model);
                    news1 = logPlus(news1, forAPSEI[tNK+nK+k][1] + s1 + sc + backAPSEI[tNK+NK+nK+K+k][2]);
                    news2 = logPlus(news2, forAPSEI[tNK+nK+k][3] + s2 + sc + backAPSEI[tNK+NK+nK+K+k][2]);
                    news3 = logPlus(news3, forAPSEI[tNK+nK+k][4] + s3 + sc + backAPSEI[tNK+NK+nK+K+k][2]);
                    
                    for (int nextK=0; nextK<ALPHABET_SIZE; nextK++){
                        sucKmer = successingKmer(k, nextK);
                        sc = score(sig[t], kmer_seq[n], sucKmer, 0.05, model);
                        newa1 = logPlus(newa1, forAPSEI[tNK+nK+k][3] + a1 + sc + backAPSEI[tNK+NK+nK+K+k+sucKmer][0]);
                        newa2 = logPlus(newa2, forAPSEI[tNK+nK+k][4] + a2 + sc + backAPSEI[tNK+NK+nK+K+k+sucKmer][0]);
                    }
                }
            }

            if (t>0 && n<N-1) {
                sc = score(sig[t-1], kmer_seq[n], k, 0.05, model);
                newi1 = logPlus(newi1, forAPSEI[tNK+nK+k][3] + i1 + sc + backAPSEI[tNK+nK+K+k][4]);
                newi2 = logPlus(newi2, forAPSEI[tNK+nK+k][4] + i2 + sc + backAPSEI[tNK+nK+K+k][4]);
            }
        }
    }

    // average over the number of transitions
    newe1=exp(newe1-newe1); // Aa
    double Ae = logPlus(logPlus(logPlus(newa1, news2), logPlus(newe4, newi1)), newp2);
    newa1=exp(newa1-Ae);
    news2=exp(news2-Ae);
    newe4=exp(newe4-Ae);
    newi1=exp(newi1-Ae);
    newp2=exp(newp2-Ae);
    double As = logPlus(newe3, newp1);
    newe3=exp(newe3-As);
    newp1=exp(newp1-As);
    double Ap = logPlus(newe2, news1);
    newe2=exp(newe2-Ap);
    news1=exp(news1-Ap);
    double Ai = logPlus(logPlus(a2, i2), logPlus(p3, s3));
    newa2=exp(newa2-Ai);
    newi2=exp(newi2-Ai);
    newp3=exp(newp3-Ai);
    news3=exp(news3-Ai);

    return tuple<double, double, double, double, double, double, double, double, double, double, double, double, double, double>({newa1, newa2, newp1, newp2, newp3, news1, news2, news3, newe1, newe2, newe3, newe4, newi1, newi2});
}

/**
 * Train emission parameter with baum welch algorithm
*/
tuple<double*, double*> trainEmission(const double* sig, unordered_map<int, array<dproxy, NUMMAT>> &logAPSEI, vector<int> &allowedKeys, const int &T, const int &N, const int &K) {
    // Emission
    // https://courses.grainger.illinois.edu/ece417/fa2021/lectures/lec15.pdf
    // https://f.hubspotusercontent40.net/hubfs/8111846/Unicon_October2020/pdf/bilmes-em-algorithm.pdf
    // gamma_t(i) is the probability of being in state i at time t
    // gamma for state M - expected number of transitions of M at given time (T) for all latent states (kmers)
    // unordered_map<int, dproxy> G, s;
    // double forg, backg;

    // for(vector<int>::iterator iter = allowedKeys.begin(); iter<allowedKeys.end(); iter++){
    //     t = *iter / K;
    //     k = *iter % K;
    //     for(int n=0; n<N; n++){
    //         idx = t*(NK)+n*K+k;
    //         // forg=-INFINITY;
    //         // backg=-INFINITY;
    //         for(int i=0; i<NUMMAT; i++){
    //             // if(isnan(forAPSEI[idx][i] + backAPSEI[idx][i])) {
    //             //     cerr<<t<<","<<n<<","<<k<<", "<<idx<<", "<<i<<", "<<forAPSEI[idx][i]<<backAPSEI[idx][i]<<endl;
    //             //     exit(1);
    //             // }
    //             G[idx] = logPlus(G[idx], forAPSEI[idx][i] + backAPSEI[idx][i]);
    //         }
    //         // for(int i=0; i<NUMMAT; i++){
    //         //     forg = logPlus(forg, forAPSEI[idx][i]);
    //         //     backg = logPlus(backg, backAPSEI[idx][i]);
    //         // }
    //         // G[idx] = forg + backg;
    //         s[t] = logPlus(s[t], G[idx]);
    //     }
    // }
    // for(vector<int>::iterator iter = allowedKeys.begin(); iter<allowedKeys.end(); iter++) {
    //     t = *iter / K;
    //     k = *iter % K;
    //     for(int n=0; n<N; n++){
    //         if (!isinf(s[t])) {
    //             idx = t*(NK)+n*K+k;
    //             G[t*(NK)+n*K+k] = G[t*(NK)+n*K+k] - s[t];
    //         }
    //     }
    // }

    // decimal space
    int k, t;
    long idx;
    double* means = new double[K];
    double* stdevs = new double[K];
    double w = -INFINITY;
    double* normFactorT = new double[K];

    // int* normFactorN = new int[K];
    // for(int n=0; n<N; n++){
    //     normFactorN[kmer_seq[n]]++;
    // }

    for(vector<int>::iterator iter = allowedKeys.begin(); iter<allowedKeys.end(); iter++){
        t = *iter / K;
        k = *iter % K;
        w = -INFINITY;
        for(int n=0; n<N; n++){
            idx = t*(NK)+n*K+k;
            if(t>0) { [[likely]]
                // means[k] += exp(G[idx]) * sig[t-1];
                w = logPlus(logPlus(logPlus(logAPSEI[idx][0], logAPSEI[idx][1]), logPlus(logAPSEI[idx][2], logAPSEI[idx][3])), logPlus(logAPSEI[idx][4], w));
            }
            // normFactorT[k] += exp(G[idx]);
        }
        w = exp(w);
        means[k] += w * sig[t-1];
        normFactorT[k] += w;
    }

    for(int k=0; k<K; k++){
        means[k] = means[k] / normFactorT[k];
        // means[k] = means[k] / T;
        // cerr<<normFactorT[k]/T<<",";
    }
    // exit(1);

    // Emission (stdev of kmers)
    for(vector<int>::iterator iter = allowedKeys.begin(); iter<allowedKeys.end(); iter++){
        t = *iter / K;
        k = *iter % K;
        // skip kmers that have a very low weight, they will not be updated
        if(normFactorT[k]/T < TRAIN_THRESHOLD){
            continue;
        }
        w = -INFINITY;
        for(int n=0; n<N; n++){
            idx = t*(NK)+n*K+k;
            if (t>0) { [[likely]]
                w = logPlus(logPlus(logPlus(logAPSEI[idx][0], logAPSEI[idx][1]), logPlus(logAPSEI[idx][2], logAPSEI[idx][3])), logPlus(logAPSEI[idx][4], w));
                // stdevs[k] += exp(logPlus(logPlus(logPlus(logAPSEI[idx][0], logAPSEI[idx][1]), logPlus(logAPSEI[idx][2], logAPSEI[idx][3])), logAPSEI[idx][4])) * pow(sig[t-1] - means[k], 2.);
            }
        }
        stdevs[k] += exp(w) * pow(sig[t-1] - means[k], 2.);
    }
    for(int k=0; k<K; k++){
        // stdevs[k] = sqrt(stdevs[k]);
        stdevs[k] = sqrt(stdevs[k] / normFactorT[k]);
    }

    // delete[] normFactorT;
    return tuple<double*, double*>({means, stdevs});
}

void printTrainedTransitionParams(const double *sig, const int *kmer_seq, unordered_map<int, array<dproxy, NUMMAT>> &forAPSEI, unordered_map<int, array<dproxy, NUMMAT>> &backAPSEI, unordered_map<int, array<dproxy, NUMMAT>> &logAPSEI, vector<int> &allowedKeys, const int &T, const int &N, const int &K, vector<tuple<double, double>> &model) {

    auto [a1, a2, p1, p2, p3, s1, s2, s3, e1, e2, e3, e4, i1, i2] = trainTransition(sig, kmer_seq, forAPSEI, backAPSEI, allowedKeys, T, N, K, model);
    cout<<"a1:"<<a1<<";a2:"<<a2<<";p1:"<<p1<<";p2:"<<p2<<";p2:"<<p3<<";s1:"<<s1<<";s2:"<<s2<<";s3:"<<s3<<";e1:"<<e1<<";e2:"<<e2<<";e3:"<<e3<<";e4:"<<e4<<";i1:"<<i1<<";i2:"<<i2<<endl;

    auto [newMeans, newStdevs] = trainEmission(sig, logAPSEI, allowedKeys, T, N, K);
    for (int k=0; k<K; k++){
        if (newStdevs[k]!=0.0){
            cout<<itoa(k)<<":"<<newMeans[k]<<","<<newStdevs[k]<<";";
        }
    }
    cout<<endl;
    delete[] newMeans;
    delete[] newStdevs;
}

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
    catch (const runtime_error& err) {
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
    // polishing dimension K = number of possible kmers
    K = int(pow(ALPHABET_SIZE, kmerSize)); // currently acceptable A, C, G, T, N
    vector<tuple<double, double>> model(K, make_tuple(-INFINITY, -INFINITY));
    readKmerModel(modelpath, model);
    string signal;
    string read;

    while(1) {
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
        // cerr<<signal.length()<<endl;
        T = count(signal.begin(), signal.end(), ',')+2; // len(sig) + 1
        double* sig = new double[T-1];
        fill_n(sig, T-1, -INFINITY);
        string value;
        stringstream ss(signal);
        int i = 0;
        while(getline(ss, value, ',')) {
            sig[i++] = stod(value);
        }
        // process read N: convert string to int array
        N = read.size() + 1; // operate on base transitions
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

        int* kmer_seq = seq2kmer(seq, N-1);
        S = S;
        NK = N*K;

        // just for time measurement
        // typedef chrono::high_resolution_clock Clock;
        // typedef std::chrono::milliseconds milliseconds;

        cerr<<"T: "<<T<<endl;
        cerr<<"N: "<<N<<endl;
        cerr<<"K: "<<K<<endl;
        cerr<<"inputsize: "<<S<<endl;

        vector<int> allowedKeys;
        preProcTK(sig, kmer_seq, allowedKeys, T, K, model);
        cerr<<"dense: "<<allowedKeys.size()/float(T*K)<<", sparse: "<<1-(allowedKeys.size()/float(T*K))<<endl;
        unordered_map<int, array<dproxy, NUMMAT>> forAPSEI;
        
        // Clock::time_point t0 = Clock::now();
        logF(sig, kmer_seq, forAPSEI, allowedKeys, T, N, K, model);
        // Clock::time_point t1 = Clock::now();
        // cerr<<"Done Forward! "<<std::chrono::duration_cast<milliseconds>(t1 - t0).count()<<endl;

        unordered_map<int, array<dproxy, NUMMAT>> backAPSEI;
        logB(sig, kmer_seq, backAPSEI, allowedKeys, T, N, K, model);
        // cerr<<"Done Backward! "<<std::chrono::duration_cast<milliseconds>(Clock::now() - t1).count()<<endl;

        double Zf = -INFINITY;
        double Zb = -INFINITY;

        // init, log(1) for any k
        for(int k=0; k<K; k++){
            // Zf = logPlus(Zf, forAPSEI[(T-1)*(NK)+(N-1)*K+k][3]);
            Zf = logPlus(Zf, forAPSEI[S-1-k][3]);
            Zb = logPlus(Zb, backAPSEI[k][3]);
        }

        // Numeric error is scaled by input size, Z in forward and backward should match by some numeric error EPSILON
        if ((isinf(Zf) || isinf(Zb) || isnan(Zf) || isnan(Zb) || abs(Zf-Zb)/(S)>=EPSILON)) {
            cerr << fixed << showpoint;
            cerr << setprecision(20);
            cerr<<"Z values between matrices do not match! forZ: "<<Zf<<", backZ: "<<Zb<<", "<<abs(Zf-Zb)/(S)<<" > "<<EPSILON<<endl;
            cerr.flush();
            exit(11);
        }

        // cerr<<"forZ: "<<Zf<<", backZ: "<<Zb<<", "<<abs(Zf-Zb)/(S)<<" <! "<<EPSILON<<endl;

        if (calcZ){
            cout<<Zf<<"\n";
            cout.flush();
        } else {
            unordered_map<int, array<dproxy, NUMMAT>> logAPSEI;
            logP(logAPSEI, forAPSEI, backAPSEI, Zf, allowedKeys, N, K);

            // train both Transitions and Emissions
            if (train) {
                // t0 = Clock::now();
                printTrainedTransitionParams(sig, kmer_seq, forAPSEI, backAPSEI, logAPSEI, allowedKeys, T, N, K, model);
                // cerr<<"Done Training! "<<std::chrono::duration_cast<milliseconds>(Clock::now() - t0).count()<<endl;
                cout<<"Z:"<<Zf<<endl;
                cout.flush();

            } else {
                list<string> segString = getBorders(logAPSEI, allowedKeys, T, N, K);

                // testing stuff
                // int numCells = S;
                // int magnitude = 20;
                // int numDenseCells[magnitude] = {0}; // 10^-0, ..., 10^-(magnitude-1)
                // double maxVal = 0;
                // long idx = 0;

                // for(int t=0; t<T; t++){
                //     for(int n=0; n<N; n++){
                //         for(int k=0; k<K; k++){
                //             idx = t*(NK)+n*K+k;
                //             maxVal = exp(max(LPA[idx], max(max(LPE[idx], LPI[idx]), max(LPP[idx], LPS[idx]))));
                //             // maxVal = max(LPA[idx], max(LPP[idx], max(LPS[idx], max(LPE[idx], LPI[idx]))));
                            
                //             for(int i=0; i<magnitude; i++){
                //                 if(maxVal > pow(10, -i)) {
                //                     for(int j=i; j<magnitude; j++) {
                //                         numDenseCells[j]+=1;
                //                     }
                //                     break;
                //                 }
                //             }
                //         }
                //     }
                // }
                // for(int i=0; i<magnitude; i++){
                //     cout<<"10^-"<<i<<": "<<numDenseCells[i]<<", "<<numDenseCells[i]/double(numCells)<<endl;
                // }
                // cout<<numCells<<endl;

                for (auto const& seg : segString) {
                    cout<<seg;
                }
                cout<<endl;
                cout.flush();
            }
        }
        delete[] sig;
        delete[] seq;
        delete[] kmer_seq;
    }
    return 0;
}
