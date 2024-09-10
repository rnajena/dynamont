// author: Jannes Spangenberg
// e-mail: jannes.spangenberg@uni-jena.de
// github: https://github.com/JannesSP
// website: https://jannessp.github.io

#include <iostream>
#include <iomanip>
#include <fstream> // file io
#include <sstream> // file io
#include <string>
#include <algorithm> // std::sort, std::stable_sort
#include <numeric>   // std::iota
#include <vector>
#include <array>
#include <unordered_map>
#include <tuple>
#include <bits/stdc++.h> // reverse strings
#include <cmath> // exp, pow, log1p, INFINITY
#include <assert.h>
#include <stdlib.h>
#include "argparse.hpp"
#include "utils.hpp"

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

const int NUMMAT = 5;
int kmerSize, stepSize;
inline constexpr double EPSILON = 1e-8; // chose by eye just to distinguish real errors from numeric errors
double a1, a2, p1, p2, p3, s1, s2, s3, e1, e2, e3, e4, i1, i2; // transition parameters
size_t T, N, K, TNK, NK;

// Asserts doubleing point compatibility at compile time
// necessary for INFINITY usage
static_assert(numeric_limits<double>::is_iec559, "IEEE 754 required");

void funcA(const size_t t, const size_t n, const size_t k, unordered_map<size_t, array<dproxy, NUMMAT>> &APSEI, unordered_map<size_t, array<dproxy, NUMMAT>> &logAPSEI, list<string>* segString, const size_t K);
void funcP(const size_t t, const size_t n, const size_t k, unordered_map<size_t, array<dproxy, NUMMAT>> &APSEI, unordered_map<size_t, array<dproxy, NUMMAT>> &logAPSEI, list<string>* segString, const size_t K);
void funcS(const size_t t, const size_t n, const size_t k, unordered_map<size_t, array<dproxy, NUMMAT>> &APSEI, unordered_map<size_t, array<dproxy, NUMMAT>> &logAPSEI, list<string>* segString, const size_t K);
void funcE(const size_t t, const size_t n, const size_t k, unordered_map<size_t, array<dproxy, NUMMAT>> &APSEI, unordered_map<size_t, array<dproxy, NUMMAT>> &logAPSEI, list<string>* segString, const size_t K);
void funcI(const size_t t, const size_t n, const size_t k, unordered_map<size_t, array<dproxy, NUMMAT>> &APSEI, unordered_map<size_t, array<dproxy, NUMMAT>> &logAPSEI, list<string>* segString, const size_t K);

// ===============================================================
// ===============================================================
// ===================== Scoring calculations ====================
// ===============================================================
// ===============================================================

/**
 * Calculates the Hamming-Distance between two given kmers in their integer base representation
 */
inline int distanceSequenceKmer(const size_t kmer_N, const size_t kmer_K) {
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
 * @param signal_T point to calculate probability density
 * @param kmer_N key for kmer N the model kmer:(mean, stdev) map
 * @param kmer_K key for kmer K the model kmer:(mean, stdev) map
 * @param affineScale affince cost for NK comparison
 * @param model map containing kmers as keys and (mean, stdev) tuples as values
 * @return log probability density value for x in the given normal distribution
 */
inline double score(const double signal_T, const size_t kmer_N, const size_t kmer_K, const double affineScale, const vector<tuple<double, double>> &model) {
    // Access elements of the model tuple directly to avoid redundant tuple creation and overhead
    const auto &[mean_N, stdev_N] = model[kmer_N];
    const auto &[mean_K, stdev_K] = model[kmer_K];

    // Precompute the scores for the individual kmers and their distance
    const double scoreNT = log_normal_pdf(signal_T, mean_N, stdev_N);
    const double scoreKT = log_normal_pdf(signal_T, mean_K, stdev_K);
    const double scoreNK = -distanceSequenceKmer(kmer_N, kmer_K) * affineScale; // log(exp(-HD * affineCost)) = -HD * affineCost

    return scoreNT + scoreKT + scoreNK;
}

/**
 * Return combined log probability density for a given value and a given normal distribution
 *
 * @param signal_T point to calculate probability density
 * @param kmer_N key for kmer N the model kmer:(mean, stdev) map
 * @param kmer_K key for kmer K the model kmer:(mean, stdev) map
 * @param model map containing kmers as keys and (mean, stdev) tuples as values
 * @return log probability density value for x in the given normal distribution
 */
inline double score(const double signal_T, const size_t kmer_N, const size_t kmer_K, const vector<tuple<double, double>> &model) {
    // Access elements of the model tuple directly to avoid redundant tuple creation and overhead
    const auto &[mean_N, stdev_N] = model[kmer_N];
    const auto &[mean_K, stdev_K] = model[kmer_K];

    // Precompute the scores for the individual kmers and their distance
    double scoreNT = log_normal_pdf(signal_T, mean_N, stdev_N);
    double scoreKT = log_normal_pdf(signal_T, mean_K, stdev_K);
    double scoreNK = -distanceSequenceKmer(kmer_N, kmer_K); // log(exp(-HD)) = -HD

    return scoreNT + scoreKT + scoreNK;
}

/**
 * Calculate the logarithmic probability matrix
 *
 * @return matrix containing logarithmic probabilities for segment borders
 */
void logP(unordered_map<size_t, array<dproxy, NUMMAT>> &logAPSEI, const unordered_map<size_t, array<dproxy, NUMMAT>> &forAPSEI, const unordered_map<size_t, array<dproxy, NUMMAT>> &backAPSEI, const double Z, const vector<size_t> &allowedKeys, const size_t N, const size_t K) {
    for(size_t tnk : allowedKeys){
        // Lookup in unordered_map is done once per index to avoid repetitive lookups
        auto &logAPSEIRef = logAPSEI[tnk]; // index must not exist
        const auto &forAPSEIRef = forAPSEI.at(tnk); // .at(idx) : index must exist
        const auto &backAPSEIRef = backAPSEI.at(tnk); // .at(idx) : index must exist
        for(int mat=0; mat<NUMMAT; ++mat) {
            logAPSEIRef[mat] = forAPSEIRef[mat] + backAPSEIRef[mat] - Z;
        }
    }
}

void logP(double* LP, const double* forM, const double* backM, const double* forE, const double* backE, const size_t N, const double Z){
    for(size_t t=0; t<T; ++t){
        const size_t tN = t*N;
        for(size_t n=0; n<N; ++n){
            const size_t idx = tN+n;
            const double valM = forM[idx] + backM[idx] - Z;
            const double valE = forE[idx] + backE[idx] - Z;
            LP[idx] = logPlus(valM, valE);
        }
    }
}

// ===============================================================
// ===============================================================
// ======================== PREPROCESSING ========================
// ===============================================================
// ===============================================================

/**
 * Calculate forward matrices using logarithmic values
 *
 * @param sig ONT raw signal with pA values
 * @param seq nucleotide sequence represented by the ONT signal
 * @param T length of the ONT raw signal + 1
 * @param N length of nucleotide sequence + 1
 * @param model map containing kmers as keys and (mean, stdev) tuples as values
 */
void ppForTN(const double* sig, const int* kmer_seq, double* M, double* E, const size_t T, const size_t N, const vector<tuple<double, double>> &model){
    constexpr double m = 0.032, e1 = 1.0, e2 = 0.968;
    for(size_t t=0; t<T; ++t){
        const size_t tN = t*N;
        for(size_t n=0; n<N; ++n){
            double mat=-INFINITY;
            double ext=-INFINITY;
            if(t==0 && n==0) [[unlikely]] {
                ext = 0;
            } else if (t>0 && n>0) [[likely]] {
                const double score = scoreKmer(sig[t-1], kmer_seq[n-1], model);
                mat=logPlus(mat, E[tN-N+(n-1)] + score + m);
                ext=logPlus(ext, M[tN-N+n] + score + e1); // e1 first extend
                ext=logPlus(ext, E[tN-N+n] + score + e2); // e2 extend further
            }
            M[tN+n] = mat;
            E[tN+n] = ext;
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
void ppBackTN(const double* sig, const int* kmer_seq, double* M, double* E, const size_t T, const size_t N, const vector<tuple<double, double>> &model) {
    constexpr double m = 0.032, e1 = 1.0, e2 = 0.968;
    for(size_t t=T; t-->0;){ // iterates from T-1 to 0, first evaluates t>0? then calculates t=t-1
        const size_t tN = t*N;
        for(size_t n=N; n-->0;){
            double mat=-INFINITY;
            double ext=-INFINITY;
            if(t==T-1 && n==N-1) [[unlikely]] {
                ext = 0;
            }
            if (t+1<T) [[likely]] {
                if (n+1<N) [[likely]] {
                    ext=logPlus(ext, M[tN+N+n+1] + scoreKmer(sig[t], kmer_seq[n], model) + m);
                }
                if (n>0) [[likely]] {
                    const double score = scoreKmer(sig[t], kmer_seq[n-1], model);
                    mat=logPlus(mat, E[tN+N+n] + score + e1); // e1 first extend
                    ext=logPlus(ext, E[tN+N+n] + score + e2); // e2 extend further
                }
            }
            M[tN+n] = mat;
            E[tN+n] = ext;
        }
    }
}

/**
 * Calculate forward matrices using logarithmic values
 *
 * @param sig ONT raw signal with pA values
 * @param T length of the ONT raw signal + 1
 * @param K number of allowed kmers
 * @param model map containing kmers as keys and (mean, stdev) tuples as values
 */
void ppForTK(const double* sig, double* M, double* E, const size_t T, const size_t K, const vector<tuple<double, double>> &model){
    constexpr double m = 0.03189915859979101, e1 = 1.0, e2 = 0.9681008434763126;
    for(size_t t=0; t<T; ++t){
        const size_t tK = t*K;
        const size_t prevTK = tK-K; // (t-1)*K
        for(size_t k=0; k<K; ++k){
            double mat=-INFINITY;
            double ext=-INFINITY;
            if (t==0) [[unlikely]] {
                ext = 0.0;
            // t>0
            } else [[likely]] {
                const double score = scoreKmer(sig[t-1], k, model);
                for(size_t preKmer=precessingKmer(k, 0, stepSize); preKmer<K; preKmer+=stepSize){
                    // mat=logPlus(mat, E[prevTK+preKmer] + scoreKmer(sig[t-1], preKmer, model) + m);
                    mat=logPlus(mat, E[prevTK+preKmer] + score + m);
                }
                ext=logPlus(ext, M[prevTK+k] + score + e1); // e1 first extend
                ext=logPlus(ext, E[prevTK+k] + score + e2); // e2 extend further
            }
            M[tK+k] = mat;
            E[tK+k] = ext;
        }
    }
}


/**
 * Calculate backward matrices using logarithmic values
 *
 * @param sig ONT raw signal with pA values
 * @param T length of the ONT raw signal + 1
 * @param K number of allowed kmers
 * @param model map containing kmers as keys and (mean, stdev) tuples as values
 */
void ppBackTK(const double* sig, double* M, double* E, const size_t T, const size_t K, const vector<tuple<double, double>> &model) {
    constexpr double m = 0.03189915859979101, e1 = 1.0, e2 = 0.9681008434763126;
    for(size_t t=T; t-->0;){ // T-1, ..., 0
        const size_t tK = t*K;
        const size_t nexttK = tK+K; // (t+1)*K
        for(size_t k=K; k-->0;){ // K-1, ..., 0
            double mat=-INFINITY;
            double ext=-INFINITY;
            if (t==T-1) [[unlikely]] {
                ext = 0.0;
            // t<T-1
            } else [[likely]] {
                const size_t startKmer = successingKmer(k, 0, stepSize);
                const size_t endKmer = startKmer + ALPHABET_SIZE;
                for(size_t sucKmer = startKmer; sucKmer<endKmer; ++sucKmer){
                    ext=logPlus(ext, M[nexttK+sucKmer] + scoreKmer(sig[t], sucKmer, model) + m);
                }
                const double score = scoreKmer(sig[t], k, model);
                mat=logPlus(mat, E[nexttK+k] + score + e1); // e1 first extend
                ext=logPlus(ext, E[nexttK+k] + score + e2); // e2 extend further

            }
            M[tK+k] = mat;
            E[tK+k] = ext;
        }
    }
}

/**
 * Collect TN key pairs with highest probability
 */
void preProcTN(const double *sig, const int *kmer_seq, unordered_map<size_t, unordered_set<size_t>> &tnMap, const size_t T, const size_t N, const vector<tuple<double, double>> &model) {
    // Allocate memory in one go to reduce overhead
    const size_t TN = T*N;
    std::vector<double> forM(TN, -INFINITY);
    std::vector<double> forE(TN, -INFINITY);
    std::vector<double> backM(TN, -INFINITY);
    std::vector<double> backE(TN, -INFINITY);
    std::vector<double> LP(TN, -INFINITY);

    // Calculate forward and backward matrices    
    ppForTN(sig, kmer_seq, forM.data(), forE.data(), T, N, model);
    ppBackTN(sig, kmer_seq, backM.data(), backE.data(), T, N, model);
    
    const double Zf = forE[TN-1];
    const double Zb = backE[0];

    // Numeric error is scaled by input size, Z in forward and backward should match by some numeric error EPSILON
    if (abs(Zf-Zb)/(TN) > EPSILON) {
        cerr<<"Z values of preProcTN matrices do not match! Zf: "<<Zf<<", Zb: "<<Zb<<", "<<abs(Zf-Zb)/(TN)<<" > "<<EPSILON<<endl;
        cerr.flush();
        exit(11);
    }
    
    // cerr<<"preProcTN: Zf: "<<Zf<<", Zb: "<<Zb<<", "<<abs(Zf-Zb)/(TN)<<" <! "<<EPSILON<<endl;

    logP(LP.data(), forM.data(), backM.data(), forE.data(), backE.data(), N, Zf);
    constexpr double SPARSE_THRESHOLD = -0.04575749056; // log(0.9)

    // extract indices with highest probability per column, until SPARSE_THRESHOLD is reached
    size_t c = 0;
    for(size_t t=0; t<T; ++t){
        double s = -INFINITY;
        // get indices of values in descending order
        for(size_t n : column_argsort(LP.data(), N, t)) {
            // collect key pairs with highest value per column t
            // allowedKeys.push_back(t*N+n);
            ++c;
            tnMap[t].insert(n);
            s = logPlus(s, LP.at(t*N+n));
            // stop if threshold is reached
            if (s > SPARSE_THRESHOLD) {
                break;
            }
        }
    }
    cerr<<"TN dense: "<<c/double(TN)<<endl;
}

/**
 * Collect TK key pairs with highest probability
 */
void preProcTK(const double *sig, unordered_map<size_t, unordered_set<size_t>> &tkMap, const size_t T, const size_t K, const vector<tuple<double, double>> &model) {
    // Allocate memory in one go to reduce overhead
    const size_t TK = T * K;
    std::vector<double> forM(TK, -INFINITY);
    std::vector<double> forE(TK, -INFINITY);
    std::vector<double> backM(TK, -INFINITY);
    std::vector<double> backE(TK, -INFINITY);
    std::vector<double> LP(TK, -INFINITY);

    ppForTK(sig, forM.data(), forE.data(), T, K, model);
    ppBackTK(sig, backM.data(), backE.data(), T, K, model);

    double Zf = -INFINITY;
    double Zb = -INFINITY;
    for(size_t k=0; k<K; ++k){
        Zf = logPlus(Zf, forE[TK-1-k]);
        Zb = logPlus(Zb, backE[k]);
    }

    // Numeric error is scaled by input size, Z in forward and backward should match by some numeric error EPSILON
    if (abs(Zf-Zb)/(TK) > EPSILON) {
        cerr<<"Z values of preProcTK matrices do not match! Zf: "<<Zf<<", Zb: "<<Zb<<", "<<abs(Zf-Zb)/(TK)<<" > "<<EPSILON<<endl;
        cerr.flush();
        exit(12);
    }

    // cerr<<"preProcTK: Zf: "<<Zf<<", Zb: "<<Zb<<", "<<abs(Zf-Zb)/(TK)<<" <! "<<EPSILON<<endl;

    logP(LP.data(), forM.data(), backM.data(), forE.data(), backE.data(), K, Zb);
    constexpr double SPARSE_THRESHOLD = -0.04575749056; // log(0.9)

    // extract indices with highest probability per column, until SPARSE_THRESHOLD is reached
    size_t c = 0;
    for(size_t t=0; t<T; ++t){
        double s = -INFINITY;
        // get indices of values in descending order
        for(size_t k : column_argsort(LP.data(), K, t)) {
            // collect key pairs with highest value per column t
            // allowedKeys.push_back(t*K+k);
            ++c;
            tkMap[t].insert(k);
            s = logPlus(s, LP.at(t*K+k));
            // stop if threshold is reached
            if (s >= SPARSE_THRESHOLD) {
                break;
            }
        }
    }
    cerr<<"TK dense: "<<c/double(TK)<<endl;
}

void preProcTNK(const double *sig, const int *kmer_seq, vector<size_t> &allowedKeys, const size_t T, const size_t N, const size_t K, const vector<tuple<double, double>> &model) {
    // perform preprocessing on partial 2d problems
    unordered_map<size_t, unordered_set<size_t>> tnMap;
    preProcTN(sig, kmer_seq, tnMap, T, N, model);
    unordered_map<size_t, unordered_set<size_t>> tkMap;
    preProcTK(sig, tkMap, T, K, model);

    for (size_t t=0; t<T; t++) {
        for(size_t n : tnMap[t]) {
            for(size_t k: tkMap[t]) {
                allowedKeys.push_back(t * NK + n * K + k);
            }
        }
    }

    // erase duplicates
    sort(allowedKeys.begin(), allowedKeys.end());
}

// ===============================================================
// ========================== ALGORITHM ==========================
// ===============================================================

/**
 * Calculate forward matrices using logarithmic values
 *
 * @param sig ONT raw signal with pA values
 * @param seq nucleotide sequence represented by the ONT signal
 * @param T length of the ONT raw signal + 1
 * @param N length of nucleotide sequence + 1
 * @param model map containing kmers as keys and (mean, stdev) tuples as values
 */
void logF(const double *sig, const int *kmer_seq, unordered_map<size_t, array<dproxy, NUMMAT>> &forAPSEI, const vector<size_t> &allowedKeys, const size_t T, const size_t N, const size_t K, const vector<tuple<double, double>> &model){
    double a, p, s, e, i;
    array<dproxy, NUMMAT> forAPSEIRef;
    size_t t, n, k;
    for(size_t tnk : allowedKeys){
        // tnk = t*NK+n*K+k
        t = tnk/NK;
        n = (tnk % NK) / K;
        k = tnk % K;
        a=-INFINITY;
        p=-INFINITY;
        s=-INFINITY;
        e=-INFINITY;
        i=-INFINITY;
        if(t==0 && n==0) [[unlikely]] {
            e=0;
        } else if(t>0 && n>0) [[likely]] {
            // Precompute expensive values
            const double extendScore = score(sig[t-1], kmer_seq[n-1], k, model);
            const double openScore = score(sig[t-1], kmer_seq[n-1], k, 0.05, model);
            const size_t baseIdx1 = tnk - NK - K - k;
            const size_t baseIdx2 = tnk - NK - k;
            const size_t baseIdx3 = tnk - NK - K;
            const size_t baseIdx4 = tnk - NK;
            const size_t baseIdx5 = tnk - K;

            // non-consecutive, differs by 5^4
            for(size_t preKmer = precessingKmer(k, 0, stepSize); preKmer<K; preKmer+=stepSize) {
                // (t-1)*NK+(n-1)*K+k'
                forAPSEIRef = forAPSEI[baseIdx1+preKmer];
                a=logPlus(a, forAPSEIRef[3] + a1 + openScore);
                a=logPlus(a, forAPSEIRef[4] + a2 + openScore);
                
                // (t-1)*NK+n*K+k'
                forAPSEIRef = forAPSEI[baseIdx2+preKmer];
                p=logPlus(p, forAPSEIRef[2] + p1 + openScore);
                p=logPlus(p, forAPSEIRef[3] + p2 + openScore);
                p=logPlus(p, forAPSEIRef[4] + p3 + openScore);
            }

            // (t-1)*NK+(n-1)*K+k
            forAPSEIRef = forAPSEI[baseIdx3];
            s=logPlus(s, forAPSEIRef[1] + s1 + openScore);
            s=logPlus(s, forAPSEIRef[3] + s2 + openScore);
            s=logPlus(s, forAPSEIRef[4] + s3 + openScore);

            // (t-1)*NK+n*K+k
            forAPSEIRef = forAPSEI[baseIdx4];
            e=logPlus(e, forAPSEIRef[0] + e1 + extendScore);
            e=logPlus(e, forAPSEIRef[1] + e2 + extendScore);
            e=logPlus(e, forAPSEIRef[2] + e3 + extendScore);
            e=logPlus(e, forAPSEIRef[3] + e4 + extendScore);

            // t*NK+(n-1)*K+k
            forAPSEIRef = forAPSEI[baseIdx5];
            i=logPlus(i, forAPSEIRef[3] + i1 + openScore);
            i=logPlus(i, forAPSEIRef[4] + i2 + openScore);
        }
        forAPSEI[tnk] = {a, p, s, e, i};
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
void logB(const double *sig, const int *kmer_seq, unordered_map<size_t, array<dproxy, NUMMAT>> &backAPSEI, vector<size_t> &allowedKeys, const size_t T, const size_t N, const size_t K, const vector<tuple<double, double>> &model){
    double a, p, s, e, i, pv, sc;
    size_t t, n, k;
    for(auto tnk = allowedKeys.rbegin(); tnk != allowedKeys.rend(); ++tnk){
        t = *tnk/NK;
        n = (*tnk % NK) / K;
        k = *tnk % K;
        a=-INFINITY;
        p=-INFINITY;
        s=-INFINITY;
        e=-INFINITY;
        i=-INFINITY;
        if(t==T-1 && n==N-1) [[unlikely]] {
            e=0;
        }
        if (t<T-1) [[likely]] { 
            // (t+1)*NK+nK+k
            const size_t next_tnk = *tnk + NK;
            // (t+1)*NK+(n+1)*K+k
            const size_t next_tnk_n = next_tnk + K;
            // Cache results of successingKmer to avoid recomputation
            const size_t sucKmerBase = successingKmer(k, 0, stepSize);
            const size_t sucKmerEnd = sucKmerBase + ALPHABET_SIZE;

            if (n>0) [[likely]] {
                // Precompute score for efficiency
                sc=score(sig[t], kmer_seq[n-1], k, model);
                
                // (t+1)*NK+n*K+k
                pv = backAPSEI[next_tnk][3];
                a = logPlus(a, pv + e1 + sc);
                p = logPlus(p, pv + e2 + sc);
                s = logPlus(s, pv + e3 + sc);
                e = logPlus(e, pv + e4 + sc);
                
                // kmer int representation is consecutive
                for(size_t sucKmer = sucKmerBase; sucKmer<sucKmerEnd; ++sucKmer){
                    sc=score(sig[t], kmer_seq[n-1], sucKmer, 0.05, model);
                    // (t+1)*NK+n*K+sucKmer
                    pv=backAPSEI[next_tnk-k+sucKmer][1];
                    s=logPlus(s, pv + p1 + sc);
                    e=logPlus(e, pv + p2 + sc);
                    i=logPlus(i, pv + p3 + sc);
                }
            }

            if (n<N-1) [[likely]] {
                sc=score(sig[t], kmer_seq[n], k, 0.05, model);
                // (t+1)*(NK)+(n+1)*K+k
                pv=backAPSEI[next_tnk_n][2];
                p=logPlus(p, pv + s1 + sc);
                e=logPlus(e, pv + s2 + sc);
                i=logPlus(i, pv + s3 + sc);

                // kmer int representation is consecutive
                for(size_t sucKmer=sucKmerBase; sucKmer<sucKmerEnd; ++sucKmer){
                    sc=score(sig[t], kmer_seq[n], sucKmer, 0.05, model);
                    // (t+1)*NK+(n+1)*K+sucKmer
                    pv=backAPSEI[next_tnk_n-k+sucKmer][0];
                    e=logPlus(e, pv + a1 + sc);
                    i=logPlus(i, pv + a2 + sc);
                }
            }
        }

        if (t>0 && n<N-1) [[likely]] {
            sc=score(sig[t-1], kmer_seq[n], k, 0.05, model);
            // t*NK+(n+1)*K+k
            pv=backAPSEI[*tnk+K][4];
            e=logPlus(e, pv + i1 + sc);
            i=logPlus(i, pv + i2 + sc);
        }

        backAPSEI[*tnk] = {a, p, s, e, i};
    }
}

/**
 * Calculate the maximum a posteriori path through LP
 *
 */
list<string> getBorders(unordered_map<size_t, array<dproxy, NUMMAT>> &logAPSEI, const vector<size_t> &allowedKeys, const size_t T, const size_t N, const size_t K){
    unordered_map<size_t, array<dproxy, NUMMAT>> APSEI;
    array<dproxy, NUMMAT> logAPSEIRef;
    double a, p, s, e, i;
    size_t t, n, k, idx;
    for(size_t tnk : allowedKeys){
        // tnk = t*NK+n*K+k
        t = tnk/NK;
        n = (tnk % NK) / K;
        k = tnk % K;
        a=-INFINITY;
        p=-INFINITY;
        s=-INFINITY;
        e=-INFINITY;
        i=-INFINITY;
        if(t==0 && n==0) [[unlikely]] {
            e=0;
        } 
        if(t>0 && n>0) [[likely]] {
            // Precompute expensive values
            const size_t baseIdx1 = tnk - NK - K - k;   // (t-1)*NK+(n-1)*K+k'
            const size_t baseIdx2 = tnk - NK - k;       // (t-1)*NK+ n   *K+k'
            const size_t baseIdx3 = tnk - NK - K;       // (t-1)*NK+(n-1)*K+k
            const size_t baseIdx4 = tnk - NK;           // (t-1)*NK+ n   *K+k
            const size_t baseIdx5 = tnk - K;            //  t   *NK+(n-1)*K+k
            const array<dproxy, NUMMAT> logAPSEIRef = logAPSEI.at(tnk); // .at(idx) : index must exist
            for(size_t preKmer = precessingKmer(k, 0, stepSize); preKmer<K; preKmer+=stepSize) {
                // (t-1)*NK+(n-1)*K+k'
                idx = baseIdx1+preKmer;
                a=max(a, APSEI[idx][3] + logAPSEIRef[0]);
                a=max(a, APSEI[idx][4] + logAPSEIRef[0]);
                
                // (t-1)*NK+n*K+k'
                idx = baseIdx2+preKmer;
                p=max(p, APSEI[idx][2] + logAPSEIRef[1]);
                p=max(p, APSEI[idx][3] + logAPSEIRef[1]);
                p=max(p, APSEI[idx][4] + logAPSEIRef[1]);
            }
            // (t-1)*NK+(n-1)*K+k
            s=max(s, APSEI[baseIdx3][1] + logAPSEIRef[2]);
            s=max(s, APSEI[baseIdx3][3] + logAPSEIRef[2]);
            s=max(s, APSEI[baseIdx3][4] + logAPSEIRef[2]);

            // (t-1)*NK+n*K+k
            e=max(e, APSEI[baseIdx4][0] + logAPSEIRef[3]);
            e=max(e, APSEI[baseIdx4][1] + logAPSEIRef[3]);
            e=max(e, APSEI[baseIdx4][2] + logAPSEIRef[3]);
            e=max(e, APSEI[baseIdx4][3] + logAPSEIRef[3]);

            // t*NK+(n-1)*K+k
            i=max(i, APSEI[baseIdx5][3] + logAPSEIRef[4]);
            i=max(i, APSEI[baseIdx5][4] + logAPSEIRef[4]);
        }
        APSEI[tnk] = {a, p, s, e, i};
    }

    // get highest k for last t and n?
    double mv = -INFINITY;
    size_t hk = 3125, ld = TNK-K;
    for(size_t k=0; k<K; ++k){
        if (APSEI[ld+k][3] >= mv) {
            mv = APSEI[ld+k][3];
            hk = k;
        }
    }

    list<string> segString;
    funcE(T-1, N-1, hk, APSEI, logAPSEI, &segString, K);
    return segString;
}

void funcA(const size_t t, const size_t n, const size_t k, unordered_map<size_t, array<dproxy, NUMMAT>> &APSEI, unordered_map<size_t, array<dproxy, NUMMAT>> &logAPSEI, list<string>* segString, const size_t K){
    if (t<=1 && n<=1){ // Start value
        segString->push_front("M"+to_string(0)+","+to_string(0)+","+itoa(k, kmerSize)+";"); // n-1 because N is 1 larger than the sequences
        return;
    }
    const size_t currentIdx = t*NK+n*K+k;
    const size_t prevIdx  = (t-1)*NK+(n-1)*K;
    // Cache the score value to avoid redundant lookups
    const double score = APSEI[currentIdx][0];
    const double logScore = logAPSEI[currentIdx][0];
    for(size_t preKmer=precessingKmer(k, 0, stepSize); preKmer<K; preKmer+=stepSize) {
        if (t>1 && n>1) {
            // Check match with E state
            if (score == APSEI[prevIdx+preKmer][3] + logScore){
                segString->push_front("M"+to_string(n-1)+","+to_string(t-1)+","+itoa(k, kmerSize)+";");
                return funcE(t-1, n-1, preKmer, APSEI, logAPSEI, segString, K);
            }
            // Check match with I state
            if (score == APSEI[prevIdx+preKmer][4] + logScore){
                segString->push_front("M"+to_string(n-1)+","+to_string(t-1)+","+itoa(k, kmerSize)+";");
                return funcI(t-1, n-1, preKmer, APSEI, logAPSEI, segString, K);
            }
        }
    }
    // If no match is found, output an error
    cerr<<"Error in backtracing funcA!: t: "<<t<<", n: "<<n<<", k: "<<k<<endl;
}

void funcE(const size_t t, const size_t n, const size_t k, unordered_map<size_t, array<dproxy, NUMMAT>> &APSEI, unordered_map<size_t, array<dproxy, NUMMAT>> &logAPSEI, list<string>* segString, const size_t K){
    const size_t currentIdx = NK*t+n*K+k;
    const size_t prevIdx  = NK*(t-1)+n*K+k;
    // Cache the score value to avoid redundant lookups
    const double score = APSEI[currentIdx][3];
    const double logScore = logAPSEI[currentIdx][3];
    if (t>0 && n>0) {
        // Check match with A state
        if (score == APSEI[prevIdx][0] + logScore){
            return funcA(t-1, n, k, APSEI, logAPSEI, segString, K);
        }
        // Check match with E state
        if (score == APSEI[prevIdx][3] + logScore){
            return funcE(t-1, n, k, APSEI, logAPSEI, segString, K);
        }
        // Check match with S state
        if (score == APSEI[prevIdx][2] + logScore){
            return funcS(t-1, n, k, APSEI, logAPSEI, segString, K);
        }
        // Check match with P state
        if (score == APSEI[prevIdx][1] + logScore){
            return funcP(t-1, n, k, APSEI, logAPSEI, segString, K);
        }
    }
    else [[unlikely]] { // Start value with t==0 && n==0
        return;
    }
    // If no match is found, output an error
    cerr<<"Error in backtracing funcE!: t: "<<t<<", n: "<<n<<", k: "<<k<<endl;
}

void funcP(const size_t t, const size_t n, const size_t k, unordered_map<size_t, array<dproxy, NUMMAT>> &APSEI, unordered_map<size_t, array<dproxy, NUMMAT>> &logAPSEI, list<string>* segString, const size_t K) {
    // Precompute commonly used indices and values
    const size_t currentIdx = t*NK+n*K+k;
    const double score = APSEI[currentIdx][1];
    const double logScore = logAPSEI[currentIdx][1];
    const size_t prevBaseIdx = NK*(t-1)+n*K;  // Common base index for previous time step
    if (t>0 && n>0) {
        for (size_t preKmer = precessingKmer(k, 0, stepSize); preKmer<K; preKmer+=stepSize) {
            const size_t prevIdx = prevBaseIdx + preKmer;
            // Check match with E state
            if (score == APSEI[prevIdx][3] + logScore) {
                segString->push_front("P" + to_string(n) + "," + to_string(t-1) + "," + itoa(k, kmerSize) + ";");
                return funcE(t-1, n, preKmer, APSEI, logAPSEI, segString, K);
            }
            // Check match with S state
            if (score == APSEI[prevIdx][2] + logScore) {
                segString->push_front("P" + to_string(n) + "," + to_string(t-1) + "," + itoa(k, kmerSize) + ";");
                return funcS(t-1, n, preKmer, APSEI, logAPSEI, segString, K);
            }
            // Check match with I state
            if (score == APSEI[prevIdx][4] + logScore) {
                segString->push_front("P" + to_string(n) + "," + to_string(t-1) + "," + itoa(k, kmerSize) + ";");
                return funcI(t-1, n, preKmer, APSEI, logAPSEI, segString, K);
            }
        }
    }
    // If no match is found, output an error
    cerr << "Error in backtracing funcP!: t: "<<t<<", n: "<<n<<", k: "<<k<<endl;
}

void funcS(const size_t t, const size_t n, const size_t k, unordered_map<size_t, array<dproxy, NUMMAT>> &APSEI, unordered_map<size_t, array<dproxy, NUMMAT>> &logAPSEI, list<string>* segString, const size_t K) {
    const size_t currentIdx = NK*t+n*K+k;
    const size_t prevIdx = NK*(t-1)+(n-1)*K+k;
    // Cache score and logScore to avoid repeated map lookups
    const double score = APSEI[currentIdx][2];
    const double logScore = logAPSEI[currentIdx][2];
    if (t>0 && n>0) {
        // Check match with E state
        if (score == APSEI[prevIdx][3] + logScore) {
            segString->push_front("S" + to_string(n-1) + "," + to_string(t-1) + "," + itoa(k, kmerSize) + ";");
            return funcE(t-1, n-1, k, APSEI, logAPSEI, segString, K);
        }
        // Check match with P state
        if (score == APSEI[prevIdx][1] + logScore) {
            segString->push_front("S" + to_string(n-1) + "," + to_string(t-1) + "," + itoa(k, kmerSize) + ";");
            return funcP(t-1, n-1, k, APSEI, logAPSEI, segString, K);
        }
        // Check match with I state
        if (score == APSEI[prevIdx][4] + logScore) {
            segString->push_front("S" + to_string(n-1) + "," + to_string(t-1) + "," + itoa(k, kmerSize) + ";");
            return funcI(t-1, n-1, k, APSEI, logAPSEI, segString, K);
        }
    }
    // If no match is found, output an error
    cerr << "Error in backtracing funcS!: t: "<<t<<", n: "<<n<<", k: "<<k<<endl;
}

void funcI(const size_t t, const size_t n, const size_t k, unordered_map<size_t, array<dproxy, NUMMAT>> &APSEI, unordered_map<size_t, array<dproxy, NUMMAT>> &logAPSEI, list<string>* segString, const size_t K) {
    const size_t currentIdx = NK*t+n*K+k;
    const size_t prevIdx = NK*t+(n-1)*K+k;
    // Cache the score and logScore to avoid repeated lookups
    const double score = APSEI[currentIdx][4];
    const double logScore = logAPSEI[currentIdx][4];
    if (t>0 && n>0) {
        // Check match with I state
        if (score == APSEI[prevIdx][4] + logScore) {
            segString->push_front("I" + to_string(n-1) + "," + to_string(t) + "," + itoa(k, kmerSize) + ";");
            return funcI(t, n-1, k, APSEI, logAPSEI, segString, K);
        } 
        // Check match with E state
        if (score == APSEI[prevIdx][3] + logScore) {
            segString->push_front("I" + to_string(n-1) + "," + to_string(t) + "," + itoa(k, kmerSize) + ";");
            return funcE(t, n-1, k, APSEI, logAPSEI, segString, K);
        }
    }
    // If no match is found, output an error
    cerr << "Error in backtracing funcI!: t: "<<t<<", n: "<<n<<", k: "<<k<<endl;
}

/**
 * Train transition parameter with baum welch algorithm
*/
tuple<double, double, double, double, double, double, double, double, double, double, double, double, double, double> trainTransition(const double *sig, const int *kmer_seq, unordered_map<size_t, array<dproxy, NUMMAT>> &forAPSEI, unordered_map<size_t, array<dproxy, NUMMAT>> &backAPSEI, vector<size_t> &allowedKeys, const size_t T, const size_t N, const size_t K, const vector<tuple<double, double>> &model) {
    // Transition parameters
    double newa1=-INFINITY, newa2=-INFINITY;
    double newp1=-INFINITY, newp2=-INFINITY, newp3=-INFINITY;
    double news1=-INFINITY, news2=-INFINITY, news3=-INFINITY;
    double newe1=1, newe2=-INFINITY, newe3=-INFINITY, newe4=-INFINITY;
    double newi1=-INFINITY, newi2=-INFINITY;
    double sc, pv;
    size_t t, n, k;

    for(size_t tnk : allowedKeys){ //<int>::iterator iter = allowedKeys.begin(); iter<allowedKeys.end(); iter++){
        // tnk = t*NK+n*K+k
        t = tnk/NK;
        n = (tnk % NK) / K;
        k = tnk % K;
        // Cache results to avoid recomputation
        auto &forAPSEI_tnk = forAPSEI[tnk];

        if (t<T-1) [[likely]] {
            const size_t sucKmerBase = successingKmer(k, 0, stepSize);
            const size_t sucKmerEnd = sucKmerBase + ALPHABET_SIZE;

            if (n>0) [[likely]] {
                sc = score(sig[t], kmer_seq[n-1], k, model);
                // (t+1)*NK+n*K+k
                pv = backAPSEI[tnk+NK][3];
                // newe1 = logPlus(newe1, forAPSEI_tnk[0] + e1 + sc + pv);
                newe2 = logPlus(newe2, forAPSEI_tnk[1] + e2 + sc + pv);
                newe3 = logPlus(newe3, forAPSEI_tnk[2] + e3 + sc + pv);
                newe4 = logPlus(newe4, forAPSEI_tnk[3] + e4 + sc + pv);

                for(size_t sucKmer=sucKmerBase; sucKmer<sucKmerEnd; ++sucKmer){
                    sc = score(sig[t], kmer_seq[n-1], sucKmer, 0.05, model);
                    // (t+1)*NK+n*K+sucKmer
                    pv = backAPSEI[tnk+NK-k+sucKmer][1];
                    newp1 = logPlus(newp1, forAPSEI_tnk[2] + p1 + sc + pv);
                    newp2 = logPlus(newp2, forAPSEI_tnk[3] + p2 + sc + pv);
                    newp3 = logPlus(newp3, forAPSEI_tnk[4] + p3 + sc + pv);
                }
            }

            if (n<N-1) [[likely]] {
                sc = score(sig[t], kmer_seq[n], k, 0.05, model);
                // (t+1)*(NK)+(n+1)*K+k
                pv = backAPSEI[tnk+NK+K][2];
                news1 = logPlus(news1, forAPSEI_tnk[1] + s1 + sc + pv);
                news2 = logPlus(news2, forAPSEI_tnk[3] + s2 + sc + pv);
                news3 = logPlus(news3, forAPSEI_tnk[4] + s3 + sc + pv);
                
                for(size_t sucKmer=sucKmerBase; sucKmer<sucKmerEnd; ++sucKmer){
                    sc = score(sig[t], kmer_seq[n], sucKmer, 0.05, model);
                    // (t+1)*NK+(n+1)*K+sucKmer
                    pv = backAPSEI[tnk+NK+K-k+sucKmer][0];
                    newa1 = logPlus(newa1, forAPSEI_tnk[3] + a1 + sc + pv);
                    newa2 = logPlus(newa2, forAPSEI_tnk[4] + a2 + sc + pv);
                }
            }
        }

        if (t>0 && n<N-1) [[likely]] {
            sc = score(sig[t-1], kmer_seq[n], k, 0.05, model);
            // t*NK+(n+1)*K+k
            pv = backAPSEI[tnk+K][4];
            newi1 = logPlus(newi1, forAPSEI_tnk[3] + i1 + sc + pv);
            newi2 = logPlus(newi2, forAPSEI_tnk[4] + i2 + sc + pv);
        }
    }

    // Final normalization and averaging of transition parameters
    // newe1=exp(newe1-newe1);
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
    double Ai = logPlus(logPlus(newa2, newi2), logPlus(newp3, news3));
    newa2=exp(newa2-Ai);
    newi2=exp(newi2-Ai);
    newp3=exp(newp3-Ai);
    news3=exp(news3-Ai);

    return tuple<double, double, double, double, double, double, double, double, double, double, double, double, double, double>({newa1, newa2, newp1, newp2, newp3, news1, news2, news3, newe1, newe2, newe3, newe4, newi1, newi2});
}

/**
 * Train emission parameter with baum welch algorithm
*/
tuple<double*, double*> trainEmission(const double* sig, unordered_map<size_t, array<dproxy, NUMMAT>> &logAPSEI, vector<size_t> &allowedKeys, const size_t T, const size_t N, const size_t K) {
    // Emission
    // https://courses.grainger.illinois.edu/ece417/fa2021/lectures/lec15.pdf
    // https://f.hubspotusercontent40.net/hubfs/8111846/Unicon_October2020/pdf/bilmes-em-algorithm.pdf
    size_t k, t;
    double* means = new double[K];
    double* stdevs = new double[K];
    double* normFactorT = new double[K];
    double w;
    constexpr double TRAIN_THRESHOLD = 1e-9; // arbitrarily set

    // First pass: Compute means and normalization factors
    for(size_t tnk : allowedKeys){
        // Skip if t == 0 (equivalent to tnk < NK)
        if (tnk<NK) [[unlikely]] {
            continue;
        }
        // tnk = t*NK+n*K+k
        t = tnk/NK;
        k = tnk % K;
        // Cache logAPSEI access
        const auto& logValues = logAPSEI.at(tnk);
        w = logPlus(logPlus(logPlus(logValues[0], logValues[1]), logPlus(logValues[2], logValues[3])), logValues[4]);
        w = exp(w); // Convert log probability to normal probability
        means[k] += w * sig[t-1];
        normFactorT[k] += w;
    }

    // Normalize the means
    for(size_t k=0; k<K; ++k){
        means[k] = means[k] / normFactorT[k];
    }

    // Emission (stdev of kmers)
    for(size_t tnk : allowedKeys){
        // tnk = t*NK+n*K+k
        k = tnk % K;
        // Skip kmers with low weight or if t == 0
        if(normFactorT[k]/T < TRAIN_THRESHOLD || tnk<NK) [[likely]] {
            continue;
        }
        t = tnk/NK;
        // Cache logAPSEI access
        const auto& logValues = logAPSEI.at(tnk);
        w = logPlus(logPlus(logPlus(logValues[0], logValues[1]), logPlus(logValues[2], logValues[3])), logValues[4]);
        // Update standard deviation
        double diff = sig[t - 1] - means[k];
        stdevs[k] += exp(w) * diff * diff;
    }

    // Normalize the standard deviations
    for(size_t k=0; k<K; ++k){
        stdevs[k] = sqrt(stdevs[k] / normFactorT[k]);
    }

    delete[] normFactorT;
    return tuple<double*, double*>({means, stdevs});
}

void trainParams(const double *sig, const int *kmer_seq, unordered_map<size_t, array<dproxy, NUMMAT>> &forAPSEI, unordered_map<size_t, array<dproxy, NUMMAT>> &backAPSEI, unordered_map<size_t, array<dproxy, NUMMAT>> &logAPSEI, vector<size_t> &allowedKeys, const size_t T, const size_t N, const size_t K, vector<tuple<double, double>> &model) {

    // newa1, newa2, newp1, newp2, newp3, news1, news2, news3, newe1, newe2, newe3, newe4, newi1, newi2
    auto [a1, a2, p1, p2, p3, s1, s2, s3, e1, e2, e3, e4, i1, i2] = trainTransition(sig, kmer_seq, forAPSEI, backAPSEI, allowedKeys, T, N, K, model);
    cout<<"a1:"<<a1<<";a2:"<<a2<<";p1:"<<p1<<";p2:"<<p2<<";p3:"<<p3<<";s1:"<<s1<<";s2:"<<s2<<";s3:"<<s3<<";e1:"<<e1<<";e2:"<<e2<<";e3:"<<e3<<";e4:"<<e4<<";i1:"<<i1<<";i2:"<<i2<<endl;

    auto [newMeans, newStdevs] = trainEmission(sig, logAPSEI, allowedKeys, T, N, K);
    for(size_t k=0; k<K; ++k){
        if (newStdevs[k]!=0.0){
            cout<<itoa(k, kmerSize)<<":"<<newMeans[k]<<","<<newStdevs[k]<<";";
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
    int pore;
    bool train, calcZ; // atrain
    string modelpath;
    const string TERM_STRING = "$";

    cerr << fixed << showpoint;
    cerr << setprecision(20);

    // Argparser
    argparse::ArgumentParser program("dynamont 3d sparsed", "0.1");
    // parameters for DP

    // double a1, a2, p1, p2, p3, s1, s2, s3, e1, e2, e3, e4, i1, i2; // transition parameters
    program.add_argument("-e1", "--extendscore1").help("Transition parameter").default_value(1.00).scan<'g', double>().store_into(e1); // e1

    program.add_argument("-e2", "--extendscore2").help("Transition parameter").default_value(0.944).scan<'g', double>().store_into(e2); // e2
    program.add_argument("-s1", "--sequencescore1").help("Transition parameter").default_value(0.056).scan<'g', double>().store_into(s1); // s1

    program.add_argument("-e3", "--extendscore3").help("Transition parameter").default_value(0.086).scan<'g', double>().store_into(e3); // e3
    program.add_argument("-p1", "--polishscore1").help("Transition parameter").default_value(0.914).scan<'g', double>().store_into(p1); // p1

    program.add_argument("-a1", "--alignscore1").help("Transition parameter").default_value(0.007).scan<'g', double>().store_into(a1); // a1
    program.add_argument("-p2", "--polishscore2").help("Transition parameter").default_value(0.758).scan<'g', double>().store_into(p2); // p2
    program.add_argument("-e4", "--extendscore4").help("Transition parameter").default_value(0.231).scan<'g', double>().store_into(e4); // e4
    program.add_argument("-s2", "--sequencescore2").help("Transition parameter").default_value(0.003).scan<'g', double>().store_into(s2); // s2
    program.add_argument("-i1", "--insertionscore1").help("Transition parameter").default_value(0.001).scan<'g', double>().store_into(i1); // i1

    program.add_argument("-a2", "--alignscore2").help("Transition parameter").default_value(0.022).scan<'g', double>().store_into(a2); // a2
    program.add_argument("-p3", "--polishscore3").help("Transition parameter").default_value(0.031).scan<'g', double>().store_into(p3); // p3
    program.add_argument("-s3", "--sequencescore3").help("Transition parameter").default_value(0.654).scan<'g', double>().store_into(s3); // s3
    program.add_argument("-i2", "--insertionscore2").help("Transition parameter").default_value(0.293).scan<'g', double>().store_into(i2); // i2

    program.add_argument("-t", "--train").help("Switch algorithm to transition and emission parameter training mode").default_value(false).implicit_value(true).store_into(train);
    program.add_argument("-z", "--calcZ").help("Switch algorithm to only calculate Z").default_value(false).implicit_value(true).store_into(calcZ);
    program.add_argument("-m", "--model").help("Path to kmer model table").default_value("/home/yi98suv/projects/dynamont/data/norm_models/rna_r9.4_180mv_70bps_extended_stdev0_25.model").store_into(modelpath);
    program.add_argument("-r", "--pore").help("Pore generation used to sequence the data").default_value(9).choices(9, 10).store_into(pore);
    // unused, just here to match the other modes
    program.add_argument("-c", "--minSegLen").help("MinSegLen + 1 is the minimal segment length").default_value(0); //.store_into(C);
    program.add_argument("-p", "--probabilty").help("Print out the segment border probability").default_value(false).implicit_value(true); //.store_into(prob);

    try {
        program.parse_args(argc, argv);
    }
    catch (const runtime_error& err) {
        cerr << err.what() << std::endl;
        cerr << program;
        return 1;
    }
    
    if (pore == 9) {
        kmerSize = 5;
    } else if (pore == 10) {
        kmerSize = 9;
    }

    // polishing dimension K = number of possible kmers
    K = pow(ALPHABET_SIZE, kmerSize); // currently acceptable A, C, G, T, N
    stepSize = pow(ALPHABET_SIZE, kmerSize-1);
    vector<tuple<double, double>> model(K, make_tuple(-INFINITY, -INFINITY));
    readKmerModel(modelpath, model, kmerSize);
    string signal, read;

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
        size_t seq_size = read.size() + (kmerSize-1);
        int* seq = new int[seq_size];
        fill_n(seq, seq_size, 4); // default: fill with N add 2 Ns to 3' of read
        i = floor(kmerSize/2);
        for(const char &c: read) {
            seq[i] = BASE2ID.at(c);
            i++;
        }
        // add NN to end of sequence
        seq[i] = 4;
        seq[i+1] = 4;

        int* kmer_seq = seq2kmer(seq, N-1, kmerSize);
        NK = N*K;
        TNK = T*NK;

        cerr<<"T: "<<T<<", "<<"N: "<<N<<", "<<"K: "<<K<<", "<<"inputsize: "<<TNK<<endl;

        vector<size_t> allowedKeys;
        preProcTNK(sig, kmer_seq, allowedKeys, T, N, K, model);
        cerr<<"dense: "<<allowedKeys.size()/double(TNK)<<" ("<<allowedKeys.size()<<")"<<", sparse: "<<1-(allowedKeys.size()/double(TNK))<<" ("<<TNK-allowedKeys.size()<<")"<<endl;
        unordered_map<size_t, array<dproxy, NUMMAT>> forAPSEI;
        
        cerr<<"forward"<<endl;
        logF(sig, kmer_seq, forAPSEI, allowedKeys, T, N, K, model);
        cerr<<"backward"<<endl;
        unordered_map<size_t, array<dproxy, NUMMAT>> backAPSEI;
        logB(sig, kmer_seq, backAPSEI, allowedKeys, T, N, K, model);

        double Zf = -INFINITY;
        double Zb = -INFINITY;
        for(size_t k=0; k<K; ++k){
            Zf = logPlus(Zf, forAPSEI[TNK-1-k][3]);
            Zb = logPlus(Zb, backAPSEI[k][3]);
        }

        // Numeric error is scaled by input size, Z in forward and backward should match by some numeric error EPSILON
        if (abs(Zf-Zb)/TNK >= EPSILON) {
            cerr<<"Z values between matrices do not match! forZ: "<<Zf<<", backZ: "<<Zb<<", "<<abs(Zf-Zb)/TNK<<" > "<<EPSILON<<endl;
            cerr.flush();
            exit(13);
        }

        // cerr<<"Zf: "<<Zf<<", Zb: "<<Zb<<", "<<abs(Zf-Zb)/TNK<<" <! "<<EPSILON<<endl;

        if (calcZ){
            cout<<Zf<<endl;
            cout.flush();
        } else {
            unordered_map<size_t, array<dproxy, NUMMAT>> logAPSEI;
            logP(logAPSEI, forAPSEI, backAPSEI, Zf, allowedKeys, N, K);

            // train both Transitions and Emissions
            if (train) {
                trainParams(sig, kmer_seq, forAPSEI, backAPSEI, logAPSEI, allowedKeys, T, N, K, model);
                cout<<"Z:"<<Zf<<endl;
                cout.flush();

            // print out segmentation string
            } else {
                list<string> segString = getBorders(logAPSEI, allowedKeys, T, N, K);

                for(auto const& seg : segString) {
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
