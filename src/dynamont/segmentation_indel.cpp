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
#include <assert.h>
#include <stdlib.h>
#include <algorithm>
#include "argparse.hpp"
#include "utils.hpp"

using namespace std;

void funcM(const size_t t, const size_t n, const double* M, const double* E, const double* D, const double* I, const double* LPM, const double* LPE, const double* LPD, const double* LPI, list<string>* segString, const size_t N);
void funcE(const size_t t, const size_t n, const double* M, const double* E, const double* D, const double* I, const double* LPM, const double* LPE, const double* LPD, const double* LPI, list<string>* segString, const size_t N);
void funcD(const size_t t, const size_t n, const double* M, const double* E, const double* D, const double* I, const double* LPM, const double* LPE, const double* LPD, const double* LPI, list<string>* segString, const size_t N);
void funcI(const size_t t, const size_t n, const double* M, const double* E, const double* D, const double* I, const double* LPM, const double* LPE, const double* LPD, const double* LPI, list<string>* segString, const size_t N);

int numKmers, kmerSize; // our model works with this kmer size
inline constexpr double EPSILON = 1e-8; // chose by eye just to distinguish real errors from numeric errors
double m1, m2, m3, e1, e2, d1, d2, i1, i2; // transition parameters
size_t T, N, TN, C;

// Asserts doubleing point compatibility at compile time
// necessary for INFINITY usage
static_assert(numeric_limits<double>::is_iec559, "IEEE 754 required");
inline double deletion(const double &signal_dp, const int &kmer, const int &sucKmer, vector<tuple<double, double>>& model) {
    const auto &[mean, stddev] = model[kmer];
    const auto &[mean_s, stddev_s] = model[sucKmer];
    const bool aboveCurDist = (signal_dp > mean + 2*stddev);
    const bool belowCurDist = (signal_dp < mean - 2*stddev);
    const bool aboveNextDist = (signal_dp > mean_s + 2*stddev_s);
    const bool belowNextDist = (signal_dp < mean_s - 2*stddev_s);
    
    if ((aboveCurDist && aboveNextDist) || (belowCurDist || belowNextDist)) {
        return min(scoreKmer(mean, kmer, model), scoreKmer(mean_s, sucKmer, model));
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
void logF(const double* sig, const int* kmer_seq, double* M, double* E, double* D, double* I, const size_t T, const size_t N, vector<tuple<double, double>>& model){
    double mat, ext, del, ins, tmp, score;
    E[0] = 0;
    for(size_t t=1; t<T; ++t){
        for(size_t n=1; n<N; ++n){
            mat=-INFINITY;
            ext=-INFINITY;
            del=-INFINITY;
            ins=-INFINITY;
            score = scoreKmer(sig[t-1], kmer_seq[n-1], model);
            
            if (t>C) [[likely]] {
                tmp=E[(t-C-1)*N+(n-1)] + score + m1;
                for(size_t l=1; l<=C; ++l){
                    tmp+=scoreKmer(sig[t-l-1], kmer_seq[n-1], model);
                }
                mat=logPlus(mat, tmp);
            }
            mat=logPlus(mat, D[(t-1)*N+(n-1)] + score + m2);
            mat=logPlus(mat, I[(t-1)*N+(n-1)] + score + m3);
            ext=logPlus(ext, M[(t-1)*N+n]     + score + e1);
            ext=logPlus(ext, E[(t-1)*N+n]     + score + e2);
            ins=logPlus(ins, E[t*N+(n-1)]     + score + i1);
            ins=logPlus(ins, I[t*N+(n-1)]     + score + i2);
            if (n<N-1) [[likely]] {
                score = deletion(sig[t-1], kmer_seq[n-1], kmer_seq[n], model);
                del=logPlus(del, E[(t-1)*N+n] + score + d1);
                del=logPlus(del, D[(t-1)*N+n] + score + d2);
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
void logB(double* sig, int* kmer_seq, double* M, double* E, double* D, double* I, const size_t T, const size_t N, vector<tuple<double, double>>& model) {
    double mat, ext, del, ins, tmp, score;
    for(size_t t=T; t-->0;){
        for(size_t n=N; n-->0;){
            mat=-INFINITY;
            ext=-INFINITY;
            del=-INFINITY;
            ins=-INFINITY;
            if(t==T-1 && n==N-1) {
                ext = 0;
            }

            // m with minimum length C
            if (t+1<T && n+1<N) [[likely]] {
                score = scoreKmer(sig[t], kmer_seq[n], model);

                if (t+C+1<T) [[likely]] {
                    tmp=M[(t+1+C)*N+(n+1)] + score + m1;
                    for (size_t l=1; l<=C; ++l){
                        tmp+=scoreKmer(sig[t+l], kmer_seq[n], model);
                    }
                    ext=logPlus(ext, tmp);
                }

                del=logPlus(del, M[(t+1)*N+(n+1)] + score + m2);
                ins=logPlus(ins, M[(t+1)*N+(n+1)] + score + m3);

                if (t>0) [[likely]] {
                    score = scoreKmer(sig[t-1], kmer_seq[n], model);
                    ext=logPlus(ext, I[t*N+(n+1)] + score + i1);
                    ins=logPlus(ins, I[t*N+(n+1)] + score + i2);
                }
            }

            if (t+1<T && n>0) [[likely]] {
                score = scoreKmer(sig[t], kmer_seq[n-1], model);
                mat=logPlus(mat, E[(t+1)*N+n] + score + e1);
                ext=logPlus(ext, E[(t+1)*N+n] + score + e2);

                if (n+1<N) [[likely]] {
                    score = deletion(sig[t], kmer_seq[n-1], kmer_seq[n], model);
                    ext=logPlus(ext, D[(t+1)*N+n] + score + d1);
                    del=logPlus(del, D[(t+1)*N+n] + score + d2);
                }
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
double* logP(const double* FOR, const double* BACK, const double Z, const size_t T, const size_t N) {
    double* LP = new double[TN];
    for(size_t t=0; t<T; ++t){
        for(size_t n=0; n<N; ++n){
            size_t x = t*N+n;
            LP[x] = FOR[x] + BACK[x] - Z;
        }
    }
    return LP;
}

/**
 * Calculate the maximum a posteriori path (MAP) through LP
 *
 */
list<string> getBorders(double* LPM, double* LPE, double* LPD, double* LPI, const size_t T, const size_t N){
    double* M = new double[TN];
    double* E = new double[TN];
    double* D = new double[TN];
    double* I = new double[TN];
    for (size_t i = 0; i<TN; ++i) {
        M[i] = -INFINITY;
        E[i] = -INFINITY;
        D[i] = -INFINITY;
        I[i] = -INFINITY;
    }
    E[0] = 0;
    double mat, ext, del, ins;
    for(size_t t=1; t<T; ++t){
        for(size_t n=1; n<N; ++n){
            mat=-INFINITY;
            ext=-INFINITY;
            del=-INFINITY;
            ins=-INFINITY;
            const size_t idx = t*N+n;
            if (t>C){
                mat=max(mat, E[(t-C-1)*N+(n-1)] + LPM[idx]); // m1
            }
            ext=max(ext, M[(t-1)*N+n]     + LPE[idx]); // e1
            ext=max(ext, E[(t-1)*N+n]     + LPE[idx]); // e2
            mat=max(mat, D[(t-1)*N+(n-1)] + LPM[idx]); // m2
            mat=max(mat, I[(t-1)*N+(n-1)] + LPM[idx]); // m3
            ins=max(ins, E[t*N+(n-1)]     + LPI[idx]); // i1
            ins=max(ins, I[t*N+(n-1)]     + LPI[idx]); // i2
            if (n<N-1) {
                del=max(del, E[(t-1)*N+n] + LPD[idx]); // d1
                del=max(del, D[(t-1)*N+n] + LPD[idx]); // d2
            }
            M[idx]=mat;
            E[idx]=ext;
            D[idx]=del;
            I[idx]=ins;
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

void funcM(const size_t t, const size_t n, const double* M, const double* E, const double* D, const double* I, const double* LPM, const double* LPE, const double* LPD, const double* LPI, list<string>* segString, const size_t N){
    const double score = M[t*N+n];
    const double logscore = LPM[t*N+n];
    // start value
    if (t<=1 && n<=1){
        segString->push_front("M"+to_string(0)+","+to_string(0)); // n-1 because N is 1 larger than the sequences
        return;
    }
    if (t>0 && n>0) {
        if (t>C && score == E[(t-C-1)*N+(n-1)] + logscore){
            segString->push_front("M"+to_string(n-1)+","+to_string(t-C-1));
            return funcE(t-C-1, n-1, M, E, D, I, LPM, LPE, LPD, LPI, segString, N);
        }
        if (score == D[(t-1)*N+(n-1)] + logscore) {
            segString->push_front("M"+to_string(n-1)+","+to_string(t-1));
            return funcD(t-1, n-1, M, E, D, I, LPM, LPE, LPD, LPI, segString, N);
        }
        if (score == I[(t-1)*N+(n-1)] + logscore) {
            segString->push_front("M"+to_string(n-1)+","+to_string(t-1));
            return funcI(t-1, n-1, M, E, D, I, LPM, LPE, LPD, LPI, segString, N);
        }
    }
}

void funcE(const size_t t, const size_t n, const double* M, const double* E, const double* D, const double* I, const double* LPM, const double* LPE, const double* LPD, const double* LPI, list<string>* segString, const size_t N){
    const double score = E[t*N+n];
    const double logscore = LPE[t*N+n];
    if (t>0) {
        if (score == M[(t-1)*N+n] + logscore){
            return funcM(t-1, n, M, E, D, I, LPM, LPE, LPD, LPI, segString, N);
        }
        if (score == E[(t-1)*N+n] + logscore){
            return funcE(t-1, n, M, E, D, I, LPM, LPE, LPD, LPI, segString, N);
        }
    }
}

void funcI(const size_t t, const size_t n, const double* M, const double* E, const double* D, const double* I, const double* LPM, const double* LPE, const double* LPD, const double* LPI, list<string>* segString, const size_t N){
    const double score = I[t*N+n];
    const double logscore = LPI[t*N+n];
    if (n>0) {
        if (score == I[t*N+(n-1)] + logscore) {
            segString->push_front("I"+to_string(n-1)+","+to_string(t));
            return funcI(t, n-1, M, E, D, I, LPM, LPE, LPD, LPI, segString, N);
        }
        if (score == E[t*N+(n-1)] + logscore) {
            segString->push_front("I"+to_string(n-1)+","+to_string(t));
            return funcE(t, n-1, M, E, D, I, LPM, LPE, LPD, LPI, segString, N);
        }
    }
}

void funcD(const size_t t, const size_t n, const double* M, const double* E, const double* D, const double* I, const double* LPM, const double* LPE, const double* LPD, const double* LPI, list<string>* segString, const size_t N){
    const double score = D[t*N+n];
    const double logscore = LPD[t*N+n];
    if (t>0) {
        if (score == D[(t-1)*N+n] + logscore) {
            segString->push_front("D"+to_string(n)+","+to_string(t-1));
            return funcI(t-1, n, M, E, D, I, LPM, LPE, LPD, LPI, segString, N);
        }
        if (score == E[(t-1)*N+n] + logscore) {
            segString->push_front("D"+to_string(n)+","+to_string(t-1));
            return funcE(t-1, n, M, E, D, I, LPM, LPE, LPD, LPI, segString, N);
        }
    }
}

/**
 * Train transition parameter with baum welch algorithm
*/
tuple<double, double, double, double, double, double, double, double, double> trainTransition(const double* sig, const int* kmer_seq, const double* forM, const double* forE, const double* forD, const double* forI, const double* backM, const double* backE, const double* backD, const double* backI, const size_t T, const size_t N, vector<tuple<double, double>>& model) {
    // Transition parameters
    double newM1 = -INFINITY, newM2 = -INFINITY, newM3 = -INFINITY;
    double newE1 = -INFINITY, newE2 = -INFINITY;
    double newI1 = -INFINITY, newI2 = -INFINITY;
    double newD1 = -INFINITY, newD2 = -INFINITY;
    double tempM = -INFINITY;
    double score;

    for(size_t t=0; t<T; ++t){
        for(size_t n=0; n<N; ++n){
            if (t+1<T && n+1<N) {
                score = scoreKmer(sig[t], kmer_seq[n], model);
                if (t+C+1<T) {
                    // m1:  forward(i)    a    e(i+1)  backward(i+1)
                    tempM = forE[t*N+n] + m1 + score + backM[(t+C+1)*N+(n+1)];
                    for(size_t l=1; l<=C; ++l){
                        tempM+=scoreKmer(sig[t+l], kmer_seq[n], model);
                    }
                    newM1 = logPlus(newM1, tempM);
                }
                newM2 = logPlus(newM2, forD[t*N+n] + m2 + score + backM[(t+1)*N+(n+1)]);
                newM3 = logPlus(newM3, forI[t*N+n] + m3 + score + backM[(t+1)*N+(n+1)]);
            }

            if (t+1<T && n>0) {
                score = scoreKmer(sig[t], kmer_seq[n-1], model);
                newE1 = logPlus(newE1, forM[t*N+n] + e1 + score + backE[(t+1)*N+n]);
                newE2 = logPlus(newE2, forE[t*N+n] + e2 + score + backE[(t+1)*N+n]);
            }

            if (n+1<N && t>0) {
                score = scoreKmer(sig[t-1], kmer_seq[n], model);
                newI1 = logPlus(newI1, forE[t*N+n] + i1 + score + backI[t*N+(n+1)]);
                newI2 = logPlus(newI2, forI[t*N+n] + i2 + score + backI[t*N+(n+1)]);
            }

            if (t+1<T && n>0 && n+1<N) {
                score = deletion(sig[t], kmer_seq[n-1], kmer_seq[n], model);
                newD1 = logPlus(newD1, forE[t*N+n] + d1 + score + backD[(t+1)*N+n]);
                newD2 = logPlus(newD2, forD[t*N+n] + d2 + score + backD[(t+1)*N+n]);
            }
        }
    }
    // average over the number of transitions
    double Am = newE1;
    newE1 -= Am;
    double Ae = logPlus(newI1, logPlus(newD1, logPlus(newE2, newM1)));
    newM1 -= Ae;
    newE2 -= Ae;
    newI1 -= Ae;
    newD1 -= Ae;
    double Ad = logPlus(newD2, newM2);
    newD2 -= Ad;
    newM2 -= Ad;
    double Ai = logPlus(newI2, newM3);
    newI2 -= Ai;
    newM3 -= Ai;

    return tuple<double, double, double, double, double, double, double, double, double>({exp(newM1), exp(newE1), exp(newE2), exp(newI1), exp(newD1), exp(newD2), exp(newM2), exp(newI2), exp(newM3)});
}

/**
 * Train emission parameter with baum welch algorithm
*/
// tuple<double*, double*> trainEmission(double* sig, int* kmer_seq, double* forM, double* forE, double* forD, double* forI, double* backM, double* backE, double* backD, double* backI, const size_t T, const size_t N, vector<tuple<double, double>>& model) {
tuple<double*, double*> trainEmission(const double* sig, const int* kmer_seq, const double* LPM, const double* LPE, const double* LPI, const size_t T, const size_t N, vector<tuple<double, double>>& model) {
    // Emission
    // https://courses.grainger.illinois.edu/ece417/fa2021/lectures/lec15.pdf
    // https://f.hubspotusercontent40.net/hubfs/8111846/Unicon_October2020/pdf/bilmes-em-algorithm.pdf
    // gamma_t(i) is the probability of being in state i at time t
    // gamma for state M - expected number of transitions of M at given time (T) for all latent states (kmers)

    double* kmers = new double[N];
    double* normFactorT = new double[N];
    double* means = new double[numKmers];
    double* stdevs = new double[numKmers];
    int* counts = new int[numKmers];

    for (size_t i = 0; i<N; ++i){
        kmers[i] = 0.0;
        normFactorT[i] = 0.0;
    }
    for (int i = 0; i<numKmers; ++i){
        means[i] = 0.0;
        stdevs[i] = 0.0;
        counts[i] = 0;
    }

    for (size_t n=1; n<N; ++n) {
        counts[kmer_seq[n-1]]++;
        for (size_t t=1; t<T; ++t) {
            double g = exp(logPlus(logPlus(LPM[t*N+n], LPE[t*N+n]), LPI[t*N+n]));
            kmers[n] += g * sig[t-1];
            normFactorT[n] += g;
        }
        kmers[n] = kmers[n] / normFactorT[n];
    }

    for (size_t n=1; n<N; ++n) {
        means[kmer_seq[n-1]] += kmers[n] / counts[kmer_seq[n-1]];
    }

    // Emission (stdev of kmers)
    fill_n(kmers, N, 0.0);
    for (size_t n=1; n<N; ++n) {
        for (size_t t=1; t<T; ++t) {
            double diff = sig[t-1] - means[kmer_seq[n-1]];  // Compute difference from mean
            kmers[n] += exp(logPlus(logPlus(LPM[t*N+n], LPE[t*N+n]), LPI[t*N+n])) * diff * diff;
        }
        kmers[n] = kmers[n] / normFactorT[n];
        stdevs[kmer_seq[n-1]] += kmers[n] / counts[kmer_seq[n-1]];
    }

    for (size_t n=1; n<N; ++n) {
        // transform vars to stdevs
        stdevs[kmer_seq[n-1]] = sqrt(stdevs[kmer_seq[n-1]]);
    }

    delete[] kmers;
    delete[] counts;
    delete[] normFactorT;
    return tuple<double*, double*>({means, stdevs});
}

void trainParams(const double* sig, const int* kmer_seq, const double* forM, const double* forE, const double* forD, const double* forI, const double* backM, const double* backE, const double* backD, const double* backI, const double* LPM, const double* LPE, const double* LPI, const size_t T, const size_t N, vector<tuple<double, double>>& model) {

    auto [newM1, newE1, newE2, newI1, newD1, newD2, newM2, newI2, newM3] = trainTransition(sig, kmer_seq, forM, forE, forD, forI, backM, backE, backD, backI, T, N, model);

    cout<<"m1:"<<newM1<<";e1:"<<newE1<<";e2:"<<newE2<<";i1:"<<newI1<<";d1:"<<newD1<<";d2:"<<newD2<<";m2:"<<newM2<<";i2:"<<newI2<<";m3:"<<newM3<<endl;

    auto [newMeans, newStdevs] = trainEmission(sig, kmer_seq, LPM, LPE, LPI, T, N, model);
    for (int i=0; i<numKmers; i++){
        if (newStdevs[i]!=0.0){
            cout<<itoa(i, kmerSize)<<":"<<newMeans[i]<<","<<newStdevs[i]<<";";
        }
    }
    cout<<endl;
}

/**
 * Read signal and read from stdin until the TERM_STRING is seen
*/
int main(int argc, char* argv[]) {
    bool train, calcZ; // atrain
    int pore;
    const string TERM_STRING = "$";
    string modelpath;

    // Argparser
    argparse::ArgumentParser program("dynamont indel", "0.1");
    // parameters for DP
    // {'e1': 1.0, 'm1': 0.10319949594779426, 'd1': 0.10610745507008552, 'e2': 0.7884055549172662, 'e3': 0.00017942642216108796, 'i1': 0.002108069659237201, 'm2': 0.02937758310704649, 'i2': 0.08849246645562929, 'm3': 0.9115075346587537, 'd2': 0.9706224218081904}
    program.add_argument("-m1", "--matchscore1").help("Transition probability for match rule 1").default_value(.033).scan<'g', double>().store_into(m1); // m1
    program.add_argument("-m2", "--matchscore2").help("Transition probability for match rule 2").default_value(.99).scan<'g', double>().store_into(m2); // m2
    program.add_argument("-m3", "--matchscore3").help("Transition probability for match rule 3").default_value(.99).scan<'g', double>().store_into(m3); // m3
    program.add_argument("-e1", "--extendscore1").help("Transition probability for extend rule 1").default_value(1.00).scan<'g', double>().store_into(e1); // e1
    program.add_argument("-e2", "--extendscore2").help("Transition probability for extend rule 2").default_value(.968).scan<'g', double>().store_into(e2); // e2
    program.add_argument("-d1", "--deletionscore1").help("Transition probability for deletion rule 1").default_value(.001).scan<'g', double>().store_into(d1); // d1
    program.add_argument("-d2", "--deletionscore2").help("Transition probability for deletion rule 2").default_value(.01).scan<'g', double>().store_into(d2); // d2
    program.add_argument("-i1", "--insertionscore1").help("Transition probability for insertion rule 1").default_value(.001).scan<'g', double>().store_into(i1); // i1
    program.add_argument("-i2", "--insertionscore2").help("Transition probability for insertion rule 2").default_value(.01).scan<'g', double>().store_into(i2); // i2
    program.add_argument("-t", "--train").help("Switch algorithm to transition and emission parameter training mode").default_value(false).implicit_value(true).store_into(train);
    program.add_argument("-z", "--calcZ").help("Switch algorithm to only calculate Z").default_value(false).implicit_value(true).store_into(calcZ);
    program.add_argument("-m", "--model").help("Path to kmer model table").default_value("/home/yi98suv/projects/dynamont/data/norm_models/rna_r9.4_180mv_70bps_extended_stdev0_25.model").store_into(modelpath);
    program.add_argument("-r", "--pore").help("Pore generation used to sequence the data").default_value(9).choices(9, 10).store_into(pore).scan<'i', int>();
    program.add_argument("-c", "--minSegLen").help("MinSegLen + 1 is the minimal segment length").default_value(0).scan<'i', int>();

    try {
        program.parse_args(argc, argv);
    }
    catch (const std::runtime_error& err) {
        cerr << err.what() << std::endl;
        cerr << program;
        return 1;
    }
    
    C = program.get<int>("minSegLen");

    if (pore == 9) {
        kmerSize = 5;
    } else if (pore == 10) {
        kmerSize = 9;
    }
    numKmers = pow(ALPHABET_SIZE, kmerSize);
    vector<tuple<double, double>> model(numKmers, make_tuple(-INFINITY, -INFINITY));
    readKmerModel(modelpath, model, kmerSize);
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
        
        // process signal: convert string to double array
        T = count(signal.begin(), signal.end(), ',')+2; // len(sig) + 1
        double* sig = new double[T-1];
        fill_n(sig, T-1, -INFINITY);
        string value;
        stringstream ss(signal);
        int i = 0;
        while(getline(ss, value, ',')) {
            sig[i++] = stod(value);
        }
        // process read: convert string to int array
        N = read.size() + 1; // operate on base transitions
        int seq_size = read.size() + (kmerSize-1); 
        int* seq = new int[seq_size];
        fill_n(seq, seq_size, 4); // default: fill with N add 2 Ns to 3' of read
        i = floor(kmerSize/2);
        for (const char &c: read) {
            seq[i] = BASE2ID.at(c);
            i++;
        }
        // add NN to end of sequence
        seq[i] = 4;
        seq[i+1] = 4;

        int* kmer_seq = seq2kmer(seq, N-1, kmerSize);
        TN = T*N;

        // cerr<<"T: "<<T<<", "<<"N: "<<N<<", "<<"inputsize: "<<TN<<endl;

        // initialize for matrices
        double* forM = new double[TN];
        double* forE = new double[TN];
        double* forD = new double[TN];
        double* forI = new double[TN];
        // initialize back matrices
        double* backM = new double[TN];
        double* backE = new double[TN];
        double* backD = new double[TN];
        double* backI = new double[TN];
        for (size_t i = 0; i<TN; ++i) {
            forM[i] = -INFINITY;
            forE[i] = -INFINITY;
            forD[i] = -INFINITY;
            forI[i] = -INFINITY;
            backM[i] = -INFINITY;
            backE[i] = -INFINITY;
            backD[i] = -INFINITY;
            backI[i] = -INFINITY;
        }
        // calculate segmentation probabilities, fill forward matrices
        logF(sig, kmer_seq, forM, forE, forD, forI, T, N, model);
        // calculate segmentation probabilities, fill backward matrices
        logB(sig, kmer_seq, backM, backE, backD, backI, T, N, model);

        double Zf = forE[TN-1];
        double Zb = backE[0];

        // Numeric error is scaled by input size, Z in forward and backward should match by some numeric error EPSILON
        if (abs(Zf-Zb)/TN > EPSILON) {
            cerr << fixed << showpoint;
            cerr << setprecision(20);
            cerr<<"Z values between matrices do not match! forE[T*N-1]: "<<Zf<<", Zb: "<<Zb<<", "<<abs(Zf-Zb)/TN<<" > "<<EPSILON<<endl;
            cerr.flush();
            exit(11);
        }

        // cerr<<"Zf: "<<Zf<<", Zb: "<<Zb<<", "<<abs(Zf-Zb)/TN<<" <! "<<EPSILON<<endl;

        if (calcZ){
            cout<<Zf<<"\n";
            cout.flush();
        } else {
            double* LPM = logP(forM, backM, Zf, T, N); // log probs
            double* LPE = logP(forE, backE, Zf, T, N); // log probs
            double* LPD = logP(forD, backD, Zf, T, N); // log probs
            double* LPI = logP(forI, backI, Zf, T, N); // log probs
            
            // train both Transitions and Emissions
            if (train) {
                trainParams(sig, kmer_seq, forM, forE, forD, forI, backM, backE, backD, backI, LPM, LPE, LPI, T, N, model);
                cout<<"Z:"<<Zb<<endl;
                cout.flush();
            } else {
                list<string> segString = getBorders(LPM, LPE, LPD, LPI, T, N);

                for (auto const& seg : segString) {
                    cout<<seg;
                }
                cout<<endl;
                cout.flush();

                // Clean up
                delete[] LPM;
                delete[] LPE;
                delete[] LPD;
                delete[] LPI;
            }
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
