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

inline constexpr int ALPHABET_SIZE = 4;
int stepSize, kmerSize; // our model works with this kmer size
inline constexpr double EPSILON = 1e-8; // chose by eye just to distinguish real errors from numeric errors
double m, e; // transition parameters
size_t K, T, TK;

// Asserts doubleing point compatibility at compile time
// necessary for INFINITY usage
static_assert(numeric_limits<double>::is_iec559, "IEEE 754 required");

/**
 * Calculate forward matrices using logarithmic values
 *
 * @param sig ONT raw signal with pA values
 * @param seq nucleotide sequence represented by the ONT signal
 * @param T length of the ONT raw signal + 1
 * @param N length of nucleotide sequence + 1
 * @param model map containing kmers as keys and (mean, stdev) tuples as values
 */
void logF(const double* sig, double* M, double* E, vector<tuple<double, double>> &model){
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
                for(size_t preKmer=precessingKmer(k, 0, stepSize, ALPHABET_SIZE); preKmer<K; preKmer+=stepSize){
                    // mat=logPlus(mat, E[prevTK+preKmer] + scoreKmer(sig[t-1], preKmer, model) + m);
                    mat=logPlus(mat, E[prevTK+preKmer] + score + m);
                }
                ext=logPlus(ext, M[prevTK+k] + score); // transition probability is always 1, e1 first extend
                ext=logPlus(ext, E[prevTK+k] + score + e); // e2 extend further
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
 * @param seq nucleotide sequence represented by the ONT signal
 * @param T length of the ONT raw signal + 1
 * @param N length of nucleotide sequence + 1
 * @param model map containing kmers as keys and (mean, stdev) tuples as values
 */
void logB(const double* sig, double* M, double* E, vector<tuple<double, double>> &model) {
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
                const size_t startKmer = successingKmer(k, 0, stepSize, ALPHABET_SIZE);
                const size_t endKmer = startKmer + ALPHABET_SIZE;
                for(size_t sucKmer = startKmer; sucKmer<endKmer; ++sucKmer){
                    ext=logPlus(ext, M[nexttK+sucKmer] + scoreKmer(sig[t], sucKmer, model) + m);
                }
                const double score = scoreKmer(sig[t], k, model);
                mat=logPlus(mat, E[nexttK+k] + score); // transition probability is always 1, e1 first extend
                ext=logPlus(ext, E[nexttK+k] + score + e); // e2 extend further
            }
            M[tK+k] = mat;
            E[tK+k] = ext;
        }
    }
}

/**
 * Train transition parameter with baum welch algorithm
*/
tuple<double, double, double> trainTransition(const double* sig, double* forM, double* forE, double* backM, double* backE, vector<tuple<double, double>> &model) {
    // Transition parameters
    double newM1 = -INFINITY, newE1 = 0, newE2 = -INFINITY;

    for (size_t t=0; t<T-1; ++t){
        const size_t tK = t*K;
        const size_t nexttK = tK+K; // (t+1)*K
        for (size_t k=0; k<K; ++k){
            const size_t startKmer = successingKmer(k, 0, stepSize, ALPHABET_SIZE);
            const size_t endKmer = startKmer + ALPHABET_SIZE;
            for(size_t sucKmer = startKmer; sucKmer<endKmer; ++sucKmer){
                newM1 = logPlus(newM1, forE[tK+k] + m + scoreKmer(sig[t], sucKmer, model) + backM[nexttK+sucKmer]);
            }
            const double score = scoreKmer(sig[t], k, model);
            // newE1=logPlus(newE1, forM[tK+k] + e1 + score + backE[nexttK+k]); // e1 first extend
            newE2=logPlus(newE2, forE[tK+k] + e + score + backE[nexttK+k]); // e2 extend further
        }
    }
    // newM1 -= log(ALPHABET_SIZE);
    // average over the number of transitions
    // double Am = newE1;
    // newE1 -= Am;
    double Ae = logPlus(newE2, newM1);
    newM1 -= Ae;
    newE2 -= Ae;

    return tuple<double, double, double>({exp(newM1), exp(newE1), exp(newE2)});
}

void trainParams(const double* sig, double* forM, double* forE, double* backM, double* backE, vector<tuple<double, double>> &model) {
    auto [newM, newE1, newE2] = trainTransition(sig, forM, forE, backM, backE, model);

    // cout<<"m1:"<<newM<<";e1:"<<newE1<<";e2:"<<newE2<<endl;
    cout<<"m1:"<<newM<<endl;

    // placeholder for Emission params, is empty here cause I did not train them
    cout<<endl;
}

/**
 * Read signal and read from stdin until the TERM_STRING is seen
*/
int main(int argc, char* argv[]) {
    int pore;
    bool train, calcZ, prob;
    string modelpath;
    const string TERM_STRING = "$";

    // Argparser
    argparse::ArgumentParser program("dynamont TK", "0.1");
    // parameters for DP
    // program.add_argument("-e1", "--extendscore1").help("Transition probability for extend rule 1").default_value(1.00).scan<'g', double>().store_into(e1); // e1
    program.add_argument("-m1", "--matchscore").help("Transition probability for match rule").default_value(.03).scan<'g', double>().store_into(m); // m
    // program.add_argument("-e2", "--extendscore2").help("Transition probability for extend rule 2").default_value(.97).scan<'g', double>().store_into(e2); // e2
    program.add_argument("-t", "--train").help("Switch algorithm to transition and emission parameter training mode").default_value(false).implicit_value(true).store_into(train);
    program.add_argument("-z", "--calcZ").help("Switch algorithm to only calculate Z").default_value(false).implicit_value(true).store_into(calcZ);
    program.add_argument("-m", "--model").help("Path to kmer model table").default_value("/home/yi98suv/projects/dynamont/data/norm_models/rna_r9.4_180mv_70bps.model").store_into(modelpath);
    program.add_argument("-r", "--pore").help("Pore generation used to sequence the data").default_value(9).choices(9, 10).scan<'i', int>().store_into(pore);
    program.add_argument("-p", "--probabilty").help("Print out the segment border probability").default_value(false).implicit_value(true).store_into(prob);
    // only for compatibility
    program.add_argument("-c", "--minSegLen").help("MinSegLen + 1 is the minimal segment length").default_value(0).scan<'i', int>();
    
    e = log(1-m);
    m = log(m);

    try {
        program.parse_args(argc, argv);
    }
    catch (const std::runtime_error& err) {
        cerr << err.what() << std::endl;
        cerr << program;
        return 1;
    }

    if (pore == 9) {
        kmerSize = 5;
    } else if (pore == 10) {
        kmerSize = 9;
    }
    K = pow(ALPHABET_SIZE, kmerSize); // currently acceptable A, C, G, T, N
    stepSize = pow(ALPHABET_SIZE, kmerSize-1);
    vector<tuple<double, double>> model(K, make_tuple(-INFINITY, -INFINITY));
    readKmerModel(modelpath, model, ALPHABET_SIZE);
    string signal, read;
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
        // ... or signal is missing
        } else if (signal.empty()) {
            cout<<"Signal missing!\n";
            cout.flush();
            return 1;
        }

        // cerr<<"DEBUG 1"<<endl;
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
        TK = T*K;

        // cerr<<"T: "<<T<<", "<<"N: "<<N<<", "<<"inputsize: "<<TK<<endl;

        // initialize matrices
        double* forM = new double[TK];
        double* forE = new double[TK];
        double* backM = new double[TK];
        double* backE = new double[TK];
        for (size_t i = 0; i<TK; ++i) {
            forM[i] = -INFINITY;
            forE[i] = -INFINITY;
            backM[i] = -INFINITY;
            backE[i] = -INFINITY;
        }
        // calculate segmentation probabilities, fill forward matrices
        logF(sig, forM, forE, model);
        // calculate segmentation probabilities, fill backward matrices
        logB(sig, backM, backE, model);

        double Zf = -INFINITY;
        double Zb = -INFINITY;
        for(size_t k=0; k<K; ++k){
            Zf = logPlus(Zf, forE[TK-1-k]);
            Zb = logPlus(Zb, backE[k]);
        }

        // Numeric error is scaled by input size, Z in forward and backward should match by some numeric error EPSILON
        if (abs(Zf-Zb)/TK > EPSILON) {
            cerr << fixed << showpoint;
            cerr << setprecision(20);
            cerr<<"Z values between matrices do not match! Zf: "<<Zf<<", Zb: "<<Zb<<", "<<abs(Zf-Zb)/(TK)<<" > "<<EPSILON<<endl;
            cerr.flush();
            exit(11);
        }

        if (calcZ){
            cout<<Zb<<endl;
            cout.flush();
        } else {
            // train both Transitions and Emissions
            trainParams(sig, forM, forE, backM, backE, model);
            cout<<"Z:"<<Zb<<endl;
            cout.flush();
        }
        delete[] forM;
        delete[] forE;
        delete[] backM;
        delete[] backE;
        delete[] sig;
    }
    return 0;
}
