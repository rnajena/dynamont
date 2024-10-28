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

void funcM(const size_t t, const size_t n, const double* M, const double* E, const double* LPM, const double* LPE, list<string>* segString, const size_t N, vector<double> &segProb);
void funcE(const size_t t, const size_t n, const double* M, const double* E, const double* LPM, const double* LPE, list<string>* segString, const size_t N, vector<double> &segProb);

inline constexpr double EPSILON = 1e-2; // chose by eye just to distinguish real errors from numeric errors
inline constexpr size_t BANDEDDPWINDOW = 150; // Banded DP approach: Window size around the maximum, adjust as needed

size_t N, T, TN;
int alphabet_size, numKmers, kmerSize; // our model works with this kmer size

unordered_map<string, double> transitions = {
    {"m1", -1.0},
    {"e1", -1.0},
    {"e2", -1.0}
};

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
void logF(const double* sig, const int* kmer_seq, double* M, double* E, const size_t T, const size_t N, vector<tuple<double, double>>& model){
    E[0] = 0;
    // Initialize max indices for the first iteration
    size_t max_idx = 0;

    for (size_t t=1; t<T; ++t){
        double max_val = -INFINITY;
        const size_t tN = t*N;

        // Determine the range around the max index
        size_t start_n = (max_idx > BANDEDDPWINDOW) ? max_idx - BANDEDDPWINDOW : 1;
        size_t end_n = (max_idx + BANDEDDPWINDOW < N) ? max_idx + BANDEDDPWINDOW : N - 1;

        for (size_t n=start_n; n<=end_n; ++n){  // Use restricted range
            const double score = scoreKmer(sig[t-1], kmer_seq[n-1], model);  // Cache scoreKmer for (t-1, n-1)

            // Update matrices
            M[tN+n] = E[(t-1)*N+(n-1)] + score + transitions["m1"];
            E[tN+n] = logPlus(M[tN-N+n] + score, E[tN-N+n] + score + transitions["e2"]);

            // Track maximum values for the next iteration
            // Update max indices for the next iteration
            if (M[tN+n] > max_val) {
                max_val = M[tN+n];
                max_idx = n;
            }
            if (E[tN+n] > max_val) {
                max_val = E[tN+n];
                max_idx = n;
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
void logB(const double* sig, const int* kmer_seq, double* M, double* E, const size_t T, const size_t N, vector<tuple<double, double>>& model) {
    E[TN-1] = 0;
    // Initialize max indices for the first iteration (starting at t=T-1)
    size_t max_idx = N-1;

    // Start backward iteration
    for (size_t t=T-1; t-->0;){
        double max_val = -INFINITY;
        const size_t tN = t*N;

        // Define the range of n based on the previous max_idx
        size_t start_n = (max_idx > BANDEDDPWINDOW) ? max_idx - BANDEDDPWINDOW : 0;
        size_t end_n = (max_idx + BANDEDDPWINDOW < N) ? max_idx + BANDEDDPWINDOW + 1 : N;

        for (size_t n = end_n; n-->start_n;){  // Use restricted range and iterate backwards
            double mat = -INFINITY;
            double ext = -INFINITY;

            if (t+1 < T && n+1 < N) [[likely]] {
                ext = logPlus(ext, M[(t+1)*N+(n+1)] + scoreKmer(sig[t], kmer_seq[n], model) + transitions["m1"]);
            }

            if (t+1<T && n > 0) [[likely]] {
                const double score = scoreKmer(sig[t], kmer_seq[n-1], model);
                mat = logPlus(mat, E[tN+N+n] + score);  // e1 first extend
                ext = logPlus(ext, E[tN+N+n] + score + transitions["e2"]);  // e2 extend further
            }

            // Store computed values in matrices
            M[tN+n] = mat;
            E[tN+n] = ext;

            // Update the maximum values for the next iteration
            if (mat > max_val) {
                max_val = mat;
                max_idx = n;
            }
            if (ext > max_val) {
                max_val = ext;
                max_idx = n;
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
 * @param Z alignment score
 * @return matrix containing logarithmic probabilities for segment borders
 */
double* logP(const double* FOR, const double* BACK, const double Z, const size_t N) {
    double* LP = new double[TN];
    for (size_t t=0; t<T; ++t){
        for (size_t n=0; n<N; ++n){
            const size_t x = t*N+n;
            LP[x] = FOR[x] + BACK[x] - Z;
        }
    }
    return LP;
}

/**
 * Calculate the maximum a posteriori path through LP
 *
 */
list<string> getBorders(const double* LPM, const double* LPE, const size_t T, const size_t N){
    double* M = new double[TN];
    double* E = new double[TN];
    // Initialize M and E in one step, no need for fill_n
    for (size_t i = 0; i < TN; ++i) {
        M[i] = -INFINITY;
        E[i] = -INFINITY;
    }

    E[0] = 0;
    for (size_t t=1; t<T; ++t){
        for (size_t n=1; n<=t && n<N; ++n){ // speed up, due to rules no need to look at upper triangle of matrices
            M[t*N+n]=E[(t-1)*N+(n-1)] + LPM[t*N+n]; //m1
            E[t*N+n]=max(M[(t-1)*N+n] + LPE[t*N+n], E[(t-1)*N+n] + LPE[t*N+n]); //e1, e2
        }
    }
    list<string> segString;
    vector<double> segProb;
    funcE(T-1, N-1, M, E, LPM, LPE, &segString, N, segProb);
    delete[] M;
    delete[] E;
    return segString;
}

void funcM(const size_t t, const size_t n, const double* M, const double* E, const double* LPM, const double* LPE, list<string>* segString, const size_t N, vector<double> &segProb){
    double score = M[t*N+n];
    double logScore = LPM[t*N+n];
    segProb.push_back(exp(logScore));
    if (t<=1 && n<=1){ // Start value
        segString->push_front("M0,0," + to_string(calculateMedian(segProb)) + ";"); // n-1 because N is 1 larger than the sequences
        return;
    }
    if (t>0 && n>0 && score == E[(t-1)*N+(n-1)] + logScore){
        segString->push_front("M" + to_string(n-1+kmerSize/2) + "," + to_string(t-1) + "," + to_string(calculateMedian(segProb)) + ";");
        segProb.clear();
        return funcE(t-1, n-1, M, E, LPM, LPE, segString, N, segProb);
    }
}

void funcE(const size_t t, const size_t n, const double* M, const double* E, const double* LPM, const double* LPE, list<string>* segString, const size_t N, vector<double> &segProb){
    double score = E[t*N+n];
    double logScore = LPE[t*N+n];
    segProb.push_back(exp(logScore));
    if (t>0 && n>0) {
        if (score == M[(t-1)*N+n] + logScore){
            return funcM(t-1, n, M, E, LPM, LPE, segString, N, segProb);
        }
        if (score == E[(t-1)*N+n] + logScore){
            return funcE(t-1, n, M, E, LPM, LPE, segString, N, segProb);
        }
    }
}

/**
 * Train transition parameter with baum welch algorithm
*/
tuple<double, double, double> trainTransition(const double* sig, const int* kmer_seq, double* forM, double* forE, double* backM, double* backE, const size_t T, const size_t N, vector<tuple<double, double>>& model) {
    // Transition parameters
    double newM1 = -INFINITY, newE1 = 0, newE2 = -INFINITY;

    for (size_t t=0; t<T; ++t){
        for (size_t n=0; n<=t && n<N; ++n){ // speed up, due to rules no need to look at upper triangle of matrices
            if (n+1<N && t+1<T) {
                // m1:                 forward(i)    a                   e(i+1)                                  backward(i+1)
                newM1 = logPlus(newM1, forE[t*N+n] + transitions["m1"] + scoreKmer(sig[t], kmer_seq[n], model) + backM[(t+1)*N+(n+1)]);
            }

            if (t+1<T && n>0) {
                double score = scoreKmer(sig[t], kmer_seq[n-1], model);
                // newE1 = logPlus(newE1, forM[t*N+n] + score + backE[(t+1)*N+n]);
                newE2 = logPlus(newE2, forE[t*N+n] + transitions["e2"] + score + backE[(t+1)*N+n]);
            }
        }
    }
    // average over the number of transitions
    // double Am = newE1;
    // newE1 = newE1 - Am;
    double Ae = logPlus(newE2, newM1);
    newM1 = newM1 - Ae;
    newE2 = newE2 - Ae;

    return tuple<double, double, double>({exp(newM1), exp(newE1), exp(newE2)});
}

/**
 * Train emission parameter with baum welch algorithm
*/
tuple<double*, double*> trainEmission(const double* sig, const int* kmer_seq, double* forM, double* forE, double* backM, double* backE, const size_t T, const size_t N, vector<tuple<double, double>>& model) {

    // TODO calculate this with LP matrices
    // see notes

    // Emission
    // https://courses.grainger.illinois.edu/ece417/fa2021/lectures/lec15.pdf
    // https://f.hubspotusercontent40.net/hubfs/8111846/Unicon_October2020/pdf/bilmes-em-algorithm.pdf
    // gamma_t(i) is the probability of being in state i at time t
    // gamma for state M - expected number of transitions of M at given time (T) for all latent states (kmers)
    // Initialize memory
    double* G = new double[TN];
    double* kmers = new double[N];
    double* d = new double[N];
    double* means = new double[numKmers];
    double* stdevs = new double[numKmers];
    int* counts = new int[numKmers];

    // init everything with zero
    for (size_t i = 0; i < N; i++) {
        kmers[i] = 0.0;
        d[i] = 0.0;
    }
    for (int i = 0; i < numKmers; i++) {
        means[i] = 0.0;
        stdevs[i] = 0.0;
        counts[i] = 0;
    }

    for (size_t t=0; t<T; ++t){
        // calibrate with the sum of transitions
        double s = -INFINITY;
        for (size_t n=0; n<=t && n<N; ++n){ // speed up, due to rules no need to look at upper triangle of matrices
            G[t*N+n] = logPlus(forM[t*N+n] + backM[t*N+n], forE[t*N+n] + backE[t*N+n]);
            s = logPlus(s, G[t*N+n]);
        }
        for (size_t n=0; n<=t && n<N; ++n){
            if (!isinf(s)) {
                G[t*N+n] -= s;
            }
        }
    }
    
    for (size_t n=1; n<N; ++n) {
        counts[kmer_seq[n-1]]++;
        for (size_t t=n; t<T; ++t) {        // speed up, due to rules no need to look at upper triangle of matrices
            double g = exp(G[t * N + n]);   // Cache the exp(G[t * N + n])
            kmers[n] += g * sig[t-1];       // Accumulate for kmers
            d[n] += g;                      // Accumulate for normalizer
        }
        if (d[n] > 0) {
            kmers[n] = kmers[n] / d[n];         // Normalize
        }
    }

    for (size_t n=1; n<N; ++n) {
        means[kmer_seq[n-1]] += kmers[n] / counts[kmer_seq[n-1]]; // Update means
    }

    // Emission (stdev of kmers)
    fill_n(kmers, N, 0);
    for (size_t n=1; n<N; ++n) {
        for (size_t t=n; t<T; ++t) { // n<=t, speed up, due to rules no need to look at upper triangle of matrices
            double diff = sig[t-1] - means[kmer_seq[n-1]];  // Compute difference from mean
            kmers[n] += exp(G[t * N + n]) * diff * diff;  // Accumulate for variance
        }
        if (d[n]>0) {
            kmers[n] = kmers[n] / d[n];
        }
        stdevs[kmer_seq[n-1]] += kmers[n] / counts[kmer_seq[n-1]];
    }

    for (size_t n=1; n<N; ++n) {
        // transform vars to stdevs
        stdevs[kmer_seq[n-1]] = sqrt(stdevs[kmer_seq[n-1]]);
    }

    delete[] G;
    delete[] kmers;
    delete[] counts;
    delete[] d;
    return tuple<double*, double*>({means, stdevs});
}

void trainParams(const double* sig, const int* kmer_seq, double* forM, double* forE, double* backM, double* backE, const size_t T, const size_t N, vector<tuple<double, double>>& model) {
    auto [newM, newE1, newE2] = trainTransition(sig, kmer_seq, forM, forE, backM, backE, T, N, model);

    cout<<"m1:"<<newM<<";e1:"<<newE1<<";e2:"<<newE2<<endl;

    auto [newMeans, newStdevs] = trainEmission(sig, kmer_seq, forM, forE, backM, backE, T, N, model);
    for (int i=0; i<numKmers; i++){
        if (newStdevs[i]!=0.0){
            cout<<itoa(i, alphabet_size, kmerSize)<<":"<<newMeans[i]<<","<<newStdevs[i]<<";";
        }
    }
    cout<<endl;
}

/**
 * Read signal and read from stdin until the TERM_STRING is seen
*/
int main(int argc, char* argv[]) {
    // speedup for I/O
    std::ios_base::sync_with_stdio(0);
    std::cin.tie(0);
    std::cout.tie(0);

    bool train, calcZ, prob;
    string pore, modelpath;

    // std::cerr precisions
    std::cerr << std::fixed << std::showpoint << std::setprecision(5);
    // std::cerr precisions
    std::cout << std::fixed << std::showpoint << std::setprecision(5);

    // Argparser
    argparse::ArgumentParser program("dynamont basic", "0.1");
    // parameters for DP
    program.add_argument("-m1", "--matchscore1").help("Segment transition probability, should be close to (expected number of nucleotdes)/(signal length). Leave at -1 if unset.").default_value(-1.0).scan<'g', double>().store_into(transitions["m1"]);
    program.add_argument("-e1", "--extendscore1").help("First extend probability.").default_value(-1.0).scan<'g', double>().store_into(transitions["e1"]);
    program.add_argument("-e2", "--extendscore2").help("Further extend probability.").default_value(-1.0).scan<'g', double>().store_into(transitions["e2"]);
    program.add_argument("-t", "--train").help("Switch algorithm to transition and emission parameter training mode").default_value(false).implicit_value(true).store_into(train);
    program.add_argument("-z", "--calcZ").help("Switch algorithm to only calculate Z").default_value(false).implicit_value(true).store_into(calcZ);
    program.add_argument("-m", "--model").help("Path to kmer model table").default_value("/home/yi98suv/projects/dynamont/data/norm_models/rna_r9.4_180mv_70bps.model").store_into(modelpath);
    program.add_argument("-r", "--pore").help("Pore used to sequence the data").default_value("rna_r9").choices("rna_r9", "dna_r9", "rna_rp4", "dna_r10").store_into(pore);
    program.add_argument("-p", "--probabilty").help("Print out the segment border probability").default_value(false).implicit_value(true).store_into(prob);

    try {
        program.parse_args(argc, argv);
    }
    catch (const std::runtime_error& err) {
        cerr << err.what() << std::endl;
        cerr << program;
        return 1;
    }

    // load default and set parameters
    if (pore == "rna_r9") {
        kmerSize = 5;
        // taken from the trained NT version of dynamont
        updateTransitions(NT_rna_r9_transitions, transitions);
    } else if (pore == "dna_r9") {
        kmerSize = 5;
        // taken from the trained NT version of dynamont
        updateTransitions(NT_dna_r9_transitions, transitions);
    } else if (pore == "rna_rp4") {
        kmerSize = 9;
        // taken from the trained NT version of dynamont
        updateTransitions(NT_rna_rp4_transitions, transitions);
    } else if (pore == "dna_r10_260bps") {
        kmerSize = 9;
        updateTransitions(NT_dna_r10_260bps_transitions, transitions);
    } else if (pore == "dna_r10_400bps") {
        kmerSize = 9;
        updateTransitions(NT_dna_r10_400bps_transitions, transitions);
    }

    // check that outgoing transitions sum up to 1
    // assert(fabs(exp(transitions["m1"]) + exp(transitions["e2"]) - 1.0) < 1e-2 && "The sum of the outgoing transitions of state E: m1 and e1 must approximately 1.0");

    // vector<tuple<double, double>> model(numKmers, make_tuple(-INFINITY, -INFINITY));
    assert(!modelpath.empty() && "Please provide a modelpath!");
    auto result = readKmerModel(modelpath, kmerSize);
    vector<tuple<double, double>> model = get<0>(result);
    alphabet_size = get<1>(result);
    numKmers = get<2>(result);
    string signal, read;

    // echo 107,107,107.2,108.0,108.9,111.2,105.7,104.3,107.1,105.7,105,105 CAAAAA| src\segment.exe
    // read input, signal and read whitespace separated in single line
    getline(cin, signal);
    getline(cin, read);

    if (signal.empty()) {
        cout<<"Signal missing!"<<endl;
        return 1;
    } else if (read.empty()) {
        cout<<"Read missing!"<<endl;
        return 2;
    }

    // cerr<<"DEBUG 1"<<endl;
    // process signal: convert string to double array
    T = count(signal.begin(), signal.end(), ',')+2; // len(sig) + 1
    double* sig = new double[T-1];
    string value;
    stringstream ss(signal);
    int i = 0;
    while(getline(ss, value, ',')) {
        sig[i++] = stof(value);
    }
    
    // process read N: convert string to int array
    N = read.size() - kmerSize + 1 + 1; // N is number of kmers in sequence + 1
    int* kmer_seq = new int[N-1];
    for (size_t n=0; n<N-1; ++n) {
        kmer_seq[n] = kmer2int(read.substr(n, kmerSize), alphabet_size);
    }

    // deallocate memory
    ss.clear();
    signal.erase();
    read.erase();
    value.erase();

    TN = T*N;

    // cerr<<"T: "<<T<<", "<<"N: "<<N<<", "<<"inputsize: "<<TN<<endl;

    // initialize matrices
    double* forM = new double[TN];
    double* forE = new double[TN];
    double* backM = new double[TN];
    double* backE = new double[TN];
    for (size_t i = 0; i<TN; ++i) {
        forM[i] = -INFINITY;
        forE[i] = -INFINITY;
        backM[i] = -INFINITY;
        backE[i] = -INFINITY;
    }
    // calculate segmentation probabilities, fill forward matrices
    logF(sig, kmer_seq, forM, forE, T, N, model);
    // calculate segmentation probabilities, fill backward matrices
    logB(sig, kmer_seq, backM, backE, T, N, model);

    const double Zf = forE[TN-1];
    const double Zb = backE[0];

    // Numeric error is scaled by input size, Z in forward and backward should match by some numeric error EPSILON
    if (abs(Zf-Zb)/TN > EPSILON || isinf(Zf) || isinf(Zb)) {
        cerr<<"Z values between matrices do not match! Zf: "<<Zf<<", Zb: "<<Zb<<", "<<abs(Zf-Zb)/(TN)<<" > "<<EPSILON<<", T: "<<T<<", N: "<<N<<endl;
        exit(11);
    }
    
    // cerr<<"Zf: "<<Zf<<", Zb: "<<Zb<<", "<<abs(Zf-Zb)/TN<<" <! "<<EPSILON<<endl;
    const double Z = max(Zf, Zb);

    if (calcZ){
        cout<<Z<<endl;
    } else {
        
        // train both Transitions and Emissions
        if (train) {
            trainParams(sig, kmer_seq, forM, forE, backM, backE, T, N, model);
            cout<<"Z:"<<Z<<endl;

        } else {
            const double* LPM = logP(forM, backM, Z, N); // log probs for segmentation
            const double* LPE = logP(forE, backE, Z, N); // log probs for extension
            const list<string> segString = getBorders(LPM, LPE, T, N);

            for (auto const& seg : segString) {
                cout<<seg;
            }
            cout<<endl;

            // calculate sum of segment probabilities
            if (prob) {
                double sum;
                for (size_t t=0; t<T; ++t) {
                    sum = -INFINITY;
                    for (size_t n=0; n<N; ++n) {
                        sum = logPlus(sum, LPM[t*N+n]);
                    }
                    cout<<sum<<",";
                }
                cout<<endl;
            }

            // Clean up
            delete[] LPM;
            delete[] LPE;
        }
    }
    delete[] forM;
    delete[] forE;
    delete[] backM;
    delete[] backE;
    delete[] sig;
    // delete[] seq;
    delete[] kmer_seq;
    // }
    return 0;
}
