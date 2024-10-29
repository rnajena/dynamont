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

void funcM(const size_t t, const size_t n, const dproxy* M, const dproxy* E, const double* LPM, const double* LPE, list<string>* segString, const size_t N, vector<double> &segProb, const int kmerSize);
void funcE(const size_t t, const size_t n, const dproxy* M, const dproxy* E, const double* LPM, const double* LPE, list<string>* segString, const size_t N, vector<double> &segProb, const int kmerSize);

constexpr double EPSILON = 1e-8; // chose by eye just to distinguish real errors from numeric errors

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
tuple<dproxy*, dproxy*> logF(const double* sig, const int* kmer_seq, const size_t T, const size_t N, const vector<tuple<double, double>>& model){
    dproxy* M = new dproxy[T*N];
    dproxy* E = new dproxy[T*N];
    const double m1 = transitions.at("m1");
    const double e2 = transitions.at("e2");
    E[0] = 0;
    for (size_t t=1; t<T; ++t){
        const size_t tN = t*N;
        for (size_t n=1; n<=t && n<N; ++n){ // speed up, due to rules no need to look at upper triangle of matrices
            const double score = scoreKmer(sig[t-1], kmer_seq[n-1], model); // Cache scoreKmer for (t-1, n-1)
            M[tN+n] = E[tN-N+(n-1)] + score + m1;
            E[tN+n] = logPlus(M[tN-N+n] + score, E[tN-N+n] + score + e2);
        }
    }
    return make_tuple(M, E);
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
tuple<dproxy*, dproxy*> logB(const double* sig, const int* kmer_seq, const size_t T, const size_t N, const vector<tuple<double, double>>& model) {
    dproxy* M = new dproxy[T*N];
    dproxy* E = new dproxy[T*N];
    const double m1 = transitions.at("m1");
    const double e2 = transitions.at("e2");
    E[T*N-1] = 0;
    for (size_t t=T-1; t-->0;){
        const size_t tN = t*N;
        for (size_t n=min(N, t+1); n-->0;){ // speed up, due to rules no need to look at upper triangle of matrices
            double mat=-INFINITY, ext=-INFINITY;
            if (n+1<N) [[likely]] {
                ext = logPlus(ext, M[tN+N+(n+1)] + scoreKmer(sig[t], kmer_seq[n], model) + m1);
            }
            if (n>0) [[likely]] {
                const double score = scoreKmer(sig[t], kmer_seq[n-1], model);
                mat=logPlus(mat, E[tN+N+n] + score); // transition probability is always 1, e1 first extend
                ext=logPlus(ext, E[tN+N+n] + score + e2); // e2 extend further
            }

            M[tN+n] = mat;
            E[tN+n] = ext;
        }
    }
    return make_tuple(M, E);
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
double* logP(const dproxy* FOR, const dproxy* BACK, const double Z, const size_t TN) {
    double* LP = new double[TN];
    for (size_t tn=0; tn<TN; ++tn){
        LP[tn] = FOR[tn] + BACK[tn] - Z;
    }
    return LP;
}

/**
 * Calculate the maximum a posteriori path through LP
 *
 */
list<string> getBorders(const double* LPM, const double* LPE, const size_t T, const size_t N, const int kmerSize){
    const size_t TN = T*N;
    dproxy* M = new dproxy[TN];
    dproxy* E = new dproxy[TN];
    E[0] = 0;
    for (size_t t=1; t<T; ++t){
        for (size_t n=1; n<=t && n<N; ++n){ // speed up, due to rules no need to look at upper triangle of matrices
            M[t*N+n]=E[(t-1)*N+(n-1)] + LPM[t*N+n]; // m1
            E[t*N+n]=max(M[(t-1)*N+n] + LPE[t*N+n], E[(t-1)*N+n] + LPE[t*N+n]); //e1, e2
        }
    }
    list<string> segString;
    vector<double> segProb;
    funcE(T-1, N-1, M, E, LPM, LPE, &segString, N, segProb, kmerSize);
    delete[] M;
    delete[] E;
    return segString;
}

void funcM(const size_t t, const size_t n, const dproxy* M, const dproxy* E, const double* LPM, const double* LPE, list<string>* segString, const size_t N, vector<double>& segProb, const int kmerSize){
    const double score = M[t*N+n];
    const double logScore = LPM[t*N+n];
    segProb.push_back(exp(logScore));
    if (t<=1 && n<=1) [[unlikely]] { // Start value
        segString->push_front("M0,0," + to_string(calculateMedian(segProb)) + ";"); // n-1 because N is 1 larger than the sequences
        return;
    }
    if (t>0 && n>0 && score == E[(t-1)*N+(n-1)] + logScore) [[likely]] {
        segString->push_front("M" + to_string(n-1+kmerSize/2) + "," + to_string(t-1) + "," + to_string(calculateMedian(segProb)) + ";");
        segProb.clear();
        return funcE(t-1, n-1, M, E, LPM, LPE, segString, N, segProb, kmerSize);
    }
}

void funcE(const size_t t, const size_t n, const dproxy* M, const dproxy* E, const double* LPM, const double* LPE, list<string>* segString, const size_t N, vector<double>& segProb, const int kmerSize){
    const double score = E[t*N+n];
    const double logScore = LPE[t*N+n];
    segProb.push_back(exp(logScore));
    if (t>0 && n>0) [[likely]] {
        if (score == M[(t-1)*N+n] + logScore){
            return funcM(t-1, n, M, E, LPM, LPE, segString, N, segProb, kmerSize);
        }
        if (score == E[(t-1)*N+n] + logScore){
            return funcE(t-1, n, M, E, LPM, LPE, segString, N, segProb, kmerSize);
        }
    }
}

/**
 * Train transition parameter with baum welch algorithm
*/
tuple<double, double, double> trainTransition(const double* sig, const int* kmer_seq, const dproxy* forM, const dproxy* forE, const dproxy* backM, const dproxy* backE, const size_t T, const size_t N, const vector<tuple<double, double>>& model) {
    // Transition parameters
    double newM1 = -INFINITY, newE1 = 1, newE2 = -INFINITY;

    for (size_t t=0; t<T-1; ++t){
        for (size_t n=0; n<=t && n<N; ++n){ // speed up, due to rules no need to look at upper triangle of matrices
            if (n+1<N) [[likely]] {
                // m1:                 forward(i)    a                      e(i+1)                                  backward(i+1)
                newM1 = logPlus(newM1, forE[t*N+n] + transitions.at("m1") + scoreKmer(sig[t], kmer_seq[n], model) + backM[(t+1)*N+(n+1)]);
            }

            if (n>0) [[likely]] {
                double score = scoreKmer(sig[t], kmer_seq[n-1], model);
                // newE1 = logPlus(newE1, forM[t*N+n] + transitions.at("e1") + score + backE[(t+1)*N+n]);
                newE2 = logPlus(newE2, forE[t*N+n] + transitions.at("e2") + score + backE[(t+1)*N+n]);
            }
        }
    }
    // average over the number of transitions
    // double Am = newE1;
    // newE1 = newE1 - Am;
    const double Ae = logPlus(newE2, newM1);
    newM1 = exp(newM1 - Ae);
    newE2 = exp(newE2 - Ae);
    return make_tuple(newM1, newE1, newE2);
}

/**
 * Train emission parameter with baum welch algorithm
*/
tuple<double*, double*> trainEmission(const double* sig, const int* kmer_seq, const dproxy* forM, const dproxy* forE, const dproxy* backM, const dproxy* backE, const size_t T, const size_t N, const vector<tuple<double, double>>& model, const int numKmers) {

    // TODO: calculate this with LP matrices
    // see notes

    // Emission
    // https://courses.grainger.illinois.edu/ece417/fa2021/lectures/lec15.pdf
    // https://f.hubspotusercontent40.net/hubfs/8111846/Unicon_October2020/pdf/bilmes-em-algorithm.pdf
    // gamma_t(i) is the probability of being in state i at time t
    // gamma for state M - expected number of transitions of M at given time (T) for all latent states (kmers)
    // Initialize memory
    const size_t TN = T*N;
    double* G = new double[TN];
    double* kmers = new double[N];
    double* d = new double[N];
    double* means = new double[numKmers];
    double* stdevs = new double[numKmers];
    int* counts = new int[numKmers];

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
        // if (d[n] > 0) {
        kmers[n] = kmers[n] / d[n];         // Normalize
        // }
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
        // if (d[n]>0) {
        kmers[n] = kmers[n] / d[n];
        // }
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
    return make_tuple(means, stdevs);
}

void trainParams(const double* sig, const int* kmer_seq, const dproxy* forM, const dproxy* forE, const dproxy* backM, const dproxy* backE, const size_t T, const size_t N, vector<tuple<double, double>>& model, const int alphabet_size, const int numKmers, const int kmerSize) {
    const auto [newM, newE1, newE2] = trainTransition(sig, kmer_seq, forM, forE, backM, backE, T, N, model);
    cout<<"m1:"<<newM<<";e1:"<<newE1<<";e2:"<<newE2<<endl;
    const auto [newMeans, newStdevs] = trainEmission(sig, kmer_seq, forM, forE, backM, backE, T, N, model, numKmers);
    for (int i=0; i<numKmers; i++){
        if (newStdevs[i]!=0.0) [[likely]] {
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

    std::cerr << std::fixed << std::showpoint << std::setprecision(5);
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
    program.add_argument("-r", "--pore").help("Pore used to sequence the data").default_value("rna_r9").choices("rna_r9", "dna_r9", "rna_rp4", "dna_r10_260bps", "dna_r10_400bps").store_into(pore);
    program.add_argument("-p", "--probabilty").help("Print out the segment border probability").default_value(false).implicit_value(true).store_into(prob);

    try {
        program.parse_args(argc, argv);
    }
    catch (const std::runtime_error& err) {
        cerr << err.what() << endl;
        cerr << program;
        return 1;
    }

    int kmerSize = 0;
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
    } else {
        cerr << "Unsupported pore: " << pore << endl;
        return 3;
    }

    assert(!modelpath.empty() && "Please provide a modelpath!");
    auto result = readKmerModel(modelpath, kmerSize);
    vector<tuple<double, double>> model = get<0>(result);
    const int alphabet_size = get<1>(result);
    const int numKmers = get<2>(result);
    
    // example 
    // 107,107,107.2,108.0,108.9,111.2,105.7,104.3,107.1,105.7,105,105
    // CAAAAA
    // read input, signal and read whitespace separated in single line
    string signal, read;
    getline(cin, signal);
    getline(cin, read);

    if (signal.empty()) {
        cout<<"Signal missing!"<<endl;
        return 1;
    } else if (read.empty()) {
        cout<<"Read missing!"<<endl;
        return 2;
    }

    // process signal: convert string to double array
    const size_t T = count(signal.begin(), signal.end(), ',')+2; // len(sig) + 1
    const size_t N = read.size() - kmerSize + 1 + 1; // N is number of kmers in sequence + 1
    const size_t TN = T*N;

    double* sig = new double[T-1];
    string value;
    stringstream ss(signal);
    int i = 0;
    while(getline(ss, value, ',')) {
        sig[i++] = stof(value);
    }

    // process read N: convert string to int array
    int* kmer_seq = new int[N-1];
    for (size_t n=0; n<N-1; ++n) {
        kmer_seq[n] = kmer2int(read.substr(n, kmerSize), alphabet_size);
    }

    // deallocate memory
    ss.clear();
    signal.erase();
    read.erase();
    value.erase();

    // cerr<<"T: "<<T<<", "<<"N: "<<N<<", "<<"inputsize: "<<TN<<endl;
    // calculate segmentation probabilities, fill forward matrices
    const auto [forM, forE] = logF(sig, kmer_seq, T, N, model);    
    // calculate segmentation probabilities, fill backward matrices
    const auto [backM, backE] = logB(sig, kmer_seq, T, N, model);    
    const double Zf = forE[TN-1];
    const double Zb = backE[0];

    // Numeric error is scaled by input size, Z in forward and backward should match by some numeric error EPSILON
    if (abs(Zf-Zb)/TN > EPSILON || isinf(Zf) || isinf(Zb)) {
        cerr<<"Z values between matrices do not match! Zf: "<<Zf<<", Zb: "<<Zb<<", "<<abs(Zf-Zb)/(TN)<<" > "<<EPSILON<<endl;
        exit(11);
    }
    
    // cerr<<"Zf: "<<Zf<<", Zb: "<<Zb<<", "<<abs(Zf-Zb)/TN<<" <! "<<EPSILON<<endl;

    if (calcZ){
        cout<<Zb<<endl;
    } else {
        
        // train both Transitions and Emissions
        if (train) {
            trainParams(sig, kmer_seq, forM, forE, backM, backE, T, N, model, alphabet_size, numKmers, kmerSize);
            cout<<"Z:"<<Zb<<endl;

        } else {
            const double* LPM = logP(forM, backM, Zb, TN); // log probs for segmentation
            const double* LPE = logP(forE, backE, Zb, TN); // log probs for extension
            const list<string> segString = getBorders(LPM, LPE, T, N, kmerSize);

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
    delete[] kmer_seq;
    return 0;
}
