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
#include <vector>
#include <cmath> // exp
#include <cassert>
#include <cstddef>
#include <algorithm>
#include <filesystem> // std::filesystem::exists
#include "argparse.hpp"
#include "utils.hpp"

void funcM(const std::size_t t, const std::size_t n, const dproxy *M, const dproxy *E, const double *LPM, const double *LPE, std::list<std::string> &segString, const std::size_t N, std::vector<double> &segProb, const int kmerSize);
void funcE(const std::size_t t, const std::size_t n, const dproxy *M, const dproxy *E, const double *LPM, const double *LPE, std::list<std::string> &segString, const std::size_t N, std::vector<double> &segProb, const int kmerSize);

constexpr double EPSILON = 1e-8; // chose by eye just to distinguish real errors from numeric errors
bool rna;

std::unordered_map<std::string, double> transitions = {
    {"m1", -1.0},
    {"e1", -1.0},
    {"e2", -1.0}};

// Asserts doubleing point compatibility at compile time
// necessary for INFINITY usage
static_assert(std::numeric_limits<double>::is_iec559, "IEEE 754 required");

/**
 * Computes forward matrices using logarithmic values for signal processing.
 *
 * This function calculates the forward matrices M and E for a given signal
 * and kmer sequence using logarithmic probabilities. It iterates over the
 * time steps and latent states, updating the matrices based on the transition
 * probabilities and emission scores.
 *
 * @param sig Pointer to the ONT raw signal array with pA values.
 * @param kmerSeq Pointer to the nucleotide sequence represented by the ONT signal.
 * @param M Matrix to store the match scores.
 * @param E Matrix to store the extend scores.
 * @param T Length of the ONT raw signal + 1.
 * @param N Length of nucleotide sequence + 1.
 * @param model Vector containing kmers as keys and (mean, stdev) tuples as values.
 */
void logF(const double *sig, const int *kmerSeq, dproxy *M, dproxy *E, const std::size_t T, const std::size_t N, const std::tuple<double, double> *model)
{
    const double m1 = transitions.at("m1");
    const double e2 = transitions.at("e2");
    E[0] = 0;
    for (std::size_t t = 1; t < T; ++t)
    {
        const std::size_t tN = t * N;
        for (std::size_t n = 1; n <= t && n < N; ++n)
        {                                                                      // speed up, due to rules no need to look at upper triangle of matrices
            const double score = scoreKmer(sig[t - 1], kmerSeq[n - 1], model); // Cache scoreKmer for (t-1, n-1)
            M[tN + n] = E[tN - N + (n - 1)] + score + m1;
            E[tN + n] = logPlus(M[tN - N + n] + score, E[tN - N + n] + score + e2);
        }
    }
}

/**
 * Computes backward matrices using logarithmic values for signal processing.
 *
 * This function calculates the backward matrices M and E for a given signal
 * and kmer sequence using logarithmic probabilities. It iterates over the
 * time steps and latent states, updating the matrices based on the transition
 * probabilities and emission scores.
 *
 * @param sig Pointer to the ONT raw signal array with pA values.
 * @param kmerSeq Pointer to the nucleotide sequence represented by the ONT signal.
 * @param M Matrix to store the match scores.
 * @param E Matrix to store the extend scores.
 * @param T Length of the ONT raw signal + 1.
 * @param N Length of nucleotide sequence + 1.
 * @param model Vector containing kmers as keys and (mean, stdev) tuples as values.
 */
void logB(const double *sig, const int *kmerSeq, dproxy *M, dproxy *E, const std::size_t T, const std::size_t N, const std::tuple<double, double> *model)
{
    const double m1 = transitions.at("m1");
    const double e2 = transitions.at("e2");
    E[T * N - 1] = 0;
    for (std::size_t t = T - 1; t-- > 0;)
    {
        const std::size_t tN = t * N;
        for (std::size_t n = std::min(N, t + 1); n-- > 0;)
        { // speed up, due to rules no need to look at upper triangle of matrices
            double mat = -INFINITY, ext = -INFINITY;
            if (n + 1 < N) [[likely]]
            {
                ext = logPlus(ext, M[tN + N + n + 1] + scoreKmer(sig[t], kmerSeq[n], model) + m1);
            }
            if (n > 0) [[likely]]
            {
                const double score = scoreKmer(sig[t], kmerSeq[n - 1], model);
                mat = logPlus(mat, E[tN + N + n] + score);      // transition probability is always 1, e1 first extend
                ext = logPlus(ext, E[tN + N + n] + score + e2); // e2 extend further
            }
            M[tN + n] = mat;
            E[tN + n] = ext;
        }
    }
}

/**
 * Calculate the logarithmic probability matrix
 *
 * @param LP Matrix to store the logarithmic probabilities.
 * @param FOR Matrix containing forward-values for segment borders.
 * @param BACK Matrix containing backward-values for extending segment.
 * @param Z Alignment score.
 * @param S Size of matrix.
 *
 * This function calculates the logarithmic probability matrix for a given forward and backward
 * matrix by adding the values of the forward and backward matrix and subtracting the alignment score.
 */
void logP(double *LP, const dproxy *FOR, const dproxy *BACK, const double Z, const std::size_t S)
{
    for (std::size_t s = 0; s < S; ++s)
    {
        LP[s] = FOR[s] + BACK[s] - Z;
    }
}

/**
 * Calculate the maximum a posteriori path through logarithmic probability matrices.
 *
 * This function computes the most probable path through given matrices of match and extend scores,
 * storing segment information in a list of strings. It utilizes dynamic programming to fill matrices
 * for match (M) and extend (E) states and performs backtracking to extract the path.
 *
 * @param segString List to store segment information as strings.
 * @param LPM Matrix containing logarithmic probabilities for match state.
 * @param LPE Matrix containing logarithmic probabilities for extend state.
 * @param T Number of time steps, representing the length of the ONT raw signal + 1.
 * @param N Number of latent states, representing the length of the nucleotide sequence + 1.
 * @param kmerSize Size of k-mer used in the analysis.
 */
void getBorders(std::list<std::string> &segString, const double *LPM, const double *LPE, const std::size_t T, const std::size_t N, const int kmerSize)
{
    const std::size_t TN = T * N;
    dproxy *M = new dproxy[TN];
    dproxy *E = new dproxy[TN];
    E[0] = 0;
    for (std::size_t t = 1; t < T; ++t)
    {
        for (std::size_t n = 1; n <= t && n < N; ++n)
        {                                                                                                      // speed up, due to rules no need to look at upper triangle of matrices
            M[t * N + n] = E[(t - 1) * N + (n - 1)] + LPM[t * N + n];                                          // m1
            E[t * N + n] = std::max(M[(t - 1) * N + n] + LPE[t * N + n], E[(t - 1) * N + n] + LPE[t * N + n]); // e1, e2
        }
    }
    std::vector<double> segProb;
    funcE(T - 1, N - 1, M, E, LPM, LPE, segString, N, segProb, kmerSize);
    delete[] M;
    delete[] E;
}

/**
 * Backtracing function for M state
 *
 * @param t current time step
 * @param n current position in sequence
 * @param M matrix containing match scores
 * @param E matrix containing extend scores
 * @param LPM matrix containing logarithmic probabilities for M state
 * @param LPE matrix containing logarithmic probabilities for E state
 * @param segString list to store segment information
 * @param N size of sequence
 * @param segProb vector to store segment probabilities
 * @param kmerSize size of k-mer
 */
void funcM(const std::size_t t, const std::size_t n, const dproxy *M, const dproxy *E, const double *LPM, const double *LPE, std::list<std::string> &segString, const std::size_t N, std::vector<double> &segProb, const int kmerSize)
{
    const double score = M[t * N + n];
    const double logScore = LPM[t * N + n];
    segProb.push_back(exp(logScore));
    // if (t<=1 && n<=1) [[unlikely]] { // Start value
    //     segString.push_front("M0,0," + formattedMedian(segProb) + ";"); // n-1 because N is 1 larger than the sequences
    //     return;
    // }
    if (t > 0 && n > 0 && score == E[(t - 1) * N + (n - 1)] + logScore) [[likely]]
    {
        segString.push_front("M" + std::to_string(n - 1 + kmerSize / 2) + "," + std::to_string(t - 1) + "," + formattedMedian(segProb) + ";");
        segProb.clear();
        return funcE(t - 1, n - 1, M, E, LPM, LPE, segString, N, segProb, kmerSize);
    }
}

/**
 * Backtracing function for E state
 *
 * @param t current time step
 * @param n current position in sequence
 * @param M matrix containing match scores
 * @param E matrix containing extend scores
 * @param LPM matrix containing logarithmic probabilities for M state
 * @param LPE matrix containing logarithmic probabilities for E state
 * @param segString list to store segment information
 * @param N size of sequence
 * @param segProb vector to store segment probabilities
 * @param kmerSize size of k-mer
 */
void funcE(const std::size_t t, const std::size_t n, const dproxy *M, const dproxy *E, const double *LPM, const double *LPE, std::list<std::string> &segString, const std::size_t N, std::vector<double> &segProb, const int kmerSize)
{
    const double score = E[t * N + n];
    const double logScore = LPE[t * N + n];
    segProb.push_back(exp(logScore));
    if (t > 0 && n > 0) [[likely]]
    {
        if (score == M[(t - 1) * N + n] + logScore)
        {
            return funcM(t - 1, n, M, E, LPM, LPE, segString, N, segProb, kmerSize);
        }
        if (score == E[(t - 1) * N + n] + logScore)
        {
            return funcE(t - 1, n, M, E, LPM, LPE, segString, N, segProb, kmerSize);
        }
    }
}

/**
 * Train transition parameter with baum welch algorithm
 *
 * @param sig ONT raw signal with pA values
 * @param kmerSeq nucleotide sequence represented by the ONT signal
 * @param forM matrix containing forward probabilities for M state
 * @param forE matrix containing forward probabilities for E state
 * @param backM matrix containing backward probabilities for M state
 * @param backE matrix containing backward probabilities for E state
 * @param T length of the ONT raw signal
 * @param N length of nucleotide sequence
 * @param model map containing kmers as keys and (mean, stdev) tuples as values
 * @return tuple containing new transition probabilities
 */
std::tuple<double, double, double> trainTransition(const double *sig, const int *kmerSeq, const dproxy *forM, const dproxy *forE, const dproxy *backM, const dproxy *backE, const std::size_t T, const std::size_t N, const std::tuple<double, double> *model)
{
    // Transition parameters
    double newM1 = -INFINITY, newE1 = 0, newE2 = -INFINITY;

    for (std::size_t t = 0; t < T - 1; ++t)
    {
        for (std::size_t n = 0; n <= t && n < N; ++n)
        { // speed up, due to rules no need to look at upper triangle of matrices
            if (n + 1 < N) [[likely]]
            {
                // m1:                 forward(i)        a                      e(i+1)                                 backward(i+1)
                newM1 = logPlus(newM1, forE[t * N + n] + transitions.at("m1") + scoreKmer(sig[t], kmerSeq[n], model) + backM[(t + 1) * N + (n + 1)]);
            }

            if (n > 0) [[likely]]
            {
                double score = scoreKmer(sig[t], kmerSeq[n - 1], model);
                // newE1 = logPlus(newE1, forM[t*N+n] + transitions.at("e1") + score + backE[(t+1)*N+n]);
                newE2 = logPlus(newE2, forE[t * N + n] + transitions.at("e2") + score + backE[(t + 1) * N + n]);
            }
        }
    }
    // average over the number of transitions
    // double Am = newE1;
    // newE1 = newE1 - Am;
    const double Ae = logPlus(newE2, newM1);
    if (!std::isinf(Ae))
    {
        newM1 = newM1 - Ae;
        newE2 = newE2 - Ae;
    }
    return std::make_tuple(
        exp(newM1),
        exp(newE1),
        exp(newE2));
}

/**
 * Train emission parameter with baum welch algorithm
 * @param sig signal
 * @param kmerSeq sequence of kmers
 * @param forM forward probabilities for state M
 * @param forE forward probabilities for state E
 * @param backM backward probabilities for state M
 * @param backE backward probabilities for state E
 * @param T number of time steps
 * @param N number of latent states (kmers)
 * @param model emission model
 * @param numKmers number of kmers
 * @return tuple of emission parameter means and stdevs
 */
std::tuple<double *, double *> trainEmission(const double *sig, const int *kmerSeq, const dproxy *forM, const dproxy *forE, const dproxy *backM, const dproxy *backE, const std::size_t T, const std::size_t N, const std::tuple<double, double> *model, const int numKmers)
{

    // TODO: calculate this with LP matrices
    // see notes

    // Emission
    // https://courses.grainger.illinois.edu/ece417/fa2021/lectures/lec15.pdf
    // https://f.hubspotusercontent40.net/hubfs/8111846/Unicon_October2020/pdf/bilmes-em-algorithm.pdf
    // gamma_t(i) is the probability of being in state i at time t
    // gamma for state M - expected number of transitions of M at given time (T) for all latent states (kmers)
    // Initialize memory
    const std::size_t TN = T * N;
    double *G = new double[TN];
    double *kmers = new double[N];
    double *d = new double[N];
    double *means = new double[numKmers];
    double *stdevs = new double[numKmers];
    int *counts = new int[numKmers];

    // init everything with zero
    for (std::size_t i = 0; i < N; i++)
    {
        kmers[i] = 0.0;
        d[i] = 0.0;
    }
    for (int i = 0; i < numKmers; i++)
    {
        means[i] = 0.0;
        stdevs[i] = 0.0;
        counts[i] = 0;
    }

    // initialize G
    for (std::size_t t = 0; t < T; ++t)
    {
        // calibrate with the sum of transitions
        double s = -INFINITY;
        for (std::size_t n = 0; n <= t && n < N; ++n)
        { // speed up, due to rules no need to look at upper triangle of matrices
            G[t * N + n] = logPlus(forM[t * N + n] + backM[t * N + n], forE[t * N + n] + backE[t * N + n]);
            s = logPlus(s, G[t * N + n]);
        }
        for (std::size_t n = 0; n <= t && n < N; ++n)
        {
            if (!std::isinf(s))
            {
                G[t * N + n] -= s;
            }
        }
    }

    for (std::size_t n = 1; n < N; ++n)
    {
        counts[kmerSeq[n - 1]]++;
        for (std::size_t t = n; t < T; ++t)
        {                                 // speed up, due to rules no need to look at upper triangle of matrices
            double g = exp(G[t * N + n]); // Cache the exp(G[t * N + n])
            kmers[n] += g * sig[t - 1];   // Accumulate for kmers
            d[n] += g;                    // Accumulate for normalizer
        }
        if (d[n] != 0)
        {
            kmers[n] = kmers[n] / d[n]; // Normalize
        }
    }

    for (std::size_t n = 1; n < N; ++n)
    {
        means[kmerSeq[n - 1]] += kmers[n] / counts[kmerSeq[n - 1]]; // Update means
    }

    // Emission (stdev of kmers)
    std::fill_n(kmers, N, 0);
    for (std::size_t n = 1; n < N; ++n)
    {
        for (std::size_t t = n; t < T; ++t)
        {                                                     // n<=t, speed up, due to rules no need to look at upper triangle of matrices
            double diff = sig[t - 1] - means[kmerSeq[n - 1]]; // Compute difference from mean
            kmers[n] += exp(G[t * N + n]) * diff * diff;      // Accumulate for variance
        }
        if (d[n] > 0)
        {
            kmers[n] = kmers[n] / d[n];
        }
        stdevs[kmerSeq[n - 1]] += kmers[n] / counts[kmerSeq[n - 1]];
    }

    for (std::size_t n = 1; n < N; ++n)
    {
        // transform vars to stdevs
        stdevs[kmerSeq[n - 1]] = sqrt(stdevs[kmerSeq[n - 1]]);
    }

    delete[] G;
    delete[] kmers;
    delete[] counts;
    delete[] d;
    return std::tuple<double *, double *>({means, stdevs});
}

/**
 * Trains transition and emission parameters using the Baum-Welch algorithm.
 *
 * @param sig Pointer to the ONT raw signal array with pA values.
 * @param kmerSeq Pointer to the nucleotide sequence represented by the ONT signal.
 * @param forM Pointer to the forward matrix for match states.
 * @param forE Pointer to the forward matrix for extend states.
 * @param backM Pointer to the backward matrix for match states.
 * @param backE Pointer to the backward matrix for extend states.
 * @param T Length of the ONT raw signal + 1.
 * @param N Length of nucleotide sequence + 1.
 * @param model Reference to a vector containing kmers as keys and (mean, stdev) tuples as values.
 * @param alphabetSize The size of the nucleotide alphabet.
 * @param numKmers The number of kmers in the sequence.
 * @param kmerSize The size of each kmer.
 */
void trainParams(const double *sig, const int *kmerSeq, const dproxy *forM, const dproxy *forE, const dproxy *backM, const dproxy *backE, const std::size_t T, const std::size_t N, std::tuple<double, double> *model, const int alphabetSize, const int numKmers, const int kmerSize)
{
    const auto [newM, newE1, newE2] = trainTransition(sig, kmerSeq, forM, forE, backM, backE, T, N, model);
    std::cout << "m1:" << newM << ";e1:" << newE1 << ";e2:" << newE2 << std::endl;
    const auto [newMeans, newStdevs] = trainEmission(sig, kmerSeq, forM, forE, backM, backE, T, N, model, numKmers);
    for (int i = 0; i < numKmers; i++)
    {
        if (newStdevs[i] != 0.0) [[likely]]
        {
            std::cout << itoa(i, alphabetSize, kmerSize, rna) << ":" << newMeans[i] << "," << newStdevs[i] << ";";
        }
    }
    std::cout << std::endl;
    delete[] newMeans;
    delete[] newStdevs;
}

/**
 * Main entry point for the dynamont command line tool.
 *
 * @param argc Number of arguments passed to the program.
 * @param argv Array of strings containing the arguments passed to the program.
 *
 * @return 0 if the program succeeded, 1 if there was an error parsing the arguments, 2 if the input was invalid, 3 if the pore was not supported, 11 if the Z values did not match between forward and backward matrices.
 *
 * @details
 * This function parses the command line arguments, loads the model, reads the input, calculates the segmentation probabilities and prints the result.
 */
#ifndef UNIT_TESTING
int main(int argc, char *argv[])
{
    // speedup for I/O
    std::ios_base::sync_with_stdio(0);
    std::cin.tie(0);
    std::cout.tie(0);

    bool train, calcZ, prob;
    std::string pore, modelpath;

    std::cerr << std::fixed << std::showpoint << std::setprecision(5);
    std::cout << std::fixed << std::showpoint << std::setprecision(5);

    // Argparser
    argparse::ArgumentParser program("dynamont basic", "0.1");
    // parameters for DP
    program.add_argument("-m", "--model").help("Path to kmer model table").required().store_into(modelpath);
    program.add_argument("-r", "--pore").help("Pore used to sequence the data").required().choices("rna_r9", "dna_r9", "rna_rp4", "dna_r10_260bps", "dna_r10_400bps").store_into(pore);

    program.add_argument("-m1", "--matchscore1").help("Segment transition probability, should be close to (expected number of nucleotdes)/(signal length). Leave at -1 if unset.").default_value(-1.0).scan<'g', double>().store_into(transitions["m1"]);
    program.add_argument("-e1", "--extendscore1").help("First extend probability.").default_value(-1.0).scan<'g', double>().store_into(transitions["e1"]);
    program.add_argument("-e2", "--extendscore2").help("Further extend probability.").default_value(-1.0).scan<'g', double>().store_into(transitions["e2"]);
    program.add_argument("-t", "--train").help("Switch algorithm to transition and emission parameter training mode").default_value(false).implicit_value(true).store_into(train);
    program.add_argument("-z", "--calcZ").help("Switch algorithm to only calculate Z").default_value(false).implicit_value(true).store_into(calcZ);
    program.add_argument("-p", "--probabilty").help("Print out the segment border probability").default_value(false).implicit_value(true).store_into(prob);

    program.parse_args(argc, argv);

    int kmerSize = 0;
    // load default and set parameters
    if (pore == "rna_r9")
    {
        kmerSize = 5;
        rna = true;
        // taken from the trained NT version of dynamont
        updateTransitions(NT_rna_r9_transitions, transitions);
    }
    else if (pore == "dna_r9")
    {
        kmerSize = 5;
        rna = false;
        // taken from the trained NT version of dynamont
        updateTransitions(NT_dna_r9_transitions, transitions);
    }
    else if (pore == "rna_rp4")
    {
        kmerSize = 9;
        rna = true;
        // taken from the trained NT version of dynamont
        updateTransitions(NT_rna_rp4_transitions, transitions);
    }
    else if (pore == "dna_r10_260bps")
    {
        kmerSize = 9;
        rna = false;
        updateTransitions(NT_dna_r10_260bps_transitions, transitions);
    }
    else if (pore == "dna_r10_400bps")
    {
        kmerSize = 9;
        rna = false;
        updateTransitions(NT_dna_r10_400bps_transitions, transitions);
    }

    assert(!modelpath.empty() && std::filesystem::exists(modelpath) && "Please provide a valid modelpath!");
    auto result = readKmerModel(modelpath, kmerSize, rna);
    std::tuple<double, double> *model = std::get<0>(result);
    const int alphabetSize = std::get<1>(result);
    const int numKmers = std::get<2>(result);

    // example
    // 107,107,107.2,108.0,108.9,111.2,105.7,104.3,107.1,105.7,105,105
    // CAAAAA
    // read input, signal and read whitespace separated in single line
    std::string signal, read;
    getline(std::cin, signal);
    getline(std::cin, read);

    if (signal.empty())
    {
        std::cout << "Signal missing!" << std::endl;
        exit(4);
    }
    else if (read.empty())
    {
        std::cout << "Read missing!" << std::endl;
        exit(5);
    }

    // process signal: convert std::string to double std::array
    const std::size_t T = count(signal.begin(), signal.end(), ',') + 2; // len(sig) + 1
    const std::size_t N = read.size() - kmerSize + 1 + 1;               // N is number of kmers in sequence + 1
    assert(T > 1 && "Signal must contain more than one value!");
    assert(N > 1 && "Read must contain more than one value!");
    const std::size_t TN = T * N;

    double *sig = new double[T - 1];
    std::string value;
    std::stringstream ss(signal);
    int i = 0;
    while (getline(ss, value, ','))
    {
        sig[i++] = stof(value);
    }

    // process read N: convert std::string to int std::array
    int *kmerSeq = new int[N - 1];
    for (std::size_t n = 0; n < N - 1; ++n)
    {
        kmerSeq[n] = kmer2int(read.substr(n, kmerSize), alphabetSize);
    }

    // deallocate memory
    ss.clear();
    signal.erase();
    read.erase();
    value.erase();

    // std::cerr<<"T: "<<T<<", "<<"N: "<<N<<", "<<"inputsize: "<<TN<<std::endl;
    // calculate segmentation probabilities, fill forward matrices
    dproxy *forM = new dproxy[TN];
    dproxy *forE = new dproxy[TN];
    logF(sig, kmerSeq, forM, forE, T, N, model);
    // calculate segmentation probabilities, fill backward matrices
    dproxy *backM = new dproxy[TN];
    dproxy *backE = new dproxy[TN];
    logB(sig, kmerSeq, backM, backE, T, N, model);
    const double Zf = forE[TN - 1];
    const double Zb = backE[0];

    // Numeric error is scaled by input size, Z in forward and backward should match by some numeric error EPSILON
    if (abs(Zf - Zb) / TN > EPSILON || std::isinf(Zf) || std::isinf(Zb))
    {
        std::cerr << "Z values between matrices do not match! Zf: " << Zf << ", Zb: " << Zb << ", " << abs(Zf - Zb) / (TN) << " > " << EPSILON << std::endl;
        delete[] sig;
        delete[] kmerSeq;
        delete[] forM;
        delete[] forE;
        delete[] backM;
        delete[] backE;
        delete[] model;
        exit(3);
    }

    // std::cerr<<"Zf: "<<Zf<<", Zb: "<<Zb<<", "<<abs(Zf-Zb)/TN<<" <! "<<EPSILON<<std::endl;

    if (calcZ)
    {
        std::cout << Zb << std::endl;
    }
    else
    {

        // train both Transitions and Emissions
        if (train)
        {
            trainParams(sig, kmerSeq, forM, forE, backM, backE, T, N, model, alphabetSize, numKmers, kmerSize);
            std::cout << "Z:" << Zb << std::endl;
        }
        else
        {
            double *LPM = new double[TN];
            logP(LPM, forM, backM, Zb, TN); // log probs for segmentation
            double *LPE = new double[TN];
            logP(LPE, forE, backE, Zb, TN); // log probs for extension
            std::list<std::string> segString;
            getBorders(segString, LPM, LPE, T, N, kmerSize);

            for (const auto &seg : segString)
            {
                std::cout << seg;
            }
            std::cout << std::endl;

            // calculate sum of segment probabilities
            if (prob)
            {
                double sum;
                for (std::size_t t = 0; t < T; ++t)
                {
                    sum = -INFINITY;
                    for (std::size_t n = 0; n < N; ++n)
                    {
                        sum = logPlus(sum, LPM[t * N + n]);
                    }
                    std::cout << sum << ",";
                }
                std::cout << std::endl;
            }
            // Clean up
            delete[] LPM;
            delete[] LPE;
        }
    }
    delete[] sig;
    delete[] kmerSeq;
    delete[] forM;
    delete[] forE;
    delete[] backM;
    delete[] backE;
    delete[] model;
}
#endif
