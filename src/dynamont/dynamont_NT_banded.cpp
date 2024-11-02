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

// using namespace std;

void funcM(const std::size_t t, const std::size_t n, const dproxy *M, const dproxy *E, const double *LPM, const double *LPE, std::list<std::string> &segString, const std::size_t N, std::vector<double> &segProb, const int kmerSize);
void funcE(const std::size_t t, const std::size_t n, const dproxy *M, const dproxy *E, const double *LPM, const double *LPE, std::list<std::string> &segString, const std::size_t N, std::vector<double> &segProb, const int kmerSize);

inline constexpr double EPSILON = 1e-2;            // chose by eye just to distinguish real errors from numeric errors
inline constexpr std::size_t BANDEDDPWINDOW = 150; // Banded DP approach: Window size around the maximum, adjust as needed
bool rna;

std::unordered_map<std::string, double> transitions = {
    {"m1", -1.0},
    {"e1", -1.0},
    {"e2", -1.0}};

// Asserts doubleing point compatibility at compile time
// necessary for INFINITY usage
static_assert(std::numeric_limits<double>::is_iec559, "IEEE 754 required");

/**
 * Calculate forward matrices using logarithmic values
 *
 * @param sig ONT raw signal with pA values
 * @param seq nucleotide sequence represented by the ONT signal
 * @param T length of the ONT raw signal + 1
 * @param N length of nucleotide sequence + 1
 * @param model map containing kmers as keys and (mean, stdev) tuples as values
 */
void logF(const double *sig, const int *kmer_seq, dproxy *M, dproxy *E, const std::size_t T, const std::size_t N, const std::vector<std::tuple<double, double>> &model)
{
    E[0] = 0;
    // Initialize max indices for the first iteration
    std::size_t max_idx = 0;

    for (std::size_t t = 1; t < T; ++t)
    {
        double max_val = -INFINITY;
        const std::size_t tN = t * N;

        // Determine the range around the max index
        std::size_t start_n = (max_idx > BANDEDDPWINDOW) ? max_idx - BANDEDDPWINDOW : 1;
        std::size_t end_n = (max_idx + BANDEDDPWINDOW < N) ? max_idx + BANDEDDPWINDOW : N - 1;

        for (std::size_t n = start_n; n <= end_n; ++n)
        {                                                                       // Use restricted range
            const double score = scoreKmer(sig[t - 1], kmer_seq[n - 1], model); // Cache scoreKmer for (t-1, n-1)

            // Update matrices
            M[tN + n] = E[(t - 1) * N + (n - 1)] + score + transitions["m1"];
            E[tN + n] = logPlus(M[tN - N + n] + score, E[tN - N + n] + score + transitions["e2"]);

            // Track maximum values for the next iteration
            // Update max indices for the next iteration
            if (M[tN + n] > max_val)
            {
                max_val = M[tN + n];
                max_idx = n;
            }
            if (E[tN + n] > max_val)
            {
                max_val = E[tN + n];
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
void logB(const double *sig, const int *kmer_seq, dproxy *M, dproxy *E, const std::size_t T, const std::size_t N, const std::vector<std::tuple<double, double>> &model)
{
    E[T * N - 1] = 0;
    // Initialize max indices for the first iteration (starting at t=T-1)
    std::size_t max_idx = N - 1;

    // Start backward iteration
    for (std::size_t t = T - 1; t-- > 0;)
    {
        double max_val = -INFINITY;
        const std::size_t tN = t * N;

        // Define the range of n based on the previous max_idx
        std::size_t start_n = (max_idx > BANDEDDPWINDOW) ? max_idx - BANDEDDPWINDOW : 0;
        std::size_t end_n = (max_idx + BANDEDDPWINDOW < N) ? max_idx + BANDEDDPWINDOW + 1 : N;

        for (std::size_t n = end_n; n-- > start_n;)
        { // Use restricted range and iterate backwards
            double mat = -INFINITY;
            double ext = -INFINITY;

            if (n + 1 < N) [[likely]]
            {
                ext = logPlus(ext, M[tN + N + n + 1] + scoreKmer(sig[t], kmer_seq[n], model) + transitions["m1"]);
            }

            if (n > 0) [[likely]]
            {
                const double score = scoreKmer(sig[t], kmer_seq[n - 1], model);
                mat = logPlus(mat, E[tN + N + n] + score);                     // e1 first extend
                ext = logPlus(ext, E[tN + N + n] + score + transitions["e2"]); // e2 extend further
            }

            // Store computed values in matrices
            M[tN + n] = mat;
            E[tN + n] = ext;

            // Update the maximum values for the next iteration
            if (mat > max_val)
            {
                max_val = mat;
                max_idx = n;
            }
            if (ext > max_val)
            {
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
 * @param S size of matrix
 * @param Z alignment score
 * @return matrix containing logarithmic probabilities for segment borders
 */
void logP(double *LP, const dproxy *FOR, const dproxy *BACK, const double Z, const std::size_t S)
{
    for (std::size_t s = 0; s < S; ++s)
    {
        LP[s] = FOR[s] + BACK[s] - Z;
    }
}

/**
 * Calculate the maximum a posteriori path through LP
 *
 */
void getBorders(std::list<std::string> &segString, const double *LPM, const double *LPE, const std::size_t T, const std::size_t N, const int kmerSize)
{
    dproxy *M = new dproxy[T * N];
    dproxy *E = new dproxy[T * N];

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

void funcM(const std::size_t t, const std::size_t n, const dproxy *M, const dproxy *E, const double *LPM, const double *LPE, std::list<std::string> &segString, const std::size_t N, std::vector<double> &segProb, const int kmerSize)
{
    const double score = M[t * N + n];
    const double logScore = LPM[t * N + n];
    segProb.push_back(exp(logScore));
    // if (t<=0 && n<=0){ // Start value
    //     segString.push_front("M0,0," + std::to_string(calculateMedian(segProb)) + ";"); // n-1 because N is 1 larger than the sequences
    //     return;
    // }
    if (t > 0 && n > 0 && score == E[(t - 1) * N + (n - 1)] + logScore)
    {
        segString.push_front("M" + std::to_string(n - 1 + kmerSize / 2) + "," + std::to_string(t - 1) + "," + std::to_string(calculateMedian(segProb)) + ";");
        segProb.clear();
        return funcE(t - 1, n - 1, M, E, LPM, LPE, segString, N, segProb, kmerSize);
    }
}

void funcE(const std::size_t t, const std::size_t n, const dproxy *M, const dproxy *E, const double *LPM, const double *LPE, std::list<std::string> &segString, const std::size_t N, std::vector<double> &segProb, const int kmerSize)
{
    const double score = E[t * N + n];
    const double logScore = LPE[t * N + n];
    segProb.push_back(exp(logScore));
    if (t > 0 && n > 0)
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
 */
std::tuple<double, double, double> trainTransition(const double *sig, const int *kmer_seq, const dproxy *forM, const dproxy *forE, const dproxy *backM, const dproxy *backE, const std::size_t T, const std::size_t N, const std::vector<std::tuple<double, double>> &model)
{
    // Transition parameters
    double newM1 = -INFINITY, newE1 = 0, newE2 = -INFINITY;

    for (std::size_t t = 0; t < T - 1; ++t)
    {
        for (std::size_t n = 0; n <= t && n < N; ++n)
        { // speed up, due to rules no need to look at upper triangle of matrices
            if (n + 1 < N)
            {
                // m1:                 forward(i)    a                   e(i+1)                                  backward(i+1)
                newM1 = logPlus(newM1, forE[t * N + n] + transitions["m1"] + scoreKmer(sig[t], kmer_seq[n], model) + backM[(t + 1) * N + (n + 1)]);
            }

            if (n > 0)
            {
                double score = scoreKmer(sig[t], kmer_seq[n - 1], model);
                // newE1 = logPlus(newE1, forM[t*N+n] + score + backE[(t+1)*N+n]);
                newE2 = logPlus(newE2, forE[t * N + n] + transitions["e2"] + score + backE[(t + 1) * N + n]);
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
 */
std::tuple<double *, double *> trainEmission(const double *sig, const int *kmer_seq, const dproxy *forM, const dproxy *forE, const dproxy *backM, const dproxy *backE, const std::size_t T, const std::size_t N, const std::vector<std::tuple<double, double>> &model, const int numKmers)
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
        counts[kmer_seq[n - 1]]++;
        for (std::size_t t = n; t < T; ++t)
        {                                 // speed up, due to rules no need to look at upper triangle of matrices
            double g = exp(G[t * N + n]); // Cache the exp(G[t * N + n])
            kmers[n] += g * sig[t - 1];   // Accumulate for kmers
            d[n] += g;                    // Accumulate for normalizer
        }
        if (d[n] > 0)
        {
            kmers[n] = kmers[n] / d[n]; // Normalize
        }
    }

    for (std::size_t n = 1; n < N; ++n)
    {
        means[kmer_seq[n - 1]] += kmers[n] / counts[kmer_seq[n - 1]]; // Update means
    }

    // Emission (stdev of kmers)
    std::fill_n(kmers, N, 0);
    for (std::size_t n = 1; n < N; ++n)
    {
        for (std::size_t t = n; t < T; ++t)
        {                                                      // n<=t, speed up, due to rules no need to look at upper triangle of matrices
            double diff = sig[t - 1] - means[kmer_seq[n - 1]]; // Compute difference from mean
            kmers[n] += exp(G[t * N + n]) * diff * diff;       // Accumulate for variance
        }
        if (d[n] > 0)
        {
            kmers[n] = kmers[n] / d[n];
        }
        stdevs[kmer_seq[n - 1]] += kmers[n] / counts[kmer_seq[n - 1]];
    }

    for (std::size_t n = 1; n < N; ++n)
    {
        // transform vars to stdevs
        stdevs[kmer_seq[n - 1]] = sqrt(stdevs[kmer_seq[n - 1]]);
    }

    delete[] G;
    delete[] kmers;
    delete[] counts;
    delete[] d;
    return std::tuple<double *, double *>({means, stdevs});
}

void trainParams(const double *sig, const int *kmer_seq, const dproxy *forM, const dproxy *forE, const dproxy *backM, const dproxy *backE, const std::size_t T, const std::size_t N, const std::vector<std::tuple<double, double>> &model, const int alphabet_size, const int numKmers, const int kmerSize)
{
    const auto [newM, newE1, newE2] = trainTransition(sig, kmer_seq, forM, forE, backM, backE, T, N, model);
    std::cout << "m1:" << newM << ";e1:" << newE1 << ";e2:" << newE2 << std::endl;

    const auto [newMeans, newStdevs] = trainEmission(sig, kmer_seq, forM, forE, backM, backE, T, N, model, numKmers);
    for (int i = 0; i < numKmers; i++)
    {
        if (newStdevs[i] != 0.0)
        {
            std::cout << itoa(i, alphabet_size, kmerSize, rna) << ":" << newMeans[i] << "," << newStdevs[i] << ";";
        }
    }
    std::cout << std::endl;
}

/**
 * Read signal and read from stdin until the TERM_STRING is seen
 */
int main(int argc, char *argv[])
{
    // speedup for I/O
    std::ios_base::sync_with_stdio(0);
    std::cin.tie(0);
    std::cout.tie(0);

    bool train, calcZ, prob;
    std::string pore, modelpath;

    // std::std::cerr precisions
    std::cerr << std::fixed << std::showpoint << std::setprecision(5);
    // std::std::cerr precisions
    std::cout << std::fixed << std::showpoint << std::setprecision(5);

    // Argparser
    argparse::ArgumentParser program("dynamont basic", "0.1");
    // parameters for DP
    program.add_argument("-m1", "--matchscore1").help("Segment transition probability, should be close to (expected number of nucleotdes)/(signal length). Leave at -1 if unset.").default_value(-1.0).scan<'g', double>().store_into(transitions["m1"]);
    program.add_argument("-e1", "--extendscore1").help("First extend probability.").default_value(-1.0).scan<'g', double>().store_into(transitions["e1"]);
    program.add_argument("-e2", "--extendscore2").help("Further extend probability.").default_value(-1.0).scan<'g', double>().store_into(transitions["e2"]);
    program.add_argument("-t", "--train").help("Switch algorithm to transition and emission parameter training mode").default_value(false).implicit_value(true).store_into(train);
    program.add_argument("-z", "--calcZ").help("Switch algorithm to only calculate Z").default_value(false).implicit_value(true).store_into(calcZ);
    program.add_argument("-m", "--model").help("Path to kmer model table").required().store_into(modelpath);
    program.add_argument("-r", "--pore").help("Pore used to sequence the data").default_value("rna_r9").choices("rna_r9", "dna_r9", "rna_rp4", "dna_r10").store_into(pore);
    program.add_argument("-p", "--probabilty").help("Print out the segment border probability").default_value(false).implicit_value(true).store_into(prob);

    try
    {
        program.parse_args(argc, argv);
    }
    catch (const std::runtime_error &err)
    {
        std::cerr << err.what() << std::endl;
        std::cerr << program;
        return 1;
    }

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

    // check that outgoing transitions sum up to 1
    // assert(fabs(exp(transitions["m1"]) + exp(transitions["e2"]) - 1.0) < 1e-2 && "The sum of the outgoing transitions of state E: m1 and e1 must approximately 1.0");

    // std::vector<std::tuple<double, double>> model(numKmers, make_tuple(-INFINITY, -INFINITY));
    assert(!modelpath.empty() && std::filesystem::exists(modelpath) && "Please provide valid a modelpath!");
    const auto result = readKmerModel(modelpath, kmerSize, rna);
    const std::vector<std::tuple<double, double>> model = std::get<0>(result);
    const int alphabet_size = std::get<1>(result);
    const int numKmers = std::get<2>(result);

    // echo 107,107,107.2,108.0,108.9,111.2,105.7,104.3,107.1,105.7,105,105 CAAAAA| src\segment.exe
    // read input, signal and read whitespace separated in single line
    std::string signal, read;
    getline(std::cin, signal);
    getline(std::cin, read);

    if (signal.empty())
    {
        std::cout << "Signal missing!" << std::endl;
        return 1;
    }
    else if (read.empty())
    {
        std::cout << "Read missing!" << std::endl;
        return 2;
    }

    // std::cerr<<"DEBUG 1"<<std::endl;
    // process signal: convert std::string to double array
    const size_t T = count(signal.begin(), signal.end(), ',') + 2; // len(sig) + 1
    const size_t N = read.size() - kmerSize + 1 + 1;               // N is number of kmers in sequence + 1
    assert(T > 1 && "Signal must contain more than one value!");
    assert(N > 1 && "Read must contain more than one value!");
    const size_t TN = T * N;

    double *sig = new double[T - 1];
    std::string value;
    std::stringstream ss(signal);
    int i = 0;
    while (getline(ss, value, ','))
    {
        sig[i++] = stof(value);
    }

    // process read N: convert std::string to int array
    int *kmer_seq = new int[N - 1];
    for (std::size_t n = 0; n < N - 1; ++n)
    {
        kmer_seq[n] = kmer2int(read.substr(n, kmerSize), alphabet_size);
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
    logF(sig, kmer_seq, forM, forE, T, N, model);
    // calculate segmentation probabilities, fill backward matrices
    dproxy *backM = new dproxy[TN];
    dproxy *backE = new dproxy[TN];
    logB(sig, kmer_seq, backM, backE, T, N, model);
    const double Zf = forE[TN - 1];
    const double Zb = backE[0];
    const double Z = std::max(Zf, Zb);

    // Numeric error is scaled by input size, Z in forward and backward should match by some numeric error EPSILON
    if (abs(Zf - Zb) / TN > EPSILON || std::isinf(Zf) || std::isinf(Zb))
    {
        std::cerr << "Z values between matrices do not match! Zf: " << Zf << ", Zb: " << Zb << ", " << abs(Zf - Zb) / (TN) << " > " << EPSILON << ", T: " << T << ", N: " << N << std::endl;
        exit(11);
    }

    // std::cerr<<"Zf: "<<Zf<<", Zb: "<<Zb<<", "<<abs(Zf-Zb)/TN<<" <! "<<EPSILON<<std::endl;

    if (calcZ)
    {
        std::cout << Z << std::endl;
    }
    else
    {

        // train both Transitions and Emissions
        if (train)
        {
            trainParams(sig, kmer_seq, forM, forE, backM, backE, T, N, model, alphabet_size, numKmers, kmerSize);
            std::cout << "Z:" << Z << std::endl;
        }
        else
        {
            double *LPM = new double[TN];
            logP(LPM, forM, backM, Z, TN); // log probs for segmentation
            double *LPE = new double[TN];
            logP(LPE, forE, backE, Z, TN); // log probs for extension
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
    delete[] forM;
    delete[] forE;
    delete[] backM;
    delete[] backE;
    delete[] sig;
    delete[] kmer_seq;
    return 0;
}
