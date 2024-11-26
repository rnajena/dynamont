// author: Jannes Spangenberg
// e-mail: jannes.spangenberg@uni-jena.de
// github: https://github.com/JannesSP
// website: https://jannessp.github.io

#pragma once

#include <iostream>
#include <iomanip> // std::setprecision
#include <fstream> // file io
#include <sstream> // file io
#include <string>
#include <algorithm> // std::sort, std::stable_sort
#include <numeric>   // std::iota
#include <vector>
#include <array>
#include <unordered_map>
#include <unordered_set>
#include <tuple>
#include <cmath> // exp, pow, log1p, INFINITY
#include <cassert>
#include <cstddef>
#include "argparse.hpp"
#include "utils.hpp"

inline constexpr int NUMMAT = 5;
inline constexpr double SPARSETHRESHOLD = log(0.95); // using paths with top X% of probability per T
inline constexpr double EPSILON = 1e-8;              // chose by eye just to distinguish real errors from numeric errors

extern std::size_t TNK, NK;
extern double ppTNm, ppTNe, ppTKm, ppTKe;
extern int alphabetSize, kmerSize, halfKmerSize, stepSize;
extern bool rna;
extern std::unordered_map<std::string, double> transitions_NTK;

// Asserts double point compatibility at compile time necessary for INFINITY usage
static_assert(std::numeric_limits<double>::is_iec559, "IEEE 754 required");

// ===============================================================
// ===============================================================
// ===================== Scoring calculations ====================
// ===============================================================
// ===============================================================

/**
 * Calculates the Hamming-Distance between two given kmers in their integer base representation
 *
 * @param kmerN
 * @param kmerK
 * @returns log(e^(−2×HD))
 */
inline int scoreHD(const std::size_t kmerN, const std::size_t kmerK)
{
    // TODO: maybe creating an alphabet class with a precalculated distance lookup table would be faster here, depending on how much memory this class/lookup table consumes, only need half the table and maybe a pattern exists that can be used to faster calculate d = -2 * HD(x, y)
    // TODO: lookup table can be combined with model object storing the mean and stddev values

    // no calculation needed if same kmer, should happen often with sparse matrix, accurate basecaller and accurate kmer model
    if (kmerN == kmerK)
        return 0;

    // https://www.codeproject.com/Tips/1274380/Cplusplus11-std-div-Benchmark
    // modulus and integer division is faster than std div
    int acc = 0;
    std::size_t n = kmerN;
    std::size_t k = kmerK;
    for (int i = 0; i < kmerSize; ++i)
    {
        // Perform modulus operation
        const int remN = n % alphabetSize;
        const int remK = k % alphabetSize;
        acc += (remN != remK);
        // Perform integer division
        n /= alphabetSize;
        k /= alphabetSize;
    }
    return -2 * acc; // log(e^(−k×HD)), maybe use k=10 for r9 RNA error rate of roughly 10 %
}

/**
 * Return combined log probability density for a given value and a given normal distribution
 *
 * @param signal point to calculate probability density
 * @param kmerN key for kmer N the model kmer:(mean, stdev) map
 * @param kmerK key for kmer K the model kmer:(mean, stdev) map
 * @param model map containing kmers as keys and (mean, stdev) tuples as values
 * @return log probability density value for x in the given normal distribution
 */
inline double score(const double signal, const std::size_t kmerN, const std::size_t kmerK, const std::tuple<double, double> *model)
{
    // Access elements of the model std::tuple directly to avoid redundant std::tuple creation and overhead
    const auto &[meanN, stdN] = model[kmerN];
    const auto &[meanK, stdK] = model[kmerK];

    // Precompute the scores for the individual kmers and their distance
    double scoreNT = logNormalPdf(signal, meanN, stdN);
    double scoreKT = logNormalPdf(signal, meanK, stdK);
    double scoreNK = scoreHD(kmerN, kmerK);

    return scoreNT + scoreKT + scoreNK;
}

/**
 * Return combined log probability density for a given value and a given normal distribution
 *
 * @param signal point to calculate probability density
 * @param kmerN key for kmer N the model kmer:(mean, stdev) map
 * @param kmerK key for kmer K the model kmer:(mean, stdev) map
 * @param model map containing kmers as keys and (mean, stdev) tuples as values
 * @return log probability density value for x in the given normal distribution
 */
inline double score(const double signal, const std::size_t kmerN, const std::size_t kmerK, const std::tuple<double, double> *model);

/**
 * Computes the posterior probabilities for a hidden Markov model (HMM) using
 * forward and backward probabilities.
 *
 * @param logAPSEI Posterior probabilities stored in an unordered map.
 * @param forAPSEI Forward probabilities stored in an unordered map.
 * @param backAPSEI Backward probabilities stored in an unordered map.
 * @param Z Normalization constant for the forward and backward probabilities.
 * @param allowedKeys Array of allowed indices for processing.
 *
 * @return Posterior probabilities for each state in the HMM.
 */
void logP(std::unordered_map<std::size_t, std::array<dproxy, NUMMAT>> &logAPSEI, const std::unordered_map<std::size_t, std::array<dproxy, NUMMAT>> &forAPSEI, const std::unordered_map<std::size_t, std::array<dproxy, NUMMAT>> &backAPSEI, const double Z, const std::vector<std::size_t> &allowedKeys);

/**
 * Calculate the log probability for a given signal index s of being in one of two states, M or E.
 *
 * @param LP pointer to the output array of log probabilities of size S
 * @param forM pointer to the array of forward probabilities of state M of size S
 * @param backM pointer to the array of backward probabilities of state M of size S
 * @param forE pointer to the array of forward probabilities of state E of size S
 * @param backE pointer to the array of backward probabilities of state E of size S
 * @param S the size of the arrays
 * @param Z the alignment score
 */
void logP(double *LP, const dproxy *forM, const dproxy *backM, const dproxy *forE, const dproxy *backE, const std::size_t S, const double Z);

// ===============================================================
// ===============================================================
// ======================== PREPROCESSING ========================
// ===============================================================
// ===============================================================

/**
 * Computes forward probabilities for a hidden Markov model (HMM) using
 * signal data, k-mer sequence, and transition probabilities.
 *
 * @param sig Pointer to the signal data array.
 * @param kmerSeq Pointer to the k-mer sequence array.
 * @param M Matrix of match probabilities.
 * @param E Matrix of extension probabilities.
 * @param T Length of the ONT raw signal + 1.
 * @param N Length of nucleotide sequence + 1.
 * @param model Array containing tuples of model parameters.
 */
void ppForTN(const double *sig, const int *kmerSeq, dproxy *M, dproxy *E, const std::size_t T, const std::size_t N, const std::tuple<double, double> *model);

/**
 * Computes backward probabilities for a hidden Markov model (HMM) using
 * signal data and transition probabilities.
 *
 * @param sig Pointer to the signal data array.
 * @param M Matrix of match probabilities.
 * @param E Matrix of extension probabilities.
 * @param T Length of the ONT raw signal + 1.
 * @param N Length of nucleotide sequence + 1.
 * @param model Array containing tuples of model parameters.
 */
void ppBackTN(const double *sig, const int *kmerSeq, dproxy *M, dproxy *E, const std::size_t T, const std::size_t N, const std::tuple<double, double> *model);

/**
 * Computes forward probabilities for a hidden Markov model (HMM) using signal
 * data and transition probabilities.
 *
 * @param sig Pointer to the signal data array.
 * @param M Matrix of match probabilities.
 * @param E Matrix of extension probabilities.
 * @param T Length of the ONT raw signal + 1.
 * @param K The number of k-mers.
 * @param model Array containing tuples of model parameters.
 */
void ppForTK(const double *sig, dproxy *M, dproxy *E, const std::size_t T, const std::size_t K, const std::tuple<double, double> *model);

/**
 * Computes backward probabilities for a hidden Markov model (HMM) using
 * signal data and transition probabilities.
 *
 * This function iteratively calculates the backward matrices M and E,
 * which represent match and extend probabilities, respectively. It uses
 * logarithmic probabilities to efficiently process the signal data and
 * k-mer sequences, updating the probabilities based on transition
 * parameters and emission scores.
 *
 * @param sig Pointer to the signal data array.
 * @param M Matrix of match probabilities.
 * @param E Matrix of extension probabilities.
 * @param T Length of the ONT raw signal + 1.
 * @param K The number of k-mers.
 * @param model Array containing tuples of model parameters.
 */
void ppBackTK(const double *sig, dproxy *M, dproxy *E, const std::size_t T, const std::size_t K, const std::tuple<double, double> *model);

/**
 * Compute a sparse representation of the posterior probability matrix for a given signal and k-mer sequence.
 *
 * @param sig Pointer to the ONT raw signal array with pA values.
 * @param kmerSeq Pointer to the k-mer sequence array.
 * @param tnMap Unordered map to store the sparse representation of the posterior probability matrix.
 * @param T Length of the ONT raw signal + 1.
 * @param N Length of nucleotide sequence + 1.
 * @param model Array containing kmers as keys and (mean, stdev) tuples as values.
 */
void preProcTN(const double *sig, const int *kmerSeq, std::unordered_map<std::size_t, std::unordered_set<std::size_t>> &tnMap, const std::size_t T, const std::size_t N, const std::tuple<double, double> *model);

/**
 * Compute a sparse representation of the posterior probability matrix for a given signal and kmer sequence.
 *
 * @param sig Pointer to the ONT raw signal array with pA values.
 * @param tkMap Unordered map to store the sparse representation of the posterior probability matrix.
 * @param T Length of the ONT raw signal + 1.
 * @param K Length of kmer sequence + 1.
 * @param model Array containing kmers as keys and (mean, stdev) tuples as values.
 */
void preProcTK(const double *sig, std::unordered_map<std::size_t, std::unordered_set<std::size_t>> &tkMap, const std::size_t T, const std::size_t K, const std::tuple<double, double> *model);

/**
 * @brief Preprocesses the signal and k-mer sequence to compute allowed k-mer keys.
 *
 * This function performs preprocessing on the ONT raw signal and k-mer sequence
 * to compute a set of allowed k-mer keys. It combines preprocessing results
 * from two partial 2D problems: signal vs nucleotide sequence (TN) and signal
 * vs k-mer sequence (TK). The allowed keys are useful for further analysis
 * such as resquiggles or error corrections.
 *
 * @param sig Pointer to the ONT raw signal array with pA values.
 * @param kmerSeq Pointer to the k-mer sequence array.
 * @param T Length of the ONT raw signal + 1.
 * @param N Length of nucleotide sequence + 1.
 * @param K Number of k-mers.
 * @param model Array containing tuples of model parameters for each k-mer.
 * @return Pointer to the array of allowed k-mer keys.
 */
std::vector<std::size_t> preProcTNK(const double *sig, const int *kmerSeq, const std::size_t T, const std::size_t N, const std::size_t K, const std::tuple<double, double> *model);

// ===============================================================
// ========================== ALGORITHM ==========================
// ===============================================================

/**
 * Computes forward probabilities for a hidden Markov model (HMM) using signal
 * data, k-mer sequences, and transition probabilities.
 *
 * @param sig Pointer to the signal data array.
 * @param kmerSeq Pointer to the k-mer sequence array.
 * @param forAPSEI Forward probabilities stored in an unordered map.
 * @param allowedKeys Array of allowed indices for processing.
 * @param T Length of the ONT raw signal + 1.
 * @param N Length of nucleotide sequence + 1.
 * @param K The number of k-mers.
 * @param model array containing tuples of model parameters.
 */
void logF(const double *sig, const int *kmerSeq, std::unordered_map<std::size_t, std::array<dproxy, NUMMAT>> &forAPSEI, const std::vector<std::size_t> &allowedKeys, const std::size_t T, const std::size_t N, const std::size_t K, const std::tuple<double, double> *model);

/**
 * Computes backward probabilities for a hidden Markov model (HMM) using signal
 * data, k-mer sequences, and transition probabilities.
 *
 * @param sig Pointer to the signal data array.
 * @param kmerSeq Pointer to the k-mer sequence array.
 * @param backAPSEI Backward probabilities stored in an unordered map.
 * @param allowedKeys Array of allowed indices for processing.
 * @param T Length of the ONT raw signal + 1.
 * @param N Length of nucleotide sequence + 1.
 * @param K The number of k-mers.
 * @param model Array containing tuples of model parameters.
 */
void logB(const double *sig, const int *kmerSeq, std::unordered_map<std::size_t, std::array<dproxy, NUMMAT>> &backAPSEI, std::vector<std::size_t> &allowedKeys, const std::size_t T, const std::size_t N, const std::size_t K, const std::tuple<double, double> *model);

/**
 * @brief Computes the segmentation borders for a sequence.
 *
 * This function calculates the segmentation borders based on a list of allowed
 * k-mer indices, forward probabilities, and log probabilities. It updates the
 * segment strings with the highest probability k-mer sequence for the given
 * time and sequence index.
 *
 * @param segString List to store the resulting segment strings.
 * @param logAPSEI Unordered map of log probabilities for each state.
 * @param allowedKeys Array of allowed k-mer indices for processing.
 * @param T Length of the ONT raw signal + 1.
 * @param N Length of nucleotide sequence + 1.
 * @param K Number of k-mers.
 */
void getBorders(std::list<std::string> &segString, std::unordered_map<std::size_t, std::array<dproxy, NUMMAT>> &logAPSEI, const std::vector<std::size_t> &allowedKeys, const std::size_t T, const std::size_t N, const std::size_t K);

/**
 * @brief Backtracing function for A state.
 *
 * This function backtraces the A state.
 * It checks which state was the previous state and calls the corresponding function.
 * If no match is found, it outputs an error message.
 *
 * @param t Current time index.
 * @param n Current segment index.
 * @param k Current kmer index.
 * @param APSEI The 2D array of scores.
 * @param logAPSEI The 2D array of log scores.
 * @param segString The list to store the segment strings.
 * @param K number of kemrs
 * @param segProb The vector to store the segment probabilities.
 *
 * @details
 * This function backtraces the A state.
 * It checks which state was the previous state and calls the corresponding function.
 * If no match is found, it outputs an error message.
 */
void funcA(const std::size_t t, const std::size_t n, const std::size_t k, std::unordered_map<std::size_t, std::array<dproxy, NUMMAT>> &APSEI, std::unordered_map<std::size_t, std::array<dproxy, NUMMAT>> &logAPSEI, std::list<std::string> &segString, const std::size_t K, std::vector<double> &segProb);

/**
 * Backtracing function for E state.
 *
 * @param t Current time index.
 * @param n Current segment index.
 * @param k Current kmer index.
 * @param APSEI The 2D array of scores.
 * @param logAPSEI The 2D array of log scores.
 * @param segString The list to store the segment strings.
 * @param K number of k-mers
 * @param segProb The vector to store the segment probabilities.
 *
 * @details
 * This function backtraces the E state.
 * It checks which state was the previous state and calls the corresponding function.
 * If no match is found, it outputs an error message.
 */
void funcE(const std::size_t t, const std::size_t n, const std::size_t k, std::unordered_map<std::size_t, std::array<dproxy, NUMMAT>> &APSEI, std::unordered_map<std::size_t, std::array<dproxy, NUMMAT>> &logAPSEI, std::list<std::string> &segString, const std::size_t K, std::vector<double> &segProb);

/**
 * @brief Backtraces and constructs segmentation for the P state.
 *
 * This function performs backtracing to construct the segmentation string
 * for the P (probationary) state, using scores from the forward iteration
 * stored in APSEI and logAPSEI. It also updates the probabilities for each
 * segment.
 *
 * @param t Current time step.
 * @param n Current sequence index.
 * @param k Current kmer index.
 * @param APSEI Map containing forward scores.
 * @param logAPSEI Map containing logarithmic forward scores.
 * @param segString List to store the segmentation strings.
 * @param K Number of kmers.
 * @param segProb Vector to store probabilities for each segment.
 */
void funcP(const std::size_t t, const std::size_t n, const std::size_t k, std::unordered_map<std::size_t, std::array<dproxy, NUMMAT>> &APSEI, std::unordered_map<std::size_t, std::array<dproxy, NUMMAT>> &logAPSEI, std::list<std::string> &segString, const std::size_t K, std::vector<double> &segProb);

/**
 * @brief Backtracing for S state
 *
 * @param[in] t time step
 * @param[in] n sequence index
 * @param[in] k kmer index
 * @param[in] APSEI map of scores from forward iteration
 * @param[in] logAPSEI map of log scores from forward iteration
 * @param[out] segString list of strings representing the segmentation
 * @param[in] K number of kmers
 * @param[out] segProb vector of probabilities for each segment
 *
 * This function is used to construct the segmentation string by backtracing
 * the most likely state path from the forward iteration scores.
 */
void funcS(const std::size_t t, const std::size_t n, const std::size_t k, std::unordered_map<std::size_t, std::array<dproxy, NUMMAT>> &APSEI, std::unordered_map<std::size_t, std::array<dproxy, NUMMAT>> &logAPSEI, std::list<std::string> &segString, const std::size_t K, std::vector<double> &segProb);

/**
 * @brief Backtracing function for the Insertion (I) state
 *
 * @param t The current time step
 * @param n The current position
 * @param k The current kmer
 * @param APSEI The forward probabilities table
 * @param logAPSEI The forward log probabilities table
 * @param segString The output segment string
 * @param K number of kmers
 * @param segProb The output segment probabilities
 *
 * @return void
 */
void funcI(const std::size_t t, const std::size_t n, const std::size_t k, std::unordered_map<std::size_t, std::array<dproxy, NUMMAT>> &APSEI, std::unordered_map<std::size_t, std::array<dproxy, NUMMAT>> &logAPSEI, std::list<std::string> &segString, const std::size_t K, std::vector<double> &segProb);

/**
 * Computes transition probabilities for a hidden Markov model (HMM) using
 * forward and backward probabilities, signal data, and k-mer sequences.
 *
 * @param sig Pointer to the signal data array.
 * @param kmerSeq Pointer to the k-mer sequence array.
 * @param forAPSEI Forward probabilities stored in an unordered map.
 * @param backAPSEI Backward probabilities stored in an unordered map.
 * @param allowedKeys Array of allowed indices for processing.
 * @param T Length of the ONT raw signal + 1.
 * @param N Length of nucleotide sequence + 1.
 * @param K The number of k-mers.
 * @param model Array containing tuples of model parameters.
 *
 * @return A tuple containing the transition probabilities for different states
 *         in the HMM.
 */
std::tuple<double, double, double, double, double, double, double, double, double, double, double, double, double, double> trainTransition(const double *sig, const int *kmerSeq, std::unordered_map<std::size_t, std::array<dproxy, NUMMAT>> &forAPSEI, std::unordered_map<std::size_t, std::array<dproxy, NUMMAT>> &backAPSEI, std::vector<std::size_t> &allowedKeys, const std::size_t T, const std::size_t N, const std::size_t K, const std::tuple<double, double> *model);

/**
 * @brief Trains emission parameters using the Baum-Welch algorithm.
 *
 * @param[in] sig Pointer to the ONT raw signal array with pA values.
 * @param[in] logAPSEI Log probability matrix for match states.
 * @param[in] allowedKeys Array of allowed transitions_NTK.
 * @param[in] T Length of the ONT raw signal + 1.
 * @param[in] N Length of nucleotide sequence + 1.
 * @param[in] K Number of kmers in the sequence.
 *
 * Outputs the trained emission parameters to the console.
 */
std::tuple<double *, double *> trainEmission(const double *sig, std::unordered_map<std::size_t, std::array<dproxy, NUMMAT>> &logAPSEI, std::vector<std::size_t> &allowedKeys, const std::size_t T, const std::size_t N, const std::size_t K);

/**
 * @brief Trains transition and emission parameters using the Baum-Welch algorithm.
 *
 * @param sig Pointer to the ONT raw signal array with pA values.
 * @param kmerSeq Pointer to the nucleotide sequence represented by the ONT signal.
 * @param forAPSEI Forward matrix for match states.
 * @param backAPSEI Backward matrix for match states.
 * @param logAPSEI Log probability matrix for match states.
 * @param allowedKeys Array of allowed transitions_NTK.
 * @param T Length of the ONT raw signal + 1.
 * @param N Length of nucleotide sequence + 1.
 * @param K Number of kmers in the sequence.
 * @param model Array containing kmers as keys and (mean, stdev) tuples as values.
 *
 * Outputs the trained parameters to the console.
 */
void trainParams(const double *sig, const int *kmerSeq, std::unordered_map<std::size_t, std::array<dproxy, NUMMAT>> &forAPSEI, std::unordered_map<std::size_t, std::array<dproxy, NUMMAT>> &backAPSEI, std::unordered_map<std::size_t, std::array<dproxy, NUMMAT>> &logAPSEI, std::vector<std::size_t> &allowedKeys, const std::size_t T, const std::size_t N, const std::size_t K, std::tuple<double, double> *model);