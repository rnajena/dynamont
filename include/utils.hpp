// author: Jannes Spangenberg
// e-mail: jannes.spangenberg@uni-jena.de
// github: https://github.com/JannesSP
// website: https://jannessp.github.io

#pragma once

#include <map>
#include <vector>
#include <array>
#include <fstream>   // file io
#include <sstream>   // file io
#include <cmath>     //std::log1p
#include <algorithm> //std::stable_sort, std::reverse
#include <numeric>   //iota
#include <set>
#include <iomanip>    //setprecision
#include <filesystem> // std::filesystem::exists
#include <iostream>

extern const double EPSILON; // chose by eye just to distinguish real errors from numeric errors
extern bool rna;

// default params for NTK
extern const std::map<std::string, double> NTK_rna_r9_transitions;
extern const std::map<std::string, double> NTK_rna_rp4_transitions;
extern const std::map<std::string, double> NTK_dna_r9_transitions;
extern const std::map<std::string, double> NTK_dna_r10_260bps_transitions;
extern const std::map<std::string, double> NTK_dna_r10_400bps_transitions;

// default params for NT
extern const std::map<std::string, double> NT_rna_r9_transitions;
extern const std::map<std::string, double> NT_rna_rp4_transitions;
extern const std::map<std::string, double> NT_dna_r9_transitions;
extern const std::map<std::string, double> NT_dna_r10_260bps_transitions;
extern const std::map<std::string, double> NT_dna_r10_400bps_transitions;

// updatable maps
extern std::map<std::string, double> transitions_NT;
extern std::map<std::string, double> transitions_NTK;

extern const std::map<char, int> BASE2ID;
extern const std::map<int, char> ID2BASE;

// https://stackoverflow.com/questions/72807569/set-default-value-of-unordered-map-if-key-doesnt-exist/72807851#72807851
// workaround to change default double value in map from 0 to -INFINITY
class dproxy
{
    double value_;

public:
    /**
     * A proxy class for setting default value of double to -INFINITY
     *
     * @param value the default value for the double
     */
    dproxy(double value = -INFINITY)
        : value_{value} {}
    operator double() { return value_; }
    operator double() const { return value_; }
};

/**
 * Sorts the column indices of a row-major-indexed double matrix.
 * Complexity is O(C * log(C)), see https://en.cppreference.com/w/cpp/algorithm/stable_sort.
 *
 * @param matrix a double matrix in row major order
 * @param C column size
 * @param t the column to sort for
 *
 * @return std::size_t std::vector with the sorted index of column in descending order
 */
std::vector<std::size_t> columnArgsort(const double *matrix, const std::size_t C, const std::size_t t);

/**
 * C++ version 0.4 std::std::string style "itoa":
 * Contributions from Stuart Lowe, Ray-Yuan Sheu,
 * Rodrigo de Salvo Braz, Luc Gallant, John Maloney
 * and Brian Hunt
 *
 * Converts a decimal to number to a number of base alphabetSize.
 * TODO: Works for base between 2 and 16 (included)
 *
 * Returns kmer in reversed direction!
 *
 * @param value input number in decimal to convert to base
 * @param alphabetSize number of allowed characters in alphabet
 * @param kmerSize length of kmer
 * @param rna true if input is RNA sequence, false if DNA sequence
 * @returns kmer as reversed std::string, should be 5' - 3' direction
 */
std::string itoa(const std::size_t value, const int alphabetSize, const int kmerSize, const bool rna);

/**
 * Converts the kmers of the model file to the integer representation using the BASE2ID map
 *
 * @param s kmer containing nucleotides
 * @param BASE2ID base to id map
 * @param alphabetSize
 * @returns integer representation of the given kmer
 */
int kmer2int(const std::string &s, const int alphabetSize);

/**
 * Reads the normal distribution parameters from a given TSV file,
 * and returns the kmer model and alphabet size.
 *
 * @param file       Path to the TSV file containing kmer parameters (mean, stdev).
 * @param kmerSize   The size of the kmers (length of the kmers in the file).
 * @param rna        True if input is RNA sequence, false if DNA sequence
 * @returns          A std::tuple containing:
 *                   1. An array of tuples, where each std::tuple holds (mean, stdev) for each kmer.
 *                   2. The alphabet size (number of unique nucleotide characters from the kmer set).
 *                   3. The total number of possible kmers (calculated as alphabetSize^kmerSize).
 */
std::tuple<std::tuple<double, double> *, int, std::size_t> readKmerModel(const std::string &file, const std::size_t kmerSize, const bool rna);

/**
 * Reads the normal distribution parameters from a given TSV file,
 * and returns the kmer model and alphabet size.
 *
 * @param file       Path to the TSV file containing kmer parameters (mean, stdev).
 * @param rna        True if input is RNA sequence, false if DNA sequence
 * @returns          A std::tuple containing:
 *                   1. An array of tuples, where each std::tuple holds (mean, stdev) for each kmer.
 *                   2. The alphabet size (number of unique nucleotide characters from the kmer set).
 *                   3. The total number of possible kmers (calculated as alphabetSize^kmerSize).
 *                   4. The number of bases in the model kmers: kmer size
 */
std::tuple<std::tuple<double, double> *, int, std::size_t, std::size_t> readKmerModel(const std::string &file, const bool rna);

// https://en.wikipedia.org/wiki/Log_probability
/**
 * Calculate addition of a+b in log space as efficiently as possible
 * with x + std::log1p(exp(y-x)) : x>y
 *
 * @param a first value
 * @param b second value
 * @return log(exp(a) + exp(b))
 */
inline double logPlus(const double x, const double y)
{
    if (std::isinf(x) && std::isinf(y))
    {
        return x;
    }
    if (x >= y)
    {
        return x + std::log1p(exp(y - x));
    }
    return y + std::log1p(exp(x - y));
}

/**
 * Calculates the integer representation of the successing kmer given the current kmer and the upcoming nucleotide
 * k_i+1 = (k_i mod base^(kmerSize-1)) * base + value(nextNt, base)
 *
 * @param currentKmer current kmer in decimal representation
 * @param nextNt successing nucleotide as a token
 * @param alphabetSize number of accepted characters
 * @param stepSize equals alphabetSize ^ (kmerSize - 1)
 * @return successing Kmer as integer representation in the current base
 */
inline std::size_t successingKmer(const std::size_t currentKmer, const int nextNt, const int stepSize, const int alphabetSize)
{
    return (currentKmer % stepSize) * alphabetSize + nextNt;
}

/**
 * Calculates the integer representation of the precessor kmer given the current kmer and the precessing nucleotide
 * k_i-1 = int(k_i/base) + value(priorNt, base) * base^(kmerSize-1)
 *
 * @param currentKmer current kmer in decimal representation
 * @param priorNt precessing nucleotide as a token
 * @param alphabetSize number of accepted characters
 * @param stepSize equals alphabetSize ^ (kmerSize - 1)
 * @return precessing Kmer as integer representation in the current base
 */
inline std::size_t precessingKmer(const std::size_t currentKmer, const int priorNt, const int stepSize, const int alphabetSize)
{
    return (currentKmer / alphabetSize) + (priorNt * stepSize);
}

// ===============================================================
// ===============================================================
// ===================== Scoring calculations ====================
// ===============================================================
// ===============================================================

// https://ethz.ch/content/dam/ethz/special-interest/mavt/dynamic-systems-n-control/idsc-dam/Lectures/Stochastic-Systems/Statistical_Methods.pdf
/**
 * Calculate log pdf for a given x, mean and standard deviation
 *
 * @param x value
 * @param m mean
 * @param s standard deviation
 * @return probability density at position x for N~(m, s²)
 */
inline double logNormalPdf(const double x, const double m, const double s)
{
    // if (s <= 0) return -INFINITY; // Handle invalid standard deviation

    constexpr double log2Pi = 1.8378770664093453; // Precomputed log(2 * M_PI)

    // Compute s² and its inverse once for repeated calls with similar `s`
    const double s_inv = 1.0 / s;
    const double invVar = s_inv * s_inv;   // 1 / variance = 1 / (s * s)
    const double logVar = 2 * std::log(s); // log(s^2) = 2 * log(s)

    // Difference
    const double diff = x - m;
    const double diffSq = diff * diff;

    // Return log PDF
    return -0.5 * (log2Pi + logVar + diffSq * invVar);
}

/**
 * Return log probability density for a given value and a given normal distribution
 *
 * @param signal point to calculate probability density
 * @param kmer key for the model kmer:(mean, stdev) map
 * @param model map containing kmers as keys and (mean, stdev) tuples as values
 * @return log probability density value for x in the given normal distribution
 */
inline double scoreKmer(const double signal, const std::size_t kmer, const std::tuple<double, double> *model)
{
    // Access elements of the model std::tuple directly to avoid redundant std::tuple creation and overhead
    const auto &[mean, stddev] = model[kmer];
    return logNormalPdf(signal, mean, stddev);
}

/**
 * @brief Updates the transition probabilities by applying logarithmic values.
 *
 * This function checks each transition in the `transitions` map. If a transition value is `-1`, it updates
 * the transition value with the logarithmic value of the corresponding entry from the `defaultVals` map.
 * Otherwise, it applies the logarithm directly to the existing transition value.
 *
 * @param defaultVals A map containing default transition values (std::string keys and double values).
 *                                 These default values are used when a transition value is set to `-1`.
 * @param newVals A map containing current transition values (std::string keys and double values).
 *                    This map is updated with logarithmic values during the function execution.
 */
void updateTransitions(const std::map<std::string, double> &defaultVals, std::map<std::string, double> &newVals);

/**
 * @brief Calculates the median of a given std::vector of double values.
 *
 * This function takes a std::vector of double values as input and returns the median value.
 * If the input std::vector is empty, it throws a std::out_of_range exception.
 *
 * @param vec A reference to the std::vector of double values.
 * @return The median value of the input std::vector.
 *
 * @throws std::out_of_range If the input std::vector is empty.
 *
 * @note The function sorts the input std::vector before calculating the median.
 *
 * @example
 * std::vector<double> values = {1.0, 2.0, 3.0, 4.0, 5.0};
 * double median = calculateMedian(values);
 * # median will be 3.0
 */
double calculateMedian(std::vector<double> &vec);

/**
 * @brief Calculates the median of a given std::vector of double values and formats it as a string.
 *
 * This function takes a std::vector of double values as input, calculates the median value using
 * the `calculateMedian` function, and formats the median as a string with a fixed precision of 5.
 *
 * @param vec A reference to the std::vector of double values.
 * @return A std::string containing the median value formatted with a fixed precision of 5.
 *
 * @throws std::out_of_range If the input std::vector is empty.
 *
 * @example
 * std::vector<double> values = {1.0, 2.0, 3.0, 4.0, 5.0};
 * std::string medianStr = formattedMedian(values);
 * // medianStr will be "3.00000"
 */
std::string formattedMedian(std::vector<double> &vec);

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
void logP(double *LP, const double *FOR, const double *BACK, const double Z, const std::size_t S);

// ===============================================================
// ===============================================================
// ========================= IO Checks ===========================
// ===============================================================
// ===============================================================

void checkModelpath(std::string modelpath);

void checkInput(const std::size_t signalSize, const std::size_t readSize, const std::size_t kmerSize);