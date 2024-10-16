// ===============================================================
// ===============================================================
// =========================== Utility ===========================
// ===============================================================
// ===============================================================

#pragma once

#include <unordered_map>
#include <vector>
#include <fstream> // file io
#include <sstream> // file io
#include <cmath> //log1p
#include <algorithm> //stable_sort
#include <numeric> //iota

using namespace std;

// default params for NTK
const unordered_map<string, double> NTK_rna_r9_transitions = {
    {"a1", 0.030570760232128354},
    {"a2", 0.6101330423858251},
    {"p1", 0.0383208584222504},
    {"p2", 0.08685413259505585},
    {"p3", 0.015545502978626788},
    {"s1", 0.0010061921517067177},
    {"s2", 0.001534865889565432},
    {"s3", 0.13149698483529376},
    {"e1", 1.0},
    {"e2", 0.9989938127816652},
    {"e3", 0.9616791335752949},
    {"e4", 0.8801771554834845},
    {"i1", 0.0008630902783022252},
    {"i2", 0.24282445342753348}
};
const unordered_map<string, double> NTK_rna_rp4_transitions = {
    {"a1", 1.0},
    {"a2", 1.0},
    {"p1", 1.0},
    {"p2", 1.0},
    {"p3", 1.0},
    {"s1", 1.0},
    {"s2", 1.0},
    {"s3", 1.0},
    {"e1", 1.0},
    {"e2", 1.0},
    {"e3", 1.0},
    {"e4", 1.0},
    {"i1", 1.0},
    {"i2", 1.0}
};
const unordered_map<string, double> NTK_dna_r9_transitions = {
    {"a1", 1.0},
    {"a2", 1.0},
    {"p1", 1.0},
    {"p2", 1.0},
    {"p3", 1.0},
    {"s1", 1.0},
    {"s2", 1.0},
    {"s3", 1.0},
    {"e1", 1.0},
    {"e2", 1.0},
    {"e3", 1.0},
    {"e4", 1.0},
    {"i1", 1.0},
    {"i2", 1.0}
};
const unordered_map<string, double> NTK_dna_r10_260bps_transitions = {
    {"a1", 1.0},
    {"a2", 1.0},
    {"p1", 1.0},
    {"p2", 1.0},
    {"p3", 1.0},
    {"s1", 1.0},
    {"s2", 1.0},
    {"s3", 1.0},
    {"e1", 1.0},
    {"e2", 1.0},
    {"e3", 1.0},
    {"e4", 1.0},
    {"i1", 1.0},
    {"i2", 1.0}
};
const unordered_map<string, double> NTK_dna_r10_400bps_transitions = {
    {"a1", 1.0},
    {"a2", 1.0},
    {"p1", 1.0},
    {"p2", 1.0},
    {"p3", 1.0},
    {"s1", 1.0},
    {"s2", 1.0},
    {"s3", 1.0},
    {"e1", 1.0},
    {"e2", 1.0},
    {"e3", 1.0},
    {"e4", 1.0},
    {"i1", 1.0},
    {"i2", 1.0}
};

// default params for NT
const unordered_map<string, double> NT_rna_r9_transitions = {
    {"m1", 0.030387},
    {"e1", 1.0},
    {"e2", 0.969613}
};
const unordered_map<string, double> NT_rna_rp4_transitions = {
    {"m1", 1.0},
    {"e1", 1.0},
    {"e2", 1.0}
};
const unordered_map<string, double> NT_dna_r9_transitions = {
    {"m1", 1.0},
    {"e1", 1.0},
    {"e2", 1.0}
};
const unordered_map<string, double> NT_dna_r10_260bps_transitions = {
    {"m1", 1.0},
    {"e1", 1.0},
    {"e2", 1.0}
};
const unordered_map<string, double> NT_dna_r10_400bps_transitions = {
    {"m1", 1.0},
    {"e1", 1.0},
    {"e2", 1.0}
};

const unordered_map<char, int> BASE2ID = {
    {'A', 0},
    {'a', 0},
    {'C', 1},
    {'c', 1},
    {'G', 2},
    {'g', 2},
    {'T', 3},
    {'t', 3},
    {'U', 3},
    {'u', 3},
    {'N', 4},
    {'n', 4}
}; // Nucleotide : Token map
const unordered_map<int, char> ID2BASE = {
    {'0', 'A'},
    {'1', 'C'},
    {'2', 'G'},
    {'3', 'T'},
    {'4', 'N'}
}; // Token : Nucleotide map

/**
 * Sorts the column indices of a row-major-indexed double matrix.
 * Complexity is O(C * log(C)), see https://en.cppreference.com/w/cpp/algorithm/stable_sort.
 * 
 * @param matrix a double matrix in row major order
 * @param C column size
 * @param t the column to sort for
 * 
 * @return size_t vector with the sorted index of column in descending order
 */
vector<size_t> column_argsort(const double* matrix, const size_t C, const size_t t);

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
 * @param ALPHABET_SIZE number of allowed characters in alphabet
 * @param kmerSize length of kmer
 * @returns kmer as reversed string, should be 5' - 3' direction
*/
string itoa(const size_t value, const int ALPHABET_SIZE, const int kmerSize);

/**
 * Converts the kmers of the model file to the integer representation using the BASE2ID map
 *
 * @param s kmer containing nucleotides 
 * @param BASE2ID base to id map 
 * @param ALPHABET_SIZE
 * @returns integer representation of the given kmer
 */
int kmer2int(const string &s, const int ALPHABET_SIZE);

/**
 * Read the normal distribution parameters from a given TSV file
 *
 * @param file path to the TSV file containing the parameters
 * @param model kmer model to fill
 * @param ALPHABET_SIZE
 */
void readKmerModel(const string &file, vector<tuple<double, double>> &model, const int ALPHABET_SIZE);

// https://en.wikipedia.org/wiki/Log_probability
/**
 * Calculate addition of a+b in log space as efficiently as possible
 * with x + log1p(exp(y-x)) : x>y
 * 
 * @param a first value
 * @param b second value
 * @return log(exp(a) + exp(b))
 */
inline double logPlus(const double x, const double y) {
    if (isinf(x) && isinf(y)) {
        return x;
    }
    if (x>=y){
        return x + log1p(exp(y-x));
    }
    return y + log1p(exp(x-y));
}

/**
 * Calculates the integer representation of the successing kmer given the current kmer and the upcoming nucleotide
 * k_i+1 = (k_i mod base^(kmerSize-1)) * base + value(nextNt, base)
 * 
 * @param currentKmer current kmer in decimal representation
 * @param nextNt successing nucleotide as a token
 * @param ALPHABET_SIZE number of accepted characters
 * @param stepSize equals ALPHABET_SIZE ^ (kmerSize - 1)
 * @return successing Kmer as integer representation in the current base
 */
inline size_t successingKmer(const size_t currentKmer, const int nextNt, const int stepSize, const int ALPHABET_SIZE) {
    return (currentKmer % stepSize) * ALPHABET_SIZE + nextNt;
}

/**
 * Calculates the integer representation of the precessor kmer given the current kmer and the precessing nucleotide
 * k_i-1 = int(k_i/base) + value(priorNt, base) * base^(kmerSize-1)
 * 
 * @param currentKmer current kmer in decimal representation
 * @param priorNt precessing nucleotide as a token
 * @param ALPHABET_SIZE number of accepted characters
 * @param stepSize equals ALPHABET_SIZE ^ (kmerSize - 1)
 * @return precessing Kmer as integer representation in the current base
 */
inline size_t precessingKmer(const size_t currentKmer, const int priorNt, const int stepSize, const int ALPHABET_SIZE) {
    return (currentKmer/ALPHABET_SIZE) + (priorNt * stepSize);
}

// ===============================================================
// ===============================================================
// ===================== Scoring calculations ====================
// ===============================================================
// ===============================================================

inline constexpr double log2Pi = 1.8378770664093453; // Precomputed log(2 * M_PI)

// https://ethz.ch/content/dam/ethz/special-interest/mavt/dynamic-systems-n-control/idsc-dam/Lectures/Stochastic-Systems/Statistical_Methods.pdf
/**
 * Calculate log pdf for a given x, mean and standard deviation
 * 
 * @param x value
 * @param m mean
 * @param s standard deviation 
 * @return probabily density at position x for N~(m, s²)
*/
inline double log_normal_pdf(const double x, const double m, const double s) {
    if (s == 0.0) {
        return -INFINITY; // Handling edge case where standard deviation is 0
    }
    
    const double variance = s * s;
    const double diff = x - m;
    
    return -0.5 * (log2Pi + log(variance) + (diff * diff) / variance);
}

/**
 * Return log probability density for a given value and a given normal distribution
 *
 * @param signal point to calculate probability density
 * @param kmer key for the model kmer:(mean, stdev) map
 * @param model map containing kmers as keys and (mean, stdev) tuples as values
 * @return log probability density value for x in the given normal distribution
 */
inline double scoreKmer(const double signal, const size_t kmer, const vector<tuple<double, double>> &model) {
    // Access elements of the model tuple directly to avoid redundant tuple creation and overhead
    const auto &[mean, stddev] = model[kmer];
    return log_normal_pdf(signal, mean, stddev);
}

/**
 * @brief Updates the transition probabilities by applying logarithmic values.
 *
 * This function checks each transition in the `transitions` map. If a transition value is `-1`, it updates
 * the transition value with the logarithmic value of the corresponding entry from the `default_transitions_vals` map.
 * Otherwise, it applies the logarithm directly to the existing transition value.
 * 
 * @param default_transitions_vals A map containing default transition values (string keys and double values).
 *                                 These default values are used when a transition value is set to `-1`.
 * @param transitions A map containing current transition values (string keys and double values).
 *                    This map is updated with logarithmic values during the function execution.
 */
void updateTransitions(const unordered_map<string, double>& default_transitions_vals, unordered_map<string, double>& transitions);