// author: Jannes Spangenberg
// e-mail: jannes.spangenberg@uni-jena.de
// github: https://github.com/JannesSP
// website: https://jannessp.github.io

#pragma once

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
#include "argparse.hpp"
#include "utils.hpp"

extern constexpr double EPSILON = 1e-8; // chose by eye just to distinguish real errors from numeric errors

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
void logF(const double *sig, const int *kmerSeq, dproxy *M, dproxy *E, const std::size_t T, const std::size_t N, const std::tuple<double, double> *model);

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
void logB(const double *sig, const int *kmerSeq, dproxy *M, dproxy *E, const std::size_t T, const std::size_t N, const std::tuple<double, double> *model);

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
void logP(double *LP, const dproxy *FOR, const dproxy *BACK, const double Z, const std::size_t S);

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
void getBorders(std::list<std::string> &segString, const double *LPM, const double *LPE, const std::size_t T, const std::size_t N, const int kmerSize);

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
 * @param kmerSize size of k-mer
 */
void traceback(std::size_t t, std::size_t n, const dproxy *M, const dproxy *E, const double *LPM, const double *LPE, std::list<std::string> &segString, const std::size_t N, const int kmerSize);

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
// void funcM(const std::size_t t, const std::size_t n, const dproxy *M, const dproxy *E, const double *LPM, const double *LPE, std::list<std::string> &segString, const std::size_t N, std::vector<double> &segProb, const int kmerSize);

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
// void funcE(const std::size_t t, const std::size_t n, const dproxy *M, const dproxy *E, const double *LPM, const double *LPE, std::list<std::string> &segString, const std::size_t N, std::vector<double> &segProb, const int kmerSize);

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
std::tuple<double, double, double> trainTransition(const double *sig, const int *kmerSeq, const dproxy *forM, const dproxy *forE, const dproxy *backM, const dproxy *backE, const std::size_t T, const std::size_t N, const std::tuple<double, double> *model);

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
std::tuple<double *, double *> trainEmission(const double *sig, const int *kmerSeq, const dproxy *forM, const dproxy *forE, const dproxy *backM, const dproxy *backE, const std::size_t T, const std::size_t N, const std::tuple<double, double> *model, const int numKmers);

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
void trainParams(const double *sig, const int *kmerSeq, const dproxy *forM, const dproxy *forE, const dproxy *backM, const dproxy *backE, const std::size_t T, const std::size_t N, std::tuple<double, double> *model, const int alphabetSize, const int numKmers, const int kmerSize);