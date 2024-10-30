// author: Jannes Spangenberg
// e-mail: jannes.spangenberg@uni-jena.de
// github: https://github.com/JannesSP
// website: https://jannessp.github.io

// ===============================================================
// ===============================================================
// =========================== Utility ===========================
// ===============================================================
// ===============================================================

#include "utils.hpp"

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
std::vector<std::size_t> column_argsort(const double *matrix, const std::size_t C, const std::size_t t)
{
    // Initialize original index locations (indices correspond to C)
    std::vector<std::size_t> idx(C);
    iota(idx.begin(), idx.end(), 0);

    // Sort indexes based on comparing values in the given column 'c'
    stable_sort(idx.begin(), idx.end(),
                [matrix, C, t](std::size_t i1, std::size_t i2)
                {
                    return matrix[t * C + i1] > matrix[t * C + i2];
                });

    return idx;
}

/**
 * C++ version 0.4 std::std::string style "itoa":
 * Contributions from Stuart Lowe, Ray-Yuan Sheu,
 * Rodrigo de Salvo Braz, Luc Gallant, John Maloney
 * and Brian Hunt
 *
 * Converts a decimal to number to a number of base alphabet_size.
 * TODO Works for base between 2 and 16 (included)
 *
 * Returns kmer in reversed direction!
 *
 * @param value input number in decimal to convert to base
 * @returns kmer as reversed std::string, should be 5' - 3' direction
 */
std::string itoa(const std::size_t value, const int alphabet_size, const int kmerSize)
{
    std::string buf;
    const int base = alphabet_size;

    // check that the base if valid
    if (base < 2 || base > 16)
        return std::to_string(value);

    enum
    {
        kMaxDigits = 35
    };
    buf.reserve(kMaxDigits); // Pre-allocate enough space.
    int quotient = value;

    // Translating number to std::string with base:
    do
    {
        buf += ID2BASE.at("0123456789abcdef"[abs(quotient % base)]);
        quotient /= base;
    } while (quotient);

    // Append the negative sign
    if (value < 0)
        buf += '-';

    while ((int)buf.length() < kmerSize)
    {
        buf += ID2BASE.at('0');
    }

    // skip this so kmer is in 5' - 3' direction for output
    // reverse( buf.begin(), buf.end() );
    return buf;
}

/**
 * Converts the kmers of the model file to the integer representation using the BASE2ID map
 *
 * @param s kmer containing nucleotides
 * @param BASE2ID base to id map
 * @param alphabet_size
 * @returns integer representation of the given kmer
 */
int kmer2int(const std::string &s, const int alphabet_size)
{
    int ret = 0;
    for (const char &c : s)
    {
        // assert (BASE2ID.at(c)>=0); // check if nucleotide is known
        ret *= alphabet_size; // move the number in base to the left
        ret += BASE2ID.at(c);
    }
    return ret;
}

#include <iostream>
/**
 * Reads the normal distribution parameters from a given TSV file,
 * and returns the kmer model and alphabet size.
 *
 * @param file       Path to the TSV file containing kmer parameters (mean, stdev).
 * @param kmerSize   The size of the kmers (length of the kmers in the file).
 * @returns          A std::tuple containing:
 *                   1. An array of tuples, where each std::tuple holds (mean, stdev) for each kmer.
 *                   2. The alphabet size (number of unique nucleotide characters from the kmer std::set).
 *                   3. The total number of possible kmers (calculated as alphabet_size^kmerSize).
 */
std::tuple<std::vector<std::tuple<double, double>>, int, std::size_t> readKmerModel(const std::string &file, const int kmerSize)
{
    std::string line, kmer, tmp;

    std::set<char> uniqueChars; // std::set to store unique characters from kmers to determine the alphabet size
    std::ifstream inputFile(file);

    // First pass: read file to collect unique characters from kmer sequences
    // Skip the header line
    getline(inputFile, line);
    while (getline(inputFile, line))
    {                                   // read line
        std::stringstream buffer(line); // parse line to std::stringstream for getline
        getline(buffer, kmer, '\t');
        // Add all unique characters in the kmer to the std::set
        for (char c : kmer)
        {
            uniqueChars.insert(c);
        }
    }
    inputFile.close();

    const int alphabet_size = (int)uniqueChars.size();
    const std::size_t numKmers = pow(alphabet_size, kmerSize);
    uniqueChars.clear(); // Clear the unique character std::set (no longer needed) to free up memory

    std::vector<std::tuple<double, double>> model(numKmers);
    inputFile.open(file); // Reopen the file for the second pass

    // Read through the file to populate the model with (mean, stdev) for each kmer
    // Skip the header line
    getline(inputFile, line);
    while (getline(inputFile, line))
    {                                   // read line
        std::stringstream buffer(line); // parse line to std::stringstream for getline
        getline(buffer, kmer, '\t');
        // legacy models are stored from 3' - 5'
        // https://github.com/nanoporetech/kmer_models
        // new models are stored in 5' - 3'
        reverse(kmer.begin(), kmer.end()); // 5-3 -> 3-5 orientation
        getline(buffer, tmp, '\t');        // level_mean
        const double mean = stod(tmp);
        getline(buffer, tmp, '\t'); // level_stdv
        const double stdev = stod(tmp);
        // model.push_back(std::make_tuple(mean, stdev));
        model[kmer2int(kmer, alphabet_size)] = std::make_tuple(mean, stdev);
    }
    inputFile.close();

    // Return a std::tuple containing:
    // 1. The model array (containing kmers and their (mean, stdev) tuples)
    // 2. The alphabet size (number of unique characters in the kmers)
    // 3. The total number of kmers (alphabet_size^kmerSize)
    return std::make_tuple(model, alphabet_size, numKmers);
}

/**
 * @brief Updates the transition probabilities by applying logarithmic values.
 *
 * This function checks each transition in the `transitions` map. If a transition value is `-1`, it updates
 * the transition value with the logarithmic value of the corresponding entry from the `default_transitions_vals` map.
 * Otherwise, it applies the logarithm directly to the existing transition value.
 *
 * @param default_transitions_vals A map containing default transition values (std::string keys and double values).
 *                                 These default values are used when a transition value is std::set to `-1`.
 * @param transitions A map containing current transition values (std::string keys and double values).
 *                    This map is updated with logarithmic values during the function execution.
 */
void updateTransitions(const std::unordered_map<std::string, double> &default_transitions_vals, std::unordered_map<std::string, double> &transitions)
{
    // Iterate over each transition in the 'transitions' std::unordered_map
    for (const auto &[param, value] : transitions)
    {
        // Check if the transition value is -1, which indicates it should use the default value
        if (value == -1.0)
        {
            // Fetch the default value for this key, and store the value in the transition map
            transitions[param] = default_transitions_vals.at(param);
        }
        // Apply the logarithmic value to the current transition value
        transitions[param] = log(value);
    }
}

// Function to calculate the median of a std::vector
double calculateMedian(std::vector<double> &vec)
{
    const std::size_t size = vec.size();

    // Sort the std::vector
    std::sort(vec.begin(), vec.end());

    // If the std::vector size is odd, return the middle element
    if (size % 2 == 1)
    {
        return vec[size / 2];
    }
    // If the std::vector size is even, return the average of the two middle elements
    else
    {
        return (vec[size / 2 - 1] + vec[size / 2]) / 2.0;
    }
}