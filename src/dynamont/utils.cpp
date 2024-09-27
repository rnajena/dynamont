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
 * @return size_t vector with the sorted index of column in descending order
 */
vector<size_t> column_argsort(const double* matrix, const size_t C, const size_t t) {
    // Initialize original index locations (indices correspond to C)
    vector<size_t> idx(C);
    iota(idx.begin(), idx.end(), 0);

    // Sort indexes based on comparing values in the given column 'c'
    stable_sort(idx.begin(), idx.end(),
        [matrix, C, t](size_t i1, size_t i2) {
            return matrix[t * C + i1] > matrix[t * C + i2];
        });

    return idx;
}

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
 * @returns kmer as reversed string, should be 5' - 3' direction
*/
string itoa(const size_t value, const int ALPHABET_SIZE, const int kmerSize) {
    string buf;
    int base = ALPHABET_SIZE;

    // check that the base if valid
    if (base < 2 || base > 16) return to_string(value);

    enum { kMaxDigits = 35 };
    buf.reserve( kMaxDigits ); // Pre-allocate enough space.
    int quotient = value;

    // Translating number to string with base:
    do {
        buf += ID2BASE.at("0123456789abcdef"[ abs( quotient % base ) ]);
        quotient /= base;
    } while ( quotient );

    // Append the negative sign
    if ( value < 0) buf += '-';

    while ((int) buf.length() < kmerSize) {
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
 * @param ALPHABET_SIZE
 * @returns integer representation of the given kmer
 */
int kmer2int(const string &s, const int ALPHABET_SIZE) {
    int ret = 0;
    for(char const &c:s){
        // assert (BASE2ID.at(c)>=0); // check if nucleotide is known
        ret*=ALPHABET_SIZE; // move the number in base to the left
        ret+=BASE2ID.at(c);
    }
    return ret;
}

#include <iostream>
/**
 * Read the normal distribution parameters from a given TSV file
 *
 * @param file path to the TSV file containing the parameters
 * @param model kmer model to fill
 * @param ALPHABET_SIZE
 */
void readKmerModel(const string &file, vector<tuple<double, double>> &model, const int ALPHABET_SIZE) {
    ifstream inputFile(file);
    string line, kmer, tmp;
    double mean, stdev;
    getline(inputFile, line);
    while(getline(inputFile, line)) { // read line
        stringstream buffer(line); // parse line to stringstream for getline
        getline(buffer, kmer, '\t');
        // legacy models are stored from 3' - 5'
        // https://github.com/nanoporetech/kmer_models
        // new models are stored in 5' - 3'
        reverse(kmer.begin(), kmer.end()); // 5-3 -> 3-5 orientation
        getline(buffer, tmp, '\t'); // level_mean
        mean = stod(tmp);
        getline(buffer, tmp, '\t'); // level_stdv
        stdev = stod(tmp);
        model[kmer2int(kmer, ALPHABET_SIZE)]=make_tuple(mean, stdev);
    }
    inputFile.close();
}