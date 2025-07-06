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
string itoa(const size_t value, const int kmerSize) {
    string buf;
    int base = kmerSize;

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
    // if ( value < 0) buf += '-';

    while ((int) buf.length() < base) {
        buf += ID2BASE.at('0');
    }

    // skip this so kmer is in 5' - 3' direction for output
    // reverse( buf.begin(), buf.end() );
    return buf;
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
string itoa(const int value, const int kmerSize) {
    string buf;
    int base = kmerSize;

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

    while ((int) buf.length() < base) {
        buf += ID2BASE.at('0');
    }

    // skip this so kmer is in 5' - 3' direction for output
    // reverse( buf.begin(), buf.end() );
    return buf;
}

/**
 * Converts a number of base ALPHABET_SIZE to a decimal number.
 * Works ONLY if ALPHABET_SIZE is smaller or equal to 10!
 * 
 * @param i input number in the given base as an array
*/
int toDeci(const int* i, const int kmerSize) {
    int ret = 0;
    int m = 1;
    for(int r = kmerSize - 1; r >= 0; r--) {
        ret += m*i[r];
        m *= ALPHABET_SIZE;
    }
    return ret;
}

/**
 * Converts the kmers of the model file to the integer representation using the BASE2ID map
 *
 * @param s kmer containing nucleotides 
 * @param BASE2ID base to id map 
 * @param kmerSize kmer size 
 * @returns integer representation of the given kmer
 */
int kmer2int(const string &s, const int kmerSize) {
    int ret = 0;
    for(char const &c:s){
        // assert (BASE2ID.at(c)>=0); // check if nucleotide is known
        ret*=kmerSize; // move the number in base to the left
        ret+=BASE2ID.at(c);
    }
    return ret;
}

/**
 * Convert the read sequence to a kmer sequence which is represented by integers.
 * 
 * @param seq read sequence
 * @param N length of the read sequence, number of nucleotides
 * @param kmerSize kmer size 
 * @return kmer sequence in integer representation
*/
int* seq2kmer(const int* seq, const size_t N, const int kmerSize) {
    int* kmer_seq = new int[N];
    int* tempKmer = new int[kmerSize];
    for(size_t n=0; n<N; n++){ // extend loop to ad 2 Ns at start and end of read
        copy(seq+n, seq+n+kmerSize, tempKmer);
        kmer_seq[n] = toDeci(tempKmer, kmerSize);
    }
    delete[] tempKmer;
    return kmer_seq;
}

/**
 * Read the normal distribution parameters from a given TSV file
 *
 * @param file path to the TSV file containing the parameters
 * @param model kmer model to fill
 */
void readKmerModel(const string &file, vector<tuple<double, double>> &model, const int kmerSize) {
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
        mean = atof(tmp.c_str());
        getline(buffer, tmp, '\t'); // level_stdv
        stdev = atof(tmp.c_str());
        model[kmer2int(kmer, kmerSize)]=make_tuple(mean, stdev);
    }
    inputFile.close();
}