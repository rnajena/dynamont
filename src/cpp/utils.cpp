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
std::vector<std::size_t> columnArgsort(const double *matrix, const std::size_t C, const std::size_t t)
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
std::string itoa(const std::size_t value, const int alphabetSize, const int kmerSize, const bool rna)
{
    std::string buf;
    const int base = alphabetSize;

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

    // skip this for RNA so kmer is in 5' - 3' direction for output
    if (!rna)
        // for DNA this must be reversed to output 5' - 3' direction
        std::reverse(buf.begin(), buf.end());

    return buf;
}

/**
 * Converts the kmers of the model file to the integer representation using the BASE2ID map
 *
 * @param s kmer containing nucleotides
 * @param BASE2ID base to id map
 * @param alphabetSize
 * @returns integer representation of the given kmer
 */
int kmer2int(const std::string &s, const int alphabetSize)
{
    int ret = 0;
    for (const char &c : s)
    {
        ret *= alphabetSize; // move the number in base to the left
        ret += BASE2ID.at(c);
    }
    return ret;
}

/**
 * Reads the normal distribution parameters from a given TSV file,
 * and returns the kmer model and alphabet size.
 *
 * @param file       Path to the TSV file containing kmer parameters (mean, stdev).
 * @param kmerSize   The size of the kmers (length of the kmers in the file).
 * @param rna        True if input is RNA sequence, false if DNA sequence
 * @returns          A std::tuple containing:
 *                   1. An array of tuples, where each std::tuple holds (mean, stdev) for each kmer.
 *                   2. The alphabet size (number of unique nucleotide characters from the kmer std::set).
 *                   3. The total number of possible kmers (calculated as alphabetSize^kmerSize).
 */
std::tuple<std::tuple<double, double> *, int, std::size_t> readKmerModel(const std::string &file, const std::size_t kmerSize, const bool rna)
{
    std::string line, kmer, tmp;

    std::set<char> uniqueChars; // std::set to store unique characters from kmers to determine the alphabet size
    std::ifstream inputFile(file);

    // First pass: read file to collect unique characters from kmer sequences
    // Skip the header line
    getline(inputFile, line);
    while (getline(inputFile, line)) // read line
    {
        std::stringstream buffer(line); // parse line to std::stringstream for getline
        getline(buffer, kmer, '\t');
        // Add all unique characters in the kmer to the std::set
        for (char c : kmer)
        {
            uniqueChars.insert(c);
        }
        if (kmer.length() != kmerSize)
        {
            std::cerr << kmer << " kmer length in model " << file << " does not match kmerSize " << kmerSize << " of pore given pore type" << std::endl;
            exit(6);
        }
    }
    inputFile.close();

    const int alphabetSize = (int)uniqueChars.size();
    const std::size_t numKmers = pow(alphabetSize, kmerSize);
    uniqueChars.clear(); // Clear the unique character std::set (no longer needed) to free up memory

    std::tuple<double, double> *model = new std::tuple<double, double>[numKmers];
    inputFile.open(file); // Reopen the file for the second pass

    // Read through the file to populate the model with (mean, stdev) for each kmer
    // Skip the header line
    getline(inputFile, line);
    while (getline(inputFile, line))
    {                                   // read line
        std::stringstream buffer(line); // parse line to std::stringstream for getline
        getline(buffer, kmer, '\t');
        // models are stored in 5' - 3'
        if (rna)
            std::reverse(kmer.begin(), kmer.end()); // 5-3 -> 3-5 orientation
        getline(buffer, tmp, '\t');                 // level_mean
        const double mean = stod(tmp);
        getline(buffer, tmp, '\t'); // level_stdv
        const double stdev = stod(tmp);
        model[kmer2int(kmer, alphabetSize)] = std::make_tuple(mean, stdev);
    }
    inputFile.close();

    // Return a std::tuple containing:
    // 1. The model array (containing kmers and their (mean, stdev) tuples)
    // 2. The alphabet size (number of unique characters in the kmers)
    // 3. The total number of kmers (alphabetSize^kmerSize)
    return std::make_tuple(model, alphabetSize, numKmers);
}

/**
 * @brief Updates the transition probabilities by applying logarithmic values.
 *
 * This function checks each transition in the `transitions` map. If a transition value is `-1`, it updates
 * the transition value with the logarithmic value of the corresponding entry from the `defaultVals` map.
 * Otherwise, it applies the logarithm directly to the existing transition value.
 *
 * @param defaultVals A map containing default transition values (std::string keys and double values).
 *                                 These default values are used when a transition value is std::set to `-1`.
 * @param newVals A map containing current transition values (std::string keys and double values).
 *                    This map is updated with logarithmic values during the function execution.
 */
void updateTransitions(const std::unordered_map<std::string, double> &defaultVals, std::unordered_map<std::string, double> &newVals)
{
    // Iterate over each transition in the 'newVals' std::unordered_map
    for (const auto &[param, value] : newVals)
    {
        // Check if the transition value is -1, which indicates it should use the default value
        if (value == -1.0)
        {
            // Fetch the default value for this key, and store the value in the transition map
            newVals[param] = defaultVals.at(param);
        }
        // Apply the logarithmic value to the current transition value
        newVals[param] = std::log(value);
    }
}

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
 * // median will be 3.0
 */
double calculateMedian(std::vector<double> &vec)
{
    const std::size_t size = vec.size();

    if (!size)
    {
        throw std::out_of_range("Could not calculate median from empty vector");
    }

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
std::string formattedMedian(std::vector<double> &vec)
{
    double median = calculateMedian(vec);
    const int precision = 5;

    // Create a string stream to format the double
    std::ostringstream out;
    out << std::fixed << std::setprecision(precision) << median;

    return out.str();
}

void checkModelpath(std::string modelpath)
{
    if (modelpath.empty() || !std::filesystem::exists(modelpath))
    {
        std::cerr << "Please provide a valid modelpath: " << modelpath << std::endl;
        exit(7);
    }
}

void checkInput(const std::size_t signalSize, const std::size_t readSize, const std::size_t kmerSize)
{
    if (signalSize < 1)
    {
        std::cerr << "Signal: " << signalSize << " smaller than 1" << std::endl;
        exit(8);
    }
    if (readSize < 1)
    {
        std::cerr << "Read: " << readSize << " smaller than 1" << std::endl;
        exit(9);
    }
    if (signalSize < 2 * readSize) // each segment has at least length 2
    {
        std::cerr << "Signal: " << signalSize << " smaller than read: " << readSize << std::endl;
        exit(10);
    }
    if (readSize < kmerSize)
    {
        std::cerr << "Read: " << readSize << " smaller than kmerSize of the pore type: " << kmerSize << std::endl;
        exit(11);
    }
}