#include <iostream> // cout
#include <fstream> // file io
#include <sstream> // file io
#include <string>
#include <map> // dictionary
#include <tuple>
#include <algorithm> // reverse strings
#include <filesystem> // reading files in directory (fast5)
#include <vector>
#include <cmath> // exp
#include <limits> // for inifinity
#include <math.h> // log1p

// #include "dp.h"

using namespace std;
namespace fs = filesystem;

// Asserts floating point compatibility at compile time
// necessary for INFINITY usage
static_assert(numeric_limits<float>::is_iec559, "IEEE 754 required");

// https://stackoverflow.com/questions/10847007/using-the-gaussian-probability-density-function-in-c
/**
 * Calculate probabilty value for a given x, mean and standard deviation of a normal distribution
 *
 * @param x value
 * @param m mean
 * @param s standard deviation 
 */
float normal_pdf(float x, float m, float s) {
    static const float inv_sqrt_2pi = 0.3989422804014327;
    float a = (x - m) / s;

    return inv_sqrt_2pi / s * exp(-0.5f * a * a);
}

/**
 * Return probability density for a given value and a given normal distribution
 *
 * @param x point to calculate probability density
 * @param kmer key for the model kmer:(mean, stdev) map
 * @param model map containing kmers as keys and (mean, stdev) tuples as values
 * @return probability density value for x in the given normal distribution
 */
float score5mer(float x, string kmer, map<string,tuple<float, float>>* model) {
    tuple<float, float> kmerModel = model->find(kmer)->second;
    return normal_pdf(x, get<0>(kmerModel), get<1>(kmerModel));
}

/**
 * Calculate forward matrices
 *
 * @param signal ONT raw signal with pA values
 * @param seq nucleotide sequence represented by the ONT signal
 * @param MM matrix containing values for segment borders 
 * @param MC matrix containing values for extending segment
 * @param T length of the ONT raw signal
 * @param N length of nucleotide sequence
 * @param model map containing kmers as keys and (mean, stdev) tuples as values
 */
void forward(float* signal, string* seq, float* MM, float* MC, int T, int N, map<string,tuple<float, float>>* model){
    // int T = sizeof(signal)/sizeof(*signal);
    // int N = seq->size();
    // float MM[T*N] = {-1};
    // float MC[T*N] = {-1};
    // TODO start at j = 2? 5mers?
    for(int i=0;i<T;i++){
        for(int j=0;j<N;j++){
            float mm=0;
            if(i>0 && j>0){
                mm+=MC[(i-1)*N+(j-1)]*score5mer(signal[i], seq->substr(j-2, 5), model);
            }
            MM[i*N+j] = mm;
            float mc=0;
            if(i>0){
                mc+=MC[(i-1)*N+j]*score5mer(signal[i], seq->substr(j-2, 5), model);
                mc+=MM[(i-1)*N+j]*score5mer(signal[i], seq->substr(j-2, 5), model);
            }
            if(i==0 && j==0){
                mc+=1;
            }
            MC[i*N+j]=mc;
        }
    }
}

/**
 * Calculate addition of a+b in log space as efficiently as possible
 *
 * @param a first value
 * @param b second value
 * @returns log(exp(a) + exp(b))
 */
float logplus(float a, float b) {
    // safety check, should not occure
    // if(a==b && isinf(a) && isinf(b)) {
    //     return a;
    // }
    if(a>=b){
        return a + log1p(b-a);
    }
    return b + log1p(a-b);
}

/**
 * Calculate forward matrices using logarithmic values
 *
 * @param signal ONT raw signal with pA values
 * @param seq nucleotide sequence represented by the ONT signal
 * @param MM matrix containing values for segment borders 
 * @param MC matrix containing values for extending segment
 * @param T length of the ONT raw signal
 * @param N length of nucleotide sequence
 * @param model map containing kmers as keys and (mean, stdev) tuples as values
 */
void logForward(float* signal, string* seq, float* MM, float* MC, int T, int N, map<string,tuple<float, float>>* model){
    // int T = sizeof(signal)/sizeof(*signal);
    // int N = seq->size();
    // float MM[T*N] = {-1};
    // float MC[T*N] = {-1};
    for(int i=0;i<T;i++){
        for(int j=2;j<N-2;j++){
            float mm=-INFINITY;
            if(i>0 && j>2){
                mm = logplus(mm, MC[(i-1)*N+(j-1)] + log(score5mer(signal[i], seq->substr(j-2, 5), model)));
            }
            MM[i*N+j] = mm;
            float mc=-INFINITY;
            if(i>0){
                mc = logplus(mc, MC[(i-1)*N+j] + log(score5mer(signal[i], seq->substr(j-2, 5), model)));
                mc = logplus(mc, MM[(i-1)*N+j] + log(score5mer(signal[i], seq->substr(j-2, 5), model)));
            }
            if(i==0 && j==2){
                mc+=0;
            }
            MC[i*N+j]=mc;
        }
    }
}

/**
 * Calculate backward matrices using logarithmic values
 *
 * @param signal ONT raw signal with pA values
 * @param seq nucleotide sequence represented by the ONT signal
 * @param MM matrix containing values for segment borders 
 * @param MC matrix containing values for extending segment
 * @param T length of the ONT raw signal
 * @param N length of nucleotide sequence
 * @param model map containing kmers as keys and (mean, stdev) tuples as values
 */
void logBackward(float* signal, string* seq, float* MM, float* MC, int T, int N, map<string,tuple<float, float>>* model) {
    for(int i=T-1;i>=0;i--){
        for(int j=N-3;j>=2;j--){
            float mm=-INFINITY;
            if(i<T-1 && j<N-3){
                mm = logplus(mm, MC[(i+1)*N+(j+1)] + log(score5mer(signal[i], seq->substr(j-2, 5), model)));
            }
            MM[i*N+j] = mm;
            float mc=-INFINITY;
            if(i<T-1){
                mc = logplus(mc, MC[(i+1)*N+j] + log(score5mer(signal[i], seq->substr(j-2, 5), model)));
                mc = logplus(mc, MM[(i+1)*N+j] + log(score5mer(signal[i], seq->substr(j-2, 5), model)));
            }
            if(i==T-1 && j==N-3){
                mc+=0;
            }
            MC[i*N+j]=mc;
        }
    }
}

// https://stackoverflow.com/questions/11140483/how-to-get-list-of-files-with-a-specific-extension-in-a-given-folder
/**
 * Collect all fast5 files in the given directory
 *
 * @param directory path to the parent directory containing files with the given extension
 * @param ext collecting files containing this extension
 * @return list of all files with the given extension in the directory
 */
vector<fs::path> extractFiles(fs::path const & directory, string const & ext = ".fast5") {
    vector<fs::path> paths;
    if (fs::exists(directory) && fs::is_directory(directory)) {
        for (auto const & entry : fs::recursive_directory_iterator(directory)) {
            if (fs::is_regular_file(entry) && entry.path().extension() == ext)
                paths.emplace_back(entry.path().filename());
        }
    }
    // else throw exception?
    return paths;
}  

/**
 * Read the normal distribution parameters from a given TSV file
 *
 * @param file path to the TSV file containing the parameters
 * @param return map containing kmers as keys and (mean, stdev) tuples as values
 */
map<string,tuple<float, float>> readKmerModel(string file) {
    ifstream inputFile;
    inputFile.open(file);
    string line, kmer, tmp;
    float mean, stdev;
    map<string,tuple<float, float>> model;

    while(getline(inputFile, line)) { // read line
        stringstream inputString(line); // parse line to inputString
        getline(inputString, kmer, '\t');
        reverse(kmer.begin(), kmer.end()); // 3-5 -> 5-3 orientation
        getline(inputString, tmp, '\t');
        mean = atof(tmp.c_str());
        getline(inputString, tmp, '\t');
        stdev = atof(tmp.c_str());
        model[kmer]=make_tuple(mean, stdev);
    }
    return model;
}


// https://stackoverflow.com/questions/865668/parsing-command-line-arguments-in-c
/**
 * Return the value of the command line argument
 *
 * @param begin command line arguments
 * @param end end of command line arguments
 * @param option command line argument to search for
 * @return value of the given argument if there is one
 */
char* getCmdOption(char ** begin, char ** end, const string & option) {

    char ** itr = find(begin, end, option);
    if (itr != end && ++itr != end)
    {
        return *itr;
    }
    return 0;
}

// https://stackoverflow.com/questions/865668/parsing-command-line-arguments-in-c
/**
 * Checks if command line argument exists
 * 
 * @param begin command line arguments
 * @param end end of command line arguments
 * @param option command line argument to search for
 * @return True or False
 */
bool cmdOptionExists(char** begin, char** end, const string& option) {
    return find(begin, end, option) != end;
}

int main(int argc, char* argv[]) {

    if(cmdOptionExists(argv, argv+argc, "-h")) {
        // TODO print help message
    }

    // input fast5 parent directory
    char * fast5_directory = getCmdOption(argv, argv + argc, "-f5");

    if (fast5_directory) {
        vector<fs::path> fast5s = extractFiles(fast5_directory);
    }
}