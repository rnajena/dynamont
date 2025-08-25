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
#include <vector>
#include <cmath> // exp
#include <assert.h>
#include <stdlib.h>
#include <algorithm>
#include <unistd.h>
#include "../include/argparse.hpp"
#include "../include/utils.hpp"

// TODO split main and rest of functions
// TODO create and use polyA.hpp
// TODO for later, runtime improvement possible by writing explicit emission functions with loc, scale and df as constexpr, only sig_val is variable as parameter

// TODO do not use this, write explicit std::cin, std::cout, std::exp, std::... -> DONE 
//using namespace std;

inline constexpr double EPSILON = 1e-5; // chose by eye just to distinguish real errors from numeric errors

// Asserts doubleing point compatibility at compile time  // ?
// necessary for INFINITY usage
static_assert(std::numeric_limits<double>::is_iec559, "IEEE 754 required");

//! ------------------------------------------ PDFs, Forward, Backward & Posterior Probability ----------------------------------------------

/**
 * DIST & PARAM IN -> 60 READS :
 * adapter t, df: 5.612094 loc: -0.759701 scale: 0.535895
 * polyA t, df: 6.022091, loc: 0.839093, scale: 0.217290
 * leader gumbel l, loc: 0.927918 , scale: 0.398849
 * transcript gumbel r, loc: -0.341699 , scale: 0.890093
 * start gumbel r, loc: -1.552134, scale: 0.415937
 */

// TODO make function inline
/**
 * logarithm t distribution PDF  : checked the correctness with scipy.stats.t
 */
double log_t_pdf(const double sig_val, const double loc, const double scale, const double df)
{

    const double pi = 3.14159265358979323846;
    const double diff = (sig_val - loc) / scale;
    const double logGammaNuPlusOneHalf = lgamma((df + 1.0) / 2.0);
    const double logGammaNuHalf = lgamma(df / 2.0);

    return logGammaNuPlusOneHalf - logGammaNuHalf - 0.5 * log(df * pi * scale * scale) - (df + 1.0) / 2.0 * log(1.0 + (diff * diff) / df);
}

// TODO make function inline
/**
 * logarithm gumbel left skewed PDF : checked the correctness with scipy.stats.gumbel_l
 */
double log_gumbel_l_pdf(const double sig_val, const double loc, const double scale)
{
    // if (scale == 0.0) {
    //     return -INFINITY; // Handling edge case where beta (scale) is 0
    // }

    const double z = -(sig_val - loc) / scale;

    return -z - exp(-z);
}

// TODO make function inline
/**
 * logarithm gumbel right skewed PDF : checked with scipy.stats.gumbel_r, ->  //! around 0.92 different with scipy.stat.gumbel_r
 */
double log_gumbel_r_pdf(const double sig_val, const double loc, const double scale)
{
    // if (scale == 0.0) {
    //     return -INFINITY; // Handling edge case where beta (scale) is 0
    // }

    const double z = (sig_val - loc) / scale;

    return -z - exp(-z);
}

/**
 * Calculate forward matrices using logarithmic values
 * 1D array for each state : 5 1D arrays
 * S L A PA TR : initialized matrices for each state
 */
void logF(double *sig, double *S, double *L, double *A, double *PA, double *TR, size_t T,
          double s, double l1, double l2, double a1, double a2, double pa1, double pa2, double tr1, double tr2)
{
    double start, leader, adapter, polya, transcript;

    S[0] = 0;

    for (size_t t = 1; t < T; ++t)
    {
        // init state accumulators
        start = -INFINITY;
        leader = -INFINITY;
        adapter = -INFINITY;
        polya = -INFINITY;
        transcript = -INFINITY;

        // calculate probabilities
        //       accumulator + (prevV *                 emission                * transition)
        start = logPlus(start, S[t - 1] + log_gumbel_r_pdf(sig[t - 1], -1.552134, 0.415937) + s);

        leader = logPlus(leader, S[t - 1] + log_gumbel_l_pdf(sig[t - 1], 0.927918, 0.398849) + l1);
        leader = logPlus(leader, L[t - 1] + log_gumbel_l_pdf(sig[t - 1], 0.927918, 0.398849) + l2);

        adapter = logPlus(adapter, L[t - 1] + log_t_pdf(sig[t - 1], -0.759701, 0.535895, 5.612094) + a1);
        adapter = logPlus(adapter, A[t - 1] + log_t_pdf(sig[t - 1], -0.759701, 0.535895, 5.612094) + a2);

        polya = logPlus(polya, A[t - 1] + log_t_pdf(sig[t - 1], 0.839093, 0.217290, 6.022091) + pa1);
        polya = logPlus(polya, PA[t - 1] + log_t_pdf(sig[t - 1], 0.839093, 0.217290, 6.022091) + pa2);

        transcript = logPlus(transcript, PA[t - 1] + log_gumbel_r_pdf(sig[t - 1], -0.341699, 0.890093) + tr1);
        transcript = logPlus(transcript, TR[t - 1] + log_gumbel_r_pdf(sig[t - 1], -0.341699, 0.890093) + tr2);

        // start = logPlus(start, TR[t-1] + log_gumbel_r_pdf(sig[t-1], -1.552134, 0.415937) + s0);

        // update state matrices
        S[t] = start;
        L[t] = leader;
        A[t] = adapter;
        PA[t] = polya;
        TR[t] = transcript;

        // TODO compress code
        // S[t] = S[t - 1] + log_gumbel_r_pdf(sig[t - 1], -1.552134, 0.415937) + s;
        // L[t] = logPlus(S[t - 1] + log_gumbel_l_pdf(sig[t - 1], 0.927918, 0.398849) + l1, L[t - 1] + log_gumbel_l_pdf(sig[t - 1], 0.927918, 0.398849) + l2);
        // ...
    }
}

/**
 * Calculate backward matrices using logarithmic values
 */
void logB(double *sig, double *S, double *L, double *A, double *PA, double *TR, size_t T,
          double s, double l1, double l2, double a1, double a2, double pa1, double pa2, double tr1, double tr2)
{

    double start, leader, adapter, polya, transcript;

    TR[T - 1] = 0;

    for (size_t t = T - 1; t-- > 0;)
    { // T-2, ..., 1, 0
        // init state accumulators
        start = -INFINITY;
        leader = -INFINITY;
        adapter = -INFINITY;
        polya = -INFINITY;
        transcript = -INFINITY;

        // calculate probabilities
        //       accumulator + (prevV *         emission(t)                   * transition)
        start = logPlus(start, S[t + 1] + log_gumbel_r_pdf(sig[t], -1.552134, 0.415937) + s);
        start = logPlus(start, L[t + 1] + log_gumbel_l_pdf(sig[t], 0.927918, 0.398849) + l1);

        leader = logPlus(leader, L[t + 1] + log_gumbel_l_pdf(sig[t], 0.927918, 0.398849) + l2);
        leader = logPlus(leader, A[t + 1] + log_t_pdf(sig[t], -0.759701, 0.535895, 5.612094) + a1);

        adapter = logPlus(adapter, A[t + 1] + log_t_pdf(sig[t], -0.759701, 0.535895, 5.612094) + a2);
        adapter = logPlus(adapter, PA[t + 1] + log_t_pdf(sig[t], 0.839093, 0.217290, 6.022091) + pa1);

        polya = logPlus(polya, PA[t + 1] + log_t_pdf(sig[t], 0.839093, 0.217290, 6.022091) + pa2);
        polya = logPlus(polya, TR[t + 1] + log_gumbel_r_pdf(sig[t], -0.341699, 0.890093) + tr1);

        transcript = logPlus(transcript, TR[t + 1] + log_gumbel_r_pdf(sig[t], -0.341699, 0.890093) + tr2);

        // transcript = logPlus(transcript, S[t+1] + log_gumbel_r_pdf(sig[t], -1.552134, 0.415937) + s0);

        // update state matrices
        S[t] = start;
        L[t] = leader;
        A[t] = adapter;
        PA[t] = polya;
        TR[t] = transcript;

        // TODO compress code
    }
}

/**
 * Calculate the logarithmic probability matrix - posterior probability
 */
double *logP(const double *F, const double *B, const double Z, const size_t T)
{
    double *LP = new double[T];
    for (size_t t = 0; t < T; ++t)
    {
        LP[t] = F[t] + B[t] - Z;
    }
    return LP;
}

//! --------------------------------------------------------- BACKTRACING SECTION ------------------------------------------------------

/**
 * define backtracing function after each state
 */

// Backtracking Funcs Declaration
void funcTR(const size_t t, const double *S, const double *L, const double *A, const double *PA, const double *TR,
            const double *LPS, const double *LPL, const double *LPA, const double *LPPA, const double *LPTR,
            std::list<std::string> &segString, std::vector<size_t> &borders, std::string prevState);

void funcS(const size_t t, const double *S, const double *L, const double *A, const double *PA, const double *TR,
           const double *LPS, const double *LPL, const double *LPA, const double *LPPA, const double *LPTR,
           std::list<std::string> &segString, std::vector<size_t> &borders, std::string prevState);

void funcL(const size_t t, const double *S, const double *L, const double *A, const double *PA, const double *TR,
           const double *LPS, const double *LPL, const double *LPA, const double *LPPA, const double *LPTR,
           std::list<std::string> &segString, std::vector<size_t> &borders, std::string prevState);

void funcA(const size_t t, const double *S, const double *L, const double *A, const double *PA, const double *TR,
           const double *LPS, const double *LPL, const double *LPA, const double *LPPA, const double *LPTR,
           std::list<std::string> &segString, std::vector<size_t> &borders, std::string prevState);

void funcPA(const size_t t, const double *S, const double *L, const double *A, const double *PA, const double *TR,
            const double *LPS, const double *LPL, const double *LPA, const double *LPPA, const double *LPTR,
            std::list<std::string> &segString, std::vector<size_t> &borders, std::string prevState);

void funcS(const size_t t, const double *S, const double *L, const double *A, const double *PA, const double *TR,
           const double *LPS, const double *LPL, const double *LPA, const double *LPPA, const double *LPTR, std::list<std::string> &segString, std::vector<size_t> &borders, std::string prevState)
{

    // base case only in S as last region
    if (t == 0)
    {
        return;
    }

    if (S[t] == S[t - 1] + LPS[t])
    {
        prevState = "START";
        segString.push_back(prevState);
        funcS(t - 1, S, L, A, PA, TR, LPS, LPL, LPA, LPPA, LPTR, segString, borders, prevState);
    }

    /*
     */
    if (S[t] == TR[t - 1] + LPS[t])
    {
        const size_t border_start = t;
        borders.push_back(border_start);
        prevState = "TRANSCRIPT";
        segString.push_back(prevState);
        funcTR(t - 1, S, L, A, PA, TR, LPS, LPL, LPA, LPPA, LPTR, segString, borders, prevState);
    }
}

void funcL(const size_t t, const double *S, const double *L, const double *A, const double *PA, const double *TR,
           const double *LPS, const double *LPL, const double *LPA, const double *LPPA, const double *LPTR, std::list<std::string> &segString, std::vector<size_t> &borders, std::string prevState)
{

    if (L[t] == S[t - 1] + LPL[t])
    {
        const size_t border_start = t;
        borders.push_back(border_start);
        prevState = "START";
        segString.push_back(prevState);
        funcS(t - 1, S, L, A, PA, TR, LPS, LPL, LPA, LPPA, LPTR, segString, borders, prevState);
    }

    if (L[t] == L[t - 1] + LPL[t])
    {
        prevState = "LEADER";
        segString.push_back(prevState);
        funcL(t - 1, S, L, A, PA, TR, LPS, LPL, LPA, LPPA, LPTR, segString, borders, prevState);
    }
}

void funcA(const size_t t, const double *S, const double *L, const double *A, const double *PA, const double *TR,
           const double *LPS, const double *LPL, const double *LPA, const double *LPPA, const double *LPTR, std::list<std::string> &segString, std::vector<size_t> &borders, std::string prevState)
{

    if (A[t] == L[t - 1] + LPA[t])
    {
        const size_t border_leader = t;
        borders.push_back(border_leader);
        prevState = "LEADER";
        segString.push_back(prevState);
        funcL(t - 1, S, L, A, PA, TR, LPS, LPL, LPA, LPPA, LPTR, segString, borders, prevState);
    }

    if (A[t] == A[t - 1] + LPA[t])
    {
        prevState = "ADAPTOR";
        segString.push_back(prevState);
        funcA(t - 1, S, L, A, PA, TR, LPS, LPL, LPA, LPPA, LPTR, segString, borders, prevState);
    }
}

void funcPA(const size_t t, const double *S, const double *L, const double *A, const double *PA, const double *TR,
            const double *LPS, const double *LPL, const double *LPA, const double *LPPA, const double *LPTR, std::list<std::string> &segString, std::vector<size_t> &borders, std::string prevState)
{

    if (PA[t] == A[t - 1] + LPPA[t])
    {
        const size_t border_adaptor = t;
        borders.push_back(border_adaptor);
        prevState = "ADAPTOR";
        segString.push_back(prevState);
        funcA(t - 1, S, L, A, PA, TR, LPS, LPL, LPA, LPPA, LPTR, segString, borders, prevState);
    }

    if (PA[t] == PA[t - 1] + LPPA[t])
    {
        prevState = "POLYA";
        segString.push_back(prevState);
        funcPA(t - 1, S, L, A, PA, TR, LPS, LPL, LPA, LPPA, LPTR, segString, borders, prevState);
    }
}

void funcTR(const size_t t, const double *S, const double *L, const double *A, const double *PA, const double *TR,
            const double *LPS, const double *LPL, const double *LPA, const double *LPPA, const double *LPTR, std::list<std::string> &segString, std::vector<size_t> &borders, std::string prevState)
{
    if (TR[t] == PA[t - 1] + LPTR[t])
    {

        const size_t border_polyA = t;
        borders.push_back(border_polyA);
        prevState = "POLYA";
        segString.push_back(prevState);
        funcPA(t - 1, S, L, A, PA, TR, LPS, LPL, LPA, LPPA, LPTR, segString, borders, prevState);
    }

    if (TR[t] == TR[t - 1] + LPTR[t])
    {
        prevState = "TRANSCRIPT";
        segString.push_back(prevState);
        funcTR(t - 1, S, L, A, PA, TR, LPS, LPL, LPA, LPPA, LPTR, segString, borders, prevState);
    }
}

/**
 * Calculate the maximum a posteriori path (backtracing) - posterioir decoding
 */
std::string getBorders(const double *LPS, const double *LPL, const double *LPA, const double *LPPA, const double *LPTR, const size_t T)
{

    double *S = new double[T];
    double *L = new double[T];
    double *A = new double[T];
    double *PA = new double[T];
    double *TR = new double[T];

    // Initialize M and E in one step, no need for fill_n
    for (size_t t = 0; t < T; ++t)
    {
        S[t] = -INFINITY;
        L[t] = -INFINITY;
        A[t] = -INFINITY;
        PA[t] = -INFINITY;
        TR[t] = -INFINITY;
    }

    double start, leader, adapter, polya, transcript;
    S[0] = 0;

    for (size_t t = 1; t < T; ++t)
    {

        // TODO compress code -> compressed code (?) : 
        start, leader, adapter, polya, transcript = -INFINITY; 
        
        start = std::max(start, S[t - 1] + LPS[t]);             // s
        leader = std::max(leader, S[t - 1] + LPL[t]);           // l1 : leave start
        leader = std::max(leader, L[t - 1] + LPL[t]);           // l2 : stay in leader
        adapter = std::max(adapter, L[t - 1] + LPA[t]);         // a1 : leave leader
        adapter = std::max(adapter, A[t - 1] + LPA[t]);         // a2 : stay in adapter
        polya = std::max(polya, A[t - 1] + LPPA[t]);            // pa1 : leader adapter
        polya = std::max(polya, PA[t - 1] + LPPA[t]);           // pa2 : stay in polyA
        transcript = std::max(transcript, PA[t - 1] + LPTR[t]); // tr1 : leave polyA
        transcript = std::max(transcript, TR[t - 1] + LPTR[t]); // tr2 : stay in trancript

        S[t] = start;
        L[t] = leader;
        A[t] = adapter;
        PA[t] = polya;
        TR[t] = transcript;
    }

    std::list<std::string> segString; // define string of most probabale states at T-1 backward
    std::vector<size_t> borders;
    segString.push_back("TRANSCRIPT"); // signal value at T - 1 pos. 100% in transcript region -> beginn recursion T - 2 onward

    funcTR(T - 1, S, L, A, PA, TR, LPS, LPL, LPA, LPPA, LPTR, segString, borders, "TRANSCRIPT");

    std::ostringstream oss;
    for (size_t i = 0; i < borders.size(); ++i)
    {
        oss << borders[i];
        if (i < borders.size() - 1)
        {
            oss << ",";
        }
    }

    delete[] S;
    delete[] L;
    delete[] A;
    delete[] PA;
    delete[] TR;

    return oss.str();
}

template <typename T> // TODO is this necessary?
void writeBorders(const std::string &save_file, const std::string &read_id, const std::vector<T> &borders)
{

    ofstream output_file(save_file, ios::app);

    if (!output_file.is_open())
    {
        cerr << "Error: Unable to open file";
        exit(EXIT_FAILURE);
    }

    output_file << read_id << ",";

    for (size_t i = 0; i < borders.size(); ++i)
    {

        output_file << borders[i];

        if (i < borders.size() - 1)
        {
            output_file << ",";
        }
        else
        {
            output_file << "\n";
        }
    }
    output_file.close();
}

/**
 * Get the signal Value from python script
 */
int main()
{

    std::cout << std::fixed << std::showpoint;
    std::cout << std::setprecision(20);

    // TODO use argparse, see dynamont

    // transition parameters :
    double s = log(0.996943171897388);
    double l1 = log(0.0030568281026119044);
    double l2 = log(0.9963280807270234);
    double a1 = log(0.003671919272976708);
    double a2 = log(0.99980542449089);
    double pa1 = log(0.0001945755091038867);
    double pa2 = log(0.9996311333837735);
    double tr1 = log(0.0003688666162265902);
    double tr2 = log(1.0);

    std::string signal_values;

    std::getline(std::cin, signal_values);

    // due to buffer error while piping
    if (signal_values.empty())
    {   
        //! diff ? : std:c* vs. printf 
        std::cerr << "no signal value are provided!";
        //printf("Error: no signal value provided.");
        return 1; // non-zero value to indicate something is wrong!
    }

    // PROCESS SIGNAL : convert string to double array
    // How many signal values are there ?  T values
    const size_t T = std::count(signal_values.begin(), signal_values.end(), ',') + 2; // len(sig) + 1

    // init a double array of T-1 elements for signal values
    double *sig = new double[T - 1];

    // put each signal value in i-position of sig
    std::string value;
    std::stringstream ss(signal_values);
    int i = 0;

    while (getline(ss, value, ','))
    {
        sig[i++] = std::stod(value);
    }

    // so far we have the signal as an array of double values in sig variable 
    // initialize Forward Backward algorithm calculation
    double *forS = new double[T];
    double *forL = new double[T];
    double *forA = new double[T];
    double *forPA = new double[T];
    double *forTR = new double[T];
    double *backS = new double[T];
    double *backL = new double[T];
    double *backA = new double[T];
    double *backPA = new double[T];
    double *backTR = new double[T];

    for (size_t t = 0; t < T; ++t)
    {
        forS[t] = -INFINITY;
        backS[t] = -INFINITY;
        forL[t] = -INFINITY;
        backL[t] = -INFINITY;
        forA[t] = -INFINITY;
        backA[t] = -INFINITY;
        forPA[t] = -INFINITY;
        backPA[t] = -INFINITY;
        forTR[t] = -INFINITY;
        backTR[t] = -INFINITY;
    }

    // calculate segmentation probabilities, fill forward matrices
    logF(sig, forS, forL, forA, forPA, forTR, T, s, l1, l2, a1, a2, pa1, pa2, tr1, tr2);
    
    // calculate segmentation probabilities, fill backward matrices
    logB(sig, backS, backL, backA, backPA, backTR, T, s, l1, l2, a1, a2, pa1, pa2, tr1, tr2);
    
    // where both values should meet each other
    const double Zf = forTR[T - 1]; // end of trancript for Forward
    const double Zb = backS[0];     // is same as beginning of start for Backward

    //! ----------------------------------------------- THE START OF MAIN CALCULATATION -----------------------------------------------

    const double *LPS = logP(forS, backS, Zf, T);
    const double *LPL = logP(forL, backL, Zf, T);
    const double *LPA = logP(forA, backA, Zf, T);
    const double *LPPA = logP(forPA, backPA, Zf, T);
    const double *LPTR = logP(forTR, backTR, Zf, T);

    std::string borders = getBorders(LPS, LPL, LPA, LPPA, LPTR, T);

    if (borders.empty())
    {
        printf("segmentation failed!");
        // always clean up before return!!
        delete[] LPS;
        delete[] LPL;
        delete[] LPA;
        delete[] LPPA;
        delete[] LPTR;
        delete[] forS;
        delete[] forL;
        delete[] forA;
        delete[] forPA;
        delete[] forTR;
        delete[] backS;
        delete[] backL;
        delete[] backA;
        delete[] backPA;
        delete[] backTR;
        delete[] sig;

        return 1;
    }

    // writes for python 
    std::cout << borders << std::endl;

    // Clean up
    delete[] LPS;
    delete[] LPL;
    delete[] LPA;
    delete[] LPPA;
    delete[] LPTR;

    delete[] forS;
    delete[] forL;
    delete[] forA;
    delete[] forPA;
    delete[] forTR;
    delete[] backS;
    delete[] backL;
    delete[] backA;
    delete[] backPA;
    delete[] backTR;
    delete[] sig;

    return 0; // if prints the border then exits with 0 !
}
