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
#include <bits/stdc++.h> // reverse strings
#include <vector>
#include <cmath> // exp
#include <assert.h>
#include <stdlib.h>
#include <algorithm>
#include "argparse.hpp"
#include "utils.hpp"

using namespace std;

// Write the backtrace functions here
void funcS(const size_t t, const double *S, const double *L, const double *A, const double *PA, const double *TR, const double *LPS, const double *LPL, const double *LPA, const double *LPPA, const double *LPTR);
void funcL(const size_t t, const double *S, const double *L, const double *A, const double *PA, const double *TR, const double *LPS, const double *LPL, const double *LPA, const double *LPPA, const double *LPTR);
void funcA(const size_t t, const double *S, const double *L, const double *A, const double *PA, const double *TR, const double *LPS, const double *LPL, const double *LPA, const double *LPPA, const double *LPTR);
void funcPA(const size_t t, const double *S, const double *L, const double *A, const double *PA, const double *TR, const double *LPS, const double *LPL, const double *LPA, const double *LPPA, const double *LPTR);
void funcTR(const size_t t, const double *S, const double *L, const double *A, const double *PA, const double *TR, const double *LPS, const double *LPL, const double *LPA, const double *LPPA, const double *LPTR);

inline constexpr double EPSILON = 1e-8;       // chose by eye just to distinguish real errors from numeric errors
double s, l1, l2, a1, a2, pa1, pa2, tr1, tr2; // transition parameters

// Asserts doubleing point compatibility at compile time
// necessary for INFINITY usage
static_assert(numeric_limits<double>::is_iec559, "IEEE 754 required");

/**
 * Calculate forward matrices using logarithmic values
 */
void logF(const double *sig, const double *S, const double *L, const double *A, const double *PA, const double *TR, const size_t T)
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
        // TODO add other PDF functions for other distributions and add model parameters
        //       accumulator + (prevV * emission                                        * transition)
        start = logPlus(start, S[t - 1] + logNormalPdf(sig[t - 1], TODO_mean, TODO_stdev) + s);

        leader = logPlus(leader, S[t - 1] + logNormalPdf(sig[t - 1], TODO_mean, TODO_stdev) + l1);
        leader = logPlus(leader, L[t - 1] + logNormalPdf(sig[t - 1], TODO_mean, TODO_stdev) + l2);

        adapter = logPlus(adapter, L[t - 1] + logNormalPdf(sig[t - 1], TODO_mean, TODO_stdev) + a1);
        adapter = logPlus(adapter, A[t - 1] + logNormalPdf(sig[t - 1], TODO_mean, TODO_stdev) + a2);

        polya = logPlus(polya, A[t - 1] + logNormalPdf(sig[t - 1], TODO_mean, TODO_stdev) + pa1);
        polya = logPlus(polya, PA[t - 1] + logNormalPdf(sig[t - 1], TODO_mean, TODO_stdev) + pa2);

        transcript = logPlus(transcript, PA[t - 1] + logNormalPdf(sig[t - 1], TODO_mean, TODO_stdev) + tr1);
        transcript = logPlus(transcript, TR[t - 1] + logNormalPdf(sig[t - 1], TODO_mean, TODO_stdev) + tr2);

        // update state matrices
        S[t] = start
            L[t] = leader;
        A[t] = adapter;
        PA[t] = polya;
        TR[t] = transcript;
    }
}

/**
 * Calculate backward matrices using logarithmic values
 */
void logB(const double *sig, const double *S, const double *L, const double *A, const double *PA, const double *TR, const size_t T)
{
    double start, leader, adapter, polya, transcript;
    TR[T - 1] = 0;
    for (size_t t = T; t-- > 0;)
    { // T-1, ..., 1, 0
        // init state accumulators
        start = -INFINITY;
        leader = -INFINITY;
        adapter = -INFINITY;
        polya = -INFINITY;
        transcript = -INFINITY;

        // calculate probabilities
        //       accumulator + (prevV * emission(t)                                     * transition)
        start = logPlus(start, S[t + 1] + logNormalPdf(sig[t], TODO_mean, TODO_stdev) + s);
        start = logPlus(start, L[t + 1] + logNormalPdf(sig[t], TODO_mean, TODO_stdev) + l1);

        leader = logPlus(leader, S[t + 1] + logNormalPdf(sig[t], TODO_mean, TODO_stdev) + l2);
        leader = logPlus(leader, A[t + 1] + logNormalPdf(sig[t], TODO_mean, TODO_stdev) + a1);

        adapter = logPlus(adapter, A[t + 1] + logNormalPdf(sig[t], TODO_mean, TODO_stdev) + a2);
        adapter = logPlus(adapter, PA[t + 1] + logNormalPdf(sig[t], TODO_mean, TODO_stdev) + pa1);

        polya = logPlus(polya, PA[t + 1] + logNormalPdf(sig[t], TODO_mean, TODO_stdev) + pa2);
        polya = logPlus(polya, TR[t + 1] + logNormalPdf(sig[t], TODO_mean, TODO_stdev) + tr1);

        transcript = logPlus(transcript, TR[t + 1] + logNormalPdf(sig[t], TODO_mean, TODO_stdev) + tr2);

        // update state matrices
        S[t] = start
            L[t] = leader;
        A[t] = adapter;
        PA[t] = polya;
        TR[t] = transcript;
    }
}

/**
 * Calculate the logarithmic probability matrix
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

/**
 * Calculate the maximum a posteriori path (backtracing)
 */
list<string> getBorders(const double *LPS, const double *LPL, const double *LPA, const double *LPPA, const double *LPTR, const size_t T)
{
    double *S = new double[T];
    double *L = new double[T];
    double *A = new double[T];
    double *PA = new double[T];
    double *TR = new double[T];

    // Initialize M and E in one step, no need for fill_n
    for (size_t t = 0; t < T; ++i)
    {
        S[i] = -INFINITY;
        L[i] = -INFINITY;
        A[i] = -INFINITY;
        PA[i] = -INFINITY;
        TR[i] = -INFINITY;
    }

    double start, leader, adapter, polya, transcript;
    S[0] = 0;
    for (size_t t = 1; t < T; ++t)
    {
        start = -INFINITY;
        leader = -INFINITY;
        adapter = -INFINITY;
        polya = -INFINITY;
        transcript = -INFINITY;

        start = max(start, S[t - 1] + LPS[t]);             // s
        leader = max(leader, S[t - 1] + LPL[t]);           // l1
        leader = max(leader, L[t - 1] + LPL[t]);           // l2
        adapter = max(adapter, L[t - 1] + LPA[t]);         // a1
        adapter = max(adapter, A[t - 1] + LPA[t]);         // a2
        polya = max(polya, A[t - 1] + LPPA[t]);            // pa1
        polya = max(polya, PA[t - 1] + LPPA[t]);           // pa2
        transcript = max(transcript, PA[t - 1] + LPTR[t]); // tr1
        transcript = max(transcript, TR[t - 1] + LPTR[t]); // tr2

        S[t] = start;
        L[t] = leader;
        A[t] = adapter;
        PA[t] = polya;
        TR[t] = transcript;
    }
    list<string> segString;
    funcTR(T - 1, S, L, A, PA, TR, LPS, LPL, LPA, LPPA, LPTR);
    delete[] S;
    delete[] L;
    delete[] A;
    delete[] PA;
    delete[] TR;
    return segString;
}

void funcS(const size_t t, const double *S, const double *L, const double *A, const double *PA, const double *TR, const double *LPS, const double *LPL, const double *LPA, const double *LPPA, const double *LPTR)
{
    // TODO
}

void funcL(const size_t t, const double *S, const double *L, const double *A, const double *PA, const double *TR, const double *LPS, const double *LPL, const double *LPA, const double *LPPA, const double *LPTR)
{
    // TODO
}

void funcA(const size_t t, const double *S, const double *L, const double *A, const double *PA, const double *TR, const double *LPS, const double *LPL, const double *LPA, const double *LPPA, const double *LPTR)
{
    // TODO
}

void funcPA(const size_t t, const double *S, const double *L, const double *A, const double *PA, const double *TR, const double *LPS, const double *LPL, const double *LPA, const double *LPPA, const double *LPTR)
{
    // TODO
}

void funcTR(const size_t t, const double *S, const double *L, const double *A, const double *PA, const double *TR, const double *LPS, const double *LPL, const double *LPA, const double *LPPA, const double *LPTR)
{
    // TODO
}

/**
 * Train transition parameter with baum welch algorithm
 */
tuple<double, double, double, double, double, double, double, double, double> trainTransition(const double *sig, const double *forS, const double *forL, const double *forA, const double *forPA, const double *forTR, const double *backS, const double *backL, const double *backA, const double *backPA, const double *backTR, const size_t T)
{
    // Transition parameters
    double newS = -INFINITY, newL1 = -INFINITY, newL2 = -INFINITY, newA1 = -INFINITY, newA2 = -INFINITY, newPA1 = -INFINITY, newPA2 = -INFINITY, newTR1 = -INFINITY, newTR2 = -INFINITY;

    for (size_t t = 0; t < T; ++t)
    {
        // TODO
    }
    // TODO average over the number of transitions

    return tuple<double, double, double, double, double, double, double, double, double>({exp(newS), exp(newL1), exp(newL2), exp(newA1), exp(newA2), exp(newPA1), exp(newPA2), exp(newTR1), exp(newTR2)});
}

void trainParams(const double *sig, double *forS, double *forL, double *forA, double *forPA, double *forTR, double *backS, double *backL, double *backA, double *backPA, double *backTR, const size_t T)
{
    auto [newS, newL1, newL2, newA1, newA2, newPA1, newPA2, newTR1, newTR2] = trainTransition(sig, forS, forL, forA, forPA, forTR, backS, backL, backA, backPA, backTR, T);

    cout << "s:" << newS << ";l1:" << newL1 << ";l2:" << newL2 << ";a1:" << newA1 << ";a2:" << newA2 << ";pa1:" << newPA1 << ";pa2:" << newPA2 << ";tr1:" << newTR1 << ";tr2:" << newTR2 << endl;
}

/**
 * Read signal and read from stdin until the TERM_STRING is seen
 */
int main(int argc, char *argv[])
{
    int pore;
    bool train, calcZ, prob;
    string modelpath;
    const string TERM_STRING = "$";

    // Argparser
    argparse::ArgumentParser program("dynamont basic", "0.1");
    // TODO parameters for DP
    // double s, l1, l2, a1, a2, pa1, pa2, tr1, tr2; // transition parameters
    program.add_argument("-s", "--startscore").help("Transition probability for start rule").default_value(1.00).scan<'g', double>().store_into(s);        // s
    program.add_argument("-l1", "--leaderscore1").help("Transition probability for leader rule 1").default_value(1.00).scan<'g', double>().store_into(l1); // l1
    // TODO rest

    program.add_argument("-t", "--train").help("Switch algorithm to transition and emission parameter training mode").default_value(false).implicit_value(true).store_into(train);
    program.add_argument("-z", "--calcZ").help("Switch algorithm to only calculate Z").default_value(false).implicit_value(true).store_into(calcZ);
    program.add_argument("-r", "--pore").help("Pore generation used to sequence the data").default_value(9).choices(9, 10).scan<'i', int>().store_into(pore);
    program.add_argument("-p", "--probabilty").help("Print out the segment border probability").default_value(false).implicit_value(true).store_into(prob);

    try
    {
        program.parse_args(argc, argv);
    }
    catch (const std::runtime_error &err)
    {
        cerr << err.what() << std::endl;
        cerr << program;
        return 1;
    }

    C = program.get<int>("minSegLen");

    if (pore == 9)
    {
        kmerSize = 5;
    }
    else if (pore == 10)
    {
        kmerSize = 9;
    }
    string signal;
    int truish = 1;

    while (truish)
    {
        // truish = 0;
        // echo 107,107,107.2,108.0,108.9,111.2,105.7,104.3,107.1,105.7,105,105 CAAAAA| src\segment.exe
        // read input, signal and read whitespace separated in single line
        getline(cin, signal);

        // break loop if termination character ...
        if (signal.find(TERM_STRING) != string::npos)
        {
            return 0;
            // ... or signal or read is missing
        }
        else if (signal.empty())
        {
            cout << "Signal missing!\n";
            cout.flush();
            return 1;
        }

        // cerr<<"DEBUG 1"<<endl;
        // process signal: convert string to double array
        T = count(signal.begin(), signal.end(), ',') + 2; // len(sig) + 1
        double *sig = new double[T - 1];
        fill_n(sig, T - 1, -INFINITY);
        string value;
        stringstream ss(signal);
        int i = 0;
        while (getline(ss, value, ','))
        {
            sig[i++] = stod(value);
        }

        // initialize matrices
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
            forS[i] = -INFINITY;
            backS[i] = -INFINITY;
            forL[i] = -INFINITY;
            backL[i] = -INFINITY;
            forA[i] = -INFINITY;
            backA[i] = -INFINITY;
            forPA[i] = -INFINITY;
            backPA[i] = -INFINITY;
            forTR[i] = -INFINITY;
            backTR[i] = -INFINITY;
        }

        // calculate segmentation probabilities, fill forward matrices
        logF(sig, forS, forL, forA, forPA, forTR, T);
        // calculate segmentation probabilities, fill backward matrices
        logB(sig, backS, backL, backA, backPA, backTR, T);

        const double Zf = forTR[T - 1];
        const double Zb = backS[0];

        // Numeric error is scaled by input size, Z in forward and backward should match by some numeric error EPSILON
        if (abs(Zf - Zb) / T > EPSILON)
        {
            cerr << fixed << showpoint;
            cerr << setprecision(20);
            cerr << "Z values between matrices do not match! Zf: " << Zf << ", Zb: " << Zb << ", " << abs(Zf - Zb) / T << " > " << EPSILON << endl;
            cerr.flush();
            exit(11);
        }

        if (calcZ)
        {
            cout << Zf << "\n";
            cout.flush();
        }
        else
        {

            // train both Transitions and Emissions
            if (train)
            {
                trainParams(sig, forS, forL, forA, forPA, forTR, backS, backL, backA, backPA, backTR, T);
                cout << "Z:" << Zb << endl;
                cout.flush();
            }
            else
            {
                const double *LPS = logP(forS, backS, Zf, T);
                const double *LPL = logP(forL, backL, Zf, T);
                const double *LPA = logP(forA, backA, Zf, T);
                const double *LPPA = logP(forPA, backPA, Zf, T);
                const double *LPTR = logP(forTR, backTR, Zf, T);
                list<string> segString = getBorders(LPS, LPL, LPA, LPPA, LPTR, T);

                for (auto const &seg : segString)
                {
                    cout << seg;
                }
                cout << endl;
                cout.flush();

                // calculate sum of segment probabilities
                if (prob)
                {
                    for (size_t t = 0; t < T; ++t)
                    {
                        cout << LPPA[t] << ",";
                    }
                    cout << endl;
                    cout.flush();
                }

                // Clean up
                delete[] LPS;
                delete[] LPL;
                delete[] LPA;
                delete[] LPPA;
                delete[] LPTR;
            }
        }
        // Clean up
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
    }
    return 0;
}
