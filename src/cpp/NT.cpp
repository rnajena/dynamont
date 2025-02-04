// author: Jannes Spangenberg
// e-mail: jannes.spangenberg@uni-jena.de
// github: https://github.com/JannesSP
// website: https://jannessp.github.io

#include "NT.hpp"

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
void logF_NT(const double *sig, const int *kmerSeq, double *M, double *E, const std::size_t T, const std::size_t N, const std::tuple<double, double> *model)
{
    const double m1 = transitions_NT.at("m1");
    const double e2 = transitions_NT.at("e2");
    E[0] = 0;
    std::size_t tN = 0;
    for (std::size_t t = 1; t < T; ++t)
    {
        tN += N;
#pragma omp parallel for
        for (std::size_t n = 1; n < N; ++n)
        {                                                                      // speed up, due to rules no need to look at upper triangle of matrices
            const double score = scoreKmer(sig[t - 1], kmerSeq[n - 1], model); // Cache scoreKmer for (t-1, n-1)
            M[tN + n] = E[tN - N + n - 1] + score + m1;
            E[tN + n] = logPlus(M[tN - N + n] + score, E[tN - N + n] + score + e2);
        }
    }
}

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
void logB_NT(const double *sig, const int *kmerSeq, double *M, double *E, const std::size_t T, const std::size_t N, const std::tuple<double, double> *model)
{
    const double m1 = transitions_NT.at("m1");
    const double e2 = transitions_NT.at("e2");
    E[T * N - 1] = 0;
    std::size_t tN = (T - 1) * N;
    for (std::size_t t = T - 1; t-- > 0;)
    {
        tN -= N;
#pragma omp parallel for
        for (std::size_t n = 0; n < N; ++n)
        { // speed up, due to rules no need to look at upper triangle of matrices
            double ext = -INFINITY;
            if (n + 1 < N) [[likely]]
            {
                ext = M[tN + N + n + 1] + scoreKmer(sig[t], kmerSeq[n], model) + m1;
            }
            if (n > 0) [[likely]]
            {
                const double score = scoreKmer(sig[t], kmerSeq[n - 1], model);
                M[tN + n] = E[tN + N + n] + score;              // transition probability is always 1, e1 first extend
                ext = logPlus(ext, E[tN + N + n] + score + e2); // e2 extend further
            }
            E[tN + n] = ext;
        }
    }
}

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
void getBorders(std::list<std::string> &segString, const double *LPM, const double *LPE, const std::size_t T, const std::size_t N, const int kmerSize)
{
    const std::size_t TN = T * N;
    double *M = new double[TN];
    double *E = new double[TN];
#pragma omp parallel for
    for (std::size_t i = 0; i < TN; ++i)
    {
        M[i] = -INFINITY;
        E[i] = -INFINITY;
    }
    E[0] = 0;
    std::size_t tN = 0;
    for (std::size_t t = 1; t < T; ++t)
    {
        tN += N;
#pragma omp parallel for
        for (std::size_t n = 1; n < N; ++n)
        {
            const std::size_t prevIdx = tN - N + n;
            const std::size_t currentIdx = tN + n;
            // Compute M and E values
            M[currentIdx] = E[prevIdx - 1] + LPM[currentIdx];                                     // m1
            E[currentIdx] = std::max(M[prevIdx] + LPE[currentIdx], E[prevIdx] + LPE[currentIdx]); // e1, e2
        }
    }
    // std::vector<double> segProb;
    // funcE(T - 1, N - 1, M, E, LPM, LPE, segString, N, segProb, kmerSize);
    traceback(T - 1, N - 1, M, E, LPM, LPE, segString, N, kmerSize);
    delete[] M;
    delete[] E;
}

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
void traceback(std::size_t t, std::size_t n, const double *M, const double *E, const double *LPM, const double *LPE, std::list<std::string> &segString, const std::size_t N, const int kmerSize)
{
    std::vector<double> segProb;
    bool isFuncM = false; // start in E

    // Iterate until the stack is empty
    while (t && n)
    {
        std::size_t curIdx = t * N + n;
        // A
        if (isFuncM)
        {
            // In funcM
            segProb.push_back(exp(LPM[curIdx]));
            segString.push_front("M" + std::to_string(n - 1 + kmerSize / 2) + "," + std::to_string(t - 1) + "," + formattedMedian(segProb) + ";");
            segProb.clear();
            // Transition to E (with t-1, n-1)
            --t;
            --n;
            isFuncM = false;
        }
        // E
        else
        {
            // In funcE
            const double logScore = LPE[curIdx];
            segProb.push_back(exp(logScore));
            isFuncM = (E[curIdx] == M[curIdx - N] + logScore);
            --t;
        }
    }
}

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
std::tuple<double, double, double> trainTransition(const double *sig, const int *kmerSeq, const double *forE, const double *backM, const double *backE, const std::size_t T, const std::size_t N, const std::tuple<double, double> *model)
{
    // Transition parameters
    double newM1 = -INFINITY, newE1 = 0, newE2 = -INFINITY;

    for (std::size_t t = 0; t < T - 1; ++t)
    {
        for (std::size_t n = 0; n <= t && n < N; ++n)
        { // speed up, due to rules no need to look at upper triangle of matrices
            if (n + 1 < N) [[likely]]
            {
                // m1:                 forward(i)        a                      e(i+1)                                 backward(i+1)
                newM1 = logPlus(newM1, forE[t * N + n] + transitions_NT.at("m1") + scoreKmer(sig[t], kmerSeq[n], model) + backM[(t + 1) * N + (n + 1)]);
            }

            if (n > 0) [[likely]]
            {
                double score = scoreKmer(sig[t], kmerSeq[n - 1], model);
                // newE1 = logPlus(newE1, forM[t*N+n] + transitions_NT.at("e1") + score + backE[(t+1)*N+n]);
                newE2 = logPlus(newE2, forE[t * N + n] + transitions_NT.at("e2") + score + backE[(t + 1) * N + n]);
            }
        }
    }
    // average over the number of transitions
    // double Am = newE1;
    // newE1 = newE1 - Am;
    const double Ae = logPlus(newE2, newM1);
    if (!std::isinf(Ae))
    {
        newM1 = newM1 - Ae;
        newE2 = newE2 - Ae;
    }
    return std::make_tuple(
        exp(newM1),
        exp(newE1),
        exp(newE2));
}

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
std::tuple<double *, double *> trainEmission(const double *sig, const int *kmerSeq, const double *forM, const double *forE, const double *backM, const double *backE, const std::size_t T, const std::size_t N, const int numKmers)
{

    // TODO: calculate this with LP matrices
    // see notes

    // Emission
    // https://courses.grainger.illinois.edu/ece417/fa2021/lectures/lec15.pdf
    // https://f.hubspotusercontent40.net/hubfs/8111846/Unicon_October2020/pdf/bilmes-em-algorithm.pdf
    // gamma_t(i) is the probability of being in state i at time t
    // gamma for state M - expected number of transitions of M at given time (T) for all latent states (kmers)
    // 0-Initialize memory
    const std::size_t TN = T * N;
    double *G = new double[TN]();
    double *kmers = new double[N]();
    double *d = new double[N]();
    double *means = new double[numKmers]();
    double *stdevs = new double[numKmers]();
    int *counts = new int[numKmers]();

    // initialize G
    for (std::size_t t = 0; t < T; ++t)
    {
        // calibrate with the sum of transitions
        double s = -INFINITY;
        for (std::size_t n = 0; n < N; ++n)
        { // speed up, due to rules no need to look at upper triangle of matrices
            G[t * N + n] = logPlus(forM[t * N + n] + backM[t * N + n], forE[t * N + n] + backE[t * N + n]);
            s = logPlus(s, G[t * N + n]);
        }
        for (std::size_t n = 0; n < N; ++n)
        {
            if (!std::isinf(s))
            {
                G[t * N + n] -= s;
            }
        }
    }

    for (std::size_t n = 1; n < N; ++n)
    {
        counts[kmerSeq[n - 1]]++;
        for (std::size_t t = n; t < T; ++t)
        {                                 // speed up, due to rules no need to look at upper triangle of matrices
            double g = exp(G[t * N + n]); // Cache the exp(G[t * N + n])
            kmers[n] += g * sig[t - 1];   // Accumulate for kmers
            d[n] += g;                    // Accumulate for normalizer
        }
        if (d[n] != 0)
        {
            kmers[n] = kmers[n] / d[n]; // Normalize
        }
    }

    for (std::size_t n = 1; n < N; ++n)
    {
        means[kmerSeq[n - 1]] += kmers[n] / counts[kmerSeq[n - 1]]; // Update means
    }

    // Emission (stdev of kmers)
    std::fill_n(kmers, N, 0);
    for (std::size_t n = 1; n < N; ++n)
    {
        for (std::size_t t = n; t < T; ++t)
        {                                                     // n<=t, speed up, due to rules no need to look at upper triangle of matrices
            double diff = sig[t - 1] - means[kmerSeq[n - 1]]; // Compute difference from mean
            kmers[n] += exp(G[t * N + n]) * diff * diff;      // Accumulate for variance
        }
        if (d[n] > 0)
        {
            kmers[n] = kmers[n] / d[n];
        }
        stdevs[kmerSeq[n - 1]] += kmers[n] / counts[kmerSeq[n - 1]];
    }

#pragma omp parallel for
    for (std::size_t i = 0; i < numKmers; ++i)
    {
        // transform vars to stdevs
        stdevs[i] = sqrt(stdevs[i]);
    }

    delete[] G;
    delete[] kmers;
    delete[] counts;
    delete[] d;
    return std::tuple<double *, double *>({means, stdevs});
}

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
void trainParams(const double *sig, const int *kmerSeq, const double *forM, const double *forE, const double *backM, const double *backE, const std::size_t T, const std::size_t N, std::tuple<double, double> *model, const int alphabetSize, const int numKmers, const int kmerSize)
{
    const auto [newM, newE1, newE2] = trainTransition(sig, kmerSeq, forE, backM, backE, T, N, model);
    std::cout << "m1:" << newM << ";e1:" << newE1 << ";e2:" << newE2 << std::endl;
    const auto [newMeans, newStdevs] = trainEmission(sig, kmerSeq, forM, forE, backM, backE, T, N, numKmers);
    for (int i = 0; i < numKmers; i++)
    {
        if (newStdevs[i] != 0.0) [[likely]]
        {
            std::cout << itoa(i, alphabetSize, kmerSize, rna) << ":" << newMeans[i] << "," << newStdevs[i] << ";";
        }
    }
    std::cout << std::endl;
    delete[] newMeans;
    delete[] newStdevs;
}