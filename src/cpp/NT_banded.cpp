// author: Jannes Spangenberg
// e-mail: jannes.spangenberg@uni-jena.de
// github: https://github.com/JannesSP
// website: https://jannessp.github.io

#include "NT_banded.hpp"

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
void logF_banded(const double *sig, const int *kmerSeq, double *M, double *E, const std::size_t T, const std::size_t B, const std::tuple<double, double> *model, const std::vector<std::tuple<long, std::size_t, std::size_t>> &bounds, const std::size_t BANDWIDTH)
{
    const double m1 = transitions_NT.at("m1");
    const double e2 = transitions_NT.at("e2");
    E[BANDWIDTH + 1] = 0;
    std::size_t tB = 0;
    long oldBStart = std::get<0>(bounds[0]);
    for (std::size_t t = 1; t < T; ++t)
    {
        tB += B;
        auto [bStart, nStart, nEnd] = bounds[t];
        // forward specific: if start == 0 move it to 1 to prevent 'if's in n-loop
        if (!nStart) [[unlikely]]
        {
            nStart = 1;
        }

        const long offset = tB - bStart + 1; // TN to TB conversion offsert
        // normal lookup: no band shift
        long mShift = -B - 1;
        long eShift = -B;

        if (oldBStart != bStart) // shift index and t != 0
        {
            // different lookup: band shifts with diagonal
            ++mShift; // const std::size_t mShift = -B;
            ++eShift; // const std::size_t e2Shift = -B + 1;
            oldBStart = bStart;
        }

#pragma omp parallel for
        for (std::size_t n = nStart; n < nEnd; ++n)
        {
            const std::size_t idx = n + offset;
            const double score = scoreKmer(sig[t - 1], kmerSeq[n - 1], model);       // Cache scoreKmer for (t-1, n-1)
            M[idx] = E[idx + mShift] + score + m1;                                   // tB - B + n - 1 : (t - 1),(n - 1)
            E[idx] = logPlus(M[idx + eShift] + score, E[idx + eShift] + score + e2); // tB - B + n (t - 1),n
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
void logB_banded(const double *sig, const int *kmerSeq, double *M, double *E, const std::size_t T, const std::size_t N, const std::size_t B, const std::tuple<double, double> *model, const std::vector<std::tuple<long, std::size_t, std::size_t>> &bounds, const std::size_t BANDWIDTH)
{
    const double m1 = transitions_NT.at("m1");
    const double e2 = transitions_NT.at("e2");
    std::size_t tB = (T - 1) * B;
    E[tB + BANDWIDTH + 1] = 0;
    long oldBStart = std::get<0>(bounds[T - 1]);
    for (std::size_t t = T - 1; t-- > 0;)
    {
        tB -= B;
        auto [bStart, nStart, nEnd] = bounds[t];

        const long offset = tB - bStart + 1; // TN to TB conversion offsert
        // normal lookup: no band shift
        long mShift = B + 1;
        long eShift = B;

        if (oldBStart != bStart) // shift index and t != 0
        {
            // different lookup: band shifts with diagonal
            --mShift; // const std::size_t mShift = +B;
            --eShift; // const std::size_t e2Shift = +B - 1;
            oldBStart = bStart;
        }

#pragma omp parallel for
        for (std::size_t n = nStart; n < nEnd; ++n)
        { // speed up, due to rules no need to look at upper triangle of matrices
            double ext = -INFINITY;
            const std::size_t idx = n + offset;
            if (n + 1 < N) [[likely]]
            {
                ext = M[idx + mShift] + scoreKmer(sig[t], kmerSeq[n], model) + m1; // tB + B + n + 1; (t+1), (n+1)
            }
            if (n > 0) [[likely]]
            {
                const double score = scoreKmer(sig[t], kmerSeq[n - 1], model);
                M[idx] = E[idx + eShift] + score;                 // transition probability is always 1, e1 first extend, tB + B + n; (t+1), n
                ext = logPlus(ext, E[idx + eShift] + score + e2); // e2 extend further
            }
            E[idx] = ext;
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
void getBorders(std::list<std::string> &segString, const double *LPM, const double *LPE, const std::size_t T, const std::size_t N, const std::size_t B, const int kmerSize, const std::vector<std::tuple<long, std::size_t, std::size_t>> &bounds, const std::size_t BANDWIDTH)
{
    const std::size_t TB = T * B;
    double *M = new double[TB];
    double *E = new double[TB];

#pragma omp parallel for
    for (std::size_t tb = 0; tb < TB; ++tb)
    {
        M[tb] = -INFINITY;
        E[tb] = -INFINITY;
    }

    E[BANDWIDTH + 1] = 0;
    std::size_t tB = 0;
    long oldBStart = std::get<0>(bounds[0]);
    for (std::size_t t = 1; t < T; ++t)
    {
        tB += B;
        auto [bStart, nStart, nEnd] = bounds[t];
        // forward specific: if start == 0 move it to 1 to prevent 'if's in n-loop
        if (!nStart) [[unlikely]]
        {
            nStart = 1;
        }
        const long offset = tB - bStart + 1; // TN to TB conversion offsert
        // normal lookup: no band shift
        long mShift = -B - 1;
        long eShift = -B;

        if (oldBStart != bStart) // shift index and t != 0
        {
            // different lookup: band shifts with diagonal
            ++mShift; // const std::size_t mShift = -B;
            ++eShift; // const std::size_t e2Shift = -B + 1;
            oldBStart = bStart;
        }

#pragma omp parallel for
        for (std::size_t n = nStart; n < nEnd; ++n)
        {
            const std::size_t idx = n + offset;
            // Compute M and E values
            M[idx] = E[idx + mShift] + LPM[idx];
            E[idx] = std::max(M[idx + eShift], E[idx + eShift]) + LPE[idx];
        }
    }
    traceback(M, E, LPM, LPE, segString, T, N, B, kmerSize, bounds, BANDWIDTH);
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
void traceback(const double *M, const double *E, const double *LPM, const double *LPE, std::list<std::string> &segString, const std::size_t T, const std::size_t N, const std::size_t B, const int kmerSize, const std::vector<std::tuple<long, std::size_t, std::size_t>> &bounds, const std::size_t BANDWIDTH)
{
    std::vector<double> segProb;
    bool isFuncM = false; // start in E

    // Iterate until the stack is empty
    std::size_t t = T - 1;
    std::size_t n = N - 1;
    std::size_t tBb = (T - 1) * B + BANDWIDTH + 1;
    long oldBStart = std::get<0>(bounds[T - 1]);
    while (t && n)
    {
        long mShift = -B - 1;
        long eShift = -B;
        if (t && oldBStart != std::get<0>(bounds[t - 1])) // shift index and t != 0
        {
            // different lookup: band shifts with diagonal
            ++mShift; // const std::size_t mShift = -B;
            ++eShift; // const std::size_t e2Shift = -B + 1;
            oldBStart = std::get<0>(bounds[t - 1]);
        }

        // A
        if (isFuncM)
        {
            // In funcM
            segProb.push_back(exp(LPM[tBb]));
            segString.push_front("M" + std::to_string(n - 1 + kmerSize / 2) + "," + std::to_string(t - 1) + "," + formattedMedian(segProb) + ";");
            segProb.clear();
            // Transition to E (with t-1, n-1)
            --t;
            --n;
            tBb += mShift;
            isFuncM = false;
        }
        // E
        else
        {
            // In funcE
            const double logScore = LPE[tBb];
            segProb.push_back(exp(logScore));
            isFuncM = (E[tBb] == M[tBb + eShift] + logScore);
            --t;
            tBb += eShift;
        }
    }
}

/**
 * @brief Computes bounds for dynamic programming in the banded HMM.
 *
 * For each time step t, this function calculates the bounds of the band
 * (lower and upper) for the dynamic programming in the banded HMM.
 * The bounds are calculated as follows:
 *   lower = midpoint - BANDWIDTH
 *   upper = midpoint + BANDWIDTH + 1
 * where midpoint = t * N / T.
 * The bounds are then adjusted so that they are within the valid range
 * of the sequence, i.e. lower >= 0 and upper <= N.
 *
 * @param T Number of time steps.
 * @param N Length of the sequence.
 * @param TB Size of the matrix used for dynamic programming.
 * @return A vector of std::tuple containing the lower and upper bounds for each time step.
 */
std::vector<std::tuple<long, std::size_t, std::size_t>> getBounds(const std::size_t T, const std::size_t N, const std::size_t TB, const std::size_t BANDWIDTH)
{
    const double NTRATIO = (double)N / T;

    std::vector<std::tuple<long, std::size_t, std::size_t>> bounds(T);

// init bounds
#pragma omp parallel for
    for (std::size_t t = 0; t < T; ++t)
    {
        const std::size_t midpoint = t * NTRATIO;
        bounds[t] = std::make_tuple(
            midpoint - BANDWIDTH,
            (midpoint >= BANDWIDTH) ? midpoint - BANDWIDTH : 0,
            (midpoint + BANDWIDTH + 1 <= N) ? midpoint + BANDWIDTH + 1 : N);
    }

    return bounds;
}