// author: Jannes Spangenberg
// e-mail: jannes.spangenberg@uni-jena.de
// github: https://github.com/JannesSP
// website: https://jannessp.github.io

#include "NTC.hpp"

std::size_t TNK, NK;
double ppTNm, ppTNe, ppTKm, ppTKe;
int alphabetSize, kmerSize, halfKmerSize, stepSize;

/**
 * Computes the posterior probabilities for a hidden Markov model (HMM) using
 * forward and backward probabilities.
 *
 * @param logAPSEI Posterior probabilities stored in an unordered map.
 * @param forAPSEI Forward probabilities stored in an unordered map.
 * @param backAPSEI Backward probabilities stored in an unordered map.
 * @param Z Normalization constant for the forward and backward probabilities.
 * @param allowedKeys Array of allowed indices for processing.
 *
 * @return Posterior probabilities for each state in the HMM.
 */
void logP(std::unordered_map<std::size_t, std::array<dproxy, NUMMAT>> &logAPSEI, const std::unordered_map<std::size_t, std::array<dproxy, NUMMAT>> &forAPSEI, const std::unordered_map<std::size_t, std::array<dproxy, NUMMAT>> &backAPSEI, const double Z, const std::vector<std::size_t> &allowedKeys)
{
    // #pragma omp parallel for
    for (const std::size_t &tnk : allowedKeys)
    {
        // Lookup in std::unordered_map is done once per index to avoid repetitive lookups
        auto &logAPSEIRef = logAPSEI[tnk];            // index must not exist
        const auto &forAPSEIRef = forAPSEI.at(tnk);   // .at(idx) : index must exist
        const auto &backAPSEIRef = backAPSEI.at(tnk); // .at(idx) : index must exist
#pragma omp parallel for
        for (int mat = 0; mat < NUMMAT; ++mat)
        {
            logAPSEIRef[mat] = forAPSEIRef[mat] + backAPSEIRef[mat] - Z;
        }
    }
}

/**
 * Calculate the log probability for a given signal index s of being in one of two states, M or E.
 *
 * @param LP pointer to the output array of log probabilities of size S
 * @param forM pointer to the array of forward probabilities of state M of size S
 * @param backM pointer to the array of backward probabilities of state M of size S
 * @param forE pointer to the array of forward probabilities of state E of size S
 * @param backE pointer to the array of backward probabilities of state E of size S
 * @param S the size of the arrays
 * @param Z the alignment score
 */
void logP(double *LP, const dproxy *forM, const dproxy *backM, const dproxy *forE, const dproxy *backE, const std::size_t S, const double Z)
{
#pragma omp parallel for
    for (std::size_t s = 0; s < S; ++s)
    {
        const double valM = forM[s] + backM[s] - Z;
        const double valE = forE[s] + backE[s] - Z;
        LP[s] = logPlus(valM, valE);
    }
}

// ===============================================================
// ===============================================================
// ======================== PREPROCESSING ========================
// ===============================================================
// ===============================================================

/**
 * Computes forward probabilities for a hidden Markov model (HMM) using
 * signal data, k-mer sequence, and transition probabilities.
 *
 * @param sig Pointer to the signal data array.
 * @param kmerSeq Pointer to the k-mer sequence array.
 * @param M Matrix of match probabilities.
 * @param E Matrix of extension probabilities.
 * @param T Length of the ONT raw signal + 1.
 * @param N Length of nucleotide sequence + 1.
 * @param model Array containing tuples of model parameters.
 */
void ppForTN(const double *sig, const int *kmerSeq, dproxy *M, dproxy *E, const std::size_t T, const std::size_t N, const std::tuple<double, double> *model)
{
    E[0] = 0;
    std::size_t tN = 0;
    for (std::size_t t = 1; t < T; ++t)
    {
        tN += N;
#pragma omp parallel for
        for (std::size_t n = 1; n < N; ++n)
        {
            const double sc = scoreKmer(sig[t - 1], kmerSeq[n - 1], model);
            M[tN + n] = E[tN - N + n - 1] + sc + ppTNm;
            E[tN + n] = logPlus(M[tN - N + n] + sc, E[tN - N + n] + sc + ppTNe);
        }
    }
}

/**
 * Computes backward probabilities for a hidden Markov model (HMM) using
 * signal data and transition probabilities.
 *
 * @param sig Pointer to the signal data array.
 * @param M Matrix of match probabilities.
 * @param E Matrix of extension probabilities.
 * @param T Length of the ONT raw signal + 1.
 * @param N Length of nucleotide sequence + 1.
 * @param model Array containing tuples of model parameters.
 */
void ppBackTN(const double *sig, const int *kmerSeq, dproxy *M, dproxy *E, const std::size_t T, const std::size_t N, const std::tuple<double, double> *model)
{
    E[T * N - 1] = 0;
    std::size_t tN = (T - 1) * N;
    for (std::size_t t = T - 1; t-- > 0;)
    { // iterates from T-2 to 0
        tN -= N;
#pragma omp parallel for
        for (std::size_t n = 0; n < N; ++n)
        {
            double ext = -INFINITY;
            if (n + 1 < N) [[likely]]
            {
                ext = M[tN + N + n + 1] + scoreKmer(sig[t], kmerSeq[n], model) + ppTNm;
            }
            if (n > 0) [[likely]]
            {
                const double sc = scoreKmer(sig[t], kmerSeq[n - 1], model);
                M[tN + n] = E[tN + N + n] + sc;                 // e1 first extend
                ext = logPlus(ext, E[tN + N + n] + sc + ppTNe); // e2 extend further
            }
            E[tN + n] = ext;
        }
    }
}

/**
 * Computes forward probabilities for a hidden Markov model (HMM) using signal
 * data and transition probabilities.
 *
 * @param sig Pointer to the signal data array.
 * @param M Matrix of match probabilities.
 * @param E Matrix of extension probabilities.
 * @param T Length of the ONT raw signal + 1.
 * @param K The number of k-mers.
 * @param model Array containing tuples of model parameters.
 */
void ppForTK(const double *sig, dproxy *M, dproxy *E, const std::size_t T, const std::size_t K, const std::tuple<double, double> *model)
{
    // init first column with log(1.0)
    for (std::size_t k = 0; k < K; ++k)
    {
        E[k] = 0.0;
    }
    std::size_t tK = 0;
    for (std::size_t t = 1; t < T; ++t)
    {
        const std::size_t prevTK = tK; // (t-1)*K
        tK += K;
        // #pragma omp parallel for
        for (std::size_t k = 0; k < K; ++k)
        {
            double mat = -INFINITY;
            const double sc = scoreKmer(sig[t - 1], k, model);
            for (std::size_t preKmer = precessingKmer(k, 0, stepSize, alphabetSize); preKmer < K; preKmer += stepSize)
            {
                mat = logPlus(mat, E[prevTK + preKmer] + sc + ppTKm);
            }
            M[tK + k] = mat;
            E[tK + k] = logPlus(M[prevTK + k] + sc, E[prevTK + k] + sc + ppTKe);
        }
    }
}

/**
 * Computes backward probabilities for a hidden Markov model (HMM) using
 * signal data and transition probabilities.
 *
 * This function iteratively calculates the backward matrices M and E,
 * which represent match and extend probabilities, respectively. It uses
 * logarithmic probabilities to efficiently process the signal data and
 * k-mer sequences, updating the probabilities based on transition
 * parameters and emission scores.
 *
 * @param sig Pointer to the signal data array.
 * @param M Matrix of match probabilities.
 * @param E Matrix of extension probabilities.
 * @param T Length of the ONT raw signal + 1.
 * @param K The number of k-mers.
 * @param model Array containing tuples of model parameters.
 */
void ppBackTK(const double *sig, dproxy *M, dproxy *E, const std::size_t T, const std::size_t K, const std::tuple<double, double> *model)
{
    // init last column with log(1.0)
    for (std::size_t k = 0; k < K; ++k)
    {
        E[(T - 1) * K + k] = 0.0;
    }
    std::size_t tK = (T - 1) * K;
    for (std::size_t t = T - 1; t-- > 0;)
    {
        const std::size_t nexttK = tK; // (t+1)*K
        tK -= K;
#pragma omp parallel for
        for (std::size_t k = 0; k < K; ++k)
        {
            double ext = -INFINITY;

            const std::size_t startKmer = successingKmer(k, 0, stepSize, alphabetSize);
            const std::size_t endKmer = startKmer + alphabetSize;
            for (std::size_t sucKmer = startKmer; sucKmer < endKmer; ++sucKmer)
            {
                ext = logPlus(ext, M[nexttK + sucKmer] + scoreKmer(sig[t], sucKmer, model) + ppTKm);
            }
            const double sc = scoreKmer(sig[t], k, model);
            M[tK + k] = E[nexttK + k] + sc;
            E[tK + k] = logPlus(ext, E[nexttK + k] + sc + ppTKe); // e2 extend further
        }
    }
}

/**
 * Compute a sparse representation of the posterior probability matrix for a given signal and k-mer sequence.
 *
 * @param sig Pointer to the ONT raw signal array with pA values.
 * @param kmerSeq Pointer to the k-mer sequence array.
 * @param tnMap Unordered map to store the sparse representation of the posterior probability matrix.
 * @param T Length of the ONT raw signal + 1.
 * @param N Length of nucleotide sequence + 1.
 * @param model Array containing kmers as keys and (mean, stdev) tuples as values.
 */
void preProcTN(const double *sig, const int *kmerSeq, std::unordered_map<std::size_t, std::unordered_set<std::size_t>> &tnMap, const std::size_t T, const std::size_t N, const std::tuple<double, double> *model)
{
    // Calculate forward and backward matrices
    const std::size_t TN = T * N;
    dproxy *forM = new dproxy[TN];
    dproxy *forE = new dproxy[TN];
    ppForTN(sig, kmerSeq, forM, forE, T, N, model);
    dproxy *backM = new dproxy[TN];
    dproxy *backE = new dproxy[TN];
    ppBackTN(sig, kmerSeq, backM, backE, T, N, model);

    const double Zf = forE[TN - 1];
    const double Zb = backE[0];

    // Numeric error is scaled by input size, Z in forward and backward should match by some numeric error EPSILON
    if (abs(Zf - Zb) / TN > EPSILON || std::isinf(Zf) || std::isinf(Zb))
    {
        std::cerr << "Z values of preProcTN matrices do not match! Zf: " << Zf << ", Zb: " << Zb << ", " << abs(Zf - Zb) / TN << " > " << EPSILON << std::endl;
        exit(1);
    }

    double *LP = new double[TN];
    logP(LP, forM, backM, forE, backE, TN, Zf);

    // extract indices with highest probability per column, until SPARSETHRESHOLD is reached
    std::size_t tN = -N;
    for (std::size_t t = 0; t < T; ++t)
    {
        tN += N;
        double s = -INFINITY;
        // get indices of values in descending order
        for (const std::size_t &n : columnArgsort(LP, N, t))
        {
            // collect key pairs with highest value per column t
            tnMap[t].insert(n);
            s = logPlus(s, LP[tN + n]);
            // stop if threshold is reached
            if (s > SPARSETHRESHOLD)
            {
                break;
            }
        }
    }
    // std::cerr<<"TN dense: "<<c/double(TN)<<"\n";

    // deallocate memory
    delete[] forM;
    delete[] forE;
    delete[] backM;
    delete[] backE;
    delete[] LP;
}

/**
 * Compute a sparse representation of the posterior probability matrix for a given signal and kmer sequence.
 *
 * @param sig Pointer to the ONT raw signal array with pA values.
 * @param tkMap Unordered map to store the sparse representation of the posterior probability matrix.
 * @param T Length of the ONT raw signal + 1.
 * @param K Length of kmer sequence + 1.
 * @param model Array containing kmers as keys and (mean, stdev) tuples as values.
 */
void preProcTK(const double *sig, std::unordered_map<std::size_t, std::unordered_set<std::size_t>> &tkMap, const std::size_t T, const std::size_t K, const std::tuple<double, double> *model)
{
    // Allocate memory in one go to reduce overhead
    const std::size_t TK = T * K;
    dproxy *forM = new dproxy[TK];
    dproxy *forE = new dproxy[TK];
    ppForTK(sig, forM, forE, T, K, model);
    dproxy *backM = new dproxy[TK];
    dproxy *backE = new dproxy[TK];
    ppBackTK(sig, backM, backE, T, K, model);

    double Zf = -INFINITY;
    double Zb = -INFINITY;
    for (std::size_t k = 0; k < K; ++k)
    {
        Zf = logPlus(Zf, forE[TK - 1 - k]);
        Zb = logPlus(Zb, backE[k]);
    }

    // Numeric error is scaled by input size, Z in forward and backward should match by some numeric error EPSILON
    if (abs(Zf - Zb) / TK > EPSILON || std::isinf(Zf) || std::isinf(Zb))
    {
        std::cerr << "Z values of preProcTK matrices do not match! Zf: " << Zf << ", Zb: " << Zb << ", " << abs(Zf - Zb) / TK << " > " << EPSILON << std::endl;
        exit(2);
    }

    // std::cerr<<"preProcTK: Zf: "<<Zf<<", Zb: "<<Zb<<", "<<abs(Zf-Zb)/(TK)<<" <! "<<EPSILON<<"\n";
    double *LP = new double[TK];
    logP(LP, forM, backM, forE, backE, TK, Zb);

    // extract indices with highest probability per column, until SPARSETHRESHOLD is reached
    std::size_t tK = -K;
    for (std::size_t t = 0; t < T; ++t)
    {
        tK += K;
        double s = -INFINITY;
        // get indices of values in descending order
        for (const std::size_t &k : columnArgsort(LP, K, t))
        {
            // collect key pairs with highest value per column t
            // allowedKeys.push_back(t*K+k);
            tkMap[t].insert(k);
            s = logPlus(s, LP[tK + k]);
            // stop if threshold is reached
            if (s >= SPARSETHRESHOLD) [[likely]] // majority of probability should be distributed to only a few keys
            {
                break;
            }
        }
    }
    // std::cerr<<"TK dense: "<<c/double(TK)<<"\n";

    // deallocate memory
    delete[] forM;
    delete[] forE;
    delete[] backM;
    delete[] backE;
    delete[] LP;
}

/**
 * @brief Preprocesses the signal and k-mer sequence to compute allowed k-mer keys.
 *
 * This function performs preprocessing on the ONT raw signal and k-mer sequence
 * to compute a set of allowed k-mer keys. It combines preprocessing results
 * from two partial 2D problems: signal vs nucleotide sequence (TN) and signal
 * vs k-mer sequence (TK). The allowed keys are useful for further analysis
 * such as resquiggles or error corrections.
 *
 * @param sig Pointer to the ONT raw signal array with pA values.
 * @param kmerSeq Pointer to the k-mer sequence array.
 * @param T Length of the ONT raw signal + 1.
 * @param N Length of nucleotide sequence + 1.
 * @param K Number of k-mers.
 * @param model Array containing tuples of model parameters for each k-mer.
 * @return Pointer to the array of allowed k-mer keys.
 */
std::vector<std::size_t> preProcTNK(const double *sig, const int *kmerSeq, const std::size_t T, const std::size_t N, const std::size_t K, const std::tuple<double, double> *model)
{
    std::vector<std::size_t> allowedKeys;
    // perform preprocessing on partial 2d problems
    std::unordered_map<std::size_t, std::unordered_set<std::size_t>> tnMap;
    preProcTN(sig, kmerSeq, tnMap, T, N, model);
    std::unordered_map<std::size_t, std::unordered_set<std::size_t>> tkMap;
    preProcTK(sig, tkMap, T, K, model);

    // combine both preprocessings using AND
    for (std::size_t t = 0; t < T; ++t)
    {
        const std::size_t tNK = t * NK;
        for (const std::size_t &n : tnMap.at(t))
        {
            // always enable the possibility to just use the read as a baseline for segmentation
            const std::size_t tNKnK = tNK + n * K;
            allowedKeys.push_back(tNKnK + kmerSeq[n - 1]);
            for (const std::size_t &k : tkMap.at(t))
            {
                // these positions will be possible resquiggles / error corrections
                allowedKeys.push_back(tNKnK + k);
            }
        }
    }

    // erase duplicates
    std::sort(allowedKeys.begin(), allowedKeys.end());
    allowedKeys.erase(std::unique(allowedKeys.begin(), allowedKeys.end()), allowedKeys.end());
    return allowedKeys;
}

// ===============================================================
// ========================== ALGORITHM ==========================
// ===============================================================

/**
 * Computes forward probabilities for a hidden Markov model (HMM) using signal
 * data, k-mer sequences, and transition probabilities.
 *
 * @param sig Pointer to the signal data array.
 * @param kmerSeq Pointer to the k-mer sequence array.
 * @param forAPSEI Forward probabilities stored in an unordered map.
 * @param allowedKeys Array of allowed indices for processing.
 * @param T Length of the ONT raw signal + 1.
 * @param N Length of nucleotide sequence + 1.
 * @param K The number of k-mers.
 * @param model array containing tuples of model parameters.
 */
void logF(const double *sig, const int *kmerSeq, std::unordered_map<std::size_t, std::array<dproxy, NUMMAT>> &forAPSEI, const std::vector<std::size_t> &allowedKeys, const std::size_t K, const std::tuple<double, double> *model)
{
    std::array<dproxy, NUMMAT> forAPSEIRef;
    for (const std::size_t &tnk : allowedKeys)
    {
        // tnk = t*NK+n*K+k
        const std::size_t t = tnk / NK;
        const std::size_t n = (tnk % NK) / K;
        const std::size_t k = tnk % K;
        double a = -INFINITY;
        double p = -INFINITY;
        double s = -INFINITY;
        double e = -INFINITY;
        double i = -INFINITY;
        if (t == 0 && n == 0) [[unlikely]]
        {
            e = 0;
        }
        else if (t > 0 && n > 0) [[likely]]
        {
            // Precompute expensive values
            const double sc = score(sig[t - 1], kmerSeq[n - 1], k, model);
            const std::size_t baseIdx1 = tnk - NK - K - k;
            const std::size_t baseIdx2 = tnk - NK - k;
            const std::size_t baseIdx3 = tnk - NK - K;
            const std::size_t baseIdx4 = tnk - NK;
            const std::size_t baseIdx5 = tnk - K;

            // non-consecutive, differs by 5^4
            for (std::size_t preKmer = precessingKmer(k, 0, stepSize, alphabetSize); preKmer < K; preKmer += stepSize)
            {
                // (t-1)*NK+(n-1)*K+k'
                forAPSEIRef = forAPSEI[baseIdx1 + preKmer];
                a = logPlus(a, forAPSEIRef[3] + transitions_NTK.at("a1") + sc);
                a = logPlus(a, forAPSEIRef[4] + transitions_NTK.at("a2") + sc);

                // (t-1)*NK+n*K+k'
                forAPSEIRef = forAPSEI[baseIdx2 + preKmer];
                p = logPlus(p, forAPSEIRef[2] + transitions_NTK.at("p1") + sc);
                p = logPlus(p, forAPSEIRef[3] + transitions_NTK.at("p2") + sc);
                p = logPlus(p, forAPSEIRef[4] + transitions_NTK.at("p3") + sc);
            }

            // (t-1)*NK+(n-1)*K+k
            forAPSEIRef = forAPSEI[baseIdx3];
            s = logPlus(s, forAPSEIRef[1] + transitions_NTK.at("s1") + sc);
            s = logPlus(s, forAPSEIRef[3] + transitions_NTK.at("s2") + sc);
            s = logPlus(s, forAPSEIRef[4] + transitions_NTK.at("s3") + sc);

            // (t-1)*NK+n*K+k
            forAPSEIRef = forAPSEI[baseIdx4];
            e = logPlus(e, forAPSEIRef[0] + sc); // e1 always 1
            e = logPlus(e, forAPSEIRef[1] + transitions_NTK.at("e2") + sc);
            e = logPlus(e, forAPSEIRef[2] + transitions_NTK.at("e3") + sc);
            e = logPlus(e, forAPSEIRef[3] + transitions_NTK.at("e4") + sc);

            // t*NK+(n-1)*K+k
            forAPSEIRef = forAPSEI[baseIdx5];
            i = logPlus(i, forAPSEIRef[3] + transitions_NTK.at("i1") + sc);
            i = logPlus(i, forAPSEIRef[4] + transitions_NTK.at("i2") + sc);
        }
        forAPSEI[tnk] = {a, p, s, e, i};
    }
}

/**
 * Computes backward probabilities for a hidden Markov model (HMM) using signal
 * data, k-mer sequences, and transition probabilities.
 *
 * @param sig Pointer to the signal data array.
 * @param kmerSeq Pointer to the k-mer sequence array.
 * @param backAPSEI Backward probabilities stored in an unordered map.
 * @param allowedKeys Array of allowed indices for processing.
 * @param T Length of the ONT raw signal + 1.
 * @param N Length of nucleotide sequence + 1.
 * @param K The number of k-mers.
 * @param model Array containing tuples of model parameters.
 */
void logB(const double *sig, const int *kmerSeq, std::unordered_map<std::size_t, std::array<dproxy, NUMMAT>> &backAPSEI, std::vector<std::size_t> &allowedKeys, const std::size_t T, const std::size_t N, const std::size_t K, const std::tuple<double, double> *model)
{
    double pv, sc;
    for (auto tnk = allowedKeys.rbegin(); tnk != allowedKeys.rend(); ++tnk)
    {
        const std::size_t t = *tnk / NK;
        const std::size_t n = (*tnk % NK) / K;
        const std::size_t k = *tnk % K;
        double a = -INFINITY;
        double p = -INFINITY;
        double s = -INFINITY;
        double e = -INFINITY;
        double i = -INFINITY;
        if (t == T - 1 && n == N - 1) [[unlikely]]
        {
            e = 0;
        }
        if (t < T - 1) [[likely]]
        {
            // (t+1)*NK+nK+k
            const std::size_t tnkNK = *tnk + NK;
            // (t+1)*NK+(n+1)*K+k
            const std::size_t tnkNKK = tnkNK + K;
            // Cache results of successingKmer to avoid recomputation
            const std::size_t sucKmerBase = successingKmer(k, 0, stepSize, alphabetSize);
            const std::size_t sucKmerEnd = sucKmerBase + alphabetSize;

            if (n > 0) [[likely]]
            {
                // Precompute score for efficiency
                sc = score(sig[t], kmerSeq[n - 1], k, model);

                // (t+1)*NK+n*K+k
                pv = backAPSEI[tnkNK][3];
                a = logPlus(a, pv + sc); // e1 always 1
                p = logPlus(p, pv + transitions_NTK.at("e2") + sc);
                s = logPlus(s, pv + transitions_NTK.at("e3") + sc);
                e = logPlus(e, pv + transitions_NTK.at("e4") + sc);

                // kmer int representation is consecutive
                for (std::size_t sucKmer = sucKmerBase; sucKmer < sucKmerEnd; ++sucKmer)
                {
                    sc = score(sig[t], kmerSeq[n - 1], sucKmer, model);
                    // (t+1)*NK+n*K+sucKmer
                    pv = backAPSEI[tnkNK - k + sucKmer][1];
                    s = logPlus(s, pv + transitions_NTK.at("p1") + sc);
                    e = logPlus(e, pv + transitions_NTK.at("p2") + sc);
                    i = logPlus(i, pv + transitions_NTK.at("p3") + sc);
                }
            }

            if (n < N - 1) [[likely]]
            {
                sc = score(sig[t], kmerSeq[n], k, model);
                // (t+1)*(NK)+(n+1)*K+k
                pv = backAPSEI[tnkNKK][2];
                p = logPlus(p, pv + transitions_NTK.at("s1") + sc);
                e = logPlus(e, pv + transitions_NTK.at("s2") + sc);
                i = logPlus(i, pv + transitions_NTK.at("s3") + sc);

                // kmer int representation is consecutive
                for (std::size_t sucKmer = sucKmerBase; sucKmer < sucKmerEnd; ++sucKmer)
                {
                    sc = score(sig[t], kmerSeq[n], sucKmer, model);
                    // (t+1)*NK+(n+1)*K+sucKmer
                    pv = backAPSEI[tnkNKK - k + sucKmer][0];
                    e = logPlus(e, pv + transitions_NTK.at("a1") + sc);
                    i = logPlus(i, pv + transitions_NTK.at("a2") + sc);
                }
            }
        }

        if (t > 0 && n < N - 1) [[likely]]
        {
            sc = score(sig[t - 1], kmerSeq[n], k, model);
            // t*NK+(n+1)*K+k
            pv = backAPSEI[*tnk + K][4];
            e = logPlus(e, pv + transitions_NTK.at("i1") + sc);
            i = logPlus(i, pv + transitions_NTK.at("i2") + sc);
        }

        backAPSEI[*tnk] = {a, p, s, e, i};
    }
}

/**
 * @brief Computes the segmentation borders for a sequence.
 *
 * This function calculates the segmentation borders based on a list of allowed
 * k-mer indices, forward probabilities, and log probabilities. It updates the
 * segment strings with the highest probability k-mer sequence for the given
 * time and sequence index.
 *
 * @param segString List to store the resulting segment strings.
 * @param logAPSEI Unordered map of log probabilities for each state.
 * @param allowedKeys Array of allowed k-mer indices for processing.
 * @param T Length of the ONT raw signal + 1.
 * @param N Length of nucleotide sequence + 1.
 * @param K Number of k-mers.
 */
void getBorders(std::list<std::string> &segString, std::unordered_map<std::size_t, std::array<dproxy, NUMMAT>> &logAPSEI, const std::vector<std::size_t> &allowedKeys, const std::size_t T, const std::size_t N, const std::size_t K)
{
    std::unordered_map<std::size_t, std::array<dproxy, NUMMAT>> APSEI;
    std::size_t idx;
    for (const std::size_t &tnk : allowedKeys)
    {
        // tnk = t*NK+n*K+k
        const std::size_t t = tnk / NK;
        const std::size_t n = (tnk % NK) / K;
        const std::size_t k = tnk % K;
        double a = -INFINITY;
        double p = -INFINITY;
        double s = -INFINITY;
        double e = -INFINITY;
        double i = -INFINITY;
        if (t == 0 && n == 0) [[unlikely]]
        {
            e = 0;
        }
        if (t > 0 && n > 0) [[likely]]
        {
            // Precompute expensive values
            const std::size_t baseIdx1 = tnk - NK - K - k;                   // (t-1)*NK+(n-1)*K+k'
            const std::size_t baseIdx2 = tnk - NK - k;                       // (t-1)*NK+ n   *K+k'
            const std::size_t baseIdx3 = tnk - NK - K;                       // (t-1)*NK+(n-1)*K+k
            const std::size_t baseIdx4 = tnk - NK;                           // (t-1)*NK+ n   *K+k
            const std::size_t baseIdx5 = tnk - K;                            //  t   *NK+(n-1)*K+k
            const std::array<dproxy, NUMMAT> logAPSEIRef = logAPSEI.at(tnk); // .at(idx) : index must exist
            for (std::size_t preKmer = precessingKmer(k, 0, stepSize, alphabetSize); preKmer < K; preKmer += stepSize)
            {
                // (t-1)*NK+(n-1)*K+k'
                idx = baseIdx1 + preKmer;
                a = std::max(a, APSEI[idx][3] + logAPSEIRef[0]);
                a = std::max(a, APSEI[idx][4] + logAPSEIRef[0]);

                // (t-1)*NK+n*K+k'
                idx = baseIdx2 + preKmer;
                p = std::max(p, APSEI[idx][2] + logAPSEIRef[1]);
                p = std::max(p, APSEI[idx][3] + logAPSEIRef[1]);
                p = std::max(p, APSEI[idx][4] + logAPSEIRef[1]);
            }
            // (t-1)*NK+(n-1)*K+k
            s = std::max(s, APSEI[baseIdx3][1] + logAPSEIRef[2]);
            s = std::max(s, APSEI[baseIdx3][3] + logAPSEIRef[2]);
            s = std::max(s, APSEI[baseIdx3][4] + logAPSEIRef[2]);

            // (t-1)*NK+n*K+k
            e = std::max(e, APSEI[baseIdx4][0] + logAPSEIRef[3]);
            e = std::max(e, APSEI[baseIdx4][1] + logAPSEIRef[3]);
            e = std::max(e, APSEI[baseIdx4][2] + logAPSEIRef[3]);
            e = std::max(e, APSEI[baseIdx4][3] + logAPSEIRef[3]);

            // t*NK+(n-1)*K+k
            i = std::max(i, APSEI[baseIdx5][3] + logAPSEIRef[4]);
            i = std::max(i, APSEI[baseIdx5][4] + logAPSEIRef[4]);
        }
        APSEI[tnk] = {a, p, s, e, i};
    }

    // std::get highest k for last t and n?
    double mv = -INFINITY;
    std::size_t hk = 3125, lastDim = TNK - K;
    for (std::size_t k = 0; k < K; ++k)
    {
        if (APSEI[lastDim + k][3] >= mv)
        {
            mv = APSEI[lastDim + k][3];
            hk = k;
        }
    }

    // std::vector<double> segProb;
    // funcE(T - 1, N - 1, hk, APSEI, logAPSEI, segString, K, segProb);
    traceback(T - 1, N - 1, hk, APSEI, logAPSEI, segString, K);
}

/**
 * @brief Backtracing function for A state.
 *
 * This function backtraces the A state.
 * It checks which state was the previous state and calls the corresponding function.
 * If no match is found, it outputs an error message.
 *
 * @param t Current time index.
 * @param n Current segment index.
 * @param k Current kmer index.
 * @param APSEI The 2D array of scores.
 * @param logAPSEI The 2D array of log scores.
 * @param segString The list to store the segment strings.
 * @param K number of kemrs
 *
 * @details
 * This function backtraces the A state.
 * It checks which state was the previous state and calls the corresponding function.
 * If no match is found, it outputs an error message.
 */
void traceback(std::size_t t, std::size_t n, std::size_t k,
               std::unordered_map<std::size_t, std::array<dproxy, NUMMAT>> &APSEI,
               std::unordered_map<std::size_t, std::array<dproxy, NUMMAT>> &logAPSEI,
               std::list<std::string> &segString, const std::size_t K)
{
    std::vector<double> segProb;

    std::size_t state = 3;
    // 0 : A
    // 1 : P
    // 2 : S
    // 3 : E
    // 4 : I

    std::size_t currentIdx = t * NK + n * K + k;
    while (t) // t>0
    {
        // state E
        if (state == 3)
        {
            if (t == 1)
            {
                // Start value
                segString.push_front("M" + std::to_string(halfKmerSize) + "," + std::to_string(0) + "," + formattedMedian(segProb) + "," + itoa(k, alphabetSize, kmerSize, rna) + ";");
                segProb.clear();
                break; // stop traceback
            }
            // const std::size_t prevIdx = NK * (t - 1) + n * K + k;
            const std::size_t prevIdx = currentIdx - NK;
            const double sc = APSEI[currentIdx][3];
            const double logScore = logAPSEI[currentIdx][3];
            segProb.push_back(exp(logScore));
            // Check match with E state
            if (sc == APSEI[prevIdx][3] + logScore)
            {
                state = 3;
            }
            // Check match with A state
            else if (sc == APSEI[prevIdx][0] + logScore)
            {
                state = 0;
            }
            // Check match with S state
            else if (sc == APSEI[prevIdx][2] + logScore)
            {
                state = 2;
            }
            // Check match with P state
            else if (sc == APSEI[prevIdx][1] + logScore)
            {
                state = 1;
            }
            else
            {
                std::cerr << "Error in backtracing funcE!: t: " << t << ", n: " << n << ", k: " << k << "\n";
            }
            --t;
            currentIdx = prevIdx;
        }
        // state A
        else if (state == 0)
        {
            if (t == 1 && n == 1)
            {
                // Start value
                segString.push_front("M" + std::to_string(halfKmerSize) + "," + std::to_string(0) + "," + formattedMedian(segProb) + "," + itoa(k, alphabetSize, kmerSize, rna) + ";");
                segProb.clear();
                break; // stop traceback
            }
            // const std::size_t prevIdx = (t - 1) * NK + (n - 1) * K + 0;
            const std::size_t prevIdx = currentIdx - NK - K - k;
            const double sc = APSEI[currentIdx][0];
            const double logScore = logAPSEI[currentIdx][0];
            segProb.push_back(exp(logScore));
            for (std::size_t preKmer = precessingKmer(k, 0, stepSize, alphabetSize); preKmer < K; preKmer += stepSize)
            {
                // Check match with E state
                if (sc == APSEI[prevIdx + preKmer][3] + logScore)
                {
                    segString.push_front("M" + std::to_string(n - 1 + halfKmerSize) + "," + std::to_string(t - 1) + "," + formattedMedian(segProb) + "," + itoa(k, alphabetSize, kmerSize, rna) + ";");
                    segProb.clear();
                    state = 3;
                    --t;
                    --n;
                    k = preKmer;
                    currentIdx = prevIdx + preKmer;
                    break;
                }
                // Check match with I state
                else if (sc == APSEI[prevIdx + preKmer][4] + logScore)
                {
                    segString.push_front("M" + std::to_string(n - 1 + halfKmerSize) + "," + std::to_string(t - 1) + "," + formattedMedian(segProb) + "," + itoa(k, alphabetSize, kmerSize, rna) + ";");
                    segProb.clear();
                    state = 4;
                    --t;
                    --n;
                    k = preKmer;
                    currentIdx = prevIdx + preKmer;
                    break;
                }
            }
        }
        // state P
        else if (state == 1)
        {
            if (t == 1)
            {
                // Start value
                segString.push_front("P" + std::to_string(halfKmerSize) + "," + std::to_string(0) + "," + formattedMedian(segProb) + "," + itoa(k, alphabetSize, kmerSize, rna) + ";");
                segProb.clear();
                break;
            }
            // const std::size_t prevBaseIdx = NK * (t - 1) + n * K + 0;
            const std::size_t prevIdx = currentIdx - NK - k;
            const double sc = APSEI[currentIdx][1];
            const double logScore = logAPSEI[currentIdx][1];
            segProb.push_back(exp(logScore));
            for (std::size_t preKmer = precessingKmer(k, 0, stepSize, alphabetSize); preKmer < K; preKmer += stepSize)
            {
                // Check match with E state
                if (sc == APSEI[prevIdx + preKmer][3] + logScore)
                {
                    segString.push_front("P" + std::to_string(n - 1 + halfKmerSize) + "," + std::to_string(t - 1) + "," + formattedMedian(segProb) + "," + itoa(k, alphabetSize, kmerSize, rna) + ";");
                    segProb.clear();
                    state = 3;
                    --t;
                    k = preKmer;
                    currentIdx = prevIdx + preKmer;
                    break;
                }
                // Check match with S state
                if (sc == APSEI[prevIdx + preKmer][2] + logScore)
                {
                    segString.push_front("P" + std::to_string(n - 1 + halfKmerSize) + "," + std::to_string(t - 1) + "," + formattedMedian(segProb) + "," + itoa(k, alphabetSize, kmerSize, rna) + ";");
                    segProb.clear();
                    state = 2;
                    --t;
                    k = preKmer;
                    currentIdx = prevIdx + preKmer;
                    break;
                }
                // Check match with I state
                if (sc == APSEI[prevIdx + preKmer][4] + logScore)
                {
                    segString.push_front("P" + std::to_string(n - 1 + halfKmerSize) + "," + std::to_string(t - 1) + "," + formattedMedian(segProb) + "," + itoa(k, alphabetSize, kmerSize, rna) + ";");
                    segProb.clear();
                    state = 4;
                    --t;
                    k = preKmer;
                    currentIdx = prevIdx + preKmer;
                    break;
                }
            }
        }
        // state S
        else if (state == 2)
        {
            if (t == 1 && n == 1)
            {
                break; // stop traceback
            }
            // const std::size_t prevIdx = NK * (t - 1) + (n - 1) * K + k;
            const std::size_t prevIdx = currentIdx - NK - K;
            const double sc = APSEI[currentIdx][2];
            const double logScore = logAPSEI[currentIdx][2];
            segProb.push_back(exp(logScore));

            // Check match with E state
            if (sc == APSEI[prevIdx][3] + logScore)
            {
                state = 3;
            }
            // Check match with P state
            else if (sc == APSEI[prevIdx][1] + logScore)
            {
                state = 1;
            }
            // Check match with I state
            else if (sc == APSEI[prevIdx][4] + logScore)
            {
                state = 4;
            }
            --t;
            --n;
            currentIdx = prevIdx;
        }
        // state I
        else if (state == 4)
        {
            if (n == 1)
            {
                break; // stop traceback
            }
            // const std::size_t prevIdx = NK * t + (n - 1) * K + k;
            const std::size_t prevIdx = currentIdx - K;
            const double sc = APSEI[currentIdx][4];
            const double logScore = logAPSEI[currentIdx][4];
            segProb.push_back(exp(logScore));

            // Check match with I state
            if (sc == APSEI[prevIdx][4] + logScore)
            {
                state = 4;
            }
            // Check match with E state
            if (sc == APSEI[prevIdx][3] + logScore)
            {
                state = 3;
            }
            --n;
            currentIdx = prevIdx;
        }
    }
}

/**
 * Computes transition probabilities for a hidden Markov model (HMM) using
 * forward and backward probabilities, signal data, and k-mer sequences.
 *
 * @param sig Pointer to the signal data array.
 * @param kmerSeq Pointer to the k-mer sequence array.
 * @param forAPSEI Forward probabilities stored in an unordered map.
 * @param backAPSEI Backward probabilities stored in an unordered map.
 * @param allowedKeys Array of allowed indices for processing.
 * @param T Length of the ONT raw signal + 1.
 * @param N Length of nucleotide sequence + 1.
 * @param K The number of k-mers.
 * @param model Array containing tuples of model parameters.
 *
 * @return A tuple containing the transition probabilities for different states
 *         in the HMM.
 */
std::tuple<double, double, double, double, double, double, double, double, double, double, double, double, double, double> trainTransition(const double *sig, const int *kmerSeq, std::unordered_map<std::size_t, std::array<dproxy, NUMMAT>> &forAPSEI, std::unordered_map<std::size_t, std::array<dproxy, NUMMAT>> &backAPSEI, std::vector<std::size_t> &allowedKeys, const std::size_t T, const std::size_t N, const std::size_t K, const std::tuple<double, double> *model)
{
    // Transition parameters
    double newa1 = -INFINITY, newa2 = -INFINITY;
    double newp1 = -INFINITY, newp2 = -INFINITY, newp3 = -INFINITY;
    double news1 = -INFINITY, news2 = -INFINITY, news3 = -INFINITY;
    double newe1 = 0, newe2 = -INFINITY, newe3 = -INFINITY, newe4 = -INFINITY;
    double newi1 = -INFINITY, newi2 = -INFINITY;
    double sc, pv;

    for (const std::size_t &tnk : allowedKeys)
    {
        // tnk = t*NK+n*K+k
        const std::size_t t = tnk / NK;
        const std::size_t n = (tnk % NK) / K;
        const std::size_t k = tnk % K;
        // Cache results to avoid recomputation
        auto &forAPSEI_tnk = forAPSEI[tnk];

        if (t < T - 1) [[likely]]
        {
            const std::size_t sucKmerBase = successingKmer(k, 0, stepSize, alphabetSize);
            const std::size_t sucKmerEnd = sucKmerBase + alphabetSize;

            if (n > 0) [[likely]]
            {
                sc = score(sig[t], kmerSeq[n - 1], k, model);
                // (t+1)*NK+n*K+k
                pv = backAPSEI[tnk + NK][3];
                // newe1 = logPlus(newe1, forAPSEI_tnk[0] + e1 + sc + pv);
                newe2 = logPlus(newe2, forAPSEI_tnk[1] + transitions_NTK.at("e2") + sc + pv);
                newe3 = logPlus(newe3, forAPSEI_tnk[2] + transitions_NTK.at("e3") + sc + pv);
                newe4 = logPlus(newe4, forAPSEI_tnk[3] + transitions_NTK.at("e4") + sc + pv);

                for (std::size_t sucKmer = sucKmerBase; sucKmer < sucKmerEnd; ++sucKmer)
                {
                    sc = score(sig[t], kmerSeq[n - 1], sucKmer, model);
                    // (t+1)*NK+n*K+sucKmer
                    pv = backAPSEI[tnk + NK - k + sucKmer][1];
                    newp1 = logPlus(newp1, forAPSEI_tnk[2] + transitions_NTK.at("p1") + sc + pv);
                    newp2 = logPlus(newp2, forAPSEI_tnk[3] + transitions_NTK.at("p2") + sc + pv);
                    newp3 = logPlus(newp3, forAPSEI_tnk[4] + transitions_NTK.at("p3") + sc + pv);
                }
            }

            if (n < N - 1) [[likely]]
            {
                sc = score(sig[t], kmerSeq[n], k, model);
                // (t+1)*(NK)+(n+1)*K+k
                pv = backAPSEI[tnk + NK + K][2];
                news1 = logPlus(news1, forAPSEI_tnk[1] + transitions_NTK.at("s1") + sc + pv);
                news2 = logPlus(news2, forAPSEI_tnk[3] + transitions_NTK.at("s2") + sc + pv);
                news3 = logPlus(news3, forAPSEI_tnk[4] + transitions_NTK.at("s3") + sc + pv);

                for (std::size_t sucKmer = sucKmerBase; sucKmer < sucKmerEnd; ++sucKmer)
                {
                    sc = score(sig[t], kmerSeq[n], sucKmer, model);
                    // (t+1)*NK+(n+1)*K+sucKmer
                    pv = backAPSEI[tnk + NK + K - k + sucKmer][0];
                    newa1 = logPlus(newa1, forAPSEI_tnk[3] + transitions_NTK.at("a1") + sc + pv);
                    newa2 = logPlus(newa2, forAPSEI_tnk[4] + transitions_NTK.at("a2") + sc + pv);
                }
            }
        }

        if (t > 0 && n < N - 1) [[likely]]
        {
            sc = score(sig[t - 1], kmerSeq[n], k, model);
            // t*NK+(n+1)*K+k
            pv = backAPSEI[tnk + K][4];
            newi1 = logPlus(newi1, forAPSEI_tnk[3] + transitions_NTK.at("i1") + sc + pv);
            newi2 = logPlus(newi2, forAPSEI_tnk[4] + transitions_NTK.at("i2") + sc + pv);
        }
    }

    // Final normalization and averaging of transition parameters
    // newe1=exp(newe1-newe1);
    const double Ae = logPlus(logPlus(logPlus(newa1, news2), logPlus(newe4, newi1)), newp2);
    if (!std::isinf(Ae))
    {
        newa1 = newa1 - Ae;
        news2 = news2 - Ae;
        newe4 = newe4 - Ae;
        newi1 = newi1 - Ae;
        newp2 = newp2 - Ae;
    }
    const double As = logPlus(newe3, newp1);
    if (!std::isinf(As))
    {
        newe3 = newe3 - As;
        newp1 = newp1 - As;
    }
    const double Ap = logPlus(newe2, news1);
    if (!std::isinf(Ap))
    {
        newe2 = newe2 - Ap;
        news1 = news1 - Ap;
    }
    const double Ai = logPlus(logPlus(newa2, newi2), logPlus(newp3, news3));
    if (!std::isinf(Ai))
    {
        newa2 = newa2 - Ai;
        newi2 = newi2 - Ai;
        newp3 = newp3 - Ai;
        news3 = news3 - Ai;
    }

    return std::make_tuple(
        exp(newa1),
        exp(newa2),
        exp(newp1),
        exp(newp2),
        exp(newp3),
        exp(news1),
        exp(news2),
        exp(news3),
        exp(newe1),
        exp(newe2),
        exp(newe3),
        exp(newe4),
        exp(newi1),
        exp(newi2));
}

/**
 * @brief Trains emission parameters using the Baum-Welch algorithm.
 *
 * @param[in] sig Pointer to the ONT raw signal array with pA values.
 * @param[in] logAPSEI Log probability matrix for match states.
 * @param[in] allowedKeys Array of allowed transitions_NTK.
 * @param[in] T Length of the ONT raw signal + 1.
 * @param[in] N Length of nucleotide sequence + 1.
 * @param[in] K Number of kmers in the sequence.
 *
 * Outputs the trained emission parameters to the console.
 */
std::tuple<double *, double *> trainEmission(const double *sig, std::unordered_map<std::size_t, std::array<dproxy, NUMMAT>> &logAPSEI, std::vector<std::size_t> &allowedKeys, const std::size_t K)
{
    // Emission
    // https://courses.grainger.illinois.edu/ece417/fa2021/lectures/lec15.pdf
    // https://f.hubspotusercontent40.net/hubfs/8111846/Unicon_October2020/pdf/bilmes-em-algorithm.pdf
    // 0-Initialize memory
    double *means = new double[K]();
    double *stdevs = new double[K]();
    double *normFactorT = new double[K]();
    double w;

    // First pass: Compute means and normalization factors
    for (const std::size_t &tnk : allowedKeys)
    {
        if (tnk < NK) [[unlikely]]
        {
            continue;
        }
        // tnk = t*NK+n*K+k
        const std::size_t t = tnk / NK;
        const std::size_t k = tnk % K;
        // Cache logAPSEI access
        const auto &logValues = logAPSEI.at(tnk);
        w = logPlus(logPlus(logPlus(logValues[0], logValues[1]), logPlus(logValues[2], logValues[3])), logValues[4]);
        w = exp(w); // Convert log probability to normal probability
        means[k] += w * sig[t - 1];
        normFactorT[k] += w;
    }

    // Normalize the means
    for (std::size_t k = 0; k < K; ++k)
    {
        if (normFactorT[k] != 0.0)
        {
            means[k] = means[k] / normFactorT[k];
        }
    }

    // Emission (stdev of kmers)
    // assuming a flat prior and integrating over N, every cell (in T x K) has on avg a prob. of 1/K, integrating over T results in T/K
    // used kmers should exceed the threshold easily
    // const double TRAINTHRESHOLD = T/K;
    const double TRAINTHRESHOLD = 1e-7; // set by eye
    for (const std::size_t &tnk : allowedKeys)
    {
        // tnk = t*NK+n*K+k
        const std::size_t k = tnk % K;
        // Skip kmers with low weight or if t == 0
        if (normFactorT[k] < TRAINTHRESHOLD || tnk < NK) [[unlikely]]
        {
            continue;
        }
        const std::size_t t = tnk / NK;
        // Cache logAPSEI access
        const auto &logValues = logAPSEI.at(tnk);
        w = logPlus(logPlus(logPlus(logValues[0], logValues[1]), logPlus(logValues[2], logValues[3])), logValues[4]);
        // Update standard deviation
        const double diff = sig[t - 1] - means[k];
        stdevs[k] += exp(w) * diff * diff;
    }

    // Normalize the standard deviations
    for (std::size_t k = 0; k < K; ++k)
    {
        if (normFactorT[k] != 0.0)
        {
            stdevs[k] = sqrt(stdevs[k] / normFactorT[k]);
        }
    }
    delete[] normFactorT;
    return std::tuple<double *, double *>({means, stdevs});
}

/**
 * @brief Trains transition and emission parameters using the Baum-Welch algorithm.
 *
 * @param sig Pointer to the ONT raw signal array with pA values.
 * @param kmerSeq Pointer to the nucleotide sequence represented by the ONT signal.
 * @param forAPSEI Forward matrix for match states.
 * @param backAPSEI Backward matrix for match states.
 * @param logAPSEI Log probability matrix for match states.
 * @param allowedKeys Array of allowed transitions_NTK.
 * @param T Length of the ONT raw signal + 1.
 * @param N Length of nucleotide sequence + 1.
 * @param K Number of kmers in the sequence.
 * @param model Array containing kmers as keys and (mean, stdev) tuples as values.
 *
 * Outputs the trained parameters to the console.
 */
void trainParams(const double *sig, const int *kmerSeq, std::unordered_map<std::size_t, std::array<dproxy, NUMMAT>> &forAPSEI, std::unordered_map<std::size_t, std::array<dproxy, NUMMAT>> &backAPSEI, std::unordered_map<std::size_t, std::array<dproxy, NUMMAT>> &logAPSEI, std::vector<std::size_t> &allowedKeys, const std::size_t T, const std::size_t N, const std::size_t K, std::tuple<double, double> *model)
{

    // newa1, newa2, newp1, newp2, newp3, news1, news2, news3, newe1, newe2, newe3, newe4, newi1, newi2
    auto [a1, a2, p1, p2, p3, s1, s2, s3, e1, e2, e3, e4, i1, i2] = trainTransition(sig, kmerSeq, forAPSEI, backAPSEI, allowedKeys, T, N, K, model);
    std::cout << "a1:" << a1 << ";a2:" << a2 << ";p1:" << p1 << ";p2:" << p2 << ";p3:" << p3 << ";s1:" << s1 << ";s2:" << s2 << ";s3:" << s3 << ";e1:" << e1 << ";e2:" << e2 << ";e3:" << e3 << ";e4:" << e4 << ";i1:" << i1 << ";i2:" << i2 << "\n";

    auto [newMeans, newStdevs] = trainEmission(sig, logAPSEI, allowedKeys, K);
    for (std::size_t k = 0; k < K; ++k)
    {
        // if (newStdevs[k] != 0.0 && !std::isnan(newStdevs[k]))
        if (newStdevs[k] != 0.0)
        {
            std::cout << itoa(k, alphabetSize, kmerSize, rna) << ":" << newMeans[k] << "," << newStdevs[k] << ";";
        }
    }
    std::cout << "\n";
    delete[] newMeans;
    delete[] newStdevs;
}