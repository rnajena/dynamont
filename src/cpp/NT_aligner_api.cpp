#include "NT_aligner_api.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace dynamont
{

NTAligner::NTAligner(
    const std::string& modelFile,
    PoreType pore,
    int threads,
    std::size_t band):
    Aligner(modelFile, pore, threads, band)
{
    initializeTransitions();
}

void NTAligner::initializeTransitions()
{
    /*
        Values are stored in log-space.

        Defaults correspond to the original NT model:
        m1 = probability of match transition
        e1 = first extension
        e2 = further extension
    */

    double m1;
    double e1 = 1.0;
    double e2;

    switch (pore_)
    {
        case PoreType::RNA002:
        {
            rna_ = true;
            kmerSize_ = 5;
            m1 = 0.019889650396799997;
            e2 = 0.9801103496029998;
            break;
        }
        case PoreType::RNA004:
        {
            rna_ = true;
            kmerSize_ = 9;
            m1 = 0.031111753637096777;
            e2 = 0.9688882463622581;
            break;
        }
        case PoreType::DNA_R9:
        {
            rna_ = false;
            kmerSize_ = 5;
            m1 = 1.0;
            e2 = 1.0;
            break;
        }
        case PoreType::DNA_R10_260:
        {
            rna_ = false;
            kmerSize_ = 9;
            // currently same parameters as original implementation
            m1 = 0.031111753637096777;
            e2 = 0.9688882463622581;
            break;
        }
        case PoreType::DNA_R10_400:
        {
            rna_ = false;
            kmerSize_ = 9;
            // currently same parameters as original implementation
            m1 = 0.031111753637096777;
            e2 = 0.9688882463622581;
            break;
        }
        default:
            throw std::runtime_error("Unknown pore type");
    }

    transitions_.m1 = std::log(m1);
    transitions_.e1 = std::log(e1);
    transitions_.e2 = std::log(e2);
}


std::vector<Bound> NTAligner::computeBounds(
    std::size_t T,
    std::size_t N,
    std::size_t bandwidth) const
{
    std::vector<Bound> bounds(T);
    const double RATIO = static_cast<double>(N) / static_cast<double>(T);

    for (std::size_t t = 0; t < T; ++t)
    {
        std::size_t midpoint = static_cast<std::size_t>(t * RATIO);
        long start = static_cast<long>(midpoint) - static_cast<long>(bandwidth);
        std::size_t nStart = midpoint >= bandwidth ? midpoint - bandwidth : 0;
        std::size_t nEnd = midpoint + bandwidth + 1 <= N ? midpoint + bandwidth + 1 : N;
        bounds[t] = Bound{start, nStart, nEnd};
    }

    return bounds;
}

void NTAligner::forward(
    const double* signal,
    const int* kmers,
    double* M,
    double* E,
    std::size_t T,
    std::size_t B,
    std::size_t bandwidth,
    const std::vector<Bound>& bounds) const
{
    E[bandwidth + 1] = 0.0;

    for (std::size_t t = 1; t < T; ++t)
    {
        auto [bandStart, nStart, nEnd] = bounds[t];

        if (!nStart) [[unlikely]]
            nStart = 1;

        const std::size_t tB = t * B;
        const long offset = static_cast<long>(tB) - bandStart + 1;
        long mShift = -static_cast<long>(B) - 1;
        long eShift = -static_cast<long>(B);

        if (bounds[t - 1].start != bandStart)
        {
            ++mShift;
            ++eShift;
        }

        for (std::size_t n = nStart; n < nEnd; ++n)
        {
            const std::size_t idx = n + offset;
            const double score = scoreKmer(signal[t - 1], kmers[n - 1]);

            M[idx] = E[idx + mShift] + score + transitions_.m1;
            E[idx] = logPlus(
                M[idx + eShift] + score + transitions_.e1,
                E[idx + eShift] + score + transitions_.e2);
        }
    }
}

// ============================================================
// Backward DP
// ============================================================

void NTAligner::backward(
    const double* signal,
    const int* kmers,
    double* M,
    double* E,
    std::size_t T,
    std::size_t N,
    std::size_t B,
    std::size_t bandwidth,
    const std::vector<Bound>& bounds) const
{
    const double NEG_INF = -std::numeric_limits<double>::infinity();
    E[(T - 1) * B + bandwidth + 1] = 0.0;

    for (std::size_t t = T - 1; t-- > 0;)
    {
        auto [bandStart, nStart, nEnd] = bounds[t];
        const std::size_t tB = t * B;
        const long offset = static_cast<long>(tB) - bandStart + 1;
        long mShift = static_cast<long>(B) + 1;
        long eShift = static_cast<long>(B);

        if (bounds[t + 1].start != bandStart)
        {
            --mShift;
            --eShift;
        }

        for (std::size_t n = nStart; n < nEnd; ++n)
        {
            const std::size_t idx = n + offset;
            double extension = NEG_INF;

            if (n + 1 < N) [[likely]]
            {
                extension = M[idx + mShift] + scoreKmer(signal[t], kmers[n]) + transitions_.m1;
            }

            if (n > 0) [[likely]]
            {
                double score = scoreKmer(signal[t], kmers[n - 1]);
                M[idx] = E[idx + eShift] + score;
                extension = logPlus(extension, E[idx + eShift] + score + transitions_.e2);
            }

            E[idx] = extension;
        }
    }
}

// ============================================================
// Log posterior matrix
// ============================================================

void NTAligner::calculatePosterior(
    double* output,
    const double* forward,
    const double* backward,
    double Z,
    std::size_t size) const
{
    for (std::size_t i = 0; i < size; ++i)
    {
        output[i] = forward[i] + backward[i] - Z;
    }
}

// ============================================================
// Public alignment API
// ============================================================

Result NTAligner::align(
    const double* signal,
    std::size_t signalLength,
    const std::string& sequence,
    bool calcProbabilities)
{
    validateInput(signalLength, sequence.size());
    Result result;
    std::vector<int> kmers = sequenceToKmers(sequence);

    // Matrix dimensions (T x N) = (signal length + 1) x (kmer count + 1)
    const std::size_t T = signalLength + 1;
    const std::size_t N = kmers.size() + 1;
    const std::size_t bandwidth = std::min(bandwidth_, N / 2);
    const std::size_t B = 2 * bandwidth + 3;
    const std::size_t matrixSize = T * B;

    auto bounds = computeBounds(T, N, bandwidth);

    std::vector<double> forwardE(matrixSize);
    std::vector<double> forwardM(matrixSize);
    std::vector<double> backwardM(matrixSize);
    std::vector<double> backwardE(matrixSize);

    const double NEG_INF = -std::numeric_limits<double>::infinity();

    for(std::size_t i = 0; i < matrixSize; ++i)
    {
        forwardE[i] = NEG_INF;
        forwardM[i] = NEG_INF;
        backwardM[i] = NEG_INF;
        backwardE[i] = NEG_INF;
    }
    
    forward(
        signal,
        kmers.data(),
        forwardM.data(),
        forwardE.data(),
        T,
        B,
        bandwidth,
        bounds);

    backward(
        signal,
        kmers.data(),
        backwardM.data(),
        backwardE.data(),
        T,
        N,
        B,
        bandwidth,
        bounds);

    double Zf = forwardE[T * B - bandwidth - 2];
    double Zb = backwardE[bandwidth + 1];

    constexpr double epsilon = 1e-8;
    if (std::isinf(Zf) || std::isinf(Zb) ||
        std::abs(Zf - Zb) / static_cast<double>(matrixSize) > epsilon)
        throw std::runtime_error("Alignment failed: alignment scores do not match");

    result.Z = Zb;
    if (calcProbabilities == false)
        return result;

    std::vector<double> logProbM(matrixSize);
    std::vector<double> logProbE(matrixSize);
    calculatePosterior(logProbM.data(), forwardM.data(), backwardM.data(), Zb, matrixSize);
    calculatePosterior(logProbE.data(), forwardE.data(), backwardE.data(), Zb, matrixSize);

    result.segments = calculateSegments(
        logProbM.data(),
        logProbE.data(),
        T,
        N,
        B,
        bandwidth,
        bounds);

    return result;
}

std::vector<Segment> NTAligner::calculateSegments(
    const double* LPM,
    const double* LPE,
    std::size_t T,
    std::size_t N,
    std::size_t B,
    std::size_t bandwidth,
    const std::vector<Bound>& bounds
) const
{
    std::vector<Segment> result;
    result.reserve(N - 1);

    // Reconstruct traceback matrices.
    // These are not posterior matrices.
    // They contain the Viterbi path scores.
    
    const double NEG_INF = -std::numeric_limits<double>::infinity();
    const std::size_t size = T * B;
    std::vector<double> M(size, NEG_INF);
    std::vector<double> E(size, NEG_INF);

    E[bandwidth + 1] = 0.0;

    for (std::size_t t = 1; t < T; ++t)
    {
        auto [bandStart, nStart, nEnd] = bounds[t];
        
        // forward specific: if start == 0 move it to 1 to prevent 'if's in n-loop
        if (!nStart) [[unlikely]]
            nStart = 1;
 
        const std::size_t tB = t * B;
        const long offset = static_cast<long>(tB) - bandStart + 1; // TN to TB conversion offset
        long mShift = -static_cast<long>(B) - 1;
        long eShift = -static_cast<long>(B);

        if (bounds[t - 1].start != bandStart)
        {
            ++mShift;
            ++eShift;
        }

        for (std::size_t n = nStart; n < nEnd; ++n)
        {
            const std::size_t idx = n + offset;
            M[idx] = E[idx + mShift] + LPM[idx];
            E[idx] = std::max(M[idx + eShift], E[idx + eShift]) + LPE[idx];
        }
    }

    decodeMAP(
        M.data(),
        E.data(),
        LPM,
        LPE,
        T,
        N,
        B,
        bandwidth,
        bounds,
        result);
    return result;
}

// ============================================================
// MAP decoding
// ============================================================

void NTAligner::decodeMAP(
    const double* M,
    const double* E,
    const double* LPM,
    const double* LPE,
    std::size_t T,
    std::size_t N,
    std::size_t B,
    std::size_t bandwidth,
    const std::vector<Bound>& bounds,
    std::vector<Segment>& segments) const
{
    std::vector<double> segProb;
    bool isFuncM = false; // start in E

    std::size_t t = T - 1;
    std::size_t n = N - 1;
    // long bandIndex = static_cast<long>(N - 1) - bounds[T - 1].start;
    // std::size_t tBb = (T - 1) * B + static_cast<std::size_t>(bandIndex);
    std::size_t tBb = (T - 1) * B + bandwidth + 1;

    while (t && n)
    {
        long mShift = -static_cast<long>(B) - 1;
        long eShift = -static_cast<long>(B);

        if (bounds[t].start != bounds[t - 1].start)
        {
            ++mShift;
            ++eShift;
        }

        // Match state
        if (isFuncM)
        {
            segProb.push_back(std::exp(LPM[tBb]));

            // n - 1 + kmerSize / 2
            std::size_t position = n - 1 + kmerSize_ / 2;
            double probability = formattedMedian(segProb);

            segments.push_back(
                Segment{
                    'M',
                    position,
                    t - 1,
                    probability
                });


            segProb.clear();

            // transition M(t,n) -> E(t-1,n-1)
            --t;
            --n;
            tBb += mShift;
            isFuncM = false;
        }

        // Extension state
        else
        {
            double logScore = LPE[tBb];
            segProb.push_back(std::exp(logScore));
            // Decide whether previous state was M or E
            isFuncM = (E[tBb] == M[tBb + eShift] + logScore);
            --t;
            tBb += eShift;
        }
    }

    // because traceback runs backwards and vector only supports push_back, we need to reverse the segments vector to get the correct order
    std::reverse(segments.begin(), segments.end());
}

// ============================================================
// Training parameter estimation (one EM iteration)
// ============================================================

TrainingResult NTAligner::runTraining(
    const double* signal,
    const int* kmers,
    const double* forM,
    const double* forE,
    const double* backM,
    const double* backE,
    std::size_t T,
    std::size_t N,
    std::size_t B,
    const std::vector<Bound>& bounds,
    double Z)
{
    TrainingResult result;
    result.Z = Z;

    /*
        Expected emission statistics.

        For each kmer:
            weight
            weighted mean
            weighted variance

        Both Match and Extension states emit the observation,
        so their posterior probabilities are accumulated.
    */

    std::vector<double> weight(numKmers_, 0.0);
    std::vector<double> sum(numKmers_, 0.0);
    std::vector<double> sumSq(numKmers_, 0.0);

    for (std::size_t t = 1; t < T; ++t)
    {
        auto [bandStart, nStart, nEnd] = bounds[t];

        if (!nStart)
            nStart = 1;

        for (std::size_t n = nStart; n < nEnd; ++n)
        {
            const long offset = static_cast<long>(t * B) - bandStart + 1;
            const std::size_t idx = n + offset;
            const double posterior =
                std::exp(forM[idx] + backM[idx] - Z) +
                std::exp(forE[idx] + backE[idx] - Z);
            const int kmer = kmers[n - 1];
            const double observation = signal[t - 1];
            weight[kmer] += posterior;
            sum[kmer] += posterior * observation;
            sumSq[kmer] += posterior * observation * observation;
        }
    }

    // Update emission model
    result.emissionModel.resize(numKmers_);

    for (std::size_t k = 0; k < numKmers_; ++k)
    {
        if (weight[k] > 0.0)
        {
            double mean = sum[k] / weight[k];
            double variance = sumSq[k] / weight[k] - mean * mean;

            if (variance < 1e-12)
                variance = 1e-12;

            result.emissionModel[k] = EmissionParams{mean, std::sqrt(variance)};
        }
        else // no observations for this kmer
        {
            result.emissionModel[k] = model_[k];
        }
    }

    /*
        Transition update.

        At this stage we keep the existing
        transitions. A full Baum-Welch update
        requires explicit transition posterior
        matrices.

        The API is prepared for replacing this.
    */

    result.transitions =
        trainTransition(
            signal,
            kmers,
            forE,
            backM,
            backE,
            T,
            N,
            B,
            bounds);

    return result;
}

// ============================================================
// Public training API
// ============================================================

TrainingResult NTAligner::train(
    const double* signal,
    std::size_t signalLength,
    const std::string& sequence)
{
    validateInput(signalLength, sequence.size());

    std::vector<int> kmers = sequenceToKmers(sequence);

    const std::size_t T = signalLength + 1;
    const std::size_t N = kmers.size() + 1;
    const std::size_t bandwidth = std::min(bandwidth_, N / 2);
    const std::size_t B = 2 * bandwidth + 3;

    auto bounds = computeBounds(T, N, bandwidth);

    std::size_t matrixSize = T * B;
    std::vector<double> forwardM(matrixSize);
    std::vector<double> forwardE(matrixSize);
    std::vector<double> backwardM(matrixSize);
    std::vector<double> backwardE(matrixSize);

    const double NEG_INF = -std::numeric_limits<double>::infinity();
    for(std::size_t i = 0; i < matrixSize; ++i)
    {
        forwardE[i] = NEG_INF;
        forwardM[i] = NEG_INF;
        backwardM[i] = NEG_INF;
        backwardE[i] = NEG_INF;
    }

    forward(
        signal,
        kmers.data(),
        forwardM.data(),
        forwardE.data(),
        T,
        B,
        bandwidth,
        bounds);

    backward(
        signal,
        kmers.data(),
        backwardM.data(),
        backwardE.data(),
        T,
        N,
        B,
        bandwidth,
        bounds);

    const double Zf = forwardE[T * B - bandwidth - 2];
    const double Zb = backwardE[bandwidth + 1];

    constexpr double epsilon = 1e-8;
    if (std::isinf(Zf) || std::isinf(Zb) ||
        std::abs(Zf - Zb) / static_cast<double>(matrixSize) > epsilon)
        throw std::runtime_error("Training failed: alignment scores do not match");

    return runTraining(
        signal,
        kmers.data(),
        forwardM.data(),
        forwardE.data(),
        backwardM.data(),
        backwardE.data(),
        T,
        N,
        B,
        bounds,
        Zb);
}

TransitionParams NTAligner::trainTransition(
    const double* signal,
    const int* kmers,
    const double* forE,
    const double* backM,
    const double* backE,
    std::size_t T,
    std::size_t N,
    std::size_t B,
    const std::vector<Bound>& bounds) const
{
    double newM1 = -std::numeric_limits<double>::infinity();
    double newE2 = -std::numeric_limits<double>::infinity();
    const double m1 = transitions_.m1;
    const double e2 = transitions_.e2;
    long oldBandStart = bounds[T - 1].start;

    std::size_t tB = (T - 1) * B;

    for (std::size_t t = T - 1; t-- > 0;)
    {
        tB -= B;
        auto [bandStart, nStart, nEnd] = bounds[t];
        const long offset = static_cast<long>(tB) - bandStart + 1;

        long mShift = static_cast<long>(B) + 1;
        long eShift = static_cast<long>(B);

        if (oldBandStart != bandStart)
        {
            --mShift;
            --eShift;
            oldBandStart = bandStart;
        }

        for (std::size_t n = nStart; n < nEnd; ++n)
        {
            const std::size_t idx = n + offset;

            // Match transition: E(t,n) -> M(t+1,n+1)
            if (n + 1 < N)
            {
                newM1 = logPlus(newM1, forE[idx] + m1 + scoreKmer(signal[t], kmers[n]) + backM[idx + mShift]);
            }

            // Extension transition: E(t,n) -> E(t+1,n)
            if (n > 0)
            {
                const double score = scoreKmer(signal[t], kmers[n - 1]);
                newE2 = logPlus(newE2, forE[idx] + e2 + score + backE[idx + eShift]);
            }
        }
    }

    /*
        Normalize transition probabilities.

        Original code normalizes only m1/e2.
        e1 remains unused because current model
        always uses deterministic first extension.
    */

    const double norm = logPlus(newM1, newE2);
    TransitionParams result{};

    if (!std::isinf(norm))
    {
        result.m1 = std::exp(newM1 - norm);
        result.e2 = std::exp(newE2 - norm);
    }
    else
    {
        result.m1 = 0.0;
        result.e2 = 0.0;
    }

    /*
        Keep compatibility with original model.
        e1 was commented out in trainTransition().
    */

    result.e1 = std::exp(transitions_.e1);

    return result;
}

}