#pragma once

#include <vector>

#include "aligner.hpp"

namespace dynamont
{

struct Bound
{
    long start;
    std::size_t nStart;
    std::size_t nEnd;
};

class NTAligner : public Aligner
{
public:
    NTAligner(
        const std::string& modelFile,
        PoreType pore,
        int threads = 1,
        std::size_t band = 400);

    Result align(
        const double* signal,
        std::size_t signalLength,
        const std::string& sequence,
        bool calcProbabilities = false) override;

    TrainingResult train(
        const double* signal,
        std::size_t signalLength,
        const std::string& sequence) override;

private:
    TransitionParams transitions_;

    void initializeTransitions();

    std::vector<Bound> computeBounds(
        std::size_t T,
        std::size_t N,
        std::size_t bandwidth) const;

    void forward(
        const double* signal,
        const int* kmers,
        double* M,
        double* E,
        std::size_t T,
        std::size_t B,
        std::size_t bandwidth,
        const std::vector<Bound>& bounds) const;

    void backward(
        const double* signal,
        const int* kmers,
        double* M,
        double* E,
        std::size_t T,
        std::size_t N,
        std::size_t B,
        std::size_t bandwidth,
        const std::vector<Bound>& bounds) const;

    std::vector<Segment> calculateSegments(
        const double* LPM,
        const double* LPE,
        std::size_t T,
        std::size_t N,
        std::size_t B,
        std::size_t bandwidth,
        const std::vector<Bound>& bounds) const;

    void decodeMAP(
        const double* M,
        const double* E,
        const double* LPM,
        const double* LPE,
        std::size_t T,
        std::size_t N,
        std::size_t B,
        std::size_t bandwidth,
        const std::vector<Bound>& bounds,
        std::vector<Segment>& segments) const;

    void calculatePosterior(
        double* output,
        const double* forward,
        const double* backward,
        double Z,
        std::size_t size) const;

    TrainingResult runTraining(
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
        double Z);

    TransitionParams trainTransition(
        const double* signal,
        const int* kmers,
        const double* forE,
        const double* backM,
        const double* backE,
        std::size_t T,
        std::size_t N,
        std::size_t B,
        const std::vector<Bound>& bounds) const;
};

} // namespace dynamont