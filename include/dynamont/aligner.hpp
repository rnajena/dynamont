#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace dynamont
{

struct TransitionParams
{
	double m1;
	double e1;
	double e2;
};

struct EmissionParams
{
	double mean;
	double stdev;
};

enum class PoreType
{
	RNA002,
	RNA004,
	DNA_R9,
	DNA_R10_260,
	DNA_R10_400
};

struct Segment
{
	char state;
	std::size_t sequencePosition;
	std::size_t signalPosition;
	double probability;
	std::string polish = "";
};

struct Result
{
	double Z = 0.0;
	std::vector<Segment> segments;
};

struct TrainingResult
{
	double Z = 0.0;
	TransitionParams transitions;
	std::vector<EmissionParams> emissionModel;
};

double log_normal_pdf(double x, double mean, double stdev);

class Aligner
{
public:
	Aligner(
		const std::string& modelFile,
		PoreType pore,
		int threads = 1,
		std::size_t band = 400);

	virtual ~Aligner() = default;

	Aligner(const Aligner&) = delete;
	Aligner& operator=(const Aligner&) = delete;

	virtual Result align(
		const double* signal,
		std::size_t signalLength,
		const std::string& sequence,
		bool calcProbabilities = false) = 0;

	virtual TrainingResult train(
		const double* signal,
		std::size_t signalLength,
		const std::string& sequence);

protected:
	void validateInput(
		std::size_t signalLength,
		std::size_t readLength) const;

	std::vector<int> sequenceToKmers(const std::string& sequence) const;
	int kmerToInt(const std::string& sequence) const;
	std::string intToKmer(std::size_t value) const;

	double scoreKmer(double signalValue, std::size_t kmer) const;

	std::size_t successorKmer(std::size_t currentKmer, int nextNt) const;
	std::size_t predecessorKmer(std::size_t currentKmer, int priorNt) const;
	double formattedMedian(std::vector<double>& values) const;

	static double logPlus(double x, double y);

	PoreType pore_;
	bool rna_;
	int threads_;
	std::size_t bandwidth_;

	std::array<int8_t, 256> baseLookup_;
	std::vector<EmissionParams> model_;
	int alphabetSize_;
	std::size_t highestPower_;
	std::size_t numKmers_;
	std::size_t kmerSize_;

private:
	void createTable();
	void initializePore();
	void loadModel(const std::string& modelPath);
};

} // namespace dynamont
