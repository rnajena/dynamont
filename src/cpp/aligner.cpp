#include "dynamont/aligner.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <set>
#include <sstream>
#include <stdexcept>

namespace dynamont
{

Aligner::Aligner(
	const std::string& modelFile,
	PoreType pore,
	int threads,
	std::size_t band):
	pore_(pore),
	rna_(false),
	threads_(threads),
	bandwidth_(band / 2),
	alphabetSize_(0),
	highestPower_(0),
	numKmers_(0),
	kmerSize_(0)
{
	createTable();
	initializePore();
	loadModel(modelFile);

	highestPower_ = 1;
	for (std::size_t i = 1; i < kmerSize_; ++i)
	{
		highestPower_ *= static_cast<std::size_t>(alphabetSize_);
	}
}

TrainingResult Aligner::train(
	const double*,
	std::size_t,
	const std::string&)
{
	throw std::runtime_error("Training is not implemented for this aligner");
}

void Aligner::createTable()
{
	baseLookup_ = [] {
		std::array<int8_t, 256> lookup{};
		for (auto& value : lookup)
			value = -1;

		lookup['A'] = lookup['a'] = 0;
		lookup['C'] = lookup['c'] = 1;
		lookup['G'] = lookup['g'] = 2;
		lookup['T'] = lookup['t'] = lookup['U'] = lookup['u'] = 3;
		lookup['N'] = lookup['n'] = 4;
		return lookup;
	}();
}

void Aligner::initializePore()
{
	switch (pore_)
	{
		case PoreType::RNA002:
			rna_ = true;
			kmerSize_ = 5;
			break;
		case PoreType::RNA004:
			rna_ = true;
			kmerSize_ = 9;
			break;
		case PoreType::DNA_R9:
			rna_ = false;
			kmerSize_ = 5;
			break;
		case PoreType::DNA_R10_260:
		case PoreType::DNA_R10_400:
			rna_ = false;
			kmerSize_ = 9;
			break;
		default:
			throw std::runtime_error("Unknown pore type");
	}
}

void Aligner::loadModel(const std::string& modelPath)
{
	std::ifstream file(modelPath);
	if (!file)
	{
		throw std::runtime_error(
			"Could not open model file, please prove a valid model path " + modelPath);
	}

	std::string line;
	std::getline(file, line);
	std::set<char> alphabet;

	while (std::getline(file, line))
	{
		std::stringstream stream(line);
		std::string kmer;
		std::getline(stream, kmer, '\t');

		if (kmer.size() != kmerSize_)
		{
			throw std::runtime_error("Inconsistent kmer size in model");
		}

		for (char base : kmer)
		{
			alphabet.insert(base);
		}
	}

	alphabetSize_ = static_cast<int>(alphabet.size());
	numKmers_ = static_cast<std::size_t>(std::pow(alphabetSize_, kmerSize_));
	model_.resize(numKmers_);

	file.clear();
	file.seekg(0);

	std::getline(file, line);
	while (std::getline(file, line))
	{
		std::stringstream stream(line);
		std::string kmer;
		std::string mean;
		std::string stdev;
		std::getline(stream, kmer, '\t');
		std::getline(stream, mean, '\t');
		std::getline(stream, stdev, '\t');

		if (rna_)
		{
			std::reverse(kmer.begin(), kmer.end());
		}

		model_[kmerToInt(kmer)] = {std::stod(mean), std::stod(stdev)};
	}
}

void Aligner::validateInput(
	std::size_t signalLength,
	std::size_t readLength) const
{
	if (signalLength < 1)
	{
		throw std::runtime_error("Signal is empty");
	}

	if (readLength < kmerSize_)
	{
		throw std::runtime_error("Sequence shorter than model kmer size");
	}

	const std::size_t kmerCount = readLength - kmerSize_ + 1;
	if (signalLength < 2 * kmerCount)
	{
		throw std::runtime_error("Signal too short compared to sequence");
	}
}

std::vector<int> Aligner::sequenceToKmers(const std::string& sequence) const
{
	if (sequence.size() < kmerSize_)
	{
		throw std::runtime_error("Sequence shorter than kmer size");
	}

	const std::size_t lastKmerIdx = sequence.size() - kmerSize_ + 1;
	std::vector<int> kmers(lastKmerIdx);

	int value = 0;
	for (std::size_t i = 0; i < kmerSize_; ++i)
	{
		const int digit = baseLookup_[static_cast<unsigned char>(sequence[i])];
		if (digit < 0 || digit >= alphabetSize_)
		{
			throw std::runtime_error("Invalid nucleotide: " + std::string(1, sequence[i]));
		}
		value = value * alphabetSize_ + digit;
	}
	kmers[0] = value;

	for (std::size_t i = 1; i < lastKmerIdx; ++i)
	{
		const int left = baseLookup_[static_cast<unsigned char>(sequence[i - 1])];
		const int right = baseLookup_[static_cast<unsigned char>(sequence[i + kmerSize_ - 1])];
		if (right < 0 || right >= alphabetSize_)
		{
			throw std::runtime_error(
				"Invalid nucleotide: " + std::string(1, sequence[i + kmerSize_ - 1]));
		}

		value -= static_cast<int>(left * highestPower_);
		value *= alphabetSize_;
		value += right;
		kmers[i] = value;
	}

	return kmers;
}

int Aligner::kmerToInt(const std::string& sequence) const
{
	int value = 0;
	for (std::size_t i = 0; i < kmerSize_; ++i)
	{
		const int digit = baseLookup_[static_cast<unsigned char>(sequence[i])];
		if (digit < 0 || digit >= alphabetSize_)
		{
			throw std::runtime_error("Invalid nucleotide in k-mer: " + sequence);
		}
		value = value * alphabetSize_ + digit;
	}
	return value;
}

std::string Aligner::intToKmer(std::size_t value) const
{
	const char bases[] = "ACGTN";
	std::string result(kmerSize_, 'A');

	for (std::size_t i = 0; i < kmerSize_; ++i)
	{
		result[kmerSize_ - i - 1] = bases[value % static_cast<std::size_t>(alphabetSize_)];
		value /= static_cast<std::size_t>(alphabetSize_);
	}

	if (rna_)
	{
		std::reverse(result.begin(), result.end());
	}

	return result;
}

double Aligner::scoreKmer(double signalValue, std::size_t kmer) const
{
	const auto& model = model_[kmer];
	return log_normal_pdf(signalValue, model.mean, model.stdev);
}

double Aligner::formattedMedian(std::vector<double>& values) const
{
	if (values.empty())
		return 0.0;

	const std::size_t size = values.size();
	const std::size_t mid = size / 2;
	auto middle = values.begin() + static_cast<std::ptrdiff_t>(mid);
	std::nth_element(values.begin(), middle, values.end());

	if (size % 2 == 1)
		return *middle;

	const double upper = *middle;
	std::nth_element(values.begin(), middle - 1, middle);
	return (*(middle - 1) + upper) / 2.0;
}

std::size_t Aligner::successorKmer(std::size_t currentKmer, int nextNt) const
{
	return (currentKmer % highestPower_) * static_cast<std::size_t>(alphabetSize_) + nextNt;
}

std::size_t Aligner::predecessorKmer(std::size_t currentKmer, int priorNt) const
{
	return (currentKmer / static_cast<std::size_t>(alphabetSize_)) +
		(static_cast<std::size_t>(priorNt) * highestPower_);
}

double Aligner::logPlus(double x, double y)
{
	if (std::isinf(x))
		return y;
	if (std::isinf(y))
		return x;
	if (x < y)
		std::swap(x, y);
	return x + std::log1p(std::exp(y - x));
}

double log_normal_pdf(double x, double mean, double stdev)
{
	const double diff = x - mean;
	const double z = diff / stdev;
	return -0.5 * z * z - std::log(stdev) - 0.5 * std::log(2.0 * M_PI);
}

} // namespace dynamont
