#pragma once

#include <array>
#include <cstddef>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "aligner.hpp"

namespace dynamont
{

class NTKAligner : public Aligner
{
public:
	NTKAligner(
		const std::string& modelFile,
		PoreType pore,
		int threads = 1,
		std::size_t band = 400);

	Result align(
		const double* signal,
		std::size_t signalLength,
		const std::string& sequence,
		bool calcProbabilities = false) override;

private:
	static constexpr int NUMMAT = 5;
	using SparseScore = std::array<double, NUMMAT>;
	using SparseMatrix = std::unordered_map<std::size_t, SparseScore>;
	using ColumnMask = std::vector<std::unordered_set<std::size_t>>;

	std::unordered_map<std::string, double> transitions_;
	double ppTNm_;
	double ppTNe_;
	double ppTKm_;
	double ppTKe_;

	void initializeTransitions();

	int scoreHD(std::size_t kmerN, std::size_t kmerK) const;
	double score(double signal, std::size_t kmerN, std::size_t kmerK) const;

	static SparseScore makeSparseScore();
	const SparseScore& stateArrayOrDefault(const SparseMatrix& matrix, std::size_t idx) const;

	static std::vector<std::size_t> columnArgsort(
		const double* matrix,
		std::size_t columns,
		std::size_t row);

	void logP(
		SparseMatrix& logAPSEI,
		const SparseMatrix& forAPSEI,
		const SparseMatrix& backAPSEI,
		double Z,
		const std::vector<std::size_t>& allowedKeys) const;

	static void logP(
		double* LP,
		const double* forM,
		const double* backM,
		const double* forE,
		const double* backE,
		std::size_t size,
		double Z);

	void ppForTN(
		const double* signal,
		const int* kmerSeq,
		double* M,
		double* E,
		std::size_t T,
		std::size_t N) const;

	void ppBackTN(
		const double* signal,
		const int* kmerSeq,
		double* M,
		double* E,
		std::size_t T,
		std::size_t N) const;

	void ppForTK(
		const double* signal,
		double* M,
		double* E,
		std::size_t T,
		std::size_t K) const;

	void ppBackTK(
		const double* signal,
		double* M,
		double* E,
		std::size_t T,
		std::size_t K) const;

	void preProcTN(
		const double* signal,
		const int* kmerSeq,
		ColumnMask& tnMap,
		std::size_t T,
		std::size_t N) const;

	void preProcTK(
		const double* signal,
		ColumnMask& tkMap,
		std::size_t T,
		std::size_t K) const;

	std::vector<std::size_t> preProcTNK(
		const double* signal,
		const int* kmerSeq,
		std::size_t T,
		std::size_t N,
		std::size_t K) const;

	void logF(
		const double* signal,
		const int* kmerSeq,
		SparseMatrix& forAPSEI,
		const std::vector<std::size_t>& allowedKeys,
		std::size_t N,
		std::size_t K) const;

	void logB(
		const double* signal,
		const int* kmerSeq,
		SparseMatrix& backAPSEI,
		const std::vector<std::size_t>& allowedKeys,
		std::size_t T,
		std::size_t N,
		std::size_t K) const;

	void pushTraceSegment(
		std::vector<Segment>& segments,
		char state,
		std::size_t sequencePosition,
		std::size_t signalPosition,
		std::vector<double>& segProb,
		std::size_t k) const;

	void traceback(
		std::size_t t,
		std::size_t n,
		std::size_t k,
		SparseMatrix& APSEI,
		SparseMatrix& logAPSEI,
		std::vector<Segment>& segments,
		std::size_t N,
		std::size_t K) const;

	void decodeMAP(
		Result& result,
		SparseMatrix& logAPSEI,
		const std::vector<std::size_t>& allowedKeys,
		std::size_t T,
		std::size_t N,
		std::size_t K) const;
};

} // namespace dynamont
