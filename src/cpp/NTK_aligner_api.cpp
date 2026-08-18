#include "NTK_aligner_api.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace dynamont
{

namespace
{

constexpr double NEG_INF = -std::numeric_limits<double>::infinity();
constexpr double EPSILON = 1e-8;
constexpr double SPARSETHRESHOLD = -0.02227639471; // std::log(0.95);

} // namespace

NTKAligner::NTKAligner(
	const std::string& modelFile,
	PoreType pore,
	int threads,
	std::size_t band):
	Aligner(modelFile, pore, threads, band),
	ppTNm_(NEG_INF),
	ppTNe_(NEG_INF),
	ppTKm_(NEG_INF),
	ppTKe_(NEG_INF)
{
	initializeTransitions();
}

void NTKAligner::initializeTransitions()
{
	double ntMatch;
	double ntExtend;

	switch (pore_)
	{
		case PoreType::RNA002:
			ntMatch = 0.019889650396799997;
			ntExtend = 0.9801103496029998;
			transitions_ = {
				{"a1", 0.019326040280789637}, {"a2", 0.19725479693713352},
				{"p1", 0.1979799841413514}, {"p2", 0.0006135538271005425},
				{"p3", 0.7669801909288386}, {"s1", 0.27034500789657623},
				{"s2", 0.00032463686748883153}, {"s3", 0.02916688206070035},
				{"e1", 1.0}, {"e2", 0.7296549921055607}, {"e3", 0.8020200158564497},
				{"e4", 0.9797333838008437}, {"i1", 2.3852272324574183e-06},
				{"i2", 0.006598130068516047},
			};
			break;
		case PoreType::RNA004:
			ntMatch = 0.031111753637096777;
			ntExtend = 0.9688882463622581;
			transitions_ = {
				{"a1", 0.029709838889618322}, {"a2", 0.2837864344979079},
				{"p1", 0.15353628902814298}, {"p2", 0.0041495012884881655},
				{"p3", 0.47456322874771467}, {"s1", 0.05012685122100474},
				{"s2", 0.0006112333189296363}, {"s3", 0.13506593503589423},
				{"e1", 1.0}, {"e2", 0.949873148779652}, {"e3", 0.8464637109688202},
				{"e4", 0.9654529072452087}, {"i1", 7.651926003806137e-05},
				{"i2", 0.10658440170772512},
			};
			break;
		case PoreType::DNA_R9:
			ntMatch = 1.0;
			ntExtend = 1.0;
			transitions_ = {
				{"a1", 1.0}, {"a2", 1.0}, {"p1", 1.0}, {"p2", 1.0}, {"p3", 1.0},
				{"s1", 1.0}, {"s2", 1.0}, {"s3", 1.0}, {"e1", 1.0}, {"e2", 1.0},
				{"e3", 1.0}, {"e4", 1.0}, {"i1", 1.0}, {"i2", 1.0},
			};
			break;
		case PoreType::DNA_R10_260:
		case PoreType::DNA_R10_400:
			ntMatch = 0.031111753637096777;
			ntExtend = 0.9688882463622581;
			transitions_ = {
				{"a1", 0.029709838889618322}, {"a2", 0.2837864344979079},
				{"p1", 0.15353628902814298}, {"p2", 0.0041495012884881655},
				{"p3", 0.47456322874771467}, {"s1", 0.05012685122100474},
				{"s2", 0.0006112333189296363}, {"s3", 0.13506593503589423},
				{"e1", 1.0}, {"e2", 0.949873148779652}, {"e3", 0.8464637109688202},
				{"e4", 0.9654529072452087}, {"i1", 7.651926003806137e-05},
				{"i2", 0.10658440170772512},
			};
			break;
		default:
			throw std::runtime_error("Unknown pore type");
	}

	ppTNm_ = std::log(ntMatch);
	ppTNe_ = std::log(ntExtend);
	ppTKm_ = std::log(ntMatch);
	ppTKe_ = std::log(ntExtend);

	for (auto& [name, value] : transitions_)
	{
		value = std::log(value);
	}
}

NTKAligner::SparseScore NTKAligner::makeSparseScore()
{
	return {NEG_INF, NEG_INF, NEG_INF, NEG_INF, NEG_INF};
}

const NTKAligner::SparseScore& NTKAligner::stateArrayOrDefault(
	const SparseMatrix& matrix,
	std::size_t idx) const
{
	static const SparseScore empty = makeSparseScore();
	auto it = matrix.find(idx);
	return it == matrix.end() ? empty : it->second;
}

int NTKAligner::scoreHD(std::size_t kmerN, std::size_t kmerK) const
{
	if (kmerN == kmerK)
		return 0;

	int distance = 0;
	for (std::size_t i = 0; i < kmerSize_; ++i)
	{
		const int remN = static_cast<int>(kmerN % static_cast<std::size_t>(alphabetSize_));
		const int remK = static_cast<int>(kmerK % static_cast<std::size_t>(alphabetSize_));
		distance += (remN != remK);
		kmerN /= static_cast<std::size_t>(alphabetSize_);
		kmerK /= static_cast<std::size_t>(alphabetSize_);
	}
	return -2 * distance;
}

double NTKAligner::score(double signal, std::size_t kmerN, std::size_t kmerK) const
{
	return scoreKmer(signal, kmerN) + scoreKmer(signal, kmerK) + scoreHD(kmerN, kmerK);
}

std::vector<std::size_t> NTKAligner::columnArgsort(
	const double* matrix,
	std::size_t columns,
	std::size_t row)
{
	std::vector<std::size_t> idx(columns);
	std::iota(idx.begin(), idx.end(), 0);
	std::stable_sort(
		idx.begin(),
		idx.end(),
		[matrix, columns, row](std::size_t left, std::size_t right)
		{
			return matrix[row * columns + left] > matrix[row * columns + right];
		});
	return idx;
}

void NTKAligner::logP(
	SparseMatrix& logAPSEI,
	const SparseMatrix& forAPSEI,
	const SparseMatrix& backAPSEI,
	double Z,
	const std::vector<std::size_t>& allowedKeys) const
{
	for (const std::size_t tnk : allowedKeys)
	{
		SparseScore values = makeSparseScore();
		const auto& forwardValues = stateArrayOrDefault(forAPSEI, tnk);
		const auto& backwardValues = stateArrayOrDefault(backAPSEI, tnk);
		for (int mat = 0; mat < NUMMAT; ++mat)
		{
			values[mat] = forwardValues[mat] + backwardValues[mat] - Z;
		}
		logAPSEI[tnk] = values;
	}
}

void NTKAligner::logP(
	double* LP,
	const double* forM,
	const double* backM,
	const double* forE,
	const double* backE,
	std::size_t size,
	double Z)
{
	#pragma omp parallel for schedule(static)
	for (std::size_t i = 0; i < size; ++i)
	{
		LP[i] = Aligner::logPlus(
			forM[i] + backM[i] - Z,
			forE[i] + backE[i] - Z);
	}
}

void NTKAligner::ppForTN(
	const double* signal,
	const int* kmerSeq,
	double* M,
	double* E,
	std::size_t T,
	std::size_t N) const
{
	E[0] = 0.0;
	std::size_t tN = 0;
	for (std::size_t t = 1; t < T; ++t)
	{
		tN += N;
		#pragma omp parallel for if (threads_ > 1) schedule(static)
		for (std::size_t n = 1; n < N; ++n)
		{
			const double sc = scoreKmer(signal[t - 1], static_cast<std::size_t>(kmerSeq[n - 1]));
			M[tN + n] = E[tN - N + n - 1] + sc + ppTNm_;
			E[tN + n] = logPlus(M[tN - N + n] + sc, E[tN - N + n] + sc + ppTNe_);
		}
	}
}

void NTKAligner::ppBackTN(
	const double* signal,
	const int* kmerSeq,
	double* M,
	double* E,
	std::size_t T,
	std::size_t N) const
{
	E[T * N - 1] = 0.0;
	std::size_t tN = (T - 1) * N;
	for (std::size_t t = T - 1; t-- > 0;)
	{
		tN -= N;
		#pragma omp parallel for if (threads_ > 1) schedule(static)
		for (std::size_t n = 0; n < N; ++n)
		{
			double ext = NEG_INF;
			if (n + 1 < N)
			{
				ext = M[tN + N + n + 1] +
					scoreKmer(signal[t], static_cast<std::size_t>(kmerSeq[n])) + ppTNm_;
			}
			if (n > 0)
			{
				const double sc = scoreKmer(signal[t], static_cast<std::size_t>(kmerSeq[n - 1]));
				M[tN + n] = E[tN + N + n] + sc;
				ext = logPlus(ext, E[tN + N + n] + sc + ppTNe_);
			}
			E[tN + n] = ext;
		}
	}
}

void NTKAligner::ppForTK(
	const double* signal,
	double* M,
	double* E,
	std::size_t T,
	std::size_t K) const
{
	for (std::size_t k = 0; k < K; ++k)
		E[k] = 0.0;

	std::size_t tK = 0;
	for (std::size_t t = 1; t < T; ++t)
	{
		const std::size_t prevTK = tK;
		tK += K;
		#pragma omp parallel for if (threads_ > 1) schedule(static)
		for (std::size_t k = 0; k < K; ++k)
		{
			double mat = NEG_INF;
			const double sc = scoreKmer(signal[t - 1], k);
			for (std::size_t preKmer = predecessorKmer(k, 0); preKmer < K; preKmer += highestPower_)
			{
				mat = logPlus(mat, E[prevTK + preKmer] + sc + ppTKm_);
			}
			M[tK + k] = mat;
			E[tK + k] = logPlus(M[prevTK + k] + sc, E[prevTK + k] + sc + ppTKe_);
		}
	}
}

void NTKAligner::ppBackTK(
	const double* signal,
	double* M,
	double* E,
	std::size_t T,
	std::size_t K) const
{
	for (std::size_t k = 0; k < K; ++k)
		E[(T - 1) * K + k] = 0.0;

	std::size_t tK = (T - 1) * K;
	for (std::size_t t = T - 1; t-- > 0;)
	{
		const std::size_t nextTK = tK;
		tK -= K;
		#pragma omp parallel for if (threads_ > 1) schedule(static)
		for (std::size_t k = 0; k < K; ++k)
		{
			double ext = NEG_INF;
			const std::size_t startKmer = successorKmer(k, 0);
			const std::size_t endKmer = startKmer + static_cast<std::size_t>(alphabetSize_);
			for (std::size_t sucKmer = startKmer; sucKmer < endKmer; ++sucKmer)
			{
				ext = logPlus(ext, M[nextTK + sucKmer] + scoreKmer(signal[t], sucKmer) + ppTKm_);
			}
			const double sc = scoreKmer(signal[t], k);
			M[tK + k] = E[nextTK + k] + sc;
			E[tK + k] = logPlus(ext, E[nextTK + k] + sc + ppTKe_);
		}
	}
}

void NTKAligner::preProcTN(
	const double* signal,
	const int* kmerSeq,
	ColumnMask& tnMap,
	std::size_t T,
	std::size_t N) const
{
	const std::size_t TN = T * N;
	std::vector<double> forM(TN, NEG_INF);
	std::vector<double> forE(TN, NEG_INF);
	std::vector<double> backM(TN, NEG_INF);
	std::vector<double> backE(TN, NEG_INF);

	ppForTN(signal, kmerSeq, forM.data(), forE.data(), T, N);
	ppBackTN(signal, kmerSeq, backM.data(), backE.data(), T, N);

	const double Zf = forE[TN - 1];
	const double Zb = backE[0];
	if (std::abs(Zf - Zb) / static_cast<double>(TN) > EPSILON || std::isinf(Zf) || std::isinf(Zb))
	{
		throw std::runtime_error("NTK preprocessing TN failed: alignment scores do not match");
	}

	std::vector<double> LP(TN, NEG_INF);
	logP(LP.data(), forM.data(), backM.data(), forE.data(), backE.data(), TN, Zf);

	tnMap.assign(T, {});
	#pragma omp parallel for if (threads_ > 1) schedule(static)
	for (std::size_t t = 0; t < T; ++t)
	{
		double sum = NEG_INF;
		for (const std::size_t n : columnArgsort(LP.data(), N, t))
		{
			tnMap[t].insert(n);
			sum = logPlus(sum, LP[t * N + n]);
			if (sum >= SPARSETHRESHOLD)
				break;
		}
	}
}

void NTKAligner::preProcTK(
	const double* signal,
	ColumnMask& tkMap,
	std::size_t T,
	std::size_t K) const
{
	const std::size_t TK = T * K;
	std::vector<double> forM(TK, NEG_INF);
	std::vector<double> forE(TK, NEG_INF);
	std::vector<double> backM(TK, NEG_INF);
	std::vector<double> backE(TK, NEG_INF);

	ppForTK(signal, forM.data(), forE.data(), T, K);
	ppBackTK(signal, backM.data(), backE.data(), T, K);

	double Zf = NEG_INF;
	double Zb = NEG_INF;
	for (std::size_t k = 0; k < K; ++k)
	{
		Zf = logPlus(Zf, forE[TK - 1 - k]);
		Zb = logPlus(Zb, backE[k]);
	}

	if (std::abs(Zf - Zb) / static_cast<double>(TK) > EPSILON || std::isinf(Zf) || std::isinf(Zb))
	{
		throw std::runtime_error("NTK preprocessing TK failed: alignment scores do not match");
	}

	std::vector<double> LP(TK, NEG_INF);
	logP(LP.data(), forM.data(), backM.data(), forE.data(), backE.data(), TK, Zb);

	tkMap.assign(T, {});
	#pragma omp parallel for if (threads_ > 1) schedule(static)
	for (std::size_t t = 0; t < T; ++t)
	{
		double sum = NEG_INF;
		for (const std::size_t k : columnArgsort(LP.data(), K, t))
		{
			tkMap[t].insert(k);
			sum = logPlus(sum, LP[t * K + k]);
			if (sum >= SPARSETHRESHOLD)
				break;
		}
	}
}

std::vector<std::size_t> NTKAligner::preProcTNK(
	const double* signal,
	const int* kmerSeq,
	std::size_t T,
	std::size_t N,
	std::size_t K) const
{
	ColumnMask tnMap;
	preProcTN(signal, kmerSeq, tnMap, T, N);

	ColumnMask tkMap;
	preProcTK(signal, tkMap, T, K);

	const std::size_t NK = N * K;
	std::vector<std::size_t> allowedKeys;
	for (std::size_t t = 0; t < T; ++t)
	{
		const std::size_t tNK = t * NK;
		for (const std::size_t n : tnMap.at(t))
		{
			if (n == 0)
			{
				if (t == 0)
					allowedKeys.push_back(0);
				continue;
			}

			const std::size_t base = tNK + n * K;
			allowedKeys.push_back(base + static_cast<std::size_t>(kmerSeq[n - 1]));
			for (const std::size_t k : tkMap.at(t))
			{
				allowedKeys.push_back(base + k);
			}
		}
	}

	std::sort(allowedKeys.begin(), allowedKeys.end());
	allowedKeys.erase(std::unique(allowedKeys.begin(), allowedKeys.end()), allowedKeys.end());
	return allowedKeys;
}

void NTKAligner::logF(
	const double* signal,
	const int* kmerSeq,
	SparseMatrix& forAPSEI,
	const std::vector<std::size_t>& allowedKeys,
	std::size_t N,
	std::size_t K) const
{
	const std::size_t NK = N * K;
	std::vector<SparseScore> computed(allowedKeys.size(), makeSparseScore());

	for (std::size_t idx = 0; idx < allowedKeys.size(); ++idx)
	{
		const std::size_t tnk = allowedKeys[idx];
		const std::size_t t = tnk / NK;
		const std::size_t n = (tnk % NK) / K;
		const std::size_t k = tnk % K;

		double a = NEG_INF;
		double p = NEG_INF;
		double s = NEG_INF;
		double e = NEG_INF;
		double i = NEG_INF;

		if (t == 0 && n == 0)
		{
			e = 0.0;
		}
		else if (t > 0 && n > 0)
		{
			const double sc = score(signal[t - 1], static_cast<std::size_t>(kmerSeq[n - 1]), k);
			const std::size_t baseIdx1 = tnk - NK - K - k;
			const std::size_t baseIdx2 = tnk - NK - k;
			const std::size_t baseIdx3 = tnk - NK - K;
			const std::size_t baseIdx4 = tnk - NK;
			const std::size_t baseIdx5 = tnk - K;

			for (std::size_t preKmer = predecessorKmer(k, 0); preKmer < K; preKmer += highestPower_)
			{
				const auto& fromA = stateArrayOrDefault(forAPSEI, baseIdx1 + preKmer);
				a = logPlus(a, fromA[3] + transitions_.at("a1") + sc);
				a = logPlus(a, fromA[4] + transitions_.at("a2") + sc);

				const auto& fromP = stateArrayOrDefault(forAPSEI, baseIdx2 + preKmer);
				p = logPlus(p, fromP[2] + transitions_.at("p1") + sc);
				p = logPlus(p, fromP[3] + transitions_.at("p2") + sc);
				p = logPlus(p, fromP[4] + transitions_.at("p3") + sc);
			}

			const auto& fromS = stateArrayOrDefault(forAPSEI, baseIdx3);
			s = logPlus(s, fromS[1] + transitions_.at("s1") + sc);
			s = logPlus(s, fromS[3] + transitions_.at("s2") + sc);
			s = logPlus(s, fromS[4] + transitions_.at("s3") + sc);

			const auto& fromE = stateArrayOrDefault(forAPSEI, baseIdx4);
			e = logPlus(e, fromE[0] + sc);
			e = logPlus(e, fromE[1] + transitions_.at("e2") + sc);
			e = logPlus(e, fromE[2] + transitions_.at("e3") + sc);
			e = logPlus(e, fromE[3] + transitions_.at("e4") + sc);

			const auto& fromI = stateArrayOrDefault(forAPSEI, baseIdx5);
			i = logPlus(i, fromI[3] + transitions_.at("i1") + sc);
			i = logPlus(i, fromI[4] + transitions_.at("i2") + sc);
		}

		computed[idx] = {a, p, s, e, i};
	}

	for (std::size_t idx = 0; idx < allowedKeys.size(); ++idx)
	{
		forAPSEI[allowedKeys[idx]] = computed[idx];
	}
}

void NTKAligner::logB(
	const double* signal,
	const int* kmerSeq,
	SparseMatrix& backAPSEI,
	const std::vector<std::size_t>& allowedKeys,
	std::size_t T,
	std::size_t N,
	std::size_t K) const
{
	const std::size_t NK = N * K;
	std::vector<SparseScore> computed(allowedKeys.size(), makeSparseScore());

	for (std::size_t reverseIndex = 0; reverseIndex < allowedKeys.size(); ++reverseIndex)
	{
		const std::size_t idx = allowedKeys.size() - 1 - reverseIndex;
		const std::size_t tnk = allowedKeys[idx];
		const std::size_t t = tnk / NK;
		const std::size_t n = (tnk % NK) / K;
		const std::size_t k = tnk % K;

		double a = NEG_INF;
		double p = NEG_INF;
		double s = NEG_INF;
		double e = NEG_INF;
		double i = NEG_INF;

		if (t == T - 1 && n == N - 1)
		{
			e = 0.0;
		}

		if (t < T - 1)
		{
			const std::size_t tnkNK = tnk + NK;
			const std::size_t tnkNKK = tnkNK + K;
			const std::size_t sucKmerBase = successorKmer(k, 0);
			const std::size_t sucKmerEnd = sucKmerBase + static_cast<std::size_t>(alphabetSize_);

			if (n > 0)
			{
				double sc = score(signal[t], static_cast<std::size_t>(kmerSeq[n - 1]), k);
				const auto& fromE = stateArrayOrDefault(backAPSEI, tnkNK);
				a = logPlus(a, fromE[3] + sc);
				p = logPlus(p, fromE[3] + transitions_.at("e2") + sc);
				s = logPlus(s, fromE[3] + transitions_.at("e3") + sc);
				e = logPlus(e, fromE[3] + transitions_.at("e4") + sc);

				for (std::size_t sucKmer = sucKmerBase; sucKmer < sucKmerEnd; ++sucKmer)
				{
					sc = score(signal[t], static_cast<std::size_t>(kmerSeq[n - 1]), sucKmer);
					const auto& fromP = stateArrayOrDefault(backAPSEI, tnkNK - k + sucKmer);
					s = logPlus(s, fromP[1] + transitions_.at("p1") + sc);
					e = logPlus(e, fromP[1] + transitions_.at("p2") + sc);
					i = logPlus(i, fromP[1] + transitions_.at("p3") + sc);
				}
			}

			if (n < N - 1)
			{
				double sc = score(signal[t], static_cast<std::size_t>(kmerSeq[n]), k);
				const auto& fromS = stateArrayOrDefault(backAPSEI, tnkNKK);
				p = logPlus(p, fromS[2] + transitions_.at("s1") + sc);
				e = logPlus(e, fromS[2] + transitions_.at("s2") + sc);
				i = logPlus(i, fromS[2] + transitions_.at("s3") + sc);

				for (std::size_t sucKmer = sucKmerBase; sucKmer < sucKmerEnd; ++sucKmer)
				{
					sc = score(signal[t], static_cast<std::size_t>(kmerSeq[n]), sucKmer);
					const auto& fromA = stateArrayOrDefault(backAPSEI, tnkNKK - k + sucKmer);
					e = logPlus(e, fromA[0] + transitions_.at("a1") + sc);
					i = logPlus(i, fromA[0] + transitions_.at("a2") + sc);
				}
			}
		}

		if (t > 0 && n < N - 1)
		{
			const double sc = score(signal[t - 1], static_cast<std::size_t>(kmerSeq[n]), k);
			const auto& fromI = stateArrayOrDefault(backAPSEI, tnk + K);
			e = logPlus(e, fromI[4] + transitions_.at("i1") + sc);
			i = logPlus(i, fromI[4] + transitions_.at("i2") + sc);
		}

		computed[idx] = {a, p, s, e, i};
	}

	for (std::size_t idx = 0; idx < allowedKeys.size(); ++idx)
	{
		backAPSEI[allowedKeys[idx]] = computed[idx];
	}
}

void NTKAligner::pushTraceSegment(
	std::vector<Segment>& segments,
	char state,
	std::size_t sequencePosition,
	std::size_t signalPosition,
	std::vector<double>& segProb,
	std::size_t k) const
{
	segments.push_back(
		Segment{
			state,
			sequencePosition,
			signalPosition,
			formattedMedian(segProb),
			intToKmer(k),
		});
	segProb.clear();
}

void NTKAligner::traceback(
	std::size_t t,
	std::size_t n,
	std::size_t k,
	SparseMatrix& APSEI,
	SparseMatrix& logAPSEI,
	std::vector<Segment>& segments,
	std::size_t N,
	std::size_t K) const
{
	const std::size_t NK = N * K;
	std::vector<double> segProb;
	std::size_t state = 3;
	std::size_t currentIdx = t * NK + n * K + k;

	while (t)
	{
		const auto& current = stateArrayOrDefault(APSEI, currentIdx);
		const auto& currentLog = stateArrayOrDefault(logAPSEI, currentIdx);

		if (state == 3)
		{
			if (t == 1)
			{
				pushTraceSegment(segments, 'M', kmerSize_ / 2, 0, segProb, k);
				break;
			}

			const std::size_t prevIdx = currentIdx - NK;
			const double sc = current[3];
			const double logScore = currentLog[3];
			segProb.push_back(std::exp(logScore));

			const auto& prev = stateArrayOrDefault(APSEI, prevIdx);
			if (sc == prev[3] + logScore)
				state = 3;
			else if (sc == prev[0] + logScore)
				state = 0;
			else if (sc == prev[2] + logScore)
				state = 2;
			else if (sc == prev[1] + logScore)
				state = 1;

			--t;
			currentIdx = prevIdx;
		}
		else if (state == 0)
		{
			if (t == 1 && n == 1)
			{
				pushTraceSegment(segments, 'M', kmerSize_ / 2, 0, segProb, k);
				break;
			}

			const std::size_t prevIdx = currentIdx - NK - K - k;
			const double sc = current[0];
			const double logScore = currentLog[0];
			segProb.push_back(std::exp(logScore));

			for (std::size_t preKmer = predecessorKmer(k, 0); preKmer < K; preKmer += highestPower_)
			{
				const auto& prev = stateArrayOrDefault(APSEI, prevIdx + preKmer);
				if (sc == prev[3] + logScore)
				{
					pushTraceSegment(segments, 'M', n - 1 + kmerSize_ / 2, t - 1, segProb, k);
					state = 3;
					--t;
					--n;
					k = preKmer;
					currentIdx = prevIdx + preKmer;
					break;
				}
				if (sc == prev[4] + logScore)
				{
					pushTraceSegment(segments, 'M', n - 1 + kmerSize_ / 2, t - 1, segProb, k);
					state = 4;
					--t;
					--n;
					k = preKmer;
					currentIdx = prevIdx + preKmer;
					break;
				}
			}
		}
		else if (state == 1)
		{
			if (t == 1)
			{
				pushTraceSegment(segments, 'P', kmerSize_ / 2, 0, segProb, k);
				break;
			}

			const std::size_t prevIdx = currentIdx - NK - k;
			const double sc = current[1];
			const double logScore = currentLog[1];
			segProb.push_back(std::exp(logScore));

			for (std::size_t preKmer = predecessorKmer(k, 0); preKmer < K; preKmer += highestPower_)
			{
				const auto& prev = stateArrayOrDefault(APSEI, prevIdx + preKmer);
				if (sc == prev[3] + logScore)
				{
					pushTraceSegment(segments, 'P', n - 1 + kmerSize_ / 2, t - 1, segProb, k);
					state = 3;
					--t;
					k = preKmer;
					currentIdx = prevIdx + preKmer;
					break;
				}
				if (sc == prev[2] + logScore)
				{
					pushTraceSegment(segments, 'P', n - 1 + kmerSize_ / 2, t - 1, segProb, k);
					state = 2;
					--t;
					k = preKmer;
					currentIdx = prevIdx + preKmer;
					break;
				}
				if (sc == prev[4] + logScore)
				{
					pushTraceSegment(segments, 'P', n - 1 + kmerSize_ / 2, t - 1, segProb, k);
					state = 4;
					--t;
					k = preKmer;
					currentIdx = prevIdx + preKmer;
					break;
				}
			}
		}
		else if (state == 2)
		{
			if (t == 1 && n == 1)
				break;

			const std::size_t prevIdx = currentIdx - NK - K;
			const double sc = current[2];
			const double logScore = currentLog[2];
			segProb.push_back(std::exp(logScore));
			const auto& prev = stateArrayOrDefault(APSEI, prevIdx);

			if (sc == prev[3] + logScore)
				state = 3;
			else if (sc == prev[1] + logScore)
				state = 1;
			else if (sc == prev[4] + logScore)
				state = 4;

			--t;
			--n;
			currentIdx = prevIdx;
		}
		else if (state == 4)
		{
			if (n == 1)
				break;

			const std::size_t prevIdx = currentIdx - K;
			const double sc = current[4];
			const double logScore = currentLog[4];
			segProb.push_back(std::exp(logScore));
			const auto& prev = stateArrayOrDefault(APSEI, prevIdx);

			if (sc == prev[4] + logScore)
				state = 4;
			else if (sc == prev[3] + logScore)
				state = 3;

			--n;
			currentIdx = prevIdx;
		}
	}

	std::reverse(segments.begin(), segments.end());
}

void NTKAligner::decodeMAP(
	Result& result,
	SparseMatrix& logAPSEI,
	const std::vector<std::size_t>& allowedKeys,
	std::size_t T,
	std::size_t N,
	std::size_t K) const
{
	const std::size_t NK = N * K;
	SparseMatrix APSEI;

	for (const std::size_t tnk : allowedKeys)
	{
		const std::size_t t = tnk / NK;
		const std::size_t n = (tnk % NK) / K;
		const std::size_t k = tnk % K;

		double a = NEG_INF;
		double p = NEG_INF;
		double s = NEG_INF;
		double e = NEG_INF;
		double i = NEG_INF;

		if (t == 0 && n == 0)
		{
			e = 0.0;
		}
		if (t > 0 && n > 0)
		{
			const std::size_t baseIdx1 = tnk - NK - K - k;
			const std::size_t baseIdx2 = tnk - NK - k;
			const std::size_t baseIdx3 = tnk - NK - K;
			const std::size_t baseIdx4 = tnk - NK;
			const std::size_t baseIdx5 = tnk - K;
			const auto& logScore = stateArrayOrDefault(logAPSEI, tnk);

			for (std::size_t preKmer = predecessorKmer(k, 0); preKmer < K; preKmer += highestPower_)
			{
				a = std::max(a, stateArrayOrDefault(APSEI, baseIdx1 + preKmer)[3] + logScore[0]);
				a = std::max(a, stateArrayOrDefault(APSEI, baseIdx1 + preKmer)[4] + logScore[0]);

				p = std::max(p, stateArrayOrDefault(APSEI, baseIdx2 + preKmer)[2] + logScore[1]);
				p = std::max(p, stateArrayOrDefault(APSEI, baseIdx2 + preKmer)[3] + logScore[1]);
				p = std::max(p, stateArrayOrDefault(APSEI, baseIdx2 + preKmer)[4] + logScore[1]);
			}

			s = std::max(s, stateArrayOrDefault(APSEI, baseIdx3)[1] + logScore[2]);
			s = std::max(s, stateArrayOrDefault(APSEI, baseIdx3)[3] + logScore[2]);
			s = std::max(s, stateArrayOrDefault(APSEI, baseIdx3)[4] + logScore[2]);

			e = std::max(e, stateArrayOrDefault(APSEI, baseIdx4)[0] + logScore[3]);
			e = std::max(e, stateArrayOrDefault(APSEI, baseIdx4)[1] + logScore[3]);
			e = std::max(e, stateArrayOrDefault(APSEI, baseIdx4)[2] + logScore[3]);
			e = std::max(e, stateArrayOrDefault(APSEI, baseIdx4)[3] + logScore[3]);

			i = std::max(i, stateArrayOrDefault(APSEI, baseIdx5)[3] + logScore[4]);
			i = std::max(i, stateArrayOrDefault(APSEI, baseIdx5)[4] + logScore[4]);
		}

		APSEI[tnk] = {a, p, s, e, i};
	}

	double best = NEG_INF;
	std::size_t bestK = 0;
	const std::size_t lastDim = (T - 1) * NK + (N - 1) * K;
	for (std::size_t k = 0; k < K; ++k)
	{
		const double candidate = stateArrayOrDefault(APSEI, lastDim + k)[3];
		if (candidate >= best)
		{
			best = candidate;
			bestK = k;
		}
	}

	traceback(T - 1, N - 1, bestK, APSEI, logAPSEI, result.segments, N, K);
}

Result NTKAligner::align(
	const double* signal,
	std::size_t signalLength,
	const std::string& sequence,
	bool calcProbabilities)
{
	validateInput(signalLength, sequence.size());

	Result result;
	std::vector<int> kmers = sequenceToKmers(sequence);
	const std::size_t T = signalLength + 1;
	const std::size_t N = kmers.size() + 1;
	const std::size_t K = numKmers_;
	const std::size_t NK = N * K;

	std::vector<std::size_t> allowedKeys = preProcTNK(signal, kmers.data(), T, N, K);

	SparseMatrix forAPSEI;
	logF(signal, kmers.data(), forAPSEI, allowedKeys, N, K);

	SparseMatrix backAPSEI;
	logB(signal, kmers.data(), backAPSEI, allowedKeys, T, N, K);

	double Zf = NEG_INF;
	double Zb = NEG_INF;
	const std::size_t lastDim = (T - 1) * NK + (N - 1) * K;
	for (std::size_t k = 0; k < K; ++k)
	{
		Zf = logPlus(Zf, stateArrayOrDefault(forAPSEI, lastDim + k)[3]);
		Zb = logPlus(Zb, stateArrayOrDefault(backAPSEI, k)[3]);
	}

	if (std::abs(Zf - Zb) / static_cast<double>(T * NK) >= EPSILON ||
		std::isinf(Zf) || std::isinf(Zb))
	{
		throw std::runtime_error("NTK alignment failed: alignment scores do not match");
	}

	result.Z = Zb;
	if (!calcProbabilities)
		return result;

	SparseMatrix logAPSEI;
	logP(logAPSEI, forAPSEI, backAPSEI, Zb, allowedKeys);
	decodeMAP(result, logAPSEI, allowedKeys, T, N, K);
	return result;
}

} // namespace dynamont
