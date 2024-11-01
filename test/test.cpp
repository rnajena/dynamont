#include <iostream>
#include <gtest/gtest.h>

#include <vector>
#include <algorithm>

#include "utils.hpp"
#include "dynamont_NTK.cpp"

TEST(Scorehd, shouldComputeCorrectScoreHDFor5merValues)
{
    alphabet_size = 4; // Example alphabet size
    kmerSize = 5;      // Example kmer size

    std::size_t kmer_N = 94;  // 94_base10 = 1132_base4 = ACCTG
    std::size_t kmer_K = 594; // 594_base10 = 21102_base4 = GCCAG

    int expected_score = -2 * 2; // Expected score based on manual calculation
    int actual_score = scoreHD(kmer_N, kmer_K);

    ASSERT_EQ(expected_score, actual_score);
}

TEST(Scorehd, shouldComputeCorrectScoreHDFor9merValues)
{
    alphabet_size = 4; // Example alphabet size
    kmerSize = 9;      // Example kmer size

    std::size_t kmer_N = 23438;  // 23438_base10  = 011232032_base4 = ACCGTGATG
    std::size_t kmer_K = 180132; // 180132_base10 = 223332210_base4 = GGTTTGGCA

    int expected_score = -2 * 7; // Expected score based on manual calculation
    int actual_score = scoreHD(kmer_N, kmer_K);

    ASSERT_EQ(expected_score, actual_score);
}

TEST(Scorehd, shouldHandleZeroKmerValues)
{
    alphabet_size = 4; // Example alphabet size
    kmerSize = 5;      // Example kmer size

    std::size_t kmer_N = 0;
    std::size_t kmer_K = 0;

    int expected_score = 0; // Both kmers are zero, so no differences
    int actual_score = scoreHD(kmer_N, kmer_K);

    ASSERT_EQ(expected_score, actual_score);
}

TEST(ColumnArgsort, shouldSortIndicesInDescendingOrderWhenGivenAValidMatrixAndColumn)
{
    const std::size_t C = 3;
    const std::size_t t = 1;
    const double matrix[] = {1.0, 3.0, 2.0, 4.0, 6.0, 5.0};

    std::vector<std::size_t> expected = {1, 2, 0};
    std::vector<std::size_t> result = column_argsort(matrix, C, t);

    EXPECT_EQ(result, expected);
}

// Calculate logPlus for positive finite numbers
TEST(Logplus, shouldReturnCorrectLogPlusWhenBothNumbersArePositiveAndFinite)
{
    // Arrr matey, let's see if this logPlus be workin' with positive numbers!
    const double x = 2.0;
    const double y = 3.0;
    const double expected = y + log1p(exp(x - y));
    EXPECT_DOUBLE_EQ(logPlus(x, y), expected);
}

// Handle both x and y as positive infinity
TEST(Logplus, BothPositiveInfinity)
{
    // Infinity ho! Both numbers sail into the endless!
    double x = std::numeric_limits<double>::infinity();
    double y = std::numeric_limits<double>::infinity();
    double result = logPlus(x, y);
    EXPECT_TRUE(std::isinf(result));
    EXPECT_GT(result, 0);
}

// Handle both x and y as negative infinity
TEST(Logplus, BothNegativeInfinity)
{
    // Into the abyss! Both numbers vanish into the deep!
    double x = -std::numeric_limits<double>::infinity();
    double y = -std::numeric_limits<double>::infinity();
    double result = logPlus(x, y);
    EXPECT_TRUE(std::isinf(result));
    EXPECT_LT(result, 0);
}

TEST(Logplus, XPositiveInfinityYNegativeInfinity)
{
    // X sails high while Y sinks low!
    double x = std::numeric_limits<double>::infinity();
    double y = -std::numeric_limits<double>::infinity();
    double result = logPlus(x, y);
    EXPECT_TRUE(std::isinf(result));
    EXPECT_GT(result, 0);
}

TEST(Logplus, YPositiveInfinityXNegativeInfinity)
{
    // X sails high while Y sinks low!
    double y = std::numeric_limits<double>::infinity();
    double x = -std::numeric_limits<double>::infinity();
    double result = logPlus(x, y);
    EXPECT_TRUE(std::isinf(result));
    EXPECT_LT(result, 0);
}

TEST(Successingkmer, ValidInputs)
{
    // Given
    std::size_t currentKmer = 5; // AAACC
    int nextNt = 2;              // + G = AACCG
    int stepSize = 256;
    int alphabet_size = 4;

    // When
    std::size_t result = successingKmer(currentKmer, nextNt, stepSize, alphabet_size);

    // Then
    EXPECT_EQ(result, 22);
}

// Computes new k-mer correctly with valid inputs
TEST(Precessingkmer, ValidInputs)
{
    // Given
    std::size_t currentKmer = 10; // AAAGG
    int priorNt = 2;              // + G = GAAAG
    int stepSize = 256;
    int alphabet_size = 4;

    // When
    std::size_t result = precessingKmer(currentKmer, priorNt, stepSize, alphabet_size);

    // Then
    EXPECT_EQ(result, 514);
}

// Use log2Pi constant correctly in calculations
TEST(LogNormalPdf, shouldUseLog2PiConstantCorrectlyInCalculations)
{
    // Given
    const double x = 1.0;
    const double m = 0.0;
    const double s = 1.0;

    // When
    double result = log_normal_pdf(x, m, s);

    // Then
    EXPECT_NEAR(result, -1.4189385332046727, 1e-9);
}

TEST(LogNormalPdf, shouldReturnNegativeInfinityWhenStandardDeviationIsZero)
{
    // Given
    const double x = 1.0;
    const double m = 0.0;
    const double s = 0.0;

    // When
    double result = log_normal_pdf(x, m, s);

    // Then
    EXPECT_EQ(result, -INFINITY);
}

TEST(LogNormalPdf, shouldHandleVerySmallPositiveStandardDeviationValues)
{
    // Given
    const double x = 1.0;
    const double m = 0.0;
    const double s = 1e-10;

    // When
    double result = log_normal_pdf(x, m, s);

    // Then
    EXPECT_NE(result, -INFINITY);
}

// Converts a positive integer to a string representation using a valid base
TEST(Itoa, shouldConvertPositiveIntegerToStringWhenBaseIsValid)
{
    // Given
    std::size_t value = 10; // AAAGG
    int alphabet_size = 4;
    int kmerSize = 5;

    // When
    std::string result = itoa(value, alphabet_size, kmerSize);

    // Then
    EXPECT_EQ(result, "GGAAA");

    value = 0;
    result = itoa(value, alphabet_size, kmerSize);
    EXPECT_EQ(result, "AAAAA");

    value = pow(alphabet_size, kmerSize) - 1;
    result = itoa(value, alphabet_size, kmerSize);
    EXPECT_EQ(result, "TTTTT");

    kmerSize = 9;

    value = 10;
    result = itoa(value, alphabet_size, kmerSize);
    EXPECT_EQ(result, "GGAAAAAAA");

    value = 0;
    result = itoa(value, alphabet_size, kmerSize);
    EXPECT_EQ(result, "AAAAAAAAA");

    value = pow(alphabet_size, kmerSize) - 1;
    result = itoa(value, alphabet_size, kmerSize);
    EXPECT_EQ(result, "TTTTTTTTT");
}

// Converts a valid 5-mer string to an integer using the BASE2ID map
TEST(Kmer2int, shouldConvertValid5merStringToInteger)
{
    // Given
    std::string kmer = "ACGTG";
    int alphabet_size = 4;

    // When
    int result = kmer2int(kmer, alphabet_size);

    // Then
    EXPECT_EQ(result,
              0 * pow(alphabet_size, 4) +
                  1 * pow(alphabet_size, 3) +
                  2 * pow(alphabet_size, 2) +
                  3 * pow(alphabet_size, 1) +
                  2 * pow(alphabet_size, 0));
}

// Converts a valid 9-mer string to an integer using the BASE2ID map
TEST(Kmer2int, shouldConvertValid9merStringToInteger)
{
    // Given
    int alphabet_size = 9;
    std::string kmer = "ACGTTTGCA";

    // When
    int result = kmer2int(kmer, alphabet_size);

    // Then
    EXPECT_EQ(result,
              0 * pow(alphabet_size, 8) +
                  1 * pow(alphabet_size, 7) +
                  2 * pow(alphabet_size, 6) +
                  3 * pow(alphabet_size, 5) +
                  3 * pow(alphabet_size, 4) +
                  3 * pow(alphabet_size, 3) +
                  2 * pow(alphabet_size, 2) +
                  1 * pow(alphabet_size, 1) +
                  0 * pow(alphabet_size, 0));
}

// Calculates median for an odd-sized vector
TEST(Calculatemedian, OddSizedVector)
{
    // Given
    std::vector<double> vec = {3.0, 1.0, 2.0};

    // When
    double median = calculateMedian(vec);

    // Then
    EXPECT_EQ(median, 2.0);
}

// Calculates median for an even-sized vector
TEST(Calculatemedian, EvenSizedVector)
{
    // Given
    std::vector<double> vec = {4.0, 1.0, 3.0, 2.0};

    // When
    double median = calculateMedian(vec);

    // Then
    EXPECT_EQ(median, 2.5);
}

// Handles an empty vector gracefully
TEST(Calculatemedian, EmptyVector)
{
    // Given
    std::vector<double> vec;

    // When & Then
    EXPECT_THROW(calculateMedian(vec), std::out_of_range);
}

// Correctly fetches default values from default_transitions_vals
TEST(Updatetransitions, ShouldFetchDefaultValuesCorrectlyFromDefaultTransitionsVals)
{
    // Given
    std::unordered_map<std::string, double> default_transitions_vals = {{"a", 5.0}, {"b", 6.0}};
    std::unordered_map<std::string, double> transitions = {{"a", -1.0}, {"b", -1.0}};

    // When
    updateTransitions(default_transitions_vals, transitions);

    // Then
    EXPECT_DOUBLE_EQ(transitions["a"], log(5.0));
    EXPECT_DOUBLE_EQ(transitions["b"], log(6.0));
}

// Applies logarithmic transformation to all transition values
TEST(Updatetransitions, ShouldApplyLogarithmicTransformationToAllTransitionValues)
{
    // Given
    std::unordered_map<std::string, double> default_transitions_vals = {{"a", 2.0}, {"b", 3.0}};
    std::unordered_map<std::string, double> transitions = {{"a", 4.0}, {"b", 9.0}};

    // When
    updateTransitions(default_transitions_vals, transitions);

    // Then
    EXPECT_DOUBLE_EQ(transitions["a"], log(4.0));
    EXPECT_DOUBLE_EQ(transitions["b"], log(9.0));
}

// Handles empty transitions map without errors
TEST(Updatetransitions, ShouldHandleEmptyTransitionsMapWithoutErrors)
{
    // Given
    std::unordered_map<std::string, double> default_transitions_vals = {{"a", 2.0}, {"b", 3.0}};
    std::unordered_map<std::string, double> transitions;

    // When
    updateTransitions(default_transitions_vals, transitions);

    // Then
    EXPECT_TRUE(transitions.empty());
}

// Computes log probabilities correctly for given input arrays
TEST(LogpNTK, shouldComputeLogProbabilitiesCorrectlyWhenGivenValidInputArrays)
{
    // Given
    std::size_t S = 3;
    double Z = 1.0;
    double LP[3];
    dproxy forM[] = {2.0, 3.0, 4.0};
    dproxy backM[] = {1.0, 1.0, 1.0};
    dproxy forE[] = {0.5, 0.5, 0.5};
    dproxy backE[] = {0.5, 0.5, 0.5};

    // When
    logP(LP, forM, backM, forE, backE, S, Z);

    // Then
    EXPECT_NEAR(LP[0], logPlus(2.0, 0.0), 1e-9);
    EXPECT_NEAR(LP[1], logPlus(3.0, 0.0), 1e-9);
    EXPECT_NEAR(LP[2], logPlus(4.0, 0.0), 1e-9);
}

// Computes M and E arrays correctly for valid inputs
TEST(Ppfortn, ValidInputs)
{
    // Given: A valid signal array, kmer sequence, and model vector
    const double sig[] = {0.5, 1.0, 1.5};
    const int kmer_seq[] = {0, 1, 2};
    dproxy M[16];
    dproxy E[16];
    std::vector<std::tuple<double, double>> model = {{0.5, 0.1}, {1.0, 0.2}, {1.5, 0.3}};
    // When : ppForTN is called
    ppForTN(sig, kmer_seq, M, E, 4, 4, model);
    double expected_value_M4 = 1.3836465597893728;
    double expected_value_E4 = -INFINITY;
    double tolerance = 1e-9;
    // Then : M and E arrays should be computed correctly
    EXPECT_NEAR(static_cast<double>(M[5]), expected_value_M4, tolerance);
    EXPECT_EQ(static_cast<double>(E[4]), expected_value_E4);
}

// Correctly initializes E[0] to 0
TEST(Ppfortn, InitializeEZero)
{
    // Given: Any valid input parameters
    const double sig[] = {0.5};
    const int kmer_seq[] = {0};
    dproxy M[4];
    dproxy E[4];
    std::vector<std::tuple<double, double>> model = {{0.5, 0.1}};
    // When : ppForTN is called
    ppForTN(sig, kmer_seq, M, E, 2, 2, model);
    // Then : E[0] should be initialized to zero
    EXPECT_EQ(static_cast<double>(E[0]), 0);
}

// Correctly initializes E[0] to 0
TEST(Ppbacktn, InitializeEZero)
{
    // Given: Any valid input parameters
    const double sig[] = {0.5};
    const int kmer_seq[] = {0};
    dproxy M[4];
    dproxy E[4];
    std::vector<std::tuple<double, double>> model = {{0.5, 0.1}};
    // When : ppForTN is called
    ppBackTN(sig, kmer_seq, M, E, 2, 2, model);
    // Then : E[0] should be initialized to zero
    EXPECT_EQ(static_cast<double>(E[3]), 0);
}

TEST(Ppforbacktn, Zmatches)
{
    // Given: A valid signal array, kmer sequence, and model vector
    const double sig[] = {0.5, 0.5, 1.0, 1.0, 1.5, 1.5};
    const int kmer_seq[] = {0, 1, 2};
    const int T = 7;
    const int N = 4;
    std::vector<std::tuple<double, double>> model = {{0.5, 0.1}, {1.0, 0.2}, {1.5, 0.3}};
    dproxy fM[T * N];
    dproxy fE[T * N];
    dproxy bM[T * N];
    dproxy bE[T * N];
    // When : ppForTN is called
    ppForTN(sig, kmer_seq, fM, fE, 7, 4, model);
    ppBackTN(sig, kmer_seq, bM, bE, 7, 4, model);
    // Then : M and E arrays should be computed correctly
    EXPECT_EQ(static_cast<double>(fE[T * N - 1]), static_cast<double>(bE[0]));
}

TEST(Ppforbacktk, Zmatches)
{
    // Given: A valid signal array, kmer sequence, and model vector
    const double sig[] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
    const int T = 7;
    const int K = 1024;
    std::vector<std::tuple<double, double>> model(K, {0.5, 0.1});
    dproxy fM[T * K];
    dproxy fE[T * K];
    dproxy bM[T * K];
    dproxy bE[T * K];
    alphabet_size = 4;
    kmerSize = 5;
    stepSize = pow(alphabet_size, kmerSize - 1);
    // When : ppForTN is called
    ppForTK(sig, fM, fE, T, K, model);
    ppBackTK(sig, bM, bE, T, K, model);
    // Then : M and E arrays should be computed correctly
    double Zf = -INFINITY;
    double Zb = -INFINITY;
    for (std::size_t k = 0; k < K; ++k)
    {
        Zf = logPlus(Zf, fE[T * K - 1 - k]);
        Zb = logPlus(Zb, bE[k]);
    }

    EXPECT_EQ(Zf, Zb);
}

int main(int argc, char **argv)
{

    testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}