#include <iostream>
#include <gtest/gtest.h>

#include <vector>
#include <algorithm>

#include "utils.hpp"
#include "NTC.hpp"

TEST(Scorehd, shouldComputeCorrectScoreHDFor5merValues)
{
    alphabetSize = 4; // Example alphabet size
    kmerSize = 5;     // Example kmer size

    std::size_t kmerN = 94;  // 94_base10 = 1132_base4 = ACCTG
    std::size_t kmerK = 594; // 594_base10 = 21102_base4 = GCCAG

    int expectedScore = -2 * 2; // Expected score based on manual calculation
    int actualScore = scoreHD(kmerN, kmerK);

    ASSERT_EQ(expectedScore, actualScore);
}

TEST(Scorehd, shouldComputeCorrectScoreHDFor9merValues)
{
    alphabetSize = 4; // Example alphabet size
    kmerSize = 9;     // Example kmer size

    std::size_t kmerN = 23438;  // 23438_base10  = 011232032_base4 = ACCGTGATG
    std::size_t kmerK = 180132; // 180132_base10 = 223332210_base4 = GGTTTGGCA

    int expectedScore = -2 * 7; // Expected score based on manual calculation
    int actualScore = scoreHD(kmerN, kmerK);

    ASSERT_EQ(expectedScore, actualScore);
}

TEST(Scorehd, shouldHandleZeroKmerValues)
{
    alphabetSize = 4; // Example alphabet size
    kmerSize = 5;     // Example kmer size

    std::size_t kmerN = 0;
    std::size_t kmerK = 0;

    int expectedScore = 0; // Both kmers are zero, so no differences
    int actualScore = scoreHD(kmerN, kmerK);

    ASSERT_EQ(expectedScore, actualScore);
}

TEST(ColumnArgsort, shouldSortIndicesInDescendingOrderWhenGivenAValidMatrixAndColumn)
{
    const std::size_t C = 3;
    const std::size_t t = 1;
    const double matrix[] = {1.0, 3.0, 2.0, 4.0, 6.0, 5.0};

    std::vector<std::size_t> expected = {1, 2, 0};
    std::vector<std::size_t> result = columnArgsort(matrix, C, t);

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
    std::size_t nextNt = 2;      // + G = AACCG
    stepSize = 256;
    alphabetSize = 4;

    // When
    std::size_t result = successingKmer(currentKmer, nextNt, stepSize, alphabetSize);

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
    int alphabetSize = 4;

    // When
    std::size_t result = precessingKmer(currentKmer, priorNt, stepSize, alphabetSize);

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
    double result = logNormalPdf(x, m, s);

    // Then
    EXPECT_NEAR(result, -1.4189385332046727, 1e-9);
}

// deactivated check for stddev == 0
// TEST(LogNormalPdf, shouldReturnNegativeInfinityWhenStandardDeviationIsZero)
// {
//     // Given
//     const double x = 1.0;
//     const double m = 0.0;
//     const double s = 0.0;

//     // When
//     double result = logNormalPdf(x, m, s);

//     // Then
//     EXPECT_EQ(result, -INFINITY);
// }

TEST(LogNormalPdf, shouldHandleVerySmallPositiveStandardDeviationValues)
{
    // Given
    const double x = 1.0;
    const double m = 0.0;
    const double s = 1e-10;

    // When
    double result = logNormalPdf(x, m, s);

    // Then
    EXPECT_NE(result, -INFINITY);
}

// Converts a positive integer to a string representation using a valid base
TEST(Itoa, check5merRNA)
{
    // Given
    std::size_t value = 10; // AAAGG
    alphabetSize = 4;
    kmerSize = 5;

    // When
    std::string result = itoa(value, alphabetSize, kmerSize, true);

    // Then
    EXPECT_EQ(result, "GGAAA");

    value = 0;
    result = itoa(value, alphabetSize, kmerSize, true);
    EXPECT_EQ(result, "AAAAA");

    value = pow(alphabetSize, kmerSize) - 1;
    result = itoa(value, alphabetSize, kmerSize, true);
    EXPECT_EQ(result, "TTTTT");
}

TEST(Itoa, check9merRNA)
{
    // Given
    std::size_t value = 10; // AAAGG
    alphabetSize = 4;
    kmerSize = 9;

    // When
    std::string result = itoa(value, alphabetSize, kmerSize, true);
    // Then
    EXPECT_EQ(result, "GGAAAAAAA");

    value = 0;
    result = itoa(value, alphabetSize, kmerSize, true);
    EXPECT_EQ(result, "AAAAAAAAA");

    value = pow(alphabetSize, kmerSize) - 1;
    result = itoa(value, alphabetSize, kmerSize, true);
    EXPECT_EQ(result, "TTTTTTTTT");
}

TEST(Itoa, check5merDNA)
{
    // Given
    std::size_t value = 10; // AAAGG
    alphabetSize = 4;
    kmerSize = 5;

    // When
    std::string result = itoa(value, alphabetSize, kmerSize, false);

    // Then
    EXPECT_EQ(result, "AAAGG");

    value = 0;
    result = itoa(value, alphabetSize, kmerSize, false);
    EXPECT_EQ(result, "AAAAA");

    value = pow(alphabetSize, kmerSize) - 1;
    result = itoa(value, alphabetSize, kmerSize, false);
    EXPECT_EQ(result, "TTTTT");
}

TEST(Itoa, check9merDNA)
{
    // Given
    std::size_t value = 10; // AAAGG
    alphabetSize = 4;
    kmerSize = 9;

    // When
    std::string result = itoa(value, alphabetSize, kmerSize, false);

    // Then
    EXPECT_EQ(result, "AAAAAAAGG");

    value = 0;
    result = itoa(value, alphabetSize, kmerSize, false);
    EXPECT_EQ(result, "AAAAAAAAA");

    value = pow(alphabetSize, kmerSize) - 1;
    result = itoa(value, alphabetSize, kmerSize, false);
    EXPECT_EQ(result, "TTTTTTTTT");
}

// Converts a valid 5-mer string to an integer using the BASE2ID map
TEST(Kmer2int, shouldConvertValid5merStringToInteger)
{
    // Given
    std::string kmer = "ACGTG";
    alphabetSize = 4;

    // When
    int result = kmer2int(kmer, alphabetSize);

    // Then
    EXPECT_EQ(result,
              0 * pow(alphabetSize, 4) +
                  1 * pow(alphabetSize, 3) +
                  2 * pow(alphabetSize, 2) +
                  3 * pow(alphabetSize, 1) +
                  2 * pow(alphabetSize, 0));
}

// Converts a valid 9-mer string to an integer using the BASE2ID map
TEST(Kmer2int, shouldConvertValid9merStringToInteger)
{
    // Given
    std::string kmer = "ACGTTTGCA";
    alphabetSize = 9;

    // When
    int result = kmer2int(kmer, alphabetSize);

    // Then
    EXPECT_EQ(result,
              0 * pow(alphabetSize, 8) +
                  1 * pow(alphabetSize, 7) +
                  2 * pow(alphabetSize, 6) +
                  3 * pow(alphabetSize, 5) +
                  3 * pow(alphabetSize, 4) +
                  3 * pow(alphabetSize, 3) +
                  2 * pow(alphabetSize, 2) +
                  1 * pow(alphabetSize, 1) +
                  0 * pow(alphabetSize, 0));
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

// Correctly fetches default values from defaultVals
TEST(Updatetransitions, ShouldFetchDefaultValuesCorrectlyFromDefaultTransitionsVals)
{
    // Given
    std::map<std::string, double> defaultVals = {{"a", 5.0}, {"b", 6.0}};
    std::map<std::string, double> newVals = {{"a", -1.0}, {"b", -1.0}};

    // When
    updateTransitions(defaultVals, newVals);

    // Then
    EXPECT_DOUBLE_EQ(newVals["a"], log(5.0));
    EXPECT_DOUBLE_EQ(newVals["b"], log(6.0));
}

// Applies logarithmic transformation to all transition values
TEST(Updatetransitions, ShouldApplyLogarithmicTransformationToAllTransitionValues)
{
    // Given
    std::map<std::string, double> defaultVals = {{"a", 2.0}, {"b", 3.0}};
    std::map<std::string, double> newVals = {{"a", 4.0}, {"b", 9.0}};

    // When
    updateTransitions(defaultVals, newVals);

    // Then
    EXPECT_DOUBLE_EQ(newVals["a"], log(4.0));
    EXPECT_DOUBLE_EQ(newVals["b"], log(9.0));
}

// Handles empty transitions map without errors
TEST(Updatetransitions, ShouldHandleEmptyTransitionsMapWithoutErrors)
{
    // Given
    std::map<std::string, double> defaultVals = {{"a", 2.0}, {"b", 3.0}};
    std::map<std::string, double> newVals;

    // When
    updateTransitions(defaultVals, newVals);

    // Then
    EXPECT_TRUE(newVals.empty());
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
    const int kmerSeq[] = {0, 1, 2};
    dproxy M[16];
    dproxy E[16];
    std::tuple<double, double> model[] = {{0.5, 0.1}, {1.0, 0.2}, {1.5, 0.3}};
    // When : ppForTN is called
    ppForTN(sig, kmerSeq, M, E, 4, 4, model);
    double expectedValM4 = 1.3836465597893728;
    double expectedValE4 = -INFINITY;
    double tolerance = 1e-9;
    // Then : M and E arrays should be computed correctly
    EXPECT_NEAR(static_cast<double>(M[5]), expectedValM4, tolerance);
    EXPECT_EQ(static_cast<double>(E[4]), expectedValE4);
}

// Correctly initializes E[0] to 0
TEST(Ppfortn, InitializeEZero)
{
    // Given: Any valid input parameters
    const double sig[] = {0.5};
    const int kmerSeq[] = {0};
    dproxy M[4];
    dproxy E[4];
    std::tuple<double, double> model[] = {{0.5, 0.1}};
    // When : ppForTN is called
    ppForTN(sig, kmerSeq, M, E, 2, 2, model);
    // Then : E[0] should be initialized to zero
    EXPECT_EQ(static_cast<double>(E[0]), 0);
}

// Correctly initializes E[0] to 0
TEST(Ppbacktn, InitializeEZero)
{
    // Given: Any valid input parameters
    const double sig[] = {0.5};
    const int kmerSeq[] = {0};
    dproxy M[4];
    dproxy E[4];
    std::tuple<double, double> model[] = {{0.5, 0.1}};
    // When : ppForTN is called
    ppBackTN(sig, kmerSeq, M, E, 2, 2, model);
    // Then : E[0] should be initialized to zero
    EXPECT_EQ(static_cast<double>(E[3]), 0);
}

TEST(Ppforbacktn, Zmatches)
{
    // Given: A valid signal array, kmer sequence, and model vector
    const double sig[] = {0.5, 0.5, 1.0, 1.0, 1.5, 1.5};
    const int kmerSeq[] = {0, 1, 2};
    const int T = 7;
    const int N = 4;
    std::tuple<double, double> model[] = {{0.5, 0.1}, {1.0, 0.2}, {1.5, 0.3}};
    dproxy fM[T * N];
    dproxy fE[T * N];
    dproxy bM[T * N];
    dproxy bE[T * N];
    // When : ppForTN is called
    ppForTN(sig, kmerSeq, fM, fE, 7, 4, model);
    ppBackTN(sig, kmerSeq, bM, bE, 7, 4, model);
    // Then : M and E arrays should be computed correctly
    EXPECT_EQ(static_cast<double>(fE[T * N - 1]), static_cast<double>(bE[0]));
}

TEST(Ppforbacktk, Zmatches)
{
    // Given: A valid signal array, kmer sequence, and model vector
    const double sig[] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
    const int T = 7;
    const int K = 1024;
    std::tuple<double, double> *model = new std::tuple<double, double>[K];
    for (int k = 0; k < K; ++k)
    {
        model[k] = {0.5, 0.1};
    }
    dproxy fM[T * K];
    dproxy fE[T * K];
    dproxy bM[T * K];
    dproxy bE[T * K];
    alphabetSize = 4;
    kmerSize = 5;
    stepSize = pow(alphabetSize, kmerSize - 1);
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