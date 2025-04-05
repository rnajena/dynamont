// author: Jannes Spangenberg
// e-mail: jannes.spangenberg@uni-jena.de
// github: https://github.com/JannesSP
// website: https://jannessp.github.io

#include "NT.hpp"
#include "version.hpp"

int main(int argc, char *argv[])
{
  // speedup for I/O
  std::ios_base::sync_with_stdio(0);
  std::cin.tie(0);
  std::cout.tie(0);

  bool train, calcZ, prob;
  std::string pore, modelpath;

  std::cerr << std::fixed << std::showpoint << std::setprecision(11);
  std::cout << std::fixed << std::showpoint << std::setprecision(11);

  // Argparser
  argparse::ArgumentParser program("dynamont basic", PROJECT_VERSION);
  // parameters for DP
  program.add_argument("-m", "--model").help("Path to kmer model table").required().store_into(modelpath);
  program.add_argument("-r", "--pore").help("Pore used to sequence the data").required().choices("rna_r9", "dna_r9", "rna_rp4", "dna_r10_260bps", "dna_r10_400bps").store_into(pore);

  program.add_argument("-m1", "--matchscore1").help("Segment transition probability, should be close to (expected number of nucleotdes)/(signal length). Leave at -1 if unset.").default_value(-1.0).scan<'g', double>().store_into(transitions_NT["m1"]);
  program.add_argument("-e1", "--extendscore1").help("First extend probability.").default_value(-1.0).scan<'g', double>().store_into(transitions_NT["e1"]);
  program.add_argument("-e2", "--extendscore2").help("Further extend probability.").default_value(-1.0).scan<'g', double>().store_into(transitions_NT["e2"]);
  program.add_argument("--train").help("Switch algorithm to transition and emission parameter training mode").default_value(false).implicit_value(true).store_into(train);
  program.add_argument("-z", "--calcZ").help("Switch algorithm to only calculate Z").default_value(false).implicit_value(true).store_into(calcZ);
  program.add_argument("-p", "--probabilty").help("Print out the segment border probability").default_value(false).implicit_value(true).store_into(prob);
  program.add_argument("-t").help("Number of threads to use").default_value(1).scan<'i', int>();

  program.parse_args(argc, argv);

  omp_set_dynamic(1);
  omp_set_num_threads(program.get<int>("-t"));

  // load default and set parameters
  if (pore == "rna_r9")
  {
    rna = true;
    // taken from the trained NT version of dynamont
    updateTransitions(NT_rna_r9_transitions, transitions_NT);
  }
  else if (pore == "dna_r9")
  {
    rna = false;
    // taken from the trained NT version of dynamont
    updateTransitions(NT_dna_r9_transitions, transitions_NT);
  }
  else if (pore == "rna_rp4")
  {
    rna = true;
    // taken from the trained NT version of dynamont
    updateTransitions(NT_rna_rp4_transitions, transitions_NT);
  }
  else if (pore == "dna_r10_260bps")
  {
    rna = false;
    updateTransitions(NT_dna_r10_260bps_transitions, transitions_NT);
  }
  else if (pore == "dna_r10_400bps")
  {
    rna = false;
    updateTransitions(NT_dna_r10_400bps_transitions, transitions_NT);
  }

  checkModelpath(modelpath);
  auto result = readKmerModel(modelpath, rna);
  std::tuple<double, double> *model = std::get<0>(result);
  const int alphabetSize = std::get<1>(result);
  const int numKmers = std::get<2>(result);
  const int kmerSize = std::get<3>(result);

  // example
  // 107,107,107.2,108.0,108.9,111.2,105.7,104.3,107.1,105.7,105,105
  // CAAAAA
  // read input, signal and read whitespace separated in single line
  std::string signal, read;
  getline(std::cin, signal);
  getline(std::cin, read);

  if (signal.empty())
  {
    std::cerr << "Signal missing!" << std::endl;
    exit(4);
  }
  else if (read.empty())
  {
    std::cerr << "Read missing!" << std::endl;
    exit(5);
  }

  // process signal: convert std::string to double std::array
  const std::size_t T = count(signal.begin(), signal.end(), ',') + 2; // len(sig) + 1
  checkInput(T, read.size(), kmerSize);
  const std::size_t N = read.size() - kmerSize + 1 + 1; // N is number of kmers in sequence + 1
  const std::size_t TN = T * N;

  double *sig = new double[T - 1];
  std::string value;
  std::stringstream ss(signal);
  int i = 0;
  while (getline(ss, value, ','))
  {
    sig[i++] = stod(value);
  }

  // process read N: convert std::string to int std::array
  int *kmerSeq = new int[N - 1];
#pragma omp parallel for
  for (std::size_t n = 0; n < N - 1; ++n)
  {
    kmerSeq[n] = kmer2int(read.substr(n, kmerSize), alphabetSize);
  }

  // deallocate memory
  ss.clear();
  signal.erase();
  read.erase();
  value.erase();

  // std::cerr<<"T: "<<T<<", "<<"N: "<<N<<", "<<"inputsize: "<<TN<<std::endl;
  // calculate segmentation probabilities, fill forward matrices
  double *forM = new double[TN];
  double *forE = new double[TN];
  // calculate segmentation probabilities, fill backward matrices
  double *backM = new double[TN];
  double *backE = new double[TN];
#pragma omp parallel for
  for (std::size_t i = 0; i < TN; ++i)
  {
    forM[i] = -INFINITY;
    forE[i] = -INFINITY;
    backM[i] = -INFINITY;
    backE[i] = -INFINITY;
  }
  logF_NT(sig, kmerSeq, forM, forE, T, N, model);
  logB_NT(sig, kmerSeq, backM, backE, T, N, model);
  const double Zf = forE[TN - 1];
  const double Zb = backE[0];

  // Numeric error is scaled by input size, Z in forward and backward should match by some numeric error EPSILON
  if (abs(Zf - Zb) / TN > EPSILON || std::isinf(Zf) || std::isinf(Zb))
  {
    // for (std::size_t t = 0; t < T; ++t)
    // {
    //     for (std::size_t n = 0; n < N; ++n)
    //     {
    //         std::cout << logPlus(forM[t * N + n], forE[t * N + n]) << ",";
    //     }
    //     std::cout << "\n";

    //     for (std::size_t n = 0; n < N; ++n)
    //     {
    //         std::cout << logPlus(backM[t * N + n], backE[t * N + n]) << ",";
    //     }
    //     std::cout << "\n";
    // }
    // std::cout << std::flush;

    std::cerr << "Z values between matrices do not match! Zf: " << Zf << ", Zb: " << Zb << ", " << abs(Zf - Zb) / TN << " > " << EPSILON << std::endl;
    delete[] sig;
    delete[] kmerSeq;
    delete[] forM;
    delete[] forE;
    delete[] backM;
    delete[] backE;
    delete[] model;
    exit(3);
  }

  // std::cerr<<"Zf: "<<Zf<<", Zb: "<<Zb<<", "<<abs(Zf-Zb)/TN<<" <! "<<EPSILON<<std::endl;

  if (calcZ)
  {
    std::cout << Zb << std::endl;
    delete[] sig;
    delete[] kmerSeq;
    delete[] forM;
    delete[] forE;
    delete[] backM;
    delete[] backE;
    delete[] model;
  }
  else
  {
    // train both Transitions and Emissions
    if (train)
    {
      trainParams(sig, kmerSeq, forM, forE, backM, backE, T, N, model, alphabetSize, numKmers, kmerSize);
      std::cout << "Z:" << Zb << std::endl;
      delete[] sig;
      delete[] kmerSeq;
      delete[] forM;
      delete[] forE;
      delete[] backM;
      delete[] backE;
      delete[] model;
    }
    else
    {
      delete[] sig;
      delete[] kmerSeq;
      delete[] model;
      double *LPM = new double[TN];
      logP(LPM, forM, backM, Zb, TN); // log probs for segmentation
      delete[] forM;
      delete[] backM;
      double *LPE = new double[TN];
      logP(LPE, forE, backE, Zb, TN); // log probs for extension
      delete[] forE;
      delete[] backE;
      std::list<std::string> segString;
      getBorders(segString, LPM, LPE, T, N, kmerSize);

      for (const auto &seg : segString)
      {
        std::cout << seg;
      }
      std::cout << std::endl;

      // calculate sum of segment probabilities
      if (prob)
      {
        double sum;
        for (std::size_t t = 0; t < T; ++t)
        {
          sum = -INFINITY;
          for (std::size_t n = 0; n < N; ++n)
          {
            sum = logPlus(sum, LPM[t * N + n]);
          }
          std::cout << sum << ",";
        }
        std::cout << std::endl;

        // to plot probability heatmap
        // for (std::size_t t = 0; t < T; ++t)
        // {
        //     for (std::size_t n = 0; n < N; ++n)
        //     {
        //         std::cout << exp(logPlus(LPM[t * N + n], LPE[t * N + n])) << ",";
        //     }
        //     std::cout << "\n";

        //     for (std::size_t n = 0; n < N; ++n)
        //     {
        //         std::cout << logPlus(forM[t * N + n], forE[t * N + n]) << ",";
        //     }
        //     std::cout << "\n";

        //     for (std::size_t n = 0; n < N; ++n)
        //     {
        //         std::cout << logPlus(backM[t * N + n], backE[t * N + n]) << ",";
        //     }
        //     std::cout << "\n";
        // }
        // std::cout << std::flush;
      }
      // Clean up
      delete[] LPM;
      delete[] LPE;
    }
  }
}