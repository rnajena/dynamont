// author: Jannes Spangenberg
// e-mail: jannes.spangenberg@uni-jena.de
// github: https://github.com/JannesSP
// website: https://jannessp.github.io

#include "NTC.hpp"
#include "version.hpp"

int main(int argc, char *argv[])
{
    // speedup for I/O
    std::ios_base::sync_with_stdio(0);
    std::cin.tie(0);
    std::cout.tie(0);

    bool train, calcZ, prob; // atrain
    std::string pore, modelpath;

    std::cerr << std::fixed << std::showpoint << std::setprecision(11);
    std::cout << std::fixed << std::showpoint << std::setprecision(11);

    argparse::ArgumentParser program("dynamont 3d sparsed", PROJECT_VERSION);
    program.add_argument("-m", "--model").help("Path to kmer model table").required().store_into(modelpath);
    program.add_argument("-r", "--pore").help("Pore used to sequence the data").required().choices("rna_r9", "dna_r9", "rna_rp4", "dna_r10_260bps", "dna_r10_400bps").store_into(pore);

    program.add_argument("-a1", "--alignscore1").help("Transition parameter").default_value(-1.0).scan<'g', double>().store_into(transitions_NTK["a1"]);     // a1
    program.add_argument("-a2", "--alignscore2").help("Transition parameter").default_value(-1.0).scan<'g', double>().store_into(transitions_NTK["a2"]);     // a2
    program.add_argument("-e1", "--extendscore1").help("Transition parameter").default_value(-1.0).scan<'g', double>().store_into(transitions_NTK["e1"]);    // e1
    program.add_argument("-e2", "--extendscore2").help("Transition parameter").default_value(-1.0).scan<'g', double>().store_into(transitions_NTK["e2"]);    // e2
    program.add_argument("-e3", "--extendscore3").help("Transition parameter").default_value(-1.0).scan<'g', double>().store_into(transitions_NTK["e3"]);    // e3
    program.add_argument("-e4", "--extendscore4").help("Transition parameter").default_value(-1.0).scan<'g', double>().store_into(transitions_NTK["e4"]);    // e4
    program.add_argument("-s1", "--sequencescore1").help("Transition parameter").default_value(-1.0).scan<'g', double>().store_into(transitions_NTK["s1"]);  // s1
    program.add_argument("-s2", "--sequencescore2").help("Transition parameter").default_value(-1.0).scan<'g', double>().store_into(transitions_NTK["s2"]);  // s2
    program.add_argument("-s3", "--sequencescore3").help("Transition parameter").default_value(-1.0).scan<'g', double>().store_into(transitions_NTK["s3"]);  // s2
    program.add_argument("-p1", "--polishscore1").help("Transition parameter").default_value(-1.0).scan<'g', double>().store_into(transitions_NTK["p1"]);    // p1
    program.add_argument("-p2", "--polishscore2").help("Transition parameter").default_value(-1.0).scan<'g', double>().store_into(transitions_NTK["p2"]);    // p2
    program.add_argument("-p3", "--polishscore3").help("Transition parameter").default_value(-1.0).scan<'g', double>().store_into(transitions_NTK["p3"]);    // p2
    program.add_argument("-i1", "--insertionscore1").help("Transition parameter").default_value(-1.0).scan<'g', double>().store_into(transitions_NTK["i1"]); // i1
    program.add_argument("-i2", "--insertionscore2").help("Transition parameter").default_value(-1.0).scan<'g', double>().store_into(transitions_NTK["i2"]); // i2
    program.add_argument("--train").help("Switch algorithm to transition and emission parameter training mode").default_value(false).implicit_value(true).store_into(train);
    program.add_argument("-z", "--calcZ").help("Switch algorithm to only calculate Z").default_value(false).implicit_value(true).store_into(calcZ);
    // unused, just here to match the other modes
    program.add_argument("-p", "--probabilty").help("Print out the segment border probability").default_value(false).implicit_value(true).store_into(prob); //.store_into(prob);
    program.add_argument("-t").help("Number of threads to use").default_value(1).scan<'i', int>();

    program.parse_args(argc, argv);

    omp_set_dynamic(1);
    omp_set_num_threads(program.get<int>("-t"));

    // load default and set parameters
    if (pore == "rna_r9")
    {
        rna = true;
        // taken from the trained NT version of dynamont
        ppTNm = log(NT_rna_r9_transitions.at("m1")), ppTNe = log(NT_rna_r9_transitions.at("e2"));
        ppTKm = log(NT_rna_r9_transitions.at("m1")), ppTKe = log(NT_rna_r9_transitions.at("e2"));
        updateTransitions(NTK_rna_r9_transitions, transitions_NTK);
    }
    else if (pore == "dna_r9")
    {
        rna = false;
        // taken from the trained NT version of dynamont
        ppTNm = log(NT_dna_r9_transitions.at("m1")), ppTNe = log(NT_dna_r9_transitions.at("e2"));
        ppTKm = log(NT_dna_r9_transitions.at("m1")), ppTKe = log(NT_dna_r9_transitions.at("e2"));
        updateTransitions(NTK_dna_r9_transitions, transitions_NTK);
    }
    else if (pore == "rna_rp4")
    {
        rna = true;
        // taken from the trained NT version of dynamont
        ppTNm = log(NT_rna_rp4_transitions.at("m1")), ppTNe = log(NT_rna_rp4_transitions.at("e2"));
        ppTKm = log(NT_rna_rp4_transitions.at("m1")), ppTKe = log(NT_rna_rp4_transitions.at("e2"));
        updateTransitions(NTK_rna_rp4_transitions, transitions_NTK);
    }
    else if (pore == "dna_r10_260bps")
    {
        rna = false;
        ppTNm = log(NT_dna_r10_260bps_transitions.at("m1")), ppTNe = log(NT_dna_r10_260bps_transitions.at("e2"));
        ppTKm = log(NT_dna_r10_260bps_transitions.at("m1")), ppTKe = log(NT_dna_r10_260bps_transitions.at("e2"));
        updateTransitions(NTK_dna_r10_260bps_transitions, transitions_NTK);
    }
    else if (pore == "dna_r10_400bps")
    {
        rna = false;
        ppTNm = log(NT_dna_r10_400bps_transitions.at("m1")), ppTNe = log(NT_dna_r10_400bps_transitions.at("e2"));
        ppTKm = log(NT_dna_r10_400bps_transitions.at("m1")), ppTKe = log(NT_dna_r10_400bps_transitions.at("e2"));
        updateTransitions(NTK_dna_r10_400bps_transitions, transitions_NTK);
    }

    checkModelpath(modelpath);
    auto result = readKmerModel(modelpath, rna);
    std::tuple<double, double> *model = std::get<0>(result);
    alphabetSize = std::get<1>(result);
    // polishing dimension K = number of possible kmers
    const std::size_t K = std::get<2>(result);
    kmerSize = std::get<3>(result);
    halfKmerSize = kmerSize / 2;

    stepSize = pow(alphabetSize, kmerSize - 1);
    std::string signal, read;

    // echo 107,107,107.2,108.0,108.9,111.2,105.7,104.3,107.1,105.7,105,105 CAAAAA| src\segment.exe
    // read input, signal and read whitespace separated in single line
    getline(std::cin, signal);
    getline(std::cin, read);

    // exit if wrong input ...
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

    // process signal T: convert std::string to double std::array
    const std::size_t T = count(signal.begin(), signal.end(), ',') + 2; // len(sig) + 1
    checkInput(T, read.size(), kmerSize);
    const std::size_t N = read.size() - kmerSize + 1 + 1; // N is number of kmers in sequence + 1
    // std::cerr << "T: " << T << ", " << "N: " << N << ", " << "K: " << K << ", " << "inputsize: " << TNK << "\n";
    NK = N * K;
    TNK = T * NK;

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
    for (std::size_t n = 0; n < N - 1; ++n)
    {
        kmerSeq[n] = kmer2int(read.substr(n, kmerSize), alphabetSize);
    }

    // deallocate memory
    ss.clear();
    signal.erase();
    read.erase();
    value.erase();

    std::vector<std::size_t> allowedKeys = preProcTNK(sig, kmerSeq, T, N, K, model);
    // std::cerr<<"dense: "<<allowedKeys.size()/double(TNK)<<" ("<<allowedKeys.size()<<" / "<<TNK-allowedKeys.size()<<")"<<"\n"; //", sparse: "<<1-(allowedKeys.size()/double(TNK))<<" ("<<TNK-allowedKeys.size()<<")"<<"\n";
    std::unordered_map<std::size_t, std::array<dproxy, NUMMAT>> forAPSEI;

    // std::cerr<<"forward"<<"\n";
    logF(sig, kmerSeq, forAPSEI, allowedKeys, K, model);
    // std::cerr<<"backward"<<"\n";
    std::unordered_map<std::size_t, std::array<dproxy, NUMMAT>> backAPSEI;
    logB(sig, kmerSeq, backAPSEI, allowedKeys, T, N, K, model);

    double Zf = -INFINITY;
    double Zb = -INFINITY;
    for (std::size_t k = 0; k < K; ++k)
    {
        Zf = logPlus(Zf, forAPSEI[TNK - 1 - k][3]);
        Zb = logPlus(Zb, backAPSEI[k][3]);
    }

    // Numeric error is scaled by input size, Z in forward and backward should match by some numeric error EPSILON
    if (abs(Zf - Zb) / TNK >= EPSILON || std::isinf(Zf) || std::isinf(Zb))
    {
        std::cerr << "Z values between matrices do not match! forZ: " << Zf << ", backZ: " << Zb << ", " << abs(Zf - Zb) / TNK << " > " << EPSILON << std::endl;
        delete[] sig;
        delete[] kmerSeq;
        delete[] model;
        exit(3);
    }

    // std::cerr<<"Zf: "<<Zf<<", Zb: "<<Zb<<", "<<abs(Zf-Zb)/TNK<<" <! "<<EPSILON<<"\n";

    if (calcZ)
    {
        std::cout << Zf << std::endl;
    }
    else
    {
        std::unordered_map<std::size_t, std::array<dproxy, NUMMAT>> logAPSEI;
        logP(logAPSEI, forAPSEI, backAPSEI, Zf, allowedKeys);

        // train both Transitions and Emissions
        if (train)
        {
            trainParams(sig, kmerSeq, forAPSEI, backAPSEI, logAPSEI, allowedKeys, T, N, K, model);
            std::cout << "Z:" << Zf << std::endl;

            // print out segmentation std::string
        }
        else
        {
            // Alignment output
            std::list<std::string> segString;
            getBorders(segString, logAPSEI, allowedKeys, T, N, K);

            for (const auto &seg : segString)
            {
                std::cout << seg;
            }
            std::cout << std::endl;

            // calculate sum of segment border probabilities
            if (prob)
            {
                double sum = -INFINITY;
                std::size_t lastT = T, t;
                for (const std::size_t &tnk : allowedKeys)
                {
                    t = tnk / NK;
                    if (t != lastT)
                    {
                        lastT = t;
                        std::cout << sum << ",";
                        sum = -INFINITY;
                    }
                    for (const int &i : {0, 1})
                    { // sum up prob for new segment in dimension K
                        sum = logPlus(sum, logAPSEI.at(tnk)[i]);
                    }
                }
                std::cout << std::endl;
            }
        }
    }
    delete[] sig;
    delete[] kmerSeq;
    delete[] model;
    return 0;
}