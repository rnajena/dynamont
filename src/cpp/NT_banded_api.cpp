namespace dynamont {

struct NTOptions
{
    std::string pore;
    std::string modelPath;

    size_t bandwidth = 400;
    int threads = 1;

    bool calcZ = false;
    bool train = false;
    bool probability = false;

    double m1 = -1;
    double e1 = -1;
    double e2 = -1;
};

struct Segment
{
    int base;
    int sample; // signal index
    double probability;
};

struct Result
{
    double Z;

    std::vector<Segment> segments;

    std::vector<double> borderProbabilities;

    std::vector<double> trainedMeans;
    std::vector<double> trainedStddevs;
};

Result runNTBanded(
    const std::vector<double>& signal,
    const std::string& read,
    const NTOptions& options);

}