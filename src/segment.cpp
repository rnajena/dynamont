#include <iostream> // cout
#include <fstream> // file io
#include <sstream> // file io
#include <string>
#include <map> // dictionary
#include <tuple>
#include <algorithm> // reverse strings
#include <cmath> // exp
#include <limits> // for inifinity
#include <math.h> // log1p

using namespace std;

//Asserts floating point compatibility at compile time
static_assert(numeric_limits<float>::is_iec559, "IEEE 754 required");

map<string,tuple<float, float>> readKmerModel(string file) {
    ifstream inputFile;
    inputFile.open(file);
    string line, kmer, tmp;
    float mean, stdev;
    map<string,tuple<float, float>> model;

    while(getline(inputFile, line)) { // read line
        stringstream inputString(line); // parse line to inputString
        getline(inputString, kmer, '\t');
        reverse(kmer.begin(), kmer.end()); // 3-5 -> 5-3 orientation
        getline(inputString, tmp, '\t');
        mean = atof(tmp.c_str());
        getline(inputString, tmp, '\t');
        stdev = atof(tmp.c_str());
        model[kmer]=make_tuple(mean, stdev);
    }
    return model;
}

float normal_pdf(float x, float m, float s)
{
    static const float inv_sqrt_2pi = 0.3989422804014327;
    float a = (x - m) / s;

    return inv_sqrt_2pi / s * exp(-0.5f * a * a);
}

float score5mer(float signalpoint, string kmer, map<string,tuple<float, float>>* model) {
    tuple<float, float> kmerModel = model->find(kmer)->second;
    return normal_pdf(signalpoint, get<0>(kmerModel), get<1>(kmerModel));
}

void forward(float* signal, string* seq, float* MM, float* MC, int T, int N, map<string,tuple<float, float>>* model){
    // int T = sizeof(signal)/sizeof(*signal);
    // int N = seq->size();
    // float MM[T*N] = {-1};
    // float MC[T*N] = {-1};
    // TODO start at j = 2? 5mers?
    for(int i=0;i<T;i++){
        for(int j=0;j<N;j++){
            float mm=0;
            if(i>0 && j>0){
                mm+=MC[(i-1)*N+(j-1)]*score5mer(signal[i], seq->substr(j-2, 5), model);
            }
            MM[i*N+j] = mm;
            float mc=0;
            if(i>0){
                mc+=MC[(i-1)*N+j]*score5mer(signal[i], seq->substr(j-2, 5), model);
                mc+=MM[(i-1)*N+j]*score5mer(signal[i], seq->substr(j-2, 5), model);
            }
            if(i==0 && j==0){
                mc+=1;
            }
            MC[i*N+j]=mc;
        }
    }
}

float logplus(float a, float b) {
    // safety check, should not occure
    // if(a==b && isinf(a) && isinf(b)) {
    //     return a;
    // }
    if(a>=b){
        return a + log1p(b-a);
    }
    return b + log1p(a-b);
}

void logForward(float* signal, string* seq, float* MM, float* MC, int T, int N, map<string,tuple<float, float>>* model){
    // int T = sizeof(signal)/sizeof(*signal);
    // int N = seq->size();
    // float MM[T*N] = {-1};
    // float MC[T*N] = {-1};
    // TODO start at j = 2? 5mers?
    for(int i=0;i<T;i++){
        for(int j=0;j<N;j++){
            float mm=-INFINITY;
            if(i>0 && j>0){
                mm = logplus(mm, MC[(i-1)*N+(j-1)]) + score5mer(signal[i], seq->substr(j-2, 5), model);
            }
            MM[i*N+j] = mm;
            float mc=-INFINITY;
            if(i>0){
                mc = logplus(mc, MC[(i-1)*N+j]) + score5mer(signal[i], seq->substr(j-2, 5), model);
                mc = logplus(mc, MM[(i-1)*N+j]) + score5mer(signal[i], seq->substr(j-2, 5), model);
            }
            if(i==0 && j==0){
                mc+=0;
            }
            MC[i*N+j]=mc;
        }
    }
}

int main() {
    // float sig[10] = {1.0,2.0,3.0,2.0,3.0,8.0,7.0,9.0,7.0,9.0};
    // int seq[2] = {0,1};
    // int T = 5 + 1;
    // int N = 2 + 1;
    // forward(MM, MC, T, N);
    // for(int i=0;i<T;i++){
    //     for(int j=0;j<N;j++){
    //         cout<<MM[i*N+j]<<' ';
    //     }
    //     cout<<endl;
    // }
    // map<string,tuple<float, float>> model = readKmerModel();
    // cout<<get<0>(model.find("AAAAA")->second)<<"     "<<get<1>(models.find("AAAAA")->second);
    // map<string,tuple<float, float>> :: iterator i;
    // cout<<"Keys"<<"  &  "<<"Value"<<endl;
    // for (i = models.begin(); i!= models.end(); i++) {
    //     cout<<(*i).first<<"    "<<get<0>((*i).second)<<"    "<<get<1>((*i).second)<<"\n";
    // }
    // int* a;
    // int* b;
    // tuple<int*, int*> t = test();
    // a = get<0>(t);
    // cout<<a[0];
}