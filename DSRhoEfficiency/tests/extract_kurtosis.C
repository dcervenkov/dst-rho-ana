#include <fstream>

int extract_kurtosis(const char* input, const char* output) {
    TFile file(input);
    TH3F* histo = (TH3F*)file.Get("eff_histo");
    ofstream myfile;
    myfile.open(output);
    myfile << "Histo kurtosis: " << histo->GetKurtosis() << std::endl;
    return 0;
}
