#ifndef WAVEFORM_UTILS_H
#define WAVEFORM_UTILS_H
#include <iostream>
#include "TTree.h"
#include "TH1F.h"
#include <map>
#include <string>
#include <TDatime.h>

struct WaveformData {
    Double_t peakHeight;
    Double_t peak_time;
    Double_t avgBaseline;
    Double_t pulse_start;
    Double_t pulse_end;
    Long64_t waveformID;
};

// Function prototypes
std::vector<TH1D*> createHistogramsFromTree(TTree* tree);
void createAveragedHistogram(std::map<Long64_t, TH1F*>& histograms, const std::string& outputFileName);
void ApplyCutToTree(const char* input_filename, const char* tree_name, const char* output_filename, const char* cut);
TH1D* createAverageHistogram_Induction(std::vector<TH1D*> histograms);
#endif // WAVEFORM_UTILS_H
