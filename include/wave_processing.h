#ifndef HISTOGRAM_FUNCTIONS_H
#define HISTOGRAM_FUNCTIONS_H

#include "TH1F.h"

#include <cmath>
#include <iostream>

#include "TFile.h" 
#include "TRandom3.h" 
#include "TTree.h" 
#include "TCanvas.h" 
#include "TH2D.h" 
#include "TLine.h" 
#include "TF1.h"
#include "TPaveLabel.h" 
#include "TPaveText.h"
#include "TROOT.h"
#include "TStyle.h"
#include <TPaveStats.h>
#include "TSpectrum.h"
#include <fstream>
#include "wave_processing.h"
#include <fstream>
#include <vector>
#include <sys/stat.h>
#include <sstream>
#include <string>
#include <TLeaf.h>
#include "Minuit2/Minuit2Minimizer.h"
#include "TH1D.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "TCanvas.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include "TH1.h"
#include <TGraph.h>



class WaveformManager_1 {
public:
    WaveformManager_1(Long64_t iWaveform, TH1D* histo) 
        : iWaveform(iWaveform), histo(histo), peak_x(0.0), spectrum(new TSpectrum()), nfound(0), xpeaks(nullptr) {}

    std::pair<double, double> FindMaxNegativePeak(TH1D* histogram);
    void FindMaxPositivePeak(TH1D*& histogram, double cut_value);

private:
    Long64_t iWaveform;
    double peak_x;
    double peakHeight;
    TH1D* histo;
    TSpectrum* spectrum;
    Int_t nfound;
    Double_t* xpeaks;
};

class WaveformManager_3 {
public:
    WaveformManager_3(Long64_t iWaveform, TH1D* histo) 
        : iWaveform(iWaveform), histo(histo), peak_x(0.0), spectrum(new TSpectrum()), nfound(0), xpeaks(nullptr) {}

    std::pair<double, double> DisplayPeak();
    void ZoomedWaveform(TH1D*& histogram);
    void DoublePeakCheck();
    void CutWaveform(TH1D*& histogram, double peakPosition, double cutValue);
    int FindStartBin(TH1D* histogram, int peakBin, double cutValue);
    int FindEndBin(TH1D* histogram, int peakBin, double cutValue);

private:
    Long64_t iWaveform;
    double peak_x;
    double peakHeight;
    TH1D* histo;
    TSpectrum* spectrum;
    Int_t nfound;
    Double_t* xpeaks;
};

void CleaningSlice(TH1D*& histo, Long64_t iWaveform);
void AdjustBaseline(TH1D* histo);
double CalculateFWHM(TH1D* histo);
double CalculateFWZM(TH1D* histo);
double CalculateBaseline(TH1D* histo, int start, int end);
char GetPlane(TH1D* histo);
void Normalize(TH1D* histo);
TH1D* CalculateAverageHistogram(const std::vector<TH1D*>& histograms);
void WriteHistogramToPDF(const std::vector<TH1D*>& histograms, const char* pdfName);
void WriteDataToCSV(const std::vector<int>& waveformIDs, const std::vector<char>& planeLetters, const std::string& fileName);
bool FileExists(const std::string& filename);
std::vector<int> GetWaveformIDsByLetter(const std::string& fileName, char letter);
std::vector<Long64_t> extractWaveformIDs(const std::string& file_name, const std::string& tree_name, const std::string& branch_name);
double chiSquare(const double *par, TH1D* h1, TH1D* h2);
void fitHistograms(TH1D* template_histogram, TH1D* fitting_histogram, double x_guess, double scale_guess);
void WriteSingleHistogramToPDF(TH1D* histogram, const char* fileName);
void drawInterpolation(TH1D *h, const char* filename);
#endif // HISTOGRAM_FUNCTIONS_H
