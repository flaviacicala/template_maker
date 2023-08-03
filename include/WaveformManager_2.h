#ifndef WAVEFORMMANAGER_2_H
#define WAVEFORMMANAGER_2_H

#include <utility> // for std::pair
#include "TH1D.h"

class WaveformManager_2 {
public:
    WaveformManager_2();
    ~WaveformManager_2();
    
    void SetHistogram(TH1D* histogram);
    std::pair<double, double> FindMaxNegativePeak();
    void CutWaveform(double peakPosition, double cutValue);
    std::pair<double, double> FindMaxPositivePeak(double peak_start_x);
    int FindStartBin(int peakBin, double cutValue);
    int FindEndBin(int peakBin, double cutValue);

private:
    TH1D* histogram;
};

#endif /* WAVEFORMMANAGER_2_H */
