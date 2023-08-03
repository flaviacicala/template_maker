#include "WaveformManager_2.h"

#include <utility> // for std::pair
#include <limits>  // for std::numeric_limits

#include "TH1D.h"
#include "TSpectrum.h"


WaveformManager_2::WaveformManager_2() {
    // constructor implementation
    // initialize your members here if needed, for example
    histogram = nullptr;
}

WaveformManager_2::~WaveformManager_2() {
    // destructor implementation
    // delete or free your resources if needed
    // be careful with deleting or freeing resources
    // as doing it wrongly can cause crashes or memory leaks
}

void WaveformManager_2::SetHistogram(TH1D* hist) {
    histogram = hist;
}

std::pair<double, double> WaveformManager_2::FindMaxNegativePeak() {
    // Initialize peak parameters
    double peak_x = 0.0;
    double peakHeight = 0.0;
    Int_t nfound = 0;
    Double_t* xpeaks = nullptr;

    // Find peaks using TSpectrum
    TSpectrum *spectrum = new TSpectrum();
    histogram->Scale(-1); // Make all values negative for the search
    nfound = spectrum->Search(histogram,3,"",0.2); //increase second argument to make it identify less peaks

    // Get the positions of the found peaks
    xpeaks = spectrum->GetPositionX();

    // If peaks were found, take the highest one (which is the most negative in our case)
    if (nfound > 0) {
        peak_x = xpeaks[0];
        peakHeight = histogram->GetBinContent(histogram->GetXaxis()->FindBin(peak_x));

        for (int p = 1; p < nfound; p++) {
            double tempHeight = histogram->GetBinContent(histogram->GetXaxis()->FindBin(xpeaks[p]));
            if (tempHeight > peakHeight) {
                peak_x = xpeaks[p];
                peakHeight = tempHeight;
            }
        }
    }

    histogram->Scale(-1); // Convert the values back to positive for further processing

    // Cleanup
    delete spectrum;

    return std::make_pair(peak_x, peakHeight);
}

void WaveformManager_2::CutWaveform(double peakPosition, double cutValue) {
    int peakBin = histogram->FindBin(peakPosition);
    int startBin = peakBin, endBin = peakBin;

    // find the first bin to the left that is smaller than cutValue
    for(int i = peakBin; i >= 1; --i) {
        if(histogram->GetBinContent(i) < cutValue) {
            startBin = i;
            break;
        }
    }

    // find the first bin to the right that is smaller than cutValue
    for(int i = peakBin; i <= histogram->GetNbinsX(); ++i) {
        if(histogram->GetBinContent(i) < cutValue) {
            endBin = i;
            break;
        }
    }

    TH1D* cutHisto = (TH1D*)histogram->Clone(Form("%s_cut_%lld", histogram->GetName()));
    cutHisto->GetXaxis()->SetRange(startBin, endBin);

    // update the histogram pointer
    histogram = cutHisto;
}

std::pair<double, double> WaveformManager_2::FindMaxPositivePeak(double peak_start_x) {
    const int nBinsTotal = histogram->GetNbinsX();
    int endBin = histogram->GetXaxis()->FindBin(peak_start_x);
    int startBin = endBin - 30;

    // Adjust startBin if it goes below the first bin
    if (startBin < 1) {
        startBin = 1;
    }

    // Create a new histogram for the search window
    TH1D* zoomHisto = new TH1D(Form("%s_zoom_%lld", histogram->GetName()), Form("%s_zoom_%lld", histogram->GetTitle()), endBin-startBin+1, histogram->GetBinLowEdge(startBin), histogram->GetBinLowEdge(endBin+1));

    // Copy the content from the original histogram to the zoomed histogram
    for (int iBin = startBin; iBin <= endBin; ++iBin) {
        zoomHisto->SetBinContent(iBin - startBin + 1, histogram->GetBinContent(iBin));
    }

    // Variables to hold the highest peak value and its bin number
    double maxPeakValue = -std::numeric_limits<double>::max();
    int maxPeakBin = -1;

    // Use TSpectrum to find the peaks
    TSpectrum *s = new TSpectrum();
    Int_t nfound = s->Search(zoomHisto);

    // Get the positions of the found peaks
    double *xpeaks = s->GetPositionX();

    // If peaks were found, take the highest one
    if (nfound > 0) {
        maxPeakBin = zoomHisto->GetXaxis()->FindBin(xpeaks[0]);
        maxPeakValue = zoomHisto->GetBinContent(maxPeakBin);

        for (int p = 1; p < nfound; p++) {
            int tempBin = zoomHisto->GetXaxis()->FindBin(xpeaks[p]);
            double tempHeight = zoomHisto->GetBinContent(tempBin);
            if (tempHeight > maxPeakValue) {
                maxPeakBin = tempBin;
                maxPeakValue = tempHeight;
            }
        }
    }

    delete s;
    delete zoomHisto;

    // Check if a peak was found
    if (maxPeakBin != -1) {
        double peakPosition = histogram->GetBinCenter(startBin + maxPeakBin - 1);

        // Return the peak value and position as a pair
        return std::make_pair(peakPosition, maxPeakValue);
    } else {
        // No peak was found, return a pair of zeros
        return std::make_pair(0.0, 0.0);
    }
}

int WaveformManager_2::FindStartBin(int peakBin, double cutValue) {
    // find the first bin to the left that is smaller than cutValue
    for(int i = peakBin; i >= 1; --i) {
        if(histogram->GetBinContent(i) < cutValue) {
            return i;
        }
    }
    return 1; // return the first bin if no bin is found that satisfies the condition
}

int WaveformManager_2::FindEndBin(int peakBin, double cutValue) {
    // find the first bin to the right that is smaller than cutValue
    for(int i = peakBin; i <= histogram->GetNbinsX(); ++i) {
        if(histogram->GetBinContent(i) < cutValue) {
            return i;
        }
    }
    return histogram->GetNbinsX(); // return the last bin if no bin is found that satisfies the condition
}
