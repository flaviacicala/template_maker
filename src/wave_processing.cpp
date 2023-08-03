
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
#include "TMarker.h"
#include "wave_processing.h"
#include "TH1F.h"
#include <fstream>
#include <vector>
#include <sys/stat.h> 
#include <sstream>
#include <string>
#include <TLeaf.h>
#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include "TH1.h"
#include <TGraph.h>





std::pair<double, double> WaveformManager_3::DisplayPeak() {
    double peakHeight = 0.0;
    // Find peaks using TSpectrum
    if (histo == nullptr) {
        std::cout << "Invalid histogram!" << std::endl;
    }
    std::cout << "Display peak IN " << std::endl;
    TSpectrum *spectrum = new TSpectrum();
    std::cout << "TSPECTURM FORMEDDDD " << std::endl;
    nfound = spectrum->Search(histo,2,"",0.2); // Adjust parameters according to your needs
    std::cout << "Display peak1 " << std::endl;
    // Get the positions of the found peaks
    xpeaks = spectrum->GetPositionX();

    // If peaks were found, take the highest one
    if (nfound > 0) 
    {
        peak_x = xpeaks[0];
        peakHeight = histo->GetBinContent(histo->GetXaxis()->FindBin(peak_x));

        for (int p = 1; p < nfound; p++)
        {
            double tempHeight = histo->GetBinContent(histo->GetXaxis()->FindBin(xpeaks[p]));
            if (tempHeight > peakHeight) 
            {
                peak_x = xpeaks[p];
                peakHeight = tempHeight;
            }
        }

        // Now mark the highest peak on the histogram
        //TMarker *marker = new TMarker(peak_x, histo->GetBinContent(histo->GetXaxis()->FindBin(peak_x)), 20);
        //marker->SetMarkerColor(kRed);
        //marker->Draw();
        //std::cout << "Display peak 2 " << std::endl;
    }
    //std::cout << "Display peak 3" << std::endl;

    return std::make_pair(peak_x, peakHeight);
}
void WaveformManager_3::ZoomedWaveform(TH1D*& histogram) {
    const int zoomBins = 30; // 100 bins on each side of the peak

    int nBinsTotal = histogram->GetNbinsX();
    int peakBin = histogram->GetXaxis()->FindBin(peak_x);
    int startBin = std::max(1, peakBin - 15);
    int endBin = std::min(nBinsTotal, peakBin + 15);

    // Create a new histogram for the zoomed peak
    TH1D* zoomHisto = new TH1D(Form("%s_zoom_%lld", histogram->GetName(), iWaveform), Form("%s_zoom_%lld", histogram->GetTitle(), iWaveform), zoomBins, startBin, endBin);

    // Copy the content from the original histogram to the zoomed histogram
    for (int iBin = startBin; iBin <= endBin; ++iBin) {
        zoomHisto->SetBinContent(iBin - startBin + 1, histogram->GetBinContent(iBin));
    }

    // Update the histogram pointer outside the function
    histogram = zoomHisto;
}
void WaveformManager_3::DoublePeakCheck() {
    if (nfound > 1) { // more than one peak found
        for (int i = 0; i < nfound; i++) {
            double peakPos = xpeaks[i];
            double peak_x = histo->GetBinContent(histo->GetXaxis()->FindBin(peakPos));
            // Rest of your code here
        }
    }
}


std::pair<double, double> WaveformManager_1::FindMaxNegativePeak(TH1D* histogram) {
    // Initialize peak parameters
    double peak_x = 0.0;
    double peakHeight = 0.0;
    Int_t nfound = 0;
    Double_t* xpeaks = nullptr;
    
    // Find peaks using TSpectrum
    TSpectrum *spectrum = new TSpectrum();
    histogram->Scale(-1); // Make all values negative for the search
    nfound = spectrum->Search(histogram,2,"",0.2);
    
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

void WaveformManager_3::CutWaveform(TH1D*& histogram, double peakPosition, double cutValue) {
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

    TH1D* cutHisto = (TH1D*)histogram->Clone(Form("%s_cut_%lld", histogram->GetName(), iWaveform));
    cutHisto->GetXaxis()->SetRange(startBin, endBin);

    // update the histogram pointer
    histogram = cutHisto;
}


void WaveformManager_1::FindMaxPositivePeak(TH1D*& histogram, double cut_value) {
    const int nBinsTotal = histogram->GetNbinsX();
    const int peakBin = histogram->GetXaxis()->FindBin(peak_x);
    int startBin = peakBin;
    int endBin = peakBin;

    // Find the start bin (closest to the left of the peak that is less than the cut value)
    while (startBin > 1 && histogram->GetBinContent(startBin) >= cut_value) {
        --startBin;
    }

    // Find the end bin (closest to the right of the peak that is less than the cut value)
    while (endBin < nBinsTotal && histogram->GetBinContent(endBin) >= cut_value) {
        ++endBin;
    }

    const int zoomBins = endBin - startBin + 1;

    // Create a new histogram for the zoomed peak
    TH1D* zoomHisto = new TH1D(Form("%s_zoom_%lld", histogram->GetName(), iWaveform), Form("%s_zoom_%lld", histogram->GetTitle(), iWaveform), zoomBins, startBin, endBin);

    // Copy the content from the original histogram to the zoomed histogram
    for (int iBin = startBin; iBin <= endBin; ++iBin) {
        zoomHisto->SetBinContent(iBin - startBin + 1, histogram->GetBinContent(iBin));
    }

    // Update the histogram pointer outside the function
    histogram = zoomHisto;
}

int WaveformManager_3::FindStartBin(TH1D* histogram, int peakBin, double cutValue) {
    // find the first bin to the left that is smaller than cutValue
    for(int i = peakBin; i >= 1; --i) {
        if(histogram->GetBinContent(i) < cutValue) {
            return i;
        }
    }
    return 1; // return the first bin if no bin is found that satisfies the condition
}

// Function to find the bin where the cut should end
int WaveformManager_3::FindEndBin(TH1D* histogram, int peakBin, double cutValue) {
    // find the first bin to the right that is smaller than cutValue
    for(int i = peakBin; i <= histogram->GetNbinsX(); ++i) {
        if(histogram->GetBinContent(i) < cutValue) {
            return i;
        }
    }
    return histogram->GetNbinsX(); // return the last bin if no bin is found that satisfies the condition
}







void CleaningSlice(TH1D*& histo, Long64_t iWaveform) {
    // Slice the waveform to remove the last portion
    const int nBinsToSlice = 10; // Number of bins to slice from the end
    const int nBinsTotal = histo->GetNbinsX();
    const int nBinsToKeep = nBinsTotal - nBinsToSlice;

    // Create a new histogram to store the sliced waveform
    TH1D* slicedHisto = new TH1D(Form("%s_sliced_%lld", histo->GetName(), iWaveform), Form("%s_sliced_%lld", histo->GetTitle(), iWaveform), nBinsToKeep, 0., nBinsToKeep);

    // Copy the content from the original histogram to the sliced histogram
    for (int iBin = 1; iBin <= nBinsToKeep; ++iBin) {
        slicedHisto->SetBinContent(iBin, histo->GetBinContent(iBin));
    }

    // Delete the original histogram
    delete histo;

    // Replace the original histogram pointer with the sliced histogram pointer
    histo = slicedHisto;

    // Save the content of the original histogram to a temporary variable
    std::vector<double> originalContent;
    originalContent.reserve(nBinsToKeep);
    for (int iBin = 1; iBin <= nBinsToKeep; ++iBin) {
        originalContent.push_back(histo->GetBinContent(iBin));
    }

    // Compare the content of the original histogram to the content of the temporary variable
    for (int iBin = 1; iBin <= nBinsToKeep; ++iBin) {
        if (originalContent[iBin - 1] != slicedHisto->GetBinContent(iBin)) {
            std::cout << "The CleaningSlice() function is corrupting the histogram!" << std::endl;
            break;
        }
    }
}

void AdjustBaseline(TH1D* histo) {

    // Calculate the baseline for the first and last 20 bins
    double baselineStart = CalculateBaseline(histo, 1, 20);
    double baselineEnd = CalculateBaseline(histo, histo->GetNbinsX() - 19, histo->GetNbinsX());

    // Calculate the overall baseline as the average of the two
    double baseline = (baselineStart + baselineEnd) / 2;

    // Subtract the baseline from all bins
    for (int iBin = 1; iBin <= histo->GetNbinsX(); ++iBin) {
        double binContent = histo->GetBinContent(iBin);
        histo->SetBinContent(iBin, binContent - baseline);
    }
}

double CalculateFWHM(TH1D* histo) {
    // calculateFWHM function definition


    int max_bin = histo->GetMaximumBin();
    double max_val = histo->GetBinContent(max_bin);
    double max_x = histo->GetBinCenter(max_bin);

    std::cout << "Max bin: " << max_bin << std::endl;
    std::cout << "Max value: " << max_val << std::endl;
    std::cout << "Max x: " << max_x << std::endl;

    // Find the two bins with half the maximum value
    double fwhm_left = histo->FindFirstBinAbove(max_val / 2.0);
    double fwhm_right = histo->FindLastBinAbove(max_val / 2.0);

    std::cout << "FWHM Left: " << fwhm_left << std::endl;
    std::cout << "FWHM Right: " << fwhm_right << std::endl;

    // Calculate the FWHM
    double fwhm = histo->GetBinCenter(fwhm_right) - histo->GetBinCenter(fwhm_left);

    std::cout << "FWHM: " << fwhm << std::endl;

    return fwhm;
}

double CalculateFWZM(TH1D* histo) {
    // calculateFWZM function definition
      // Find the maximum bin
  int max_bin = histo->GetMaximumBin();
  double max_val = histo->GetBinContent(max_bin);
  double max_x = histo->GetBinCenter(max_bin);

  std::cout << "Max bin: " << max_bin << std::endl;
  std::cout << "Max value: " << max_val << std::endl;
  std::cout << "Max x: " << max_x << std::endl;

  // Initialize left and right bin indices
  int fwzm_left = max_bin;
  int fwzm_right = max_bin;

  // Find the left bin where the value falls to zero
  while (fwzm_left >= 1 && histo->GetBinContent(fwzm_left) > 0.0) {
    --fwzm_left;
  }

  // Find the right bin where the value falls to zero
  while (fwzm_right <= histo->GetNbinsX() && histo->GetBinContent(fwzm_right) > 0.0) {
    ++fwzm_right;
  }

  std::cout << "FWZM Left: " << fwzm_left << std::endl;
  std::cout << "FWZM Right: " << fwzm_right << std::endl;

  // Calculate the FWZM
  double fwzm = histo->GetBinCenter(fwzm_right) - histo->GetBinCenter(fwzm_left);

  std::cout << "FWZM: " << fwzm << std::endl;

  return fwzm;
}

double CalculateBaseline(TH1D* histo, int startBin, int endBin) {
    // calculateBaseline function definition
    double sum = 0;
    int count = 0;
  //  std::cout << "Start bin: " << startBin << " End bin: " << endBin << std::endl;
  //  std::cout << "Sum: " << sum << " Baseline: " << sum / (endBin - startBin + 1) << std::endl;

    for (int i = startBin; i <= endBin; ++i) {
        sum += histo->GetBinContent(i);
        ++count;
    }
    return sum / count;
}

void Normalize(TH1D* histo) {
    Double_t peakAmplitude = histo->GetMaximum();
    int nBins = histo->GetNbinsX();

    for (int iBin = 1; iBin <= nBins; ++iBin) {
        Double_t binContent = histo->GetBinContent(iBin);
        histo->SetBinContent(iBin, binContent / peakAmplitude);
    }

}

TH1D* CalculateAverageHistogram(const std::vector<TH1D*>& histograms) {
    if (histograms.empty()) {
        return nullptr;  // Return null pointer if vector is empty
    }

    int nBins = histograms[0]->GetNbinsX();
    double xMin = histograms[0]->GetXaxis()->GetXmin();
    double xMax = histograms[0]->GetXaxis()->GetXmax();

    TH1D* averageHistogram = new TH1D("AverageHistogram", "Average Histogram", nBins, xMin, xMax);
    int numHistograms = histograms.size();

    for (int iBin = 1; iBin <= nBins; ++iBin) {
        double sumBinContents = 0.0;
        for (const auto& histo : histograms) {
            sumBinContents += histo->GetBinContent(iBin);
        }
        double averageBinContent = sumBinContents / numHistograms;
        averageHistogram->SetBinContent(iBin, averageBinContent);
    }

    return averageHistogram;
}

void WriteHistogramToPDF(const std::vector<TH1D*>& histograms, const char* pdfName) {
  

    // Check if the PDF file exists
    std::ifstream ifile(pdfName);
    bool fileExists = ifile.is_open();
    ifile.close();

    // Create a TCanvas
    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);

    // Open the PDF file in "append" mode
    c1->Print(Form("%s[", pdfName), "pdf"); // First page

    // Draw each histogram on the canvas and save it to the PDF file
    for (size_t i = 0; i < histograms.size(); i++) {
        // Clear the canvas
        c1->Clear();

        // Update the histogram title to include the result
        TString title = histograms[i]->GetTitle();
  
        histograms[i]->SetTitle(title);

        // Draw the histogram
        histograms[i]->Draw();

        // Save the canvas to the PDF file
        c1->Print(pdfName, "pdf");
    }

    // Close the PDF file
    c1->Print(Form("%s]", pdfName), "pdf"); // Last page

    // Clean up
    delete c1;
}

char GetPlane(TH1D* h) {
    if (h->GetEntries() == 0) {
        return 'N';
    }

    // Search for positive peaks
    TSpectrum s;
    Int_t nfoundPositive = s.Search(h);
    Double_t* xposPositive = s.GetPositionX();
    Double_t maxPositivePeak = nfoundPositive > 0 ? h->GetBinContent(h->GetXaxis()->FindBin(xposPositive[0])) : 0;

    // Now, create a copy of histogram and multiply by -1 for negative peaks
    TH1D* hNegative = (TH1D*)h->Clone();
    hNegative->Scale(-1);
    Int_t nfoundNegative = s.Search(hNegative);
    Double_t* xposNegative = s.GetPositionX();
    Double_t maxNegativePeak = nfoundNegative > 0 ? hNegative->GetBinContent(hNegative->GetXaxis()->FindBin(xposNegative[0])) : 0;

    delete hNegative; // Clean up cloned histogram

        // Print out the peaks
    std::cout << "The value of the highest positive peak is approximately: " << maxPositivePeak << std::endl;
    std::cout << "The value of the highest negative peak is approximately: " << maxNegativePeak << std::endl;


    // Now determine the plane
    if (maxPositivePeak > 2 * maxNegativePeak) {
        return 'Y';
    } else if (maxNegativePeak > 2 * maxPositivePeak) {
        return 'U';
    } else if (maxPositivePeak > 0.5 * maxNegativePeak && maxNegativePeak > 0.5 * maxPositivePeak) {
        return 'V';
    } else {
        return 'N';
    }
}

bool FileExists(const std::string& filename) {
    struct stat buf;
    return (stat(filename.c_str(), &buf) == 0);
}

void WriteDataToCSV(const std::vector<int>& waveformIDs, const std::vector<char>& planeLetters, const std::string& fileName) {
    bool exists = FileExists(fileName);

    // Open the file in append mode
    std::ofstream file(fileName, std::ios::app);

    // Write the header line if the file did not exist
    if (!exists) {
        file << "WaveformID,PlaneLetter\n";
    }

    // Write the data
    for (size_t i = 0; i < waveformIDs.size(); ++i) {
        file << waveformIDs[i] << ',' << planeLetters[i] << '\n';
    }
}

std::vector<int> GetWaveformIDsByLetter(const std::string& fileName, char letter) {
    std::vector<int> waveformIDs;
    std::ifstream file(fileName);

    // Check if file was opened successfully
    if (!file.is_open()) {
        std::cerr << "Unable to open file " << fileName << std::endl;
        return waveformIDs;
    }

    std::string line;
    // Skip the header line
    std::getline(file, line);

    // Read data line by line
    while (std::getline(file, line)) {
        std::istringstream ss(line);
        std::string field;

        // Get WaveformID
        std::getline(ss, field, ',');
        int waveformID = std::stoi(field);

        // Get PlaneLetter
        std::getline(ss, field, ',');
        char planeLetter = field[0];

        if (planeLetter == letter) {
            waveformIDs.push_back(waveformID);
        }
    }

    return waveformIDs;
}

 std::vector<Long64_t> extractWaveformIDs(const std::string& file_name, const std::string& tree_name, const std::string& branch_name) {
    std::vector<Long64_t> waveformIDs;

    // Open the ROOT file
    TFile* file = TFile::Open(file_name.c_str());
    if (!file || file->IsZombie()) {
        printf("Error: Cannot open file\n");
        return waveformIDs;
    }

    // Get the tree
    TTree* tree = dynamic_cast<TTree*>(file->Get(tree_name.c_str()));
    if (!tree) {
        printf("Tree not found\n");
        return waveformIDs;
    }

    // Get the branch
    TBranch* branch = tree->GetBranch(branch_name.c_str());
    if (!branch) {
        printf("Branch not found\n");
        return waveformIDs;
    }

    // Get the number of entries
    Long64_t nEntries = branch->GetEntries();

    // Loop over entries and extract waveformIDs
    for (Long64_t i = 0; i < nEntries; ++i) {
        branch->GetEntry(i);
        TLeaf* leaf = branch->GetLeaf(branch_name.c_str());
        Long64_t waveformID = leaf->GetValue();
        waveformIDs.push_back(waveformID);
    }

    // Clean up
    delete file;

    return waveformIDs;
}

double chiSquare(const double *par, TH1D* h1, TH1D* h2) {
    // Create copies of the input histograms so we don't modify the originals
    TH1D* h1_copy = (TH1D*)h1->Clone();
    TH1D* h2_copy = (TH1D*)h2->Clone();

    // Apply scale factor and shift
    h1_copy->Scale(par[1]);
    h1_copy->GetXaxis()->Set(h1_copy->GetXaxis()->GetNbins(), h1_copy->GetXaxis()->GetXmin() + par[0], h1_copy->GetXaxis()->GetXmax() + par[0]);

    // Compute Chi2 value
    double chi2 = h1_copy->Chi2Test(h2_copy, "UU CHI2/NDF");

    // Clean up
    delete h1_copy;
    delete h2_copy;

    return chi2;
}


void fitHistograms(TH1D* template_histogram, TH1D* fitting_histogram, double x_guess, double scale_guess) {
    // Create a new histogram for the modified template
    TH1D* modified_template_histogram = new TH1D(*template_histogram);

    // Initialize minimizer
    ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "");

    // Set tolerance and print level
    minimizer->SetTolerance(0.001);
    minimizer->SetPrintLevel(1);

    // Create functor for chi2
    ROOT::Math::Functor functor([&](const double *par) { return chiSquare(par, modified_template_histogram, fitting_histogram); }, 2);
    minimizer->SetFunction(functor);

    // Initial values: x_shift = 0, scale = 1
    minimizer->SetVariable(0, "x_shift", x_guess+1, 1);
    minimizer->SetVariable(1, "scale", scale_guess, 0.1);
    minimizer->SetVariableLowerLimit(1, 0.0);

    // Perform the fit
    minimizer->Minimize();

    // Print results
    const double *xs = minimizer->X();
    double chi_sq = minimizer->MinValue();
    int N = fitting_histogram->GetNbinsX(); // Number of data points
    int k = 2; // Number of free parameters (x_shift and scale)
    int dof = N - k; // Degrees of freedom
    double reduced_chi_sq = chi_sq / dof;
    std::cout << "SCALE GUESS: " << scale_guess << std::endl; 
    std::cout << "Best fit: x_shift = " << xs[0] << ", scale = " << xs[1]<< std::endl;
    std::cout << "Chi square: " << chi_sq << std::endl;
    std::cout << "Reduced chi square: " << reduced_chi_sq << std::endl;


    TH1F *new_histogram = new TH1F("new_histogram","New Histogram",modified_template_histogram->GetNbinsX(), modified_template_histogram->GetXaxis()->GetXmin()+int(xs[0]), modified_template_histogram->GetXaxis()->GetXmax()+int(xs[0]));


    // Modify the new histogram and apply the shift and scale -- something is going wrong here!!!!!!!!!!!
    for (int i = 1; i <= modified_template_histogram->GetNbinsX(); ++i) {
      //  double x = modified_template_histogram->GetBinCenter(i);
        double bin_height = modified_template_histogram->GetBinContent(i);
       // int bin_number_shifted = modified_template_histogram->GetXaxis()->FindBin(x) + 30;
        std::cout <<"new bin " << i << "bin height " << bin_height << std::endl;
        new_histogram->SetBinContent(i, bin_height * xs[1]);
        
    }
    TCanvas* c = new TCanvas("c", "c", 800, 600);
    // Drawing the histograms
        // Get min and max for x and y across both histograms to ensure full visibility


    // Draw an empty frame with the range covering both histograms
    
    fitting_histogram->SetLineColor(kGreen);
    fitting_histogram->Draw("SAME");


    new_histogram->SetLineColor(kBlue);
    new_histogram->Draw("SAME");

    c->Update(); // Make sure to update the canvas

    c->SaveAs("fit_result2.pdf");


    // Cleanup
    delete c;
    delete minimizer;
    delete modified_template_histogram;
}

void WriteSingleHistogramToPDF(TH1D* histogram, const char* fileName) {
  TCanvas canvas("canvas", "canvas");
  histogram->Draw();

  canvas.Print(Form("%s.pdf", fileName));
}

void drawInterpolation(TH1D *h, const char* filename) {

    int nBins = h->GetNbinsX();
    double *x = new double[nBins];
    double *y = new double[nBins];

    for (int i = 0; i < nBins; ++i) {
        x[i] = h->GetBinCenter(i+1);
        y[i] = h->GetBinContent(i+1);
    }

    TGraph *gr = new TGraph(nBins, x, y);

    // Define a new function which uses a polynomial of degree 5
    TF1 *f1 = new TF1("f1", "pol5", h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());

    // Fit the histogram
    gr->Fit(f1);

    // Draw histogram and function
    TCanvas *c1 = new TCanvas("c1","c1",800,600);
    h->Draw();
    gr->SetLineColor(kRed);
    gr->Draw("SAME C");  // 'C' draws a smooth curve
    f1->SetLineColor(kBlue);
    f1->Draw("SAME");

    // Save as PDF
    c1->SaveAs(filename);

    delete [] x;
    delete [] y;
    delete c1;
    delete gr;
    delete f1;
}



