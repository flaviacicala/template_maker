#include "waveform_utils.h"
#include "wave_processing.h"
#include "TFile.h"
#include "TCanvas.h"
#include <iostream>
#include <TDatime.h>

std::vector<TH1D*> createHistogramsFromTree(TTree* tree) {
    const std::string m_file_name = "/sbnd/app/users/moriarty/workdir_2/my_larsoft/template_maker_final/datafile/good_muon_hitdumper.root";

    // Open root file
    TFile* root_file = TFile::Open( m_file_name.c_str() );
    TTree* root_tree = dynamic_cast<TTree*>( root_file->Get("hitdumper/hitdumpertree" ) );

    // Declare vector to access all waveforms, in all events, from root tree
    std::vector<std::vector<double>> rawadcs;
    auto* rawadcs_ptr = &rawadcs;
    root_tree->SetBranchStatus("*", 0); // disable all branches
    root_tree->SetBranchStatus("rawadcs", 1); // enable specific branch
    root_tree->SetBranchAddress("rawadcs", &rawadcs_ptr);

    root_tree->GetEntry(0);

    std::vector<TH1D*> histograms;
  
    Double_t pulse_start;
    tree->SetBranchAddress("pulse_start", &pulse_start);

    Double_t pulse_end;
    tree->SetBranchAddress("pulse_end", &pulse_end);

    Long64_t waveformID;
    tree->SetBranchAddress("waveformID", &waveformID);

    Double_t peak_time;
    tree->SetBranchAddress("peak_time", &peak_time);

    Long64_t nEntries = tree->GetEntries();
    Long64_t nHistograms = 0;

    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);

        if (pulse_end <= pulse_start) {
            std::cerr << "Warning: pulse_end is not greater than pulse_start for entry " << i << std::endl;
            continue;
        }

        // Compute the number of bins with width 1.
        Long64_t nBins = static_cast<Long64_t>(pulse_end - pulse_start);

        // Check if waveformID is within the range of rawadcs indices
        if (waveformID < 0 || waveformID >= static_cast<Long64_t>(rawadcs.size())) {
            std::cerr << "Warning: waveformID " << waveformID << " is out of rawadcs range" << std::endl;
            continue;
        }

        TH1D* histo { new TH1D(("waveform" + std::to_string(i)).c_str(), ("waveform " + std::to_string(i)).c_str(), nBins, 0.,  nBins) };

        // Fill the histogram with the rawadc values within pulse_start and pulse_end.
        for (Long64_t bin = static_cast<Long64_t>(pulse_start); bin < static_cast<Long64_t>(pulse_end) && bin < static_cast<Long64_t>(rawadcs[waveformID].size()); ++bin) {
            histo->SetBinContent(bin - static_cast<Long64_t>(pulse_start), rawadcs[waveformID][bin]);
        }
        AdjustBaseline(histo);
        CleaningSlice(histo, waveformID);
        Normalize(histo);
        histograms.push_back(histo);
        nHistograms++;
        
        std::cout << "Entry " << i << ": pulse_start = " << pulse_start 
                  << ", pulse_end = " << pulse_end 
                  << ", waveformID = " << waveformID 
                  << ", peak_time = " << peak_time << std::endl;
    }

    std::cout << "Number of created histograms: " << nHistograms << std::endl;

    return histograms;
}




void ApplyCutToTree(const char* input_filename, const char* tree_name, const char* output_filename, const char* cut) {
    // Open the input ROOT file
    TFile *inputFile = TFile::Open(input_filename);
    if (!inputFile || inputFile->IsZombie()) {
        std::cout << "Error: Cannot open file " << input_filename << std::endl;
        return;
    }

    // Get the TTree
    TTree *tree = (TTree*)inputFile->Get(tree_name);
    if (!tree) {
        std::cout << "Error: Cannot find TTree " << tree_name << std::endl;
        return;
    }

    // Open the output ROOT file
    TFile *outputFile = new TFile(output_filename, "RECREATE");
    if (!outputFile || outputFile->IsZombie()) {
        std::cout << "Error: Cannot open file " << output_filename << std::endl;
        return;
    }

    // Copy the tree applying the cut
    TTree *newTree = tree->CopyTree(cut);

    // Write the new TTree to the file
    newTree->Write();

    // Clean up
    delete newTree;
    delete tree;
    delete outputFile;
    delete inputFile;
}

TH1D* createAverageHistogram_Induction(std::vector<TH1D*> histograms) {
    // Find the histogram with the most bins
    TH1D* maxBinHistogram = nullptr;
    int maxBins = 0;
    for (auto histo : histograms) {
        int nBins = histo->GetXaxis()->GetNbins();
        if (nBins > maxBins) {
            maxBins = nBins;
            maxBinHistogram = histo;
        }
    }

    if (!maxBinHistogram) {
        std::cerr << "No input histograms provided." << std::endl;
        return nullptr;
    }

    // Create the resulting histogram with the same bin structure as maxBinHistogram
    TH1D* avgHisto = (TH1D*)maxBinHistogram->Clone("avgHisto");
    avgHisto->Reset();

    // Vector to hold the sum of bin contents for each bin across all histograms
    std::vector<double> sumBinContents(maxBins, 0.0);
    
    for (auto histo : histograms) {
        // Center the histogram on its largest absolute peak
        int peakBin = histo->GetMaximumBin();
        if (std::abs(histo->GetBinContent(histo->GetMinimumBin())) > histo->GetBinContent(peakBin))
            peakBin = histo->GetMinimumBin();

        histo->GetXaxis()->SetRangeUser(peakBin - maxBins/2, peakBin + maxBins/2);

        // Add bin contents to the corresponding sum
        for (int i = 1; i <= maxBins; i++) {
            sumBinContents[i-1] += histo->GetBinContent(i);
        }
    }

    // Set the resulting histogram bins to the average bin contents
    for (int i = 1; i <= maxBins; i++) {
        avgHisto->SetBinContent(i, sumBinContents[i-1]/histograms.size());
    }

    return avgHisto;
}
