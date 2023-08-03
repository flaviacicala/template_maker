#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TColor.h"
#include <vector>
#include <string>

void heatmap_multi(std::vector<std::string> filenames, std::vector<std::string> treenames, const char* branch1_name, const char* branch2_name) {
    if (filenames.size() != treenames.size()) {
        printf("Error: The size of filenames and treenames should be the same.\n");
        return;
    }
    
    for (size_t i = 0; i < filenames.size(); ++i) {
        // Open the ROOT file
        TFile *file = TFile::Open(filenames[i].c_str());

        // Get the tree
        TTree *tree = (TTree*)file->Get(treenames[i].c_str());

        if (!file || file->IsZombie()) {
            printf("Error: Cannot open file\n");
            continue;
        }

        // Check if the tree exists
        if (!tree) {
            printf("Tree not found\n");
            continue;
        }

        // Variables to hold data from the branches
        Long64_t branch1_value;
        double branch2_value;
        Long64_t min_branch1 = std::numeric_limits<Long64_t>::max();  // Initialize with maximum value
        Long64_t max_branch1 = std::numeric_limits<Long64_t>::min();  // Initialize with minimum value

        // Set the branch addresses
        tree->SetBranchAddress(branch1_name, &branch1_value);
        tree->SetBranchAddress(branch2_name, &branch2_value);

        // Find the minimum and maximum values for branch1
        Long64_t nentries = tree->GetEntries();
        for (Long64_t i = 0; i < nentries; i++) {
            tree->GetEntry(i);
            if (branch1_value < min_branch1) {
                min_branch1 = branch1_value;
            }
            if (branch1_value > max_branch1) {
                max_branch1 = branch1_value;
            }
        }

        // Get the minimum and maximum values of each branch
        double branch2_min = tree->GetMinimum(branch2_name);
        double branch2_max = tree->GetMaximum(branch2_name);

        // Create a 2D histogram (heat map)
        TH2D *hist = new TH2D("hist", Form("Heat Map;%s;%s", branch1_name, branch2_name), 1000, min_branch1, max_branch1, 1000, branch2_min, branch2_max);

        // Fill the histogram
        for (Long64_t i = 0; i < nentries; i++) {
            tree->GetEntry(i);
            hist->Fill(branch1_value, branch2_value);
        }

        // Create a canvas
        TCanvas *c1 = new TCanvas("c1", "Heat Map", 700, 700);
        gStyle->SetOptStat(0); // Hide the statistics box

        // Draw the histogram
        hist->Draw("COLZ"); // Z axis colored

        // Save the canvas as a PDF file
        std::string pdf_filename = filenames[i].substr(0, filenames[i].size() - 5) + ".pdf"; // Change the extension from .root to .pdf
        c1->SaveAs(pdf_filename.c_str());

        // Clean up
        delete c1;
        delete hist;
        delete tree;
        delete file;
    }
}
