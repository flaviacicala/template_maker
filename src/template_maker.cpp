#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "waveform_utils.h"
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <unistd.h>  // Required for sleep() function.
#include "wave_processing.h"

int main() {
    // define the cuts here
    std::string cutCase1 = "peakHeight < 50";
    std::string cutCase2 = "peakHeight < 50";

    std::array<std::string, 2> cuts = {cutCase1, cutCase2};

    // define the filenames here
    std::vector<std::string> filenames = {
        "/sbnd/app/users/moriarty/workdir_2/my_larsoft/template_maker_final/build/dataCase1.root",
        "/sbnd/app/users/moriarty/workdir_2/my_larsoft/template_maker_final/build/dataCase2.root",
    };

    for (size_t i = 0; i < filenames.size(); ++i) {
        std::string inputFile = filenames[i];
        std::string treeName = "Tree" + std::to_string(i+1);
        std::string cutOutputFile = "cutData" + std::to_string(i+1) + ".root";

        // Apply the cut to the tree and get the new tree.
        ApplyCutToTree(inputFile.c_str(), treeName.c_str(), cutOutputFile.c_str(), cuts[i].c_str());

        TFile::Open(cutOutputFile.c_str())->Close();  // Close the file.

        sleep(2);  // Wait for 10 seconds.

        TFile* inFileCut = TFile::Open(cutOutputFile.c_str(), "READ");

          TTree* cutTree = dynamic_cast<TTree*>( inFileCut->Get(treeName.c_str()) );



    
        // Create histograms from the cut tree.
        auto histograms = createHistogramsFromTree(cutTree);
        

        // Print the number of bins for each histogram
        for (auto histogram : histograms) {
            std::cout << "Histogram has " << histogram->GetNbinsX() << " bins and max value: " << histogram->GetMaximum() << std::endl;
        }
        auto average_histogram = CalculateAverageHistogram(histograms);
        std::string pdfName = "average_histo_" + std::to_string(i+1) + ".pdf";

        WriteSingleHistogramToPDF(average_histogram, pdfName.c_str());


        for (auto histogram : histograms) {
    delete histogram;
}


        inFileCut->Close();
        delete inFileCut;
    }

    std::cout << "Finished processing all files." << std::endl;

    return 0;
}
