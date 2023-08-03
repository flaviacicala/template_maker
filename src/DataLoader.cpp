#include "DataLoader.h"
#include "WaveformManager_2.h"
#include "TFile.h"
#include "TTree.h"
#include <iostream> 


DataLoader::DataLoader(const std::string& filename)
    : m_file_name(filename) {
    // Original file and tree
    root_file = TFile::Open(m_file_name.c_str());
    root_tree = dynamic_cast<TTree*>(root_file->Get("hitdumper/hitdumpertree"));

    // New file and tree
    new_root_file = TFile::Open("/sbnd/app/users/moriarty/workdir_2/my_larsoft/template_maker_final/datafile/treeFile.root");
    if (!new_root_file || new_root_file->IsZombie()) {
        //std::cerr << "Error opening new file: treeFile.root" << std::endl;
        return;
    }

    new_root_tree = dynamic_cast<TTree*>(new_root_file->Get("T"));
    if (!new_root_tree) {
        //std::cerr << "Error getting tree from new file: treeFile.root" << std::endl;
        return;
    }

    // Initialize mhit_tpc_vec and mhit_plane_vec
    getTPCValues();
    getPlaneValues();




  //  root_tree->Print();
  //  new_root_tree->Print();
}

DataLoader::~DataLoader() {
    root_file->Close();
    new_root_file->Close();
}
std::vector<TH1D*> DataLoader::loadWaveforms(size_t eventIndex) {
    if (waveforms.empty()) {
        root_tree->SetBranchStatus("rawadcs", 1); // enable rawadcs branch

        if (eventIndex >= root_tree->GetEntries()) {
            return waveforms;
        }

        // Declare vector to access all waveforms, in all events, from root tree
        std::vector<std::vector<double>> rawadcs;
        root_tree->SetBranchAddress("rawadcs", &rawadcs);
        
        root_tree->GetEntry(0);

        for(size_t waveformIndex = 0; waveformIndex < rawadcs.size(); ++waveformIndex) {
            if (mhit_tpc_vec[waveformIndex] != 0 || mhit_plane_vec[waveformIndex] != 0) {
                TH1D* histo = new TH1D(("waveform_histo" + std::to_string(waveformIndex)).c_str(),
                                        "waveform_histo", 
                                        3400, 
                                        0.,  
                                        3400. );
                int Bin_counter = 0;
                for(const auto& iter : rawadcs[waveformIndex]) {
                    histo->SetBinContent(Bin_counter, iter);
                    ++Bin_counter;
                }
                waveforms.push_back(histo);
            }
        }

        delete root_tree;
        root_tree = nullptr;
    }
    return waveforms;
}


void DataLoader::getTPCValues() {
    // Get the total number of entries in the tree
    Long64_t nEntries = new_root_tree->GetEntries();

    // Prepare a variable to hold the branch value
    int mhit_tpc;

    // Associate the branch with the variable
    new_root_tree->SetBranchAddress("mhit_tpc", &mhit_tpc);

    // Clear the vector to hold all the values
    mhit_tpc_vec.clear();
    mhit_tpc_vec.reserve(nEntries); // reserve space for efficiency

    // Loop over the entries in the tree
    for (Long64_t i = 0; i < nEntries; i++) {
        // Get the i-th entry in the tree
        new_root_tree->GetEntry(i);

        // Add the value of mhit_tpc for this entry to the vector
        mhit_tpc_vec.push_back(mhit_tpc);
    }
}

void DataLoader::getPlaneValues() {

    // Get the total number of entries in the tree
    Long64_t nEntries = new_root_tree->GetEntries();

    // Prepare a variable to hold the branch value
    int mhit_plane;

    // Associate the branch with the variable
    new_root_tree->SetBranchAddress("mhit_plane", &mhit_plane);

    // Clear the vector to hold all the values
    mhit_plane_vec.clear();
    mhit_plane_vec.reserve(nEntries); // reserve space for efficiency

    // Loop over the entries in the tree
    for (Long64_t i = 0; i < nEntries; i++) {
        // Get the i-th entry in the tree
        new_root_tree->GetEntry(i);

        // Add the value of mhit_plane for this entry to the vector
        mhit_plane_vec.push_back(mhit_plane);
    }
}




std::vector<int> DataLoader::getTPCVec() {
    return mhit_tpc_vec;
}

std::vector<int> DataLoader::getPlaneVec() {
    return mhit_plane_vec;
}

