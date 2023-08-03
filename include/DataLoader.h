#ifndef DATALOADER_H
#define DATALOADER_H

#include <string>
#include <vector>
#include "TH1D.h" 
#include "TFile.h"
#include "TTree.h"
#include <iostream>  

class DataLoader {
public:
    DataLoader(const std::string& filename);
    ~DataLoader();
    std::vector<TH1D*> loadWaveforms(size_t eventIndex);  // updated method declaration
    int getTPC(size_t eventIndex);
    int getPlane(size_t eventIndex);
    void getTPCValues();  
    void getPlaneValues();  
    std::vector<int> getTPCVec(); 
    std::vector<int> getPlaneVec(); 


private:
    std::string m_file_name;
    TFile* root_file;
    TFile* new_root_file;
    TTree* root_tree;
    TTree* new_root_tree;
    std::vector<int> mhit_tpc_vec;  
    std::vector<int> mhit_plane_vec; 
    std::vector<double> min_ranges;
    std::vector<double> max_ranges;
    std::vector<TH1D*> waveforms;
};

#endif // DATALOADER_H
