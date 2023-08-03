#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "WaveformManager_2.h"
#include "wave_processing.h"
#include "heatmap_multi.h"

#include "DataLoader.h"


#include <string>
#include <utility> // for std::pair
#include <limits>  // for std::numeric_limits

struct WaveformData {
    Double_t peakHeight;
    Double_t peak_time;
    Double_t avgBaseline;
    Double_t pulse_start;
    Double_t pulse_end;
    Long64_t waveformID;
};


int main()
{
std::vector<std::string> filenames = {
    "dataCase1.root",
    "dataCase2.root",
    "dataCase3.root",
    "dataCase4.root"
};
std::vector<std::string> treenames = {
    "Tree1", 
    "Tree2", 
    "Tree3", 
    "Tree4"
};


std::vector<TH1D*> waveforms;
TFile *fileCase1 = new TFile("dataCase1.root", "RECREATE"); 
TFile *fileCase2 = new TFile("dataCase2.root", "RECREATE"); 
TFile *fileCase3 = new TFile("dataCase3.root", "RECREATE"); 
TFile *fileCase4 = new TFile("dataCase4.root", "RECREATE"); 

WaveformData dataCase1, dataCase2, dataCase3, dataCase4;

TTree* treeCase1 = new TTree("Tree1", "Data from Case 1");
treeCase1->Branch("peakHeight", &dataCase1.peakHeight);
treeCase1->Branch("peak_time", &dataCase1.peak_time);
treeCase1->Branch("avgBaseline", &dataCase1.avgBaseline);
treeCase1->Branch("pulse_start", &dataCase1.pulse_start);
treeCase1->Branch("pulse_end", &dataCase1.pulse_end);
treeCase1->Branch("waveformID", &dataCase1.waveformID);
treeCase1->SetDirectory(fileCase1);



TTree* treeCase2 = new TTree("Tree2", "Data from Case 2");
treeCase2->Branch("peakHeight", &dataCase2.peakHeight);
treeCase2->Branch("peak_time", &dataCase2.peak_time);
treeCase2->Branch("avgBaseline", &dataCase2.avgBaseline);
treeCase2->Branch("pulse_start", &dataCase2.pulse_start);
treeCase2->Branch("pulse_end", &dataCase2.pulse_end);
treeCase2->Branch("waveformID", &dataCase2.waveformID);
treeCase2->SetDirectory(fileCase2);



TTree* treeCase3 = new TTree("Tree3", "Data from Case 3");
treeCase3->Branch("peakHeight", &dataCase3.peakHeight);
treeCase3->Branch("peak_time", &dataCase3.peak_time);
treeCase3->Branch("avgBaseline", &dataCase3.avgBaseline);
treeCase3->Branch("pulse_start", &dataCase3.pulse_start);
treeCase3->Branch("pulse_end", &dataCase3.pulse_end);
treeCase3->Branch("waveformID", &dataCase3.waveformID);
treeCase3->SetDirectory(fileCase3);



TTree* treeCase4 = new TTree("Tree4", "Data from Case 4");
treeCase4->Branch("peakHeight", &dataCase4.peakHeight);
treeCase4->Branch("peak_time", &dataCase4.peak_time);
treeCase4->Branch("avgBaseline", &dataCase4.avgBaseline);
treeCase4->Branch("pulse_start", &dataCase4.pulse_start);
treeCase4->Branch("pulse_end", &dataCase4.pulse_end);
treeCase4->Branch("waveformID", &dataCase4.waveformID);
treeCase4->SetDirectory(fileCase4);



const std::string m_file_name = "/sbnd/app/users/moriarty/workdir_2/my_larsoft/template_maker_final/datafile/good_muon_hitdumper.root";

// Open root file
  TFile* root_file = TFile::Open( m_file_name.c_str() );
    TTree* root_tree = dynamic_cast<TTree*>( root_file->Get("hitdumper/hitdumpertree" ) );

root_tree->SetBranchStatus("*", 0); // disable all branches
    root_tree->SetBranchStatus("rawadcs", 1); // enable specific branches
    root_tree->SetBranchStatus("min_ranges", 1);
    root_tree->SetBranchStatus("max_ranges", 1);



  //declare vector to access all waveforms, in all events, from root tree
  std::vector<std::vector<double>> rawadcs;
  auto* rawadcs_ptr = &rawadcs;

  //link the vector to the relevant branch in the root tree
  root_tree->SetBranchAddress("rawadcs", &rawadcs_ptr);

  //access the smallest amplitude value of each waveform
  std::vector< double > min_ranges;
  auto* min_ranges_ptr = &min_ranges;
  root_tree->SetBranchAddress("min_ranges", &min_ranges_ptr);
 
  //access the largest amplitude value of each waveform
  std::vector< double > max_ranges;
  auto* max_ranges_ptr = &max_ranges;
  root_tree->SetBranchAddress("max_ranges", &max_ranges_ptr);
  std::vector<TH1D*> histograms;
      root_tree->GetEntry(0);
//declare vector to contain all waveforms, as histograms, in for this event only
      std::vector< TH1D* > TPC_histos;
      
      //How many waveforms do we have in this event?
      const Long64_t nWaveforms { static_cast<Long64_t> (rawadcs.size() ) };






DataLoader loader(m_file_name);
std::cout << "here 0" << std::endl;
WaveformManager_2 wm;
std::cout << "here 0.5" << std::endl;

// Get the TPC and Plane values

std::vector<int> mhit_tpc_vec = loader.getTPCVec();
std::vector<int> mhit_plane_vec = loader.getPlaneVec();
std::cout << "Size of mhit_tpc_vec: " << mhit_tpc_vec.size() << std::endl;
std::cout << "Size of mhit_plane_vec: " << mhit_plane_vec.size() << std::endl;

std::vector<int> tpcVec = loader.getTPCVec();
std::vector<int> planeVec = loader.getPlaneVec();



std::cout << "here 2" << std::endl;

// Ensure that both vectors have the same size
if (mhit_tpc_vec.size() != mhit_plane_vec.size()) {
    std::cerr << "Mismatch in sizes of TPC and Plane vectors" << std::endl;
    return -1;
}
std::cout << "here 3" << std::endl;
std::cout << "here 4" << std::endl;





for (size_t i = 0; i < mhit_tpc_vec.size(); ++i) {
 //   std::cout <<"event number = " << 0 << ", waveform number = " << i << ", min_range = " << (min_ranges)[i] << ", max_range = " <<  (max_ranges)[i] << ", number of time ticks =  " << rawadcs[i].size() <<std::endl;


    TH1D* histo { new TH1D("waveform_histo", "waveform_histo", 3400, 0.,  3400. ) };

    	  int nBins = rawadcs[i].size();

	  int Bin_counter =0;






//    std::cout << "here 4" << std::endl;
    
    if (min_ranges[i] == 0 & max_ranges[i] ==0) {
        continue;
    }
    for(const auto& iter : (rawadcs)[i] )
	    {
	      //Enable or disable this printout according to how much information yout want to see on-screen as the code runs:
	      //std::cout <<"event number = " << iEntry << ", waveform number = " << iWaveform << ", bin_counter = " << Bin_counter << ", adc_value = " << iter << std::endl;
	      
	      //fill bin number Bin_counter
	      histo->SetBinContent(Bin_counter, iter );
	      ++Bin_counter;
	    }  

    CleaningSlice(histo, i);
    wm.SetHistogram(histo);

    int tpc = mhit_tpc_vec[i];
  //  std::cout << "this is TPC: " << tpc << std::endl;
    int plane = mhit_plane_vec[i];
 //   std::cout << "this is Plane: " << plane << std::endl;

    int composite = tpc * 10 + plane;





        double baseline = CalculateBaseline(histo,0,100);

        switch (composite) {
            case 0: { // tpc = 0, plane = 1
                auto negativePeak1 = wm.FindMaxNegativePeak();
                dataCase1.peakHeight = negativePeak1.second;
                dataCase1.peak_time = negativePeak1.first;
                dataCase1.avgBaseline = baseline;
                dataCase1.pulse_start = wm.FindStartBin(negativePeak1.first, baseline);
                dataCase1.pulse_end = wm.FindEndBin(negativePeak1.first, baseline);
                dataCase1.waveformID = i;
                treeCase1->Fill();

                auto positivePeak1 = wm.FindMaxPositivePeak(negativePeak1.first);
                dataCase2.peakHeight = positivePeak1.second;
                dataCase2.peak_time = positivePeak1.first;
                dataCase2.avgBaseline = baseline;
                dataCase2.pulse_start = wm.FindStartBin(positivePeak1.first, baseline);
                dataCase2.pulse_end = wm.FindEndBin(positivePeak1.first, baseline);
                dataCase2.waveformID = i;
                treeCase2->Fill();
                break;
            }
            case 10:  { // tpc = 1, plane = 1
                auto negativePeak2 = wm.FindMaxNegativePeak();
                dataCase3.peakHeight = negativePeak2.second;
                dataCase3.peak_time = negativePeak2.first;
                dataCase3.avgBaseline = baseline;
                dataCase3.pulse_start = wm.FindStartBin(negativePeak2.first, baseline);
                dataCase3.pulse_end = wm.FindEndBin(negativePeak2.first, baseline);
                dataCase3.waveformID = i;
                treeCase3->Fill();

                auto positivePeak2 = wm.FindMaxPositivePeak(negativePeak2.first);
                dataCase4.peakHeight = positivePeak2.second;
                dataCase4.peak_time = positivePeak2.first;
                dataCase4.avgBaseline = baseline;
                dataCase4.pulse_start = wm.FindStartBin(positivePeak2.first, baseline);
                dataCase4.pulse_end = wm.FindEndBin(positivePeak2.first, baseline);
                dataCase4.waveformID = i;
                treeCase4->Fill();
                break;
            }
            default:{
                // Handle unexpected cases
                break;
            }
        }
    }

    fileCase1->Write();
    fileCase2->Write();
    fileCase3->Write();
    fileCase4->Write();
    
    fileCase1->Close();
    fileCase2->Close();
    fileCase3->Close();
    fileCase4->Close();

    heatmap_multi(filenames, treenames, "waveformID", "peak_time");


    return 0;
}
