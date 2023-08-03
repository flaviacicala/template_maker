#include <iostream>
#include <TFile.h>
#include <TTree.h>

void CreateTree()
{
    std::cout << "Starting program...\n";

    // Create a new file + a clone of old tree in new file
    TFile *file = new TFile("treeFile.root","RECREATE");
    TTree *tree = new TTree("T","tree with interesting variables");

    int mhit_tpc, mhit_plane;
    tree->Branch("mhit_tpc",&mhit_tpc,"mhit_tpc/I");
    tree->Branch("mhit_plane",&mhit_plane,"mhit_plane/I");

    std::cout << "File and tree created. Starting loop...\n";

    for(int i = 0; i < 11224; i++) {
        if(i >= 0 && i < 1700) {
            mhit_tpc = 0; mhit_plane = 0;
        } else if(i >= 1700 && i < 3900) {
            mhit_tpc = 0; mhit_plane = 1;
        } else if(i >= 3900 && i < 5700) {
            mhit_tpc = 0; mhit_plane = 2;
        } else if(i >= 5700 && i < 7900) {
            mhit_tpc = 1; mhit_plane = 0;
        } else if(i >= 7900 && i < 9500) {
            mhit_tpc = 1; mhit_plane = 0;
        } else if(i >= 9500 && i < 11224) {
            mhit_tpc = 1; mhit_plane = 2;
        }

        if(i % 1000 == 0) {
            std::cout << "Processing entry: " << i << "\n";
        }

        tree->Fill();
    }

    std::cout << "Loop ended. Writing tree...\n";

    tree->Print();
    tree->Write();

    std::cout << "Tree written to file. Closing file...\n";

    delete file;

    std::cout << "File closed. Program finished.\n";
}

int main() 
{
    CreateTree();
    return 0;
}
