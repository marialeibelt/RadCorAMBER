#include <iostream>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"

using namespace std;

int main()
{
    string inputfile = "05_08_evtgen/out/05_08_evtgen.root";

    cout << "Opening ROOT file: " << inputfile << endl;

    TFile *f = TFile::Open(inputfile.c_str());

    if (!f || f->IsZombie()) {
        cerr << "Error: could not open ROOT file!" << endl;
        return 1;
    }

    TTree *tree = (TTree*)f->Get("Output");

    if (!tree) {
        cerr << "Error: TTree 'Output' not found!" << endl;
        return 1;
    }


    vector<int>* scatteredPID = nullptr;
    vector<double>* scatteredEnergy = nullptr;
    vector<double>* scatteredMomentumX = nullptr;
    vector<double>* scatteredMomentumY = nullptr;
    vector<double>* scatteredMomentumZ = nullptr;

    double weight;


    tree->SetBranchAddress("scatteredPID", &scatteredPID);
    tree->SetBranchAddress("scatteredEnergy", &scatteredEnergy);
    tree->SetBranchAddress("scatteredMomentumX", &scatteredMomentumX);
    tree->SetBranchAddress("scatteredMomentumY", &scatteredMomentumY);
    tree->SetBranchAddress("scatteredMomentumZ", &scatteredMomentumZ);
    tree->SetBranchAddress("weight", &weight);


    // Plot weight distribution
    TCanvas *c = new TCanvas("c", "Weight distribution", 800, 600);

    tree->Draw("weight");

    string outputpdf = "weights.pdf";
    c->SaveAs(outputpdf.c_str());

    cout << "Weight histogram created successfully." << endl;
    cout << "Saved as: " << outputpdf << endl;


    return 0;
}