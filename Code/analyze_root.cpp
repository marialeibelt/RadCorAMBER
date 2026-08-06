#include <iostream>
#include <vector>

#include "TFile.h"
#include "TTree.h"

using namespace std;

int main()
{
    string inputfile = "05_08_evtgen/out/05_08_evtgen.root";

    TFile *f = TFile::Open(inputfile.c_str());

    TTree *tree = (TTree*)f->Get("Output");


    vector<int>* scatteredPID = nullptr;
    vector<double>* scatteredEnergy = nullptr;
    vector<double>* scatteredMomentumX = nullptr;
    vector<double>* scatteredMomentumY = nullptr;
    vector<double>* scatteredMomentumZ = nullptr;


    tree->SetBranchAddress(
        "scatteredPID",
        &scatteredPID
    );

    tree->SetBranchAddress(
        "scatteredEnergy",
        &scatteredEnergy
    );

    tree->SetBranchAddress(
        "scatteredMomentumX",
        &scatteredMomentumX
    );

    tree->SetBranchAddress(
        "scatteredMomentumY",
        &scatteredMomentumY
    );

    tree->SetBranchAddress(
        "scatteredMomentumZ",
        &scatteredMomentumZ
    );


    tree->GetEntry(0);


    cout << "Number of scattered particles: "
         << scatteredPID->size()
         << endl;


    for(size_t i=0; i<scatteredPID->size(); i++)
    {
        cout << "Particle " << i << endl;

        cout << " PID = "
             << scatteredPID->at(i)
             << endl;

        cout << " E = "
             << scatteredEnergy->at(i)
             << endl;

        cout << " px = "
             << scatteredMomentumX->at(i)
             << endl;

        cout << " py = "
             << scatteredMomentumY->at(i)
             << endl;

        cout << " pz = "
             << scatteredMomentumZ->at(i)
             << endl;
    }


    return 0;
}