#include <iostream>
#include <vector>
#include <cmath>
#include <filesystem>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1D.h"

using namespace std;


// ============================================================
// Function for plotting and saving a histogram
// ============================================================

void saveHistogram(TH1D *hist,const string &xTitle,const string &yTitle,const string &outputFile,bool logY = true)
{
    TCanvas canvas("canvas",hist->GetTitle(),800,600);

    if (logY)
        canvas.SetLogy();

    hist->GetXaxis()->SetTitle(xTitle.c_str());
    hist->GetYaxis()->SetTitle(yTitle.c_str());

    hist->Draw("HIST");

    canvas.SaveAs(outputFile.c_str());

    cout << "Histogram created: "<< outputFile << endl;
}


// ============================================================
// Main
// ============================================================

int main()
{
    string inputfile ="05_08_evtgen/out/05_08_evtgen.root";
    string outputFolder = "analysis_output_05_08_unweighted";
    filesystem::create_directories(outputFolder);
    cout << "Output folder: "<< outputFolder << endl;

    cout << "Opening ROOT file: "<< inputfile << endl;

    TFile *f = TFile::Open(inputfile.c_str());

    if (!f || f->IsZombie()) {
        cerr << "Error: could not open ROOT file!" << endl;
        return 1;
    }

    TTree *tree =(TTree*)f->Get("Output");

    if (!tree) {
        cerr << "Error: TTree 'Output' not found!"<< endl;
        return 1;
    }


    // ========================================================
    // Branches
    // ========================================================

    vector<int>* scatteredPID = nullptr;

    vector<double>* scatteredEnergy = nullptr;

    vector<double>* scatteredMomentumX = nullptr;
    vector<double>* scatteredMomentumY = nullptr;
    vector<double>* scatteredMomentumZ = nullptr;

    double weight;


    tree->SetBranchAddress("scatteredPID",&scatteredPID);
    tree->SetBranchAddress("scatteredEnergy",&scatteredEnergy);
    tree->SetBranchAddress("scatteredMomentumX",&scatteredMomentumX);
    tree->SetBranchAddress("scatteredMomentumY",&scatteredMomentumY);
    tree->SetBranchAddress("scatteredMomentumZ",&scatteredMomentumZ);
    tree->SetBranchAddress("weight",&weight);

    Long64_t nEvents = tree->GetEntries();

    cout << "Number of events: "<< nEvents << endl;


    // ========================================================
    // 1. Weight distribution
    // ========================================================

    TH1D *hWeight =new TH1D("hWeight","Weight distribution",500,0,1.1);

    for (Long64_t i = 0; i < nEvents; i++) {
        tree->GetEntry(i);
        hWeight->Fill(weight);
    }

    saveHistogram(hWeight,"Weight","Entries",outputFolder + "/weights.pdf");


    // ========================================================
    // 2. Scattered muon energy
    // ========================================================

    TH1D *hMuonEnergy = new TH1D("hMuonEnergy","Scattered muon energy",500,0,10);

    for (Long64_t i = 0; i < nEvents; i++) {
        tree->GetEntry(i);
        for (size_t j = 0; j < scatteredPID->size(); j++) {
            if (scatteredPID->at(j) != 13)
                continue;
            hMuonEnergy->Fill(scatteredEnergy->at(j));
        }
    }

    saveHistogram(hMuonEnergy,"Scattered muon energy [GeV]","Entries",outputFolder + "/scattered_muon_energy.pdf");


    // ========================================================
    // 3. Photon energy
    // ========================================================

    TH1D *hPhotonEnergy = new TH1D("hPhotonEnergy","Photon energy",500,0,10);

    for (Long64_t i = 0; i < nEvents; i++) {
        tree->GetEntry(i);
        for (size_t j = 0; j < scatteredPID->size(); j++) {

            // Currently photons have PID 22
            if (scatteredPID->at(j) != 22)
                continue;

            hPhotonEnergy->Fill(scatteredEnergy->at(j));
        }
    }

    saveHistogram(hPhotonEnergy,"Photon energy [GeV]","Entries",outputFolder + "/photon_energy.pdf");


    // ========================================================
    // 4. Photon angle
    // ========================================================

    TH1D *hPhotonAngle = new TH1D("hPhotonAngle","Photon angle - CMS",500,0,20);

    for (Long64_t i = 0; i < nEvents; i++) {
        tree->GetEntry(i);
        for (size_t j = 0;j < scatteredPID->size();j++) {

            if (scatteredPID->at(j) != 22)
                continue;

            double px = scatteredMomentumX->at(j);
            double py = scatteredMomentumY->at(j);
            double pz = scatteredMomentumZ->at(j);

            double pPerp =sqrt(px * px +py * py);

            double theta = atan2(pPerp,pz);

            // rad -> mrad
            theta *= 1000.0;

            hPhotonAngle->Fill(theta);
        }
    }


    saveHistogram(hPhotonAngle,"Photon angle #theta_{#gamma}^{CMS} [mrad]","Entries",outputFolder + "/photon_angle_cms.pdf");


    // ========================================================
    // 5. Scattered muon angle - CMS
    // ========================================================

    TH1D *hMuonAngleCMS =new TH1D("hMuonAngleCMS","Scattered muon angle - CMS",500,0,5);

    for (Long64_t i = 0; i < nEvents; i++) {
        tree->GetEntry(i);
        for (size_t j = 0;j < scatteredPID->size();j++) {

            if (scatteredPID->at(j) != 13)
                continue;

            double px = scatteredMomentumX->at(j);
            double py = scatteredMomentumY->at(j);
            double pz = scatteredMomentumZ->at(j);

            double pPerp = sqrt(px * px + py * py);
            double theta = atan2(pPerp,pz);

            // rad -> mrad
            theta *= 1000.0;

            hMuonAngleCMS->Fill(theta);
        }
    }

    saveHistogram(hMuonAngleCMS,"Scattered muon angle #theta_{#mu}^{CMS} [mrad]","Entries",outputFolder + "/scattered_muon_angle_cms.pdf");


    // ========================================================
    // 6. Scattered muon angle - LAB
    // ========================================================

    TH1D *hMuonAngleLAB = new TH1D("hMuonAngleLAB","Scattered muon angle - LAB",500,0,5);

    // Beam energy and proton mass
    const double beamEnergyLab = 100.0;       // GeV
    const double protonMass = 0.938272088;    // GeV

    // CMS -> LAB boost
    double beta =beamEnergyLab / (beamEnergyLab + protonMass);
    double gamma =1.0 / sqrt(1.0 - beta * beta);

    for (Long64_t i = 0; i < nEvents; i++) {
        tree->GetEntry(i);
        for (size_t j = 0;j < scatteredPID->size();j++) {

            if (scatteredPID->at(j) != 13)
                continue;

            double px =scatteredMomentumX->at(j);
            double py =scatteredMomentumY->at(j);
            double pz =scatteredMomentumZ->at(j);

            double E = scatteredEnergy->at(j);
            double pzLAB = gamma * (pz + beta * E);
            double pPerp = sqrt(px * px + py * py);

            double thetaLAB = atan2(pPerp,pzLAB);

            // rad -> mrad
            thetaLAB *= 1000.0;

            hMuonAngleLAB->Fill(thetaLAB);
        }
    }

    saveHistogram(hMuonAngleLAB,"Scattered muon angle #theta_{#mu}^{LAB} [mrad]","Entries",outputFolder + "/scattered_muon_angle_lab.pdf");


    // ========================================================
    // 7. Q^2
    // ========================================================

    TH1D *hQ2 = new TH1D("hQ2","Momentum transfer Q^{2}",500,0,0.05);

    const double muonMass = 0.105658375; // GeV

    // Initial muon in CMS
    // E = 6.833951 GeV
    const double EbeamCMS = 6.833951293;
    const double pBeamCMS = sqrt(EbeamCMS * EbeamCMS - muonMass * muonMass);

    for (Long64_t i = 0; i < nEvents; i++) {
        tree->GetEntry(i);
        for (size_t j = 0;j < scatteredPID->size();j++) {

            if (scatteredPID->at(j) != 13)
                continue;

            double px = scatteredMomentumX->at(j);
            double py = scatteredMomentumY->at(j);
            double pz = scatteredMomentumZ->at(j);

            double E = scatteredEnergy->at(j);

            // Initial and final muon four-vectors
            double dE = EbeamCMS - E;
            double dpx = -px;
            double dpy = -py;
            double dpz = pBeamCMS - pz;

            double q2 = dE * dE - dpx * dpx - dpy * dpy - dpz * dpz;

            // Q^2 = -q^2
            q2 = -q2;

            hQ2->Fill(q2);
        }
    }

    saveHistogram(hQ2,"Q^{2} [GeV^{2}]","Entries",outputFolder + "/Q2.pdf");


    // ========================================================
    // 8. Energy conservation
    // ========================================================

    TH1D *hEnergyConservation = new TH1D("hEnergyConservation","Energy conservation",500,-1e-6,1e-6);

    // Initial-state energies in CMS
    const double initialMuonEnergy = 6.833951293;
    const double initialProtonEnergy = 6.897251706; 
    const double initialTotalEnergy = initialMuonEnergy + initialProtonEnergy;
 
    for (Long64_t i = 0; i < nEvents; i++) {
        tree->GetEntry(i);

        double finalTotalEnergy = 0.0;

        for (size_t j = 0;j < scatteredEnergy->size();j++) {
            finalTotalEnergy += scatteredEnergy->at(j);
        }

        double deltaE = finalTotalEnergy - initialTotalEnergy;

        hEnergyConservation->Fill(deltaE);
    }

    saveHistogram(hEnergyConservation,"#Delta E = E_{final} - E_{initial} [GeV]","Entries",outputFolder + "/energy_conservation.pdf");


    // ========================================================
    // Done
    // ========================================================

    cout << endl;
    cout << "All histograms created successfully." << endl;

    f->Close();

    return 0;
}