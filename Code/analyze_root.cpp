#include <iostream>
#include <vector>
#include <cmath>
#include <filesystem>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLegend.h"

using namespace std;


// ============================================================
// Function for plotting a single histogram
// ============================================================

void saveHistogram(
    TH1D *hist,
    const string &xTitle,
    const string &yTitle,
    const string &outputFile,
    bool logY = true)
{
    TCanvas canvas("canvas", hist->GetTitle(), 800, 600);

    if (logY)
        canvas.SetLogy();

    hist->GetXaxis()->SetTitle(xTitle.c_str());
    hist->GetYaxis()->SetTitle(yTitle.c_str());

    hist->Draw("HIST");

    canvas.SaveAs(outputFile.c_str());

    cout << "Histogram created: " << outputFile << endl;
}


// ============================================================
// Function for comparing LO and NLO
// Each histogram is normalized individually to integral = 1
// ============================================================

void saveComparison(
    TH1D *hLO,
    TH1D *hNLO,
    const string &xTitle,
    const string &outputFile,
    bool logY = true)
{
    // Normalize each histogram separately
    if (hLO->Integral() > 0)
        hLO->Scale(1.0 / hLO->Integral());

    if (hNLO->Integral() > 0)
        hNLO->Scale(1.0 / hNLO->Integral());

    hLO->SetFillColor(kGreen+2);
    hLO->SetLineColor(kGreen+2);
    hLO->SetFillStyle(3004);

    hNLO->SetLineColor(kPink+10);
    hNLO->SetLineWidth(2);

    hLO->GetXaxis()->SetTitle(xTitle.c_str());
    hLO->GetYaxis()->SetTitle("Normalized entries");

    TCanvas canvas("canvas", "LO vs NLO", 800, 600);

    if (logY)
        canvas.SetLogy();

    hLO->Draw("HIST");
    hNLO->Draw("HIST SAME");

    TLegend legend(0.70, 0.75, 0.88, 0.88);
    legend.AddEntry(hLO, "LO", "f");
    legend.AddEntry(hNLO, "NLO", "l");
    legend.Draw();

    canvas.SaveAs(outputFile.c_str());

    cout << "Comparison created: " << outputFile << endl;
}


// ============================================================
// Main
// ============================================================

int main(int argc, char* argv[])
{
    // ========================================================
    // Command-line arguments
    // ========================================================

    if (argc != 4) {
        cerr << "Usage: " << argv[0]
             << " <LO_input.root> <NLO_input.root> <output_folder>"
             << endl;
        return 1;
    }

    string loInputfile  = argv[1];
    string nloInputfile = argv[2];
    string outputFolder = argv[3];


    // ========================================================
    // Output folders
    // ========================================================

    string loFolder   = outputFolder + "/LO";
    string nloFolder  = outputFolder + "/NLO";
    string compFolder = outputFolder + "/comparison";

    filesystem::create_directories(loFolder);
    filesystem::create_directories(nloFolder);
    filesystem::create_directories(compFolder);

    cout << "LO input file:  " << loInputfile << endl;
    cout << "NLO input file: " << nloInputfile << endl;
    cout << "Output folder:  " << outputFolder << endl;


    // ========================================================
    // Open ROOT files
    // ========================================================

    cout << "Opening LO ROOT file: " << loInputfile << endl;

    TFile *fLO = TFile::Open(loInputfile.c_str());

    if (!fLO || fLO->IsZombie()) {
        cerr << "Error: could not open LO ROOT file!" << endl;
        return 1;
    }


    cout << "Opening NLO ROOT file: " << nloInputfile << endl;

    TFile *fNLO = TFile::Open(nloInputfile.c_str());

    if (!fNLO || fNLO->IsZombie()) {
        cerr << "Error: could not open NLO ROOT file!" << endl;
        fLO->Close();
        return 1;
    }

    // ========================================================
    // Trees
    // ========================================================

    TTree *treeLO = (TTree*)fLO->Get("Output");
    TTree *treeNLO = (TTree*)fNLO->Get("Output");

    if (!treeLO) {
        cerr << "Error: TTree 'Output' not found in LO file!" << endl;
        return 1;
    }

    if (!treeNLO) {
        cerr << "Error: TTree 'Output' not found in NLO file!" << endl;
        return 1;
    }


    Long64_t nEventsLO = treeLO->GetEntries();
    Long64_t nEventsNLO = treeNLO->GetEntries();

    cout << "Number of LO events:  " << nEventsLO << endl;
    cout << "Number of NLO events: " << nEventsNLO << endl;


    // ========================================================
    // Branch pointers
    // ========================================================

    vector<int>* scatteredPIDLO = nullptr;

    vector<double>* scatteredEnergyLO = nullptr;
    vector<double>* scatteredMomentumXLO = nullptr;
    vector<double>* scatteredMomentumYLO = nullptr;
    vector<double>* scatteredMomentumZLO = nullptr;

    double weightLO;


    vector<int>* scatteredPIDNLO = nullptr;

    vector<double>* scatteredEnergyNLO = nullptr;
    vector<double>* scatteredMomentumXNLO = nullptr;
    vector<double>* scatteredMomentumYNLO = nullptr;
    vector<double>* scatteredMomentumZNLO = nullptr;

    double weightNLO;


    // ========================================================
    // Set branch addresses - LO
    // ========================================================

    treeLO->SetBranchAddress("scatteredPID", &scatteredPIDLO);
    treeLO->SetBranchAddress("scatteredEnergy", &scatteredEnergyLO);
    treeLO->SetBranchAddress("scatteredMomentumX", &scatteredMomentumXLO);
    treeLO->SetBranchAddress("scatteredMomentumY", &scatteredMomentumYLO);
    treeLO->SetBranchAddress("scatteredMomentumZ", &scatteredMomentumZLO);
    treeLO->SetBranchAddress("weight", &weightLO);


    // ========================================================
    // Set branch addresses - NLO
    // ========================================================

    treeNLO->SetBranchAddress("scatteredPID", &scatteredPIDNLO);
    treeNLO->SetBranchAddress("scatteredEnergy", &scatteredEnergyNLO);
    treeNLO->SetBranchAddress("scatteredMomentumX", &scatteredMomentumXNLO);
    treeNLO->SetBranchAddress("scatteredMomentumY", &scatteredMomentumYNLO);
    treeNLO->SetBranchAddress("scatteredMomentumZ", &scatteredMomentumZNLO);
    treeNLO->SetBranchAddress("weight", &weightNLO);


    // ========================================================
    // Constants
    // ========================================================

    const double beamEnergyLab = 100.0;
    const double protonMass = 0.938272088;
    const double muonMass = 0.105658375;

    const double EbeamCMS = 6.833951293;
    const double pBeamCMS =
        sqrt(EbeamCMS * EbeamCMS - muonMass * muonMass);

    const double initialMuonEnergy = 6.833951293;
    const double initialProtonEnergy = 6.897251706;
    const double initialTotalEnergy =
        initialMuonEnergy + initialProtonEnergy;

    // CMS -> LAB boost
    double beta =
        beamEnergyLab / (beamEnergyLab + protonMass);

    double gamma =
        1.0 / sqrt(1.0 - beta * beta);


    // ========================================================
    // 1. Weight distributions
    // ========================================================

    TH1D *hWeightLO =
        new TH1D("hWeightLO",
                 "LO Weight distribution",
                 500, 0, 1.1);

    TH1D *hWeightNLO =
        new TH1D("hWeightNLO",
                 "NLO Weight distribution",
                 500, 0, 1.1);


    for (Long64_t i = 0; i < nEventsLO; i++) {

        treeLO->GetEntry(i);

        hWeightLO->Fill(weightLO);
    }


    for (Long64_t i = 0; i < nEventsNLO; i++) {

        treeNLO->GetEntry(i);

        hWeightNLO->Fill(weightNLO);
    }


    saveHistogram(
        hWeightLO,
        "Weight",
        "Entries",
        loFolder + "/weights.pdf");

    saveHistogram(
        hWeightNLO,
        "Weight",
        "Entries",
        nloFolder + "/weights.pdf");


    // ========================================================
    // 2. Scattered muon energy
    // ========================================================

    TH1D *hMuonEnergyLO =
        new TH1D("hMuonEnergyLO",
                 "LO Scattered muon energy",
                 500, 0, 10);

    TH1D *hMuonEnergyNLO =
        new TH1D("hMuonEnergyNLO",
                 "NLO Scattered muon energy",
                 500, 0, 10);


    for (Long64_t i = 0; i < nEventsLO; i++) {

        treeLO->GetEntry(i);

        for (size_t j = 0; j < scatteredPIDLO->size(); j++) {

            if (scatteredPIDLO->at(j) != 13)
                continue;

            hMuonEnergyLO->Fill(
                scatteredEnergyLO->at(j));
        }
    }


    for (Long64_t i = 0; i < nEventsNLO; i++) {

        treeNLO->GetEntry(i);

        for (size_t j = 0; j < scatteredPIDNLO->size(); j++) {

            if (scatteredPIDNLO->at(j) != 13)
                continue;

            hMuonEnergyNLO->Fill(
                scatteredEnergyNLO->at(j));
        }
    }


    saveHistogram(
        hMuonEnergyLO,
        "Scattered muon energy [GeV]",
        "Entries",
        loFolder + "/scattered_muon_energy.pdf");

    saveHistogram(
        hMuonEnergyNLO,
        "Scattered muon energy [GeV]",
        "Entries",
        nloFolder + "/scattered_muon_energy.pdf");

    saveComparison(
        hMuonEnergyLO,
        hMuonEnergyNLO,
        "Scattered muon energy [GeV]",
        compFolder + "/scattered_muon_energy.pdf");


    // ========================================================
    // 3. Photon energy
    // ========================================================
    TH1D *hPhotonEnergyNLO =
        new TH1D("hPhotonEnergyNLO",
                 "NLO Photon energy",
                 500, -1, 10);

    for (Long64_t i = 0; i < nEventsNLO; i++) {

        treeNLO->GetEntry(i);

        for (size_t j = 0; j < scatteredPIDNLO->size(); j++) {

            if (scatteredPIDNLO->at(j) != 22)
                continue;

            hPhotonEnergyNLO->Fill(
                scatteredEnergyNLO->at(j));
        }
    }

    saveHistogram(
        hPhotonEnergyNLO,
        "Photon energy [GeV]",
        "Entries",
        nloFolder + "/photon_energy.pdf");


    // ========================================================
    // 4. Photon angle - CMS
    // ========================================================
    TH1D *hPhotonAngleNLO =
        new TH1D("hPhotonAngleNLO",
                 "NLO Photon angle - CMS",
                 500, 0, 20);


    for (Long64_t i = 0; i < nEventsNLO; i++) {

        treeNLO->GetEntry(i);

        for (size_t j = 0; j < scatteredPIDNLO->size(); j++) {

            if (scatteredPIDNLO->at(j) != 22)
                continue;

            double px = scatteredMomentumXNLO->at(j);
            double py = scatteredMomentumYNLO->at(j);
            double pz = scatteredMomentumZNLO->at(j);

            double pPerp = sqrt(px * px + py * py);

            double theta = atan2(pPerp, pz);

            theta *= 1000.0;

            hPhotonAngleNLO->Fill(theta);
        }
    }

    saveHistogram(
        hPhotonAngleNLO,
        "Photon angle #theta_{#gamma}^{CMS} [mrad]",
        "Entries",
        nloFolder + "/photon_angle_cms.pdf");


    // ========================================================
    // 5. Scattered muon angle - CMS
    // ========================================================

    TH1D *hMuonAngleCMSLO =
        new TH1D("hMuonAngleCMSLO",
                 "LO Scattered muon angle - CMS",
                 500, 3, 31);

    TH1D *hMuonAngleCMSNLO =
        new TH1D("hMuonAngleCMSNLO",
                 "NLO Scattered muon angle - CMS",
                 500, 3, 31);


    for (Long64_t i = 0; i < nEventsLO; i++) {

        treeLO->GetEntry(i);

        for (size_t j = 0; j < scatteredPIDLO->size(); j++) {

            if (scatteredPIDLO->at(j) != 13)
                continue;

            double px = scatteredMomentumXLO->at(j);
            double py = scatteredMomentumYLO->at(j);
            double pz = scatteredMomentumZLO->at(j);

            double pPerp = sqrt(px * px + py * py);

            double theta = atan2(pPerp, pz);

            theta *= 1000.0;

            hMuonAngleCMSLO->Fill(theta);
        }
    }


    for (Long64_t i = 0; i < nEventsNLO; i++) {

        treeNLO->GetEntry(i);

        for (size_t j = 0; j < scatteredPIDNLO->size(); j++) {

            if (scatteredPIDNLO->at(j) != 13)
                continue;

            double px = scatteredMomentumXNLO->at(j);
            double py = scatteredMomentumYNLO->at(j);
            double pz = scatteredMomentumZNLO->at(j);

            double pPerp = sqrt(px * px + py * py);

            double theta = atan2(pPerp, pz);

            theta *= 1000.0;

            hMuonAngleCMSNLO->Fill(theta);
        }
    }


    saveHistogram(
        hMuonAngleCMSLO,
        "Scattered muon angle #theta_{#mu}^{CMS} [mrad]",
        "Entries",
        loFolder + "/scattered_muon_angle_cms.pdf");

    saveHistogram(
        hMuonAngleCMSNLO,
        "Scattered muon angle #theta_{#mu}^{CMS} [mrad]",
        "Entries",
        nloFolder + "/scattered_muon_angle_cms.pdf");

    saveComparison(
        hMuonAngleCMSLO,
        hMuonAngleCMSNLO,
        "Scattered muon angle #theta_{#mu}^{CMS} [mrad]",
        compFolder + "/scattered_muon_angle_cms.pdf");


    // ========================================================
    // 6. Scattered muon angle - LAB
    // ========================================================

    TH1D *hMuonAngleLABLO =
        new TH1D("hMuonAngleLABLO",
                 "LO Scattered muon angle - LAB",
                 500, 0, 2.5);

    TH1D *hMuonAngleLABNLO =
        new TH1D("hMuonAngleLABNLO",
                 "NLO Scattered muon angle - LAB",
                 500, 0, 2.5);


    for (Long64_t i = 0; i < nEventsLO; i++) {

        treeLO->GetEntry(i);

        for (size_t j = 0; j < scatteredPIDLO->size(); j++) {

            if (scatteredPIDLO->at(j) != 13)
                continue;

            double px = scatteredMomentumXLO->at(j);
            double py = scatteredMomentumYLO->at(j);
            double pz = scatteredMomentumZLO->at(j);

            double E = scatteredEnergyLO->at(j);

            double pzLAB =
                gamma * (pz + beta * E);

            double pPerp =
                sqrt(px * px + py * py);

            double thetaLAB =
                atan2(pPerp, pzLAB);

            thetaLAB *= 1000.0;

            hMuonAngleLABLO->Fill(thetaLAB);
        }
    }


    for (Long64_t i = 0; i < nEventsNLO; i++) {

        treeNLO->GetEntry(i);

        for (size_t j = 0; j < scatteredPIDNLO->size(); j++) {

            if (scatteredPIDNLO->at(j) != 13)
                continue;

            double px = scatteredMomentumXNLO->at(j);
            double py = scatteredMomentumYNLO->at(j);
            double pz = scatteredMomentumZNLO->at(j);

            double E = scatteredEnergyNLO->at(j);

            double pzLAB =
                gamma * (pz + beta * E);

            double pPerp =
                sqrt(px * px + py * py);

            double thetaLAB =
                atan2(pPerp, pzLAB);

            thetaLAB *= 1000.0;

            hMuonAngleLABNLO->Fill(thetaLAB);
        }
    }


    saveHistogram(
        hMuonAngleLABLO,
        "Scattered muon angle #theta_{#mu}^{LAB} [mrad]",
        "Entries",
        loFolder + "/scattered_muon_angle_lab.pdf");

    saveHistogram(
        hMuonAngleLABNLO,
        "Scattered muon angle #theta_{#mu}^{LAB} [mrad]",
        "Entries",
        nloFolder + "/scattered_muon_angle_lab.pdf");

    saveComparison(
        hMuonAngleLABLO,
        hMuonAngleLABNLO,
        "Scattered muon angle #theta_{#mu}^{LAB} [mrad]",
        compFolder + "/scattered_muon_angle_lab.pdf");


    // ========================================================
    // 7. Q^2
    // ========================================================

    TH1D *hQ2LO =
        new TH1D("hQ2LO",
                 "LO Momentum transfer Q^{2}",
                 500, 0, 0.05);

    TH1D *hQ2NLO =
        new TH1D("hQ2NLO",
                 "NLO Momentum transfer Q^{2}",
                 500, 0, 0.05);


    for (Long64_t i = 0; i < nEventsLO; i++) {

        treeLO->GetEntry(i);

        for (size_t j = 0; j < scatteredPIDLO->size(); j++) {

            if (scatteredPIDLO->at(j) != 13)
                continue;

            double px = scatteredMomentumXLO->at(j);
            double py = scatteredMomentumYLO->at(j);
            double pz = scatteredMomentumZLO->at(j);

            double E = scatteredEnergyLO->at(j);

            double dE = EbeamCMS - E;
            double dpx = -px;
            double dpy = -py;
            double dpz = pBeamCMS - pz;

            double q2 =
                dE * dE
                - dpx * dpx
                - dpy * dpy
                - dpz * dpz;

            q2 = -q2;

            hQ2LO->Fill(q2);
        }
    }


    for (Long64_t i = 0; i < nEventsNLO; i++) {

        treeNLO->GetEntry(i);

        for (size_t j = 0; j < scatteredPIDNLO->size(); j++) {

            if (scatteredPIDNLO->at(j) != 13)
                continue;

            double px = scatteredMomentumXNLO->at(j);
            double py = scatteredMomentumYNLO->at(j);
            double pz = scatteredMomentumZNLO->at(j);

            double E = scatteredEnergyNLO->at(j);

            double dE = EbeamCMS - E;
            double dpx = -px;
            double dpy = -py;
            double dpz = pBeamCMS - pz;

            double q2 =
                dE * dE
                - dpx * dpx
                - dpy * dpy
                - dpz * dpz;

            q2 = -q2;

            hQ2NLO->Fill(q2);
        }
    }


    saveHistogram(
        hQ2LO,
        "Q^{2} [GeV^{2}]",
        "Entries",
        loFolder + "/Q2.pdf");

    saveHistogram(
        hQ2NLO,
        "Q^{2} [GeV^{2}]",
        "Entries",
        nloFolder + "/Q2.pdf");

    saveComparison(
        hQ2LO,
        hQ2NLO,
        "Q^{2} [GeV^{2}]",
        compFolder + "/Q2.pdf");


    // ========================================================
    // 8. Energy conservation
    // ========================================================

    TH1D *hEnergyConservationLO =
        new TH1D("hEnergyConservationLO",
                 "LO Energy conservation",
                 500, -1e-6, 1e-6);

    TH1D *hEnergyConservationNLO =
        new TH1D("hEnergyConservationNLO",
                 "NLO Energy conservation",
                 500, -1e-6, 1e-6);


    for (Long64_t i = 0; i < nEventsLO; i++) {

        treeLO->GetEntry(i);

        double finalTotalEnergy = 0.0;

        for (size_t j = 0; j < scatteredEnergyLO->size(); j++)
            finalTotalEnergy += scatteredEnergyLO->at(j);

        double deltaE =
            finalTotalEnergy - initialTotalEnergy;

        hEnergyConservationLO->Fill(deltaE);
    }


    for (Long64_t i = 0; i < nEventsNLO; i++) {

        treeNLO->GetEntry(i);

        double finalTotalEnergy = 0.0;

        for (size_t j = 0; j < scatteredEnergyNLO->size(); j++)
            finalTotalEnergy += scatteredEnergyNLO->at(j);

        double deltaE =
            finalTotalEnergy - initialTotalEnergy;

        hEnergyConservationNLO->Fill(deltaE);
    }


    saveHistogram(
        hEnergyConservationLO,
        "#Delta E = E_{final} - E_{initial} [GeV]",
        "Entries",
        loFolder + "/energy_conservation.pdf");

    saveHistogram(
        hEnergyConservationNLO,
        "#Delta E = E_{final} - E_{initial} [GeV]",
        "Entries",
        nloFolder + "/energy_conservation.pdf");

    saveComparison(
        hEnergyConservationLO,
        hEnergyConservationNLO,
        "#Delta E = E_{final} - E_{initial} [GeV]",
        compFolder + "/energy_conservation.pdf",
        false);


    // ========================================================
    // Done
    // ========================================================

    cout << endl;
    cout << "All histograms created successfully." << endl;

    cout << "LO events:  " << nEventsLO << endl;
    cout << "NLO events: " << nEventsNLO << endl;

    fLO->Close();
    fNLO->Close();

    return 0;
}