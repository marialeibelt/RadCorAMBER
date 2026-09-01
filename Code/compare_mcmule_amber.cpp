#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <filesystem>
#include <memory>
#include <algorithm>
#include <utility>

#include <zlib.h>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"
#include "TLorentzVector.h"

using namespace std;


// ============================================================
// Configuration
// ============================================================

// In the TGEANT particle block observed so far:
//   incoming beam muon: GeV
//   final-state particles: MeV
//
// Therefore final-state particle four-vectors are converted
// with this factor before comparison to the McMule ROOT file.
constexpr double TGEANT_FINALSTATE_TO_GEV = 1.0e-3;

// Plot ranges
constexpr double Q2_MIN = 1.0e-5;
constexpr double Q2_MAX = 1.0e-1;

constexpr double EMU_MIN = 0.0;
constexpr double EMU_MAX = 105.0;

constexpr double PHOTON_E_MIN = 0.0;
constexpr double PHOTON_E_MAX = 100.0;

constexpr double PHOTON_THETA_MIN = 0.0;
constexpr double PHOTON_THETA_MAX = 15.0; // mrad


// ============================================================
// TGEANT particle representation
// ============================================================

struct Particle {
    int index = 0;
    int pid   = 0;

    double px     = 0.0;
    double py     = 0.0;
    double pz     = 0.0;
    double energy = 0.0;
    double mass   = 0.0;
};


// ============================================================
// Small helpers
// ============================================================

bool readGzLine(gzFile file, string& line)
{
    const int BUFFER_SIZE = 10000;
    char buffer[BUFFER_SIZE];

    if (!gzgets(file, buffer, BUFFER_SIZE))
        return false;

    line = buffer;

    while (!line.empty() &&
           (line.back() == '\n' || line.back() == '\r')) {
        line.pop_back();
    }

    return true;
}


string getFirstToken(const string& line)
{
    stringstream ss(line);
    string token;
    ss >> token;
    return token;
}


bool isEventStart(const string& line)
{
    return getFirstToken(line) == "#EVENT";
}


bool isEventEnd(const string& line)
{
    return getFirstToken(line) == "#EVENT_END";
}


bool getEventNumber(const string& line, long long& eventNumber)
{
    stringstream ss(line);
    string tag;

    if (!(ss >> tag >> eventNumber))
        return false;

    return tag == "#EVENT";
}


bool parseSingleInteger(const string& line, int& value)
{
    stringstream ss(line);

    if (!(ss >> value))
        return false;

    string extra;

    if (ss >> extra)
        return false;

    return true;
}


bool parseParticle(const string& line, Particle& p)
{
    stringstream ss(line);

    double dummy1 = 0.0;
    double dummy2 = 0.0;
    double dummy3 = 0.0;

    if (!(ss
          >> p.index
          >> p.pid
          >> dummy1
          >> dummy2
          >> dummy3
          >> p.px
          >> p.py
          >> p.pz
          >> p.energy
          >> p.mass)) {
        return false;
    }

    return true;
}


bool findParticleBlock(
    const vector<string>& eventLines,
    vector<Particle>& particles)
{
    particles.clear();

    for (size_t i = 0; i < eventLines.size(); ++i) {

        int nParticles = 0;

        if (!parseSingleInteger(eventLines[i], nParticles))
            continue;

        if (nParticles <= 0 || nParticles > 100)
            continue;

        if (i + static_cast<size_t>(nParticles) >= eventLines.size())
            continue;

        vector<Particle> candidate;
        candidate.reserve(nParticles);

        bool valid = true;

        for (int j = 0; j < nParticles; ++j) {

            Particle p;

            if (!parseParticle(eventLines[i + 1 + j], p)) {
                valid = false;
                break;
            }

            candidate.push_back(p);
        }

        if (valid) {
            particles = std::move(candidate);
            return true;
        }
    }

    return false;
}


// Reads exactly one real TGEANT event.
// Returns false only when no further #EVENT block exists.
bool readNextTgeantEvent(
    gzFile file,
    long long& eventNumber,
    vector<Particle>& particles,
    bool& hasParticleBlock)
{
    string line;

    while (readGzLine(file, line)) {

        if (!isEventStart(line))
            continue;

        eventNumber = -1;
        getEventNumber(line, eventNumber);

        vector<string> eventLines;
        bool foundEnd = false;

        while (readGzLine(file, line)) {

            if (isEventEnd(line)) {
                foundEnd = true;
                break;
            }

            eventLines.push_back(line);
        }

        if (!foundEnd) {
            cerr << "WARNING: no #EVENT_END found for TGEANT event "
                 << eventNumber << endl;
        }

        hasParticleBlock =
            findParticleBlock(eventLines, particles);

        return true;
    }

    return false;
}


// ============================================================
// Four-vector helpers
// ============================================================

TLorentzVector makeFourVector(
    double px,
    double py,
    double pz,
    double energy)
{
    TLorentzVector p;
    p.SetPxPyPzE(px, py, pz, energy);
    return p;
}


TLorentzVector makeTgeantFinalStateFourVector(
    const Particle& p)
{
    return makeFourVector(
        p.px     * TGEANT_FINALSTATE_TO_GEV,
        p.py     * TGEANT_FINALSTATE_TO_GEV,
        p.pz     * TGEANT_FINALSTATE_TO_GEV,
        p.energy * TGEANT_FINALSTATE_TO_GEV
    );
}


TLorentzVector makeTgeantBeamFourVector(
    const Particle& p)
{
    // Current observed TGEANT convention:
    // beam particle is already stored in GeV.
    return makeFourVector(
        p.px,
        p.py,
        p.pz,
        p.energy
    );
}


double calculateQ2(
    const TLorentzVector& incomingMuon,
    const TLorentzVector& outgoingMuon)
{
    const TLorentzVector q =
        incomingMuon - outgoingMuon;

    return -q.M2();
}


double thetaMrad(const TLorentzVector& p)
{
    return p.Theta() * 1000.0;
}


// ============================================================
// ROOT-particle extraction
// ============================================================

bool getRootParticle(
    int wantedPID,
    const vector<int>* pid,
    const vector<double>* energy,
    const vector<double>* px,
    const vector<double>* py,
    const vector<double>* pz,
    TLorentzVector& result)
{
    if (!pid || !energy || !px || !py || !pz)
        return false;

    for (size_t i = 0; i < pid->size(); ++i) {

        if (pid->at(i) != wantedPID)
            continue;

        result = makeFourVector(
            px->at(i),
            py->at(i),
            pz->at(i),
            energy->at(i)
        );

        return true;
    }

    return false;
}


// ============================================================
// TGEANT-particle extraction
// ============================================================

bool getTgeantBeamMuon(
    const vector<Particle>& particles,
    TLorentzVector& result)
{
    for (const auto& p : particles) {

        if (abs(p.pid) != 13)
            continue;

        result = makeTgeantBeamFourVector(p);
        return true;
    }

    return false;
}


bool getTgeantScatteredMuon(
    const vector<Particle>& particles,
    TLorentzVector& result)
{
    int muonCounter = 0;

    for (const auto& p : particles) {

        if (abs(p.pid) != 13)
            continue;

        ++muonCounter;

        // Current observed ordering:
        // first muon  = beam muon
        // second muon = scattered muon
        if (muonCounter == 2) {
            result = makeTgeantFinalStateFourVector(p);
            return true;
        }
    }

    return false;
}


bool getTgeantPhoton(
    const vector<Particle>& particles,
    TLorentzVector& result)
{
    for (const auto& p : particles) {

        if (p.pid != 22)
            continue;

        result = makeTgeantFinalStateFourVector(p);
        return true;
    }

    return false;
}


// ============================================================
// Logarithmic binning
// ============================================================

vector<double> makeLogBins(
    int nBins,
    double minValue,
    double maxValue)
{
    vector<double> edges(nBins + 1);

    const double logMin = log10(minValue);
    const double logMax = log10(maxValue);

    for (int i = 0; i <= nBins; ++i) {

        const double fraction =
            static_cast<double>(i) / nBins;

        edges[i] =
            pow(10.0, logMin + fraction * (logMax - logMin));
    }

    return edges;
}


// ============================================================
// Plot helpers
// ============================================================

void save2D(
    TH2D* hist,
    const string& outputFile,
    bool logX,
    bool logY,
    bool logZ,
    bool diagonal = true)
{
    TCanvas canvas("canvas2d", hist->GetTitle(), 850, 700);
    canvas.SetRightMargin(0.15);

    if (logX) canvas.SetLogx();
    if (logY) canvas.SetLogy();
    if (logZ) canvas.SetLogz();

    hist->Draw("COLZ");

    if (diagonal) {

        const double xmin = hist->GetXaxis()->GetXmin();
        const double xmax = hist->GetXaxis()->GetXmax();

        TLine line(xmin, xmin, xmax, xmax);
        line.SetLineStyle(2);
        line.SetLineWidth(2);
        line.Draw("SAME");
    }

    canvas.SaveAs(outputFile.c_str());

    cout << "Created: " << outputFile << endl;
}


void saveNormalizedComparison(
    TH1D* before,
    TH1D* after,
    const string& xTitle,
    const string& outputFile,
    bool logY = true)
{
    unique_ptr<TH1D> hBefore(
        static_cast<TH1D*>(before->Clone("before_plot"))
    );

    unique_ptr<TH1D> hAfter(
        static_cast<TH1D*>(after->Clone("after_plot"))
    );

    hBefore->SetDirectory(nullptr);
    hAfter->SetDirectory(nullptr);

    if (hBefore->Integral() > 0.0)
        hBefore->Scale(1.0 / hBefore->Integral());

    if (hAfter->Integral() > 0.0)
        hAfter->Scale(1.0 / hAfter->Integral());

    hBefore->SetLineColor(kBlue + 1);
    hBefore->SetLineWidth(2);

    hAfter->SetLineColor(kRed + 1);
    hAfter->SetLineWidth(2);

    hBefore->GetXaxis()->SetTitle(xTitle.c_str());
    hBefore->GetYaxis()->SetTitle("Normalized entries");

    const double maxY =
        max(hBefore->GetMaximum(), hAfter->GetMaximum());

    if (maxY > 0.0)
        hBefore->SetMaximum(1.25 * maxY);

    TCanvas canvas("canvas1d", "Before vs after TGEANT", 850, 650);

    if (logY)
        canvas.SetLogy();

    hBefore->Draw("HIST");
    hAfter->Draw("HIST SAME");

    TLegend legend(0.68, 0.76, 0.88, 0.88);
    legend.AddEntry(hBefore.get(), "Before TGEANT", "l");
    legend.AddEntry(hAfter.get(),  "After TGEANT",  "l");
    legend.Draw();

    canvas.SaveAs(outputFile.c_str());

    cout << "Created: " << outputFile << endl;
}


void saveSingleHistogram(
    TH1D* hist,
    const string& xTitle,
    const string& yTitle,
    const string& outputFile,
    bool logY = false)
{
    TCanvas canvas("canvas_single", hist->GetTitle(), 850, 650);

    if (logY)
        canvas.SetLogy();

    hist->GetXaxis()->SetTitle(xTitle.c_str());
    hist->GetYaxis()->SetTitle(yTitle.c_str());

    hist->Draw("HIST");

    canvas.SaveAs(outputFile.c_str());

    cout << "Created: " << outputFile << endl;
}


// ============================================================
// Main
// ============================================================

int main(int argc, char* argv[])
{
    if (argc != 4) {

        cerr
            << "Usage:\n"
            << "  " << argv[0]
            << " <McMule_input.root>"
            << " <TGEANT_output.tgeant.gz>"
            << " <output_folder>\n";

        return 1;
    }

    const string rootInputFile = argv[1];
    const string tgeantInputFile = argv[2];
    const string outputFolder = argv[3];

    filesystem::create_directories(outputFolder);


    // --------------------------------------------------------
    // ROOT input: before TGEANT
    // --------------------------------------------------------

    TFile* rootFile =
        TFile::Open(rootInputFile.c_str());

    if (!rootFile || rootFile->IsZombie()) {

        cerr
            << "ERROR: could not open ROOT input file: "
            << rootInputFile
            << endl;

        return 1;
    }


    TTree* tree =
        static_cast<TTree*>(rootFile->Get("Output"));

    if (!tree) {

        cerr
            << "ERROR: TTree 'Output' not found in "
            << rootInputFile
            << endl;

        rootFile->Close();

        return 1;
    }


    // --------------------------------------------------------
    // ROOT branches
    // --------------------------------------------------------

    int beamPID = 0;

    double beamEnergy = 0.0;
    double beamMomentumX = 0.0;
    double beamMomentumY = 0.0;
    double beamMomentumZ = 0.0;

    vector<int>* scatteredPID = nullptr;

    vector<double>* scatteredEnergy = nullptr;
    vector<double>* scatteredMomentumX = nullptr;
    vector<double>* scatteredMomentumY = nullptr;
    vector<double>* scatteredMomentumZ = nullptr;


    tree->SetBranchAddress("beamPID", &beamPID);
    tree->SetBranchAddress("beamEnergy", &beamEnergy);
    tree->SetBranchAddress("beamMomentumX", &beamMomentumX);
    tree->SetBranchAddress("beamMomentumY", &beamMomentumY);
    tree->SetBranchAddress("beamMomentumZ", &beamMomentumZ);

    tree->SetBranchAddress("scatteredPID", &scatteredPID);
    tree->SetBranchAddress("scatteredEnergy", &scatteredEnergy);
    tree->SetBranchAddress("scatteredMomentumX", &scatteredMomentumX);
    tree->SetBranchAddress("scatteredMomentumY", &scatteredMomentumY);
    tree->SetBranchAddress("scatteredMomentumZ", &scatteredMomentumZ);


    const Long64_t nRootEvents =
        tree->GetEntries();


    // --------------------------------------------------------
    // TGEANT input: after TGEANT
    // --------------------------------------------------------

    gzFile tgeantFile =
        gzopen(tgeantInputFile.c_str(), "rb");

    if (!tgeantFile) {

        cerr
            << "ERROR: could not open TGEANT file: "
            << tgeantInputFile
            << endl;

        rootFile->Close();

        return 1;
    }


    // --------------------------------------------------------
    // Histograms
    // --------------------------------------------------------

    gStyle->SetOptStat(0);

    const vector<double> q2Bins =
        makeLogBins(120, Q2_MIN, Q2_MAX);


    TH2D hQ2BeforeAfter(
        "hQ2BeforeAfter",
        "Q^{2}: before vs after TGEANT;"
        "Q^{2}_{before} [GeV^{2}];"
        "Q^{2}_{after} [GeV^{2}]",
        120, q2Bins.data(),
        120, q2Bins.data()
    );


    TH2D hMuonEnergyBeforeAfter(
        "hMuonEnergyBeforeAfter",
        "Scattered muon energy: before vs after TGEANT;"
        "E_{#mu,before} [GeV];"
        "E_{#mu,after} [GeV]",
        150, EMU_MIN, EMU_MAX,
        150, EMU_MIN, EMU_MAX
    );


    TH1D hMuonEnergyBefore(
        "hMuonEnergyBefore",
        "Scattered muon energy",
        250, EMU_MIN, EMU_MAX
    );


    TH1D hMuonEnergyAfter(
        "hMuonEnergyAfter",
        "Scattered muon energy",
        250, EMU_MIN, EMU_MAX
    );


    TH1D hMuonEnergyLoss(
        "hMuonEnergyLoss",
        "Muon energy difference",
        400, -5.0, 5.0
    );


    TH2D hPhotonEnergyBeforeAfter(
        "hPhotonEnergyBeforeAfter",
        "Photon energy: before vs after TGEANT;"
        "E_{#gamma,before} [GeV];"
        "E_{#gamma,after} [GeV]",
        150, PHOTON_E_MIN, PHOTON_E_MAX,
        150, PHOTON_E_MIN, PHOTON_E_MAX
    );


    TH1D hPhotonEnergyBefore(
        "hPhotonEnergyBefore",
        "Photon energy",
        250, PHOTON_E_MIN, PHOTON_E_MAX
    );


    TH1D hPhotonEnergyAfter(
        "hPhotonEnergyAfter",
        "Photon energy",
        250, PHOTON_E_MIN, PHOTON_E_MAX
    );


    TH2D hPhotonAngleBeforeAfter(
        "hPhotonAngleBeforeAfter",
        "Photon angle: before vs after TGEANT;"
        "#theta_{#gamma,before} [mrad];"
        "#theta_{#gamma,after} [mrad]",
        150, PHOTON_THETA_MIN, PHOTON_THETA_MAX,
        150, PHOTON_THETA_MIN, PHOTON_THETA_MAX
    );


    TH1D hPhotonAngleBefore(
        "hPhotonAngleBefore",
        "Photon angle",
        250, PHOTON_THETA_MIN, PHOTON_THETA_MAX
    );


    TH1D hPhotonAngleAfter(
        "hPhotonAngleAfter",
        "Photon angle",
        250, PHOTON_THETA_MIN, PHOTON_THETA_MAX
    );


    // --------------------------------------------------------
    // Counters
    // --------------------------------------------------------

    long long realTgeantEvents = 0;
    long long comparedEvents = 0;

    long long noParticleBlock = 0;
    long long noRootMuon = 0;
    long long noTgeantMuon = 0;

    long long rootPhotonEvents = 0;
    long long comparedPhotonEvents = 0;
    long long missingTgeantPhoton = 0;

    long long nonPositiveQ2Before = 0;
    long long nonPositiveQ2After = 0;

    long long eventNumberMismatch = 0;


    // ========================================================
    // Event-by-event comparison
    // ========================================================

    for (Long64_t i = 0; i < nRootEvents; ++i) {

        long long tgeantEventNumber = -1;
        vector<Particle> particles;
        bool hasParticleBlock = false;


        if (!readNextTgeantEvent(
                tgeantFile,
                tgeantEventNumber,
                particles,
                hasParticleBlock)) {

            cout
                << "Reached end of TGEANT file after "
                << realTgeantEvents
                << " events."
                << endl;

            break;
        }


        ++realTgeantEvents;


        if (tgeantEventNumber != i) {

            ++eventNumberMismatch;

            if (eventNumberMismatch <= 10) {
                cerr
                    << "WARNING: ROOT entry "
                    << i
                    << " is being compared to TGEANT event "
                    << tgeantEventNumber
                    << ". Event-wise matching assumes identical ordering."
                    << endl;
            }
        }


        tree->GetEntry(i);


        if (!hasParticleBlock) {

            ++noParticleBlock;
            continue;
        }


        // ----------------------------------------------------
        // BEFORE TGEANT: ROOT four-vectors
        // ----------------------------------------------------

        const TLorentzVector rootBeam =
            makeFourVector(
                beamMomentumX,
                beamMomentumY,
                beamMomentumZ,
                beamEnergy
            );


        TLorentzVector rootMuon;

        if (!getRootParticle(
                13,
                scatteredPID,
                scatteredEnergy,
                scatteredMomentumX,
                scatteredMomentumY,
                scatteredMomentumZ,
                rootMuon)) {

            ++noRootMuon;
            continue;
        }


        // ----------------------------------------------------
        // AFTER TGEANT: particle block four-vectors
        // ----------------------------------------------------

        TLorentzVector tgeantBeam;
        TLorentzVector tgeantMuon;


        if (!getTgeantBeamMuon(
                particles,
                tgeantBeam) ||
            !getTgeantScatteredMuon(
                particles,
                tgeantMuon)) {

            ++noTgeantMuon;
            continue;
        }


        // ----------------------------------------------------
        // Q^2
        // ----------------------------------------------------

        const double q2Before =
            calculateQ2(rootBeam, rootMuon);

        const double q2After =
            calculateQ2(tgeantBeam, tgeantMuon);


        if (q2Before > 0.0 && q2After > 0.0) {

            hQ2BeforeAfter.Fill(
                q2Before,
                q2After
            );
        }

        if (q2Before <= 0.0)
            ++nonPositiveQ2Before;

        if (q2After <= 0.0)
            ++nonPositiveQ2After;


        // ----------------------------------------------------
        // Muon energy
        // ----------------------------------------------------

        hMuonEnergyBefore.Fill(rootMuon.E());
        hMuonEnergyAfter.Fill(tgeantMuon.E());

        hMuonEnergyBeforeAfter.Fill(
            rootMuon.E(),
            tgeantMuon.E()
        );


        // Positive means:
        // after-TGEANT muon has less energy than before-TGEANT muon.
        const double muonEnergyLoss =
            rootMuon.E() - tgeantMuon.E();

        hMuonEnergyLoss.Fill(
            muonEnergyLoss
        );


        // ----------------------------------------------------
        // Photon
        // ----------------------------------------------------

        TLorentzVector rootPhoton;

        const bool hasRootPhoton =
            getRootParticle(
                22,
                scatteredPID,
                scatteredEnergy,
                scatteredMomentumX,
                scatteredMomentumY,
                scatteredMomentumZ,
                rootPhoton
            );


        if (hasRootPhoton) {

            ++rootPhotonEvents;

            TLorentzVector tgeantPhoton;

            const bool hasTgeantPhoton =
                getTgeantPhoton(
                    particles,
                    tgeantPhoton
                );


            hPhotonEnergyBefore.Fill(
                rootPhoton.E()
            );

            hPhotonAngleBefore.Fill(
                thetaMrad(rootPhoton)
            );


            if (hasTgeantPhoton) {

                ++comparedPhotonEvents;

                hPhotonEnergyAfter.Fill(
                    tgeantPhoton.E()
                );

                hPhotonAngleAfter.Fill(
                    thetaMrad(tgeantPhoton)
                );


                hPhotonEnergyBeforeAfter.Fill(
                    rootPhoton.E(),
                    tgeantPhoton.E()
                );


                hPhotonAngleBeforeAfter.Fill(
                    thetaMrad(rootPhoton),
                    thetaMrad(tgeantPhoton)
                );
            }

            else {

                ++missingTgeantPhoton;
            }
        }


        ++comparedEvents;
    }


    // --------------------------------------------------------
    // Check whether TGEANT contains more events than ROOT
    // --------------------------------------------------------

    long long extraEventNumber = -1;
    vector<Particle> extraParticles;
    bool extraHasParticleBlock = false;

    const bool hasExtraTgeantEvent =
        readNextTgeantEvent(
            tgeantFile,
            extraEventNumber,
            extraParticles,
            extraHasParticleBlock
        );


    if (hasExtraTgeantEvent) {

        cerr
            << "\nWARNING: TGEANT contains more events than the ROOT input."
            << "\nA one-to-one before/after comparison is then not guaranteed."
            << "\nThis can also be a sign that generator events were reused."
            << endl;
    }


    // ========================================================
    // Plots
    // ========================================================

    save2D(
        &hQ2BeforeAfter,
        outputFolder + "/Q2_before_vs_after.pdf",
        true,
        true,
        true,
        true
    );


    save2D(
        &hMuonEnergyBeforeAfter,
        outputFolder + "/muon_energy_before_vs_after.pdf",
        false,
        false,
        true,
        true
    );


    saveNormalizedComparison(
        &hMuonEnergyBefore,
        &hMuonEnergyAfter,
        "Scattered muon energy [GeV]",
        outputFolder + "/muon_energy_distribution.pdf",
        true
    );


    saveSingleHistogram(
        &hMuonEnergyLoss,
        "E_{#mu,before} - E_{#mu,after} [GeV]",
        "Events",
        outputFolder + "/muon_energy_loss.pdf",
        true
    );


    save2D(
        &hPhotonEnergyBeforeAfter,
        outputFolder + "/photon_energy_before_vs_after.pdf",
        false,
        false,
        true,
        true
    );


    saveNormalizedComparison(
        &hPhotonEnergyBefore,
        &hPhotonEnergyAfter,
        "Photon energy [GeV]",
        outputFolder + "/photon_energy_distribution.pdf",
        true
    );


    save2D(
        &hPhotonAngleBeforeAfter,
        outputFolder + "/photon_angle_before_vs_after.pdf",
        false,
        false,
        true,
        true
    );


    saveNormalizedComparison(
        &hPhotonAngleBefore,
        &hPhotonAngleAfter,
        "Photon angle #theta_{#gamma} [mrad]",
        outputFolder + "/photon_angle_distribution.pdf",
        true
    );


    // ========================================================
    // Store histograms
    // ========================================================

    const string histogramFileName =
        outputFolder + "/before_after_histograms.root";

    TFile histogramFile(
        histogramFileName.c_str(),
        "RECREATE"
    );

    hQ2BeforeAfter.Write();

    hMuonEnergyBeforeAfter.Write();
    hMuonEnergyBefore.Write();
    hMuonEnergyAfter.Write();
    hMuonEnergyLoss.Write();

    hPhotonEnergyBeforeAfter.Write();
    hPhotonEnergyBefore.Write();
    hPhotonEnergyAfter.Write();

    hPhotonAngleBeforeAfter.Write();
    hPhotonAngleBefore.Write();
    hPhotonAngleAfter.Write();

    histogramFile.Close();


    // ========================================================
    // Summary
    // ========================================================

    cout << endl;
    cout << "==============================================" << endl;
    cout << "Before / after TGEANT comparison" << endl;
    cout << "==============================================" << endl;

    cout
        << "ROOT events available:            "
        << nRootEvents
        << endl;

    cout
        << "TGEANT events read:               "
        << realTgeantEvents
        << endl;

    cout
        << "Events compared:                  "
        << comparedEvents
        << endl;

    cout
        << "No TGEANT particle block:         "
        << noParticleBlock
        << endl;

    cout
        << "No scattered muon in ROOT:        "
        << noRootMuon
        << endl;

    cout
        << "No scattered muon in TGEANT:      "
        << noTgeantMuon
        << endl;

    cout
        << "ROOT events with photon:          "
        << rootPhotonEvents
        << endl;

    cout
        << "Photon events compared:           "
        << comparedPhotonEvents
        << endl;

    cout
        << "Photon missing in TGEANT block:   "
        << missingTgeantPhoton
        << endl;

    cout
        << "Non-positive Q2 before:           "
        << nonPositiveQ2Before
        << endl;

    cout
        << "Non-positive Q2 after:            "
        << nonPositiveQ2After
        << endl;

    cout
        << "Event-number mismatches:          "
        << eventNumberMismatch
        << endl;

    cout
        << "Histograms written to:            "
        << histogramFileName
        << endl;

    cout << "==============================================" << endl;


    gzclose(tgeantFile);
    rootFile->Close();

    return 0;
}
