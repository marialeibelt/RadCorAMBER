// AMBER TGEANT ASCII -> ROOT: generator/primary-particle truth + histograms.
//
// CAUTION: This parser reads a candidate *primary particle block* in each
// #EVENT ... #EVENT_END block. It does NOT parse and archive detector hits,
// trajectories, energy deposits or all other TGEANT ASCII data.
// Never delete .tgeant.gz based solely on this ROOT output.
//
// Usage: analyze_AMBER_eventgen input.tgeant.gz output.root

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <zlib.h>

#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"
#include "TLorentzVector.h"

struct Particle {
    int index = 0, pid = 0;
    double dummy1 = 0, dummy2 = 0, dummy3 = 0;
    double px = 0, py = 0, pz = 0, energy = 0, mass = 0;
};

// Unlike a fixed buffer, this reads even unusually long TGEANT lines intact.
static bool readGzLine(gzFile file, std::string& line) {
    line.clear();
    char buffer[8192];
    while (true) {
        char* p = gzgets(file, buffer, sizeof(buffer));
        if (!p) return !line.empty();
        line += buffer;
        if (!line.empty() && line.back() == '\n') break;
    }
    while (!line.empty() && (line.back() == '\n' || line.back() == '\r'))
        line.pop_back();
    return true;
}

static std::string firstToken(const std::string& s) {
    std::istringstream ss(s);
    std::string token;
    ss >> token;
    return token;
}

static bool oneInteger(const std::string& s, int& n) {
    std::istringstream ss(s);
    std::string other;
    if (!(ss >> n)) return false;
    return !(ss >> other);
}

static bool parseParticle(const std::string& s, Particle& p) {
    std::istringstream ss(s);
    if (!(ss >> p.index >> p.pid >> p.dummy1 >> p.dummy2 >> p.dummy3
             >> p.px >> p.py >> p.pz >> p.energy >> p.mass)) return false;
    // Only check the format, not the physical interpretation of its columns.
    return std::isfinite(p.px) && std::isfinite(p.py) &&
           std::isfinite(p.pz) && std::isfinite(p.energy);
}

static bool findPrimaryBlock(const std::vector<std::string>& lines,
                             std::vector<Particle>& particles) {
    // Require at least two muons to reduce the risk of choosing a different
    // integer-headed detector block. This heuristic must be checked against
    // an actual TGEANT ASCII format example before relying on it for truth.
    for (size_t i = 0; i < lines.size(); ++i) {
        int n = 0;
        if (!oneInteger(lines[i], n) || n <= 0 || n > 1000) continue;
        if (i + static_cast<size_t>(n) >= lines.size()) continue;
        std::vector<Particle> candidate;
        candidate.reserve(n);
        bool ok = true;
        for (int j = 0; j < n; ++j) {
            Particle p;
            if (!parseParticle(lines[i + 1 + j], p)) { ok = false; break; }
            candidate.push_back(p);
        }
        if (!ok) continue;
        int nMuon = 0;
        for (const auto& p : candidate) if (std::abs(p.pid) == 13) ++nMuon;
        if (nMuon >= 2) {
            particles = std::move(candidate);
            return true;
        }
    }
    return false;
}

static double energyScaleToGeV(double e) {
    // Observed TGEANT format can mix an ~100 GeV incoming beam muon with an
    // ~100000 MeV outgoing muon. Preserve RAW values as well as conversion.
    // A value <1000 is tentatively GeV, >1000 tentatively MeV.
    // Treat these as explicitly diagnosed HEURISTICS, not certified units.
    return std::abs(e) > 1000.0 ? 0.001 : 1.0;
}

static double angleMrad(const TLorentzVector& p) {
    return 1000.0 * std::atan2(std::hypot(p.Px(), p.Py()), p.Pz());
}

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " input.tgeant.gz output.root\n";
        return 2;
    }
    const std::string inputFile = argv[1];
    const std::string outputFile = argv[2];
    gzFile gz = gzopen(inputFile.c_str(), "rb");
    if (!gz) {
        std::cerr << "ERROR: cannot open " << inputFile << "\n";
        return 2;
    }
    TFile output(outputFile.c_str(), "RECREATE");
    if (output.IsZombie()) {
        std::cerr << "ERROR: cannot create " << outputFile << "\n";
        gzclose(gz);
        return 2;
    }
    output.SetCompressionLevel(7);

    TH1D h_muon_energy("h_muon_energy", "Scattered Muon Energy;E_{#mu}^{lab} (GeV);Events / bin",
                        500, 0., 101.);
    TH1D h_muon_theta("h_muon_theta", "Muon Scattering Angle;#theta_{#mu}^{lab} (mrad);Events / bin",
                       500, 0., 8.);
    TH1D h_Q2("h_Q2", "Momentum Transfer;Q^{2} (GeV^{2});Events / bin",
               500, 0., 0.18);
    TH1D h_muon_energy_loss("h_muon_energy_loss", "Muon Energy Loss;#Delta E_{#mu} (GeV);Events / bin",
                            500, 0., 100.);
    TH1D h_photon_energy("h_photon_energy", "Real Photon Energy;E_{#gamma}^{lab} (GeV);Photons / bin",
                         1000, 0., 100.);
    TH1D h_photon_theta("h_photon_theta", "Real Photon Polar Angle;#theta_{#gamma}^{lab} (mrad);Photons / bin",
                        500, 0., 50.);
    TH1D h_photon_theta_full("h_photon_theta_full", "Real Photon Polar Angle (Full Range);#theta_{#gamma}^{lab} (mrad);Photons / bin",
                             500, 0., 3200.);
    // This histogram merges by addition with hadd, enabling a reliable
    // event-count summary without depending on ROOT's GetEntries bookkeeping.
    TH1D h_event_summary("h_event_summary", "TGEANT Event Summary;Category;Events",
                         7, 0.5, 7.5);
    h_event_summary.GetXaxis()->SetBinLabel(1, "event_blocks");
    h_event_summary.GetXaxis()->SetBinLabel(2, "primary_block_found");
    h_event_summary.GetXaxis()->SetBinLabel(3, "two_muons_found");
    h_event_summary.GetXaxis()->SetBinLabel(4, "valid_muon_kinematics");
    h_event_summary.GetXaxis()->SetBinLabel(5, "Q2_nonnegative");
    h_event_summary.GetXaxis()->SetBinLabel(6, "events_with_primary_photon");
    h_event_summary.GetXaxis()->SetBinLabel(7, "primary_photons_total");
    for (TH1D* h : {&h_muon_energy, &h_muon_theta, &h_Q2,
                    &h_muon_energy_loss, &h_photon_energy,
                    &h_photon_theta, &h_photon_theta_full, &h_event_summary})
        h->Sumw2();

    // One row per TGEANT #EVENT, including events without a valid particle
    // block. Only the identified primary particle block is preserved.
    TTree truth("TGEANTTruth", "Parsed primary particles and MC truth (NOT detector hits)");
    std::string sourceFile = inputFile, rawEventHeader;
    Long64_t eventNumber = -1;
    Bool_t primaryBlockFound = false, validMuonKinematics = false;
    Int_t nMuons = 0, nPhotons = 0;
    Double_t muonEnergyGeV = std::numeric_limits<double>::quiet_NaN();
    Double_t thetaMuonMrad = std::numeric_limits<double>::quiet_NaN();
    Double_t Q2GeV2 = std::numeric_limits<double>::quiet_NaN();
    std::vector<int> pIndex, pPID;
    std::vector<double> pDummy1Raw, pDummy2Raw, pDummy3Raw;
    std::vector<double> pPxRaw, pPyRaw, pPzRaw, pEnergyRaw, pMassRaw;
    std::vector<double> pScaleToGeV, pPxGeV, pPyGeV, pPzGeV, pEnergyGeV;

    truth.Branch("sourceFile", &sourceFile);
    truth.Branch("eventNumber", &eventNumber);
    truth.Branch("rawEventHeader", &rawEventHeader);
    truth.Branch("primaryBlockFound", &primaryBlockFound);
    truth.Branch("nMuons", &nMuons);
    truth.Branch("nPrimaryPhotons", &nPhotons);
    truth.Branch("validMuonKinematics", &validMuonKinematics);
    truth.Branch("muonEnergy_GeV", &muonEnergyGeV);
    truth.Branch("thetaMuon_mrad", &thetaMuonMrad);
    truth.Branch("Q2_GeV2", &Q2GeV2);
    truth.Branch("particleIndex", &pIndex);
    truth.Branch("particlePID", &pPID);
    truth.Branch("particleDummy1Raw", &pDummy1Raw);
    truth.Branch("particleDummy2Raw", &pDummy2Raw);
    truth.Branch("particleDummy3Raw", &pDummy3Raw);
    truth.Branch("particlePxRaw", &pPxRaw);
    truth.Branch("particlePyRaw", &pPyRaw);
    truth.Branch("particlePzRaw", &pPzRaw);
    truth.Branch("particleEnergyRaw", &pEnergyRaw);
    truth.Branch("particleMassRaw", &pMassRaw);
    truth.Branch("particleScaleToGeV", &pScaleToGeV);
    truth.Branch("particlePx_GeV", &pPxGeV);
    truth.Branch("particlePy_GeV", &pPyGeV);
    truth.Branch("particlePz_GeV", &pPzGeV);
    truth.Branch("particleEnergy_GeV", &pEnergyGeV);

    long long foundEvents = 0, noBlock = 0, fewerThanTwoMuons = 0;
    long long suspectEnergy = 0, negativeQ2 = 0, eventsWithPhoton = 0;
    long long photonRecords = 0, unitMismatch = 0, invalidEventEnd = 0;
    std::string line;
    while (readGzLine(gz, line)) {
        if (firstToken(line) != "#EVENT") continue;
        ++foundEvents;
        h_event_summary.Fill(1.);
        rawEventHeader = line;
        eventNumber = -1;
        {
            std::istringstream ss(line);
            std::string tag;
            ss >> tag >> eventNumber;
        }
        primaryBlockFound = false;
        validMuonKinematics = false;
        nMuons = nPhotons = 0;
        muonEnergyGeV = thetaMuonMrad = Q2GeV2 = std::numeric_limits<double>::quiet_NaN();
        pIndex.clear(); pPID.clear();
        pDummy1Raw.clear(); pDummy2Raw.clear(); pDummy3Raw.clear();
        pPxRaw.clear(); pPyRaw.clear(); pPzRaw.clear(); pEnergyRaw.clear(); pMassRaw.clear();
        pScaleToGeV.clear(); pPxGeV.clear(); pPyGeV.clear(); pPzGeV.clear(); pEnergyGeV.clear();

        std::vector<std::string> eventLines;
        bool ended = false;
        while (readGzLine(gz, line)) {
            if (firstToken(line) == "#EVENT_END") { ended = true; break; }
            eventLines.push_back(line);
        }
        if (!ended) {
            ++invalidEventEnd;
            std::cerr << "ERROR: truncated event (missing #EVENT_END), event="
                      << eventNumber << " in " << inputFile << "\n";
            break;
        }

        std::vector<Particle> particles;
        primaryBlockFound = findPrimaryBlock(eventLines, particles);
        if (!primaryBlockFound) {
            ++noBlock;
            truth.Fill();
            continue;
        }
        h_event_summary.Fill(2.);

        std::vector<size_t> muonIndices;
        for (size_t i = 0; i < particles.size(); ++i) {
            const Particle& p = particles[i];
            pIndex.push_back(p.index);
            pPID.push_back(p.pid);
            pDummy1Raw.push_back(p.dummy1);
            pDummy2Raw.push_back(p.dummy2);
            pDummy3Raw.push_back(p.dummy3);
            pPxRaw.push_back(p.px);
            pPyRaw.push_back(p.py);
            pPzRaw.push_back(p.pz);
            pEnergyRaw.push_back(p.energy);
            pMassRaw.push_back(p.mass);
            // Incoming muon (first record) may have different units to
            // outgoing particles in this ASCII format; use muon positions.
            if (std::abs(p.pid) == 13) muonIndices.push_back(i);
            if (p.pid == 22) ++nPhotons;
        }
        nMuons = static_cast<Int_t>(muonIndices.size());
        if (nPhotons) {
            ++eventsWithPhoton;
            photonRecords += nPhotons;
            h_event_summary.Fill(6.);
            h_event_summary.Fill(7., nPhotons);
        }
        if (muonIndices.size() < 2) {
            ++fewerThanTwoMuons;
            truth.Fill();
            continue;
        }
        h_event_summary.Fill(3.);

        const Particle& beam = particles[muonIndices[0]];
        const Particle& scattered = particles[muonIndices[1]];
        double beamScale = energyScaleToGeV(beam.energy);
        double finalScale = energyScaleToGeV(scattered.energy);
        if (beamScale != finalScale) ++unitMismatch;

        pScaleToGeV.reserve(particles.size());
        pPxGeV.reserve(particles.size());
        pPyGeV.reserve(particles.size());
        pPzGeV.reserve(particles.size());
        pEnergyGeV.reserve(particles.size());
        for (size_t i = 0; i < particles.size(); ++i) {
            const double scale = (i == muonIndices[0]) ? beamScale : finalScale;
            const Particle& p = particles[i];
            pScaleToGeV.push_back(scale);
            pPxGeV.push_back(p.px * scale);
            pPyGeV.push_back(p.py * scale);
            pPzGeV.push_back(p.pz * scale);
            pEnergyGeV.push_back(p.energy * scale);
        }
        TLorentzVector mu_in, mu_out;
        const size_t iIn = muonIndices[0], iOut = muonIndices[1];
        mu_in.SetPxPyPzE(pPxGeV[iIn], pPyGeV[iIn], pPzGeV[iIn], pEnergyGeV[iIn]);
        mu_out.SetPxPyPzE(pPxGeV[iOut], pPyGeV[iOut], pPzGeV[iOut], pEnergyGeV[iOut]);
        muonEnergyGeV = mu_out.E();
        thetaMuonMrad = angleMrad(mu_out);
        Q2GeV2 = -(mu_in - mu_out).M2();
        validMuonKinematics = std::isfinite(thetaMuonMrad) &&
                              std::isfinite(muonEnergyGeV) &&
                              muonEnergyGeV > 0 && muonEnergyGeV <= 110 &&
                              mu_in.E() > 0 && mu_in.E() <= 110;
        if (validMuonKinematics) {
            h_event_summary.Fill(4.);
            h_muon_energy.Fill(muonEnergyGeV);
            h_muon_theta.Fill(thetaMuonMrad);
            h_muon_energy_loss.Fill(mu_in.E() - mu_out.E());
            if (std::isfinite(Q2GeV2) && Q2GeV2 >= 0) {
                h_Q2.Fill(Q2GeV2);
                h_event_summary.Fill(5.);
            } else {
                ++negativeQ2;
            }
        } else {
            ++suspectEnergy;
            if (suspectEnergy <= 3) {
                std::cerr << "WARNING: suspect muon kinematics in " << inputFile
                          << " event=" << eventNumber << " E_in=" << mu_in.E()
                          << " E_out=" << mu_out.E() << " (GeV, heuristic)\n";
            }
        }
        for (size_t i = 0; i < particles.size(); ++i) {
            if (particles[i].pid != 22) continue;
            const double E = pEnergyGeV[i];
            if (std::isfinite(E) && E >= 0) h_photon_energy.Fill(E);
            TLorentzVector photon;
            photon.SetPxPyPzE(pPxGeV[i], pPyGeV[i], pPzGeV[i], pEnergyGeV[i]);
            double angle = angleMrad(photon);
            if (std::isfinite(angle)) {
                h_photon_theta.Fill(angle);
                h_photon_theta_full.Fill(angle);
            }
        }
        truth.Fill();
    }

    int status = Z_OK;
    gzerror(gz, &status);
    gzclose(gz);
    if (invalidEventEnd || (status != Z_OK && status != Z_STREAM_END)) {
        std::cerr << "ERROR: incomplete or corrupt TGEANT input " << inputFile << "\n";
        output.Close();
        return 3;
    }
    output.cd();
    for (TH1D* h : {&h_muon_energy, &h_muon_theta, &h_Q2,
                    &h_muon_energy_loss, &h_photon_energy,
                    &h_photon_theta, &h_photon_theta_full, &h_event_summary})
        h->Write();
    truth.Write();
    output.Close();
    std::cout << "TGEANT events found: " << foundEvents << "\n"
              << "Primary particle blocks: " << foundEvents - noBlock << "\n"
              << "Events with fewer than 2 muons: " << fewerThanTwoMuons << "\n"
              << "Invalid muon kinematics: " << suspectEnergy << "\n"
              << "Q2 negative/nonfinite: " << negativeQ2 << "\n"
              << "Events with primary photon: " << eventsWithPhoton << "\n"
              << "Primary photon records: " << photonRecords << "\n"
              << "Mixed beam/outgoing units (heuristic): " << unitMismatch << "\n"
              << "Output (histograms + TGEANTTruth): " << outputFile << "\n"
              << "NOTE: detector hits and full raw event text NOT archived.\n";
    return 0;
}
