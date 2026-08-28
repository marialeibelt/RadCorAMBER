#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <utility>

#include <zlib.h>

#include "TFile.h"
#include "TH1D.h"
#include "TLorentzVector.h"


// ============================================================
// Particle
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
// Read one line from gzipped file
// ============================================================

bool readGzLine(
    gzFile file,
    std::string& line
) {

    const int BUFFER_SIZE = 10000;
    char buffer[BUFFER_SIZE];

    if (!gzgets(file, buffer, BUFFER_SIZE))
        return false;

    line = buffer;

    // Remove newline characters
    while (!line.empty() &&
           (line.back() == '\n' ||
            line.back() == '\r')) {

        line.pop_back();
    }

    return true;
}


// ============================================================
// Read first token of a line
// ============================================================

std::string getFirstToken(
    const std::string& line
) {

    std::stringstream ss(line);

    std::string token;
    ss >> token;

    return token;
}


// ============================================================
// Exact event markers
//
// IMPORTANT:
//
// "#EVENT_END" also starts with "#EVENT" as a string.
// Therefore we must compare the complete first token.
// ============================================================

bool isEventStart(
    const std::string& line
) {

    return getFirstToken(line) == "#EVENT";
}


bool isEventEnd(
    const std::string& line
) {

    return getFirstToken(line) == "#EVENT_END";
}


// ============================================================
// Read TGEANT event number
//
// Example:
//
// #EVENT  4  1.00497496  1
//
// -> eventNumber = 4
// ============================================================

bool getEventNumber(
    const std::string& line,
    long long& eventNumber
) {

    std::stringstream ss(line);

    std::string tag;

    if (!(ss >> tag >> eventNumber))
        return false;

    return tag == "#EVENT";
}


// ============================================================
// Check whether a line contains exactly one integer
//
// Example:
//
// 3
//
// This is used as a possible particle-count line.
// ============================================================

bool parseSingleInteger(
    const std::string& line,
    int& value
) {

    std::stringstream ss(line);

    if (!(ss >> value))
        return false;

    // Make sure there is nothing else on the line
    std::string extra;

    if (ss >> extra)
        return false;

    return true;
}


// ============================================================
// Parse one particle line
//
// Current format observed in your TGEANT file:
//
// index
// pid
// dummy
// dummy
// dummy
// px
// py
// pz
// energy
// mass
//
// "valid" here ONLY means:
// the line can be read in this numerical format.
// It does NOT mean physically valid.
// ============================================================

bool parseParticle(
    const std::string& line,
    Particle& p
) {

    std::stringstream ss(line);

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


// ============================================================
// Find particle block inside one TGEANT event
//
// At the moment we know that the correct block has the form:
//
// N
// particle line 1
// particle line 2
// ...
// particle line N
//
// We search ONLY inside one real event:
//
// #EVENT
// ...
// #EVENT_END
//
// Therefore #EVENT_END can no longer produce a fake event.
// ============================================================

bool findParticleBlock(
    const std::vector<std::string>& eventLines,
    std::vector<Particle>& particles
) {

    particles.clear();


    for (size_t i = 0; i < eventLines.size(); ++i) {

        int nParticles = 0;


        // Is this line a single integer?
        if (!parseSingleInteger(
                eventLines[i],
                nParticles)) {

            continue;
        }


        // Reasonable sanity range
        if (nParticles <= 0 ||
            nParticles > 100) {

            continue;
        }


        // Check that enough lines remain
        if (i + static_cast<size_t>(nParticles)
            >= eventLines.size()) {

            continue;
        }


        std::vector<Particle> candidateParticles;

        candidateParticles.reserve(nParticles);

        bool validBlock = true;


        for (int j = 0; j < nParticles; ++j) {

            Particle p;

            if (!parseParticle(
                    eventLines[i + 1 + j],
                    p)) {

                validBlock = false;
                break;
            }

            candidateParticles.push_back(p);
        }


        if (validBlock) {

            particles =
                std::move(candidateParticles);

            return true;
        }
    }


    return false;
}


// ============================================================
// Main
// ============================================================

int main(
    int argc,
    char* argv[]
) {

    // --------------------------------------------------------
    // Arguments
    // --------------------------------------------------------

    if (argc != 2 && argc != 3) {

        std::cerr
            << "Usage:\n"
            << "  "
            << argv[0]
            << " input.tgeant.gz [output.root]\n";

        return 1;
    }


    const std::string inputFile =
        argv[1];


    const std::string outputFile =
        (argc == 3)
        ? argv[2]
        : "tgeant_histograms.root";


    // --------------------------------------------------------
    // Open compressed TGEANT file
    // --------------------------------------------------------

    gzFile file =
        gzopen(
            inputFile.c_str(),
            "rb"
        );


    if (!file) {

        std::cerr
            << "ERROR: Could not open "
            << inputFile
            << "\n";

        return 1;
    }


    std::cout
        << "Reading: "
        << inputFile
        << "\n";


    // --------------------------------------------------------
    // ROOT output
    // --------------------------------------------------------

    TFile output(
        outputFile.c_str(),
        "RECREATE"
    );


    if (output.IsZombie()) {

        std::cerr
            << "ERROR: Could not create "
            << outputFile
            << "\n";

        gzclose(file);

        return 1;
    }


    // ========================================================
    // Histograms
    // ========================================================

    TH1D h_muon_energy(
        "h_muon_energy",
        "Scattered muon energy;"
        "E_{#mu'} [GeV];"
        "Events",
        100,
        0.0,
        100.0
    );


    TH1D h_muon_theta(
        "h_muon_theta",
        "Scattered muon angle;"
        "#theta_{#mu'} [mrad];"
        "Events",
        100,
        0.0,
        5.0
    );


    TH1D h_Q2(
        "h_Q2",
        "Momentum transfer;"
        "Q^{2} [GeV^{2}];"
        "Events",
        100,
        0.0,
        0.05
    );


    // ========================================================
    // Event loop
    // ========================================================

    long long eventCounter = 0;

    long long analyzedEvents = 0;

    long long eventsWithoutParticleBlock = 0;

    long long eventsWithoutTwoMuons = 0;

    long long negativeQ2Events = 0;


    std::string line;


    while (readGzLine(file, line)) {


        // ----------------------------------------------------
        // Search for REAL #EVENT
        //
        // #EVENT_END is NOT accepted here.
        // ----------------------------------------------------

        if (!isEventStart(line))
            continue;


        ++eventCounter;


        // ----------------------------------------------------
        // Read real TGEANT event number
        // ----------------------------------------------------

        long long tgeantEventNumber = -1;

        getEventNumber(
            line,
            tgeantEventNumber
        );


        // ----------------------------------------------------
        // Read complete event until #EVENT_END
        // ----------------------------------------------------

        std::vector<std::string> eventLines;


        bool foundEventEnd = false;


        while (readGzLine(file, line)) {


            if (isEventEnd(line)) {

                foundEventEnd = true;
                break;
            }


            eventLines.push_back(line);
        }


        if (!foundEventEnd) {

            std::cerr
                << "WARNING: No #EVENT_END found for TGEANT event "
                << tgeantEventNumber
                << "\n";

            break;
        }


        // ----------------------------------------------------
        // Find primary particle block
        // ----------------------------------------------------

        std::vector<Particle> particles;


        if (!findParticleBlock(
                eventLines,
                particles)) {


            ++eventsWithoutParticleBlock;


            if (eventCounter <= 10) {

                std::cout
                    << "\nTGEANT Event "
                    << tgeantEventNumber
                    << ": no particle block found\n";
            }


            continue;
        }


        // ----------------------------------------------------
        // Find muons
        // ----------------------------------------------------

        std::vector<size_t> muonIndices;


        for (size_t i = 0;
             i < particles.size();
             ++i) {


            if (std::abs(particles[i].pid) == 13) {

                muonIndices.push_back(i);
            }
        }


        // ----------------------------------------------------
        // Debug print for first 10 real TGEANT events
        // ----------------------------------------------------

        if (eventCounter <= 10) {

            std::cout
                << "\nTGEANT Event "
                << tgeantEventNumber
                << ": particles = "
                << particles.size()
                << ", muons = "
                << muonIndices.size()
                << "\n";


            for (const auto& p : particles) {

                std::cout
                    << "  PID = "
                    << p.pid

                    << " E(raw) = "
                    << p.energy

                    << " px(raw) = "
                    << p.px

                    << " py(raw) = "
                    << p.py

                    << " pz(raw) = "
                    << p.pz

                    << "\n";
            }
        }


        // ----------------------------------------------------
        // Need at least incoming + scattered muon
        // ----------------------------------------------------

        if (muonIndices.size() < 2) {

            ++eventsWithoutTwoMuons;

            continue;
        }


        // ----------------------------------------------------
        // Current ordering observed in your TGEANT output:
        //
        // first muon  = incoming beam muon
        // second muon = scattered muon
        //
        // We should later verify this from the exact TGEANT
        // particle block definition.
        // ----------------------------------------------------

        const Particle& mu_in_raw =
            particles[
                muonIndices[0]
            ];


        const Particle& mu_out_raw =
            particles[
                muonIndices[1]
            ];


        // ====================================================
        // Units
        // ====================================================
        //
        // What we currently observe:
        //
        // incoming:
        //
        //     E ~ 100
        //
        // -> apparently GeV
        //
        // outgoing:
        //
        //     E ~ 100000
        //
        // -> apparently MeV
        //
        // Therefore we currently convert ONLY final-state
        // TGEANT values from MeV to GeV.
        //
        // IMPORTANT:
        // This should later be verified in the TGEANT ROOT
        // input / output code.
        // ====================================================

        // ----------------------------------------------------
        // Incoming muon
        //
        // currently already in GeV
        // ----------------------------------------------------

        TLorentzVector mu_in;

        mu_in.SetPxPyPzE(
            mu_in_raw.px,
            mu_in_raw.py,
            mu_in_raw.pz,
            mu_in_raw.energy
        );


        // ----------------------------------------------------
        // Scattered muon
        //
        // currently MeV -> GeV
        // ----------------------------------------------------

        TLorentzVector mu_out;

        mu_out.SetPxPyPzE(
            mu_out_raw.px,

            mu_out_raw.py,

            mu_out_raw.pz,

            mu_out_raw.energy
        );


        // ====================================================
        // Scattered muon energy
        // ====================================================

        const double E_mu_GeV =
            mu_out.E();


        h_muon_energy.Fill(
            E_mu_GeV
        );


        // ====================================================
        // Scattered muon angle
        //
        // TLorentzVector::Theta() -> radians
        // ====================================================

        const double theta_mrad =
            mu_out.Theta()
            * 1000.0;


        h_muon_theta.Fill(
            theta_mrad
        );


        // ====================================================
        // Q²
        //
        // q = p_in - p_out
        //
        // Q² = -q²
        //
        // Since both four-vectors are now in GeV,
        // Q² automatically comes out in GeV².
        // ====================================================

        const TLorentzVector q =
            mu_in - mu_out;


        const double Q2 =
            -q.M2();


        if (Q2 >= 0.0) {

            h_Q2.Fill(
                Q2
            );
        }

        else {

            ++negativeQ2Events;
        }


        ++analyzedEvents;


        // ====================================================
        // Debug first 10 real events
        // ====================================================

        if (eventCounter <= 10) {

            std::cout
                << "\n  Converted quantities:\n";


            std::cout
                << "    Incoming muon:"
                << " E = "
                << mu_in.E()
                << " GeV"
                << " pz = "
                << mu_in.Pz()
                << " GeV\n";


            std::cout
                << "    Scattered muon:"
                << " E = "
                << mu_out.E()
                << " GeV"
                << " pz = "
                << mu_out.Pz()
                << " GeV\n";


            std::cout
                << "    theta = "
                << theta_mrad
                << " mrad\n";


            std::cout
                << "    Q2 = "
                << Q2
                << " GeV^2\n";
        }
    }


    // ========================================================
    // Finish
    // ========================================================

    gzclose(file);


    std::cout
        << "\n==============================================\n";

    std::cout
        << "TGEANT analysis summary\n";

    std::cout
        << "==============================================\n";


    std::cout
        << "Real #EVENT blocks found:   "
        << eventCounter
        << "\n";


    std::cout
        << "Events analyzed:            "
        << analyzedEvents
        << "\n";


    std::cout
        << "No particle block:          "
        << eventsWithoutParticleBlock
        << "\n";


    std::cout
        << "Fewer than two muons:       "
        << eventsWithoutTwoMuons
        << "\n";


    std::cout
        << "Negative Q2 values:         "
        << negativeQ2Events
        << "\n";


    // ========================================================
    // Write ROOT file
    // ========================================================

    output.cd();

    h_muon_energy.Write();
    h_muon_theta.Write();
    h_Q2.Write();

    output.Close();


    std::cout
        << "Histograms written to:      "
        << outputFile
        << "\n";

    std::cout
        << "==============================================\n";


    return 0;
}