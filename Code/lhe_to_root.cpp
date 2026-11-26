#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <random>
#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"

struct Particle {

    int pid = 0;
    int status = 0;

    double px = 0.0;
    double py = 0.0;
    double pz = 0.0;
    double energy = 0.0;
};


bool readEvent(
    std::ifstream &file,

    double &weight,

    int &beamPID,
    double &beamEnergy,

    double &beamMomentumX,
    double &beamMomentumY,
    double &beamMomentumZ,

    std::vector<int> &scatteredPID,
    std::vector<double> &scatteredEnergy,

    std::vector<double> &scatteredMomentumX,
    std::vector<double> &scatteredMomentumY,
    std::vector<double> &scatteredMomentumZ
) {

    std::string line;

    // --------------------------------------------------
    // Find next <event> block
    // --------------------------------------------------

    bool foundEvent = false;

    while (std::getline(file, line)) {

        if (line.find("<event") != std::string::npos &&
            line.find("</event") == std::string::npos) {

            foundEvent = true;
            break;
        }
    }

    if (!foundEvent)
        return false;


    // --------------------------------------------------
    // Read event header
    // --------------------------------------------------

    if (!std::getline(file, line)) {
        throw std::runtime_error(
            "Unexpected end of file while reading event header."
        );
    }

    std::stringstream header(line);

    int nParticles = 0;
    int idprup = 0;

    double scale = 0.0;
    double alphaEM = 0.0;
    double alphaS = 0.0;

    if (!(header
        >> nParticles
        >> idprup
        >> weight
        >> scale
        >> alphaEM
        >> alphaS)) {

        throw std::runtime_error(
            "Could not parse LHE event header:\n" + line
        );
    }


    // --------------------------------------------------
    // Read particles
    // --------------------------------------------------

    std::vector<Particle> particles;
    particles.reserve(nParticles);

    for (int i = 0; i < nParticles; ++i) {

        if (!std::getline(file, line)) {
            throw std::runtime_error(
                "Unexpected end of file while reading particle."
            );
        }

        std::stringstream ss(line);

        Particle p;

        int mother1 = 0;
        int mother2 = 0;
        int color1 = 0;
        int color2 = 0;

        double mass = 0.0;
        double lifetime = 0.0;
        double spin = 0.0;

        if (!(ss
            >> p.pid
            >> p.status
            >> mother1
            >> mother2
            >> color1
            >> color2
            >> p.px
            >> p.py
            >> p.pz
            >> p.energy
            >> mass
            >> lifetime
            >> spin)) {

            throw std::runtime_error(
                "Could not parse LHE particle line:\n" + line
            );
        }


        // McMule photon convention
        if (p.pid == 0)
            p.pid = 22;

        particles.push_back(p);
    }


    // --------------------------------------------------
    // Move to end of event
    // --------------------------------------------------

    bool foundEventEnd = false;

    while (std::getline(file, line)) {

        if (line.find("</event>") != std::string::npos) {
            foundEventEnd = true;
            break;
        }
    }

    if (!foundEventEnd) {
        throw std::runtime_error(
            "Could not find </event>."
        );
    }


    // --------------------------------------------------
    // Reset output variables
    // --------------------------------------------------

    beamPID = 0;

    beamEnergy = 0.0;

    beamMomentumX = 0.0;
    beamMomentumY = 0.0;
    beamMomentumZ = 0.0;

    scatteredPID.clear();
    scatteredEnergy.clear();

    scatteredMomentumX.clear();
    scatteredMomentumY.clear();
    scatteredMomentumZ.clear();


    // --------------------------------------------------
    // Identify incoming muon and final-state particles
    // --------------------------------------------------

    bool foundBeam = false;

    for (const auto &p : particles) {

        // incoming muon
        if (p.status == -1 &&
            std::abs(p.pid) == 13) {

            if (foundBeam) {
                throw std::runtime_error(
                    "More than one incoming muon found in event."
                );
            }

            beamPID = p.pid;

            beamEnergy = p.energy;

            beamMomentumX = p.px;
            beamMomentumY = p.py;
            beamMomentumZ = p.pz;

            foundBeam = true;
        }


        // final-state particles
        if (p.status == 1) {

            scatteredPID.push_back(p.pid);

            scatteredEnergy.push_back(p.energy);

            scatteredMomentumX.push_back(p.px);
            scatteredMomentumY.push_back(p.py);
            scatteredMomentumZ.push_back(p.pz);
        }
    }


    if (!foundBeam) {
        throw std::runtime_error(
            "No incoming muon found in event."
        );
    }


    return true;
}


double findMaxWeight(
    const std::string &inputFile,
    long long &totalEvents,
    long long &negativeEvents
) {

    std::ifstream file(inputFile);

    if (!file.is_open()) {
        throw std::runtime_error(
            "Cannot open input file: " + inputFile
        );
    }


    double maxWeight = 0.0;
    double weight = 0.0;

    int beamPID = 0;

    double beamEnergy = 0.0;

    double beamMomentumX = 0.0;
    double beamMomentumY = 0.0;
    double beamMomentumZ = 0.0;

    std::vector<int> scatteredPID;
    std::vector<double> scatteredEnergy;

    std::vector<double> scatteredMomentumX;
    std::vector<double> scatteredMomentumY;
    std::vector<double> scatteredMomentumZ;


    totalEvents = 0;
    negativeEvents = 0;


    while (readEvent(
        file,
        weight,
        beamPID,
        beamEnergy,
        beamMomentumX,
        beamMomentumY,
        beamMomentumZ,
        scatteredPID,
        scatteredEnergy,
        scatteredMomentumX,
        scatteredMomentumY,
        scatteredMomentumZ
    )) {

        ++totalEvents;

        if (!std::isfinite(weight)) {
            throw std::runtime_error(
                "Non-finite event weight encountered."
            );
        }

        if (weight < 0.0) {

            ++negativeEvents;

            // Current strategy:
            // negative NLO events are not used.
            continue;
        }

        maxWeight = std::max(maxWeight, weight);
    }


    return maxWeight;
}


int main(int argc, char **argv) {

    try {

        // ==================================================
        // Beam / frame definition
        // ==================================================

        const double Ebeam = 100.0;       // GeV
        const double mmu   = 0.105658375; // GeV
        const double mp    = 0.938272088; // GeV

        const double pbeam =
            std::sqrt(Ebeam * Ebeam - mmu * mmu);

        // Velocity of the CMS relative to the lab frame
        const double beta =
            pbeam / (Ebeam + mp);


        if (argc != 3) {

            std::cout
                << "Usage: ./lhe_to_root input.lhe output.root\n";

            return 1;
        }


        const std::string inputFile = argv[1];
        const std::string outputFile = argv[2];


        // ==================================================
        // Log file
        // ==================================================

        std::ofstream logFile("lhe_to_root.log");

        if (!logFile.is_open()) {
            std::cerr << "Cannot open log file.\n";
            return 1;
        }


        // ==================================================
        // Random number generator
        // ==================================================

        // Fixed seed -> reproducible results while debugging
        const unsigned long long seed = 123456789ULL;

        std::mt19937_64 rng(seed);

        std::uniform_real_distribution<double> uniform(0.0, 1.0);


        // ==================================================
        // First pass: maximum positive weight
        // ==================================================

        long long totalEvents = 0;
        long long negativeEvents = 0;

        const double maxWeight =
            findMaxWeight(
                inputFile,
                totalEvents,
                negativeEvents
            );


        std::cout
            << "Total events: "
            << totalEvents
            << "\n";

        std::cout
            << "Negative events: "
            << negativeEvents
            << "\n";

        std::cout
            << "Maximum positive weight: "
            << maxWeight
            << "\n";

        std::cout
            << "CMS -> lab beta: "
            << beta
            << "\n";


        if (maxWeight <= 0.0) {

            std::cerr
                << "Error: no positive event weights found.\n";

            return 1;
        }


        // ==================================================
        // Second pass
        // ==================================================

        std::ifstream lhe(inputFile);

        if (!lhe.is_open()) {

            std::cerr
                << "Cannot open "
                << inputFile
                << "\n";

            return 1;
        }


        // ==================================================
        // ROOT output
        // ==================================================

        TFile outfile(
            outputFile.c_str(),
            "RECREATE"
        );

        if (outfile.IsZombie()) {

            std::cerr
                << "Could not create ROOT file "
                << outputFile
                << "\n";

            return 1;
        }


        TTree tree(
            "Output",
            "McMule events"
        );


        // unweighted output weight
        double weight = 1.0;

        // original McMule weight
        double generatorWeight = 0.0;


        double vertexX = 0.0;
        double vertexY = 0.0;
        double vertexZ = -3200.0;


        int beamPID = 0;

        double beamEnergy = 0.0;

        double beamMomentumX = 0.0;
        double beamMomentumY = 0.0;
        double beamMomentumZ = 0.0;


        std::vector<int> scatteredPID;
        std::vector<double> scatteredEnergy;

        std::vector<double> scatteredMomentumX;
        std::vector<double> scatteredMomentumY;
        std::vector<double> scatteredMomentumZ;


        // ==================================================
        // ROOT branches
        // ==================================================

        tree.Branch(
            "weight",
            &weight
        );

        tree.Branch(
            "generatorWeight",
            &generatorWeight
        );


        tree.Branch("vertexX", &vertexX);
        tree.Branch("vertexY", &vertexY);
        tree.Branch("vertexZ", &vertexZ);


        tree.Branch("beamPID", &beamPID);

        tree.Branch(
            "beamEnergy",
            &beamEnergy
        );

        tree.Branch(
            "beamMomentumX",
            &beamMomentumX
        );

        tree.Branch(
            "beamMomentumY",
            &beamMomentumY
        );

        tree.Branch(
            "beamMomentumZ",
            &beamMomentumZ
        );


        tree.Branch(
            "scatteredPID",
            &scatteredPID
        );

        tree.Branch(
            "scatteredEnergy",
            &scatteredEnergy
        );

        tree.Branch(
            "scatteredMomentumX",
            &scatteredMomentumX
        );

        tree.Branch(
            "scatteredMomentumY",
            &scatteredMomentumY
        );

        tree.Branch(
            "scatteredMomentumZ",
            &scatteredMomentumZ
        );


        // ==================================================
        // Derived radiative-correction observables (lab frame)
        //
        // All quantities follow uniquely from the McMule
        // four-vectors. Photon quantities refer to the hardest
        // photon; they are set to NO_VALUE if there is none.
        // ECAL2 hit = straight-line projection of the photon from
        // the vertex to the ECAL2 plane (no magnetic field acts on
        // photons). Units: GeV, GeV^2, mrad, mm.
        // ==================================================

        const double NO_VALUE = -9999.0;

        // ECAL2 position from TGEANT share/prm_2022/alignment_prm.xml (mm)
        const double ecal2X = 450.0;
        const double ecal2Y = -3.968;
        const double ecal2Z = 33252.0;

        Long64_t eventID = 0;

        int nPhotons = 0;
        int hasPhoton = 0;

        double photonEnergy = NO_VALUE;
        double photonTheta = NO_VALUE;        // mrad, w.r.t. z axis
        double photonPhi = NO_VALUE;          // rad
        double photonAngleToBeam = NO_VALUE;  // mrad, to incoming muon  (ISR-like if small)
        double photonAngleToMuon = NO_VALUE;  // mrad, to outgoing muon  (FSR-like if small)
        double photonECAL2X = NO_VALUE;       // mm, lab
        double photonECAL2Y = NO_VALUE;       // mm, lab
        double photonECAL2XRel = NO_VALUE;    // mm, relative to ECAL2 centre
        double photonECAL2YRel = NO_VALUE;    // mm, relative to ECAL2 centre

        double muonTheta = NO_VALUE;          // mrad
        double Q2mu = NO_VALUE;               // -(k - k')^2
        double Q2p = NO_VALUE;                // -(P - P')^2 = 2 m_p T_p
        double protonKineticEnergy = NO_VALUE;
        double protonTheta = NO_VALUE;        // mrad
        double energyBalance = NO_VALUE;      // E_beam - E_mu' - T_p (= E_gamma)

        tree.Branch("eventID", &eventID);
        tree.Branch("nPhotons", &nPhotons);
        tree.Branch("hasPhoton", &hasPhoton);
        tree.Branch("photonEnergy", &photonEnergy);
        tree.Branch("photonTheta", &photonTheta);
        tree.Branch("photonPhi", &photonPhi);
        tree.Branch("photonAngleToBeam", &photonAngleToBeam);
        tree.Branch("photonAngleToMuon", &photonAngleToMuon);
        tree.Branch("photonECAL2X", &photonECAL2X);
        tree.Branch("photonECAL2Y", &photonECAL2Y);
        tree.Branch("photonECAL2XRel", &photonECAL2XRel);
        tree.Branch("photonECAL2YRel", &photonECAL2YRel);
        tree.Branch("muonTheta", &muonTheta);
        tree.Branch("Q2mu", &Q2mu);
        tree.Branch("Q2p", &Q2p);
        tree.Branch("protonKineticEnergy", &protonKineticEnergy);
        tree.Branch("protonTheta", &protonTheta);
        tree.Branch("energyBalance", &energyBalance);


        // Separate small tree: one entry per photon.
        // Link back to "Output" via eventID
        // (e.g. Output->BuildIndex("eventID")).
        TTree photonTree(
            "Photons",
            "McMule photons (one entry per photon)"
        );

        photonTree.Branch("eventID", &eventID);
        photonTree.Branch("weight", &weight);
        photonTree.Branch("photonEnergy", &photonEnergy);
        photonTree.Branch("photonTheta", &photonTheta);
        photonTree.Branch("photonPhi", &photonPhi);
        photonTree.Branch("photonAngleToBeam", &photonAngleToBeam);
        photonTree.Branch("photonAngleToMuon", &photonAngleToMuon);
        photonTree.Branch("photonECAL2X", &photonECAL2X);
        photonTree.Branch("photonECAL2Y", &photonECAL2Y);
        photonTree.Branch("photonECAL2XRel", &photonECAL2XRel);
        photonTree.Branch("photonECAL2YRel", &photonECAL2YRel);
        photonTree.Branch("muonTheta", &muonTheta);
        photonTree.Branch("Q2mu", &Q2mu);
        photonTree.Branch("Q2p", &Q2p);
        photonTree.Branch("protonKineticEnergy", &protonKineticEnergy);
        photonTree.Branch("energyBalance", &energyBalance);


        // Photon quantities for one photon four-vector
        auto fillPhoton = [&](const TLorentzVector &g,
                              const TLorentzVector &k,
                              const TLorentzVector &kp,
                              bool haveMuon) {

            photonEnergy = g.E();
            photonTheta = 1000.0 * g.Theta();
            photonPhi = g.Phi();
            photonAngleToBeam = 1000.0 * g.Angle(k.Vect());
            photonAngleToMuon =
                haveMuon ? 1000.0 * g.Angle(kp.Vect()) : NO_VALUE;

            if (g.Pz() > 0.0) {
                const double dz = ecal2Z - vertexZ;
                photonECAL2X = vertexX + g.Px() / g.Pz() * dz;
                photonECAL2Y = vertexY + g.Py() / g.Pz() * dz;
                photonECAL2XRel = photonECAL2X - ecal2X;
                photonECAL2YRel = photonECAL2Y - ecal2Y;
            } else {
                photonECAL2X = photonECAL2Y = NO_VALUE;
                photonECAL2XRel = photonECAL2YRel = NO_VALUE;
            }
        };


        // ==================================================
        // Event loop
        // ==================================================

        long long acceptedEvents = 0;


        while (readEvent(
            lhe,
            generatorWeight,
            beamPID,
            beamEnergy,
            beamMomentumX,
            beamMomentumY,
            beamMomentumZ,
            scatteredPID,
            scatteredEnergy,
            scatteredMomentumX,
            scatteredMomentumY,
            scatteredMomentumZ
        )) {


            // ------------------------------------------------
            // Current treatment of negative NLO events
            // ------------------------------------------------

            if (generatorWeight < 0.0)
                continue;


            // ------------------------------------------------
            // Accept/reject unweighting
            // ------------------------------------------------

            double probability =
                generatorWeight / maxWeight;

            // Numerical safety
            if (probability > 1.0)
                probability = 1.0;


            const double randomNumber =
                uniform(rng);


            if (randomNumber >= probability)
                continue;


            // Output event is unweighted
            weight = 1.0;


            // ==================================================
            // CMS -> LAB
            //
            // IMPORTANT:
            // CMS moving in +z relative to lab.
            // Therefore CMS -> lab uses +beta.
            // ==================================================

            TLorentzVector beam;

            beam.SetPxPyPzE(
                beamMomentumX,
                beamMomentumY,
                beamMomentumZ,
                beamEnergy
            );


            beam.Boost(
                0.0,
                0.0,
                +beta
            );


            beamMomentumX = beam.Px();
            beamMomentumY = beam.Py();
            beamMomentumZ = beam.Pz();
            beamEnergy     = beam.E();


            // ------------------------------------------------
            // Boost final-state particles
            // ------------------------------------------------

            for (size_t i = 0;
                 i < scatteredPID.size();
                 ++i) {

                TLorentzVector p;

                p.SetPxPyPzE(
                    scatteredMomentumX[i],
                    scatteredMomentumY[i],
                    scatteredMomentumZ[i],
                    scatteredEnergy[i]
                );


                p.Boost(
                    0.0,
                    0.0,
                    +beta
                );


                scatteredMomentumX[i] =
                    p.Px();

                scatteredMomentumY[i] =
                    p.Py();

                scatteredMomentumZ[i] =
                    p.Pz();

                scatteredEnergy[i] =
                    p.E();
            }


            // ------------------------------------------------
            // Derived observables (lab frame)
            // ------------------------------------------------

            eventID = acceptedEvents;

            nPhotons = 0;
            hasPhoton = 0;
            photonEnergy = photonTheta = photonPhi = NO_VALUE;
            photonAngleToBeam = photonAngleToMuon = NO_VALUE;
            photonECAL2X = photonECAL2Y = NO_VALUE;
            photonECAL2XRel = photonECAL2YRel = NO_VALUE;
            muonTheta = Q2mu = Q2p = NO_VALUE;
            protonKineticEnergy = protonTheta = energyBalance = NO_VALUE;

            const TLorentzVector k(
                beamMomentumX, beamMomentumY, beamMomentumZ, beamEnergy
            );
            const TLorentzVector P(0.0, 0.0, 0.0, mp);  // proton at rest in lab

            TLorentzVector kp, Pp;
            bool haveMuon = false;
            bool haveProton = false;
            std::vector<TLorentzVector> photons;

            for (size_t i = 0; i < scatteredPID.size(); ++i) {

                const TLorentzVector v(
                    scatteredMomentumX[i],
                    scatteredMomentumY[i],
                    scatteredMomentumZ[i],
                    scatteredEnergy[i]
                );

                if (std::abs(scatteredPID[i]) == 13 && !haveMuon) {
                    kp = v;
                    haveMuon = true;
                } else if (scatteredPID[i] == 2212 && !haveProton) {
                    Pp = v;
                    haveProton = true;
                } else if (scatteredPID[i] == 22) {
                    photons.push_back(v);
                }
            }

            if (haveMuon) {
                muonTheta = 1000.0 * kp.Theta();
                Q2mu = -(k - kp).M2();
            }

            if (haveProton) {
                Q2p = -(P - Pp).M2();
                protonKineticEnergy = Pp.E() - mp;
                protonTheta = 1000.0 * Pp.Theta();
            }

            if (haveMuon && haveProton) {
                energyBalance = k.E() - kp.E() - protonKineticEnergy;
            }

            nPhotons = static_cast<int>(photons.size());
            hasPhoton = nPhotons > 0 ? 1 : 0;

            // One Photons-tree entry per photon
            for (const auto &g : photons) {
                fillPhoton(g, k, kp, haveMuon);
                photonTree.Fill();
            }

            // Output tree keeps the hardest photon
            if (hasPhoton) {
                const auto hardest = std::max_element(
                    photons.begin(), photons.end(),
                    [](const TLorentzVector &a, const TLorentzVector &b) {
                        return a.E() < b.E();
                    }
                );
                fillPhoton(*hardest, k, kp, haveMuon);
            }


            ++acceptedEvents;

            tree.Fill();
        }


        // ==================================================
        // Write ROOT file
        // ==================================================

        outfile.cd();

        tree.Write();

        photonTree.Write();

        outfile.Close();


        // ==================================================
        // Log
        // ==================================================

        logFile
            << "Finished.\n"
            << "Random seed: "
            << seed << "\n"
            << "CMS -> lab beta: "
            << beta << "\n"
            << "Total events processed: "
            << totalEvents << "\n"
            << "Negative events: "
            << negativeEvents << "\n"
            << "Maximum positive weight: "
            << maxWeight << "\n"
            << "Unweighted events written: "
            << acceptedEvents << "\n";


        if (totalEvents > 0) {

            logFile
                << "Acceptance fraction: "
                << static_cast<double>(acceptedEvents)
                   / static_cast<double>(totalEvents)
                << "\n";
        }


        return 0;
    }

    catch (const std::exception &e) {

        std::cerr
            << "ERROR: "
            << e.what()
            << "\n";

        return 1;
    }
}