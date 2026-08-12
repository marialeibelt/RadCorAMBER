#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>

#include "TFile.h"
#include "TTree.h"


struct Particle {

    int pid;
    int status;

    double px;
    double py;
    double pz;
    double energy;
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
)
{
    std::string line;

    while (std::getline(file, line)) {

        if (line.find("<event") != std::string::npos)
            break;
    }

    if (file.eof())
        return false;


    // ========================================================
    // Event header
    // ========================================================

    std::getline(file, line);

    std::stringstream header(line);

    int nParticles;
    int idprup;
    double scale;
    double alphaEM;
    double alphaS;

    header
        >> nParticles
        >> idprup
        >> weight
        >> scale
        >> alphaEM
        >> alphaS;


    // ========================================================
    // Read particles
    // ========================================================

    std::vector<Particle> particles;

    for (int i = 0; i < nParticles; i++) {

        std::getline(file, line);

        std::stringstream ss(line);

        Particle p;

        int mother1, mother2;
        int color1, color2;

        ss
            >> p.pid
            >> p.status
            >> mother1
            >> mother2
            >> color1
            >> color2
            >> p.px
            >> p.py
            >> p.pz
            >> p.energy;

        // mass, lifetime, spin
        double dummy;

        ss >> dummy;
        ss >> dummy;
        ss >> dummy;


        // McMule photon convention
        if (p.pid == 0)
            p.pid = 22;

        particles.push_back(p);
    }

    // Skip </event>
    std::getline(file, line);

    // ========================================================
    // Extract beam and scattered particles
    // ========================================================

    bool foundBeam = false;

    scatteredPID.clear();
    scatteredEnergy.clear();
    scatteredMomentumX.clear();
    scatteredMomentumY.clear();
    scatteredMomentumZ.clear();


    for (auto &p : particles) {

        // Incoming muon
        if (
            p.status == -1 &&
            std::abs(p.pid) == 13
        ) {
            beamPID = p.pid;

            beamEnergy = p.energy;

            beamMomentumX = p.px;
            beamMomentumY = p.py;
            beamMomentumZ = p.pz;

            foundBeam = true;
        }


        // Final-state particles
        if (p.status == 1) {

            scatteredPID.push_back(p.pid);

            scatteredEnergy.push_back(p.energy);

            scatteredMomentumX.push_back(p.px);
            scatteredMomentumY.push_back(p.py);
            scatteredMomentumZ.push_back(p.pz);
        }
    }

    if (!foundBeam) {
        std::cerr << "Warning: no incoming muon found\n";
    }

    return true;
}


int main(int argc, char **argv)
{
    // ========================================================
    // Command line arguments
    // ========================================================

    if (argc != 3) {
        std::cout << "Usage: ./lhe_to_root input.lhe output.root\n";
        return 1;
    }

    std::string inputFile = argv[1];
    std::string outputFile = argv[2];


    // ========================================================
    // Log file
    // ========================================================

    std::ofstream logFile("lhe_to_root.log");

    if (!logFile.is_open()) {
        std::cerr << "Cannot open log file\n";
        return 1;
    }


    // ========================================================
    // Open LHE file
    // ========================================================

    std::ifstream lhe(inputFile);

    if (!lhe.is_open()) {
        std::cerr << "Cannot open " << inputFile << "\n";
        return 1;
    }

    std::cout << "Reading LHE file: " << inputFile << std::endl;


    // ========================================================
    // Create ROOT file
    // ========================================================

    TFile outfile(outputFile.c_str(),"RECREATE");

    if (outfile.IsZombie()) {
        std::cerr << "Error: could not create ROOT file " << outputFile << "\n";
        return 1;
    }

    TTree tree("Output","McMule events");


    // ========================================================
    // Variables
    // ========================================================

    double weight;

    double vertexX = 0.0;
    double vertexY = 0.0;
    double vertexZ = -3200.0;


    int beamPID;

    double beamEnergy;

    double beamMomentumX;
    double beamMomentumY;
    double beamMomentumZ;


    std::vector<int> scatteredPID;

    std::vector<double> scatteredEnergy;

    std::vector<double> scatteredMomentumX;
    std::vector<double> scatteredMomentumY;
    std::vector<double> scatteredMomentumZ;


    // ========================================================
    // ROOT branches
    // ========================================================

    tree.Branch("weight", &weight);

    tree.Branch("vertexX", &vertexX);
    tree.Branch("vertexY", &vertexY);
    tree.Branch("vertexZ", &vertexZ);

    tree.Branch("beamPID", &beamPID);
    tree.Branch("beamEnergy", &beamEnergy);
    tree.Branch("beamMomentumX", &beamMomentumX);
    tree.Branch("beamMomentumY", &beamMomentumY);
    tree.Branch("beamMomentumZ", &beamMomentumZ);

    tree.Branch("scatteredPID", &scatteredPID);
    tree.Branch("scatteredEnergy", &scatteredEnergy);
    tree.Branch("scatteredMomentumX", &scatteredMomentumX);
    tree.Branch("scatteredMomentumY", &scatteredMomentumY);
    tree.Branch("scatteredMomentumZ", &scatteredMomentumZ);


    // ========================================================
    // Read all events
    // ========================================================

    long long totalEvents = 0;
    long long negativeEvents = 0;


    while (
        readEvent(
            lhe,
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
        )
    ) {

        totalEvents++;

        if (weight < 0.0)
            negativeEvents++;
        tree.Fill();

        if (totalEvents % 100000 == 0) {
            std::cout << totalEvents << " events written" << std::endl;
            logFile << totalEvents << " events written\n";
        }
    }

    outfile.Write();
    outfile.Close();

    lhe.close();


    logFile
        << "Finished.\n"
        << "Total events processed: " << totalEvents << "\n"
        << "Negative events: " << negativeEvents << "\n"
        << "Events written to ROOT: " << totalEvents << "\n";

    logFile.close();

    std::cout << "\nFinished." << std::endl;
    std::cout << "Total events processed: " << totalEvents << std::endl;
    std::cout << "Negative events: " << negativeEvents << std::endl;
    std::cout << "Events written to ROOT: " << totalEvents << std::endl;
    std::cout << "Output file: " << outputFile << std::endl;

    return 0;
}