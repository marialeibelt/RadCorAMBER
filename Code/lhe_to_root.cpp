#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <random>
#include <algorithm>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"

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

) {


    std::string line;


    while(std::getline(file,line)) {

        if(line.find("<event") != std::string::npos)
            break;

    }


    if(file.eof())
        return false;


    // event header
    std::getline(file,line);

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



    std::vector<Particle> particles;



    for(int i=0;i<nParticles;i++) {


        std::getline(file,line);


        std::stringstream ss(line);


        Particle p;


        int mother1,mother2;
        int color1,color2;


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
        if(p.pid==0)
            p.pid=22;


        particles.push_back(p);

    }


    // skip </event>
    std::getline(file,line);



    bool foundBeam=false;


    scatteredPID.clear();
    scatteredEnergy.clear();
    scatteredMomentumX.clear();
    scatteredMomentumY.clear();
    scatteredMomentumZ.clear();



    for(auto &p : particles) {


        if(
            p.status==-1 &&
            std::abs(p.pid)==13
        ) {

            beamPID=p.pid;

            beamEnergy=p.energy;

            beamMomentumX=p.px;
            beamMomentumY=p.py;
            beamMomentumZ=p.pz;

            foundBeam=true;
        }



        if(p.status==1) {

            scatteredPID.push_back(p.pid);

            scatteredEnergy.push_back(p.energy);

            scatteredMomentumX.push_back(p.px);
            scatteredMomentumY.push_back(p.py);
            scatteredMomentumZ.push_back(p.pz);

        }

    }


    if(!foundBeam) {

        std::cerr
        << "Warning: no incoming muon found\n";

    }


    return true;

}


double findMaxWeight(
    const std::string &inputFile,
    long long &totalEvents,
    long long &negativeEvents
)
{
    std::ifstream file(inputFile);

    if (!file.is_open()) {
        std::cerr << "Cannot open " << inputFile << "\n";
        return 0.0;
    }

    double maxWeight = 0.0;
    double weight;

    int beamPID;
    double beamEnergy;
    double beamMomentumX, beamMomentumY, beamMomentumZ;

    std::vector<int> scatteredPID;
    std::vector<double> scatteredEnergy;
    std::vector<double> scatteredMomentumX;
    std::vector<double> scatteredMomentumY;
    std::vector<double> scatteredMomentumZ;

    totalEvents = 0;
    negativeEvents = 0;

    while (readEvent(
        file, weight,
        beamPID, beamEnergy,
        beamMomentumX, beamMomentumY, beamMomentumZ,
        scatteredPID, scatteredEnergy,
        scatteredMomentumX, scatteredMomentumY, scatteredMomentumZ
    )) {
        totalEvents++;

        if (weight < 0.0) {
            negativeEvents++;
            continue;
        }

        maxWeight = std::max(maxWeight, weight);
    }

    return maxWeight;
}

int main(int argc, char **argv)
{
    const double Ebeam = 100.0;       // GeV
    const double mmu   = 0.105658375; // GeV
    const double mp    = 0.938272088; // GeV

    const double pbeam = std::sqrt(Ebeam * Ebeam - mmu * mmu);

    const double beta = pbeam / (Ebeam + mp);

    if (argc != 3) {
        std::cout << "Usage: ./lhe_to_root input.lhe output.root\n";
        return 1;
    }

    std::string inputFile = argv[1];
    std::string outputFile = argv[2];

    std::ofstream logFile("lhe_to_root.log");

    if (!logFile.is_open()) {
        std::cerr << "Cannot open log file\n";
        return 1;
    }

    // Random number generator
    std::random_device rd;
    std::mt19937_64 rng(rd());

    // First pass: find maximum positive weight
    long long totalEvents;
    long long negativeEvents;

    double maxWeight = findMaxWeight(
        inputFile,
        totalEvents,
        negativeEvents
    );
    
    std::cout << "Total events: "
              << totalEvents << std::endl;

    std::cout << "Negative events: "
              << negativeEvents << std::endl;

    std::cout << "Maximum positive weight: "
              << maxWeight << std::endl;

    if (maxWeight <= 0.0) {
        std::cerr << "Error: no positive event weights found.\n";
        return 1;
    }

    // Second pass: read events and perform accept/reject
    std::ifstream lhe(inputFile);

    if (!lhe.is_open()) {
        std::cerr << "Cannot open "
                  << inputFile << "\n";
        return 1;
    }


    TFile outfile(
        outputFile.c_str(),
        "RECREATE"
    );

    TTree tree(
        "Output",
        "McMule events"
    );


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


    // ROOT branches
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


    long long acceptedEvents = 0;


    while (readEvent(
        lhe, weight,
        beamPID, beamEnergy,
        beamMomentumX, beamMomentumY, beamMomentumZ,
        scatteredPID, scatteredEnergy,
        scatteredMomentumX, scatteredMomentumY, scatteredMomentumZ
    )) {

        // Ignore negative-weight events
        if (randomNumber < probability) {

            weight = 1.0;

            // =========================
            // CMS -> Lab boost
            // =========================

            // Boost incoming muon
            TLorentzVector beam;

            beam.SetPxPyPzE(
                beamMomentumX,
                beamMomentumY,
                beamMomentumZ,
                beamEnergy
            );

            beam.Boost(0.0, 0.0, beta);

            beamMomentumX = beam.Px();
            beamMomentumY = beam.Py();
            beamMomentumZ = beam.Pz();
            beamEnergy     = beam.E();


            // Boost final-state particles
            for (size_t i = 0; i < scatteredPID.size(); ++i) {

                TLorentzVector p;

                p.SetPxPyPzE(
                    scatteredMomentumX[i],
                    scatteredMomentumY[i],
                    scatteredMomentumZ[i],
                    scatteredEnergy[i]
                );

                p.Boost(0.0, 0.0, beta);

                scatteredMomentumX[i] = p.Px();
                scatteredMomentumY[i] = p.Py();
                scatteredMomentumZ[i] = p.Pz();
                scatteredEnergy[i]    = p.E();
            }


            // =========================
            // Write event
            // =========================

            acceptedEvents++;

            tree.Fill();
        }
    }


    outfile.Write();
    outfile.Close();


    logFile
        << "Finished.\n"
        << "Total events processed: " << totalEvents << "\n"
        << "Negative events: " << negativeEvents << "\n"
        << "Maximum positive weight: " << maxWeight << "\n"
        << "Unweighted events written: " << acceptedEvents << "\n";

    logFile.close();

    return 0;
}