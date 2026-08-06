#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>

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
            abs(p.pid)==13
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




int main(int argc,char **argv)
{

    std::ofstream logFile("lhe_to_root.log");

    if(!logFile.is_open()) {

        std::cerr << "Cannot open log file\n";
        return 1;
    }



    if(argc!=3) {

        std::cout
        << "Usage: ./lhe_to_root input.lhe output.root\n";

        return 1;
    }



    std::string inputFile=argv[1];
    std::string outputFile=argv[2];



    std::ifstream lhe(inputFile);


    if(!lhe.is_open()) {

        std::cerr
        << "Cannot open "
        << inputFile
        << "\n";

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


    double vertexX=0.0;
    double vertexY=0.0;
    double vertexZ=-3200.0;



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



    tree.Branch(
        "weight",
        &weight
    );


    tree.Branch("vertexX",&vertexX);
    tree.Branch("vertexY",&vertexY);
    tree.Branch("vertexZ",&vertexZ);


    tree.Branch(
        "beamPID",
        &beamPID
    );


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



    int counter=0;


    while(
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


        tree.Fill();

        counter++;


        if(counter%10000==0)
            logFile
            << counter
            << " events processed\n";

    }



    outfile.Write();

    outfile.Close();



    logFile
    << "Finished. Events written: "
    << counter
    << "\n";

    logFile.close();


    return 0;
}