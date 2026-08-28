#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <random>
#include <cmath>
#include <string>
#include <filesystem>

#include "TFile.h"
#include "TTree.h"

using namespace std;


// ============================================================
// Read cross sections from txt file
// ============================================================

bool readCrossSections(const string& filename,
                       double& sigmaLO,
                       double& sigmaOnlyR,
                       double& sigmaNLO,
                       double& sigmaFull)
{
    ifstream file(filename);

    if (!file.is_open()) {
        cerr << "ERROR: Could not open cross section file: "
             << filename << endl;
        return false;
    }

    // Skip comments and read the four numbers
    string line;
    vector<double> values;

    while (getline(file, line)) {

        // Ignore empty lines
        if (line.empty())
            continue;

        // Ignore comments
        if (line[0] == '#')
            continue;

        try {
            values.push_back(stod(line));
        }
        catch (...) {
            cerr << "WARNING: Could not read line: "
                 << line << endl;
        }
    }

    file.close();

    if (values.size() < 4) {
        cerr << "ERROR: Cross section file contains fewer than 4 values."
             << endl;
        return false;
    }

    sigmaLO    = values[0];
    sigmaOnlyR = values[1];
    sigmaNLO   = values[2];
    sigmaFull  = values[3];

    return true;
}


// ============================================================
// Main
// ============================================================

int main(int argc, char* argv[])
{
    // --------------------------------------------------------
    // Check command line arguments
    // --------------------------------------------------------

    if (argc != 5 && argc != 6) {
        cerr << "Usage:\n"
            << "  " << argv[0]
            << " LO.root NLO.root cross_sections.txt output_folder [suffix]\n";
        return 1;
    }

    string loFileName  = argv[1];
    string nloFileName = argv[2];
    string csFileName  = argv[3];
    string outputDir   = argv[4];

    // Optional suffix, e.g. "run0", "test", "xi01"
    string suffix = "";

    if (argc == 6) {
        suffix = "_" + string(argv[5]);
    }

    std::filesystem::create_directories(outputDir);

    // --------------------------------------------------------
    // Read cross sections
    // --------------------------------------------------------

    double sigmaLO;
    double sigmaOnlyR;
    double sigmaNLO;
    double sigmaFull;

    if (!readCrossSections(
            csFileName,
            sigmaLO,
            sigmaOnlyR,
            sigmaNLO,
            sigmaFull))
    {
        return 1;
    }

    cout << "==============================================" << endl;
    cout << "Cross sections" << endl;
    cout << "==============================================" << endl;

    cout << "sigma_LO    = " << sigmaLO    << endl;
    cout << "sigma_onlyR = " << sigmaOnlyR << endl;
    cout << "sigma_NLO   = " << sigmaNLO   << endl;
    cout << "sigma_full  = " << sigmaFull  << endl;

    double fractionLO  = sigmaLO  / sigmaFull;
    double fractionNLO = sigmaNLO / sigmaFull;

    cout << "LO fraction  = " << fractionLO  << endl;
    cout << "NLO fraction = " << fractionNLO << endl;

    double fractionSum = fractionLO + fractionNLO;

    cout << "Fraction sum = " << fractionSum << endl;

    if (std::abs(fractionSum - 1.0) > 1e-6) {

        cerr << "WARNING: LO and NLO fractions do not add up to 1." << endl;
        cerr << "Check the meaning of sigma_LO, sigma_NLO and sigma_full."
            << endl;
    }


    // --------------------------------------------------------
    // Open input ROOT files
    // --------------------------------------------------------

    TFile* loFile = TFile::Open(loFileName.c_str(), "READ");

    if (!loFile || loFile->IsZombie()) {
        cerr << "ERROR: Could not open " << loFileName << endl;
        return 1;
    }

    TFile* nloFile = TFile::Open(nloFileName.c_str(), "READ");

    if (!nloFile || nloFile->IsZombie()) {
        cerr << "ERROR: Could not open " << nloFileName << endl;
        return 1;
    }


    // --------------------------------------------------------
    // Get trees
    // --------------------------------------------------------

    TTree* loTree = nullptr;
    TTree* nloTree = nullptr;

    loFile->GetObject("Output", loTree);
    nloFile->GetObject("Output", nloTree);

    if (!loTree) {
        cerr << "ERROR: Could not find tree 'Output' in "
             << loFileName << endl;
        return 1;
    }

    if (!nloTree) {
        cerr << "ERROR: Could not find tree 'Output' in "
             << nloFileName << endl;
        return 1;
    }


    Long64_t nLOavailable  = loTree->GetEntries();
    Long64_t nNLOavailable = nloTree->GetEntries();

    cout << "Available LO events  = "
         << nLOavailable << endl;

    cout << "Available NLO events = "
         << nNLOavailable << endl;

    cout << endl;


    // --------------------------------------------------------
    // Determine maximum possible statistics
    //
    // We want:
    //
    // N_LO  / N_total = sigma_LO  / sigma_total
    // N_NLO / N_total = sigma_NLO / sigma_total
    //
    // Therefore:
    //
    // N_total <= N_LO_available  / fractionLO
    // N_total <= N_NLO_available / fractionNLO
    // --------------------------------------------------------

    double maxNtotalFromLO =
        static_cast<double>(nLOavailable) / fractionLO;

    double maxNtotalFromNLO =
        static_cast<double>(nNLOavailable) / fractionNLO;

    double maxNtotal =
        min(maxNtotalFromLO, maxNtotalFromNLO);

    Long64_t Ntotal =
        static_cast<Long64_t>(floor(maxNtotal));

    // --------------------------------------------------------
    // Determine required number of events
    // --------------------------------------------------------

    Long64_t NLOtarget =
        static_cast<Long64_t>(
            llround(Ntotal * fractionNLO)
        );

    Long64_t LOtarget =
        Ntotal - NLOtarget;


    // --------------------------------------------------------
    // Safety checks
    // --------------------------------------------------------

    if (LOtarget > nLOavailable) {

        cerr << "ERROR: Requested "
             << LOtarget
             << " LO events, but only "
             << nLOavailable
             << " are available." << endl;

        return 1;
    }

    if (NLOtarget > nNLOavailable) {

        cerr << "ERROR: Requested "
             << NLOtarget
             << " NLO events, but only "
             << nNLOavailable
             << " are available." << endl;

        return 1;
    }


    cout << "==============================================" << endl;
    cout << "Final event sample" << endl;
    cout << "==============================================" << endl;

    cout << "N_total = " << Ntotal << endl;
    cout << "N_LO    = " << LOtarget << endl;
    cout << "N_NLO   = " << NLOtarget << endl;

    cout << endl;

    cout << "Actually used LO events:  "
         << LOtarget << " / " << nLOavailable << endl;

    cout << "Actually used NLO events: "
         << NLOtarget << " / " << nNLOavailable << endl;

    cout << endl;


    // --------------------------------------------------------
    // Randomly select LO events
    // --------------------------------------------------------

    vector<Long64_t> loIndices(nLOavailable);

    for (Long64_t i = 0; i < nLOavailable; ++i)
        loIndices[i] = i;


    // Fixed seed for reproducibility
    mt19937 rng(123456789);

    shuffle(loIndices.begin(), loIndices.end(), rng);

    loIndices.resize(LOtarget);


    // --------------------------------------------------------
    // Randomly select NLO events
    //
    // If NLOtarget == nNLOavailable, this simply uses all
    // NLO events.
    // --------------------------------------------------------

    vector<Long64_t> nloIndices(nNLOavailable);

    for (Long64_t i = 0; i < nNLOavailable; ++i)
        nloIndices[i] = i;

    shuffle(nloIndices.begin(), nloIndices.end(), rng);

    nloIndices.resize(NLOtarget);


    // --------------------------------------------------------
    // Output files
    // --------------------------------------------------------

    string loOutputName =
        outputDir + "/LO_scaled" + suffix + ".root";

    string nloOutputName =
        outputDir + "/NLO_scaled" + suffix + ".root";


    cout << "Output files:" << endl;
    cout << "  " << loOutputName << endl;
    cout << "  " << nloOutputName << endl;


    TFile* loOut =
        TFile::Open(loOutputName.c_str(), "RECREATE");

    TFile* nloOut =
        TFile::Open(nloOutputName.c_str(), "RECREATE");


    if (!loOut || loOut->IsZombie()) {
        cerr << "ERROR: Could not create "
            << loOutputName << endl;
        return 1;
    }

    if (!nloOut || nloOut->IsZombie()) {
        cerr << "ERROR: Could not create "
            << nloOutputName << endl;
        return 1;
    }


    // --------------------------------------------------------
    // Create output trees
    //
    // Clone the complete tree structure but start empty.
    // --------------------------------------------------------

    loOut->cd();
    TTree* loOutTree = loTree->CloneTree(0);

    nloOut->cd();
    TTree* nloOutTree = nloTree->CloneTree(0);


    // --------------------------------------------------------
    // Fill LO output
    // --------------------------------------------------------

    cout << "Writing LO events..." << endl;

    for (Long64_t index : loIndices) {

        loTree->GetEntry(index);
        loOutTree->Fill();
    }


    // --------------------------------------------------------
    // Fill NLO output
    // --------------------------------------------------------

    cout << "Writing NLO events..." << endl;

    for (Long64_t index : nloIndices) {

        nloTree->GetEntry(index);
        nloOutTree->Fill();
    }


    // --------------------------------------------------------
    // Write files
    // --------------------------------------------------------

    loOut->cd();
    loOutTree->Write();

    nloOut->cd();
    nloOutTree->Write();

    loOut->Close();
    nloOut->Close();

    loFile->Close();
    nloFile->Close();


    // --------------------------------------------------------
    // Finished
    // --------------------------------------------------------

    cout << endl;
    cout << "==============================================" << endl;
    cout << "Finished!" << endl;
    cout << "==============================================" << endl;

    cout << "Created:" << endl;
    cout << "  " << loOutputName << endl;
    cout << "  " << nloOutputName << endl;

    return 0;
}