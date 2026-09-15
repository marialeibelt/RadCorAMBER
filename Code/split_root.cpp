#include <iostream>
#include <string>
#include <algorithm>

#include "TFile.h"
#include "TTree.h"

int main(int argc, char* argv[]) {

    if (argc != 4) {
        std::cerr << "Usage:\n"
                  << argv[0]
                  << " input.root output_folder batch_size\n";
        return 1;
    }

    std::string inputFile  = argv[1];
    std::string baseName = inputFile;

    // Pfad entfernen
    size_t slash = baseName.find_last_of("/\\");
    if (slash != std::string::npos) {
        baseName = baseName.substr(slash + 1);
    }

    // .root entfernen
    if (baseName.size() >= 5 &&
        baseName.substr(baseName.size() - 5) == ".root") {
        baseName = baseName.substr(0, baseName.size() - 5);
    }

    std::string outputDir  = argv[2];
    Long64_t batchSize     = std::stoll(argv[3]);

    // --------------------------------------------------
    // Open input
    // --------------------------------------------------

    TFile *fin = TFile::Open(inputFile.c_str(), "READ");

    if (!fin || fin->IsZombie()) {
        std::cerr << "ERROR: Could not open input file:\n"
                  << inputFile << std::endl;
        return 1;
    }

    TTree *tree = dynamic_cast<TTree*>(fin->Get("Output"));

    if (!tree) {
        std::cerr << "ERROR: TTree 'Output' not found.\n";
        fin->Close();
        return 1;
    }

    Long64_t nEntries = tree->GetEntries();

    std::cout << "Input file: " << inputFile << "\n";
    std::cout << "Entries:    " << nEntries << "\n";
    std::cout << "Batch size: " << batchSize << "\n\n";

    // --------------------------------------------------
    // Split
    // --------------------------------------------------

    Long64_t start = 0;
    int batch = 0;

    while (start < nEntries) {

        Long64_t nThisBatch =
            std::min(batchSize, nEntries - start);

        Long64_t end = start + nThisBatch;

        std::string outputFile =
            outputDir
            + "/"
            + baseName
            + "_batch"
            + std::to_string(batch)
            + "_"
            + std::to_string(nThisBatch)
            + ".root";

        std::cout
            << "Batch " << batch
            << ": entries "
            << start << " - " << end - 1
            << " (" << nThisBatch << " events)"
            << "\n -> " << outputFile
            << std::endl;

        TFile *fout =
            TFile::Open(outputFile.c_str(), "RECREATE");

        TTree *newTree =
            tree->CloneTree(0);

        for (Long64_t i = start; i < end; ++i) {
            tree->GetEntry(i);
            newTree->Fill();
        }

        newTree->Write();
        fout->Close();

        delete fout;

        start = end;
        batch++;
    }

    fin->Close();

    std::cout << "\nDone. Created "
              << batch
              << " ROOT files."
              << std::endl;

    return 0;
}