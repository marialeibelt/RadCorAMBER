#include <iostream>
#include <string>

#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TString.h"
#include "TSystem.h"


void plotHistogram(TH1D* hist, const char* xTitle, const char* outputFile,
                   const char* outputDir, const char* yTitle = "Events",
                   int lineWidth = 2, const char* drawOption = "HIST")
{
    std::cout << "Plotting " << outputFile << "..." << std::endl;

    TCanvas* canvas = new TCanvas();

    hist->SetLineWidth(lineWidth);
    hist->SetTitle("");
    hist->GetXaxis()->SetTitle(xTitle);
    hist->GetYaxis()->SetTitle(yTitle);
    hist->Draw(drawOption);

    std::cout << "Saving " << outputFile << "..." << std::endl;

    canvas->SaveAs(Form("%s/%s",outputDir,outputFile));

    std::cout << "Saved." << std::endl;

    delete canvas;
}


void plot_analyze_AMBER_eventgen(
    const char* inputFile =
        "/nfs/momos/user/mleibelt/TGEANT_runs/output/PRM_run007_histograms.root",

    const char* outputDir =
        "/nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/Figures/AMBER_eventgen"
)
{
    std::cout << "A: start" << std::endl;
    gSystem->mkdir(outputDir,true);
    std::cout << "B: directory ready" << std::endl;
    TFile* file = TFile::Open(inputFile,"READ");
    std::cout << "C: file opened" << std::endl;

    if (!file || file->IsZombie()) {
        std::cerr << "ERROR: Could not open " << inputFile << std::endl;
        return;
    }

    std::cout << "Reading: " << inputFile << std::endl;


    // --------------------------------------------------------
    // Get histograms
    // --------------------------------------------------------

    std::cout << "Getting histograms..." << std::endl;

    TH1D* h_muon_energy = dynamic_cast<TH1D*>(file->Get("h_muon_energy"));
    TH1D* h_muon_theta = dynamic_cast<TH1D*>(file->Get("h_muon_theta"));
    TH1D* h_Q2 = dynamic_cast<TH1D*>(file->Get("h_Q2"));

    std::cout << "Histograms loaded." << std::endl;

    if (!h_muon_energy || !h_muon_theta || !h_Q2) {
        std::cerr << "ERROR: Could not find all histograms!" << std::endl;
        file->Close();
        delete file;
        return;
    }


    // --------------------------------------------------------
    // General ROOT style
    // --------------------------------------------------------

    gStyle->SetOptStat(0);
    gStyle->SetTitleSize(0.045,"XYZ");
    gStyle->SetLabelSize(0.04,"XYZ");
    gStyle->SetPadLeftMargin(0.13);
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetPadBottomMargin(0.13);
    gStyle->SetPadTopMargin(0.05);


    // --------------------------------------------------------
    // Plots
    // --------------------------------------------------------

    plotHistogram(h_muon_energy,"E_{#mu'} [GeV]","muon_energy.pdf",outputDir);
    plotHistogram(h_muon_theta,"#theta_{#mu'} [mrad]","muon_theta.pdf",outputDir);
    plotHistogram(h_Q2,"Q^{2} [GeV^{2}]","Q2.pdf",outputDir);


    // --------------------------------------------------------
    // Information
    // --------------------------------------------------------

    std::cout << "\n==============================================" << std::endl;
    std::cout << "Plotting finished" << std::endl;
    std::cout << "==============================================" << std::endl;
    std::cout << "Input file:           " << inputFile << std::endl;
    std::cout << "Output directory:     " << outputDir << std::endl;
    std::cout << "Muon energy entries:  " << h_muon_energy->GetEntries() << std::endl;
    std::cout << "Muon theta entries:   " << h_muon_theta->GetEntries() << std::endl;
    std::cout << "Q2 entries:           " << h_Q2->GetEntries() << std::endl;


    file->Close();
    delete file;
}