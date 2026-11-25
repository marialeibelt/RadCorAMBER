// Plot ROOT histograms made by the PHAST UserEvent module.
#include "TCanvas.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TKey.h"
#include "TObject.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include <iostream>
#include <memory>

void plot_analyze_phast_prm(const char* histogramFile, const char* figureFolder) {
    std::unique_ptr<TFile> in(TFile::Open(histogramFile,"READ"));
    if (!in || in->IsZombie()) {
        std::cerr << "ERROR: Cannot open PHAST histogram file: " << histogramFile << "\n";
        return;
    }
    TDirectory* dir = in->GetDirectory("AMBER_PRM");
    if (!dir) {
        std::cerr << "ERROR: No AMBER_PRM directory in ROOT file.\n";
        in->ls();
        return;
    }
    gSystem->mkdir(figureFolder,true);
    gStyle->SetPalette(kBird);
    gStyle->SetNumberContours(60);
    TIter next(dir->GetListOfKeys());
    int count=0;
    while (TObject* item = next()) {
        auto* key = dynamic_cast<TKey*>(item);
        if (!key) continue;
        std::unique_ptr<TObject> obj(key->ReadObj());
        auto* hist = dynamic_cast<TH1*>(obj.get());
        if (!hist) continue;
        TCanvas c(TString::Format("canvas_%d",count),hist->GetName(),1100,750);
        c.SetLeftMargin(0.13); c.SetBottomMargin(0.13);
        const bool twoD = hist->InheritsFrom(TH2::Class());
        c.SetRightMargin(twoD ? 0.15 : 0.05);
        gStyle->SetOptStat(twoD ? 0 : 1110);
        hist->SetLineWidth(2);
        hist->Draw(twoD ? "COLZ" : "HIST");
        c.Update();
        const TString filename = TString(figureFolder)+"/"+hist->GetName();
        c.SaveAs(filename+".png");
        c.SaveAs(filename+".pdf");
        ++count;
    }
    std::cout << "Saved " << count << " PHAST histograms in " << figureFolder << std::endl;
}
