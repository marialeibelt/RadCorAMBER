// Compare identically named MC-Truth and PHAST-Reco histograms.
// ROOT: root -l -b -q 'plot_compare_tgeant_phast.C("mc.root","phast.root","plots")'
// Overlays are SHAPE comparisons only. They do NOT represent reconstruction
// efficiency or resolution without matched events and equivalent selections.
#include <algorithm>
#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TPad.h>
#include <TString.h>
#include <TSystem.h>
#include <TStyle.h>
#include <iostream>
#include <memory>
#include <vector>

void plot_compare_tgeant_phast(const char* mcfile,
                               const char* phastfile,
                               const char* outputdir) {
    std::unique_ptr<TFile> fMC(TFile::Open(mcfile,"READ"));
    std::unique_ptr<TFile> fReco(TFile::Open(phastfile,"READ"));
    if (!fMC || fMC->IsZombie() || !fReco || fReco->IsZombie()) {
        std::cerr << "ERROR: Unable to open one or both ROOT files\n";
        return;
    }
    TDirectory* mcDir = fMC->GetDirectory("MC_TRUTH");
    TDirectory* recoDir = fReco->GetDirectory("AMBER_PRM");
    if (!mcDir || !recoDir) {
        std::cerr << "ERROR: Expected MC_TRUTH in MC file, AMBER_PRM in PHAST file.\n";
        fMC->ls(); fReco->ls();
        return;
    }
    gSystem->mkdir(outputdir,true);
    gStyle->SetOptStat(0);
    // The files can have different event counts and selection efficiencies;
    // normalize *the plotted histograms* only, without altering source ROOT.
    const std::vector<const char*> oneD = {
        "h_vertex_x", "h_vertex_y", "h_vertex_z",
        "h_muon_p", "h_muon_E", "h_muon_theta", "h_muon_theta_fine", "h_muon_phi",
        "h_q2_nominal_beam_approx", "h_q2_nominal_beam_approx_fine",
        "h_nu_nominal_beam_approx"
    };
    int nplots = 0;
    for (const auto* name : oneD) {
        const auto* originalMC = dynamic_cast<TH1*>(mcDir->Get(name));
        const auto* originalReco = dynamic_cast<TH1*>(recoDir->Get(name));
        if (!originalMC || !originalReco) {
            std::cerr << "Skipping missing pair: " << name << "\n";
            continue;
        }
        if (originalMC->GetDimension()!=1 || originalReco->GetDimension()!=1) continue;
        auto hm = std::unique_ptr<TH1>(dynamic_cast<TH1*>(originalMC->Clone(TString(name)+"_MC_copy")));
        auto hr = std::unique_ptr<TH1>(dynamic_cast<TH1*>(originalReco->Clone(TString(name)+"_Reco_copy")));
        hm->SetDirectory(nullptr); hr->SetDirectory(nullptr);
        const double nm = hm->Integral(1,hm->GetNbinsX());
        const double nr = hr->Integral(1,hr->GetNbinsX());
        if (nm<=0 || nr<=0) {
            std::cerr << "Skipping empty/negative histogram pair: " << name << "\n";
            continue;
        }
        hm->Scale(1.0/nm); hr->Scale(1.0/nr);
        hm->SetLineColor(kBlue+1); hr->SetLineColor(kRed+1);
        hm->SetLineWidth(2); hr->SetLineWidth(2);
        const double ymax=std::max(hm->GetMaximum(),hr->GetMaximum());
        hm->SetMaximum(1.18*ymax);
        hm->SetTitle(TString::Format("%s; %s;Normalized entries",name,hm->GetXaxis()->GetTitle()));
        TCanvas c(TString::Format("c_%s",name),name,1000,700);
        c.SetLeftMargin(0.14); c.SetBottomMargin(0.13);
        hm->Draw("HIST"); hr->Draw("HIST SAME");
        TLegend legend(0.59,0.76,0.89,0.90);
        legend.SetBorderSize(0);
        legend.AddEntry(hm.get(),"TGEANT MC truth", "l");
        legend.AddEntry(hr.get(),"PHAST reconstructed", "l");
        legend.Draw();
        TString filename=TString(outputdir)+"/compare_"+name;
        c.SaveAs(filename+".png"); c.SaveAs(filename+".pdf");
        ++nplots;
    }
    const std::vector<const char*> twoD = {
        "h_vertex_xy", "h_muon_E_vs_theta", "h_q2_vs_theta_nominal_beam_approx"
    };
    for (const auto* name : twoD) {
        const auto* hm = dynamic_cast<TH2*>(mcDir->Get(name));
        const auto* hr = dynamic_cast<TH2*>(recoDir->Get(name));
        if (!hm || !hr) continue;
        TCanvas c(TString::Format("c2_%s",name),name,1600,650);
        c.Divide(2,1);
        c.cd(1); gPad->SetRightMargin(0.16); hm->Draw("COLZ");
        c.cd(2); gPad->SetRightMargin(0.16); hr->Draw("COLZ");
        TString filename=TString(outputdir)+"/compare_"+name;
        c.SaveAs(filename+".png"); c.SaveAs(filename+".pdf");
        ++nplots;
    }
    std::cout << "Wrote " << nplots << " MC-vs-reconstructed plots to " << outputdir << "\n";
    std::cout << "NOTE: generator/truth-level vs reconstructed-level shape comparison,"
              << " NOT a detector resolution or data/MC validation.\n";
}
