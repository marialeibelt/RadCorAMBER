#include <iostream>
#include <vector>
#include <cmath>
#include <filesystem>
#include <memory>
#include <string>
#include <algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"

using namespace std;

constexpr double MUON_MASS = 0.105658375; // GeV


// ============================================================
// Result container for one sample
// ============================================================

struct Q2Result {

    unique_ptr<TH1D> hExact;
    unique_ptr<TH1D> hApprox;

    unique_ptr<TH1D> hRelativeDifference;
    unique_ptr<TH1D> hAbsoluteDifference;

    unique_ptr<TH2D> hExactVsApprox;
    unique_ptr<TH2D> hRelativeDifferenceVsTheta;
    unique_ptr<TH2D> hRelativeDifferenceVsEnergy;

    long long analyzedMuons = 0;
    long long invalidCases = 0;

    double sumAbsRelativeDifference = 0.0;
    double maxAbsRelativeDifference = 0.0;
};


// ============================================================
// Plot helpers
// ============================================================

void save2D(
    TH2D* hist,
    const string& outputFile,
    bool logX = false,
    bool logY = false,
    bool logZ = true,
    bool diagonal = false)
{
    TCanvas canvas("canvas2d", hist->GetTitle(), 850, 700);
    canvas.SetRightMargin(0.15);

    if (logX) canvas.SetLogx();
    if (logY) canvas.SetLogy();
    if (logZ) canvas.SetLogz();

    hist->Draw("COLZ");

    if (diagonal) {

        const double xmin = hist->GetXaxis()->GetXmin();
        const double xmax = hist->GetXaxis()->GetXmax();

        TLine line(xmin, xmin, xmax, xmax);
        line.SetLineStyle(2);
        line.SetLineWidth(2);
        line.Draw("SAME");
    }

    canvas.SaveAs(outputFile.c_str());

    cout << "Created: " << outputFile << endl;
}


void save1D(
    TH1D* hist,
    const string& outputFile,
    bool logY = true)
{
    TCanvas canvas("canvas1d", hist->GetTitle(), 850, 650);

    if (logY)
        canvas.SetLogy();

    hist->Draw("HIST");

    canvas.SaveAs(outputFile.c_str());

    cout << "Created: " << outputFile << endl;
}


// Exact and approximation for ONE sample on the same plot.
void saveExactApproxComparison(
    TH1D* hExact,
    TH1D* hApprox,
    const string& title,
    const string& outputFile,
    bool logY = true)
{
    unique_ptr<TH1D> exact(
        static_cast<TH1D*>(hExact->Clone("exact_plot"))
    );

    unique_ptr<TH1D> approx(
        static_cast<TH1D*>(hApprox->Clone("approx_plot"))
    );

    exact->SetDirectory(nullptr);
    approx->SetDirectory(nullptr);

    exact->SetLineColor(kBlue + 1);
    exact->SetLineWidth(2);

    approx->SetLineColor(kRed + 1);
    approx->SetLineWidth(2);
    approx->SetLineStyle(2);

    exact->SetTitle(title.c_str());
    exact->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
    exact->GetYaxis()->SetTitle("Entries");

    const double ymax =
        max(exact->GetMaximum(), approx->GetMaximum());

    if (ymax > 0.0)
        exact->SetMaximum(1.25 * ymax);

    TCanvas canvas("canvas_exact_approx", title.c_str(), 850, 650);

    if (logY)
        canvas.SetLogy();

    exact->Draw("HIST");
    approx->Draw("HIST SAME");

    TLegend legend(0.67, 0.76, 0.88, 0.88);
    legend.AddEntry(exact.get(),  "Exact Q^{2}", "l");
    legend.AddEntry(approx.get(), "Ultrarel. Q^{2}", "l");
    legend.Draw();

    canvas.SaveAs(outputFile.c_str());

    cout << "Created: " << outputFile << endl;
}


// LO and NLO, exact and approximation, all in one plot.
// LO and NLO are normalized individually because they can contain
// different numbers of events.
void saveAllQ2Comparison(
    TH1D* hExactLO,
    TH1D* hApproxLO,
    TH1D* hExactNLO,
    TH1D* hApproxNLO,
    const string& outputFile)
{
    unique_ptr<TH1D> exactLO(
        static_cast<TH1D*>(hExactLO->Clone("exactLO_plot"))
    );

    unique_ptr<TH1D> approxLO(
        static_cast<TH1D*>(hApproxLO->Clone("approxLO_plot"))
    );

    unique_ptr<TH1D> exactNLO(
        static_cast<TH1D*>(hExactNLO->Clone("exactNLO_plot"))
    );

    unique_ptr<TH1D> approxNLO(
        static_cast<TH1D*>(hApproxNLO->Clone("approxNLO_plot"))
    );

    exactLO->SetDirectory(nullptr);
    approxLO->SetDirectory(nullptr);
    exactNLO->SetDirectory(nullptr);
    approxNLO->SetDirectory(nullptr);

    // Normalize LO pair with the same factor
    if (exactLO->Integral() > 0.0) {

        const double scaleLO =
            1.0 / exactLO->Integral();

        exactLO->Scale(scaleLO);
        approxLO->Scale(scaleLO);
    }

    // Normalize NLO pair with the same factor
    if (exactNLO->Integral() > 0.0) {

        const double scaleNLO =
            1.0 / exactNLO->Integral();

        exactNLO->Scale(scaleNLO);
        approxNLO->Scale(scaleNLO);
    }

    exactLO->SetLineColor(kBlue + 1);
    exactLO->SetLineWidth(2);

    approxLO->SetLineColor(kBlue + 1);
    approxLO->SetLineWidth(2);
    approxLO->SetLineStyle(2);

    exactNLO->SetLineColor(kRed + 1);
    exactNLO->SetLineWidth(2);

    approxNLO->SetLineColor(kRed + 1);
    approxNLO->SetLineWidth(2);
    approxNLO->SetLineStyle(2);

    exactLO->SetTitle(
        "Q^{2}: exact and ultrarelativistic approximation"
    );

    exactLO->GetXaxis()->SetTitle(
        "Q^{2} [GeV^{2}]"
    );

    exactLO->GetYaxis()->SetTitle(
        "Normalized entries"
    );

    const double ymax =
        max(
            max(exactLO->GetMaximum(), approxLO->GetMaximum()),
            max(exactNLO->GetMaximum(), approxNLO->GetMaximum())
        );

    if (ymax > 0.0)
        exactLO->SetMaximum(1.25 * ymax);

    TCanvas canvas(
        "canvas_all_q2",
        "Q2 exact and approximation",
        850,
        650
    );

    canvas.SetLogy();

    exactLO->Draw("HIST");
    approxLO->Draw("HIST SAME");
    exactNLO->Draw("HIST SAME");
    approxNLO->Draw("HIST SAME");

    TLegend legend(0.62, 0.68, 0.88, 0.88);

    legend.AddEntry(
        exactLO.get(),
        "LO exact",
        "l"
    );

    legend.AddEntry(
        approxLO.get(),
        "LO ultrarel.",
        "l"
    );

    legend.AddEntry(
        exactNLO.get(),
        "NLO exact",
        "l"
    );

    legend.AddEntry(
        approxNLO.get(),
        "NLO ultrarel.",
        "l"
    );

    legend.Draw();

    canvas.SaveAs(
        outputFile.c_str()
    );

    cout
        << "Created: "
        << outputFile
        << endl;
}


// ============================================================
// Analyze one ROOT file
// ============================================================

Q2Result analyzeFile(
    const string& inputFile,
    const string& label)
{
    Q2Result result;


    // --------------------------------------------------------
    // Histograms
    // --------------------------------------------------------

    result.hExact =
        make_unique<TH1D>(
            ("hQ2Exact_" + label).c_str(),
            ("Exact Q^{2} - " + label).c_str(),
            500,
            0.0,
            0.05
        );


    result.hApprox =
        make_unique<TH1D>(
            ("hQ2Approx_" + label).c_str(),
            ("Ultrarelativistic Q^{2} - " + label).c_str(),
            500,
            0.0,
            0.05
        );


    // Your observed values are around 1.12e-6.
    result.hRelativeDifference =
        make_unique<TH1D>(
            ("hRelativeDifference_" + label).c_str(),
            ("Relative Q^{2} difference - " + label
             + ";(Q^{2}_{UR}-Q^{2}_{exact})/Q^{2}_{exact};Events").c_str(),
            400,
            1.0e-6,
            1.2e-6
        );


    result.hAbsoluteDifference =
        make_unique<TH1D>(
            ("hAbsoluteDifference_" + label).c_str(),
            ("Absolute Q^{2} difference - " + label
             + ";Q^{2}_{UR}-Q^{2}_{exact} [GeV^{2}];Events").c_str(),
            400,
            0.0,
            6.0e-8
        );


    result.hExactVsApprox =
        make_unique<TH2D>(
            ("hQ2ExactVsApprox_" + label).c_str(),
            ("Exact vs ultrarelativistic Q^{2} - " + label
             + ";Q^{2}_{exact} [GeV^{2}];Q^{2}_{UR} [GeV^{2}]").c_str(),
            300,
            0.0,
            0.05,
            300,
            0.0,
            0.05
        );


    result.hRelativeDifferenceVsTheta =
        make_unique<TH2D>(
            ("hRelativeDifferenceVsTheta_" + label).c_str(),
            ("Relative Q^{2} difference vs muon angle - " + label
             + ";#theta_{#mu} [mrad];"
               "(Q^{2}_{UR}-Q^{2}_{exact})/Q^{2}_{exact}").c_str(),
            250,
            0.0,
            2.5,
            250,
            1.0e-6,
            1.2e-6
        );


    result.hRelativeDifferenceVsEnergy =
        make_unique<TH2D>(
            ("hRelativeDifferenceVsEnergy_" + label).c_str(),
            ("Relative Q^{2} difference vs scattered muon energy - " + label
             + ";E'_{#mu} [GeV];"
               "(Q^{2}_{UR}-Q^{2}_{exact})/Q^{2}_{exact}").c_str(),
            250,
            0.0,
            105.0,
            250,
            1.0e-6,
            1.2e-6
        );


    // --------------------------------------------------------
    // Open ROOT file
    // --------------------------------------------------------

    TFile* file =
        TFile::Open(
            inputFile.c_str()
        );


    if (!file || file->IsZombie()) {

        cerr
            << "ERROR: could not open "
            << inputFile
            << endl;

        return result;
    }


    TTree* tree =
        static_cast<TTree*>(
            file->Get("Output")
        );


    if (!tree) {

        cerr
            << "ERROR: TTree 'Output' not found in "
            << inputFile
            << endl;

        file->Close();

        return result;
    }


    // --------------------------------------------------------
    // Branches
    // --------------------------------------------------------

    double beamEnergy = 0.0;

    double beamMomentumX = 0.0;
    double beamMomentumY = 0.0;
    double beamMomentumZ = 0.0;


    vector<int>* scatteredPID =
        nullptr;

    vector<double>* scatteredEnergy =
        nullptr;

    vector<double>* scatteredMomentumX =
        nullptr;

    vector<double>* scatteredMomentumY =
        nullptr;

    vector<double>* scatteredMomentumZ =
        nullptr;


    tree->SetBranchAddress(
        "beamEnergy",
        &beamEnergy
    );

    tree->SetBranchAddress(
        "beamMomentumX",
        &beamMomentumX
    );

    tree->SetBranchAddress(
        "beamMomentumY",
        &beamMomentumY
    );

    tree->SetBranchAddress(
        "beamMomentumZ",
        &beamMomentumZ
    );


    tree->SetBranchAddress(
        "scatteredPID",
        &scatteredPID
    );

    tree->SetBranchAddress(
        "scatteredEnergy",
        &scatteredEnergy
    );

    tree->SetBranchAddress(
        "scatteredMomentumX",
        &scatteredMomentumX
    );

    tree->SetBranchAddress(
        "scatteredMomentumY",
        &scatteredMomentumY
    );

    tree->SetBranchAddress(
        "scatteredMomentumZ",
        &scatteredMomentumZ
    );


    const Long64_t nEvents =
        tree->GetEntries();


    cout << endl;
    cout
        << label
        << " input: "
        << inputFile
        << endl;

    cout
        << label
        << " events: "
        << nEvents
        << endl;


    // ========================================================
    // Event loop
    // ========================================================

    for (Long64_t i = 0;
         i < nEvents;
         ++i) {


        tree->GetEntry(i);


        if (!scatteredPID ||
            !scatteredEnergy ||
            !scatteredMomentumX ||
            !scatteredMomentumY ||
            !scatteredMomentumZ) {

            ++result.invalidCases;
            continue;
        }


        const double pMuon =
            sqrt(
                beamMomentumX * beamMomentumX
                +
                beamMomentumY * beamMomentumY
                +
                beamMomentumZ * beamMomentumZ
            );


        if (pMuon <= 0.0) {

            ++result.invalidCases;
            continue;
        }


        for (size_t j = 0;
             j < scatteredPID->size();
             ++j) {


            if (scatteredPID->at(j) != 13)
                continue;


            const double EPrime =
                scatteredEnergy->at(j);


            const double pxPrime =
                scatteredMomentumX->at(j);

            const double pyPrime =
                scatteredMomentumY->at(j);

            const double pzPrime =
                scatteredMomentumZ->at(j);


            const double pPrime =
                sqrt(
                    pxPrime * pxPrime
                    +
                    pyPrime * pyPrime
                    +
                    pzPrime * pzPrime
                );


            if (pPrime <= 0.0) {

                ++result.invalidCases;
                continue;
            }


            // Full angle relative to the actual beam direction.
            double cosTheta =
                (
                    beamMomentumX * pxPrime
                    +
                    beamMomentumY * pyPrime
                    +
                    beamMomentumZ * pzPrime
                )
                /
                (pMuon * pPrime);


            cosTheta =
                max(
                    -1.0,
                    min(1.0, cosTheta)
                );


            const double theta =
                acos(cosTheta);


            const double thetaMrad =
                theta * 1000.0;


            // =================================================
            // Exact
            //
            // Q^2 =
            // 2 ( E_mu E'_mu
            //     - p_mu p'_mu cos(theta_mu)
            //     - m_mu^2 )
            // =================================================

            const double Q2Exact =
                2.0
                *
                (
                    beamEnergy * EPrime
                    -
                    pMuon * pPrime * cosTheta
                    -
                    MUON_MASS * MUON_MASS
                );


            // =================================================
            // Ultrarelativistic approximation
            //
            // Q^2 =
            // 4 E_mu E'_mu sin^2(theta_mu / 2)
            // =================================================

            const double sinHalfTheta =
                sin(theta / 2.0);


            const double Q2Approx =
                4.0
                *
                beamEnergy
                *
                EPrime
                *
                sinHalfTheta
                *
                sinHalfTheta;


            if (Q2Exact <= 0.0) {

                ++result.invalidCases;
                continue;
            }


            const double absoluteDifference =
                Q2Approx - Q2Exact;


            const double relativeDifference =
                absoluteDifference / Q2Exact;


            result.hExact->Fill(
                Q2Exact
            );


            result.hApprox->Fill(
                Q2Approx
            );


            result.hExactVsApprox->Fill(
                Q2Exact,
                Q2Approx
            );


            result.hAbsoluteDifference->Fill(
                absoluteDifference
            );


            result.hRelativeDifference->Fill(
                relativeDifference
            );


            result.hRelativeDifferenceVsTheta->Fill(
                thetaMrad,
                relativeDifference
            );


            result.hRelativeDifferenceVsEnergy->Fill(
                EPrime,
                relativeDifference
            );


            const double absRelativeDifference =
                abs(relativeDifference);


            result.sumAbsRelativeDifference +=
                absRelativeDifference;


            result.maxAbsRelativeDifference =
                max(
                    result.maxAbsRelativeDifference,
                    absRelativeDifference
                );


            ++result.analyzedMuons;
        }
    }


    file->Close();

    return result;
}


// ============================================================
// Summary helper
// ============================================================

void printSummary(
    const string& label,
    const Q2Result& result)
{
    cout << endl;
    cout << "----------------------------------------------" << endl;
    cout << label << endl;
    cout << "----------------------------------------------" << endl;

    cout
        << "Analyzed muons:                  "
        << result.analyzedMuons
        << endl;

    cout
        << "Invalid/skipped cases:           "
        << result.invalidCases
        << endl;


    if (result.analyzedMuons > 0) {

        const double meanAbsRelativeDifference =
            result.sumAbsRelativeDifference
            /
            static_cast<double>(
                result.analyzedMuons
            );


        cout
            << "Mean |relative difference|:      "
            << meanAbsRelativeDifference
            << endl;


        cout
            << "Maximum |relative difference|:   "
            << result.maxAbsRelativeDifference
            << endl;
    }
}


// ============================================================
// Main
// ============================================================

int main(int argc, char* argv[])
{
    if (argc != 4) {

        cerr
            << "Usage:\n"
            << "  "
            << argv[0]
            << " <LO_input.root>"
            << " <NLO_input.root>"
            << " <output_folder>\n";

        return 1;
    }


    const string loInputFile =
        argv[1];

    const string nloInputFile =
        argv[2];

    const string outputFolder =
        argv[3];


    const string loFolder =
        outputFolder + "/LO";

    const string nloFolder =
        outputFolder + "/NLO";

    const string comparisonFolder =
        outputFolder + "/comparison";


    filesystem::create_directories(
        loFolder
    );

    filesystem::create_directories(
        nloFolder
    );

    filesystem::create_directories(
        comparisonFolder
    );


    gStyle->SetOptStat(0);


    // ========================================================
    // Analyze LO and NLO
    // ========================================================

    Q2Result lo =
        analyzeFile(
            loInputFile,
            "LO"
        );


    Q2Result nlo =
        analyzeFile(
            nloInputFile,
            "NLO"
        );


    // ========================================================
    // LO plots
    // ========================================================

    saveExactApproxComparison(
        lo.hExact.get(),
        lo.hApprox.get(),
        "LO: exact vs ultrarelativistic Q^{2}",
        loFolder + "/Q2_exact_and_ultrarel.pdf",
        true
    );


    save2D(
        lo.hExactVsApprox.get(),
        loFolder + "/Q2_exact_vs_ultrarel.pdf",
        false,
        false,
        true,
        true
    );


    save1D(
        lo.hRelativeDifference.get(),
        loFolder + "/Q2_relative_difference.pdf",
        true
    );


    save1D(
        lo.hAbsoluteDifference.get(),
        loFolder + "/Q2_absolute_difference.pdf",
        true
    );


    save2D(
        lo.hRelativeDifferenceVsTheta.get(),
        loFolder + "/Q2_relative_difference_vs_theta.pdf",
        false,
        false,
        true,
        false
    );


    save2D(
        lo.hRelativeDifferenceVsEnergy.get(),
        loFolder + "/Q2_relative_difference_vs_energy.pdf",
        false,
        false,
        true,
        false
    );


    // ========================================================
    // NLO plots
    // ========================================================

    saveExactApproxComparison(
        nlo.hExact.get(),
        nlo.hApprox.get(),
        "NLO: exact vs ultrarelativistic Q^{2}",
        nloFolder + "/Q2_exact_and_ultrarel.pdf",
        true
    );


    save2D(
        nlo.hExactVsApprox.get(),
        nloFolder + "/Q2_exact_vs_ultrarel.pdf",
        false,
        false,
        true,
        true
    );


    save1D(
        nlo.hRelativeDifference.get(),
        nloFolder + "/Q2_relative_difference.pdf",
        true
    );


    save1D(
        nlo.hAbsoluteDifference.get(),
        nloFolder + "/Q2_absolute_difference.pdf",
        true
    );


    save2D(
        nlo.hRelativeDifferenceVsTheta.get(),
        nloFolder + "/Q2_relative_difference_vs_theta.pdf",
        false,
        false,
        true,
        false
    );


    save2D(
        nlo.hRelativeDifferenceVsEnergy.get(),
        nloFolder + "/Q2_relative_difference_vs_energy.pdf",
        false,
        false,
        true,
        false
    );


    // ========================================================
    // Combined LO/NLO plot
    // ========================================================

    saveAllQ2Comparison(
        lo.hExact.get(),
        lo.hApprox.get(),
        nlo.hExact.get(),
        nlo.hApprox.get(),
        comparisonFolder
            + "/Q2_LO_NLO_exact_and_ultrarel.pdf"
    );


    // ========================================================
    // Save all histograms
    // ========================================================

    const string outputRootFile =
        outputFolder
        + "/Q2approx_test_histograms.root";


    TFile output(
        outputRootFile.c_str(),
        "RECREATE"
    );


    output.mkdir("LO");
    output.cd("LO");

    lo.hExact->Write();
    lo.hApprox->Write();
    lo.hRelativeDifference->Write();
    lo.hAbsoluteDifference->Write();
    lo.hExactVsApprox->Write();
    lo.hRelativeDifferenceVsTheta->Write();
    lo.hRelativeDifferenceVsEnergy->Write();


    output.cd();
    output.mkdir("NLO");
    output.cd("NLO");

    nlo.hExact->Write();
    nlo.hApprox->Write();
    nlo.hRelativeDifference->Write();
    nlo.hAbsoluteDifference->Write();
    nlo.hExactVsApprox->Write();
    nlo.hRelativeDifferenceVsTheta->Write();
    nlo.hRelativeDifferenceVsEnergy->Write();


    output.Close();


    // ========================================================
    // Summary
    // ========================================================

    cout << endl;
    cout << "==============================================" << endl;
    cout << "Q2 approximation test" << endl;
    cout << "==============================================" << endl;


    printSummary(
        "LO",
        lo
    );


    printSummary(
        "NLO",
        nlo
    );


    cout << endl;

    cout
        << "Histograms written to:           "
        << outputRootFile
        << endl;

    cout << "==============================================" << endl;


    return 0;
}
