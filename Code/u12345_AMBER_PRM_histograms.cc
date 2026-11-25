// Reconstructed distributions directly from CORAL mDST using PHAST.
// A PHAST UserEvent module (not a standalone ROOT executable).
// No detector energy loss, event-by-event muon-defined Q2, or TGEANT truth matching here.
// Candidate definition: uniquely reconstructed negative outgoing track with p>70 GeV/c.
// This is a PRM-motivated kinematic candidate, NOT validated muon PID.
// With real photons, muon-side and proton-side Q2 can differ.

#include "PaEvent.h"
#include "Phast.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include <cmath>
#include <vector>

namespace {
constexpr double muMass = 0.105658375;  // GeV
constexpr double nominalBeamE = 100.0;  // GeV, only for approximate Q2

struct Histograms {
    TH1D *cutflow, *nTracks, *nVertices, *nOut, *nMuons;
    TH1D *vx, *vy, *vz, *trackP, *trackTheta, *trackXX0, *zFirst, *zLast;
    TH1D *muP, *muE, *muTheta, *muThetaFine, *muPhi, *q2Approx, *q2ApproxFine, *nuApprox;
    TH2D *vxy, *eVsTheta, *q2VsTheta;

    void book() {
        Phast::Ref().HistFileDir("AMBER_PRM");
        cutflow = new TH1D("h_cutflow", "Selection;Step;Events", 4, 0, 4);
        cutflow->GetXaxis()->SetBinLabel(1, "All events");
        cutflow->GetXaxis()->SetBinLabel(2, "Primary vertex");
        cutflow->GetXaxis()->SetBinLabel(3, "Outgoing particle");
        cutflow->GetXaxis()->SetBinLabel(4, "1 negative track, p>70");
        nTracks = new TH1D("h_n_tracks", "Reconstructed tracks;N_{tracks};Events", 80,-0.5,79.5);
        nVertices = new TH1D("h_n_vertices", "Reconstructed vertices;N_{vertices};Events", 20,-0.5,19.5);
        nOut = new TH1D("h_n_outgoing", "Outgoing particles in best vertex;N_{out};Events",20,-0.5,19.5);
        nMuons = new TH1D("h_n_muon_candidates", "High-p negative outgoing track candidates;N;Events",8,-0.5,7.5);
        vx = new TH1D("h_vertex_x", "Best vertex x;x [cm];Events",200,-20,20);
        vy = new TH1D("h_vertex_y", "Best vertex y;y [cm];Events",200,-20,20);
        vz = new TH1D("h_vertex_z", "Best vertex z;z [cm];Events",300,-500,500);
        vxy = new TH2D("h_vertex_xy", "Vertex x/y;x [cm];y [cm]",120,-12,12,120,-12,12);
        trackP = new TH1D("h_track_p", "Outgoing tracks;Momentum [GeV/c];Tracks",240,0,120);
        trackTheta = new TH1D("h_track_theta", "Outgoing tracks;Polar angle [mrad];Tracks",200,0,100);
        trackXX0 = new TH1D("h_track_xx0", "Track material budget;X/X_{0};Tracks",200,0,100);
        zFirst = new TH1D("h_track_z_first", "Track first hit;z_{first} [cm];Tracks",260,-500,8000);
        zLast = new TH1D("h_track_z_last", "Track last hit;z_{last} [cm];Tracks",260,-500,8000);
        muP = new TH1D("h_muon_p", "Scattered #mu momentum;p_{#mu'} [GeV/c];Events",240,0,120);
        muE = new TH1D("h_muon_E", "Scattered #mu energy;E_{#mu'} [GeV];Events",240,0,120);
        muTheta = new TH1D("h_muon_theta", "Scattered #mu angle to +z;#theta_{#mu'} [mrad];Events",240,0,12);
        muThetaFine = new TH1D("h_muon_theta_fine", "Scattered #mu angle (zoom);#theta_{#mu'} [mrad];Events",250,0,2);
        muPhi = new TH1D("h_muon_phi", "Scattered #mu azimuth;#phi_{#mu'} [rad];Events",128,-TMath::Pi(),TMath::Pi());
        q2Approx = new TH1D("h_q2_nominal_beam_approx", "Q^{2}, nominal 100 GeV +z beam (approx);Q^{2}_{approx} [GeV^{2}];Events",250,0,0.05);
        q2ApproxFine = new TH1D("h_q2_nominal_beam_approx_fine", "Q^{2} nominal beam (zoom);Q^{2}_{approx} [GeV^{2}];Events",250,0,0.001);
        nuApprox = new TH1D("h_nu_nominal_beam_approx", "Approx energy transfer (NOT detector energy loss);E_{beam}^{nom}-E_{#mu'} [GeV];Events",220,-10,100);
        eVsTheta = new TH2D("h_muon_E_vs_theta", "Reconstructed scattered #mu;#theta [mrad];E_{#mu'} [GeV]",120,0,12,160,0,120);
        q2VsTheta = new TH2D("h_q2_vs_theta_nominal_beam_approx", "Q^{2}_{approx} vs angle;#theta [mrad];Q^{2}_{approx} [GeV^{2}]",120,0,12,120,0,0.05);
    }
};
} // namespace

// Change UserEvent slot to match your local PHAST setup if needed.
void UserEvent12345(PaEvent& e) {
    static Histograms h;
    static bool first = true;
    if (first) { h.book(); first = false; }

    h.cutflow->Fill(0.5);
    h.nTracks->Fill(e.NTrack());
    h.nVertices->Fill(e.NVertex());
    const int iv = e.iBestPrimaryVertex();
    if (iv < 0) return;
    h.cutflow->Fill(1.5);

    const PaVertex& v = e.vVertex(iv);
    h.vx->Fill(v.X());
    h.vy->Fill(v.Y());
    h.vz->Fill(v.Z());
    h.vxy->Fill(v.X(),v.Y());
    h.nOut->Fill(v.NOutParticles());
    if (v.NOutParticles() == 0) return;
    h.cutflow->Fill(2.5);

    std::vector<int> muonCandidates;
    for (int i=0; i<v.NOutParticles(); ++i) {
        int ip = v.iOutParticle(i);
        const PaParticle& p = e.vParticle(ip);
        const PaTPar& par = p.ParInVtx(iv);
        // IsMuPrim() is NOT used: it does not establish scattered-muon PID.
        // For PRM, use a reproducible kinematic candidate definition.
        if (p.Q() == -1 && par.Mom() > 70.0) muonCandidates.push_back(ip);
        h.trackP->Fill(par.Mom());
        const auto momentum = par.LzVec(0.0).Vect();
        if (momentum.Mag2()>0) h.trackTheta->Fill(momentum.Theta()*1000.0);
        if (p.iTrack() >= 0) {
            const PaTrack& t = e.vTrack(p.iTrack());
            h.trackXX0->Fill(t.XX0()); // radiation lengths; NOT dE/dx
            h.zFirst->Fill(t.ZFirst());
            h.zLast->Fill(t.ZLast());
        }
    }
    h.nMuons->Fill(muonCandidates.size());
    if (muonCandidates.size()!=1) return;
    h.cutflow->Fill(3.5);

    const PaParticle& mu = e.vParticle(muonCandidates.front());
    const TLorentzVector pOut = mu.ParInVtx(iv).LzVec(muMass);
    const double theta = pOut.Theta()*1000.0;
    h.muP->Fill(pOut.P());
    h.muE->Fill(pOut.E());
    h.muTheta->Fill(theta);
    h.muThetaFine->Fill(theta);
    h.muPhi->Fill(pOut.Phi());
    h.eVsTheta->Fill(theta,pOut.E());

    // Diagnostic only! Use reconstructed incoming muon for physics-level Q2.
    const double pBeam = std::sqrt(nominalBeamE*nominalBeamE-muMass*muMass);
    const TLorentzVector pIn(0.0,0.0,pBeam,nominalBeamE);
    const double q2 = -(pIn-pOut).M2();
    h.q2Approx->Fill(q2);
    h.q2ApproxFine->Fill(q2);
    h.nuApprox->Fill(nominalBeamE-pOut.E());
    h.q2VsTheta->Fill(theta,q2);
}
