/// \file CheckTracks.C
/// \brief Simple macro to check ITSU tracks

#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <array>

#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TH2F.h>
#include <TNtuple.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TString.h>

#include "SimulationDataFormat/TrackReference.h"
#include "SimulationDataFormat/MCTrack.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "DataFormatsITSMFT/Cluster.h"
#include "DataFormatsITS/TrackITS.h"
#endif

void CheckTracks(std::string tracfile = "o2trac_its.root", std::string clusfile = "o2clus_its.root", std::string hitfile = "o2sim.root")
{
  
    using namespace o2::ITSMFT;
    using namespace o2::ITS;

    TFile* f = TFile::Open("CheckTracks_continuousRO.root", "recreate");

    TNtuple* nt = new TNtuple("ntt", "track ntuple",
                                //"mcYOut:recYOut:"
                                "mcZOut:recZOut:"
                                "mcPhiOut:recPhiOut:"
                                "mcThetaOut:recThetaOut:"
                                "mcPhi:recPhi:"
                                "mcLam:recLam:"
                                "mcPt:recPt:"
                                "ipD:ipZ:label");

    // MC tracks
    TFile* hitFile = TFile::Open(hitfile.data());
    TTree* mcTree = (TTree*)hitFile->Get("o2sim");
    std::vector<o2::MCTrack>* mcArr = nullptr;
    mcTree->SetBranchAddress("MCTrack", &mcArr);
    std::vector<o2::TrackReference>* mcTrackRefs = nullptr;
    mcTree->SetBranchAddress("TrackRefs", &mcTrackRefs);

    // Clusters
    TFile* clustFile = TFile::Open(clusfile.data());
    TTree* clusTree = (TTree*)clustFile->Get("o2sim");
    std::vector<Cluster>* clusArr = nullptr;
    clusTree->SetBranchAddress("ITSCluster", &clusArr);

    // Cluster MC labels
    o2::dataformats::MCTruthContainer<o2::MCCompLabel>* clusLabArr = nullptr;
    clusTree->SetBranchAddress("ITSClusterMCTruth", &clusLabArr);


    // Reconstructed tracks
    TFile* recTrackFile = TFile::Open(tracfile.data());
    TTree* recTree = (TTree*)recTrackFile->Get("o2sim");
    std::vector<TrackITS>* recArr = nullptr;
    recTree->SetBranchAddress("ITSTrack", &recArr);

    // Track MC labels
    o2::dataformats::MCTruthContainer<o2::MCCompLabel>* trkLabArr = nullptr;
    recTree->SetBranchAddress("ITSTrackMCTruth", &trkLabArr);


    // Get track Read Out Frame arrays
    //for clusters
    TTree* clusterROFTree = (TTree*)clustFile->Get("ITSClustersROF");
    std::vector<o2::ITSMFT::ROFRecord>* ROFRecordClusterArrray = nullptr;
    clusterROFTree->SetBranchAddress("ITSClustersROF", &ROFRecordClusterArrray);
    std::vector<o2::ITSMFT::ROFRecord>& ROFRecordClusterArrrayRef = *ROFRecordClusterArrray;
    clusterROFTree->GetEntry(0);

    TTree* clusterMC2ROFTree = (TTree*)clustFile->Get("ITSClustersMC2ROF");
    std::vector<o2::ITSMFT::MC2ROFRecord>* MC2ROFRecordClusterArrray = nullptr;
    clusterMC2ROFTree->SetBranchAddress("ITSClustersMC2ROF", &MC2ROFRecordClusterArrray);
    clusterMC2ROFTree->GetEntry(0);

    // for tracks
    TTree* trackROFTree = (TTree*)recTrackFile->Get("ITSTracksROF");
    std::vector<o2::ITSMFT::ROFRecord>* ROFRecordTrackArrray = nullptr;
    trackROFTree->SetBranchAddress("ITSTracksROF", &ROFRecordTrackArrray);
    std::vector<o2::ITSMFT::ROFRecord>& ROFRecordTrackArrrayRef = *ROFRecordTrackArrray;
    trackROFTree->GetEntry(0);

    TTree* trackMC2ROFTree = (TTree*)recTrackFile->Get("ITSTracksMC2ROF");
    std::vector<o2::ITSMFT::MC2ROFRecord>* MC2ROFRecordTrackArrray = nullptr;
    trackMC2ROFTree->SetBranchAddress("ITSTracksMC2ROF", &MC2ROFRecordTrackArrray);
    trackMC2ROFTree->GetEntry(0);

    unsigned int minROF = 0; 
    unsigned int maxROF = 0;

    unsigned int rofIndex = 0;
    unsigned int rofNEntries = 0;


    Int_t tf = 0, nrec = 0;
    Int_t lastEventID = -1;
    Int_t nev = mcTree->GetEntries();


    // Efficiency histograms
    // pT
    Int_t nb = 100;
    Double_t xbins[nb + 1], ptcutl = 0.01, ptcuth = 10.;
    Double_t a = TMath::Log(ptcuth / ptcutl) / nb;

    for (Int_t i = 0; i <= nb; i++) xbins[i] = ptcutl * TMath::Exp(i * a);
    
    TH1D* num = new TH1D("num", ";#it{p}_{T} (GeV/#it{c});Efficiency (fake-track rate)", nb, xbins);
    num->Sumw2();
    TH1D* fak = new TH1D("fak", ";#it{p}_{T} (GeV/#it{c});Fak", nb, xbins);
    fak->Sumw2();
    TH1D* den = new TH1D("den", ";#it{p}_{T} (GeV/#it{c});Den", nb, xbins);
    den->Sumw2();

    // Phi
    TH1D* hGoodVsPhi = new TH1D("hGoodVsPhi", ";#it{#phi};N_{Good tracks}", 180, -TMath::Pi(), TMath::Pi());
    hGoodVsPhi->Sumw2();
    TH1D* hFakeVsPhi = new TH1D("hFakeVsPhi", ";#it{#phi};N_{Fake tracks}", 180, -TMath::Pi(), TMath::Pi());
    hFakeVsPhi->Sumw2();
    TH1D* hMCvsPhi = new TH1D("hMCvsPhi", ";#it{#phi};N_{MC tracks}", 180, -TMath::Pi(), TMath::Pi());
    hMCvsPhi->Sumw2();

    // Lambda
    TH1D* hGoodVsLambda = new TH1D("hGoodVsLambda", ";#it{#lambda};N_{Good tracks}", 180, -TMath::Pi()/2, TMath::Pi()/2);
    hGoodVsLambda->Sumw2();
    TH1D* hFakeVsLambda = new TH1D("hFakeVsLambda", ";#it{#lambda};N_{Fake tracks}", 180, -TMath::Pi()/2, TMath::Pi()/2);
    hFakeVsLambda->Sumw2();
    TH1D* hMCvsLambda = new TH1D("hMCvsLambda", ";#it{#lambda};N_{MC tracks}", 180, -TMath::Pi()/2, TMath::Pi()/2);
    hMCvsLambda->Sumw2();




    for(int iEvent=0; iEvent<nev; iEvent++){

        Int_t nGen = 0, nGoo = 0, nFak = 0;

        mcTree->GetEvent(iEvent);

        Int_t nmc = mcArr->size();
        Int_t nmcrefs = mcTrackRefs->size();

        cout <<"################" <<endl;
        cout <<"n mc track : " <<nmc <<endl;

        while (nmc--) {

            const auto& mcTrack = (*mcArr)[nmc];

            Int_t mID = mcTrack.getMotherTrackId();
            
            if (mID >= 0) continue; // Select primary particles
          
            Int_t pdg = mcTrack.GetPdgCode();
            
            if (TMath::Abs(pdg) != 211) continue; // Select pions
          
            int ok = 0;


            // Check the availability of clusters

            // LOOP on :  MC2ROFRecord array
            for(auto& it : *MC2ROFRecordClusterArrray){        

                minROF = it.minROF;
                maxROF = it.maxROF;


                // LOOP on : ROFRecord array
                for(int iROF = minROF; iROF<=maxROF; iROF++){

                    auto eventIndex = ROFRecordClusterArrrayRef[iROF].getROFEntry();

                    int eventID = eventIndex.getEvent();

                    clusTree->GetEvent(eventID);

                    rofIndex = eventIndex.getIndex();
                    rofNEntries = ROFRecordClusterArrrayRef[iROF].getNROFEntries();


                    // LOOP on : clusters array
                    for(int iCluster = rofIndex; iCluster < rofIndex+rofNEntries; iCluster++){

                        const Cluster& c = (*clusArr)[iCluster];
                        auto lab = (clusLabArr->getLabels(iCluster))[0];

                        if (lab.getEventID() != iEvent) continue;
          
                        if (lab.getTrackID() != nmc) continue;
          
                        auto r = c.getX();

                        if (TMath::Abs(r - 2.2) < 0.5) 
                            ok |= 0b1;
          
                        if (TMath::Abs(r - 3.0) < 0.5)
                            ok |= 0b10;
          
                        if (TMath::Abs(r - 3.8) < 0.5)
                            ok |= 0b100;

                        if (TMath::Abs(r - 19.5) < 0.5)
                            ok |= 0b1000;

                        if (TMath::Abs(r - 24.5) < 0.5)
                            ok |= 0b10000;

                        if (TMath::Abs(r - 34.5) < 0.5)
                            ok |= 0b100000;

                        if (TMath::Abs(r - 39.5) < 0.5)
                            ok |= 0b1000000;


                        if (ok == 0b1111111) break;
                    }

                    if (ok == 0b1111111) break;

                }

                if (ok == 0b1111111) break;

            }

            if (ok != 0b1111111) continue;
        

            nGen++; // Generated tracks for the efficiency calculation

            // Float_t mcYOut=-1., recYOut=-1.;
            Float_t mcZOut = -1., recZOut = -1.;
            Float_t mcPhiOut = -1., recPhiOut = -1.;
            Float_t mcThetaOut = -1., recThetaOut = -1.;
            Float_t mcPx = mcTrack.GetStartVertexMomentumX();
            Float_t mcPy = mcTrack.GetStartVertexMomentumY();
            Float_t mcPz = mcTrack.GetStartVertexMomentumZ();
            Float_t mcPhi = TMath::ATan2(mcPy, mcPx), recPhi = -1.;
            Float_t mcPt = mcTrack.GetPt(), recPt = -1.;
            Float_t mcLam = TMath::ATan2(mcPz, mcPt), recLam = -1.;
            Float_t ip[2]{ 0., 0. };
            Float_t label = -123456789.;

            den->Fill(mcPt);


            bool recoTrackFound = false;


            // LOOP on :  MC2ROFRecord array
            for(auto& it : *MC2ROFRecordTrackArrray){        

                minROF = it.minROF;
                maxROF = it.maxROF;

                // LOOP on : ROFRecord array
                for(int iROF = minROF; iROF<maxROF+1; iROF++){

                    auto eventIndex = ROFRecordTrackArrrayRef[iROF].getROFEntry();

                    int eventID = eventIndex.getEvent();

                    clusTree->GetEvent(eventID);
                    recTree->GetEvent(eventID);

                    rofIndex = eventIndex.getIndex();
                    rofNEntries = ROFRecordTrackArrrayRef[iROF].getNROFEntries();

                    // LOOP on : tracks array
                    for(int iTrack = rofIndex; iTrack < rofIndex+rofNEntries; iTrack++){

                        const TrackITS& recTrack = (*recArr)[iTrack];

                        auto mclab = (trkLabArr->getLabels(iTrack))[0];
                
                        auto id = mclab.getEventID();

                        if (id != iEvent) continue;
          
                        Int_t lab = mclab.getTrackID();

                        if (TMath::Abs(lab) != nmc) continue;

                        for (auto& ref : *mcTrackRefs) {
                    
                            if (ref.getUserId() != 6) continue;
            
                            if (ref.getTrackID() != nmc) continue;
            
                            // mcYOut=ref.LocalY();
                            mcZOut = ref.Z();
                            mcPhiOut = ref.Phi();
                            mcThetaOut = ref.Theta();
                            break;
                        }

                        auto out = recTrack.getParamOut();
                        // recYOut = out.getY();
                        recZOut = out.getZ();
                        recPhiOut = out.getPhi();
                        recThetaOut = out.getTheta();

                        std::array<float, 3> p;
                        recTrack.getPxPyPzGlo(p);
                        recPt = recTrack.getPt();
                        recPhi = TMath::ATan2(p[1], p[0]);
                        recLam = TMath::ATan2(p[2], recPt);

                        Float_t vx = 0., vy = 0., vz = 0.; // Assumed primary vertex
                        Float_t bz = 5.;                   // Assumed magnetic field
                        recTrack.getImpactParams(vx, vy, vz, bz, ip);
                        label = lab;

                        if (label > 0) {
                            nGoo++; // Good found tracks for the efficiency calculation
                            num->Fill(mcPt);
                            recoTrackFound = true;
                            break;
                        } 
                        else {
                            nFak++; // Fake-track rate calculation
                            fak->Fill(mcPt);
                        }
                    }

            
                    nt->Fill( // mcYOut,recYOut,
                        mcZOut, recZOut, mcPhiOut, recPhiOut, mcThetaOut, recThetaOut, mcPhi, recPhi, mcLam, recLam, mcPt, recPt, ip[0],
                        ip[1], label);

                    if(recoTrackFound == true) break;
                }

                if(recoTrackFound == true) break;
            }

        }

        
        cout <<"nGen : " <<nGen <<endl;
        cout <<"nGoo : " <<nGoo <<endl;
        cout <<"nFake : " <<nFak <<endl;
        
        if (nGen > 0) {
            Float_t eff = nGoo / Float_t(nGen);
            Float_t rat = nFak / Float_t(nGen);
            std::cout << "Good found tracks: " << nGoo << ",  efficiency: " << eff << ",  fake-track rate: " << rat << std::endl;
        }

        cout <<"################" <<endl;
    }
    

    for(int i=0; i<fak->GetNbinsX(); i++){

        if(fak->GetBinContent(i+1) != 0) cout <<fak->GetBinContent(i+1) <<"     " <<den->GetBinContent(i+1) <<endl;
    }
    
  
    // "recPt>0" means "found tracks only"
    // "label>0" means "found good tracks only"
    new TCanvas;
    nt->Draw("ipD", "recPt>0 && label>0");
    new TCanvas;
    nt->Draw("mcLam-recLam", "recPt>0 && label>0");
    new TCanvas;
    nt->Draw("mcPt-recPt", "recPt>0 && label>0");
    new TCanvas;
    nt->Draw("mcZOut-recZOut", "recPt>0 && label>0 && abs(mcZOut-recZOut)<0.025");
    new TCanvas;
    nt->Draw("mcPhiOut-recPhiOut", "recPt>0 && label>0");
    new TCanvas;
    nt->Draw("mcThetaOut-recThetaOut", "recPt>0 && label>0");

    TCanvas* c1 = new TCanvas;
    c1->SetLogx();
    c1->SetGridx();
    c1->SetGridy();
    num->Divide(num, den, 1, 1, "b");
    num->Draw("histe");
    fak->Divide(fak, den, 1, 1, "b");
    fak->SetLineColor(2);
    fak->Draw("histesame");

    f->Write();
    f->Close();

    
}








