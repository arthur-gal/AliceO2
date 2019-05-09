/// \file CheckDigits.C
/// \brief Simple macro to check ITSU clusters

#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <TCanvas.h>
#include <TFile.h>
#include <TH2F.h>
#include <TNtuple.h>
#include <TString.h>
#include <TTree.h>

#include "ITSBase/GeometryTGeo.h"
#include "DataFormatsITSMFT/Cluster.h"
#include "ITSMFTSimulation/Hit.h"
#include "MathUtils/Cartesian3D.h"
#include "MathUtils/Utils.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#endif

void CheckClusters_continuousRO(std::string clusfile = "o2clus_its.root", std::string hitfile = "o2sim.root", std::string inputGeom = "O2geometry.root", std::string paramfile = "o2sim_par.root")
{

	using namespace o2::Base;
	using namespace o2::ITS;

  	using o2::ITSMFT::Cluster;
  	using o2::ITSMFT::Hit;

  	// Outputs file
  	TFile* f = TFile::Open("CheckClusters_continuousRO.root", "recreate");
  	TNtuple* nt = new TNtuple("ntc", "cluster ntuple", "x:y:z:dx:dz:lab:rof:ev:hlx:hlz:clx:clz");

  	// Map to store nTuple per topology
  	short nx = 0; // effective cluster size in X
  	short nz = 0; // effective cluster size in Z

  	std::map<short, TNtuple*> nTupleMAP;
  	bool findID = false;

  	short nxOffset = 8;
  	short nMask = 0xff;
  	short nTupleID = 0; // effective cluster size in X (1st byte) and Z (2nd byte) directions

  	// Geometry
  	o2::Base::GeometryManager::loadGeometry(inputGeom, "FAIRGeom");
  	auto gman = o2::ITS::GeometryTGeo::Instance();
  	gman->fillMatrixCache(o2::utils::bit2Mask(o2::TransformType::T2L, o2::TransformType::T2GRot,
												o2::TransformType::L2G)); // request cached transforms

  	// Hits
	TFile* hitFile = TFile::Open(hitfile.data());
	TTree* hitTree = (TTree*)hitFile->Get("o2sim");
	std::vector<o2::ITSMFT::Hit>* hitArray = nullptr;
	hitTree->SetBranchAddress("ITSHit", &hitArray);

  	// Clusters
  	TFile* clustFile = TFile::Open(clusfile.data());
  	TTree* clusTree = (TTree*)clustFile->Get("o2sim");
  	std::vector<Cluster>* clusArr = nullptr;
  	clusTree->SetBranchAddress("ITSCluster", &clusArr);

  	// Cluster MC labels
  	o2::dataformats::MCTruthContainer<o2::MCCompLabel>* clusLabArr = nullptr;
  	clusTree->SetBranchAddress("ITSClusterMCTruth", &clusLabArr);


  	// Get Read Out Frame arrays
  	TTree* ROFTree = (TTree*)clustFile->Get("ITSClustersROF");
  	std::vector<o2::ITSMFT::ROFRecord>* ROFRecordArrray = nullptr;
  	ROFTree->SetBranchAddress("ITSClustersROF", &ROFRecordArrray);
  	std::vector<o2::ITSMFT::ROFRecord>& ROFRecordArrrayRef = *ROFRecordArrray;
  	ROFTree->GetEntry(0);

  	TTree* MC2ROFTree = (TTree*)clustFile->Get("ITSClustersMC2ROF");
  	std::vector<o2::ITSMFT::MC2ROFRecord>* MC2ROFRecordArrray = nullptr;
  	MC2ROFTree->SetBranchAddress("ITSClustersMC2ROF", &MC2ROFRecordArrray);
  	MC2ROFTree->GetEntry(0);

  	unsigned int minROF = 0; 
  	unsigned int maxROF = 0;

  	unsigned int rofIndex = 0;
  	unsigned int rofNEntries = 0;

  	int lastReadHitEv = -1;


  	cout <<"TEST ttree" <<clusTree->GetEntries() <<endl;


  	// LOOP on :  MC2ROFRecord array
  	for(auto& it : *MC2ROFRecordArrray){	  	

	  	minROF = it.minROF;
	  	maxROF = it.maxROF;


	  	// LOOP on : ROFRecord array
	  	for(int iROF = minROF; iROF<=maxROF; iROF++){

			auto eventIndex = ROFRecordArrrayRef[iROF].getROFEntry();

			clusTree->GetEvent(eventIndex.getEvent());

			rofIndex = eventIndex.getIndex();
			rofNEntries = ROFRecordArrrayRef[iROF].getNROFEntries();


			// LOOP on : clusters array
			for(int iCluster = rofIndex; iCluster < rofIndex+rofNEntries; iCluster++){

	  			// cluster is in tracking coordinates always
				Cluster& c = (*clusArr)[iCluster];

	 			Int_t chipID = c.getSensorID();

	  			const auto locC = c.getXYZLoc(*gman);    // convert from tracking to local frame
	  			const auto gloC = c.getXYZGloRot(*gman); // convert from tracking to global frame

	  			auto lab = (clusLabArr->getLabels(iCluster))[0];

	  			float dx = 0, dz = 0;
	  			int trID = lab.getTrackID();
	  			int ievH = lab.getEventID();
	  
	  			Point3D<float> locH, locHsta;

	  			if (trID >= 0) { // is this cluster from hit or noise ?
		
					Hit* p = nullptr;
		
					if (lastReadHitEv != ievH) {
		  				
		  				hitTree->GetEvent(ievH);
		  				lastReadHitEv = ievH;
					}

					for (auto& ptmp : *hitArray) {
		  
		  				if (ptmp.GetDetectorID() != chipID) continue;
			
		  				if (ptmp.GetTrackID() != trID) continue;
		  
		  				p = &ptmp;
		  				break;
					}
		
		
					if (!p) {

		  				printf("did not find hit (scanned HitEvs %d %d) for cluster of tr%d on chip %d\n", ievH, it.eventRecordID, trID, chipID);
		  				locH.SetXYZ(0.f, 0.f, 0.f);
					} 
					else {
		  	
		  				// mean local position of the hit
		  				locH = gman->getMatrixL2G(chipID) ^ (p->GetPos()); // inverse conversion from global to local
		  				locHsta = gman->getMatrixL2G(chipID) ^ (p->GetPosStart());
		  				locH.SetXYZ(0.5 * (locH.X() + locHsta.X()), 0.5 * (locH.Y() + locHsta.Y()), 0.5 * (locH.Z() + locHsta.Z()));
		  				// std::cout << "chip "<< p->GetDetectorID() << "  PposGlo " << p->GetPos() << std::endl;
		  				// std::cout << "chip "<< c->getSensorID() << "  PposLoc " << locH << std::endl;
		  				dx = locH.X() - locC.X();
		  				dz = locH.Z() - locC.Z();

					}

					nx = c.getNx();
        			nz = c.getNz();

        			nTupleID = 0;
        			nTupleID = ((nx & nMask) << nxOffset) + (nz & nMask);
        			findID = false;

        			for(auto& mapIT : nTupleMAP){ 

          				if(nTupleID == mapIT.first){

            				findID = true;
            				mapIT.second->Fill(gloC.X(), gloC.Y(), gloC.Z(), dx, dz, trID, c.getROFrame(), it.eventRecordID, locH.X(), locH.Z(), locC.X(), locC.Z());
           					break;
          				}
        			}

        			if(findID == false){
          				nTupleMAP[nTupleID] = new TNtuple(Form("Nx_%i_Nz_%i", nx, nz), Form("Nx_%i_Nz_%i", nx, nz), "x:y:z:dx:dz:lab:rof:ev:hlx:hlz:clx:clz");
          				nTupleMAP[nTupleID]->Fill(gloC.X(), gloC.Y(), gloC.Z(), dx, dz, trID, c.getROFrame(), it.eventRecordID, locH.X(), locH.Z(), locC.X(), locC.Z());
        			}

				}

				nt->Fill(gloC.X(), gloC.Y(), gloC.Z(), dx, dz, trID, c.getROFrame(), it.eventRecordID, locH.X(), locH.Z(), locC.X(), locC.Z());
	  		}

		}

	}
 

  	new TCanvas;
  	nt->Draw("y:x");
  	new TCanvas;
  	nt->Draw("dx:dz");

  	TCanvas* cnv = nullptr;
 	TCanvas* cnvBis = nullptr;

  	f->cd();
  	nt->Write();


  	for(auto& mapIT : nTupleMAP){

    	int nx = (mapIT.first >> nxOffset) & nMask;
    	int nz = mapIT.first & nMask;

    	cnv = new TCanvas(Form("%i_%i_xy", nx, nz), "", 800, 600);
    	mapIT.second->Draw("x:y");

    	cnvBis = new TCanvas(Form("%i_%i_dxdz", nx, nz), "", 800, 600);
    	mapIT.second->Draw("dx:dz");

    	mapIT.second->Write();
    	cnv->Write();
    	cnvBis->Write();

    	delete cnv;
    	delete cnvBis;

    	cnv = nullptr;
    	cnvBis = nullptr;

  	}     

  	f->Close();











}
