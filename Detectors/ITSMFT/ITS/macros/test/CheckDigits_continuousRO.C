/// \file CheckDigits.C
/// \brief Simple macro to check ITSU digits

#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <TCanvas.h>
#include <TFile.h>
#include <TH2F.h>
#include <TNtuple.h>
#include <TString.h>
#include <TTree.h>

#include <vector>
#include "ITSBase/GeometryTGeo.h"
#include "ITSMFTBase/Digit.h"
#include "ITSMFTBase/SegmentationAlpide.h"
#include "ITSMFTSimulation/Hit.h"
#include "MathUtils/Utils.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "DetectorsBase/GeometryManager.h"

#include "CommonDataFormat/EvIndex.h"

#endif

using namespace o2::Base;


void CheckDigits_continuousRO(std::string digifile = "itsdigits.root", std::string hitfile = "o2sim.root", std::string inputGeom = "O2geometry.root", std::string paramfile = "o2sim_par.root"){


	using o2::ITSMFT::Digit;
	using o2::ITSMFT::Hit;
	using o2::ITSMFT::SegmentationAlpide;
	using namespace o2::ITS;

	TFile* f = TFile::Open("CheckDigits_continuousRO.root", "recreate");

	TNtuple* nt = new TNtuple("ntd", "digit ntuple", "id:x:y:z:rowD:colD:rowH:colH:xlH:zlH:xlcH:zlcH:dx:dz");

	// Geometry
	o2::Base::GeometryManager::loadGeometry(inputGeom, "FAIRGeom");
	auto* gman = o2::ITS::GeometryTGeo::Instance();
	gman->fillMatrixCache(o2::utils::bit2Mask(o2::TransformType::L2G));

	SegmentationAlpide seg;


	// Hits
	TFile* hitFile = TFile::Open(hitfile.data());
	TTree* hitTree = (TTree*)hitFile->Get("o2sim");
	std::vector<o2::ITSMFT::Hit>* hitArray = nullptr;
	hitTree->SetBranchAddress("ITSHit", &hitArray);

	// Digits
	TFile* digFile = TFile::Open(digifile.data());
	TTree* digTree = (TTree*)digFile->Get("o2sim");

	std::vector<o2::ITSMFT::Digit>* digArr = nullptr;
	digTree->SetBranchAddress("ITSDigit", &digArr);
	std::vector<o2::ITSMFT::Digit>& digArrRef = *digArr;

	o2::dataformats::MCTruthContainer<o2::MCCompLabel>* labels = nullptr;
	digTree->SetBranchAddress("ITSDigitMCTruth", &labels);

	int nevD = digTree->GetEntries(); // digits in cont. readout may be grouped as few events per entry
	int nevH = hitTree->GetEntries(); // hits are stored as one event per entry
	int lastReadHitEv = -1;

	int nDigitRead = 0, nDigitFilled = 0;


	// Get Read Out Frame arrays
	TTree* ROFTree = (TTree*)digFile->Get("ITSDigitROF");
	std::vector<o2::ITSMFT::ROFRecord>* ROFRecordArrray = nullptr;
	ROFTree->SetBranchAddress("ITSDigitROF", &ROFRecordArrray);
	std::vector<o2::ITSMFT::ROFRecord>& ROFRecordArrrayRef = *ROFRecordArrray;
	ROFTree->GetEntry(0);

	TTree* MC2ROFTree = (TTree*)digFile->Get("ITSDigitMC2ROF");
	std::vector<o2::ITSMFT::MC2ROFRecord>* MC2ROFRecordArrray = nullptr;
	MC2ROFTree->SetBranchAddress("ITSDigitMC2ROF", &MC2ROFRecordArrray);
	MC2ROFTree->GetEntry(0);


	unsigned int minROF = 0; // 1st ROFrame MC2ROFrame contribution
	unsigned int maxROF = 0; // Last ROFrame MC2ROFrame contribution

	unsigned int rofIndex = 0;
	unsigned int rofNEntries = 0;



  	// LOOP on :  MC2ROFRecord array
  	for(auto& it : *MC2ROFRecordArrray){	  	

	  	minROF = it.minROF;
	  	maxROF = it.maxROF;


	  	// LOOP on : ROFRecord array
	  	for(int iROF = minROF; iROF<=maxROF; iROF++){

			auto eventIndex = ROFRecordArrrayRef[iROF].getROFEntry();

			digTree->GetEvent(eventIndex.getEvent());

			rofIndex = eventIndex.getIndex();
			rofNEntries = ROFRecordArrrayRef[iROF].getNROFEntries();


			// LOOP on : digits array
			for(int iDigit = rofIndex; iDigit < rofIndex+rofNEntries; iDigit++){


		  		Int_t ix = digArrRef[iDigit].getRow(), iz = digArrRef[iDigit].getColumn();
		  		Float_t x = 0.f, z = 0.f;

		  		seg.detectorToLocal(ix, iz, x, z);

		  		const Point3D<float> locD(x, 0., z);

		  		Int_t chipID = digArrRef[iDigit].getChipIndex();

		  		const auto& labs = labels->getLabels(iDigit);

		  		int trID = labs[0].getTrackID();
		  		int ievH = labs[0].getEventID();

		  		if (trID >= 0) { // not a noise

					nDigitRead++;

					const auto gloD = gman->getMatrixL2G(chipID)(locD); // convert to global
					float dx = 0., dz = 0.;


					if (lastReadHitEv != ievH) {

			  			hitTree->GetEvent(ievH);
			  			lastReadHitEv = ievH;
					}

					bool ok = false;

					for (auto& p : *hitArray) {

			  			if (p.GetDetectorID() != chipID) continue;
					
			  			if (p.GetTrackID() != trID) continue;
				

			  			auto locH = gman->getMatrixL2G(chipID) ^ (p.GetPos()); // inverse conversion from global to local
			  			auto locHsta = gman->getMatrixL2G(chipID) ^ (p.GetPosStart());

			  			locH.SetXYZ(0.5 * (locH.X() + locHsta.X()), 0.5 * (locH.Y() + locHsta.Y()), 0.5 * (locH.Z() + locHsta.Z()));

			  			int row, col;
			  			float xlc, zlc;

			  			seg.localToDetector(locH.X(), locH.Z(), row, col);
			  			seg.detectorToLocal(row, col, xlc, zlc);
				
			  			nt->Fill(chipID, gloD.X(), gloD.Y(), gloD.Z(), ix, iz, row, col, locH.X(), locH.Z(), xlc, zlc,
					   			locH.X() - locD.X(), locH.Z() - locD.Z());

			  			ok = true;
			  			nDigitFilled++;
			  			break;

			  		} // end loop on hits array

			  		if (!ok) printf("did not find hit for digit %d in ev %d: MCEv:%d MCTrack %d\n", iDigit, it.eventRecordID, ievH, trID);
				}

		  	} // end loop on digits array

	  
		} // end loop on ROFRecords array


	}// end loop on MC2ROFRecords array


  	new TCanvas;
  	nt->Draw("y:x");
	new TCanvas;
	nt->Draw("dx:dz");
	
  	f->Write();
  	f->Close();
  	printf("read %d filled %d\n", nDigitRead, nDigitFilled);


}
