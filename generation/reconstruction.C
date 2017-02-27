/* * * 

Macro for reconstruction of jets in events using FastJet.
Its input are: Kinematics(...).root and galice(...).root (generated for example by generation.C)
containing kinematic information about each particle in event.

Writes to tree info about: 
	- jets themselves 
	- particles in jets
	- particles in event - enables to run again e.g. with other algorithm or parameter

It can be run using dedeicated macros: 
(multi)run_recon.C 
or for executing also reconstruction:
(multi)run_gener_recon.C 

Parameter 'fileNameParams' defines input files used as input.

* * */

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <Riostream.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TStopwatch.h>
#include <TDatime.h>
#include <TRandom.h>
#include <TDatabasePDG.h>
#include <TParticle.h>
#include <TArrayI.h>
#include <TTree.h>
#include <TFile.h>
#include <THnSparse.h>
#include "TCanvas.h"
#include "TROOT.h"
#include "TClassTable.h"
#include "TPDGCode.h"
#include <typeinfo>

#include "AliGenerator.h"
#include "AliPDG.h"
#include "AliRunLoader.h"
#include "AliRun.h"
#include "AliStack.h"
#include "AliHeader.h"
#include "AliGenPythia.h"	//### E-TCint::AutoLoadCallback: failure loading library libAliPythia6 for class AliGenPythia
#include "AliPythia.h"

//for jets:
#include "AliFastJetHeaderV1.h"
#include "AliFastJetFinder.h"
#include "AliAODEvent.h"
#include "AliMCParticle.h"
//#include "AliJetCalTrkTrackKine.h"


#endif


void reconstruction(Int_t evNumber=3000, const char* fileNameParams="") 
{

	Int_t nEvToPrint=0;


	TString* galiceFileName = new TString("galice");
	*galiceFileName += fileNameParams;
	*galiceFileName += ".root";


	// Connect the Root Galice file containing Geometry, Kine and Hits
	AliRunLoader* rl = AliRunLoader::Open(galiceFileName->Data());

	TDatabasePDG* DataBase = new TDatabasePDG();
		

	// main parameter for anti-kt algorithm
	Float_t jetRadius = 0.4;


	//
	// Jet finder settings go via the FastJetHeader
	//
	AliFastJetHeaderV1 *header = new AliFastJetHeaderV1;
	header->SetBGMode(0);
	// header->SetRadius(0.4);
	//Float_t jetRadius=0.2;//0.3//0.2//was 0.4//for Gustavo 0.3
	header->SetRparam(jetRadius); 
	//header->SetJetEtaMin(-0.9+jetRadius); // to be added the jet eta min
	//header->SetJetEtaMax(0.9-jetRadius); // to be added the jet eta max
	// no eta limits:
	header->SetJetEtaMin(-1e10);
	header->SetJetEtaMax(1e10); 
	header->SetPtMin(5.); //5. is default

	// set reconstruction algorithm
	header->SetAlgorithm(fastjet::antikt_algorithm);

	AliFastJetFinder *FastJet = new AliFastJetFinder;
	FastJet->SetJetHeader(header);



	// Set up dummy AOD event for JETAN output
	AliAODEvent *aod = new AliAODEvent();
	aod->CreateStdContent();
	FastJet->ConnectAOD(aod);


	// Set up input structures for FastJet
	AliJetCalTrkEvent JetFinderEvent(0,1);
	TClonesArray aliplist("AliMCParticle",1000);



	// Create some histograms
	TH1F *thetaH	=	new TH1F("thetaH","Theta distribution",180,0,180);
	TH1F *phiH		= new TH1F("phiH","Phi distribution" ,180,-180,180);
	TH1F *etaH		= new TH1F("etaH","Pseudorapidity",120,-12,12);
	TH1F *yH			= new TH1F("yH","Rapidity distribution",120,-12,12);
	TH1F *eH			= new TH1F("eH","Energy distribution",100,0,220);
	TH1F *eetaH		= new TH1F("eetaH","Pseudorapidity",120,0,12);
	TH1F *ptH			= new TH1F("ptH","Pt distribution",150,0,220);

	TH1F *eInH				= new TH1F("eInH", "electron in jet", 5, -0.5,4.5);
	TH1F *ptRelH			= new TH1F("ptRelH","relative electron momentum", 100,0.,5.);
	TH1F *eJetAlfaH		= new TH1F("eJetAlfaH", "angle electron - jet", 1000,0.,180.);
	TH1F *radMomH			= new TH1F("radMomH", "radial moment", 200, 0., 1.);
	TH1F *radMom_riH	= new TH1F("radMom_riH", "radial moment r_is", 100, 0., 1.); // upLim must be > jetRadius | 
	TH1F* angularH 		= new TH1F("angularH", "jet angularities", 100, 0.,200. );

	Int_t bins1[2] = {100000};
	Double_t xmin1[2] = {0.};
	Double_t xmax1[2] = {5.};

	THnSparseD* vrtxHsparse = new THnSparseD("vrtxHsparse","particle vertex X abs sparse",1, bins1, xmin1, xmax1); 

	TH1D* NsecVertices1H = new TH1D("NsecVertices1H", "number of vertices found in jet with at least 1 particles", 20, -0.5, 19.5);
	TH1D* NsecVertices2H = new TH1D("NsecVertices2H", "number of vertices found in jet with at least 2 particles", 20, -0.5, 19.5);
	TH1D* NsecVertices3H = new TH1D("NsecVertices3H", "number of vertices found in jet with at least 3 particles", 20, -0.5, 19.5);

	TH1D* SVnH = new TH1D("SVnH", "number of vertices found in jet with at least 3 particles and sigma restriction", 20, -0.5, 19.5);
	TH1D* SVrH = new TH1D("SVrH", "r of vertices found in jet with at least 3 particles and sigma restriction", 100, 0., 20.);


	Int_t pdgRange = 8000;
	TH1F *pdgH = new TH1F("PDG", "Particles PDG (PID)", pdgRange, -pdgRange/2+0.5, pdgRange/2+0.5);

	TH1F pt_pdgH[pdgRange];
	const char* part_names[pdgRange];

	TH1F* tmpH = new TH1F("","",150,0,220);
	for(Int_t i_hist=0; i_hist<pdgRange; ++i_hist){
		pt_pdgH[i_hist] = *tmpH;
	}



	TH2F *h7=new TH2F("h7","Effective area charged vs p_{T}^{jet}, all jets",50,0,100,50,0,1);
	h7->GetXaxis()->SetTitle("p_{T}^{jet}");
	h7->GetYaxis()->SetTitle("A");

	TH1F *h16=new TH1F("h16","all jet multiplicity",100,-0.5,99.5);
	h16->GetXaxis()->SetTitle("jet multiplicity");

	TH1F *h17=new TH1F("h17","number of all jets in event",20,-0.5,19.5);
	h17->GetXaxis()->SetTitle("#jet in event");

	TH1F *h28=new TH1F("h28","number of jets after cuts",20,-0.5,19.5);
	h28->GetXaxis()->SetTitle("#jet in event");

	TH2F *hJetEtaPhiAll=new TH2F("hJetEtaPhiAll","#eta^{jet} vs #phi^{jet} all jets",100,0,6.5,100,-5,5);
	hJetEtaPhiAll->GetXaxis()->SetTitle("#phi^{jet}");
	hJetEtaPhiAll->GetYaxis()->SetTitle("#eta^{jet}");
	TH2F *hJetEtaPhi=new TH2F("hJetEtaPhi","#eta^{jet} vs #phi^{jet}",100,0,6.5,100,-1,1);
	hJetEtaPhi->GetXaxis()->SetTitle("#phi^{jet}");
	hJetEtaPhi->GetYaxis()->SetTitle("#eta^{jet}");


	rl->LoadKinematics();
	rl->LoadHeader(); 

	// tree initialization
	TString* treeFileName = new TString("treeOfJets");
	*treeFileName += fileNameParams;
	*treeFileName += ".root";


	TFile* treeFile = new TFile(treeFileName->Data(),"RECREATE");
	treeFile->cd();

	typedef struct {
		Double_t ptRel;			// e- to jet relative pt
		Double_t radMom;		// radial moment
		Double_t angular;		// angularity
		Double_t svR;				// sec vrtx displacement in xy plane
		Int_t svN;					// number of sec vtx found in jet, with sigma_r restriction
		Int_t mult;					// jet multiplicity aka counts
		Int_t eIn;					// flag if e- in jet

	} MyJetObservables;
	static MyJetObservables jetObserv;

	AliAODJet* jet = new AliAODJet();

	TClonesArray* eventParticles = new TClonesArray("TParticle",1000);
	TClonesArray &evParticles = *eventParticles;

	TClonesArray* jetParticles = new TClonesArray("TParticle",1000);
	TClonesArray &jParticles = *jetParticles;

	TClonesArray* eventMCParticles = new TClonesArray("AliMCParticle",1000);
	TClonesArray &evMCParticles = *eventMCParticles;

	TTree* tree = new TTree("jetsTree","jets tree's title", 0);
	tree->Branch("jetsBranch", "AliAODJet", &jet);	
	tree->Branch("evParticlesBranch", &eventParticles);
	tree->Branch("jetParticlesBranch", &jetParticles);
	tree->Branch("jetObservablesBranch",&jetObserv,"ptRel/D:radMom:angular:svR:svN/I:mult:eIn");




	//
	//										Events loop 
	//

	Int_t counter1 = 0;
	Int_t counter2 = 0;
	Int_t counter3 = 0;
	for (Int_t nev=0; nev< evNumber; nev++) {
		if(nev < nEvToPrint || !(nev%100)) cout<<"\n\n============ EVENT "<<nev<<" ==================";
		rl->GetEvent(nev);
		AliStack* stack = rl->Stack();

		aliplist.Clear();
		JetFinderEvent.Clear();
		evParticles.Clear();


		TProcessID::SetObjectCount(0); // Needed for TRefs in AliCalTrkTrack and AliAODJet


		//
		// Loop over primary particles 
		//

		Int_t npart = stack->GetNprimary();
		Int_t offset = 0;
		for (Int_t part=0; part<npart; part++) {
				counter1++;
				TParticle *MPart = stack->Particle(part);
				Int_t mpart		= MPart->GetPdgCode();
				Int_t child1	= MPart->GetFirstDaughter();
				Int_t child2	= MPart->GetLastDaughter();	
				Int_t mother	= MPart->GetFirstMother();

				const char* pdg_name = MPart->GetName();

				Float_t Pt 		= MPart->Pt();
				Float_t E			= MPart->Energy();
				Float_t Pz 		= MPart->Pz();
				Float_t Py		= MPart->Py();
				Float_t Px		= MPart->Px();
				Float_t pT		= TMath::Sqrt(Px*Px+Py*Py);
				Float_t theta = MPart->Theta();
				Float_t phi	 	= MPart->Phi()-TMath::Pi();
				Float_t eta		= -TMath::Log(TMath::Tan(theta/2.));
				Float_t y			= 0.5*TMath::Log((E+Pz+1.e-13)/(E-Pz+1.e-13));
				
				if (child1 >= 0) { //cout<<"\nWARNING, \'continue\' in loop due to child1 >= 0 \tit's "<<pdg_name<<endl;
				offset++;
				counter2++;
				continue;
			}
				else if ((mpart < 6 && mpart >=-6) || mpart == 21) { //cout<<"\nWARNING, \'continue\' in loop due to q, q_bar or g in event \tit's "<<pdg_name<<endl;
				offset++;
				counter2++;
				continue;
			}


				thetaH->Fill(theta*180./TMath::Pi());
				phiH->Fill(phi*180./TMath::Pi());
				etaH->Fill(eta); 
				eetaH->Fill(eta,E); 
				yH->Fill(y);
				eH->Fill(E); 
				ptH->Fill(pT);

				if (mpart < -pdgRange/2 || mpart > pdgRange/2) printf("WARNING: particle pdg=%d out of range!",mpart);
				pdgH->Fill(mpart);

				(pt_pdgH[mpart+pdgRange/2]).Fill(Pt);
				part_names[mpart+pdgRange/2] = pdg_name;


			// ghosts i.e. have momentum exactly in z direction, 
			// see WARN. in type 1 or 2 (pair of q or g)
			if (pT==0.0) {
				cout<<"WARNING: ghost-candidate: pt="<<pT<<", eta="<<eta<<", theta="<<theta<<", type="<<pdg_name<<endl;
				//continue;
			}

		
			// key to fastjet being run -- filling JetFinderEvent
			// also filling TClonesArrays, which are to written in tree 
			new (evParticles[part-offset]) TParticle(*MPart);
			JetFinderEvent.AddCalTrkTrackKine((AliMCParticle*)evMCParticles[part-offset],1,1);
			counter3++;

		} // primary loop


		if(nev < nEvToPrint){
			cout<<"\nJetFinderEvent->N: "<<JetFinderEvent.GetNCalTrkTracks();
			cout<<"\nJetFinderEvent: "; JetFinderEvent.fJetCalTrkTrack->Print(); cout<<"end of JetFinderEvent"; 
			cout<<"\nJetFinderEvent2: \n"; 
			for(int ii=0;ii<JetFinderEvent.GetNCalTrkTracks();++ii){
				AliMCParticle amcp(*((AliMCParticle*)JetFinderEvent.fJetCalTrkTrack->At(ii)));
				amcp.Print();
				cout<<"\n Px: "<<amcp.Px()<<" Py: "<<amcp.Px()<<" Pz: "<<amcp.Pz();
				cout<<"\n Pt: "<<amcp.Pt(); cout<<" Eta: "<<amcp.Eta()<<" Phi: "<<amcp.Phi();			
			}
			cout<<"end of JetFinderEvent2"; 

		}


		//
		// Fastjet part
		//

		aod->ClearStd();
		FastJet->Reset();
		FastJet->SetCalTrkEvent(JetFinderEvent);
		FastJet->ProcessEvent();



		h17->Fill(aod->GetNJets());



		for (Int_t iJet = 0; iJet < aod->GetNJets(); iJet++) { //loop over jets
			if(nev < nEvToPrint) cout<<"#"<<iJet<<"__";

			jet = aod->GetJet(iJet);
			h7->Fill(jet->Pt(),jet->EffectiveAreaCharged());
			hJetEtaPhiAll->Fill(jet->Phi(),jet->Eta());
			h16->Fill(jet->GetRefTracks()->GetEntriesFast());

			jParticles.Clear();

			if(nev < nEvToPrint){
				cout<<"\njet nr "<<iJet<<": "; jet->GetRefTracks()->Print(); 
				for(int i=0;i<jet->GetRefTracks()->GetEntriesFast();++i) {cout<<"\n--->>"; cout<<((AliMCParticle*)jet->GetRefTracks()->At(i))->Pt();}
				cout<<"\nend of jet nr\n";
			}

		
			// filling TCA jetParticles: 
			// loop over jetParticles {
			// 		loop over stackParticles{
			//	 	 if jetPart == stackPart => Fill TCA jetParticles 
			//		}
			//	}
			Int_t nIterIn 	= stack->GetNprimary();
			Int_t nIterOut	= jet->GetRefTracks()->GetEntriesFast();
			Int_t jetPartCounter = 0;
		
			Int_t eInJet = 0;
			Double_t radialMoment = 0.;
			Double_t angularity = 0.;		
			Double_t ptRelMax = -1.;
			TClonesArray secVertices("TVector3", 20);
			Int_t secVcounts[20] = {};
			Int_t isecV = 0;
			Int_t NsecV[3] = {}; // how many secV with at least 1, 2 or 3 particles

			for(int iJetPart=0; iJetPart < nIterOut; iJetPart++){ // loop over particles found in jet
				// filling TCA jetParticles
				Int_t offset = 0;
				AliMCParticle* jetp = ((AliMCParticle*)jet->GetRefTracks()->At(iJetPart));
				for(int iStack=0; iStack < nIterIn; iStack++){ // loop over particles on stack
					if (nev < nEvToPrint) cout<<"\nloop: "<<iJetPart<<" "<<iStack;
					//alleventParticles = JetFinderEvent.fJetCalTrkTrack;
					//new (allevParticles[i_track]) AliJetCalTrkTrack(JetFinderEvent.GetCalTrkTrack(i_track));
					TParticle* stackp = stack->Particle(iStack);
					Float_t deltaPhi = abs(jetp->Phi() - stackp->Phi());
					Float_t deltaEn = abs(jetp->E() - stackp->Energy());
					//cout<<"\tPt: "<<jetp->Pt()<<" "<<stackp->Pt()<<"\t E: "<<jetp->E()<<" "<<stackp->Energy();
					if (stackp->GetFirstDaughter() >= 0) {/*cout<<"\nWARNING, \'continue\' in loop due to GetFirstDaughter >= 0 ";*/ offset++; continue;}
					// Fill
					if (deltaPhi < 1e-7 && deltaEn < 1e-7) {
						new (jParticles[jetPartCounter]) TParticle(*stackp);
						// if e-
						if (stackp->GetPdgCode() == 11) {
							eInJet++;
							TVector3 eleP(stackp->Px(), stackp->Py(), stackp->Pz());
							Double_t alfa = jet->MomentumVector()->Angle(eleP);
							Double_t ptRel = sin(alfa)*stackp->Pt();
							eJetAlfaH->Fill(alfa*180./TMath::Pi());
							ptRelMax = max(ptRelMax, ptRel); 
							//cout<<"PtRel: "<<ptRel<<" alfa: "<<alfa*180./TMath::Pi()<<" sin(alfa): "<<sin(alfa)<<" jet->Pt(): "<<jet->Pt()<<" e->Pt(): "<<stackp->Pt()<<endl;
						}
						jetPartCounter++; 
						if (nev < nEvToPrint) cout<<"\n\tparticle found\t"; 				
						if (nev < nEvToPrint) cout<<"\njetPartCounter: "<< jetPartCounter; 


						if (stackp->GetPdgCode() != 22){ // exclude photons from sec vertex searches
							TVector3 v(stackp->Vx(), stackp->Vy(), stackp->Vz());
							if (v != TVector3(0.,0.,0.) && stackp->R() > 0.){
								if (nev < nEvToPrint){ cout<<"secV candidate: "; v.Print();}
								Int_t indexFound = -1;
								for(int i=0; i<secVertices.GetEntries(); ++i) if( v == *((TVector3*)secVertices.At(i)) ) {
									indexFound = i; 
									if (nev < nEvToPrint) { cout<<"found at "<<i<<endl;} 
								} 
								if (indexFound != -1) secVcounts[indexFound] += 1;
								else{ new(secVertices[isecV]) TVector3(v); secVcounts[isecV] = 1; isecV++; 
									if (nev < nEvToPrint) { cout<<"jet"<<iJet<<" NEW secVertex: "<<isecV<<" "; v.Print(); }	 
								}
							}
						}
						break; // break from stack particles loop after stack particle is equal to jet particle
					} // jet particles found at stack
				}// stack particles loop

				//calculating radial moment, angularity defined as in: http://jets.physics.harvard.edu/qvg/
				Double_t deltaPhi = max(jet->Phi(),jetp->Phi()) - min(jet->Phi(),jetp->Phi());
				Double_t deltaEta = max(jet->Eta(),jetp->Eta()) - min(jet->Eta(),jetp->Eta());
				deltaPhi = deltaPhi < TMath::Pi() ? deltaPhi : deltaPhi-2*TMath::Pi();
				deltaEta = deltaEta < TMath::Pi() ? deltaEta : deltaEta-2*TMath::Pi();

				Double_t r_i = sqrt(pow(deltaPhi,2) + pow(deltaEta,2));
				// radial moment
				radialMoment += jetp->Pt()/jet->Pt() * r_i;
				radMom_riH->Fill(r_i);
				if (r_i > jetRadius*1.4 || r_i < 0.) { cout<<"\n\nWARNING: r_i in radial moment > 1.4*jetRadius or < 0\n"; cout<<r_i<<endl<<"\tphiJet, phiPart, etaJet, etaPart: "<< jet->Phi()<<" "<<jetp->Phi()<<"\t"<< jet->Eta()<<" "<<jetp->Eta()<<endl;
				}

				// angularity
				Double_t a_angu = 1.; // parameter
				Double_t theta = TMath::Pi()*r_i/(2*jetRadius);
				Double_t f = pow(sin(theta), a_angu) * pow(1-cos(theta), 1-a_angu);
				angularity += jetp->E() * f;


				//vertex
				Double_t vx = abs(jetp->Xv());
				Double_t vy = abs(jetp->Yv()); 
				
				Double_t Vx[1] = {vx};
				vrtxHsparse->Fill(Vx);

	//			cout<<"\njet nr "<<iJet<<": "; //jet->GetRefTracks()->Print(); 
	//			for(int ii=0;ii<jet->GetRefTracks()->GetEntriesFast();++ii) {
	//				cout<<"part "<<ii;
	//				AliMCParticle* amcp = (AliMCParticle*)jet->GetRefTracks()->At(ii);
	//				cout<<" phi eta pt: "<<amcp->Phi()<<" "<<amcp->Eta()<<" "<<amcp->Pt()<<endl;
	//			}
	//			cout<<"\nend of jet nr\n";


			} // jet particles loop
			radMomH->Fill(radialMoment);
			angularH->Fill(angularity);
			eInH->Fill(eInJet);
			if (ptRelMax > 0) ptRelH->Fill(ptRelMax); //ptRelMax = -1 means no electron in jet


			//count sec vertices with min X particles
			Int_t Nsv = 0;
			Double_t RsvMax = -1;
			for(int i=0; i<secVertices.GetEntries(); ++i){
				Int_t n_tracks = secVcounts[i];
				if (n_tracks > 0) NsecV[0]+=1;
				if (n_tracks > 1) NsecV[1]+=1;
				if (n_tracks > 2) NsecV[2]+=1;


				if (secVcounts[i] < 3) continue;
				if (nev < nEvToPrint){ cout<<"secVertex with N>=3: counts and (x,y,z):	N="<<secVcounts[i]<<"\t"; (*(TVector3*)secVertices.At(i)).Print();}
				// sigma_vtx as a function of n_trakcs
				// http://personalpages.to.infn.it/~milano/presentations/LeonardoMilano_LHCBeyond.pdf s.29
				// probably better ref. https://cds.cern.ch/record/1648854/files/CERN-PH-EP-2014-031.pdf but for 7 TeV, not 0.9
				//Double_t sigma_vtx = sqrt( pow(166,2) + pow(540/sqrt(pow(n_tracks,0.9)),2) ) *1e-3; // um->mm, 0.9 TeV
				Double_t sigma_vtx = sqrt( pow(40,2) + pow(297/sqrt(pow(n_tracks,1.37)),2) ) *1e-3; // um->mm, 7 TeV
				TVector3 sv = *(TVector3*)secVertices.At(i);
				Double_t R = sqrt(sv.x()*sv.x() + sv.y()*sv.y());
				if (nev < nEvToPrint) cout<<"\t\tsigma_vtx = "<<sigma_vtx<<"\t sv->r = "<< R <<endl;
				if (R > 3*sigma_vtx*sqrt(2)){
					if (nev < nEvToPrint) cout<<"\t>>> accepted <<<\n";
					SVrH->Fill(R);
					RsvMax = max(R, RsvMax);
					Nsv++;
				}																																																																							
			}

			NsecVertices1H->Fill(NsecV[0]);
			NsecVertices2H->Fill(NsecV[1]);
			NsecVertices3H->Fill(NsecV[2]);
			SVnH->Fill(Nsv);



			// observables struct filling
			jetObserv.mult		= jet->GetRefTracks()->GetEntriesFast();
			jetObserv.eIn			= eInJet;
			jetObserv.svN			= Nsv;
			jetObserv.svR			= RsvMax;
			jetObserv.ptRel		= ptRelMax;
			jetObserv.radMom	= radialMoment;
			jetObserv.angular	= angularity;



			if (iJet != 0) {
				// write eventParticles only once per event
				// eventParticles are at once filled with all particles in event 
				// while we are now in jet loop
				evParticles.Clear(); 
				evMCParticles.Clear();
			}

			// some printing
			if (nev < nEvToPrint) {cout<<"\njet->GetRefTracks()->Print():"; jet->GetRefTracks()->Print(); cout<<"end of jet->GetRefTracks()->Print()\n";}
			if (nev < nEvToPrint) {cout<<"\neventParticles->Print():"; eventParticles->Print(); cout<<"end of eventParticles->Print()\n";}
			if (nev < nEvToPrint) {cout<<"\njetParticles->Print():"; jetParticles->Print(); cout<<"end of jetParticles->Print()\n";}


			// Filling tree
			if (nev < nEvToPrint) cout<<"\nJust before tree->Fill ...\t";
			tree->Fill();
			if (nev < nEvToPrint) cout<<"... Just after tree->Fill\n"; 

		}//loop over jets

		AliAODTrack *obj =NULL;

	} //event loop


	// Finish 
	treeFile->cd();
	//tree->Print();
	tree->Write();

	//cout<<"\n\n\n counter1: "<<counter1<<" counter2: "<<counter2<<" counter3: "<<counter3<<endl;



//Create a canvas, set the view range, show histograms

	TString* histFileName = new TString("histos"); 
	*histFileName += fileNameParams;
	*histFileName += ".root";

	TFile* histFile = new TFile(histFileName->Data(),"RECREATE");
	histFile->cd();


	// Canvases
/*
	TCanvas *c = new TCanvas("c","Canvas 1",400,10,600,700);
		c->Divide(2,2);
		c->cd(1); h7->Draw();
		c->cd(2); h16->Draw();
		c->cd(3); h17->Draw();
		c->cd(4); hJetEtaPhiAll->Draw();
		c->Write();


	TCanvas *c1 = new TCanvas("c1","Canvas 1",400,10,600,700);
		c1->Divide(2,2);
		c1->cd(1); ptH->Draw();
		c1->cd(2); etaH->Draw();
		c1->cd(3); yH->Draw();
		c1->cd(4); eH->Draw();
		c1->Write();


	TCanvas *c2 = new TCanvas("c2","Canvas 2",400,10,600,700);
		c2->Divide(2,2);
		c2->cd(1); phiH->Draw();
		c2->cd(2); thetaH->Draw();
		c2->cd(3); eetaH->Draw();
		c2->cd(4); etaH->Draw();
		c2->Write();
 

	TCanvas *c3 = new TCanvas("c3","Canvas 3 - PDG",400,10,600,700);
		pdgH->Draw();

		Int_t N_part_types = 0;
		Int_t NentriesLimit = 100;
		for (Int_t i_pt = 0; i_pt < pdgRange; ++i_pt){
				if (pt_pdgH[i_pt].GetEntries() < NentriesLimit) continue;
				N_part_types++;
		//printf("\n not empty bin: >>>\t %d \t pdg: %d <=> name: %s \t\t\t NEntries: %f <<<", i_pt, i_pt-pdgRange/2, part_names[i_pt], pt_pdgH[i_pt].GetEntries());					

		}
		for(int i=1; i<pdgH->GetNbinsX(); ++i) if(pdgH->GetBinContent(i) > 0) cout<< "found: pdg="<<i-pdgRange/2<<", name="<<part_names[i]<<endl;
		c3->Write();


	TCanvas *c4 = new TCanvas("c4","Canvas - Pt for particles species",400,10,600,700);
		int a = sqrt(N_part_types), b = sqrt(N_part_types);
		if (a*b < N_part_types) a+=1;
		c4->Divide(a,b);

		int N_plt = 0;
		for (Int_t i_pt = 0; i_pt < pdgRange; ++i_pt){
			if (pt_pdgH[i_pt].GetEntries() < NentriesLimit) continue; //
			//printf("%c \n", part_names[i_pt][0]);
			N_plt++; 

			TH1F* pt_tmp = new TH1F("","", 150, 0, 15);
			*pt_tmp = pt_pdgH[i_pt];
			pt_tmp->SetTitle(part_names[i_pt]);
			c4->cd(N_plt);
			pt_tmp->Draw();	
		}
		c4->Write();

	TCanvas *c5 = new TCanvas("c5","selected Pt",400,10,600,700);
		c5->Divide(2,2);

		TList *listPions = new TList;
		TList *listKaons = new TList;
		TList *listProtons = new TList;

		for (Int_t i_pt = 0; i_pt < pdgRange; ++i_pt){
			if (pt_pdgH[i_pt].GetEntries() < NentriesLimit) continue;
			TString pname(part_names[i_pt]);

			if (pname.BeginsWith("pi")) listPions->Add(pt_pdgH+i_pt);
			if (pname.BeginsWith("K")) listKaons->Add(pt_pdgH+i_pt);
			if (pname.BeginsWith("proton")) listProtons->Add(pt_pdgH+i_pt);
		}

		TH1F *ptPions = (TH1F*)ptH->Clone("ptPions");
			ptPions->SetTitle("ptPions");
			ptPions->Reset();
			ptPions->Merge(listPions);
			ptPions->Write();
		TH1F *ptKaons = (TH1F*)ptH->Clone("ptKaons");
			ptKaons->SetTitle("ptKaons");
			ptKaons->Reset();
			ptKaons->Merge(listKaons);
			ptKaons->Write();
		TH1F *ptProtons = (TH1F*)ptH->Clone("ptProtons");
			ptProtons->SetTitle("ptProtons");
			ptProtons->Reset();
			ptProtons->Merge(listProtons);
			ptProtons->Write();

		c5->cd(1); ptPions->Draw();
		c5->cd(2); ptKaons->Draw();
		c5->cd(3); ptProtons->Draw();
		c5->Write();
*/



	thetaH->Write();
	phiH->Write();
	etaH->Write();
	yH->Write();
	eH->Write();
	eetaH->Write();
	ptH->Write();
	pdgH->Write();

	eInH->Write();
	ptRelH->Write();
	eJetAlfaH->Write();
	radMomH->Write();
	radMom_riH->Write();
	angularH->Write();
		
	vrtxHsparse->Write();
	NsecVertices1H->Write();
	NsecVertices2H->Write();
	NsecVertices3H->Write();

	SVnH->Write();
	SVrH->Write();

	h7->Write();
	h16->Write();
	h17->Write();
	h28->Write();

	hJetEtaPhiAll->Write();

}







