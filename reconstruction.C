/* * * 

Variation of standard reconstruction.C
with additional jet tagging based on:
- presence of charmed/beauty hadron in jet 
  (refered to as Exp - experimental method)
- last parton before hadronization (TrueLast)
- very original parton (TrueFirst)

usually abs(PDG codes) are used in resulting tree, with 2 exceptions:
0 for gluons
2 for protons (they appear frequently in TrueFirst)


Macro can be run with:
$ aliroot -l -q 'run_recon.C(25,"DATA/LHC17f8g_20_255618_090/")'

where: 
- 25 is no. events 
- "DATA/LHC17f8g_20_255618_090/" is the relative path to directory containing Kinematics.root and galice.root

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
#include "TString.h"
#include <typeinfo>

#include "AliGenerator.h"
#include "AliPDG.h"
#include "AliRunLoader.h"
#include "AliRun.h"
#include "AliStack.h"
#include "AliHeader.h"
#include "AliGenPythia.h"  //### E-TCint::AutoLoadCallback: failure loading library libAliPythia6 for class AliGenPythia
#include "AliPythia.h"

//for jets:
#include "AliFastJetHeaderV1.h"
#include "AliFastJetFinder.h"
#include "AliAODEvent.h"
#include "AliMCParticle.h"
//#include "AliJetCalTrkTrackKine.h"

#endif


void print_ancestors(TParticle* part0, AliStack* stack){
  printf("%s   ", part0->GetName());
  TParticle* part = (TParticle*)part0->Clone();
  while(part->GetFirstMother() != -1){
    part = stack->Particle(part->GetFirstMother());
    printf("<-%s(%d)", part->GetName(), part->GetPdgCode());
  }

  printf("\n\t X ");
  part = (TParticle*)part0->Clone();
  while(part->GetFirstMother() != -1){
    printf("<-(%d)", part->GetFirstMother());
    part = stack->Particle(part->GetFirstMother());
  }
}


void print_children(TParticle* part0, int stack_pid, AliStack* stack, int iGeneration=0){
  if (iGeneration == 0) printf("\n*");
  printf("-> %s (%d)[%.2f GeV]", part0->GetName(), stack_pid, part0->Energy());
  TParticle* part = (TParticle*)part0->Clone();
  if(part->GetFirstDaughter() != -1){
    for(int iDaughter=part->GetFirstDaughter(); iDaughter <= part->GetLastDaughter() ; ++iDaughter){
      printf("\n");
      TParticle* part_tmp = stack->Particle(iDaughter);
      TParticle* partMother = (TParticle*) part_tmp->Clone();
      for(int i=0; i<iGeneration; ++i)  printf(" .");
      print_children(part_tmp, iDaughter, stack, iGeneration+1);
    }
  }

}



void reconstruction(Int_t evNumber=200, const char* pathToFile="") {

  TString* galiceFileName = new TString(pathToFile);
  *galiceFileName += "/galice.root";


  // Connect the Root Galice file containing Geometry, Kine and Hits
  AliRunLoader* rl = AliRunLoader::Open(galiceFileName->Data());

  TDatabasePDG* PDGbase = new TDatabasePDG();
    

  // main parameter for anti-kt algorithm
  Float_t jetRadius = 0.4;


  //
  // Jet finder settings go via the FastJetHeader
  //
  AliFastJetHeaderV1 *header = new AliFastJetHeaderV1;
  header->SetBGMode(0);
  header->SetRparam(jetRadius); 

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
  TH1F *thetaH      = new TH1F("thetaH","Theta distribution",180,0,180);
  TH1F *phiH        = new TH1F("phiH","Phi distribution" ,180,-180,180);
  TH1F *etaH        = new TH1F("etaH","Pseudorapidity",120,-12,12);
  TH1F *yH          = new TH1F("yH","Rapidity distribution",120,-12,12);
  TH1F *eH          = new TH1F("eH","Energy distribution",100,0,220);
  TH1F *eetaH       = new TH1F("eetaH","Pseudorapidity",120,0,12);
  TH1F *ptH         = new TH1F("ptH","Pt distribution",150,0,220);

  TH1F *eInH        = new TH1F("eInH", "electron in jet", 5, -0.5,4.5);
  TH1F *ptRelH      = new TH1F("ptRelH","relative electron momentum", 100,0.,5.);
  TH1F *eJetAlfaH   = new TH1F("eJetAlfaH", "angle electron - jet", 1000,0.,180.);
  TH1F *radMomH     = new TH1F("radMomH", "radial moment", 200, 0., 1.);
  TH1F *radMom_riH  = new TH1F("radMom_riH", "radial moment r_is", 100, 0., 1.); // upLim must be > jetRadius | 
  TH1F* angularH    = new TH1F("angularH", "jet angularities", 100, 0.,200. );

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
  TString* treeFileName = new TString(pathToFile);
  *treeFileName += "/jetsTree.root";


  TFile* treeFile = new TFile(treeFileName->Data(),"RECREATE");
  treeFile->cd();

  typedef struct {
    Double_t ptRel;      // e- to jet relative pt
    Double_t radMom;    // radial moment
    Double_t angular;    // angularity
    Double_t svR;        // sec vrtx displacement in xy plane
    Int_t svN;          // number of sec vtx found in jet, with sigma_r restriction
    Int_t mult;          // jet multiplicity aka counts
    Int_t eIn;          // flag if e- in jet
    Int_t tagExp;
    Int_t tagTrueLast;
    Int_t tagTrueFirst;
    Int_t tagExpC;
    Int_t tagTrueLastC;
    Int_t tagTrueFirstC;
    Int_t tagExpB;
    Int_t tagTrueLastB;
    Int_t tagTrueFirstB;

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
  tree->Branch("jetObservablesBranch",&jetObserv,"ptRel/D:radMom:angular:svR:svN/I:mult:eIn:tagExp:tagTrueLast:tagTrueFirst:tagExpC:tagTrueLastC:tagTrueFirstC:tagExpB:tagTrueLastB:tagTrueFirstB");

  //
  //                    Events loop 
  //

  Int_t counter1 = 0;
  Int_t counter2 = 0;
  Int_t counter3 = 0;

  //TDatabasePDG* base = new TDatabasePDG();
  for (Int_t nev=0; nev < evNumber; nev++) {
    cout<<"\n\n============ EVENT "<<nev<<" ==================";
    if(!(nev%100)) cout<<"\n\n============ EVENT "<<nev<<" ==================";
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
        Int_t mpart    = MPart->GetPdgCode();
        Int_t child1   = MPart->GetFirstDaughter();
        Int_t child2   = MPart->GetLastDaughter();  
        Int_t mother   = MPart->GetFirstMother();

        const char* pdg_name = MPart->GetName();

        Float_t Pt    = MPart->Pt();
        Float_t E     = MPart->Energy();
        Float_t Pz    = MPart->Pz();
        Float_t Py    = MPart->Py();
        Float_t Px    = MPart->Px();
        Float_t pT    = TMath::Sqrt(Px*Px+Py*Py);
        Float_t theta   = MPart->Theta();
        Float_t phi     = MPart->Phi()-TMath::Pi();
        Float_t eta     = -TMath::Log(TMath::Tan(theta/2.));
        Float_t y     = 0.5*TMath::Log((E+Pz+1.e-13)/(E-Pz+1.e-13));
        
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
      new (evMCParticles[part-offset]) AliMCParticle(MPart);
      JetFinderEvent.AddCalTrkTrackKine((AliMCParticle*)evMCParticles[part-offset],1,1);
      counter3++;

    } // primary loop
    //printf("\nEnd of primary loop\n");




    //
    // Fastjet part
    //

    aod->ClearStd();
    FastJet->Reset();
    FastJet->SetCalTrkEvent(JetFinderEvent);
    FastJet->ProcessEvent();

    h17->Fill(aod->GetNJets());


    for (Int_t iJet = 0; iJet < aod->GetNJets(); iJet++) { //loop over jets
      printf("\n\n\t\t# JET nr %d\n", iJet);

      jet = aod->GetJet(iJet);
      h7->Fill(jet->Pt(),jet->EffectiveAreaCharged());
      hJetEtaPhiAll->Fill(jet->Phi(),jet->Eta());
      h16->Fill(jet->GetRefTracks()->GetEntriesFast());

      jParticles.Clear();

    
      // filling TCA jetParticles: 
      // loop over jetParticles {
      //     loop over stackParticles{
      //      if jetPart == stackPart => Fill TCA jetParticles 
      //    }
      //  }
      Int_t nIterIn   = stack->GetNprimary();
      Int_t nIterOut  = jet->GetRefTracks()->GetEntriesFast();
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


            if (stackp->GetPdgCode() != 22){ // exclude photons from sec vertex searches
              TVector3 v(stackp->Vx(), stackp->Vy(), stackp->Vz());
              if (v != TVector3(0.,0.,0.) && stackp->R() > 0.){
                Int_t indexFound = -1;
                for(int i=0; i<secVertices.GetEntries(); ++i) if( v == *((TVector3*)secVertices.At(i)) ) {
                  indexFound = i; 
                } 
                if (indexFound != -1) secVcounts[indexFound] += 1;
                else{ new(secVertices[isecV]) TVector3(v); secVcounts[isecV] = 1; isecV++; }
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
        // sigma_vtx as a function of n_trakcs
        // http://personalpages.to.infn.it/~milano/presentations/LeonardoMilano_LHCBeyond.pdf s.29
        // probably better ref. https://cds.cern.ch/record/1648854/files/CERN-PH-EP-2014-031.pdf but for 7 TeV, not 0.9
        //Double_t sigma_vtx = sqrt( pow(166,2) + pow(540/sqrt(pow(n_tracks,0.9)),2) ) *1e-3; // um->mm, 0.9 TeV
        Double_t sigma_vtx = sqrt( pow(40,2) + pow(297/sqrt(pow(n_tracks,1.37)),2) ) *1e-3; // um->mm, 7 TeV
        TVector3 sv = *(TVector3*)secVertices.At(i);
        Double_t R = sqrt(sv.x()*sv.x() + sv.y()*sv.y());
        if (R > 3*sigma_vtx*sqrt(2)){
          SVrH->Fill(R);
          RsvMax = max(R, RsvMax);
          Nsv++;
        }                                                                                                                                              
      }

      NsecVertices1H->Fill(NsecV[0]);
      NsecVertices2H->Fill(NsecV[1]);
      NsecVertices3H->Fill(NsecV[2]);
      SVnH->Fill(Nsv);



      ////////////////////////////////////
      //
      // Extended TAGGING
      //
      ////////////////////////////////////


      Int_t tagExp = 0;
      Int_t tagTrueLast = 0;
      Int_t tagTrueFirst = 0;

      Int_t tagExpC = 0;
      Int_t tagTrueLastC = 0;
      Int_t tagTrueFirstC = 0;
      Int_t tagExpB = 0;
      Int_t tagTrueLastB = 0;
      Int_t tagTrueFirstB = 0;

      // loop over jet particles for tagging
      for(int iJetPart=0; iJetPart < jet->GetRefTracks()->GetEntriesFast(); iJetPart++){
          printf("\n %d.", iJetPart);

          TParticle* jetp = ((AliMCParticle*)jet->GetRefTracks()->At(iJetPart))->Particle();
          print_ancestors(jetp, stack);      
          TParticle* motherp = jetp;

          Int_t tagTrueFirstTmp = 0;
          Int_t tagTrueLastTmp = 0;
          Int_t tagExpTmp = 0;

          // tag 1: true first
          motherp = jetp;
          while(kTRUE){
            motherp = stack->Particle(motherp->GetFirstMother());
            if (motherp->GetFirstMother() == -1){ 
                printf("\n\t true first: %s (%d)", motherp->GetName(), motherp->GetPdgCode());
                if (abs(motherp->GetPdgCode()) == 4) tagTrueFirstC = 1;
                if (abs(motherp->GetPdgCode()) == 5) tagTrueFirstB = 1;

                if (abs(motherp->GetPdgCode()) == 21) tagTrueFirstTmp = 0;   // gluon
                else if (abs(motherp->GetPdgCode()) == 2212) tagTrueFirstTmp = 2;  // proton
                else tagTrueFirstTmp = abs(motherp->GetPdgCode());

                break;
            }
          }
          tagTrueFirst = max(tagTrueFirst, tagTrueFirstTmp);  // b > c > l;  (b,c) => b :
                                                              // if b and c quarks appered in parts' history - assign b

          // tag 2: true last
          motherp = jetp;
          while(kTRUE){
            motherp = stack->Particle(motherp->GetFirstMother());
            if (abs(motherp->GetPdgCode()) < 7 || motherp->GetPdgCode() == 21){ 
                printf("\n\t true last: %s (%d)", motherp->GetName(), motherp->GetPdgCode());  

                if (abs(motherp->GetPdgCode()) == 4) tagTrueLastC = 1;
                if (abs(motherp->GetPdgCode()) == 5) tagTrueLastB = 1;

                if (abs(motherp->GetPdgCode()) == 21) tagTrueLastTmp = 0;
                else if (abs(motherp->GetPdgCode()) == 2212) tagTrueLastTmp = 2;
                else tagTrueLastTmp = abs(motherp->GetPdgCode());

                break;
            }
            if (motherp->GetFirstMother() == -1){ 
                printf("\n\t true last (end): %s (%d)", motherp->GetName(), motherp->GetPdgCode());
                tagTrueLastTmp = -1;
                break;
            }
          }
          tagTrueLast = max(tagTrueLast, tagTrueLastTmp); 
          
          // tag 3: exp
          motherp = jetp;
          while(kTRUE){
              motherp = stack->Particle(motherp->GetFirstMother());
              Int_t code = abs(motherp->GetPdgCode());

              
              if (  ! TDatabasePDG::Instance()->GetParticle(code) ){  // segmentation for part. XXX(e.g. 4124, 14122)->Charm()
                  printf("\n\n\n\nW A R N I N G: XXX particle detected \n\n\n\n");
                  continue;
              }
              //else{ cout<<"$ "<< TDatabasePDG::Instance()->GetParticle(code) <<endl;} 

              Bool_t charm = (motherp->Charm() != 0  ||(code>400 && code<500));
              Bool_t beauty = (motherp->Beauty() != 0 || (code>500 && code<600));
              if (charm || beauty){
                printf("\n\tEXP %s (%d) %d", motherp->GetName(), motherp->GetPdgCode(), motherp->GetPDG()->Stable());

                if (charm)  tagExpC = 1;
                if (beauty) tagExpB = 1;

                if (abs(motherp->GetPdgCode()) == 21) tagExpTmp = 0;
                else if (abs(motherp->GetPdgCode()) == 2212) tagExpTmp = 2;
                else tagExpTmp = int(code/100);

                break;
              }
              if (motherp->GetFirstMother() == -1){ 
                printf("\n\t EXP (end): light"); 
                tagExpTmp = -1; 
                break;
              }
          }
          tagExp = max(tagExp, tagExpTmp); 

          if (iJetPart+1 == jet->GetRefTracks()->GetEntriesFast()){
              // print report for jet      //made from particles printed so far, it's updated after every particle in jet
		      printf("\n\n------- true first=%d  true last=%d  exp=%d", tagTrueFirst, tagTrueLast, tagExp);
		      printf("\n------c: true first=%d  true last=%d  exp=%d", tagTrueFirstC, tagTrueLastC, tagExpC);
		      printf("\n------b: true first=%d  true last=%d  exp=%d", tagTrueFirstB, tagTrueLastB, tagExpB);
          }

      }



      // observables struct filling
      jetObserv.mult    = jet->GetRefTracks()->GetEntriesFast();
      jetObserv.eIn      = eInJet;
      jetObserv.svN      = Nsv;
      jetObserv.svR      = RsvMax;
      jetObserv.ptRel    = ptRelMax;
      jetObserv.radMom  = radialMoment;
      jetObserv.angular  = angularity;
      jetObserv.tagTrueFirst = tagTrueFirst;
      jetObserv.tagTrueLast = tagTrueLast;
      jetObserv.tagExp = tagExp;

      jetObserv.tagTrueFirstC = tagTrueFirstC;
      jetObserv.tagTrueLastC = tagTrueLastC;
      jetObserv.tagExpC = tagExpC;
      jetObserv.tagTrueFirstB = tagTrueFirstB;
      jetObserv.tagTrueLastB = tagTrueLastB;
      jetObserv.tagExpB = tagExpB;
      




      if (iJet != 0) {
        // write eventParticles only once per event
        // eventParticles are at once filled with all particles in event 
        // while we are now in jet loop
        evParticles.Clear(); 
        evMCParticles.Clear();
      }



      // Filling tree
      tree->Fill();


    }//loop over jets

    AliAODTrack *obj =NULL;

  } //event loop


  // Finish 
  treeFile->cd();
  //tree->Print();
  tree->Write();

  //cout<<"\n\n\n counter1: "<<counter1<<" counter2: "<<counter2<<" counter3: "<<counter3<<endl;



//Create a canvas, set the view range, show histograms

  TString* histFileName = new TString(pathToFile); 
  *histFileName += "/histos.root";

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







