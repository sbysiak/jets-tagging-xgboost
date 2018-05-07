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
#include "TRandom3.h"
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

#include "AliAnalysisTaskEmcalJet.h"
#include "AliJetContainer.h"
#include "AliEmcalJet.h"
#include "AliHFJetsTaggingVertex.h"
#include "AliAnalysisTaskJetExtractorHF.h"

#endif


void calc_observables_flavour(const char* pathToFile="", const char* flavour="") {
// TODO:
// 1) convert into a useable script
// 2) vertices:
//    - old way: check if svN is always equal to no. constituents
//                -> if yes replace """ v == *((TVector3*)secVertices.At(i)) """ with v - *((TVector3*)secVertices.At(i)) < sigma
//    - fill old-way verices without sigma condition (it's already included in det. simulation)
//    - use vertices from tree
// 3) tagging

  cout<<"starting run for flavour: "<<flavour<<endl;

  Double_t jetRadius = 0.4;

  TH1F *thetaH      = new TH1F("thetaH","Theta distribution",180,0,180);
  TH1F *phiH        = new TH1F("phiH","Phi distribution" ,90,-180,180);
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

  TH1I *isPhysPrimH = new TH1I("isPhysPrimH", "isPhysicalPrimary in all particles from kinematic tree", 5, -0.5, 4.5);
  TH1I *isPhysPrimFJH = new TH1I("isPhysPrimFJH", "isPhysicalPrimary in FJ input", 5, -0.5, 4.5);
  TH1I *eMotherPdgH = new TH1I("eMotherPdgH", "pdg of electron's mother (jet electrons only)", 2000, -1000,1000);

  Int_t bins1[2] = {100000};
  Double_t xmin1[2] = {0.};
  Double_t xmax1[2] = {5.};

  THnSparseD* vrtxHsparse = new THnSparseD("vrtxHsparse","particle vertex X abs sparse",1, bins1, xmin1, xmax1);

  TH1D* NsecVertices1H = new TH1D("NsecVertices1H", "number of vertices found in jet with at least 1 particles", 20, -0.5, 19.5);
  TH1D* NsecVertices2H = new TH1D("NsecVertices2H", "number of vertices found in jet with at least 2 particles", 20, -0.5, 19.5);
  TH1D* NsecVertices3H = new TH1D("NsecVertices3H", "number of vertices found in jet with at least 3 particles", 20, -0.5, 19.5);

  TH1D* SVnH = new TH1D("SVnH", "number of vertices found in jet with at least 3 particles and sigma restriction", 20, -0.5, 19.5);
  TH1D* SVrH = new TH1D("SVrH", "r of vertices found in jet with at least 3 particles and sigma restriction", 100, 0., 20.);
  TH1D* ratioRsigmaRH = new TH1D("ratopRsigmaRH", "R / sigmaR", 100, 0., 10.);

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



  TString* inputFileName = new TString(pathToFile);
  *inputFileName += "/AnalysisResults.root";
  TFile *inputFile = TFile::Open(inputFileName->Data());
  TDirectory* dir = inputFile->GetDirectory("ChargedJetsHadronCF");
  TString listName = TString("AliAnalysisTaskJetExtractorHF_Jet_AKTChargedR040_tracks_pT0150_E_scheme_RhoR020KT_");
  listName += flavour;
  listName += "Jets_histos";
  AliEmcalList* emcallist_flavour = (AliEmcalList*) dir->Get(listName.Data());
  TTree* tree_flavour = (TTree*) emcallist_flavour->Last();

  // create tree
  TString* jetObservFileName = new TString(pathToFile);
  *jetObservFileName += "/jetObserv.root";
  TFile* jetObservFile = new TFile(jetObservFileName->Data(),"UPDATE");
  jetObservFile->cd();
  TDirectory* flavourDir = jetObservFile->mkdir(flavour);
  flavourDir->cd();

  typedef struct {
    Double_t pT;
    Double_t phi;
    Double_t eta;

    Double_t ptRel;      // e- to jet relative pt
    Double_t radMom;     // radial moment
    Double_t angular;    // angularity
    Double_t svR;        // sec vrtx displacement in xy plane
    Int_t svN;           // number of sec vtx found in jet, with sigma_r restriction
    Int_t mult;          // jet multiplicity aka counts
    Int_t eIn;           // flag if e- in jet
    //Int_t eMotherPdg;
    //Int_t tagExp;
    //Int_t tagTrueLast;
    //Int_t tagTrueFirst;
    //Int_t tagExpC;
    //Int_t tagTrueLastC;
    //Int_t tagTrueFirstC;
    //Int_t tagExpB;
    //Int_t tagTrueLastB;
    //Int_t tagTrueFirstB;


  } MyJetObservables;
  static MyJetObservables jetObserv;

  TTree* jetObservTree = new TTree("jetObservTree","jets tree's title", 0);
  jetObservTree->Branch("jetObservablesBranch",&jetObserv,"pT/D:phi:eta:ptRel:radMom:angular:svR:svN/I:mult:eIn");




  AliBasicJet *jet = new AliBasicJet();
  TBranch* jetBranch = tree_flavour->GetBranch("Jets");
  jetBranch->SetAddress(&jet);

  Int_t njets = tree_flavour->GetEntries();
  // loop over jets
  for(int i=0;i<njets;++i){
    tree_flavour->GetEntry(i);

    TVector3 jetMom = TVector3();
    jetMom.SetPtEtaPhi(jet->Pt(), jet->Eta(), jet->Phi());

    Int_t eInJet = 0;
    Double_t radialMoment = 0.;
    Double_t angularity = 0.;
    Double_t ptRelMax = -1.;
    TClonesArray secVertices("TVector3", 20);
    Int_t secVcounts[20] = {};
    Int_t isecV = 0;
    Int_t NsecV[3] = {}; // how many secV with at least 1, 2 or 3 particles

    Int_t mult = jet->GetNumbersOfConstituents();
    AliBasicJetConstituent* constit = new AliBasicJetConstituent();
    Bool_t unknown_PIDs_only = kTRUE;
    // loop over jet constituents
    for(int iConstit=0; iConstit < mult; iConstit++){
      constit = jet->GetJetConstituent(iConstit);

      TVector3 constitMom = TVector3();
      constitMom.SetPtEtaPhi(constit->Pt(), constit->Eta(), constit->Phi());

      if (constit->PID()->TruthPID() != 9) unknown_PIDs_only = kFALSE;
      // if e-
      if (constit->PID()->TruthPID() == 0){
        eInJet++;
        TVector3 eleMom(constitMom.Px(), constitMom.Py(), constitMom.Pz());
        Double_t alfa = jetMom.Angle(eleMom);
        Double_t ptRel = sin(alfa)*constit->Pt();
        ptRelMax = max(ptRelMax, ptRel);

        eJetAlfaH->Fill(alfa*180./TMath::Pi());
        //cout<<"PtRel: "<<ptRel<<" alfa: "<<alfa*180./TMath::Pi()<<" sin(alfa): "<<sin(alfa)<<" jet->Pt(): "<<jet->Pt()<<" e->Pt(): "<<constit->Pt()<<endl;
      }

      // find vertices
      /* TODO #2:
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
      */


      //calculating radial moment, angularity defined as in: http://jets.physics.harvard.edu/qvg/
      Double_t deltaPhi = max(jetMom.Phi(),constitMom.Phi()) - min(jetMom.Phi(),constitMom.Phi());
      Double_t deltaEta = max(jetMom.Eta(),constitMom.Eta()) - min(jetMom.Eta(),constitMom.Eta());
      deltaPhi = deltaPhi < TMath::Pi() ? deltaPhi : deltaPhi-2*TMath::Pi();
      deltaEta = deltaEta < TMath::Pi() ? deltaEta : deltaEta-2*TMath::Pi();

      Double_t r_i = sqrt(pow(deltaPhi,2) + pow(deltaEta,2));
      // radial moment
      radialMoment += constitMom.Pt()/jetMom.Pt() * r_i;
      radMom_riH->Fill(r_i);
      if (r_i > jetRadius*1.4 || r_i < 0.) {
        cout<<"\n\nWARNING: r_i in radial moment > 1.4*jetRadius or < 0\n"; cout<<r_i<<endl<<"\tphiJet, phiPart, etaJet, etaPart: "<< jet->Phi()<<" "<<constit->Phi()<<"\t"<< jet->Eta()<<" "<<constit->Eta()<<endl;
      }

      // angularity
      Double_t a_angu = 1.; // parameter
      Double_t theta = TMath::Pi()*r_i/(2*jetRadius);
      Double_t fac = pow(sin(theta), a_angu) * pow(1-cos(theta), 1-a_angu);
      angularity += constitMom.Mag() * fac;  // energy instead of momentum should be used !!!

      //vertex
      //Double_t vx = abs(jetp->Xv());
      //Double_t vy = abs(jetp->Yv());

      //Double_t Vx[1] = {vx};
      //vrtxHsparse->Fill(Vx);


    }  // end of loop over jet constituents

    radMomH->Fill(radialMoment);
    angularH->Fill(angularity);
    eInH->Fill(eInJet);
    if (ptRelMax > 0) ptRelH->Fill(ptRelMax); //ptRelMax = -1 means no electron in jet


    Int_t nSecV = jet->GetNumbersOfSecVertices();
    Int_t nAccSecV = 0;
    Double_t secVMaxR = -1;
    AliBasicJetSecondaryVertex* secV = new AliBasicJetSecondaryVertex();
    // loop over secondary vertices
    for(int iSecV=0; iSecV < nSecV; iSecV++){
      secV = jet->GetSecondaryVertex(iSecV);
      // cout<<"secV: (Vx,Vy,Vz)=("<<secV->Vx()<<","<<secV->Vy()<<","<<secV->Vz()<<
      //       ")\t|  Lxy="<<secV->Lxy()<<"+/-"<<secV->SigmaLxy()<<
      //       " \t|  Dispersion="<<secV->Dispersion()<<endl;
      Double_t R = abs(secV->Lxy());
      Double_t sigmaR = secV->SigmaLxy();
      ratioRsigmaRH->Fill(R/sigmaR);
      if(R > 3*sigmaR){
        SVrH->Fill(R);
        secVMaxR = max(R, secVMaxR);
        nAccSecV++;
      }
    }

    //count sec vertices with min X particles
    // Int_t Nsv = 0;
    // Double_t RsvMax = -1;
    // for(int i=0; i<secVertices.GetEntries(); ++i){
    //   Int_t n_tracks = secVcounts[i];
    //   if (n_tracks > 0) NsecV[0]+=1;
    //   if (n_tracks > 1) NsecV[1]+=1;
    //   if (n_tracks > 2) NsecV[2]+=1;
    //
    //
    //   if (secVcounts[i] < 3) continue;
    //   // sigma_vtx as a function of n_trakcs
    //   // http://personalpages.to.infn.it/~milano/presentations/LeonardoMilano_LHCBeyond.pdf s.29
    //   // probably better ref. https://cds.cern.ch/record/1648854/files/CERN-PH-EP-2014-031.pdf but for 7 TeV, not 0.9
    //   //Double_t sigma_vtx = sqrt( pow(166,2) + pow(540/sqrt(pow(n_tracks,0.9)),2) ) *1e-3; // um->mm, 0.9 TeV
    //   Double_t sigma_vtx = sqrt( pow(40,2) + pow(297/sqrt(pow(n_tracks,1.37)),2) ) *1e-3; // um->mm, 7 TeV
    //   TVector3 sv = *(TVector3*)secVertices.At(i);
    //   Double_t R = sqrt(sv.x()*sv.x() + sv.y()*sv.y());
    //   if (R > 3*sigma_vtx*sqrt(2)){
    //     SVrH->Fill(R);
    //     RsvMax = max(R, RsvMax);
    //     Nsv++;
    //   }
    // }
    //
    // //NsecVertices1H->Fill(NsecV[0]);
    // //NsecVertices2H->Fill(NsecV[1]);
    // //NsecVertices3H->Fill(NsecV[2]);
    SVnH->Fill(nAccSecV);

    // TODO #3 - tagging


    // observables struct filling
    jetObserv.pT       = jet->Pt();
    jetObserv.phi      = jet->Phi();
    jetObserv.eta      = jet->Eta();
    jetObserv.mult     = jet->GetNumbersOfConstituents();
    jetObserv.eIn      = eInJet;
    //jetObserv.eMotherPdg = eMother;
    jetObserv.svN      = nAccSecV;
    jetObserv.svR      = secVMaxR;
    jetObserv.ptRel    = ptRelMax;
    jetObserv.radMom   = radialMoment;
    jetObserv.angular  = angularity;
    // jetObserv.tagTrueFirst = tagTrueFirst;
    // jetObserv.tagTrueLast = tagTrueLast;
    // jetObserv.tagExp = tagExp;
    //
    // jetObserv.tagTrueFirstC = tagTrueFirstC;
    // jetObserv.tagTrueLastC = tagTrueLastC;
    // jetObserv.tagExpC = tagExpC;
    // jetObserv.tagTrueFirstB = tagTrueFirstB;
    // jetObserv.tagTrueLastB = tagTrueLastB;
    // jetObserv.tagExpB = tagExpB;

    // Filling tree
    jetObservTree->Fill();


  }  // end of loop over jets


  // Finish
  flavourDir->cd();
  //jetObservFile->cd();
  //tree->Print();
  jetObservTree->Write();


  TString* histFileName = new TString(pathToFile);
  *histFileName += "/histos.root";
  TFile* histFile = new TFile(histFileName->Data(),"UPDATE");
  histFile->cd();
  TDirectory* flavourHistDir = histFile->mkdir(flavour);
  flavourHistDir->cd();


  thetaH->Write();
  phiH->Write();
  etaH->Write();
  yH->Write();
  eH->Write();
  eetaH->Write();
  ptH->Write();
  //pdgH->Write();

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
  ratioRsigmaRH->Write();

  SVnH->Write();
  SVrH->Write();

  h7->Write();
  h16->Write();
  h17->Write();
  h28->Write();

  hJetEtaPhiAll->Write();

  eMotherPdgH->Write();
  isPhysPrimH->Write();
  isPhysPrimFJH->Write();



  jetObservFile->Close();
  histFile->Close();
  inputFile->Close();
}







void calc_observables(const char* pathToFile=""){
  TString* jetObservFileName = new TString(pathToFile);
  *jetObservFileName += "/jetObserv.root";
  TFile* jetObservFile = new TFile(jetObservFileName->Data(),"RECREATE");
  jetObservFile->Close();

  TString* histFileName = new TString(pathToFile);
  *histFileName += "/histos.root";
  TFile* histFile = new TFile(histFileName->Data(),"RECREATE");
  histFile->Close();

  calc_observables_flavour(pathToFile, "light");
  calc_observables_flavour(pathToFile, "b");
  calc_observables_flavour(pathToFile, "c");

}
