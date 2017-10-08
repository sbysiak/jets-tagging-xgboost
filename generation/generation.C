/* * * 

Macro for generating specific events: 
pair of gluons or quark-antiquark in back-to-back configuration in random direction
which hadronize and form jets.
Used for comparision between jets arising from different partons.

Macro uses Pythia 8 generator (accessible through AliPythia).

It can be run using dedeicated macros: 
(multi)run_gener.C 
or for executing also reconstruction:
(multi)run_gener_recon.C 

Parameters:
	* nev            - number of events
	* type           - 1 for gluon-gluon and 2 for quark-antiquark pair
	* parton_en      - energy of quark or gluon in GeV
	* q_id 			 - quark type: 1=u, 2=d, 3=s, 4=c, 5=b, 6=t
	* fileNameParams - string containing above params, handled by (multi)run_(...).C, used for distinguish outputs

* * */



#include <Riostream.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TStopwatch.h>
#include <TDatime.h>
#include <TRandom.h>
#include <TDatabasePDG.h>
#include <TParticle.h>
#include <TArrayI.h>
#include <TTree.h>
#include <TString.h>
#include <TProfile.h>
#include <THnSparse.h>
#include <TTimeStamp.h>
#include <TMath.h>
#include <TSystem.h>
#include <TRandom3.h>
#include <TCanvas.h>
#include <TFile.h>

#include "AliGenerator.h"
#include "AliPDG.h"
#include "AliRunLoader.h"
#include "AliRun.h"
#include "AliStack.h"
#include "AliHeader.h"
#include "AliGenPythia.h"	//### E-TCint::AutoLoadCallback: failure loading library libAliPythia6 for class AliGenPythia

#include "AliGenPythiaPlus.h"
#include "AliPythia.h"
#include "AliPythia8.h"

#include <string>
#include <iostream>
#include <fstream>



// must be outside function due to logs implementation
// in order to use fileNameParams, which is generation() parameter
// need to either: (1) put log functions inside generation()
// (2) write logger class
const char* outFileName = "output_generation.root";  // used for logs and plots

void initLog();
void writeLog(const char *, Double_t, Double_t, Bool_t printOnScreen=kTRUE);
void writeLog(const char *,                     Bool_t printOnScreen=kTRUE);
void writeLog(const char *, Double_t ,          Bool_t printOnScreen=kTRUE);
void writeLog(const char *, const char *,       Bool_t printOnScreen=kTRUE);
void writeLine();

AliGenPythiaPlus*  CreateGenerator();
void fillPartons(int, double, int, Pythia8::Event&, Pythia8::ParticleData&, TRandom3 rand);
Float_t qMass(Int_t pdgCode);





// main function
void generation(Int_t nev = 2000, Int_t type = 2, Double_t parton_e = 20., Int_t q_id=2, const char* fileNameParams=""){

  TString* kineFileName = new TString("Kinematics");
  *kineFileName += fileNameParams;
  *kineFileName += ".root";

  TString* galiceFileName = new TString("galice");
  *galiceFileName += fileNameParams;
  *galiceFileName += ".root";

  string processes[2] = {"q qbar", "g g"};				


  initLog();
  writeLog("Macro run with following parameters:", kFALSE);
  writeLog("Process number:     ", type, kFALSE);
  writeLog("Process type:       ", (processes[type]).c_str(), kFALSE);
  writeLog("Number of events:   ", nev, kFALSE);
  writeLog("Parton energy:      ", parton_e, kFALSE);
  writeLog("Output File Name:   ", outFileName, kFALSE);
  writeLog("File name params:   ", fileNameParams, kFALSE);


  TProfile* meanPtH = new TProfile("meanPtH", "", 100, -0.5, 99.5);

  Int_t multNbins = 152;  
  Double_t binsMult[153];
  for (Int_t i=0; i<=multNbins; i++) { binsMult[i] = -0.5 + i; }
  binsMult[152] = 1000.;
  const Int_t ptNbins = 81;
  Double_t binsPt[82] = {0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 40.0, 45.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0, 140.0, 150.0, 160.0, 180.0, 200.0};
  const Int_t etaNbins = 140;
  Double_t binsEta[141] = {-7,-6.9,-6.8,-6.7,-6.6,-6.5,-6.4,-6.3,-6.2,-6.1,-6,-5.9,-5.8,-5.7,-5.6,-5.5,-5.4,-5.3,-5.2,-5.1,-5,-4.9,-4.8,-4.7,-4.6,-4.5,-4.4,-4.3,-4.2,-4.1,-4,-3.9,-3.8,-3.7,-3.6,-3.5,-3.4,-3.3,-3.2,-3.1,-3,-2.9,-2.8,-2.7,-2.6,-2.5,-2.4,-2.3,-2.2,-2.1,-2,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9,5,5.1,5.2,5.3,5.4,5.5,5.6,5.7,5.8,5.9,6,6.1,6.2,6.3,6.4,6.5,6.6,6.7,6.8,6.9,7};
  const Int_t zvNbins = 12;
  Double_t binsZv[13] = {-30.,-25.,-20.,-15.,-10.,-5.,0.,5.,10.,15.,20.,25.,30.};  
  Int_t binsMCPrimTrackHist[4]={zvNbins,ptNbins,etaNbins,multNbins};
  
  TH1D* spectraPt = new TH1D("spectraPt","",ptNbins,binsPt);
  spectraPt->Sumw2();
  TH2D* h2PtMult = new TH2D("h2PtMult","Pt:Mult",ptNbins,binsPt,multNbins,binsMult);
  h2PtMult->GetXaxis()->SetTitle("#p_{T}");
  h2PtMult->GetYaxis()->SetTitle("mult");
  h2PtMult->Sumw2();
  TH2D* h2PtEta = new TH2D("h2PtEta","Pt:Eta",ptNbins,binsPt,etaNbins,binsEta);
  h2PtEta->GetXaxis()->SetTitle("#p_{T}");
  h2PtEta->GetYaxis()->SetTitle("#eta");
  h2PtEta->Sumw2();

  // without pt equally divided binning
  TH2D *h2EtaPhi=new TH2D("h2EtaPhi","#eta vs #phi",100,0,6.5,100,-5,5);
  h2EtaPhi->GetXaxis()->SetTitle("#phi");
  h2EtaPhi->GetYaxis()->SetTitle("#eta");
  h2EtaPhi->Sumw2();
  TH2D *h2ThetaPhi=new TH2D("h2ThetaPhi","#theta vs #phi",100,0,6.5, 100,0.,3.2);
  h2ThetaPhi->GetXaxis()->SetTitle("#phi");
  h2ThetaPhi->GetYaxis()->SetTitle("#theta");  
  h2ThetaPhi->Sumw2();

  THnSparseD* hSMCTrackHist = new THnSparseD("hSMCTrackHist","Zv:mcpT:mceta:multTrue",4,binsMCPrimTrackHist); 
  hSMCTrackHist->SetBinEdges(0,binsZv);
  hSMCTrackHist->SetBinEdges(1,binsPt);
  hSMCTrackHist->SetBinEdges(2,binsEta);
  hSMCTrackHist->SetBinEdges(3,binsMult);
  hSMCTrackHist->GetAxis(0)->SetTitle("Zv (cm)");
  hSMCTrackHist->GetAxis(1)->SetTitle("MC p_{T} (GeV/c)");
  hSMCTrackHist->GetAxis(2)->SetTitle("#eta (MC)");
  hSMCTrackHist->GetAxis(3)->SetTitle("true multiplicity (MC) in |eta|<0.8");  

  TH1D* h1MCEvents = new TH1D("h1MCEvents","h1MCEvents",multNbins,binsMult);
  h1MCEvents->GetXaxis()->SetTitle("multMC true in |eta|<1.0");
  TH1D* h1MCEvents08 = new TH1D("h1MCEvents08","h1MCEvents08",multNbins,binsMult);
  h1MCEvents08->GetXaxis()->SetTitle("multMC true in |eta|<0.8");
  TH1D* h1MCEvents_alice = new TH1D("h1MCEvents_alice","h1MCEvents_alice",multNbins,binsMult);
  h1MCEvents_alice->GetXaxis()->SetTitle("multMC true in |eta|<0.8 && p_{T}>0.15");
  TH1D* h1MCEvents_All = new TH1D("h1MCEvents_All","h1MCEvents_All",3,0,3);
  h1MCEvents_All->GetXaxis()->SetTitle("multMC true");
  
  gSystem->Load("liblhapdf.so");    
  gSystem->Load("libEGPythia6.so"); 
  gSystem->Load("libpythia6.so");   
  gSystem->Load("libAliPythia6.so");
  gSystem->Load("libpythia8.so");   
  gSystem->Load("libAliPythia8.so");
//  gSystem->Setenv("PYTHIA8DATA", gSystem->ExpandPathName("$ALICE_ROOT/PYTHIA8/pythia8175/xmldoc"));
  gSystem->Setenv("PYTHIA8DATA", gSystem->ExpandPathName("$ALICE_ROOT/PYTHIA8/pythia8/xmldoc"));
  gSystem->Setenv("LHAPDF",      gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF"));
  gSystem->Setenv("LHAPATH",     gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF/PDFsets"));

  writeLog("Environment setting                  DONE", kFALSE);



  //
  //                    Ali Classes initialization
  //

  using namespace Pythia8;

  //  RunLoader
  AliRunLoader* rl = AliRunLoader::Open(galiceFileName->Data(),"FASTRUN","recreate");
  rl->SetKineFileName(kineFileName->Data());
  rl->SetCompressionLevel(2);
  rl->SetNumberOfEventsPerFile(nev);
  rl->LoadKinematics("RECREATE");
  rl->MakeTree("E");
  gAlice->SetRunLoader(rl);
  //  Create stack
  rl->MakeStack();
  AliStack* stack = rl->Stack();

  
  //  Header
  AliHeader* header = rl->GetHeader();
  
  writeLine();
  writeLog("PYTHIA8 run with following user settings:", kFALSE);
  //  Create and Initialize Generator
  AliGenPythiaPlus *gener = CreateGenerator();

  (AliPythia8::Instance())->Pythia8()->readString("ProcessLevel:all = off");  // !!!
  writeLog("PYTHIA setting:   ProcessLevel:all = off", kFALSE);
  // Switch off automatic event listing in favour of manual.
	//  (AliPythia8::Instance())->Pythia8()->readString("Next:numberShowInfo = 0");
	//  (AliPythia8::Instance())->Pythia8()->readString("Next:numberShowProcess = 0");
	//  (AliPythia8::Instance())->Pythia8()->readString("Next:numberShowEvent = 0");
	//  writeLog("PYTHIA setting:   Next:numberShowInfo = 0", kFALSE);
	//  writeLog("PYTHIA setting:   Next:numberShowProcess = 0", kFALSE);  
	//  writeLog("PYTHIA setting:   Next:numberShowEvent = 0", kFALSE);
  gener->Init();
  gener->SetStack(stack);
  writeLog("Ali classes initialization           DONE");



  double ee = parton_e; // Set typical energy per parton.

  Event& event      = (AliPythia8::Instance())->Pythia8()->event;
  ParticleData& pdt = (AliPythia8::Instance())->Pythia8()->particleData;
  TRandom3 rand     = TRandom3(0);


  //
  //                    Events Loop
  //

  for (Int_t i_ev = 0; i_ev < nev; i_ev++) {  // >> Event Loop
    if (!(i_ev%100)) cout<<"Event number: "<<i_ev<<" of "<<nev<<endl;
    //  Initialize event
    header->Reset(0,i_ev);
    rl->SetEventNumber(i_ev);
    stack->Reset();
    rl->MakeTree("K");
    //stack->ConnectTree();

    fillPartons( type, ee, q_id, event, pdt, rand );

    //  Generate event
    gener->Generate();

    //
    // edgars' part
    //
    Int_t multAll=0; // total mult
    Int_t multEta10=0; // mult in eta +-1.0
    Int_t multEta08=0; // mult in eta +-0.8
    Int_t multEta_alice=0;
    
    //  Analysis
    Int_t npart = stack->GetNprimary();
    Float_t nch   = 0.;
    Float_t sumPt = 0.;
    // Calculate number of charged particles
    for (Int_t i_part = 0; i_part < npart; i_part++) {
      TParticle *particle = stack->Particle(i_part);
      Int_t pdg = particle->GetPdgCode();
      if (!(stack->IsPhysicalPrimary(i_part)))   continue;
      if (particle->GetPDG()->Charge() == 0)   continue;
      //if (TMath::Abs(particle->Eta()) > 1)   continue;
      if (particle->Pt() < 0.)                 continue;
      nch   += 1.;
      if ( TMath::Abs(particle->Eta()) < 1.0 ) multEta10++;
      if ( TMath::Abs(particle->Eta()) < 0.8 ) multEta08++;
      if ( TMath::Abs(particle->Eta()) < 0.8 && particle->Pt()>0.15) multEta_alice++;
      //if (particle->Energy() < 4.) stack->fParticles.RemoveAt(i_part);
      //if (particle->Energy() > 4.) stack->FlagTrack(i_part);
		
    }

    h1MCEvents_All->Fill(0); //Multiplicity of the events without any eta restriction
    if(multEta10>0) h1MCEvents_All->Fill(1);
    if(multEta08>0) h1MCEvents_All->Fill(2);
    
    h1MCEvents->Fill(multEta10);
    h1MCEvents08->Fill(multEta08);
    h1MCEvents_alice->Fill(multEta_alice);

    for (Int_t i_part = 0; i_part < npart; i_part++) {
      TParticle *particle = stack->Particle(i_part);
      Int_t pdg = particle->GetPdgCode();
      if (!(stack->IsPhysicalPrimary(i_part)))   continue;
      if (particle->GetPDG()->Charge() == 0)   continue;
      if (particle->Pt() <= 0)               continue;

      Float_t pt = particle->Pt();
      Float_t eta = particle->Eta();
      Float_t phi = particle->Phi();
      Float_t theta = particle->Theta();
      
      sumPt += pt;
      h2PtMult->Fill(pt,nch);
      spectraPt->Fill(pt);
      h2PtEta->Fill(pt,eta);
      Double_t coords[4] = {0.,pt,eta,double(multEta08)};
      hSMCTrackHist->Fill(coords);
	  h2EtaPhi->Fill(phi, eta);
	  h2ThetaPhi->Fill(phi, theta);
    }
    if (nch >= 1.) {
      sumPt /= nch;
      meanPtH->Fill(nch, sumPt);
    }

    //  Finish event
    header->SetNprimary(stack->GetNprimary());
    header->SetNtrack(stack->GetNtrack());  
    //      I/O
    //	

    stack->FinishEvent();
    header->SetStack(stack);

    rl->TreeE()->Fill();
    rl->WriteKinematics("OVERWRITE");
 
  } // event loop



  //
  //                         Termination
  //

  //  Generator
  gener->FinishRun();
  //  Write file
  rl->WriteHeader("OVERWRITE");
  gener->Write();
  rl->Write();

  // Write histos
  TFile* fout = TFile::Open(outFileName,"RECREATE");
/*
  hSMCTrackHist->Write();
  h1MCEvents->Write();
  h1MCEvents08->Write();
  h1MCEvents_alice->Write();
  h1MCEvents_All->Write();
  spectraPt->Write();
*/
  TCanvas* cc = new TCanvas("cc","Canvas cc",400,10,600,700);
  cc->Divide(1,2);
  cc->cd(1);
  h2EtaPhi->Draw();
  cc->cd(2);
  h2ThetaPhi->Draw();

  h2EtaPhi->Write();
  h2ThetaPhi->Write();
  spectraPt->Write();
  h2PtMult->Write();
  h2PtEta->Write();
  h2EtaPhi->Write();
  h2ThetaPhi->Write();


  delete rl;

  //AliPythia8::Instance()->Pythia8()->stat();

} // end of Pythia8(...)







AliGenPythiaPlus*  CreateGenerator(){
  AliGenPythiaPlus* gener = new AliGenPythiaPlus(AliPythia8::Instance());

    TDatime dt;
    UInt_t seed = dt.Get();
    cout     <<"Random seed set to "<< seed << endl;

  (AliPythia8::Instance())->ReadString("Random:setSeed = on");
  (AliPythia8::Instance())->ReadString(Form("Random:seed = %u", seed%900000000)); 

  gener->SetProcess(kPyMbDefault);

  //   Centre of mass energy 
  //   gener->SetEnergyCMS(energy); // was used in thesis

  return gener;
}




void fillPartons(int type, double ee, int q_id, Pythia8::Event& event, Pythia8::ParticleData& pdt, TRandom3 rndm) {
  // Reset event record to allow for new event.

  event.reset();
  rndm.SetSeed();

  double pi = TMath::Pi();
  double theta = 40.395+(180.-2*40.395)*rndm.Rndm(); //theta: 40.395 - 139.605  <=>  Eta: -1 - 1
  double cTheta = cos(theta);
  double sTheta = sin(theta);
  double phi = 2*pi*rndm.Rndm();
  double cPhi = cos(phi);
  double sPhi = sin(phi);

  // Information on a q qbar system, to be hadronized.
  if (type == 1) {
    int    id = q_id;
    double mm = qMass(id); // http://home.thep.lu.se/~torbjorn/pythia81html/ParticleData.html
    double pp = Pythia8::sqrtpos(ee*ee - mm*mm);

    event.append(  id, 23, 101,   0,  sPhi*sTheta*pp,  cPhi*sTheta*pp,  cTheta*pp, ee, mm);
    event.append( -id, 23,   0, 101, -sPhi*sTheta*pp, -cPhi*sTheta*pp, -cTheta*pp, ee, mm);

  // Information on a g g system, to be hadronized.
  } else if (type == 2) {

    event.append( 21, 23, 102, 103,  sPhi*sTheta*ee,  cPhi*sTheta*ee,  cTheta*ee, ee);
    event.append( 21, 23, 103, 102, -sPhi*sTheta*ee, -cPhi*sTheta*ee, -cTheta*ee, ee);
  
  } else{
	writeLog("\nERROR: unhandled type: ", double(type));
    exit (EXIT_FAILURE);
  } 
}




Float_t qMass(Int_t pdgCode){
	// implemented due to error while executing > pdt.m0(id) <
	// masses comes from http://home.thep.lu.se/~torbjorn/pythia81html/ParticleData.html
	// so it should not change execution of program

	if (pdgCode > 6 || pdgCode < 1) {writeLog("\nERROR: given pdg is not quark's: ", double(pdgCode)); exit (EXIT_FAILURE);}
	else if (pdgCode == 1) return 0.33000;
	else if (pdgCode == 2) return 0.33000;
	else if (pdgCode == 3) return 0.50000;
	else if (pdgCode == 4) return 1.50000;
	else if (pdgCode == 5) return 4.80000;
	else if (pdgCode == 6) {writeLog("\nWARNING: you are using top quark with not well definied mass ! \n", double(pdgCode)); return 171.00000;}
}




//________________________________________________________
//
//  						 L O G S
//________________________________________________________


void initLog(){   
  // cleaning file
  std::string logFile(outFileName);
  std::string stroot(".root");
  logFile.replace(logFile.find(stroot.data()), stroot.length(), ".log");
  std::ofstream myfile;
  myfile.open(logFile.data());
  myfile.close();

  TTimeStamp time;
  string datetime = time.AsString("cl");
  writeLog("         Date:            ",(datetime.substr(0,10)).c_str(), kFALSE);
  writeLog("         Starting Time:   ",(datetime.substr(11,8)).c_str(), kFALSE);
  writeLine();
}

void writeLog(const char *addToLog, Double_t number1, Double_t number2, Bool_t printOnScreen){
  std::string logFile(outFileName);
  std::string stroot(".root");
  logFile.replace(logFile.find(stroot.data()), stroot.length(), ".log");

  std::ofstream myfile;
  myfile.open(logFile.data(),std::ios_base::app); 

  if(number1==9999.){
    myfile << addToLog<< "\n";
    if (printOnScreen) cout << addToLog<< "\n";
  }else if(number2==9999.){
    myfile << addToLog<< number1 << "\n";
    if (printOnScreen) cout << addToLog<< number1 << "\n";
  }else{
    myfile << addToLog<< number1 << " " << number2 << "\n";
    if (printOnScreen) cout << addToLog<< number1 << " " << number2 << "\n";
  }

  myfile.close();
}

void writeLog(const char *addToLog, Bool_t printOnScreen){
    writeLog(addToLog, 9999., 9999., printOnScreen);
}

void writeLog(const char *addToLog, Double_t number1, Bool_t printOnScreen){
    writeLog(addToLog, number1, 9999., printOnScreen );
}

void writeLog(const char *addToLog, const char *alsoadd, Bool_t printOnScreen){
  std::string logFile(outFileName);
  std::string stroot(".root");
  logFile.replace(logFile.find(stroot.data()), stroot.length(), ".log");

  std::ofstream myfile;
  myfile.open(logFile.data(),std::ios_base::app);

  myfile << addToLog<< alsoadd << "\n";
  if (printOnScreen) cout << addToLog << alsoadd << "\n";
}

void writeLine(){
  writeLog("________________________________________________________", kFALSE);
  writeLog("", kFALSE);
}

