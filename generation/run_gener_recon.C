/* * * 
script executing generation.C and then reconstruction.C with following parameters:
	* nev           - number of events
	* type          - 1 for gluon-gluon and 2 for quark-antiquark pair
	* parton_en     - energy of quark or gluon in GeV
	* q_id 			- quark type: 1=u, 2=d, 3=s, 4=c, 5=b, 6=t
	* longFileName  - include params above in output files names

Example of usage:
	$> aliroot 'run_gener_recon.C(1000,1,50,5)'
Consider also using multirun_gener_recon.sh

* * */

void rungen_gener_recon(Int_t nev=300, Int_t type=1, Double_t parton_en=20., Int_t q_id=2, const char* longFileName=kTRUE){

  // Simulation and reconstruction
  TStopwatch timer;
  timer.Start();

  gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_ROOT -I$ALICE_ROOT/EVGEN -I/home/sebbys/alice/ali-master/AliRoot/PYTHIA8/pythia8205/include/Pythia8 -I/home/sebbys/alice/ali-master/AliRoot/PYTHIA8/pythia8205/include -I$ALICE_PHYSICS/include -I$FASTJET_ROOT/include");
  gSystem->Load("liblhapdf");      // Parton density functions
  gSystem->Load("libEGPythia6");   // TGenerator interface
  gSystem->Load("libpythia6");     // Pythia
  gSystem->Load("libAliPythia6");  // ALICE specific implementations

  gROOT->ProcessLine(".L loadlibs.C");
  loadlibs();


  TString* fileNameParams = new TString("");
  if(longFileName == "kTRUE"){
	if (type == 2){ *fileNameParams += "_g"; }  
	else if (type == 1) {*fileNameParams += "_q"; *fileNameParams += q_id;}
	else {*fileNameParams += "_q"; *fileNameParams += "X";}
	*fileNameParams += "_en";  *fileNameParams += parton_en;
	*fileNameParams += "_nev"; *fileNameParams += nev;

  }

  TString* ts = new TString("generation.C+("); 
		 *ts += nev;		*ts += ","; 
		 *ts +=	type;		*ts += ","; 
		 *ts += parton_en;	*ts += ","; 
		 *ts += q_id;		*ts += ",";	
		 *ts+="\""; *ts+=fileNameParams->Data(); *ts+="\"";
		 *ts += ")";
  gROOT->Macro(*ts);


  for (int i=0;i<3;++i) cout<<"_______________________________________________\n";
  cout<<"\tG E N E R A T I O N    D O N E\n";
  for (int i=0;i<3;++i) cout<<"_______________________________________________\n";


  TString* ts2 = new TString("reconstruction.C+("); 
		*ts2+=nev; 	*ts2+=","; 
		*ts2+="\""; *ts2+=fileNameParams->Data(); *ts2+="\"";

  gROOT->Macro(*ts2);
  timer.Stop();
  timer.Print();
}
