/* * * 
script executing generation.C with following parameters:
	* nev           - number of events
	* type          - 1 for gluon-gluon and 2 for quark-antiquark pair
	* parton_en     - energy of quark or gluon in GeV
	* q_id 			- quark type: 1=u, 2=d, 3=s, 4=c, 5=b, 6=t
	* longFileName  - include params above in output files names

Example of usage:
	$> aliroot 'run_gener.C(1000,1,50,5)'
or to run with reconstruction at once: use run_gener_recon.C.
Consider also using multirun_gener.sh

* * */

void run_gener(Int_t nev=300, Int_t type=1, Double_t parton_en=20., Int_t q_id=2, const char* longFileName="kTRUE"){
  TStopwatch timer;
  timer.Start();

  gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_ROOT -I$ALICE_ROOT/EVGEN -I$HOME/alice/ali-master/AliRoot/PYTHIA8/pythia8205/include/Pythia8 -I$HOME/alice/ali-master/AliRoot/PYTHIA8/pythia8205/include -I$ALICE_PHYSICS/include -I$FASTJET_ROOT/include");
  gSystem->Load("liblhapdf");      // Parton density functions
  gSystem->Load("libEGPythia6");   // TGenerator interface
  gSystem->Load("libpythia6");     // Pythia
  gSystem->Load("libAliPythia6");  // ALICE specific implementations


  TString* fileNameParams = new TString("");
  if(longFileName == "kTRUE"){
	if (type == 2){ *fileNameParams += "_g"; }  
	else if (type == 1) {*fileNameParams += "_q"; *fileNameParams += q_id;}
	else {*fileNameParams += "_X"; *fileNameParams += "X";}
	*fileNameParams += "_en";  *fileNameParams += parton_en;
	*fileNameParams += "_nev"; *fileNameParams += nev;
  }


  TString* command = new TString("generation.C+("); 
		 *command += nev;		*command += ","; 
		 *command += type;		*command += ","; 
		 *command += parton_en;	*command += ",";
		 *command += q_id;		*command += ",";
		 *command+="\""; *command+=fileNameParams->Data(); *command+="\"";
		 *command += ")";
  gROOT->Macro(*command);

  timer.Stop();
  timer.Print();
}
