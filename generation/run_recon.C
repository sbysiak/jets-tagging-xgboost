/* * * 
script executing reconstruction.C with following parameters:
	* nev           - number of events
	* type          - 1 for gluon-gluon and 2 for quark-antiquark pair
	* parton_en     - energy of quark or gluon in GeV
	* q_id 					- quark type: 1=u, 2=d, 3=s, 4=c, 5=b, 6=t

	Above parameters defines which input files (Kinematics_(...).root, galice_(...).root) should be used

Example of usage:
	$> aliroot 'run_recon.C(1000,1,50,5)'
or to run with generation at once: use run_gener_recon.C.
Consider also using multirun_recon.sh

* * */

void run_recon(Int_t nev=200, Int_t type=1, Double_t parton_en=20., Int_t q_id=2){
  // Simulation and reconstruction
  TStopwatch timer;
  timer.Start();
  gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_ROOT -I$ALICE_ROOT/EVGEN -I$ALICE_ROOT/EVE/macros -I$ALICE_PHYSICS/include -I$FASTJET_ROOT/include");


  gROOT->ProcessLine(".L loadlibs.C");
  loadlibs();

  TString* fileNameParams = new TString("");

	if (type == 2){ *fileNameParams += "_g"; }  
	else if (type == 1) {*fileNameParams += "_q"; *fileNameParams += q_id;}
	else {*fileNameParams += "_X"; *fileNameParams += "X";}
	*fileNameParams += "_en";  *fileNameParams += parton_en;
	*fileNameParams += "_nev"; *fileNameParams += nev;


  //n = nev;  // Use the requested number of events
  //gROOT->Macro("readPrimariesFastJet.C+(n1, n2)");
  TString* ts2 = new TString("reconstruction.C+("); 
			*ts2+=nev; 	*ts2 += ",";	
			*ts2+="\""; *ts2+=fileNameParams->Data(); *ts2+="\"";
			*ts2+=")";

  gROOT->Macro(*ts2);
  timer.Stop();
  timer.Print();
}
