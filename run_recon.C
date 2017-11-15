void run_recon(Int_t nev=200, const char* pathToFile="DATA/LHC17f8g_20_255618_001/"){
  TStopwatch timer;
  timer.Start();
  gSystem->SetIncludePath("-I$ALICE_ROOT/include -I$ALICE_ROOT -I$ALICE_PHYSICS/include -I$ALICE_PHYSICS -I$FASTJET_ROOT/include");

  gROOT->ProcessLine(".L loadlibs.C");
  loadlibs();

  cout<<"libs loaded"<<endl;
  TString* ts2 = new TString("reconstruction.C+("); 
			*ts2+=nev; 	*ts2 += ",";	
			*ts2+="\""; *ts2+=pathToFile; *ts2+="\"";
			*ts2+=")";

  cout<<"running macro"<<endl;
  gROOT->Macro(*ts2);
  timer.Stop();
  timer.Print();
}
