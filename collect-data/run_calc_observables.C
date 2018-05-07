void run_calc_observables(const char* pathToFile="."){
  TStopwatch timer;
  timer.Start();
  gSystem->SetIncludePath("-I$ALICE_ROOT/include -I$ALICE_ROOT -I$ALICE_PHYSICS/include -I$ALICE_PHYSICS -I$FASTJET_ROOT/include -I/home/sebbys/alice/ali-master/AliPhysics/PWGHF/jetsHF");

  gROOT->ProcessLine(".L loadlibs.C");
  loadlibs();

  cout<<"libs loaded"<<endl;
  TString* command = new TString("calc_observables.C+(");
    *command+="\"";
    *command+=pathToFile;
    *command+="\"";
  *command+=")";

  cout<<"running macro"<<endl;
  gROOT->Macro(*command);
  timer.Stop();
  timer.Print();
}
