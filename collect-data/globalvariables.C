{

//AliTrackContainer::SetDefTrackCutsPeriod("lhc11h"); // for correct hybrid track selection/track tagging
Double_t kTrackPtCut = 0.15;
Double_t kGhostArea = 0.005;
Double_t kJetPtCut = 1;

TString kRandomizedTracksName = "tracks_randomized";
TString kTracksName = "tracks";
TString kTpcKtJetsR020Name(Form("Jet_KTChargedR020_%s_pT0150_pt_scheme",kTracksName.Data())); 

TString kTpcJetsR020Name(Form("Jet_AKTChargedR020_%s_pT0150_pt_scheme",kTracksName.Data()));
TString kTpcJetsR030Name(Form("Jet_AKTChargedR030_%s_pT0150_pt_scheme",kTracksName.Data()));
TString kTpcJetsR040Name(Form("Jet_AKTChargedR040_%s_pT0150_pt_scheme",kTracksName.Data()));


// Physics selection
UInt_t kComPhysSel = AliVEvent::kAny;
UInt_t kPhysSel = AliVEvent::kMB + AliVEvent::kSemiCentral + AliVEvent::kCentral;

}
