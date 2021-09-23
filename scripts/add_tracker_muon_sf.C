void add_tracker_muon_sf(TString input_filename="Ntuples/2018/T3MSelectionTreeInput_preapproval_v1_combined_2018UL_calo_check.root",
								 TString output_filename="Ntuples/2018/T3MSelectionTreeInput_preapproval_v1_combined_2018UL_calo_check_muonsf.root"
								 TString weight_filename = "Efficiencies/2018/NUM_Tau3MuTracker_DEN_genTracks_pt_abseta.root"){

  float tracker_muon_sf = 0.0;
  float pt;
  float eta;
  bool threeGlobal;
  
  TFile* weight_file = new TFile(weight_filename, "read");
  TH1F* h2 = (TH1F*)weight_file->Get("NUM_Tau3MuTracker_DEN_genTracks_pt_abseta");


  TString tree_list[4] = {"TreeB", "TreeS_Ds", "TreeS_Bu", "TreeS_Bd"};

  TFile* input_file = new TFile(input_filename, "read");
  TTree* old_tree[4];

  for (int i=0; i<4; i++){
    old_tree[i] = (TTree*)input_file->Get(tree_list[i]);
  }

  TFile* output_file = new TFile(output_filename, "recreate");
  TTree* new_tree[4];
  for (int i=0; i<4; i++) {
    new_tree[i] = old_tree[i]->CloneTree(0);
    new_tree[i]->Branch("tracker_muon_sf", &tracker_muon_sf);
    }

  for (int i=0; i<4; i++){
    int nentries = old_tree[i]->GetEntriesFast();
    old_tree[i]->SetBranchAddress("var_Muon3_Pt", &pt);
    old_tree[i]->SetBranchAddress("var_Muon3_Eta", &eta);
    old_tree[i]->SetBranchAddress("threeGlobal", &threeGlobal);

    for (int j=0; j<nentries; j++){
      old_tree[i]->GetEntry(j);
      int binx = h2->GetXaxis()->FindBin(pt);
      int biny = h2->GetYaxis()->FindBin(abs(eta));
      if (threeGlobal==0){
		tracker_muon_sf = h2->GetBinContent(binx, biny);
		if (tracker_muon_sf==0) cout<<pt<<" "<<eta<<endl;
		}
      else tracker_muon_sf = 0.0;
		if (i==0) tracker_muon_sf=1.0;
      new_tree[i]->Fill();
    }

  }

  output_file->cd();
  for (int i=0; i<4; i++) new_tree[i]->Write();
  output_file->Close();

  input_file->Close();
  weight_file->Close();

  return;
}
