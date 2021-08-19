
void tempMakePlots()
{
	vector<TString> file_list = {
		"output_data/DYHists_v2p6_DYtoLL_DoubleEG_RunB_EE.root",
		"output_data/DYHists_v2p6_DYtoLL_DoubleEG_RunC_EE.root",
		"output_data/DYHists_v2p6_DYtoLL_DoubleEG_RunD_EE.root",
		"output_data/DYHists_v2p6_DYtoLL_DoubleEG_RunE_EE.root",
		"output_data/DYHists_v2p6_DYtoLL_DoubleEG_RunF_EE.root",
		"output_data/DYHists_v2p6_DYtoLL_DoubleEG_RunG_EE.root",
		"output_data/DYHists_v2p6_DYtoLL_DoubleEG_RunHver2_EE.root",
		"output_data/DYHists_v2p6_DYtoLL_DoubleEG_RunHver3_EE.root"
	};
	
	vector<TFile*> files;
	vector<TH1D*> histInvMass;
	vector<TH1D*> histRapidity;
	vector<TH1D*> histPt;
	int nFiles = file_list.size();
	for(int i=0;i<nFiles;i++){
		files.push_back(new TFile(file_list.at(i)));
		histInvMass.push_back((TH1D*)files.at(i)->Get("hMassReco"));
		histRapidity.push_back((TH1D*)files.at(i)->Get("hRapidityReco"));
		histPt.push_back((TH1D*)files.at(i)->Get("hPtReco"));
	}
	
	for(int i=1;i<nFiles;i++){
		cout << "Adding hist from: " << file_list.at(i) << endl;
		histInvMass.at(0)->Add(histInvMass.at(i));
		histRapidity.at(0)->Add(histRapidity.at(i));
		histPt.at(0)->Add(histPt.at(i));
	}	

	histInvMass.at(0)->GetXaxis()->SetTitle("invariant mass [GeV]");
	histInvMass.at(0)->GetXaxis()->SetNoExponent();
	histInvMass.at(0)->GetXaxis()->SetMoreLogLabels();
	histRapidity.at(0)->GetXaxis()->SetTitle("rapidity");
	histPt.at(0)->GetXaxis()->SetTitle("p_{T} [GeV]");

	TCanvas*c1 = new TCanvas("c1","",0,0,1000,1000);
	c1->SetGrid();
	c1->SetLogx();
	c1->SetLogy();
	histInvMass.at(0)->Draw("pe");
	c1->SaveAs("plots/DYtoEE_InvMassData.png");

	TCanvas*c2 = new TCanvas("c2","",0,0,1000,1000);
	c2->SetGrid();
	c2->SetLogy();
	histRapidity.at(0)->Draw("pe");
	c2->SaveAs("plots/DYtoEE_RapidityData.png");

	TCanvas*c3 = new TCanvas("c3","",0,0,1000,1000);
	c3->SetGrid();
	c3->SetLogy();
	histPt.at(0)->Draw("pe");
	c3->SaveAs("plots/DYtoEE_PtData.png");
}
