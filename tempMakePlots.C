
void tempMakePlots()
{
	//gROOT->SetBatch(true);
	vector<TString> file_list = {

		"output_data/DYHists_v2p6_Data_RunB_EE.root",
		"output_data/DYHists_v2p6_Data_RunC_EE.root",
		"output_data/DYHists_v2p6_Data_RunD_EE.root",
		"output_data/DYHists_v2p6_Data_RunE_EE.root",
		"output_data/DYHists_v2p6_Data_RunF_EE.root",
		"output_data/DYHists_v2p6_Data_RunG_EE.root",
		"output_data/DYHists_v2p6_Data_RunHver2_EE.root",
		"output_data/DYHists_v2p6_Data_RunHver3_EE.root"
/*
		"output_data/DYHists_v2p6_DYtoLL_M10to50_EE.root",
		"output_data/DYHists_v2p6_DYtoLL_M50to100_EE.root",
		"output_data/DYHists_v2p6_DYtoLL_M100to200_EE.root",
		"output_data/DYHists_v2p6_DYtoLL_M200to400_EE.root",
		"output_data/DYHists_v2p6_DYtoLL_M400to500_EE.root",
		"output_data/DYHists_v2p6_DYtoLL_M500to700_EE.root",
		"output_data/DYHists_v2p6_DYtoLL_M700to800_EE.root",
		"output_data/DYHists_v2p6_DYtoLL_M800to1000_EE.root",
		"output_data/DYHists_v2p6_DYtoLL_M1000to1500_EE.root",
		"output_data/DYHists_v2p6_DYtoLL_M1500to2000_EE.root",
		"output_data/DYHists_v2p6_DYtoLL_M2000to3000_EE.root"
*/
	};
	
	vector<TFile*> files;
	vector<TH1D*> histInvMass;
	vector<TH1D*> histRapidity;
	vector<TH1D*> histPt;
	int nFiles = file_list.size();
	for(int i=0;i<nFiles;i++){
		cout << "Adding hist from: " << file_list.at(i) << endl;
		files.push_back(new TFile(file_list.at(i)));
		histInvMass.push_back((TH1D*)files.at(i)->Get("hMassReco"));
		histRapidity.push_back((TH1D*)files.at(i)->Get("hRapidityReco"));
		histPt.push_back((TH1D*)files.at(i)->Get("hPtReco"));
	}
	
	for(int i=1;i<nFiles;i++){
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
	histInvMass.at(0)->Draw("hist");
//	c1->SaveAs("plots/DYtoEE_InvMassData.png");

	TCanvas*c2 = new TCanvas("c2","",0,0,1000,1000);
	c2->SetGrid();
	c2->SetLogy();
	histRapidity.at(0)->Draw("hist");
//	c2->SaveAs("plots/DYtoEE_RapidityData.png");

	TCanvas*c3 = new TCanvas("c3","",0,0,1000,1000);
	c3->SetGrid();
	c3->SetLogy();
	histPt.at(0)->Draw("hist");
//	c3->SaveAs("plots/DYtoEE_PtData.png");
}
