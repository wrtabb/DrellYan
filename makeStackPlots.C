
vector<TString> file_data = {
	"output_data/DYHists_v2p6_Data_RunB_EE.root",
	"output_data/DYHists_v2p6_Data_RunC_EE.root",
	"output_data/DYHists_v2p6_Data_RunD_EE.root",
	"output_data/DYHists_v2p6_Data_RunE_EE.root",
	"output_data/DYHists_v2p6_Data_RunF_EE.root",
	"output_data/DYHists_v2p6_Data_RunG_EE.root",
	"output_data/DYHists_v2p6_Data_RunHver2_EE.root",
	"output_data/DYHists_v2p6_Data_RunHver3_EE.root"
};

vector<TString> file_diboson = {
	"output_data/DYHists_v2p6_Diboson_WW_EE.root",
	"output_data/DYHists_v2p6_Diboson_WZ_EE.root",
	"output_data/DYHists_v2p6_Diboson_ZZ_EE.root"
};

vector<TString> file_taus = {
	"output_data/DYHists_v2p6_TauTau_DYLL_M100to200_TauTau_EE.root",
	"output_data/DYHists_v2p6_TauTau_DYLL_M200to400_TauTau_EE.root",
	"output_data/DYHists_v2p6_TauTau_DYLL_M400to500_TauTau_EE.root",
	"output_data/DYHists_v2p6_TauTau_DYLL_M500to700_TauTau_EE.root",
	"output_data/DYHists_v2p6_TauTau_DYLL_M700to800_TauTau_EE.root",
	"output_data/DYHists_v2p6_TauTau_DYLL_M1000to1500_TauTau_EE.root",
	"output_data/DYHists_v2p6_TauTau_DYLL_M1500to2000_TauTau_EE.root",
	"output_data/DYHists_v2p6_TauTau_DYLL_M2000to3000_TauTau_EE.root"
};

vector<TString> file_DYLL = {
	"output_data/DYHists_v2p6_DYtoLL_M1000to1500_EE.root",
	"output_data/DYHists_v2p6_DYtoLL_M100to200_EE.root",
	"output_data/DYHists_v2p6_DYtoLL_M10to50_EE.root",
	"output_data/DYHists_v2p6_DYtoLL_M1500to2000_EE.root",
	"output_data/DYHists_v2p6_DYtoLL_M2000to3000_EE.root",
	"output_data/DYHists_v2p6_DYtoLL_M200to400_EE.root",
	"output_data/DYHists_v2p6_DYtoLL_M400to500_EE.root",
	"output_data/DYHists_v2p6_DYtoLL_M500to700_EE.root",
	"output_data/DYHists_v2p6_DYtoLL_M50to100_EE.root",
	"output_data/DYHists_v2p6_DYtoLL_M700to800_EE.root",
	"output_data/DYHists_v2p6_DYtoLL_M800to1000_EE.root"
};

vector<TString> file_top = {
	"output_data/DYHists_v2p6_Top_ST_tW_EE.root",
	"output_data/DYHists_v2p6_Top_ST_tbarW_EE.root",
	"output_data/DYHists_v2p6_Top_ttbar_M0to700_EE.root",
	"output_data/DYHists_v2p6_Top_ttbar_M1000toInf_EE.root",
	"output_data/DYHists_v2p6_Top_ttbar_M700to1000_EE.root"
};

vector<TString> file_Fake = {
	"output_data/DYHists_v2p6_Fake_WJetsToLNuExt_EE.root",
	"output_data/DYHists_v2p6_Fake_WJetsToLNu_EE.root"
};

TH1D*GetHistogram(vector<TString> filesvector);

void makeStackPlots()
{
	gROOT->SetBatch(true);
	TCanvas*c1 = new TCanvas("c1","",0,0,1000,1000);
	c1->SetGrid();
	c1->SetLogx();
	c1->SetLogy();

	TH1D*hMassDYLL = GetHistogram(file_DYLL);
	TH1D*hMassTops = GetHistogram(file_top);
	TH1D*hMassFake = GetHistogram(file_Fake);
	TH1D*hMassDiboson = GetHistogram(file_diboson);
	//TH1D*hMassTau = GetHistogram(file_taus);

	THStack*hStack_Mass = new THStack("hStack_Mass","");
	hStack_Mass->Add(hMassTops);
	hStack_Mass->Add(hMassFake);
	hStack_Mass->Add(hMassDiboson);
	//hStack_Mass->Add(hMassTau);
	hStack_Mass->Add(hMassDYLL);
	
	TH1D*hMassData = GetHistogram(file_data);
	hStack_Mass->SetMinimum(1);;
	hStack_Mass->Draw("hist");
	hMassData->Draw("pe,same");

	c1->SaveAs("tempDataVsMCStack.png");
}

TH1D*GetHistogram(vector<TString> filesvector)
{
	TString histMassName = "hMassReco";
	int nFiles = filesvector.size();
	TFile*file_load[nFiles];
	TH1D*hMassArray[nFiles];
	TH1D*hMass;
	for(int i=0;i<nFiles;i++){
		file_load[i] = new TFile(filesvector.at(i));
		hMassArray[i] = (TH1D*)file_load[i]->Get(histMassName);
		if(i==0) hMass = (TH1D*)hMassArray[i]->Clone();
		else hMass->Add(hMassArray[i]);
	}

	return hMass;
}

