using namespace DrellYanVariables;

enum LogType{
	LOG_X,
	LOG_Y,
	LOG_BOTH,
	LINEAR
};

TCanvas*MakeCanvas(TString cname,LogType logtype);

void makePlots(LepType lepType)
{
	TFile*load_file;
	TH1D*hMass;
	TH1D*hRapidity;
	TH1D*hPt;
	TString saveNameMass = "plots/InvMass";
	TString saveNameRapidity = "plots/Rapidity";
	TString saveNamePt = "plots/Pt";
	if(lepType==ELE){
		load_file = new TFile("data/DYHists_EE.root");
		hMass = (TH1D*)load_file->Get("histInvMass_EE");
		hRapidity = (TH1D*)load_file->Get("histRapidity_EE");
		hPt = (TH1D*)load_file->Get("histPt_EE");
		
		saveNameMass += "_EE.png";
		saveNameRapidity += "_EE.png";
		saveNamePt += "_EE.png";
	}
	else if(lepType==MUON){
		load_file = new TFile("data/DYHists_MuMu.root");
		hMass = (TH1D*)load_file->Get("histInvMass_MuMu");
		hRapidity = (TH1D*)load_file->Get("histRapidity_MuMu");
		hPt = (TH1D*)load_file->Get("histPt_MuMu");

		saveNameMass += "_MuMu.png";
		saveNameRapidity += "_MuMu.png";
		saveNamePt += "_MuMu.png";
	}
	hMass->SetFillColor(kOrange-2);
	hRapidity->SetFillColor(kOrange-2);
	hPt->SetFillColor(kOrange-2);

	TCanvas*cMass = MakeCanvas("InvMass",LOG_BOTH);
	hMass->Draw("hist"); 
	TCanvas*cRapidity = MakeCanvas("Rapidity",LOG_Y);
	hRapidity->Draw("hist"); 
	TCanvas*cPt = MakeCanvas("Pt",LOG_Y);
	hPt->Draw("hist"); 

	cMass->SaveAs(saveNameMass);
	cRapidity->SaveAs(saveNameRapidity);
	cPt->SaveAs(saveNamePt);
}

TCanvas*MakeCanvas(TString cname,LogType logtype)
{
	TCanvas*canvas = new TCanvas(cname,"",0,0,1000,1000);
	canvas->SetGrid();
	if(logtype==LOG_X)canvas->SetLogx();
	else if(logtype==LOG_Y)canvas->SetLogy();
	else if(logtype==LOG_BOTH){
		canvas->SetLogx();
		canvas->SetLogy();
	}
	return canvas;
}
