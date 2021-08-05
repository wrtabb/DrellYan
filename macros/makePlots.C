using namespace DrellYanVariables;

enum LogType{
	LOG_X,
	LOG_Y,
	LOG_BOTH,
	LINEAR
};

TCanvas*MakeCanvas(TString cname,LogType logtype);

void makePlots(LepType lepType,SampleType sampleType,NtupleType ntupType)
{
	int nProcesses = 4;
	TFile*load_file;
	TH1D*hMass[nProcesses];
	TH1D*hRapidity[nProcesses];
	TH1D*hPt[nProcesses];
	TString saveNameMass;
	TString saveNameRapidity;
	TString saveNamePt;
	
	if(sampleType==SAMPLE_LL){
		saveNameMass = "plots/DYtoLL";
		saveNameRapidity = "plots/DYtoLL";
		saveNamePt = "plots/DYtoLL";
	}
	else if(sampleType==SAMPLE_TAU){
		saveNameMass = "plots/TAU";
		saveNameRapidity = "plots/TAU";
		saveNamePt = "plots/TAU";
	}
	saveNameMass += "_InvMass";
	saveNameRapidity += "_Rapidity";
	saveNamePt += "_Pt";

	vector<TString> histTag = {
		"HardProcess",
		"FSR",
		"Reco",
		"Pruned"
	};
	int l;
	int nHistTag = histTag.size();
	TString file_load= "data/DYHists";
	if(ntupType==V2P3) file_load += "_v2p3";
	else if(ntupType==V2P6) file_load += "_v2p6";
	else if(ntupType==TEST) file_load += "_TEST";
	else if(ntupType==SINGLE_TEST) file_load += "_SINGLE_TEST";
	else{
		cout << "ERROR: NtupleType not properly chosen" << endl;
		return;
	}

	if(sampleType==SAMPLE_LL) file_load += "_DYtoLL";
	else if(sampleType==SAMPLE_TOP) file_load += "_TT";
	else if(sampleType==SAMPLE_FAKE) file_load += "_Fakes";
	else if(sampleType==SAMPLE_DIBOSON) file_load += "_Dibosons";
	else if(sampleType==SAMPLE_TAU) file_load += "_TauTau";
	else if(sampleType==SAMPLE_DATA) file_load += "_Data";
	else if(sampleType==SAMPLE_ALL) file_load += "_All";
	else{
		cout << "ERROR: SampleType not properly chosen" << endl;
		return;
	}

	if(lepType==ELE) file_load += "_EE";
	else if(lepType==MUON) file_load += "_MuMu";
	else {
		cout << "ERROR: LepType not properly chosen" << endl;
		return;
	}
	
	file_load += ".root";

	if(lepType==ELE){
		load_file = new TFile(file_load);
		for(int i=0;i<nProcesses;i++){
			TString massLoad = "histInvMass_EE_DYtoLL";
			massLoad += histTag.at(i);
			TString rapLoad = "histRapidity_EE_DYtoLL";
			rapLoad += histTag.at(i);
			TString ptLoad = "histPt_EE_DYtoLL";
			ptLoad += histTag.at(i);

			hMass[i] = (TH1D*)load_file->Get(massLoad);
			hRapidity[i] = (TH1D*)load_file->Get(rapLoad);
			hPt[i] = (TH1D*)load_file->Get(ptLoad);
		}	
		saveNameMass += "_EE.png";
		saveNameRapidity += "_EE.png";
		saveNamePt += "_EE.png";
	}
	else if(lepType==MUON){
		load_file = new TFile(file_load);
		for(int i=0;i<nProcesses;i++){
			TString massLoad = "histInvMass_MuMu_DYtoLL";
			massLoad += histTag.at(i);
			TString rapLoad = "histRapidity_MuMu_DYtoLL";
			rapLoad += histTag.at(i);
			TString ptLoad = "histPt_MuMu_DYtoLL";
			ptLoad += histTag.at(i);

			hMass[i] = (TH1D*)load_file->Get(massLoad);
			hRapidity[i] = (TH1D*)load_file->Get(rapLoad);
			hPt[i] = (TH1D*)load_file->Get(ptLoad);
		}	

		saveNameMass += "_MuMu.png";
		saveNameRapidity += "_MuMu.png";
		saveNamePt += "_MuMu.png";
	}

	for(int i=0;i<nProcesses;i++){
		hMass[i]->SetFillColor(kOrange-2);
		hMass[i]->SetMinimum(1e-1);
		hRapidity[i]->SetFillColor(kOrange-2);
		hPt[i]->SetFillColor(kOrange-2);
	}

	TCanvas*cMass = MakeCanvas("InvMass",LOG_BOTH);
	TCanvas*cRapidity = MakeCanvas("Rapidity",LOG_Y);
	TCanvas*cPt = MakeCanvas("Pt",LOG_Y);
	for(int i=0;i<nProcesses;i++){
		cMass->cd(i);
		hMass[i]->Draw("hist");

		cRapidity->cd(i);
		hRapidity[i]->Draw("hist");
		
		cPt->cd(i);
		hPt[i]->Draw("hist");
	}

	cMass->SaveAs(saveNameMass);
	delete cMass;
	cRapidity->SaveAs(saveNameRapidity);
	delete cRapidity;
	cPt->SaveAs(saveNamePt);
	delete cPt;
}

TCanvas*MakeCanvas(TString cname,LogType logtype)
{
	TCanvas*canvas = new TCanvas(cname,"",0,0,1000,1000);
	canvas->Divide(2,2);
	if(logtype==LOG_X)canvas->SetLogx();
	else if(logtype==LOG_Y)canvas->SetLogy();
	else if(logtype==LOG_BOTH){
		canvas->SetLogx();
		canvas->SetLogy();
	}
	return canvas;
}
