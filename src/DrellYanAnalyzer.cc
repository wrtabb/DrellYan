#include "DrellYanAnalyzer.hh"
///////////////////////////////////////////////
//-----Initialize DrellYanAnalyzer Class-----//
///////////////////////////////////////////////
DrellYanAnalyzer::DrellYanAnalyzer(DrellYanVariables::NtupleType ntupType,
				   DrellYanVariables::LepType lepType,
				   DrellYanVariables::FileName fileName)
{
	using namespace DrellYanVariables;
	_lepType = lepType;
	_ntupType = ntupType;
	_fileName = fileName;
	if(ntupType==V2P6) _FileToLoad = base_directory_v2p6;
	else if(ntupType==TEST) _FileToLoad = base_directory_test;
	else{
		cout << "NtupleType not properly chosen in DrellYanAnalyzer()" << endl;
		cout << "See include/DrellYanVariables.h for more information" << endl;
	}

	if(_lepType==ELE) _FileToLoad+="skims_EE/";
	else if(_lepType==MUON) _FileToLoad+="skims_MuMu/";

	_FileToLoad += GetSampleName();
	_sampleType = GetSampleType();
	if(_sampleType==SAMPLE_DATA) _isMC = false;
	else _isMC = true;

	// Scale Factor files
	_fRecoSF  = new TFile("data/Reco_SF.root");
	_fMedIDSF = new TFile("data/MediumID_SF.root");
	_fLeg2SF  = new TFile("data/Leg2_SF.root");
	_hRecoSF  = (TH2F*)_fRecoSF->Get("EGamma_SF2D"); 
	_hMedIDSF = (TH2F*)_fMedIDSF->Get("EGamma_SF2D"); 
	_hLeg2SF  = (TH2F*)_fLeg2SF->Get("EGamma_SF2D");

	// Pileup weighting files
	TFile*puFile = new TFile("data/pileup.root");
        _hPileupRatio = (TH1F*)puFile->Get("hPileupRatio");

}//end function DrellYanAnalyzer()

DrellYanVariables::SampleType DrellYanAnalyzer::GetSampleType()
{
	using namespace DrellYanVariables;
	SampleType sampleType;
	if(_fileName==DYLL_M10to50     || _fileName==DYLL_M50to100   || 
	   _fileName==DYLL_M100to200   || _fileName==DYLL_M200to400  || 
	   _fileName==DYLL_M400to500   || _fileName==DYLL_M500to700  ||
	   _fileName==DYLL_M700to800   || _fileName==DYLL_M800to1000 || 
	   _fileName==DYLL_M1000to1500 || _fileName==DYLL_M1500to2000|| 
	   _fileName==DYLL_M2000to3000) 
		sampleType=SAMPLE_LL;

	else if(_fileName==ST_tbarW || _fileName==ST_tW || _fileName==	ttbar_M0to700 || 
		_fileName==ttbar_M700to1000 || _fileName==ttbar_M1000toInf)
		sampleType=SAMPLE_TOP;

	else if(_fileName==WJetsToLNu_amcatnlo_ext2v5 || 
		_fileName==WJetsToLNu_amcatnlo_ext ||
		_fileName==WJetsToLNu_amcatnlo)
		sampleType=SAMPLE_FAKE;

	else if(_fileName==WW || _fileName==WZ || _fileName==ZZ)
		sampleType=SAMPLE_DIBOSON;

	else if(_fileName==DYLL_M10to50_TauTau     || _fileName==DYLL_M50to100_TauTau   ||
		_fileName==DYLL_M100to200_TauTau   || _fileName==DYLL_M200to400_TauTau  ||
		_fileName==DYLL_M400to500_TauTau   || _fileName==DYLL_M500to700_TauTau  ||
		_fileName==DYLL_M700to800_TauTau   || _fileName==DYLL_M800to1000_TauTau ||
		_fileName==DYLL_M1000to1500_TauTau || _fileName==DYLL_M1500to2000_TauTau||
		_fileName==DYLL_M2000to3000_TauTau)
		sampleType=SAMPLE_TAU;

	else if(_fileName==DATA_RunB     || _fileName==DATA_RunC || _fileName==DATA_RunD || 
		_fileName==DATA_RunE     || _fileName==DATA_RunF || _fileName==DATA_RunG ||
		_fileName==DATA_RunHver2 || _fileName==DATA_RunHver3)
		sampleType=SAMPLE_DATA;

	else{
		cout << "ERROR in GetSampleType()" << endl;
		cout << "file name not properly defined" << endl;
		return SAMPLE_ERR;	
	}
	return sampleType;
}

///////////////////////
//-----Load data-----//
///////////////////////
int DrellYanAnalyzer::LoadData()
{
	int returnCode = LoadTrees();
	return returnCode;	
}
Long64_t DrellYanAnalyzer::GetNevents()
{
	return _nEvents;	
}

bool DrellYanAnalyzer::GetIsMC()
{
	return _isMC;
}//end GetIsMC

int DrellYanAnalyzer::LoadTrees()
{
	using namespace DrellYanVariables;
	int returnCode = 1;
	TTimeStamp ts_start;
	cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
	TStopwatch totaltime;
	totaltime.Start();

	TString load_file = _FileToLoad+".root";
	cout << "Begin loading " << load_file << endl;
	_tree = new TChain(_treeName);
	_tree->Add(load_file);
	_nEvents = _tree->GetEntries();;

	totaltime.Stop();
	Double_t TotalCPURunTime = totaltime.CpuTime();
	Double_t TotalRunTime = totaltime.RealTime();
	TTimeStamp ts_end;
	cout << endl;
	cout << "End loading trees:" << endl;
	cout << "Entries loaded: " << _nEvents << endl;
	cout << "**************************************************************************" << endl;
	cout << "Total CPU RunTime: " << TotalCPURunTime << " seconds" << endl;
	cout << "Total Real RunTime: " << TotalRunTime << " seconds" << endl;
	cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;
	cout << "**************************************************************************" << endl;
	cout << endl;

	InitializeBranches(_tree);
	
	return returnCode;
}

int DrellYanAnalyzer::InitializeBranches(TChain*chain)
{
	using namespace DrellYanVariables;
	if(_lepType==ELE){
		//Electrons
		chain->SetBranchAddress("Nelectrons",&Nelectrons,&b_Nelectrons);
		chain->SetBranchAddress("Electron_pT",&Electron_pT,&b_Electron_pT);
		chain->SetBranchAddress("Electron_eta",&Electron_eta,&b_Electron_eta);
		chain->SetBranchAddress("Electron_phi",&Electron_phi,&b_Electron_phi);
		chain->SetBranchAddress("Electron_passMediumID",&Electron_passMediumID,
					 &b_Electron_passMediumID);
	}
	if(_lepType==MUON){
	//Muons
		chain->SetBranchAddress("nMuon",&nMuon,&b_nMuon);
		chain->SetBranchAddress("Muon_pT",&Muon_pT,&b_Muon_pT);
		chain->SetBranchAddress("Muon_eta",&Muon_eta,&b_Muon_eta);
		chain->SetBranchAddress("Muon_phi",&Muon_phi,&b_Muon_phi);
		chain->SetBranchAddress("Muon_passTightID",&Muon_passTightID,
					&b_Muon_passTightID);
		chain->SetBranchAddress("Muon_charge",&Muon_charge,&b_Muon_charge);
		chain->SetBranchAddress("Muon_PfChargedHadronIsoR04",
					&Muon_PfChargedHadronIsoR04,
					&b_Muon_PfChargedHadronIsoR04);
		chain->SetBranchAddress("Muon_PfNeutralHadronIsoR04",
					&Muon_PfNeutralHadronIsoR04,
					&b_Muon_PfNeutralHadronIsoR04);
		chain->SetBranchAddress("Muon_PfGammaIsoR04",
					&Muon_PfGammaIsoR04,
					&b_Muon_PfGammaIsoR04);
		chain->SetBranchAddress("Muon_PFSumPUIsoR04",
					&Muon_PFSumPUIsoR04,
					&b_Muon_PFSumPUIsoR04);
	}
	//Trigger
	chain->SetBranchAddress("HLT_ntrig",&HLT_ntrig,&b_HLT_ntrig);
	chain->SetBranchAddress("HLT_trigType",&HLT_trigType,&b_HLT_trigType);
	chain->SetBranchAddress("HLT_trigFired",&HLT_trigFired,&b_HLT_trigFired);
	chain->SetBranchAddress("HLT_trigName",&pHLT_trigName);

	chain->SetBranchAddress("PVz",&PVz,&b_PVz);
	chain->SetBranchAddress("nVertices",&nVertices,&b_nVertices);
	chain->SetBranchAddress("nPileUp",&nPileUp,&b_nPileUp);

	chain->SetBranchAddress("vtxTrkCkt1Pt",&pvtxTrkCkt1Pt);
	chain->SetBranchAddress("vtxTrkCkt2Pt",&pvtxTrkCkt2Pt);
	chain->SetBranchAddress("vtxTrkChi2",&pvtxTrkChi2);
	chain->SetBranchAddress("vtxTrkNdof",&pvtxTrkNdof);

	if(_isMC){
		chain->SetBranchAddress("GENnPair", &GENnPair, &b_GENnPair);
		chain->SetBranchAddress("GENLepton_eta", &GENLepton_eta, &b_GENLepton_eta);
		chain->SetBranchAddress("GENLepton_phi",&GENLepton_phi, &b_GENLepton_phi);
		chain->SetBranchAddress("GENLepton_pT",&GENLepton_pT, &b_GENLepton_pT);
		chain->SetBranchAddress("GENLepton_ID",&GENLepton_ID, &b_GENLepton_ID);
		chain->SetBranchAddress("GENLepton_isHardProcess",&GENLepton_isHardProcess,
					&b_GENLepton_isHardProcess);
		chain->SetBranchAddress("GENLepton_fromHardProcessFinalState",
					&GENLepton_fromHardProcessFinalState,
					&b_GENLepton_fromHardProcessFinalState);
		chain->SetBranchAddress("nGenOthers",&nGenOthers,&b_nGenOthers);
		chain->SetBranchAddress("GenOthers_eta",&GenOthers_eta,&b_GenOthers_eta);
		chain->SetBranchAddress("GenOthers_phi",&GenOthers_phi,&b_GenOthers_phi);
		chain->SetBranchAddress("GenOthers_pT",&GenOthers_pT,&b_GenOthers_pT);
		chain->SetBranchAddress("GenOthers_ID",&GenOthers_ID,&b_GenOthers_ID);
		chain->SetBranchAddress("GenOthers_isHardProcess",&GenOthers_isHardProcess,
					&b_GenOthers_isHardProcess);
		chain->SetBranchAddress("GenOthers_isPromptFinalState",
					&GenOthers_isPromptFinalState,
					&b_GenOthers_isPromptFinalState);
		chain->SetBranchAddress("GENEvt_weight",&GENEvt_weight,&b_GENEvt_weight);
	}
	return 1;
}//end initializing trees

void DrellYanAnalyzer::InitializeHistograms()
{
	using namespace DrellYanVariables;
	TString lepTag;
	
	TString histNameMass = "hMass";
	histNameMass += "HardProcess";
	TString histNameRapidity = "hRapidity";
	histNameRapidity += "HardProcess";
	TString histNamePt = "hPt";
	histNamePt += "HardProcess";

	_hMassHardProcess = new TH1D(histNameMass,"",nMassBins,massbins);
	DefineHistogramProperties(_hMassHardProcess);
	_hRapidityHardProcess = new TH1D(histNameRapidity,"",nRapidityBins,rapiditybins);
	DefineHistogramProperties(_hRapidityHardProcess);
	_hPtHardProcess = new TH1D(histNamePt,"",nPtBins,ptbins);
	DefineHistogramProperties(_hPtHardProcess);

	histNameMass = "hMass";
	histNameMass += "Reco";
	histNameRapidity = "hRapidity";
	histNameRapidity += "Reco";
	histNamePt = "hPt";
	histNamePt += "Reco";

	_hMassReco = new TH1D(histNameMass,"",nMassBins,massbins);
	DefineHistogramProperties(_hMassReco);
	_hRapidityReco = new TH1D(histNameRapidity,"",nRapidityBins,rapiditybins);
	DefineHistogramProperties(_hRapidityReco);
	_hPtReco = new TH1D(histNamePt,"",nPtBins,ptbins);
	DefineHistogramProperties(_hPtReco);
}//end InitializeHistograms

void DrellYanAnalyzer::DefineHistogramProperties(TH1*hist)
{
	using namespace DrellYanVariables;
	if(_sampleType==SAMPLE_LL){
		hist->SetFillColor(kOrange-2);
		hist->SetLineColor(kOrange+3);
	}
	else if(_sampleType==SAMPLE_TOP){
		hist->SetFillColor(kBlue+2);
		hist->SetLineColor(kBlue+3);
	}
	else if(_sampleType==SAMPLE_FAKE){
		hist->SetFillColor(kViolet+5);
		hist->SetLineColor(kViolet+3);
	}
	else if(_sampleType==SAMPLE_DIBOSON){
		hist->SetFillColor(kRed+2);
		hist->SetLineColor(kRed+4);
	}
	else if(_sampleType==SAMPLE_TAU){
		hist->SetFillColor(kGreen+2);
		hist->SetLineColor(kGreen+3);
	}
	else if(_sampleType==SAMPLE_DATA){
		hist->SetMarkerStyle(20);
		hist->SetMarkerColor(kBlack);
		hist->SetLineColor(kBlack);
		hist->SetFillColor(kWhite);
	}
	else cout << "ERROR in DefineHistogramProperties" << endl;

	return;
}//end DefineHistogramProperties()

/////////////////////////
//-----Get Weights-----//
/////////////////////////
double DrellYanAnalyzer::GetSampleWeights()
{
	// This function returns the weights which are determeined per sample
	using namespace DrellYanVariables;
	//-----Cross section-----//
	double xSecWeight = 1.0;
	double xSec = 1.0;

	// Drell-Yan to leptons
	if(_fileName==DYLL_M10to50 ||
	   _fileName==DYLL_M10to50_TauTau) xSec = xSec_LL.at(0);
	else if(_fileName==DYLL_M50to100 ||
		_fileName==DYLL_M50to100_TauTau) xSec = xSec_LL.at(1);
	else if(_fileName==DYLL_M100to200 ||
		_fileName==DYLL_M100to200_TauTau) xSec = xSec_LL.at(2);
	else if(_fileName==DYLL_M200to400 ||
		_fileName==DYLL_M200to400_TauTau) xSec = xSec_LL.at(3);
	else if(_fileName==DYLL_M400to500 ||
		_fileName==DYLL_M400to500_TauTau) xSec = xSec_LL.at(4);
	else if(_fileName==DYLL_M500to700 ||
		_fileName==DYLL_M500to700_TauTau) xSec = xSec_LL.at(5);
	else if(_fileName==DYLL_M700to800 ||
		_fileName==DYLL_M700to800_TauTau) xSec = xSec_LL.at(6);
	else if(_fileName==DYLL_M800to1000 ||
		_fileName==DYLL_M800to1000_TauTau) xSec = xSec_LL.at(7);
	else if(_fileName==DYLL_M1000to1500 ||
		_fileName==DYLL_M1000to1500_TauTau) xSec = xSec_LL.at(8);
	else if(_fileName==DYLL_M1500to2000 ||
		_fileName==DYLL_M1500to2000_TauTau) xSec = xSec_LL.at(9);
	else if(_fileName==DYLL_M2000to3000 ||
		_fileName==DYLL_M2000to3000_TauTau) xSec = xSec_LL.at(10);

	// tops
	else if(_fileName==ST_tbarW) 	     xSec = xSec_tops.at(0);
	else if(_fileName==ST_tW) 	     xSec = xSec_tops.at(1);
	else if(_fileName==ttbar_M0to700)    xSec = xSec_tops.at(2);
	else if(_fileName==ttbar_M700to1000) xSec = xSec_tops.at(3);
	else if(_fileName==ttbar_M1000toInf) xSec = xSec_tops.at(4);
	
	// Fakes
	else if(_fileName==WJetsToLNu_amcatnlo)     xSec = xSec_fakes.at(0); 
	else if(_fileName==WJetsToLNu_amcatnlo_ext) xSec = xSec_fakes.at(1);

	// Dibosons
	else if(_fileName==WW) xSec = xSec_dibosons.at(0);
	else if(_fileName==WZ) xSec = xSec_dibosons.at(1);
	else if(_fileName==ZZ) xSec = xSec_dibosons.at(2);
	
	//When cross section weights are applied with gen weights,
	//You don't divide by the number of events
	//For the TEST cases, I don't use gen weights because
	//Gen weights are only applied on the full sample
	if(_ntupType==V2P6) xSecWeight = dataLuminosity*xSec;
	else if(_ntupType==TEST) xSecWeight = dataLuminosity*xSec/_nEvents;

	//-----Total Weight-----//
	double weight = xSecWeight;
	return weight;
}//end GetSampleWeights()

double DrellYanAnalyzer::GetEventWeights(double pt1,double pt2,double eta1,double eta2)
{
	using namespace DrellYanVariables;
	
	// Pileup weights
	double puWeight = GetPUWeight();

	// Electron Scale Factors
	double eleSF = 1.0;
	if(_isMC && _lepType==ELE){
		eleSF = GetEleScaleFactors(pt1,pt2,eta1,eta2);
	}
	return eleSF*puWeight;	
}//end GetEventWeights()

double DrellYanAnalyzer::GetGenWeightSum()
{
	using namespace DrellYanVariables;

	TTimeStamp ts_start;
	cout << endl;
	cout << "****************************************" << endl;
	cout << "Getting gen weights :" << endl;
	cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
	TStopwatch totaltime;
	totaltime.Start();

	double sumGenWeight = 0.0;
	double sumRawGenWeight = 0.0;
	double varGenWeight = 0.0;
	double genWeight;
	Long64_t localEntry;
	for(Long64_t i=0;i<_nEvents;i++){
		localEntry = _tree->LoadTree(i);
		b_GENEvt_weight->GetEntry(localEntry);
		genWeight = GENEvt_weight/fabs(GENEvt_weight);//normalized genweight
		sumGenWeight += genWeight;
		varGenWeight += GENEvt_weight*GENEvt_weight;//variance of genweights
		sumRawGenWeight += GENEvt_weight;
		Counter(i,_nEvents);
	}
	totaltime.Stop();
	Double_t TotalCPURunTime = totaltime.CpuTime();
	Double_t TotalRunTime = totaltime.RealTime();
	TTimeStamp ts_end;
	cout << endl;
	cout << "End getting gen weights:" << endl;
	cout << "**************************************************************************" << endl;
	cout << "Total CPU RunTime: " << TotalCPURunTime << " seconds" << endl;
	cout << "Total Real RunTime: " << TotalRunTime << " seconds" << endl;
	cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;
	cout << "**************************************************************************" << endl;
	cout << endl;
	return sumGenWeight;
}//end GetGenWeight

double DrellYanAnalyzer::GetPUWeight()
{
	using namespace DrellYanVariables;
	double puWeight = 1.0;	

	puWeight = _hPileupRatio->GetBinContent(_hPileupRatio->FindBin(nPileUp));
	return puWeight;
}

double DrellYanAnalyzer::GetEleScaleFactors(double pt1,double pt2,double eta1,double eta2)
{

	double sfReco1,sfReco2,sfID1,sfID2,sfHLT,sfWeight;

	double ptBinHigh = 500;
	double etaBinLow = -2.5;
	double etaBinHigh = 2.5;
 
	if(pt1>ptBinHigh) pt1 = ptBinHigh;
	if(pt2>ptBinHigh) pt2 = ptBinHigh;
	if(eta1>etaBinHigh) eta1 = etaBinHigh;
	if(eta2>etaBinHigh) eta2 = etaBinHigh;
	if(eta1<etaBinLow)  eta1 = etaBinLow;
	if(eta2<etaBinLow)  eta2 = etaBinLow;

	sfReco1 =  _hRecoSF->GetBinContent(_hRecoSF->FindBin(eta1,pt1));
	sfReco2 =  _hRecoSF->GetBinContent(_hRecoSF->FindBin(eta2,pt2));
	sfID1   = _hMedIDSF->GetBinContent(_hMedIDSF->FindBin(eta1,pt1));
	sfID2   = _hMedIDSF->GetBinContent(_hMedIDSF->FindBin(eta2,pt2));
	sfHLT   = (_hLeg2SF->GetBinContent(_hLeg2SF->FindBin(eta1,pt1)))*
		  (_hLeg2SF->GetBinContent(_hLeg2SF->FindBin(eta2,pt2)));

	sfWeight = sfReco1*sfReco2*sfID1*sfID2*sfHLT;

	return sfWeight;
}//end GetScaleFactors

///////////////////////////////////////
//-----Get Leptons From an Event-----//
///////////////////////////////////////
int DrellYanAnalyzer::EventLoop()
{
	using namespace DrellYanVariables;
	InitializeHistograms();

	//Initialize lepton mass variable
	double lMass;
	if(_lepType==ELE) lMass = eMass;
	else if(_lepType==MUON) lMass = muMass;
	else{
		cout << "Error in EventLoop(): Need to specify ELE or MUON" << endl;
		 return 0;
	}

	double xSecWeight = 1.0;
	double sumGenWeight = 1.0;
	if(_isMC){
		xSecWeight = GetSampleWeights();
		sumGenWeight = GetGenWeightSum();
	}

	//Loop over all events in the loaded tree
	TTimeStamp ts_start;
	cout << endl;
	cout << "********************************" << endl;
	cout << "Begin EventLoop:" << endl;
	cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
	TStopwatch totaltime;
	totaltime.Start();
	for(Long64_t iEntry=0;iEntry<_nEvents;iEntry++){

		//Initialize variables such that if they are not assigned a value 
		//in the loop, they will be in the underflow bin of their histogram

		//hard process variables
		double hardPt1 		= -1000;
		double hardPt2 		= -1000;
		double hardEta1 	= -1000;
		double hardEta2 	= -1000;
		double hardPhi1 	= -1000;
		double hardPhi2 	= -1000; 
		double invMassHard 	= -1000;
		double rapidityHard 	= -1000;
		double ptHard 		= -1000;
		
		//FSR variables
		double fsrPt1 		= -1000;
		double fsrPt2 		= -1000;
		double fsrEta1 		= -1000;
		double fsrEta2 		= -1000;
		double fsrPhi1 		= -1000;
		double fsrPhi2 		= -1000; 
		double invMassFSR 	= -1000;
		double rapidityFSR 	= -1000;
		double ptFSR 		= -1000;

		//Reco variables
		double recoPt1 		= -1000;
		double recoPt2 		= -1000;
		double recoEta1 	= -1000;
		double recoEta2 	= -1000;
		double recoPhi1 	= -1000;
		double recoPhi2 	= -1000; 
		double invMassReco 	= -1000;
		double rapidityReco 	= -1000;
		double ptReco 		= -1000;

		//Pruned variables
		double prunedPt1 	= -1000;
		double prunedPt2 	= -1000;
		double prunedEta1 	= -1000;
		double prunedEta2 	= -1000;
		double prunedPhi1 	= -1000;
		double prunedPhi2 	= -1000; 
		double invMassPruned 	= -1000;
		double rapidityPruned 	= -1000;
		double ptPruned 	= -1000;

		//Initialize lepton indices
		int iHard1    = -1;
		int iHard2    = -1;
		int iFSR1     = -1;
		int iFSR2     = -1;
		int iDressed1 = -1;
		int iDressed2 = -1;
		int iLep1     = -1;
		int iLep2     = -1;

		//Get event from tree
		_tree->GetEntry(iEntry);

		//Check if event passes HLT criteria
//		bool passHLT = PassHLT(); 
		bool passHLT = true;
		//Get gen level leptons
		int nDileptonsGen = GetGenLeptons(iHard1,iHard2,
						  iFSR1,iFSR2);

		//Get reco level leptons
		int nDileptonsReco = GetRecoLeptons(iLep1,iLep2);

		// CHeck if reco leptons match with gen leptons
		bool passGenToRecoMatch1 = PassGenToRecoMatch(iFSR1,iLep1);
		bool passGenToRecoMatch2 = PassGenToRecoMatch(iFSR2,iLep2);
		bool passGenToReco = true;
		if(_isMC){
			if(!(passGenToRecoMatch1 && passGenToRecoMatch2))  
				passGenToReco = false;
		}

		//reco quantities
		if(nDileptonsReco==1 && passHLT){
			if(_lepType==ELE){
				recoPt1  = Electron_pT[iLep1];
				recoEta1 = Electron_eta[iLep1];
				recoPhi1 = Electron_phi[iLep1];
				recoPt2  = Electron_pT[iLep2];
				recoEta2 = Electron_eta[iLep2];
				recoPhi2 = Electron_phi[iLep2];
			}
			else if(_lepType==MUON){
				recoPt1  = Muon_pT[iLep1];
				recoEta1 = Muon_eta[iLep1];
				recoPhi1 = Muon_phi[iLep1];
				recoPt2  = Muon_pT[iLep2];
				recoEta2 = Muon_eta[iLep2];
				recoPhi2 = Muon_phi[iLep2];
			}
			else{
				cout << "ERROR in Event Loop!" << endl;
				cout << "Lepton type not specified" << endl;
				return 0;
			}

			invMassReco   = GetInvMass(recoPt1,recoEta1,recoPhi1,
						   recoPt2,recoEta2,recoPhi2);
			rapidityReco = GetRapidity(recoPt1,recoEta1,recoPhi1,
						   recoPt2,recoEta2,recoPhi2);
			ptReco 		   = GetPt(recoPt1,recoEta1,recoPhi1,
						   recoPt2,recoEta2,recoPhi2);
		}// end if passHLT and passGenToReco

		// Gen level quantities
		if(nDileptonsGen==1 && passHLT){
			hardPt1  = GENLepton_pT[iHard1];
			hardEta1 = GENLepton_eta[iHard1];
			hardPhi1 = GENLepton_phi[iHard1];
			hardPt2  = GENLepton_pT[iHard2];
			hardEta2 = GENLepton_eta[iHard2];
			hardPhi2 = GENLepton_phi[iHard2];

			fsrPt1  = GENLepton_pT[iFSR1];
			fsrEta1 = GENLepton_eta[iFSR1];
			fsrPhi1 = GENLepton_phi[iFSR1];
			fsrPt2  = GENLepton_pT[iFSR2];
			fsrEta2 = GENLepton_eta[iFSR2];
			fsrPhi2 = GENLepton_phi[iFSR2];

			//This is temporary
			//Still need to write a function
			//To calculate pruned leptons
			prunedPt1  = GENLepton_pT[iHard1];
			prunedEta1 = GENLepton_eta[iHard1];
			prunedPhi1 = GENLepton_phi[iHard1];
			prunedPt2  = GENLepton_pT[iHard2];
			prunedEta2 = GENLepton_eta[iHard2];
			prunedPhi2 = GENLepton_phi[iHard2];

			//hard process quantities
			invMassHard   = GetInvMass(hardPt1,hardEta1,hardPhi1,
				 		   hardPt2,hardEta2,hardPhi2);
			rapidityHard = GetRapidity(hardPt1,hardEta1,hardPhi1,
						   hardPt2,hardEta2,hardPhi2);
			ptHard 	           = GetPt(hardPt1,hardEta1,hardPhi1,
				       		   hardPt2,hardEta2,hardPhi2);

			//FSR quantities
			invMassFSR   = GetInvMass(fsrPt1,fsrEta1,fsrPhi1,
						  fsrPt2,fsrEta2,fsrPhi2);
			rapidityFSR = GetRapidity(fsrPt1,fsrEta1,fsrPhi1,
						  fsrPt2,fsrEta2,fsrPhi2);
			ptFSR 		  = GetPt(fsrPt1,fsrEta1,fsrPhi1,
				      		  fsrPt2,fsrEta2,fsrPhi2);


			//Pruned quantities
			invMassPruned   = GetInvMass(prunedPt1,prunedEta1,
						     prunedPhi1,prunedPt2,
						     prunedEta2,prunedPhi2);
			rapidityPruned = GetRapidity(prunedPt1,prunedEta1,
						     prunedPhi1,prunedPt2,
						     prunedEta2,prunedPhi2);
			ptPruned 	     = GetPt(prunedPt1,prunedEta1,
					 	     prunedPhi1,prunedPt2,
					 	     prunedEta2,prunedPhi2);
		}//end 

		//Calculate genweight for MC samples
		double genWeight = 1.0;
		double eventWeight = 1.0;
		if(_isMC){
			genWeight = (GENEvt_weight/fabs(GENEvt_weight))/sumGenWeight;
			eventWeight = GetEventWeights(recoPt1,recoPt2,recoEta1,recoEta2);
		}
	
		double weights = xSecWeight*eventWeight*genWeight;

		//Fill all histograms
		_hMassHardProcess->Fill(invMassHard,weights);
		_hRapidityHardProcess->Fill(rapidityHard,weights);
		_hPtHardProcess->Fill(ptHard,weights);
		_hMassReco->Fill(invMassReco,weights);
		_hRapidityReco->Fill(rapidityReco,weights);
		_hPtReco->Fill(ptReco,weights);
		Counter(iEntry,_nEvents);
	}//end event loop

	//Save histograms to a root file
	SaveResults();

	totaltime.Stop();
	Double_t TotalCPURunTime = totaltime.CpuTime();
	Double_t TotalRunTime = totaltime.RealTime();
	TTimeStamp ts_end;
	cout << endl;
	cout << "Ending EventLoop:" << endl;
	cout << "**************************************************************************" << endl;
	cout << "Total CPU RunTime: " << TotalCPURunTime << " seconds" << endl;
	cout << "Total Real RunTime: " << TotalRunTime << " seconds" << endl;
	cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;
	cout << "**************************************************************************" << endl;
	cout << endl;

	return 1;
}//end EventLoop

int DrellYanAnalyzer::GetGenLeptons(int &iHard1,int &iHard2,
				    int &iFSR1,int &iFSR2)
{
	using namespace DrellYanVariables;
	int lepID = 0;
	if      (_lepType==ELE)  lepID = 11;
	else if (_lepType==MUON) lepID = 13;
	else {
		cout << "ERROR: Appropriate lepton not selected." << endl;
		cout << "NOTE: This analysis does not handle taus." << endl;
		return 0;
	}
	int nDileptons = 0;

	for(int iLep=0;iLep<GENnPair;iLep++){
		for(int jLep=iLep+1;jLep<GENnPair;jLep++){
			//Lepton selection
			if(!(abs(GENLepton_ID[iLep])==lepID && 
			     abs(GENLepton_ID[jLep])==lepID)) continue;
			//require opposite sign
			if(GENLepton_ID[iLep]*GENLepton_ID[jLep]>0) continue;
			//Get leptons from hard process
			if(GENLepton_isHardProcess[iLep]==1 && 
			   GENLepton_isHardProcess[jLep]==1){
				iHard1 = iLep;
				iHard2 = jLep;
				nDileptons++;
			}//end if hard process
			//Get leptons from FSR
			if(GENLepton_fromHardProcessFinalState[iLep]==1 &&
			GENLepton_fromHardProcessFinalState[jLep]==1){
				iFSR1 = iLep;
				iFSR2 = jLep;
			}//end if FSR
		}//end jLep loop
	}//end iLep loop
	return nDileptons;
}//end GetGenLeptons()

int DrellYanAnalyzer::GetRecoLeptons(int &leadLep,int &subLep)
{
	//I'm uncertain about how &leadLep will be connected with &leadEle or &leadMu
	//Need to remember to check this behavior when I can
	using namespace DrellYanVariables;

	double ptLead = -1000;
        double ptSub = -1000;
        double etaLead = -1000;
        double etaSub = -1000;
	double phiLead = -1000;
	double phiSub = -1000;
	double chargeLead = 0;
	double chargeSub = 0;
        int idxLead = -1;
        int idxSub = -1;

	if(_lepType==ELE){
		if(Nelectrons<2) return 0;
		for(int iEle = 0; iEle < Nelectrons; iEle++) {
			if(!Electron_passMediumID[iEle]) continue;
			if(Electron_pT[iEle] > ptLead){
				ptLead = Electron_pT[iEle];
				etaLead = Electron_eta[iEle];
				idxLead = iEle;
			}
		}//end lead pt loop

		for(int iEle = 0; iEle < Nelectrons; iEle++) {
			if(!Electron_passMediumID[iEle]) continue;
			if(Electron_pT[iEle] > ptSub && Electron_pT[iEle] < ptLead){
				ptSub = Electron_pT[iEle];
				etaSub = Electron_eta[iEle];
				idxSub = iEle;
			}
		}//end subleading pt loop
	}// end if ELE

	else if(_lepType==MUON){
		if(nMuon<2) return 0;
		// Need to select two muons by which two have the smallest vtx chi2
		// Temporarily I am choosing the two with highest pT
		for(int iMu=0;iMu<nMuon;iMu++){
			if(!(Muon_passTightID[iMu] && PassMuonIsolation(iMu))) continue;
			if(Muon_pT[iMu] > ptLead){
				ptLead = Muon_pT[iMu];
				etaLead = Muon_eta[iMu];
				chargeLead = Muon_charge[iMu];
				idxLead = iMu;
			}// end if muon pass tight ID and pass isolation
		}//end lead loop
		for(int iMu=0;iMu<nMuon;iMu++){
			if(!(Muon_passTightID[iMu] && PassMuonIsolation(iMu))) continue;
			if(Muon_pT[iMu] > ptSub && Muon_pT[iMu] < ptLead){
				ptSub = Muon_pT[iMu];
				etaSub = Muon_eta[iMu];
				chargeSub = Muon_charge[iMu];
				idxSub = iMu;
			}// end if muon pass tight ID and pass isolation
		}//end sub loop 
		//Other cuts not yet applied:
		//smallest dimuon vertex chi2
		//PassVertexChi2(iMu,jMu);

		if(PassMuonAngle(ptLead,etaLead,phiLead,muMass,ptSub,etaSub,phiSub,muMass)&&
		   chargeLead*chargeSub<0 && PassAcceptance(ptLead,ptSub,etaLead,etaSub)){
			leadLep = idxLead;
			subLep = idxSub;
			return 1;
		}//end pass acceptance
		else return 0;
	}// end if MUON

	// Make sure dilepton pair passes kinematic cuts
	// Muons and electrons use the same Kinematic cuts
	if(!PassAcceptance(ptLead,ptSub,etaLead,etaSub)){
		return 0;
	}
	else{
		leadLep = idxLead;
		subLep = idxSub;
		return 1;
	}
}

int DrellYanAnalyzer::GetRecoMuons(int &leadMu,int &subMu)
{
	using namespace DrellYanVariables;
	double ptLead = -1000;
	double ptSub = -1000;
	double etaLead = -1000;
	double etaSub = -1000;
	double phiLead = -1000;
	double phiSub = -1000;
	double chargeLead = 0;
	double chargeSub = 0;
	int idxLead = -1;
	int idxSub = -1;

	//For muons we want exactly two passing all selection criteria
	//The two must be chosen by which has the smallest vertex chi2
	//This is not yet implemented
	
	if(nMuon<2) return 0;
	for(int iMu=0;iMu<nMuon;iMu++){
		if(Muon_passTightID[iMu] && PassMuonIsolation(iMu)){ 
			if(Muon_pT[iMu] > ptLead){
				ptLead = Muon_pT[iMu];
				etaLead = Muon_eta[iMu];
				chargeLead = Muon_charge[iMu];
				idxLead = iMu;
			}
			if(Muon_pT[iMu] > ptSub && Muon_pT[iMu] < ptLead){
				ptSub = Muon_pT[iMu];
				etaSub = Muon_eta[iMu];
				chargeSub = Muon_charge[iMu];
				idxSub = iMu;
			}
		}// end if muon pass tight ID and pass isolation
	}//end iMu loop
	//Other cuts not yet applied:
	//smallest dimuon vertex chi2
	//PassVertexChi2(iMu,jMu);

	if(PassMuonAngle(ptLead,etaLead,phiLead,muMass,ptSub,etaSub,phiSub,muMass) && 
	   chargeLead*chargeSub<0 && PassAcceptance(ptLead,ptSub,etaLead,etaSub)){
		leadMu = idxLead;
		subMu = idxSub;
		return 1;
	}//end pass acceptance
	else return 0;
}//end GetRecoMuons()

//////////////////////////
//-----Cut criteria-----//
//////////////////////////
bool DrellYanAnalyzer::PassAcceptance(double pt1,double pt2,double eta1,double eta2)
{
	using namespace DrellYanVariables;
	//Acceptance criteria is the same for muons and electrons
	
	//I think we are not excluding the gap region in the ECAL, so I commented these out
	//But the lines can be uncommented if we want to exclude the gap
	if(abs(eta1)>etaGapLow && abs(eta1)<etaGapHigh) return false;
	if(abs(eta2)>etaGapLow && abs(eta2)<etaGapHigh) return false;
	
	//Ensure both leptons are below |eta| = 2.4
	if(abs(eta1)>etaHigh||abs(eta2)>etaHigh) return false;

	//Make sure that the leading lepton has pT greater than 28 GeV and 
	//the subleading has pT greater than 17 GeV
	if(!( (pt1>ptLow && pt2>ptHigh) || (pt1>ptHigh && pt2>ptLow) )) return false;
	return true;
}//end PassAcceptance()

bool DrellYanAnalyzer::PassHLT()
{
	using namespace DrellYanVariables;
	TString trigName;
	int trigNameSize = pHLT_trigName->size();
	bool passHLT = false;
	for(int iHLT=0;iHLT<trigNameSize;iHLT++) {
		trigName = pHLT_trigName->at(iHLT);
		if(_lepType==ELE){
			if(trigName.CompareTo(electronTrigger)==0 && HLT_trigFired[iHLT]==1)
			{
				passHLT = true;
				break;
			}//end if trigName...
		}//end if lepType
		else if(_lepType==MUON){
			if((trigName.CompareTo(muonTrigger1)==0 || 
			    trigName.CompareTo(muonTrigger2)==0 ) && 
			    HLT_trigFired[iHLT]==1){
				passHLT = true;
				break;
			}//end if trigName...
		}//end if lepType
	}//end loop over triggers
	return passHLT;
}//end PassHLT

bool DrellYanAnalyzer::PassGenToRecoMatch(int genIndex,int &recoIndex)
{
	using namespace DrellYanVariables;
	if(_lepType==ELE) return PassGenToRecoMatchEle(genIndex,recoIndex);
	//if(_lepType==MUON) return PassGenToRecoMatchMu(genIndex,recoIndex);
	if(_lepType==MUON) return true;
	else{
		cout << "Error in DrellYanAnalyzer::PassGenToRecoMatch()!" << endl;
		cout << "lepType not properly defined" << endl;
		cout << "Must be ELE or MUON" << endl;
		return false;
	}
}//end PassgenToRecoMatch

double DrellYanAnalyzer::GetVertexChi2(int index1,int index2)
{
	double vtxChi2 = 1000000;
	return vtxChi2;	
}//end PassVertexChi2

bool DrellYanAnalyzer::PassGenToRecoMatchEle(int genIndex,int &recoIndex)
{
	using namespace DrellYanVariables;
	double dR,deta,dphi;
	float dRMin = 100000;
	recoIndex=-1;
	for(int iLep=0;iLep<Nelectrons;iLep++){
		deta=Electron_eta[iLep]-GENLepton_eta[genIndex];
		dphi=abs(Electron_phi[iLep]-GENLepton_phi[genIndex]);
		if(dphi>pi) dphi=2*pi-dphi;
		dR=sqrt(deta*deta+dphi*dphi);

		if(dR<dRMin){
			recoIndex=iLep;
			dRMin=dR;
		}
	}
	bool matchFound = true;
	if(dRMin>=dRMinCut){
		recoIndex=-1;
		matchFound=false;
	}
	return matchFound;
}//end PassGenToRecoMatchEle

bool DrellYanAnalyzer::PassGenToRecoMatchMu(int genIndex,int &recoIndex)
{
	using namespace DrellYanVariables;
	double dR,deta,dphi;
	float dRMin = 100000;
	recoIndex=-1;
	for(int iLep=0;iLep<nMuon;iLep++){
		deta=Muon_eta[iLep]-GENLepton_eta[genIndex];
		dphi=abs(Muon_phi[iLep]-GENLepton_phi[genIndex]);
		if(dphi>pi) dphi=2*pi-dphi;
		dR=sqrt(deta*deta+dphi*dphi);

		if(dR<dRMin){
			recoIndex=iLep;
			dRMin=dR;
		}
	}
	bool matchFound = true;
	if(dRMin>=dRMinCut){
		recoIndex=-1;
		matchFound=false;
	}
	return matchFound;
}//end PassGenToRecomatchMu

bool DrellYanAnalyzer::PassMuonAngle(double pt1,double eta1,double phi1,double mass1,
				     double pt2,double eta2,double phi2,double mass2)
{
	using namespace DrellYanVariables;
	double angle;
	angle = CalcVariable(pt1,eta1,phi1,mass1,pt2,eta2,phi2,mass2,MUON_ANGLE);
	if(angle<pi-0.005) return true;
	else return false;
}//end PassMuonAngle

bool DrellYanAnalyzer::PassMuonIsolation(int index)
{
	using namespace DrellYanVariables;
	double chargedIso = Muon_PfChargedHadronIsoR04[index];
	double neutralIso = Muon_PfNeutralHadronIsoR04[index];
	double gammaIso	  = Muon_PfGammaIsoR04[index];
	double sumPUPt    = Muon_PFSumPUIsoR04[index];
	double pT         = Muon_pT[index];
	double iso_dBeta;
       	iso_dBeta = (chargedIso+max(0.0,neutralIso+gammaIso-0.5*sumPUPt))/pT;
	if(iso_dBeta < 0.15) return true;
	else return false;
}

///////////////////////
//-----Utilities-----//
///////////////////////
void DrellYanAnalyzer::SaveResults()
{
	using namespace DrellYanVariables;
	TString filesave = "output_data/DYHists";
	if(_ntupType==V2P6) filesave+= "_v2p6";
	else if (_ntupType==V2P3) filesave+= "_v2p3";
	else if(_ntupType==TEST) filesave += "_TEST";
	else if(_ntupType==SINGLE_FILE) filesave += "_SingleSampleTest";
	else{
		cout << "********************************************************" << endl;
		cout << "ERROR in SaveResults()" << endl;
		cout << "NtupleType not correctly defined" << endl;
		cout << "See include/DrellYanVariables.h for list of NtupleTypes" << endl;
		cout << "********************************************************" << endl;
		return;
	}

	if(_sampleType==SAMPLE_LL) 	      filesave += "_DYtoLL";
	else if(_sampleType==SAMPLE_TOP)      filesave += "_Top";
	else if(_sampleType==SAMPLE_FAKE)     filesave += "_Fake";
	else if (_sampleType==SAMPLE_DIBOSON) filesave += "_Diboson";
	else if (_sampleType==SAMPLE_TAU)     filesave += "_TauTau";
	else if (_sampleType==SAMPLE_DATA)    filesave += "_Data";
	else {
		cout << "*****************************************************" << endl;
		cout << "ERROR in SaveResults() " << endl;	
		cout << "SampleType not correctly defined" << endl;
		cout << "See include/DrellYanVariables.h for list of SampleTypes" << endl;
		cout << "*****************************************************" << endl;
	}

	//MC
	if(_fileName==DYLL_M10to50)          filesave += "_M10to50";
	else if(_fileName==DYLL_M50to100)    filesave += "_M50to100";
	else if(_fileName==DYLL_M100to200)   filesave += "_M100to200";
	else if(_fileName==DYLL_M200to400)   filesave += "_M200to400";
	else if(_fileName==DYLL_M400to500)   filesave += "_M400to500";
	else if(_fileName==DYLL_M500to700)   filesave += "_M500to700";
	else if(_fileName==DYLL_M700to800)   filesave += "_M700to800";
	else if(_fileName==DYLL_M800to1000)  filesave += "_M800to1000";
	else if(_fileName==DYLL_M1000to1500) filesave += "_M1000to1500";
	else if(_fileName==DYLL_M1500to2000) filesave += "_M1500to2000";
	else if(_fileName==DYLL_M2000to3000) filesave += "_M2000to3000";

	// Taus
	else if(_fileName==DYLL_M10to50)     filesave += "_M10to50";
	else if(_fileName==DYLL_M50to100)    filesave += "_M50to100";
	else if(_fileName==DYLL_M100to200)   filesave += "_M100to200";
	else if(_fileName==DYLL_M200to400)   filesave += "_M200to400";
	else if(_fileName==DYLL_M400to500)   filesave += "_M400to500";
	else if(_fileName==DYLL_M500to700)   filesave += "_M500to700";
	else if(_fileName==DYLL_M700to800)   filesave += "_M700to800";
	else if(_fileName==DYLL_M800to1000)  filesave += "_M800to1000";
	else if(_fileName==DYLL_M1000to1500) filesave += "_M1000to1500";
	else if(_fileName==DYLL_M1500to2000) filesave += "_M1500to2000";
	else if(_fileName==DYLL_M2000to3000) filesave += "_M2000to3000";
	
	//Data
	else if(_fileName==DATA_RunB)     filesave += "_RunB";
	else if(_fileName==DATA_RunC)     filesave += "_RunC";
        else if(_fileName==DATA_RunD)     filesave += "_RunD";
        else if(_fileName==DATA_RunE)     filesave += "_RunE";
        else if(_fileName==DATA_RunF)     filesave += "_RunF";
        else if(_fileName==DATA_RunG)     filesave += "_RunG";
        else if(_fileName==DATA_RunHver2) filesave += "_RunHver2";
        else if(_fileName==DATA_RunHver3) filesave += "_RunHver3";

	// Tops
	else if(_fileName==ST_tbarW)         filesave += "_ST_tbarW";
	else if(_fileName==ST_tW) 	     filesave += "_ST_tW";
	else if(_fileName==ttbar_M0to700)    filesave += "_ttbar_M0to700";
	else if(_fileName==ttbar_M700to1000) filesave += "_ttbar_M700to1000";
	else if(_fileName==ttbar_M1000toInf) filesave += "_ttbar_M1000toInf";

	// Fakes
	else if(_fileName==WJetsToLNu_amcatnlo)     filesave += "_WJetsToLNu";  
	else if(_fileName==WJetsToLNu_amcatnlo_ext) filesave += "_WJetsToLNuExt";

	// Dibosons
	else if(_fileName==WW) filesave += "_WW";
	else if(_fileName==WZ) filesave += "_WZ";
	else if(_fileName==ZZ) filesave += "_ZZ";

	else{
		cout << "_FileToLoad not properly defined" << endl;
		return;
	}

	if(_lepType==ELE) filesave += "_EE.root";
	else if(_lepType==MUON) filesave += "_MuMu.root";

	cout << "*********************************************************" << endl;
	cout << "Saving histograms in file: " << filesave << endl;
	cout << "*********************************************************" << endl;

	TFile*file = new TFile(filesave,"recreate");
	_hMassHardProcess->Write();
	_hRapidityHardProcess->Write();
	_hPtHardProcess->Write();
	_hMassReco->Write();
	_hRapidityReco->Write();
	_hPtReco->Write();
	file->Write();
	file->Close();
	return;
}

TString DrellYanAnalyzer::GetSampleName()
{
	using namespace DrellYanVariables;

	TString sampleName;

	// Dileptons
	if(_lepType==ELE){
		if(_fileName==DYLL_M10to50) 	     sampleName = "DYLL_M10to50_EE";
		else if(_fileName==DYLL_M50to100)    sampleName = "DYLL_M50to100_EE";
		else if(_fileName==DYLL_M100to200)   sampleName = "DYLL_M100to200_EE";
		else if(_fileName==DYLL_M200to400)   sampleName = "DYLL_M200to400_EE";
		else if(_fileName==DYLL_M400to500)   sampleName = "DYLL_M400to500_EE";
		else if(_fileName==DYLL_M500to700)   sampleName = "DYLL_M500to700_EE";
		else if(_fileName==DYLL_M700to800)   sampleName = "DYLL_M700to800_EE";
		else if(_fileName==DYLL_M800to1000)  sampleName = "DYLL_M800to1000_EE";
		else if(_fileName==DYLL_M1000to1500) sampleName = "DYLL_M1000to1500_EE";
		else if(_fileName==DYLL_M1500to2000) sampleName = "DYLL_M1500to2000_EE";
		else if(_fileName==DYLL_M2000to3000) sampleName = "DYLL_M2000to3000_EE";
		else if(_fileName==DATA_RunB)    sampleName = "crab_DoubleEG_RunB";
		else if(_fileName==DATA_RunC)    sampleName = "crab_DoubleEG_RunC";
		else if(_fileName==DATA_RunD)    sampleName = "crab_DoubleEG_RunD";
		else if(_fileName==DATA_RunE)    sampleName = "crab_DoubleEG_RunE";
		else if(_fileName==DATA_RunF)    sampleName = "crab_DoubleEG_RunF";
		else if(_fileName==DATA_RunG)    sampleName = "crab_DoubleEG_RunG";
		else if(_fileName==DATA_RunHver2)sampleName = "crab_DoubleEG_RunHver2";
		else if(_fileName==DATA_RunHver3)sampleName = "crab_DoubleEG_RunHver3";
	}//end if leptype = ELE
	else if(_lepType==MUON){
		if(_fileName==DYLL_M10to50)          sampleName = "DYLL_M10to50_MuMu";
		else if(_fileName==DYLL_M50to100)    sampleName = "DYLL_M50to100_MuMu";
		else if(_fileName==DYLL_M100to200)   sampleName = "DYLL_M100to200_MuMu";
		else if(_fileName==DYLL_M200to400)   sampleName = "DYLL_M200to400_MuMu";
		else if(_fileName==DYLL_M400to500)   sampleName = "DYLL_M400to500_MuMu";
		else if(_fileName==DYLL_M500to700)   sampleName = "DYLL_M500to700_MuMu";
		else if(_fileName==DYLL_M700to800)   sampleName = "DYLL_M700to800_MuMu";
		else if(_fileName==DYLL_M800to1000)  sampleName = "DYLL_M800to1000_MuMu";
		else if(_fileName==DYLL_M1000to1500) sampleName = "DYLL_M1000to1500_MuMu";
		else if(_fileName==DYLL_M1500to2000) sampleName = "DYLL_M1500to2000_MuMu";
		else if(_fileName==DYLL_M2000to3000) sampleName = "DYLL_M2000to3000_MuMu";
		else if(_fileName==DATA_RunB)    sampleName = "SingleMuon_Run2016B";
		else if(_fileName==DATA_RunC)    sampleName = "SingleMuon_Run2016C";
		else if(_fileName==DATA_RunD)    sampleName = "SingleMuon_Run2016D";
		else if(_fileName==DATA_RunE)    sampleName = "SingleMuon_Run2016E";
		else if(_fileName==DATA_RunF)    sampleName = "SingleMuon_Run2016F";
		else if(_fileName==DATA_RunG)    sampleName = "SingleMuon_Run2016G";
		else if(_fileName==DATA_RunHver2)sampleName = "SingleMuon_Run2016Hver2";
		else if(_fileName==DATA_RunHver3)sampleName = "SingleMuon_Run2016Hver3";
	}//end if leptype = MUON

	// Tops
	if(_fileName==ST_tbarW)		     sampleName = "ST_tbarW";
	else if(_fileName==ST_tW)	     sampleName = "ST_tW";
	else if(_fileName==ttbar_M0to700)    sampleName = "ttbar_M0to700";
	else if(_fileName==ttbar_M700to1000) sampleName = "ttbar_M700to1000";
	else if(_fileName==ttbar_M1000toInf) sampleName = "ttbar_M1000toInf";

	// Fakes
	else if(_fileName==WJetsToLNu_amcatnlo)     
		sampleName = "WJetsToLNu_amcatnlo";
	else if(_fileName==WJetsToLNu_amcatnlo_ext) 
		sampleName = "WJetsToLNu_amcatnloi_ext";

	// Dibosons
	else if(_fileName==WW) sampleName = "WW";
	else if(_fileName==WZ) sampleName = "WZ";
	else if(_fileName==ZZ) sampleName = "ZZ";

	// Taus
	if(_fileName==DYLL_M10to50_TauTau)          sampleName = "DYLL_M10to50_TauTau";
	else if(_fileName==DYLL_M50to100_TauTau)    sampleName = "DYLL_M50to100_TauTau";
	else if(_fileName==DYLL_M100to200_TauTau)   sampleName = "DYLL_M100to200_TauTau";
	else if(_fileName==DYLL_M200to400_TauTau)   sampleName = "DYLL_M200to400_TauTau";
	else if(_fileName==DYLL_M400to500_TauTau)   sampleName = "DYLL_M400to500_TauTau";
	else if(_fileName==DYLL_M500to700_TauTau)   sampleName = "DYLL_M500to700_TauTau";
	else if(_fileName==DYLL_M700to800_TauTau)   sampleName = "DYLL_M700to800_TauTau";
	else if(_fileName==DYLL_M800to1000_TauTau)  sampleName = "DYLL_M800to1000_TauTau";
	else if(_fileName==DYLL_M1000to1500_TauTau) sampleName = "DYLL_M1000to1500_TauTau";
	else if(_fileName==DYLL_M1500to2000_TauTau) sampleName = "DYLL_M1500to2000_TauTau";
	else if(_fileName==DYLL_M2000to3000_TauTau) sampleName = "DYLL_M2000to3000_TauTau";
	
	return sampleName;
}//end GetSampleName

double DrellYanAnalyzer::CalcVariable(double pt1,double eta1,double phi1,double mass1,
				      double pt2,double eta2,double phi2,double mass2,
				      DrellYanVariables::VarType varType)
{
	using namespace DrellYanVariables;

	TLorentzVector vLep1;
	TLorentzVector vLep2;
	vLep1.SetPtEtaPhiM(pt1,eta1,phi1,mass1);
	vLep2.SetPtEtaPhiM(pt2,eta2,phi2,mass2);

	if(varType==INV_MASS) return (vLep1+vLep2).M();
	else if(varType==RAPIDITY) return (vLep1+vLep2).Rapidity();
	else if(varType==PT) return (vLep1+vLep2).Pt();
	else if(varType==MUON_ANGLE) return (vLep1.Angle(vLep2.Vect()));
	else{
		cout << "ERROR in CalcVariable()!" << endl;
		cout << "VarType not properly chosen" << endl;
		cout << "See DrellYanAnalyzer.hh" << endl;
	       	return -100000000;
	}
}//end CalcVariable

double DrellYanAnalyzer::GetInvMass(double pt1,double eta1,double phi1,
		  	     	    double pt2,double eta2,double phi2)
{
	using namespace DrellYanVariables;
	double mass;
	if(_lepType==ELE) mass = eMass;
	else if(_lepType==MUON) mass = muMass;
	else{
		cout << "Lepton type not correctly chosen" << endl;
		return -10000;
	}
	double invMass = CalcVariable(pt1,eta1,phi1,mass,pt2,eta2,phi2,mass,INV_MASS);
	return invMass;
}

double DrellYanAnalyzer::GetRapidity(double pt1,double eta1,double phi1,
		  	     	     double pt2,double eta2,double phi2)
{
	using namespace DrellYanVariables;
	double mass;
	if(_lepType==ELE) mass = eMass;
	else if(_lepType==MUON) mass = muMass;
	else{
		cout << "Lepton type not correctly chosen" << endl;
		return -10000;
	}
	double rapidity = CalcVariable(pt1,eta1,phi1,mass,pt2,eta2,phi2,mass,RAPIDITY);
	return rapidity;
}

double DrellYanAnalyzer::GetPt(double pt1,double eta1,double phi1,
		  	       double pt2,double eta2,double phi2)
{
	using namespace DrellYanVariables;
	double mass;
	if(_lepType==ELE) mass = eMass;
	else if(_lepType==MUON) mass = muMass;
	else{
		cout << "Lepton type not correctly chosen" << endl;
		return -10000;
	}
	double pT = CalcVariable(pt1,eta1,phi1,mass,pt2,eta2,phi2,mass,PT);
	return pT;
}

void DrellYanAnalyzer::Counter(Long64_t i,Long64_t total){
	int P = 100*(i)/(total);  
	TTimeStamp eventTimeStamp;
	if(i%(total/100)==0)
	cout << "[Time: " << eventTimeStamp.AsString("s") << "] " << P << "%" << endl;
}
