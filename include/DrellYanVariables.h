
namespace DrellYanVariables
{
	TString treeName = "recoTree/DYTree";

	//Ntuple version
	enum NtupleType{
		V2P3,
		V2P6,
		TEST,
		SINGLE_TEST
	};
	
	//Lepton selection
	enum LepType{
		ELE,
		MUON,
		NONE
	};

	//Sample Types
	enum SampleType{
		SAMPLE_LL,
		SAMPLE_TOP,
		SAMPLE_FAKE,
		SAMPLE_DIBOSON,
		SAMPLE_TAU,
		SAMPLE_DATA,
		SAMPLE_ALL
	};

	//Variable Types for placing in histograms
	enum VarType{
		INV_MASS,
		RAPIDITY,
		PT
	};
	
	//Process types
	enum ProcessType{
		HARD,
		FSR,
		RECO,
		PRUNED
	};

	//-----V2.6-----//	
	TString base_directory_v2p6 = "root://xrootd-local.unl.edu///store/user/wtabb/DrellYan_13TeV_2016/v2p6/skims/";

	vector<TString> dy_EE_v2p6 = {
		//Drell-Yan MC signal
		//Electrons
		"DYLL_M10to50_EE.root",					
		"DYLL_M50to100_EE.root",	
		"DYLL_M100to200_EE.root",				
		"DYLL_M200to400_EE.root",				
		"DYLL_M400to500_EE.root",				
		"DYLL_M500to700_EE.root",				
		"DYLL_M700to800_EE.root",				
		"DYLL_M800to1000_EE.root",				
		"DYLL_M1000to1500_EE.root",				
		"DYLL_M1500to2000_EE.root",				
		"DYLL_M2000to3000_EE.root"				
	};
	vector<TString> dy_MuMu_v2p6 = {
		//Muons
		"DYLL_M10to50_MuMu.root",				
		"DYLL_M50to100_MuMu.root",	
		"DYLL_M100to200_MuMu.root",				
		"DYLL_M200to400_MuMu.root",				
		"DYLL_M400to500_MuMu.root",				
		"DYLL_M500to700_MuMu.root",				
		"DYLL_M700to800_MuMu.root",				
		"DYLL_M800to1000_MuMu.root",				
		"DYLL_M1000to1500_MuMu.root",				
		"DYLL_M1500to2000_MuMu.root",				
		"DYLL_M2000to3000_MuMu.root"				
	};
	vector<TString> tops_v2p6 = {
		//Backgrounds
		//tops
		"ST_tbarW.root",					
		"ST_tW.root",						
		"ttbar_truncated_M0To700.root",				
		"ttbar_M700to1000.root",				
		"ttbar_M1000toInf.root"				
	};
	vector<TString> fakes_v2p6 = {
		//W+Jets (Fakes)		
		"WJetsToLNu_amcatnlo_ext2v5.root",			
		"WJetsToLNu_amcatnlo_ext.root",				
		"WJetsToLNu_amcatnlo.root"				
	};
	vector<TString> dibosons_v2p6 = {
		//EW (bosons + DY->TauTau)
		//Bosons
		"WW.root",						
		"WZ.root",						
		"ZZ.root"						
	};
	vector<TString> taus_v2p6 = {
		//Taus
		"DYLL_M10to50_TauTau.root",				
		"DYLL_M50to100_TauTau.root",	
		"DYLL_M100to200_TauTau.root",				
		"DYLL_M200to400_TauTau.root",				
		"DYLL_M400to500_TauTau.root",				
		"DYLL_M500to700_TauTau.root",				
		"DYLL_M700to800_TauTau.root",				
		"DYLL_M800to1000_TauTau.root",				
		"DYLL_M1000to1500_TauTau.root",				
		"DYLL_M1500to2000_TauTau.root",				
		"DYLL_M2000to3000_TauTau.root"				
	};
	vector<TString> data_v2p6 = {
		//Drell-Yan Data
		"crab_DoubleEG_RunB.root",					
		"crab_DoubleEG_RunC.root",					
		"crab_DoubleEG_RunD.root",					
		"crab_DoubleEG_RunE.root",					
		"crab_DoubleEG_RunF.root",					
		"crab_DoubleEG_RunG.root",					
		"crab_DoubleEG_RunHver2.root",				
		"crab_DoubleEG_RunHver3.root"				
	};

	//-----TEST-----//
	TString base_directory_test = "root://xrootd-local.unl.edu///store/user/wtabb/DrellYan_13TeV_2016/v2p6/skims_testSkim/";

	vector<TString> dy_EE_test = {
		//Drell-Yan MC signal
		//Electrons
		"DYLL_M10to50_EE.root",					
		"DYLL_M50to100_EE.root",	
		"DYLL_M100to200_EE.root",				
		"DYLL_M200to400_EE.root",				
		"DYLL_M400to500_EE.root",				
		"DYLL_M500to700_EE.root",				
		"DYLL_M700to800_EE.root",				
		"DYLL_M800to1000_EE.root",				
		"DYLL_M1000to1500_EE.root",				
		"DYLL_M1500to2000_EE.root",				
		"DYLL_M2000to3000_EE.root"				
	};
	vector<TString> dy_MuMu_test = {
		//Muons
		"DYLL_M10to50_MuMu.root",				
		"DYLL_M50to100_MuMu.root",	
		"DYLL_M100to200_MuMu.root",				
		"DYLL_M200to400_MuMu.root",				
		"DYLL_M400to500_MuMu.root",				
		"DYLL_M500to700_MuMu.root",				
		"DYLL_M700to800_MuMu.root",				
		"DYLL_M800to1000_MuMu.root",				
		"DYLL_M1000to1500_MuMu.root",				
		"DYLL_M1500to2000_MuMu.root",				
		"DYLL_M2000to3000_MuMu.root"				
	};

	vector<TString> dy_EE_SingleTest = {
		//Drell-Yan MC signal
		//Electrons
		"DYLL_M10to50_EE.root",					
	//	"DYLL_M50to100_EE.root",	
	//	"DYLL_M100to200_EE.root",				
	//	"DYLL_M200to400_EE.root",				
	//	"DYLL_M400to500_EE.root",				
	//	"DYLL_M500to700_EE.root",				
	//	"DYLL_M700to800_EE.root",				
	//	"DYLL_M800to1000_EE.root",				
	//	"DYLL_M1000to1500_EE.root",				
	//	"DYLL_M1500to2000_EE.root",				
	//	"DYLL_M2000to3000_EE.root"				
	};
	vector<TString> dy_MuMu_SingleTest = {
		//Muons
		"DYLL_M10to50_MuMu.root",				
	//	"DYLL_M50to100_MuMu.root",	
	//	"DYLL_M100to200_MuMu.root",				
	//	"DYLL_M200to400_MuMu.root",				
	//	"DYLL_M400to500_MuMu.root",				
	//	"DYLL_M500to700_MuMu.root",				
	//	"DYLL_M700to800_MuMu.root",				
	//	"DYLL_M800to1000_MuMu.root",				
	//	"DYLL_M1000to1500_MuMu.root",				
	//	"DYLL_M1500to2000_MuMu.root",				
	//	"DYLL_M2000to3000_MuMu.root"				
	};

	vector<double> xSec_LL = {
		18610.0,	//DYLL_10to50 v1,v2,ext1v1 combined (NLO)
		1923.26*3,	//DYLL_50to100(NNLO)
		78.1258*3,	//DYLL_100to200(NNLO)
		2.73309*3,	//DYLL_200to400(NNLO) 
		0.142945*3,	//DYLL_400to500(NNLO)
		0.0809755*3,	//DYLL_500to700(NNLO)
		0.0125589*3,	//DYLL_700to800(NNLO)
		0.0105845*3,	//DYLL_800to1000(NNLO)
		0.00556507*3,	//DYLL_1000to1500(NNLO)
		0.000730495*3,	//DYLL_1500to2000(NNLO)
		0.00016844*3	//DYLL_2000to3000(NNLO)
	};
	vector<double> ntuple_xsec_tops = {
		35.85,	//ST_tbarW(NNLO)
		35.85,	//ST_tW(NNLO)
		//Need to calculate cross section of M0to700 from full range
		// xsec_ttbar = 831.76
		// nEntriesTotal = 76949785
		// nEntriesM0to700 = 
		1.0,//ttbar_truncated_M0To700				
		76.605,	//ttbar_M700to1000(NNLO)	
		20.578	//ttbar_M1000toInf(NNLO)
	};

	//Cut criteria
	const double etaGapLow = 1.4442;
	const double etaGapHigh = 1.566;
	const double etaHigh = 2.4;
	const double ptLow = 17;
	const double ptHigh = 28;
	const float dRMinCut = 0.3;
	
	const TString electronTrigger = "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*";
	const TString muonTrigger1 = "HLT_IsoMu24_v*";
	const TString muonTrigger2 = "HLT_IsoTkMu24_v*";

	const int MPSIZE = 2000;
	int nPileUp;
	double PVz;
	const int dataLuminosity = 35867; //Run2016B to Run2016H JSON. unit: /pb, Updated at 2017.07.30

	std::vector<double> vtxTrkCkt1Pt;
	std::vector<double>*pvtxTrkCkt1Pt = &vtxTrkCkt1Pt;

	std::vector<double> vtxTrkCkt2Pt;
	std::vector<double>*pvtxTrkCkt2Pt = &vtxTrkCkt2Pt;

	std::vector<double> vtxTrkChi2;
	std::vector<double>*pvtxTrkChi2 = &vtxTrkChi2;	

	std::vector<double> vtxTrkNdof;
	std::vector<double>*pvtxTrkNdof = &vtxTrkNdof;

	//Gen lepton variables
	int GENnPair;//number of gen leptons per event
	double GENEvt_weight;
	double GENLepton_phi[MPSIZE];
	double GENLepton_eta[MPSIZE];
	double GENLepton_pT[MPSIZE];
	double GENLepton_Px[MPSIZE];
	double GENLepton_Py[MPSIZE];
	double GENLepton_Pz[MPSIZE];
	double GENLepton_E[MPSIZE];
	int GENLepton_ID[MPSIZE];
	bool GENLepton_isHardProcess[MPSIZE];
	bool GENLepton_fromHardProcessFinalState[MPSIZE];

	//Gen others variables
	int nGenOthers;
	double GenOthers_phi[MPSIZE];
	double GenOthers_eta[MPSIZE];
	double GenOthers_pT[MPSIZE];
	double GenOthers_Px[MPSIZE];
	double GenOthers_Py[MPSIZE];
	double GenOthers_Pz[MPSIZE];
	double GenOthers_E[MPSIZE];
	int GenOthers_ID[MPSIZE];
	int GenOthers_isHardProcess[MPSIZE];
	int GenOthers_isPromptFinalState[MPSIZE];

	//Electron variables
	int Nelectrons;
	double Electron_Energy[MPSIZE];  //no muon
	double Electron_pT[MPSIZE];
	double Electron_Px[MPSIZE];
	double Electron_Py[MPSIZE];
	double Electron_Pz[MPSIZE];
	double Electron_eta[MPSIZE];
	double Electron_phi[MPSIZE];
	int Electron_charge[MPSIZE];
	double Electron_etaSC[MPSIZE]; //no muon
	double Electron_phiSC[MPSIZE]; //no muon
	double Electron_dxy[MPSIZE];
	double Electron_dz[MPSIZE];
	double Electron_EnergySC[MPSIZE]; //no muon
	double Electron_etSC[MPSIZE]; //no muon
	bool Electron_passMediumID[MPSIZE];
	double eMass = 0.000511;

	//Muon variables
	int nMuon;
	int Nmuons;
	double Muon_pT[MPSIZE];
	double Muon_Px[MPSIZE];
	double Muon_Py[MPSIZE];
	double Muon_Pz[MPSIZE];
	double Muon_eta[MPSIZE];
	double Muon_phi[MPSIZE];
	int Muon_charge[MPSIZE];
	double Muon_dxy[MPSIZE];
	double Muon_dz[MPSIZE];
	bool Muon_passTightID[MPSIZE];
	double Muon_PfChargedHadronIsoR04[MPSIZE];
	double Muon_PfNeutralHadronIsoR04[MPSIZE];
	double Muon_PfGammaIsoR04[MPSIZE];
	double Muon_PFSumPUIsoR04[MPSIZE];
	double Muon_trkiso[MPSIZE];
	double muMass = 0.105658;

	//HLT variables
	int HLT_ntrig;
	int HLT_trigType[MPSIZE];
	int HLT_trigFired[MPSIZE];
	std::vector<std::string> HLT_trigName;
	std::vector<std::string> *pHLT_trigName = &HLT_trigName;

	int nVertices;

	//-----Branches-----//
	TBranch*b_runNum;
	TBranch*b_evtNum;
	TBranch*b_lumiBlock;
	TBranch*b_PUweight;
	TBranch*b_Nelectrons;
	TBranch*b_nVertices;
	TBranch*b_nPileUp;

	TBranch*b__prefiringweight;
	TBranch*b__prefiringweightup;
	TBranch*b__prefiringweightdown;
	
	//Trigger
	TBranch*b_HLT_ntrig;
	TBranch*b_HLT_trigType;
	TBranch*b_HLT_trigFired;

	//Muons
	TBranch*b_nMuon;
	TBranch*b_Nmuons;
	TBranch*b_PVz;
	TBranch*b_Muon_pT;
	TBranch*b_Muon_Px;
	TBranch*b_Muon_Py;
	TBranch*b_Muon_Pz;
	TBranch*b_Muon_eta;
	TBranch*b_Muon_phi;
	TBranch*b_Muon_charge;
	TBranch*b_Muon_dxy;
	TBranch*b_Muon_dz;
	TBranch*b_Muon_passTightID;
	TBranch*b_Muon_PfChargedHadronIsoR04;
	TBranch*b_Muon_PfNeutralHadronIsoR04;
	TBranch*b_Muon_PfGammaIsoR04;
	TBranch*b_Muon_PFSumPUIsoR04;
	TBranch*b_Muon_trkiso;

	//Electrons
	TBranch*b_Electron_Energy;
	TBranch*b_Electron_pT;
	TBranch*b_Electron_Px;
	TBranch*b_Electron_Py;
	TBranch*b_Electron_Pz;
	TBranch*b_Electron_eta;
	TBranch*b_Electron_phi;
	TBranch*b_Electron_charge;
	TBranch*b_Electron_etaSC;
	TBranch*b_Electron_phiSC;
	TBranch*b_Electron_dxy;
	TBranch*b_Electron_dz;
	TBranch*b_Electron_EnergySC;
	TBranch*b_Electron_etSC;
	TBranch*b_Electron_passMediumID;

	//Gen leptons
	TBranch*b_GENnPair;
	TBranch*b_GENLepton_phi;
	TBranch*b_GENLepton_eta;
	TBranch*b_GENLepton_pT;
	TBranch*b_GENLepton_Px;
	TBranch*b_GENLepton_Py;
	TBranch*b_GENLepton_Pz;
	TBranch*b_GENLepton_E;
	TBranch*b_GENLepton_mother;
	TBranch*b_GENLepton_mother_pT;
	TBranch*b_GENLepton_charge;
	TBranch*b_GENLepton_status;
	TBranch*b_GENLepton_ID;
	TBranch*b_GENLepton_isPrompt;
	TBranch*b_GENLepton_isPromptFinalState;
	TBranch*b_GENLepton_isTauDecayProduct;
	TBranch*b_GENLepton_isPromptTauDecayProduct;
	TBranch*b_GENLepton_isDirectPromptTauDecayProductFinalState;
	TBranch*b_GENLepton_isHardProcess;
	TBranch*b_GENLepton_isLastCopy;
	TBranch*b_GENLepton_isLastCopyBeforeFSR;
	TBranch*b_GENLepton_isPromptDecayed;
	TBranch*b_GENLepton_isDecayedLeptonHadron;
	TBranch*b_GENLepton_fromHardProcessBeforeFSR;
	TBranch*b_GENLepton_fromHardProcessDecayed;
	TBranch*b_GENLepton_fromHardProcessFinalState;
	TBranch*b_GENLepton_isMostlyLikePythia6Status3;
	TBranch*b_GENEvt_weight;
	TBranch*b_GENEvt_QScale;
	TBranch*b_GENEvt_x1;
	TBranch*b_GENEvt_x2;
	TBranch*b_GENEvt_alphaQCD;
	TBranch*b_GENEvt_alphaQED;

	//Gen others
	TBranch*b_nGenOthers;
	TBranch*b_GenOthers_phi;
	TBranch*b_GenOthers_eta;
	TBranch*b_GenOthers_pT;
	TBranch*b_GenOthers_Px;
	TBranch*b_GenOthers_Py;
	TBranch*b_GenOthers_Pz;
	TBranch*b_GenOthers_E ;
	TBranch*b_GenOthers_ID;
	TBranch*b_GenOthers_isHardProcess;
	TBranch*b_GenOthers_isPromptFinalState;



	vector<TString> file_locations_v2p6 = {
		//Drell-Yan MC signal
		//Electrons
		"DYLL_M10to50_EE.root",					//0
		"DYLL_M50toInf_truncated_M50To100_EE_base.root",	//1
		"DYLL_M100to200_EE.root",				//2
		"DYLL_M200to400_EE.root",				//3
		"DYLL_M400to500_EE.root",				//4
		"DYLL_M500to700_EE.root",				//5
		"DYLL_M700to800_EE.root",				//6
		"DYLL_M800to1000_EE.root",				//7
		"DYLL_M1000to1500_EE.root",				//8
		"DYLL_M1500to2000_EE.root",				//9
		"DYLL_M2000to3000_EE.root",				//10

		//Muons
		"DYLL_M10to50_MuMu.root",				//11
		"DYLL_M50toInf_truncated_M50To100_MuMu_base.root",	//12
		"DYLL_M100to200_MuMu.root",				//13
		"DYLL_M200to400_MuMu.root",				//14
		"DYLL_M400to500_MuMu.root",				//15
		"DYLL_M500to700_MuMu.root",				//16
		"DYLL_M700to800_MuMu.root",				//17
		"DYLL_M800to1000_MuMu.root",				//18
		"DYLL_M1000to1500_MuMu.root",				//19
		"DYLL_M1500to2000_MuMu.root",				//20
		"DYLL_M2000to3000_MuMu.root",				//21
		
		//Backgrounds
		//tops
		"ST_tbarW.root",					//22
		"ST_tW.root",						//23
		"ttbar_truncated_M0To700.root",				//24
		"ttbar_M700to1000.root",				//25
		"ttbar_M1000toInf.root",				//26
		//W+Jets (Fakes)		
		"WJetsToLNu_amcatnlo_ext2v5.root",			//27
		"WJetsToLNu_amcatnlo_ext.root",				//28
		"WJetsToLNu_amcatnlo.root",				//29
		//EW (bosons + DY->TauTau)
		//Bosons
		"WW.root",						//30
		"WZ.root",						//31
		"ZZ.root",						//32
		//Taus
		"DYLL_M10to50_TauTau.root",				//33
		"DYLL_M50toInf_truncated_M50To100_TauTau_base.root",	//34
		"DYLL_M100to200_TauTau.root",				//35
		"DYLL_M200to400_TauTau.root",				//36
		"DYLL_M400to500_TauTau.root",				//37
		"DYLL_M500to700_TauTau.root",				//38
		"DYLL_M700to800_TauTau.root",				//39
		"DYLL_M800to1000_TauTau.root",				//40
		"DYLL_M1000to1500_TauTau.root",				//41
		"DYLL_M1500to2000_TauTau.root",				//42
		"DYLL_M2000to3000_TauTau.root",				//43

		//Drell-Yan Data
		"DoubleEG_RunB.root",					//44
		"DoubleEG_RunC.root",					//
		"DoubleEG_RunD.root",					//
		"DoubleEG_RunE.root",					//
		"DoubleEG_RunF.root",					//
		"DoubleEG_RunG.root",					//
		"DoubleEG_RunHver2.root",				//
		"DoubleEG_RunHver3.root",				//
	};

	enum FileType{
		DYLL_M10to50_EE,				
		DYLL_M50to100_EE,	
		DYLL_M100to200_EE,				
		DYLL_M200to400_EE,				
		DYLL_M400to500_EE,				
		DYLL_M500to700_EE,				
		DYLL_M700to800_EE,				
		DYLL_M800to1000_EE,				
		DYLL_M1000to1500_EE,				
		DYLL_M1500to2000_EE,				
		DYLL_M2000to3000_EE,				
		DYLL_M10to50_MuMu,			
		DYLL_M50to100_MuMu,	
		DYLL_M100to200_MuMu,				
		DYLL_M200to400_MuMu,				
		DYLL_M400to500_MuMu,				
		DYLL_M500to700_MuMu,				
		DYLL_M700to800_MuMu,				
		DYLL_M800to1000_MuMu,				
		DYLL_M1000to1500_MuMu,				
		DYLL_M1500to2000_MuMu,				
		DYLL_M2000to3000_MuMu,				
		ST_tbarW,					
		ST_tW,						
		ttbar_M0To700,			   
		ttbar_M700to1000,				
		ttbar_M1000toInf,				
		WJetsToLNu_amcatnlo_ext2v5,			
		WJetsToLNu_amcatnlo_ext,			   
		WJetsToLNu_amcatnlo,				
		WW,						
		WZ,						
		ZZ,						
		DYLL_M10to50_TauTau_ext1v1,			
		DYLL_M10to50_TauTau_v1,				
		DYLL_M10to50_TauTau_v2,				
		DYLL_M50toInf_truncated_M50To100_TauTau_base,	
		DYLL_M100to200_TauTau,				
		DYLL_M200to400_TauTau,				
		DYLL_M400to500_TauTau,				
		DYLL_M500to700_TauTau,				
		DYLL_M700to800_TauTau,				
		DYLL_M800to1000_TauTau,				
		DYLL_M1000to1500_TauTau,			   
		DYLL_M1500to2000_TauTau,			   
		DYLL_M2000to3000_TauTau,			   
		DoubleEG_RunB,					
		DoubleEG_RunC,					
		DoubleEG_RunD,					
		DoubleEG_RunE,					
		DoubleEG_RunF,					
		DoubleEG_RunG,					
		DoubleEG_RunHver2,				
		DoubleEG_RunHver3				
	};

	//-----Histogram Parameters-----//
	double massbins[] = {
		15,
		20,
		25,
		30,
		35,
		40,
		45,
		50,
		55,
		60,
		64,
		68,
		72,
		76,
		81,
		86,
		91,
		96,
		101,
		106,
		110,
		115,
		120,
		126,
		133,
		141,
		150,
		160,
		171,
		185,
		200,
		220,
		243,
		273,
		320,
		380,
		440,
		510,
		600,
		700,
		830,
		1000,
		1500,
		3000
	};
	int nMassBins = size(massbins)-1;//43;
	double rapiditybins[] = {
		-2.4,
		-2.3,
		-2.2,
		-2.1,
		-2.0,
		-1.9,
		-1.8,
		-1.7,
		-1.6,
		-1.5,
		-1.4,
		-1.3,
		-1.2,
		-1.1,
		-1.0,
		-0.9,
		-0.8,
		-0.7,
		-0.6,
		-0.5,
		-0.4,
		-0.3,
		-0.2,
		-0.1,
		0.0,
		0.1,
		0.2,
		0.3,
		0.4,
		0.5,
		0.6,
		0.7,
		0.8,
		0.9,
		1.0,
		1.1,
		1.2,
		1.3,
		1.4,
		1.5,
		1.6,
		1.7,
		1.8,
		1.9,
		2.0,
		2.1,
		2.2,
		2.3,
		2.4
	};
	int nRapidityBins = size(rapiditybins)-1;//49;

	double ptbins[] = {
		0,
		5,
		10,
		15,
		20,
		25,
		30,
		35,
		40,
		45,
		50,
		55,
		60,
		65,
		70,
		75,
		80,
		85,
		90,
		95,
		100,
		105,
		110,
		115,
		120,
		125,
		130,
		135,
		140,
		145,
		150,
		155,
		160,
		165,
		170,
		175,
		180,
		185,
		190,
		195,
		200,
		205,
		210,
		215,
		220,
		225,
		230,
		235,
		240,
		245,
		250,
		255,
		260,
		265,
		270,
		275,
		280,
		285,
		290,
		295,
		300,
		305,
		310,
		315,
		320,
		325,
		330,
		335,
		340,
		345,
		350,
		355,
		360,
		365,
		370,
		375,
		380,
		385,
		390,
		395,
		400,
		405,
		410,
		415,
		420,
		425,
		430,
		435,
		440,
		445,
		450,
		455,
		460,
		465,
		470,
		475,
		480,
		485,
		490,
		495,
		500,
	};
	int nPtBins = size(ptbins)-1;//100

	vector<TString> _histTypes = {
                "histInvMass",
                "histRapidity",
                "histPt"
        };
        int _nHistTypes = 4*_histTypes.size();
}
