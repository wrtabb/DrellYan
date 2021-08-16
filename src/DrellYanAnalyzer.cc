#include "DrellYanAnalyzer.hh"
///////////////////////////////////////////////
//-----Initialize DrellYanAnalyzer Class-----//
///////////////////////////////////////////////
DrellYanAnalyzer::DrellYanAnalyzer(DrellYanVariables::NtupleType ntupType,
				   DrellYanVariables::SampleType sampleType,
				   DrellYanVariables::LepType lepType,
				   TString sampleToLoad)
{
	using namespace DrellYanVariables;
	_sampleType = sampleType;
	_lepType = lepType;
	_ntupType = ntupType;
	TFile*puFile = new TFile("pileup.root");
        _hPileupRatio = (TH1F*)puFile->Get("hPileupRatio");
	_FileToLoad = base_directory_v2p6;
	_FileToLoad += sampleToLoad;
}//end function DrellYanAnalyzer()

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

int DrellYanAnalyzer::LoadTrees()
{
	using namespace DrellYanVariables;
	int returnCode = 1;
	TTimeStamp ts_start;
	cout << "Begin loading trees:" << endl;
	cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
	TStopwatch totaltime;
	totaltime.Start();

	_tree = new TChain(_treeName);
	_tree->Add(_FileToLoad);
	_nEvents = _tree->GetEvents();;

	totaltime.Stop();
	Double_t TotalCPURunTime = totaltime.CpuTime();
	Double_t TotalRunTime = totaltime.RealTime();
	TTimeStamp ts_end;
	cout << endl;
	cout << "End loading trees:" << endl;
	cout << "**************************************************************************" << endl;
	cout << "Total CPU RunTime: " << TotalCPURunTime << " seconds" << endl;
	cout << "Total Real RunTime: " << TotalRunTime << " seconds" << endl;
	cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;
	cout << "**************************************************************************" << endl;
	cout << endl;

	bool isMC = true;
	TBranch*testBranch = (TBranch*)_tree->
		GetListOfBranches()->FindObject("GENEvt_weight");
	if(!testBranch) isMC = false; 
	InitializeBranches(_tree,isMC);
	
	return returnCode;
}

int DrellYanAnalyzer::InitializeBranches(TChain*chain,bool isMC)
{
	using namespace DrellYanVariables;
	LepType lepType = _lepType;
	if(lepType!=MUON){
		//Electrons
		chain->SetBranchAddress("Nelectrons",&Nelectrons,&b_Nelectrons);
		chain->SetBranchAddress("Electron_pT",&Electron_pT,&b_Electron_pT);
		chain->SetBranchAddress("Electron_eta",&Electron_eta,&b_Electron_eta);
		chain->SetBranchAddress("Electron_phi",&Electron_phi,&b_Electron_phi);
		chain->SetBranchAddress("Electron_passMediumID",&Electron_passMediumID,
					 &b_Electron_passMediumID);
		chain->SetBranchAddress("vtxTrkCkt1Pt",&pvtxTrkCkt1Pt);
		chain->SetBranchAddress("vtxTrkCkt2Pt",&pvtxTrkCkt2Pt);
		chain->SetBranchAddress("vtxTrkChi2",&pvtxTrkChi2);
		chain->SetBranchAddress("vtxTrkNdof",&pvtxTrkNdof);
	}
	if(lepType!=ELE){
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

	if(isMC){
		//Gen leptons
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
		//Gen others
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
		//Gen weights
		chain->SetBranchAddress("GENEvt_weight",&GENEvt_weight,&b_GENEvt_weight);
	}
	return 1;
}//end initializing trees

void DrellYanAnalyzer::InitializeHistograms()
{
	using namespace DrellYanVariables;
	SampleType sampleType = _sampleType; //type of sample being processed
	LepType lepType = _lepType; //type of lepton
	int nSamples = _nSampleTypes; //number of samples being processed
	int nVars = 3; //number of variables to plot (mass, rapidity, pt, ...)
	int nProcs = 4;
	TString lepTag;
	vector<TString> sampleTag;

	if(sampleType==SAMPLE_LL) sampleTag.push_back("DYtoLL");
	else if(sampleType==SAMPLE_TOP) sampleTag.push_back("Tops");
	else if(sampleType==SAMPLE_FAKE) sampleTag.push_back("Fakes");
	else if(sampleType==SAMPLE_DIBOSON) sampleTag.push_back("Dibosons");
	else if(sampleType==SAMPLE_TAU) sampleTag.push_back("Taus");
	else if(sampleType==SAMPLE_DATA) sampleTag.push_back("Data");
	else if(sampleType==SAMPLE_ALL){
		sampleTag.push_back("DYtoLL");
		sampleTag.push_back("Tops");
		sampleTag.push_back("Fakes");
		sampleTag.push_back("Dibosons");
		sampleTag.push_back("Taus");
		sampleTag.push_back("Data");
	}
	else {
		cout << "**********************************" << endl;
		cout << "ERROR in InitializeHistograms()!!!" << endl;	
		cout << "Sample type not properly defined" << endl;
		return;
	}

	if(lepType==ELE) lepTag = "_EE_";
	else lepTag = "_MuMu_";

	vector<TString> varName = {
		"InvMass",
		"Rapidity",
		"Pt"
	};

	vector<TString> procTag = {
		"HardProcess",
		"FSR",
		"Reco",
		"Pruned"
	};

	_nHists = nSamples*nVars*nProcs;
	TH1D*hist[nSamples][nVars][nProcs];
	for(int i=0;i<nSamples;i++){
		for(int j=0;j<nVars;j++){
			for(int k=0;k<nProcs;k++){
				TString hName = "hist";
				hName += varName.at(j);
				hName += lepTag;
				hName += sampleTag.at(i);  
				hName += procTag.at(k);
				if(j==0)hist[i][j][k] = 
					new TH1D(hName,"",nMassBins,massbins);
				if(j==1)hist[i][j][k] = 
					new TH1D(hName,"",nRapidityBins,rapiditybins);
				if(j==2)hist[i][j][k] = 
					new TH1D(hName,"",nPtBins,ptbins);
				_hists.push_back(hist[i][j][k]);
			}//end loop over processes
		}//end loop over variables
	}//end loop over samples
}//end InitializeHistograms

/////////////////////////
//-----Get Weights-----//
/////////////////////////
double DrellYanAnalyzer::GetWeights(int index1,int index2)
{
	using namespace DrellYanVariables;
	Long64_t nEntries = _trees.at(index1).at(index2)->GetEntries();
	SampleType sampleType = _sampleType;
	//-----Cross section-----//
	vector<double> xSec;
	if(sampleType==SAMPLE_LL) xSec = xSec_LL;
	else if(sampleType==SAMPLE_TAU) xSec = xSec_LL;
	double xSecWeight = 1.0;
	if(_ntupType==V2P6 || SINGLE_FILE) xSecWeight = dataLuminosity*xSec.at(index2)/1.0;
	else if(_ntupType==TEST) xSecWeight = dataLuminosity*xSec.at(index2)/nEntries;

	//-----Pileup Weight-----//
	double puWeight = GetPUWeight();

	//-----Gen Weights-----//
	double sumGenWeight = GetGenWeightSum(index1,index2);
	double genWeight = 1.0;
	if(_ntupType==V2P6) genWeight = (GENEvt_weight/fabs(GENEvt_weight))/sumGenWeight;  

	//-----Total Weight-----//
	double weight = xSecWeight*genWeight*puWeight;
	return weight;
}

double DrellYanAnalyzer::GetGenWeightSum(int index1,int index2)
{
	using namespace DrellYanVariables;
	SampleType sampleType = _sampleType;

	double sumGenWeight = 0.0;
	double sumRawGenWeight = 0.0;
	double varGenWeight = 0.0;
	double genWeight;
	Long64_t nEntries = _trees.at(index1).at(index2)->GetEntries();
	Long64_t localEntry;
	for(Long64_t i=0;i<nEntries;i++){
		localEntry = _trees.at(index1).at(index2)->LoadTree(i);
		b_GENEvt_weight->GetEntry(localEntry);
		genWeight = GENEvt_weight/fabs(GENEvt_weight);//normalized genweight
		sumGenWeight += genWeight;
		varGenWeight += GENEvt_weight*GENEvt_weight;//variance of genweights
		sumRawGenWeight += GENEvt_weight;
	}
	return sumGenWeight;
}//end GetGenWeight

double DrellYanAnalyzer::GetPUWeight()
{
	using namespace DrellYanVariables;
	double puWeight = 1.0;	
	puWeight = _hPileupRatio->GetBinContent(_hPileupRatio->FindBin(nPileUp));
	return puWeight;
}

///////////////////////////////////////
//-----Get Leptons From an Event-----//
///////////////////////////////////////
int DrellYanAnalyzer::EventLoop()
{
	using namespace DrellYanVariables;
	LepType lepType = _lepType;
	SampleType sampleType = _sampleType;

	InitializeHistograms();
	int iHard1,iHard2;
	int iFSR1,iFSR2;
	int iDressed1,iDressed2;
	int iLep1,iLep2;
	double lMass;
	Long64_t eventCount = 0;
	if(lepType==ELE) lMass = eMass;
	else if(lepType==MUON) lMass = muMass;
	else{
		cout << "Error in EventLoop(): Need to specify ELE or MUON" << endl;
		 return 0;
	}
	//Loop over all samples, i.e. DYtoLL, Fakes, data, etc.
	for(int i=0;i<_nSampleTypes;i++){
		int nFiles = _trees.at(i).size();
		Long64_t totalEvents = _nEvents;
		bool isMC = true;
		vector<TString> fileList;
		if(_sampleType==SAMPLE_LL) fileList = _files_LL; 
		//Loop over files within a sample
		for(int j=0;j<nFiles;j++){
			cout << "Processing file: " << fileList.at(j) << endl;
			TBranch*testBranch = (TBranch*)_trees.at(i).at(j)->
				GetListOfBranches()->FindObject("GENEvt_weight");
			if(!testBranch) isMC = false; 
			Long64_t nentries = _trees.at(i).at(j)->GetEntries();
			
			//cross section weights
			vector<double> xSec;
                        if(sampleType==SAMPLE_LL) xSec = xSec_LL;
                        else if(sampleType==SAMPLE_TAU) xSec = xSec_LL;
			double xSecWeight = 1.0;
			if(_ntupType==TEST) 
				xSecWeight = dataLuminosity*xSec.at(j)/nentries;
			else xSecWeight = dataLuminosity*xSec.at(j)/1.0;

			//Gen weights
			double genWeight = 1.0;
			double sumGenWeight = 1.0;
			double genEvtWeight;
			double absGenEvtWeight;
			double evtWeightRatio;
			if(_ntupType!=TEST && isMC){ 
				sumGenWeight = GetGenWeightSum(i,j);
				genEvtWeight = GENEvt_weight;
				absGenEvtWeight = fabs(GENEvt_weight);
				evtWeightRatio = genEvtWeight/absGenEvtWeight;
				genWeight = evtWeightRatio/sumGenWeight;  
			}

			//Loop over all events in the loaded tree
			for(Long64_t iEntry=0;iEntry<nentries;iEntry++){
				eventCount++;

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

				//Get events from tree
				_trees.at(i).at(j)->GetEntry(iEntry);

				//Calculate pileup weight
				double puWeight = GetPUWeight();

				//Get gen level leptons
				int nDileptonsGen = GetGenLeptons(iHard1,iHard2,
					  	          	  iFSR1,iFSR2);
				//Get reco level leptons
				int nDileptonsReco = GetRecoLeptons(iLep1,iLep2);

				if(nDileptonsGen==1 && nDileptonsReco==1){
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

					if(lepType==ELE){
						recoPt1  = Electron_pT[iLep1];
						recoEta1 = Electron_eta[iLep1];
						recoPhi1 = Electron_phi[iLep1];
						recoPt2  = Electron_pT[iLep2];
						recoEta2 = Electron_eta[iLep2];
						recoPhi2 = Electron_phi[iLep2];
					}
					else if(lepType==MUON){
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

					//hard process quantities
					invMassHard = GetInvMass(hardPt1,hardEta1,hardPhi1,
							       	 hardPt2,hardEta2,hardPhi2);
					rapidityHard = GetRapidity(hardPt1,hardEta1,hardPhi1,
								   hardPt2,hardEta2,hardPhi2);
					ptHard = GetPt(hardPt1,hardEta1,hardPhi1,
						       hardPt2,hardEta2,hardPhi2);

					//FSR quantities
					invMassFSR = GetInvMass(fsrPt1,fsrEta1,fsrPhi1,
							       	fsrPt2,fsrEta2,fsrPhi2);
					rapidityFSR = GetRapidity(fsrPt1,fsrEta1,fsrPhi1,
								  fsrPt2,fsrEta2,fsrPhi2);
					ptFSR = GetPt(fsrPt1,fsrEta1,fsrPhi1,
						      fsrPt2,fsrEta2,fsrPhi2);

					//reco quantities
					invMassReco = GetInvMass(recoPt1,recoEta1,recoPhi1,
							       	 recoPt2,recoEta2,recoPhi2);
					rapidityReco = GetRapidity(recoPt1,recoEta1,recoPhi1,
								   recoPt2,recoEta2,recoPhi2);
					ptReco = GetPt(recoPt1,recoEta1,recoPhi1,
						       recoPt2,recoEta2,recoPhi2);

					//Pruned quantities
					invMassPruned = GetInvMass(prunedPt1,prunedEta1,
								   prunedPhi1,prunedPt2,
								   prunedEta2,prunedPhi2);
					rapidityPruned = GetRapidity(prunedPt1,prunedEta1,
								   prunedPhi1,prunedPt2,
								   prunedEta2,prunedPhi2);
					ptPruned = GetPt(prunedPt1,prunedEta1,
						   	 prunedPhi1,prunedPt2,
						         prunedEta2,prunedPhi2);
				}//end if nDileptons>0
			
				double weights = xSecWeight*genWeight*puWeight;
				vector<double> var = {
					invMassHard,
					invMassFSR,
					invMassReco,
					invMassPruned,
					rapidityHard,
					rapidityFSR,
					rapidityReco,
					rapidityPruned,
					ptHard,
					ptFSR,
					ptReco,
					ptPruned
				};
				int varSize = var.size();
				int l=0;

				//Loop over histograms
				for(int k=0;k<_nHists;k++){
					_hists.at(k)->Fill(var.at(l),weights);
					l++;
					if(l>varSize) l=0;
				}//end loop over histograms
				Counter(eventCount,totalEvents);
			}//end event loop
		}//end loop over files
	}//end loop over samples						
	SaveResults();
	return 1;
}//end EventLoop

int DrellYanAnalyzer::GetGenLeptons(int &iHard1,int &iHard2,
				    int &iFSR1,int &iFSR2)
{
	using namespace DrellYanVariables;
	LepType lepType = _lepType;
	int lepID = 0;
	if      (lepType==ELE)  lepID = 11;
	else if (lepType==MUON) lepID = 13;
	else {
		cout << "ERROR: Appropriate lepton not selected." << endl;
		cout << "NOTE: This analysis does not handle taus." << endl;
		return 0;
	}
	int nDileptons = 0;

	for(int iLep=0;iLep<GENnPair;iLep++){
		for(int jLep=iLep+1;jLep<GENnPair;jLep++){
			//Lepton selection
			if(!(abs(GENLepton_ID[iLep])==lepID&&abs(GENLepton_ID[jLep])==lepID)) 
				continue;
			//require opposite sign for electrons
			if(GENLepton_ID[iLep]*GENLepton_ID[jLep]>0&&lepID==11) 
				continue;
			if(GENLepton_isHardProcess[iLep]==1&&GENLepton_isHardProcess[jLep]==1){
				iHard1 = iLep;
				iHard2 = jLep;
				nDileptons++;
			}//end if hard process
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
	LepType lepType = _lepType;
	if(lepType==ELE) return GetRecoElectrons(leadLep,subLep);
	else if(lepType==MUON) return GetRecoMuons(leadLep,subLep);
	else {
		cout << "Must choose ELE or MUON for lepType" << endl;
		return 0;
	}
}

int DrellYanAnalyzer::GetRecoElectrons(int &leadEle,int &subEle)
{
	using namespace DrellYanVariables;
	int numDielectrons = 0;
	double pt1;
	double pt2;
	double eta1;
	double eta2;

	for(int iEle = 0; iEle < Nelectrons; iEle++) {
		if(!Electron_passMediumID[iEle]) continue;
		for(int jEle = iEle+1; jEle < Nelectrons; jEle++) {
			if(!Electron_passMediumID[jEle]) continue;
			pt1 = Electron_pT[iEle];
			pt2 = Electron_pT[jEle];
			eta1 = Electron_eta[iEle];
			eta2 = Electron_eta[jEle];
			if(PassAcceptance(pt1,pt2,eta1,eta2)){
				numDielectrons++;
				if(pt1>pt2){
					leadEle = iEle;
					subEle = jEle;
				}
				else {
					leadEle= jEle;
					subEle = iEle;
				}
			}
		}//end jEle loop
	}//end iEle loop

	return numDielectrons;
}//end GetRecoElecrons()

int DrellYanAnalyzer::GetRecoMuons(int &leadMu,int &subMu)
{
	using namespace DrellYanVariables;
	int numDimuons = 0;
	double pt1;
	double pt2;
	double eta1;
	double eta2;
	double phi1;
	double phi2;
	double charge1;
	double charge2;

	for(int iMu=0;iMu<nMuon;iMu++){
		if(!Muon_passTightID[iMu]) continue;
		if(!PassMuonIsolation(iMu)) continue;
		for(int jMu=iMu+1;jMu<nMuon;jMu++){
			pt1 = Muon_pT[iMu];
			pt2 = Muon_pT[jMu];
			eta1 = Muon_eta[iMu];
			eta2 = Muon_eta[jMu];
			phi1 = Muon_phi[iMu];
			phi2 = Muon_phi[jMu];
			charge1 = Muon_charge[iMu];
			charge2 = Muon_charge[jMu];

			if(!Muon_passTightID[jMu]) continue;
			if(!PassMuonIsolation(jMu)) continue;
			if(!PassMuonAngle(pt1,eta1,phi1,muMass,pt2,eta2,phi2,muMass)) 
				continue;
			if(charge1*charge2>0) continue;

			//Other cuts not yet applied:
			//smallest dimuon vertex chi2
			//PassVertexChi2(iMu,jMu);
			if(PassAcceptance(pt1,pt2,eta1,eta2)){
				numDimuons++;
				if(pt1>pt2){
					leadMu = iMu;
					subMu = jMu;
				}
				else {
					leadMu= jMu;
					subMu = iMu;
				}
			}
		}//end jMu loop
	}//end iMu loop

	return numDimuons;
}//end GetRecoLeptons()

//////////////////////////
//-----Cut criteria-----//
//////////////////////////
bool DrellYanAnalyzer::PassAcceptance(double pt1,double pt2,double eta1,double eta2)
{
	using namespace DrellYanVariables;
	//Acceptance criteria is the same for muons and electrons
	
	//I think we are not excluding the gap region in the ECAL, so I commented these out
	//But the lines can be uncommented if we want to exclude the gap
	//if(abs(eta1)>etaGapLow && abs(eta1)<etaGapHigh) return false;
	//if(abs(eta2)>etaGapLow && abs(eta2)<etaGapHigh) return false;
	
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
	LepType lepType = _lepType;
	TString trigName;
	TString triggerUsed;
	int trigNameSize = pHLT_trigName->size();
	bool passHLT = false;
	for(int iHLT=0;iHLT<trigNameSize;iHLT++) {
		trigName = pHLT_trigName->at(iHLT);
		if(lepType==ELE){
			if(trigName.CompareTo(electronTrigger)==0 && HLT_trigFired[iHLT]==1){
				passHLT = true;
				break;
			}//end if trigName...
		}//end if lepType
		else if(lepType==MUON){
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
	LepType lepType = _lepType;
	if(lepType==ELE) return PassGenToRecoMatchEle(genIndex,recoIndex);
	if(lepType==MUON) return PassGenToRecoMatchMu(genIndex,recoIndex);
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
	LepType lepType = _lepType;
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
	double pfIso_dBeta;
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
	TString filesave = "DYHists";
	if(_ntupType==V2P6) filesave+= "_v2p6";
	else if (_ntupType==V2P3) filesave+= "_v2p3";
	else if(_ntupType==TEST) filesave += "_TEST";
	else if(_ntupType==SINGLE_TEST) filesave += "_SingleSampleTest";
	else{
		cout << "********************************************************" << endl;
		cout << "ERROR in SaveResults()" << endl;
		cout << "NtupleType not correctly defined" << endl;
		cout << "See include/DrellYanVariables.h for list of NtupleTypes" << endl;
		cout << "********************************************************" << endl;
		return;
	}

	if(_sampleType==SAMPLE_LL) filesave += "_DYtoLL";
	else if(_sampleType==SAMPLE_TOP) filesave += "_TT";
	else if(_sampleType==SAMPLE_FAKE) filesave += "_Fakes";
	else if (_sampleType==SAMPLE_DIBOSON) filesave += "_Diboson";
	else if (_sampleType==SAMPLE_TAU) filesave += "_TauTau";
	else if (_sampleType==SAMPLE_DATA) filesave += "_Data";
	else if (_sampleType==SAMPLE_ALL) filesave += "_All";
	else {
		cout << "*****************************************************" << endl;
		cout << "ERROR in SaveResults() " << endl;	
		cout << "SampleType not correctly defined" << endl;
		cout << "See include/DrellYanVariables.h for list of SampleTypes" << endl;
		cout << "*****************************************************" << endl;
	}


	if(_lepType==ELE) filesave += "_EE.root";
	else if(_lepType==MUON) filesave += "_MuMu.root";

	cout << "*********************************************************" << endl;
	cout << "Saving histograms in file: " << filesave << endl;
	cout << "*********************************************************" << endl;

	TFile*file = new TFile(filesave,"recreate");
	for(int i=0;i<_nHists;i++){
		_hists.at(i)->Write();
	}
	file->Write();
	file->Close();
	return;
}

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
