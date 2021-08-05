#include "../include/DrellYanAnalyzer.hh"
///////////////////////////////////////////////
//-----Initialize DrellYanAnalyzer Class-----//
///////////////////////////////////////////////
DrellYanAnalyzer::DrellYanAnalyzer(DrellYanVariables::NtupleType ntupType,
				   DrellYanVariables::SampleType sampleType,
				   DrellYanVariables::LepType lepType)
{
	using namespace DrellYanVariables;
	_sampleType = sampleType;
	_lepType = lepType;
	_ntupType = ntupType;
	TFile*puFile = new TFile("data/pileup.root");
        _hPileupRatio = (TH1F*)puFile->Get("hPileupRatio");

	if(ntupType==V2P6){
		_base_directory = base_directory_v2p6;
		if(lepType==ELE) _files_LL = dy_EE_v2p6;
		else if(lepType==MUON) _files_LL = dy_MuMu_v2p6;

		_files_tops = tops_v2p6;
		_files_fakes = fakes_v2p6;
		_files_dibosons = dibosons_v2p6;
		_files_taus = taus_v2p6;
		_files_data = data_v2p6;

		if(sampleType==SAMPLE_LL) _files.push_back(_files_LL);
		else if(sampleType==SAMPLE_TOP) _files.push_back(_files_tops);
		else if(sampleType==SAMPLE_FAKE) _files.push_back(_files_fakes);
		else if(sampleType==SAMPLE_DIBOSON) _files.push_back(_files_dibosons);
		else if(sampleType==SAMPLE_TAU) _files.push_back(_files_taus);
		else if(sampleType==SAMPLE_DATA) _files.push_back(_files_data);
		else if(sampleType==SAMPLE_ALL){
			_files.push_back(_files_LL);
			_files.push_back(_files_tops);
			_files.push_back(_files_fakes);
			_files.push_back(_files_dibosons);
			_files.push_back(_files_taus);
			_files.push_back(_files_data);
		}
		else {
			cout << "*************************************" << endl;
			cout << "* SaupleType must be chosen!        *" << endl;
			cout << "* See ./include/DrellYanVariables.h *" << endl;
			cout << "*************************************" << endl;
			return;
		}
	}//end if ntupType
	else if(ntupType==TEST){
		_base_directory = base_directory_test;
		if(lepType==ELE) _files_LL = dy_EE_test;
		else if(lepType==MUON) _files_LL = dy_MuMu_test;
		_files.push_back(_files_LL);
	}//end if ntupType
	else if(ntupType==SINGLE_TEST){
		_base_directory = base_directory_v2p6;
		if(lepType==ELE) _files_LL = dy_EE_SingleTest;
		else if(lepType==MUON) _files_LL = dy_MuMu_SingleTest;
		_files.push_back(_files_LL);
	}//end if ntupType
	_nSampleTypes = _files.size();
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

	int nSamples = _files.size();
	cout << "Number of samples to load: " << nSamples << endl;
	int nFiles;
	Long64_t nEvents = 0;
	for(int i=0;i<nSamples;i++){
		nFiles = _files.at(i).size();
		cout << "Files in sample: " <<  nFiles << endl;

		std::vector<TChain*>tempChain;
		for(int j=0;j<nFiles;j++){
			cout << "Loading file: " << endl;
			cout << _files.at(i).at(j) << endl;
			tempChain.push_back(new TChain(DrellYanVariables::treeName));
			tempChain.at(j)->Add(_base_directory+_files.at(i).at(j));
			nEvents += tempChain.at(j)->GetEntries();
			cout << "Loaded " << tempChain.at(j)->GetEntries() << endl;
		}
		_trees.push_back(tempChain);
	}
	_nEvents = nEvents;

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

	for(int i=0;i<_nSampleTypes;i++){
		nFiles = _files.at(i).size();
		bool isMC = true;
		for(int j=0;j<nFiles;j++){
			TBranch*testBranch = (TBranch*)_trees.at(i).at(j)->
				GetListOfBranches()->FindObject("GENEvt_weight");
			if(!testBranch) isMC = false; 
			InitializeBranches(_trees.at(i).at(j),isMC);
		}
	}
	
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
	}
	if(lepType!=ELE){
	//Muons
		chain->SetBranchAddress("nMuon",&nMuon,&b_nMuon);
		chain->SetBranchAddress("Nmuons",&Nmuons,&b_Nmuons);
		chain->SetBranchAddress("Muon_pT",&Muon_pT,&b_Muon_pT);
		chain->SetBranchAddress("Muon_eta",&Muon_eta,&b_Muon_eta);
		chain->SetBranchAddress("Muon_phi",&Muon_phi,&b_Muon_phi);
		chain->SetBranchAddress("Muon_passTightID",&Muon_passTightID,
					&b_Muon_passTightID);
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
	if(_ntupType==V2P6) xSecWeight = dataLuminosity*xSec.at(index2)/1.0;
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
	for(int i=0;i<_nSampleTypes;i++){
		int nFiles = _trees.at(i).size();
		Long64_t totalEvents = _nEvents;
		bool isMC = true;
		for(int j=0;j<nFiles;j++){//files within a sample (i.e. mass ranges)
			cout << "Processing file: " << dy_EE_v2p6.at(j) << endl;
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

			for(Long64_t iEntry=0;iEntry<nentries;iEntry++){
				eventCount++;

				//hard process variables
				double hardPt1 =  -1000;
				double hardPt2 =  -1000;
				double hardEta1 = -1000;
				double hardEta2 = -1000;
				double hardPhi1 = -1000;
				double hardPhi2 = -1000; 
				double invMassHard =  -1000;
				double rapidityHard = -1000;
				double ptHard =       -1000;
				
				//FSR variables
				double fsrPt1 =  -1000;
				double fsrPt2 =  -1000;
				double fsrEta1 = -1000;
				double fsrEta2 = -1000;
				double fsrPhi1 = -1000;
				double fsrPhi2 = -1000; 
				double invMassFSR =   -1000;
				double rapidityFSR =  -1000;
				double ptFSR =        -1000;

				//Reco variables
				double recoPt1 =  -1000;
				double recoPt2 =  -1000;
				double recoEta1 = -1000;
				double recoEta2 = -1000;
				double recoPhi1 = -1000;
				double recoPhi2 = -1000; 
				double invMassReco =  -1000;
				double rapidityReco = -1000;
				double ptReco =       -1000;

				//Pruned variables
				double prunedPt1 =  -1000;
				double prunedPt2 =  -1000;
				double prunedEta1 = -1000;
				double prunedEta2 = -1000;
				double prunedPhi1 = -1000;
				double prunedPhi2 = -1000; 
				double invMassPruned =  -1000;
				double rapidityPruned = -1000;
				double ptPruned =       -1000;

				double var;	
				_trees.at(i).at(j)->GetEntry(iEntry);
				double puWeight = GetPUWeight();
				int nDileptonsGen = GetGenLeptons(iHard1,iHard2,
					  	          	  iFSR1,iFSR2);
				int nDileptonsReco = GetRecoLeptons(iLep1,iLep2);

				if(nDileptonsGen>0 && nDileptonsReco>0){
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
					for(int k=0;k<_nHists;k++){
						_hists.at(k)->Fill(var.at(l),weights);
						l++;
						if(l>varSize) l=0;
					}
				}//end if nDileptons>0
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

	for(int iMu=0;iMu<Nmuons;iMu++){
		if(!Muon_passTightID[iMu]) continue;
		if(!PassMuonIsolation(iMu)) continue;
		for(int jMu=iMu+1;jMu<Nmuons;jMu++){
			pt1 = Muon_pT[iMu];
			pt2 = Muon_pT[jMu];
			eta1 = Muon_eta[iMu];
			eta2 = Muon_eta[jMu];
			phi1 = Muon_phi[iMu];
			phi2 = Muon_phi[jMu];

			if(!Muon_passTightID[jMu]) continue;
			if(!PassMuonIsolation(jMu)) continue;
			if(!PassMuonAngle(pt1,eta1,phi1,muMass,pt2,eta2,phi2,muMass)) 
				continue;

			//Other cuts not yet applied:
			//2 muons with opposite charge
			//smallest dimuon vertex chi2
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
}

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
}

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
}

bool DrellYanAnalyzer::PassGenToRecoMatchMu(int genIndex,int &recoIndex)
{
	using namespace DrellYanVariables;
	double dR,deta,dphi;
	float dRMin = 100000;
	recoIndex=-1;
	for(int iLep=0;iLep<Nmuons;iLep++){
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
	TString filesave = "data/DYHists";
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
