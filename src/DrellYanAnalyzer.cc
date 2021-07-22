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
		_files_LL = dy_EE_test;
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
	SampleType sampleType = _sampleType;
	LepType lepType = _lepType;
	vector<TString> histTypes = _histTypes;
	int nHistTypes = _nHistTypes;
	vector<TString> sampleNames;

	if(sampleType==SAMPLE_LL && lepType==ELE) sampleNames.push_back("_EE");
	if(sampleType==SAMPLE_LL && lepType==MUON) sampleNames.push_back("_MuMu");
	else if(sampleType==SAMPLE_TOP) sampleNames.push_back("_Tops");
	else if(sampleType==SAMPLE_FAKE) sampleNames.push_back("_Fakes");
	else if(sampleType==SAMPLE_DIBOSON) sampleNames.push_back("Dibosons");
	else if(sampleType==SAMPLE_TAU)	sampleNames.push_back("_TauTau");
	else if(sampleType==SAMPLE_DATA) sampleNames.push_back("_EGData");
	else if(sampleType==SAMPLE_ALL){
		if(lepType==ELE) sampleNames.push_back("_EE");
		else if(lepType==MUON) sampleNames.push_back("_MuMu");
		sampleNames.push_back("_Tops");
		sampleNames.push_back("_Fakes");
		sampleNames.push_back("Dibosons");
		sampleNames.push_back("_TauTau");
		sampleNames.push_back("_EGData");
	}//end if SAMPLE_ALL

	int nSamples = _nSampleTypes;
	TString totalName;
	vector<TH1D*> _histVars;
	for(int i=0;i<nSamples;i++){
		TString massName = "histInvMass_";
		massName += sampleNames.at(i);	
		TH1D*hMass = new TH1D(massName,"",nMassBins,massbins);
		TString rapidityName = "histRapidity_";
		rapidityName += sampleNames.at(i);	
		TH1D*hRapidity = new TH1D(rapidityName,"",nRapidityBins,rapiditybins);
		TString ptName = "histPt_";
		ptName += sampleNames.at(i);	
		TH1D*hPt = new TH1D(ptName,"",nPtBins,ptbins);

		_histVars.push_back(hMass);
		_histVars.push_back(hRapidity);
		_histVars.push_back(hPt);
	
		_hists.push_back(_histVars);
	}//end loop over sample types
	
}//end InitializeHistograms

/////////////////////////
//-----Get Weights-----//
/////////////////////////
double DrellYanAnalyzer::GetWeights(int index1,int index2)
{
	using namespace DrellYanVariables;
	SampleType sampleType = _sampleType;
	//-----Cross section-----//
	vector<double> xSec;
	if(sampleType==SAMPLE_LL) xSec = xSec_LL;
	double xSecWeight = dataLuminosity*xSec.at(index2)/1.0;

	double sumGenWeight = GetGenWeightSum(index1,index2);
	double genWeight = (GENEvt_weight/fabs(GENEvt_weight))/sumGenWeight;  
	double weight = xSecWeight*genWeight;
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


///////////////////////////////////////
//-----Get Leptons From an Event-----//
///////////////////////////////////////
int DrellYanAnalyzer::EventLoop()
{
	using namespace DrellYanVariables;
	int nHists = _nHistTypes;
	LepType lepType = _lepType;
	InitializeHistograms();
	int iHard1,iHard2;
	int iFSR1,iFSR2;
	int iDressed1,iDressed2;
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
		for(int j=0;j<nFiles;j++){//files within a sample (i.e. mass ranges)
			Long64_t nentries = _trees.at(i).at(j)->GetEntries();
			double weights = GetWeights(i,j);
			for(Long64_t iEntry=0;iEntry<nentries;iEntry++){
				eventCount++;
				double pt1,pt2;
				double eta1,eta2;
				double phi1,phi2; 
				double var;	
				_trees.at(i).at(j)->GetEntry(iEntry);
				int nDileptons = DrellYanAnalyzer::GetGenLeptons(iHard1,iHard2,
					  	  		                 iFSR1,iFSR2);

				if(nDileptons>0){
					pt1  = GENLepton_pT[iHard1];
					eta1 = GENLepton_eta[iHard1];
					phi1 = GENLepton_phi[iHard1];
					pt2  = GENLepton_pT[iHard2];
					eta2 = GENLepton_eta[iHard2];
					phi2 = GENLepton_phi[iHard2];
					double invMass = CalcVariable(pt1,eta1,phi1,lMass,
								      pt2,eta2,phi2,lMass,
								      INV_MASS);
					double rapidity = CalcVariable(pt1,eta1,phi1,lMass,
								       pt2,eta2,phi2,lMass,
								       RAPIDITY);
					double pt = CalcVariable(pt1,eta1,phi1,lMass,
								 pt2,eta2,phi2,lMass,
								 PT);
					for(int k=0;k<nHists;k++){
						if(k==0) var = invMass;
						else if(k==1) var = rapidity;
						else if(k==2) var = pt;
						_hists.at(i).at(k)->Fill(var,weights);
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

	for(int iMu=0;iMu<Nmuons;iMu++){
		if(!Muon_passTightID[iMu]) continue;
		for(int jMu=iMu+1;jMu<Nmuons;jMu++){
			if(!Muon_passTightID[jMu]) continue;
			//Other cuts not yet applied:
			//PFIso/pt<0.15
			//2 muons with opposite charge
			//Angle between muons < pi - 0.005 rad
			//smallest dimuon vertex chi2
			pt1 = Muon_pT[iMu];
			pt2 = Muon_pT[jMu];
			eta1 = Muon_eta[iMu];
			eta2 = Muon_eta[jMu];
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
	if(lepType==ELE) triggerUsed = electronTrigger;
	else if(lepType==MUON) triggerUsed = muonTrigger1;
	else {
		cout << "Error in DrellYanAnalyzer::PassHLT()" << endl;
		cout << "lepType does not exist" << endl;
		cout << "Must choose 'ELE' or 'MUON'" << endl;
		return false;
	}
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
			    trigName.CompareTo(muonTrigger1)==0 ) && 
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
	double pi = 3.1415926;
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
	double pi = 3.1415926;
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
}

///////////////////////
//-----Utilities-----//
///////////////////////
void DrellYanAnalyzer::SaveResults()
{
	using namespace DrellYanVariables;
	TFile*file = new TFile("data/DrellYanHistograms.root","recreate");
	for(int i=0;i<_nSampleTypes;i++){
		int nFiles = _files.at(i).size();
		for(int j=0;j<_nHistTypes;j++){
			_hists.at(i).at(j)->Write();
		}
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

	TLorentzVector vGenLepton1;
	TLorentzVector vGenLepton2;
	vGenLepton1.SetPtEtaPhiM(pt1,eta1,phi1,mass1);
	vGenLepton2.SetPtEtaPhiM(pt2,eta2,phi2,mass2);
	double invMass = (vGenLepton1+vGenLepton2).M();
	double rapidity= (vGenLepton1+vGenLepton2).Rapidity();
	double pt = (vGenLepton1+vGenLepton2).Pt();

	if(varType==INV_MASS) return invMass;
	else if(varType==RAPIDITY) return rapidity;
	else if(varType==PT) return pt;
	else return -100000000;
}

void DrellYanAnalyzer::Counter(Long64_t i,Long64_t total){
	int P = 100*(i)/(total);  
	TTimeStamp eventTimeStamp;
	if(i%(total/100)==0)
	cout << "[Time: " << eventTimeStamp.AsString("s") << "] " << P << "%" << endl;
}
