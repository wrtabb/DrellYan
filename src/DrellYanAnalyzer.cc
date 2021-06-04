#include "../include/DrellYanAnalyzer.hh"

DrellYanAnalyzer::DrellYanAnalyzer(DrellYanVariables::NtupleType ntupType)
{
	if(ntupType==DrellYanVariables::V2P6){
		_base_directory = DrellYanVariables::base_directory_v2p6;
		_files_EE = DrellYanVariables::dy_EE_v2p6;
		_files_MuMu = DrellYanVariables::dy_MuMu_v2p6;
		_files_tops = DrellYanVariables::tops_v2p6;
		_files_fakes = DrellYanVariables::fakes_v2p6;
		_files_dibosons = DrellYanVariables::dibosons_v2p6;
		_files_taus = DrellYanVariables::taus_v2p6;
		_files_data = DrellYanVariables::data_v2p6;

		_files.push_back(_files_EE);
		_files.push_back(_files_MuMu);
		_files.push_back(_files_tops);
		_files.push_back(_files_fakes);
		_files.push_back(_files_dibosons);
		_files.push_back(_files_taus);
		_files.push_back(_files_data);
	}
}

//-----Load data-----//
int DrellYanAnalyzer::LoadData()
{
	int returnCode = LoadTrees();
	return returnCode;	
}

int DrellYanAnalyzer::LoadTrees()
{
	int returnCode = 1;
	TTimeStamp ts_start;
	cout << "Begin loading trees:" << endl;
	cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;
	TStopwatch totaltime;
	totaltime.Start();

	int nFiles;
	for(int i=0;i<_nSampleTypes;i++){
		nFiles = _files.at(i).size();
		cout << "Files in sample: " <<  nFiles << endl;

		std::vector<TChain*>tempChain;
		for(int j=0;j<nFiles;j++){
			cout << "Loading file: " << endl;
			cout << _files.at(i).at(j) << endl;
			tempChain.push_back(new TChain(DrellYanVariables::treeName));
			tempChain.at(j)->Add(_base_directory+_files.at(i).at(j));
			cout << "Loaded " << tempChain.at(j)->GetEntries() << endl;
		}
		_trees.push_back(tempChain);
	}

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
			if(i==6) isMC=false;
			InitializeBranches(_trees.at(i).at(j),isMC);
		}
	}
	
	return returnCode;
}

int DrellYanAnalyzer::InitializeBranches(TChain*chain,bool isMC)
{
	using namespace DrellYanVariables;

	//Electrons
	chain->SetBranchAddress("Nelectrons",&Nelectrons,&b_Nelectrons);
	chain->SetBranchAddress("Electron_pT",&Electron_pT,&b_Electron_pT);
	chain->SetBranchAddress("Electron_eta",&Electron_eta,&b_Electron_eta);
	chain->SetBranchAddress("Electron_phi",&Electron_phi,&b_Electron_phi);
	chain->SetBranchAddress("Electron_passMediumID",&Electron_passMediumID,
				 &b_Electron_passMediumID);
	//Muons
	chain->SetBranchAddress("nMuon",&nMuon,&b_nMuon);
	chain->SetBranchAddress("Nmuons",&Nmuons,&b_Nmuons);
	chain->SetBranchAddress("Muon_pT",&Muon_pT,&b_Muon_pT);
	chain->SetBranchAddress("Muon_eta",&Muon_eta,&b_Muon_eta);
	chain->SetBranchAddress("Muon_phi",&Muon_phi,&b_Muon_phi);
	chain->SetBranchAddress("Muon_passTightID",&Muon_passTightID,&b_Muon_passTightID);
	
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
}

//-----Calculations
double DrellYanAnalyzer::GetWeights()
{
	return 0.0;
}

double DrellYanAnalyzer::CalcInvariantMass()
{
	return 0.0;
}

int DrellYanAnalyzer::GetGenLeptons(DrellYanVariables::LepType lepType,int &idxHardLep1,
				    int &idxHardLep2,int &idxFSRLep1,int &idxFSRLep2)
{
	using namespace DrellYanVariables;
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
				idxHardLep1 = iLep;
				idxHardLep2 = jLep;
				nDileptons++;
			}//end if hard process
			if(GENLepton_fromHardProcessFinalState[iLep]==1 &&
			GENLepton_fromHardProcessFinalState[jLep]==1){
				idxFSRLep1 = iLep;
				idxFSRLep2 = jLep;
			}//end if FSR
		}//end jLep loop
	}//end iLep loop
	return nDileptons;
}//end GetGenLeptons()

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

	for(int iMu = 0; iMu < Nelectrons; iMu++) {
		if(!Muon_passTightID[iMu]) continue;
		for(int jMu = iMu+1; jMu < Nelectrons; jMu++) {
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

//-----Cut criteria-----//
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

bool DrellYanAnalyzer::PassHLT(DrellYanVariables::LepType lepType)
{
	using namespace DrellYanVariables;

	TString trigName;
	TString triggerUsed;
	if(lepType==ELE) triggerUsed = electronTrigger;
	else if(lepType==MUON) triggerUsed = muonTrigger;
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
		if(trigName.CompareTo(triggerUsed)==0 && HLT_trigFired[iHLT]==1){
			passHLT = true;
			break;
		}//end if trigName...
	}//end loop over triggers
	return passHLT;
}

bool DrellYanAnalyzer::PassGenToRecoMatch()
{
	return false;
}

