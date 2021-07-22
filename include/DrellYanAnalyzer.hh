#ifndef DrellYanAnalyzer_HH
#define DrellYanAnalyzer_HH

//ROOT functions
#include <TChain.h>
#include <TH1F.h>
#include <TString.h>
#include <iostream>
#include <TTimeStamp.h>
#include <TStopwatch.h>
#include <TLorentzVector.h>
#include <TFile.h>

//Global Variables for the DrellYanVariables namespace
#include "DrellYanVariables.h"

class DrellYanAnalyzer
{
	public:
		DrellYanAnalyzer(DrellYanVariables::NtupleType ntupType,
				 DrellYanVariables::SampleType sampletype,
				 DrellYanVariables::LepType lepType);
		int GetRecoLeptons(int &leadLep,int &subLep);
		int GetGenLeptons(int &iHard1,int &iHard2,
				  int &iFSR1,int &iFSR2);
		int LoadData();
		Long64_t GetNevents();
		int EventLoop();
	private:
		Long64_t _nEvents;
		int _nSampleTypes;
		DrellYanVariables::SampleType _sampleType;
		DrellYanVariables::LepType _lepType;

		//Histogram definitions
		vector<vector<TH1D*>> _hists;
	
		//Tree definitions
		vector<TChain*>_treeLL;
		vector<TChain*>_treeTops;
		vector<TChain*>_treeFakes;
		vector<TChain*>_treeDibosons;
		vector<TChain*>_treeTaus;
		vector<TChain*>_treeData;
		vector<vector<TChain*>> _trees;
	
		//Files
		TString _base_directory;
		vector<TString> _files_LL;//# of files: 11
		vector<TString> _files_tops;//# of files: 5
		vector<TString> _files_fakes;//# of files: 3
		vector<TString> _files_dibosons;//# of files: 3
		vector<TString> _files_taus;//# of files: 11
		vector<TString> _files_data;//# of files: 8
		vector<vector<TString>> _files;

		//-----Functions-----//
		//loading data
		int LoadTrees();
		int InitializeBranches(TChain*chain,bool isMC);
		void InitializeHistograms();

		//calculations
		double GetWeights(int index,int index2);
		double GetGenWeightSum(int index,int index2);
		double GetPUWeight();
		double CalcVariable(double pt1,double eta1,double phi1,double mass1,
				    double pt2,double eta2,double phi2,double mass2,
				    DrellYanVariables::VarType varType);

		//cuts
		bool PassAcceptance(double pt1,double pt2,double eta1,double eta2);
		bool PassHLT();
		bool PassGenToRecoMatch(int genIndex,int &recoIndex);
		bool PassGenToRecoMatchEle(int genIndex,int &recoIndex);
		bool PassGenToRecoMatchMu(int genIndex,int &recoIndex);

		//Get dilepton pairs
		int GetRecoElectrons(int &leadEle,int &subEle);
		int GetRecoMuons(int &leadMu,int &subMu);

		void SaveResults();
		void Counter(Long64_t i,Long64_t total);
};//end class DrellYanAnalyzer

#endif


