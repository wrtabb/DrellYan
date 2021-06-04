#ifndef DrellYanAnalyzer_HH
#define DrellYanAnalyzer_HH

//ROOT functions
#include <TChain.h>
#include <TH1F.h>
#include <TString.h>
#include <iostream>
#include <TTimeStamp.h>
#include <TStopwatch.h>

#include "DrellYanVariables.h"

class DrellYanAnalyzer
{
	public:
		DrellYanAnalyzer(DrellYanVariables::NtupleType ntupType);
		int LoadData();
	private:
		static const int _nFilesEE = 11; 
		static const int _nFilesMuMu = 11;
		static const int _nFilesTops = 5;
		static const int _nFilesFakes = 3;
		static const int _nFilesDibosons = 3;
		static const int _nFilesTaus = 11;
		static const int _nFilesData = 8;
		static const int _nSampleTypes = 7;

		//Tree definitions
		TChain*_treeDY_EE[_nFilesEE];
		TChain*_treeDY_MuMu[_nFilesMuMu];
		TChain*_treeTops[_nFilesTops];
		TChain*_treeFakes[_nFilesFakes];
		TChain*_treeDibosons[_nFilesDibosons];
		TChain*_treeTaus[_nFilesTaus];
		TChain*_treeData[_nFilesData];
		vector<vector<TChain*>> _trees;
	
		//Variables
		TString _base_directory;
		vector<TString> _files_EE;
		vector<TString> _files_MuMu;
		vector<TString> _files_tops;
		vector<TString> _files_fakes;
		vector<TString> _files_dibosons;
		vector<TString> _files_taus;
		vector<TString> _files_data;
		vector<vector<TString>> _files;

		//-----Functions-----//
		//loading data
		int LoadTrees();
		int InitializeBranches(TChain*chain,bool isMC);

		//calculations
		double GetWeights();
		double CalcInvariantMass();

		//cuts
		bool PassAcceptance(double pt1,double pt2,double eta1,double eta2);
		bool PassHLT(DrellYanVariables::LepType lepType);
		bool PassGenToRecoMatch();

		//Get dilepton pairs
		int GetGenLeptons(DrellYanVariables::LepType lepType,int &idxHardLep1,
				  int &idxHardLep2,int &idxFSRLep1,int &idxFSRLep2);
		int GetRecoElectrons(int &leadEle,int &subEle);
		int GetRecoMuons(int &leadMu,int &subMu);
};//end class DrellYanAnalyzer

#endif


