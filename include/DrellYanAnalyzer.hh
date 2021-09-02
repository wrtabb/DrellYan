#ifndef DrellYanAnalyzer_HH
#define DrellYanAnalyzer_HH

//ROOT functions
#include <TChain.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TString.h>
#include <iostream>
#include <TTimeStamp.h>
#include <TStopwatch.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <TMath.h>

//Global Variables for the DrellYanVariables namespace
#include "DrellYanVariables.h"

class DrellYanAnalyzer
{
	public:
		DrellYanAnalyzer(DrellYanVariables::NtupleType ntupType,
				 DrellYanVariables::LepType lepType,
				 DrellYanVariables::FileName fileName);
		int GetRecoLeptons(int &leadLep,int &subLep);
		int GetGenLeptons(int &iHard1,int &iHard2,
				  int &iFSR1,int &iFSR2);
		int LoadData();
		Long64_t GetNevents();
		int EventLoop();
		bool GetIsMC();
	private:
		Long64_t _nEvents;
		bool _isMC;

		DrellYanVariables::SampleType _sampleType;
		DrellYanVariables::LepType _lepType;
		DrellYanVariables::NtupleType _ntupType;
		DrellYanVariables::FileName _fileName;

		//Histogram definitions
		TH1D*_hMassHardProcess;
		TH1D*_hRapidityHardProcess;
		TH1D*_hPtHardProcess;
		TH1D*_hMassReco;
		TH1D*_hRapidityReco;
		TH1D*_hPtReco;
		TH1F*_hPileupRatio;
		
		TH2F*_hRecoSF;
		TH2F*_hMedIDSF;
		TH2F*_hLeg2SF;
		void DefineHistogramProperties(TH1*hist);
	
		//Tree definitions
		TChain*_tree;
	
		//Files
		TString _base_directory;
		TString _FileToLoad;

		TFile*_fRecoSF;
		TFile*_fMedIDSF;
		TFile*_fLeg2SF;

		//-----Functions-----//
		//loading data
		int LoadTrees();
		int InitializeBranches(TChain*chain);
		void InitializeHistograms();
		TString GetSampleName();
		DrellYanVariables::SampleType GetSampleType();

		//calculations
		double GetSampleWeights();
		double GetEventWeights(double pt1,double pt2,double eta1,double eta2);
		double GetGenWeightSum();
		double GetPUWeight();
		double GetEleScaleFactors(double pt1,double pt2,double eta1,double eta2);
		double GetInvMass(double pt1,double eta1,double phi1,
				  double pt2,double eta2,double phi2);
		double GetRapidity(double pt1,double eta1,double phi1,
				   double pt2,double eta2,double phi2);
		double GetPt(double pt1,double eta1,double phi1,
			     double pt2,double eta2,double phi2);
		double CalcVariable(double pt1,double eta1,double phi1,double mass1,
				    double pt2,double eta2,double phi2,double mass2,
				    DrellYanVariables::VarType varType);
		double GetVertexChi2(int index1,int index2);

		//cuts
		bool PassAcceptance(double pt1,double pt2,double eta1,double eta2);
		bool PassHLT();
		bool PassGenToRecoMatch(int genIndex,int &recoIndex);
		bool PassGenToRecoMatchEle(int genIndex,int &recoIndex);
		bool PassGenToRecoMatchMu(int genIndex,int &recoIndex);
		bool PassMuonIsolation(int index);
		bool PassMuonAngle(double pt1,double eta1,double phi1,double mass1,
				   double pt2,double eta2,double phi2,double mass2);

		//Get dilepton pairs
		int GetRecoElectrons(int &leadEle,int &subEle);
		int GetRecoMuons(int &leadMu,int &subMu);

		void SaveResults();
		void Counter(Long64_t i,Long64_t total);
};//end class DrellYanAnalyzer

#endif


