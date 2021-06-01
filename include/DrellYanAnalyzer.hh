#ifndef DrellYanAnalyzer_HH
#define DrellYanAnalyzer_HH

//ROOT functions
#include <TChain.h>
#include <TH1F.h>

class DrellYanAnalyzer
{
	public:
		DrellYanAnzlyer();
		int LoadData();
	private:
		//Tree definitions
		TChain*_treeDY_MC;
		TChain*_treeTops_MC;
		TChain*_treeEW_MC;
		TChain*_treeData;

		//Functions
		int LoadTrees();
		int InitializeBranches();
		double GetWeights();
		double CalcInvariantMass();
		bool PassAcceptance();
		bool PassHLT();
		bool PassGenToRecoMatch();
};//end class DrellYanAnalyzer

#endif


