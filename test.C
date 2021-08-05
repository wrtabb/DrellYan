#include "include/DrellYanAnalyzer.hh"

using namespace DrellYanVariables;

void RunLoop(NtupleType ntupType,SampleType sampleType,LepType lepType);

void test()
{
	//RunLoop(TEST,SAMPLE_LL,ELE);
	RunLoop(TEST,SAMPLE_LL,MUON);
}

void RunLoop(NtupleType ntupType,SampleType sampleType,LepType lepType)
{
	DrellYanAnalyzer*dy = new DrellYanAnalyzer(ntupType,sampleType,lepType);
	dy->LoadData();
	Long64_t totalEntries = dy->GetNevents();
	cout << "Total entries loaded: " << totalEntries << endl;
	dy->EventLoop();

	delete dy;
}
