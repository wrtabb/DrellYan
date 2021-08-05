#include "include/DrellYanAnalyzer.hh"

using namespace DrellYanVariables;

void test(LepType lepType,SampleType sampleType,NtupleType ntupType)
{
	DrellYanAnalyzer*dy = new DrellYanAnalyzer(ntupType,sampleType,lepType);
	dy->LoadData();
	Long64_t totalEntries = dy->GetNevents();
	cout << "Total entries loaded: " << totalEntries << endl;
	dy->EventLoop();

	delete dy;
}
