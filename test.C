#include "DrellYanAnalyzer.hh"

using namespace DrellYanVariables;

void test(NtupleType ntupType,SampleType sampleType,LepType lepType,TString fileName)
{
	DrellYanAnalyzer*dy = new DrellYanAnalyzer(ntupType,sampleType,lepType,fileName);
	dy->LoadData();
	Long64_t totalEntries = dy->GetNevents();
	cout << "Total entries loaded: " << totalEntries << endl;
	dy->EventLoop();

	delete dy;
}
