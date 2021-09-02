#include "DrellYanAnalyzer.hh"

using namespace DrellYanVariables;

void test(NtupleType ntupType,LepType lepType,FileName fileName)
{
	DrellYanAnalyzer*dy = new DrellYanAnalyzer(ntupType,lepType,fileName);
	dy->LoadData();
	dy->EventLoop();

	delete dy;
}
