#include "DrellYanAnalyzer.hh"

using namespace DrellYanVariables;

void test(NtupleType ntupType,LepType lepType,FileName fileName)
{
//	NtupleType ntupType = V2P6;
//	LepType lepType = ELE;
//	FileName fileName = DYLL_M700to800;
	DrellYanAnalyzer*dy = new DrellYanAnalyzer(ntupType,lepType,fileName);
	dy->LoadData();
	dy->EventLoop();

	delete dy;
}
