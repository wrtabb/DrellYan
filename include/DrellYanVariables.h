
namespace DrellYanVariables
{
	enum NtupleTypes{
		V2P3,
		V2P6,
		V2P7,
		V2P8
	};

	TString base_directory_v2p6 = "root://xrootd-local.unl.edu///store/user/wtabb/DrellYan_13TeV_2016/v2p6/skims/";
	vector<TString> file_locations_v2p6 = {
		//Drell-Yan MC signal
		//Electrons
		"DYLL_M10to50_EE_ext1v1.root",				//0
		"DYLL_M10to50_EE_v1.root",				//1
		"DYLL_M10to50_EE_v2.root",				//2
		"DYLL_M50toInf_truncated_M50To100_EE_base.root",	//3
		"DYLL_M100to200_EE.root",				//4
		"DYLL_M200to400_EE.root",				//5
		"DYLL_M400to500_EE.root",				//6
		"DYLL_M500to700_EE.root",				//7
		"DYLL_M700to800_EE.root",				//8
		"DYLL_M800to1000_EE.root",				//9
		"DYLL_M1000to1500_EE.root",				//10
		"DYLL_M1500to2000_EE.root",				//11
		"DYLL_M2000to3000_EE.root",				//12

		//Muons
		"DYLL_M10to50_MuMu_ext1v1.root",			//13
		"DYLL_M10to50_MuMu_v1.root",				//14
		"DYLL_M10to50_MuMu_v2.root",				//15
		"DYLL_M50toInf_truncated_M50To100_MuMu_base.root",	//16
		"DYLL_M100to200_MuMu.root",				//17
		"DYLL_M200to400_MuMu.root",				//18
		"DYLL_M400to500_MuMu.root",				//19
		"DYLL_M500to700_MuMu.root",				//20
		"DYLL_M700to800_MuMu.root",				//21
		"DYLL_M800to1000_MuMu.root",				//22
		"DYLL_M1000to1500_MuMu.root",				//23
		"DYLL_M1500to2000_MuMu.root",				//24
		"DYLL_M2000to3000_MuMu.root",				//25
		
		//Backgrounds
		//tops
		"ST_tbarW.root",					//26
		"ST_tW.root",						//27
		"ttbarBackup.root",					//28
		"ttbarBackup_truncated_M0To700.root",			//29
		"ttbar_M1000toInf.root",				//30
		"ttbar_M700to1000.root",				//31
		"ttbar.root",						//32
		"ttbar_truncated_M0To700.root",				//33
		//W+Jets (Fakes)		
		"WJetsToLNu_amcatnlo_ext2v5.root",			//34
		"WJetsToLNu_amcatnlo_ext.root",				//35
		"WJetsToLNu_amcatnlo.root",				//36
		//EW (bosons + DY->TauTau)
		//Bosons
		"WW.root",						//37
		"WZ.root",						//38
		"ZZ.root",						//39
		//Taus
		"DYLL_M10to50_TauTau_ext1v1.root",			//40
		"DYLL_M10to50_TauTau_v1.root",				//41
		"DYLL_M10to50_TauTau_v2.root",				//42
		"DYLL_M50toInf_truncated_M50To100_TauTau_base.root",	//43
		"DYLL_M100to200_TauTau.root",				//44
		"DYLL_M200to400_TauTau.root",				//45
		"DYLL_M400to500_TauTau.root",				//46
		"DYLL_M500to700_TauTau.root",				//47
		"DYLL_M700to800_TauTau.root",				//48
		"DYLL_M800to1000_TauTau.root",				//49
		"DYLL_M1000to1500_TauTau.root",				//50
		"DYLL_M1500to2000_TauTau.root",				//51
		"DYLL_M2000to3000_TauTau.root",				//52

		//Drell-Yan Data
		"DoubleEG_RunB.root",					//53
		"DoubleEG_RunC.root",					//54
		"DoubleEG_RunD.root",					//55
		"DoubleEG_RunE.root",					//56
		"DoubleEG_RunF.root",					//57
		"DoubleEG_RunG.root",					//58
		"DoubleEG_RunHver2.root",				//59
		"DoubleEG_RunHver3.root",				//60
	};
}
