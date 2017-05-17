#ifndef GetUnitHistogram_h
#define GetUnitHistogram_h 1

#include "TH1F.h"
class RootHistograms{
	
	public:
	
	RootHistograms();
	~RootHistograms();
//	TH1F * UnitResponseShape;
	TH1F* MakeUnitHistogram();
	
	private:
	TH1F *UnitShapeCh1;
	TH1F *BaseLineShapeCh1Raw;
	

	
	};


#endif
