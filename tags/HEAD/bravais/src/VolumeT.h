/* $Id: VolumeT.h,v 1.1.1.1 2002-02-28 02:13:08 jzimmer Exp $ */
#include <iostream>
#include "dArrayT.h"
#include "dArray2DT.h"

class ifstreamT;

class VolumeT {
protected:
	int nSD;
	double volume;
public:
        VolumeT(int n);
	~VolumeT();
	int GetDimensions();
	double GetVolume();
	void WriteFile();

	virtual void DefineBoundary(ifstreamT& in) = 0;
	virtual void CalculateVolume() = 0;
	virtual void FillVolume() = 0;
};

class BoxT : public VolumeT {
private:
        dArrayT length;
        dArray2DT surfaces;
public:
        BoxT(int n);
        ~BoxT();
        void DefineBoundary(ifstreamT& in);
        void CalculateVolume();
        void FillVolume();
};





