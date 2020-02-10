/*
Copyright (c) 2020 Richard King

The stressRefine analysis executable "SRwithMkl" is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

"SRwithMkl" is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

The terms of the GNU General Public License are explained in the file COPYING.txt,
also available at <https://www.gnu.org/licenses/>

Note that in it's current form "SRwithMkl" make use of the pardiso sparse solver
in the Intel MKL library, with which it must be linked.
Copyright (c) 2018 Intel Corporation.
You may use and redistribute the Intel MKL library, without modification, provided the conditions of
the Intel Simplified Software license are met:
https://software.intel.com/en-us/license/intel-simplified-software-license

It is perfectly permissable to replace the use of the pardiso software from the MKL library
with an equivalent open-source solver
*/

//////////////////////////////////////////////////////////////////////
//
// SRpostProcess.h: interface for the SRpostProcess class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(SRPOSTPROCESS_INCLUDED)
#define SRPOSTPROCESS_INCLUDED

#define FRDPRINT analysis.cgxFrdFile.PrintLine
#define FRDPRINTNORET analysis.cgxFrdFile.Print
#define FRDPRINTRET analysis.cgxFrdFile.PrintReturn()

class SRfile;
class SRvec3;

enum SRentityType;


enum SRstressComponent { xxComponent, yyComponent, zzComponent, xyComponent, xzComponent, yzComponent };

class SRpostStressVector
{
public:
	SRpostStressVector(){ for (int i = 0; i < 6; i++) stress[i] = 0.0; };
	void PlusAssign(double stresst[])
	{
		for (int i = 0; i < 6; i++)
			stress[i] += stresst[i];
	};
	void Assign(double stresst[])
	{
		for (int i = 0; i < 6; i++)
			stress[i] = stresst[i];
	};
	void Scale(double s)
	{
		for (int i = 0; i < 6; i++)
			stress[i] *= s;
	};
	double* Getstress(){ return stress; };
private:
	double stress[6];
};

class SRpostProcess
{
	friend class SRinput;
public:
	void PostProcess();
	double GlobalStrainSmooth();
	double ElementStrainSmoothByMaterial();
	void CalculateMaxStress();
	double StressMaxCheck(SRelement* elem, SRvec3& pos, double stress[]);
	void StressMaxCheckvsAllowable(SRelement* elem, double svmMax);
	void CalculateElementStresses(SRelement* elem, bool checkMaxOnly = false);
	void fillSacricialElementNodalStress(SRelement* elem, bool doMaxClipping);
	void PostProcessElementStresses();
	void Cleanup();
	double *getSmoothedStrainVec(int i){ return smoothedStrains[i].GetVector(); };
	int getElSmooth(int i){ return elSmooth.Get(i); };
	bool findSacrificialElementsAtKinks();
	int MeshToFrd();
	double getNodalStress(int i, int j){ return nodalStress.Get(i, j); };
	SRvec3& getNodeDisp(int i) { return nodeDisps.Get(i); };
	void freeNodalStress(){ nodalStress.Free(); };
	void allocateNodeDisps(int n){ nodeDisps.Allocate(n); };
	void putNodeDisp(int i, SRvec3& v){ nodeDisps.Put(i, v); };

private:
	SRdoubleVector smoothedStrains[6];
	SRintVector elemOnTopoList;
	SRintVector elSmooth;
	SRintVector nodalStressCount;
	SRdoubleMatrix nodalStress;
	bool outputCgx;
	SRintVector nodeElementOwnerSave;
	SRvector <SRvec3> nodeDisps;
	double svmAvg;
	SRintVector nodalPOrders;

public:

	//f06 output:
	void OutputF06();
	void OutputF06UseElementGetStress();
	void pageCheckF06(int& nlines, bool disp = false);
	void writeHeaderF06();
	void writeDispHeaderF06();
};

#endif // !defined(SRPOSTPROCESS_INCLUDED)
