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
// SRerrorCheck.h: interface for the SRerrorCheck class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(SRERRORCHECK_INCLUDED)
#define SRERRORCHECK_INCLUDED

class SRelement;
class SRface;
class SRconstraint;

class SRerrorCheck
{
public:
	SRerrorCheck();
	void Initialize();
	void CleanUp();
	double fillFaceTractionJumps();
	double FindFaceTractionJumpError(SRelement* elem);
	double FindSmoothedVsRawError(SRelement* elem);
	bool FindRequiredPOrder(SRelement* elem, int& p);
	double FindError(SRelement* elem);
	double getElementFaceTraction(int lface, SRface *face, SRelement* elem, double rf, double sf, double& r, double& s, double& t, double* tract);
	double getStressMax(){ return stressMax; };
	void setStressMax(double s){ stressMax = s; };
	double getStrainMax(){ return strainMax; };
	void setStrainMax(double e){ strainMax = e; };
	bool getSmallMaxStressDetected(){ return smallMaxStressDetected; };
	void setSmallMaxStressDetected(bool tf){ smallMaxStressDetected = tf; };
	double GetStressMaxForErrorCheck(SRelement* elem);
	double GetStrainMaxForErrorCheck(SRelement* elem);
	void SetLowStressTol(double tol){ lowStressTolFinalAdapt = tol; };
	double GetLowStressTol(){ return lowStressTol; };
	int GetnumSacrificialElements(){ return numSacrificialElements; };
	bool AutoSacrificialElements();
	//support routines for AutoSacrificialElements:
	bool checkForSymCon(SRface* face, SRconstraint* con, int connedDof);
	bool checkFacesNormal(SRface *face, int lej, SRface* face2, int lej2);
	bool checkForTangentialCon(SRface* face, SRconstraint* con, SRface* face2, SRconstraint* con2, int connedDof, int connedDof2);
	bool checkForOrthogonalCon(SRface* face, SRconstraint* con, int connedDof);
	double getElSvmAtMidNode(SRelement* elem, int mid);

	void setLowStressTolFinalAdapt(double tol){ lowStressTolFinalAdapt = tol; };
	void setCheckReentrant(bool tf){ checkReentrant = tf; };
	void SetUp(bool finalAdapt);

private:
	double strainMax;
	double stressMax;
	bool smallMaxStressDetected;
	double lowStressTol;
	double lowStressTolFinalAdapt;
	bool checkReentrant;
	int numSacrificialElements;
};

#endif // !defined(SRERRORCHECK_INCLUDED)
