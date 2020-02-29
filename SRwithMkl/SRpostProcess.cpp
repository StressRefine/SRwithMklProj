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
// SRpostProcess.cpp: implementation of the SRpostProcess class.
//
//////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <search.h>
#include "SRmodel.h"
#include "SRanalysis.h"
#include "SRelement.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#endif

extern SRmodel model;
extern SRanalysis analysis;

void SRpostProcess::PostProcess()
{
	//postprocessing:
	//output mesh, displacement results, and stress results

	analysis.zeroStressMax();
	for (int i = 0; i < model.GetNumMaterials(); i++)
	{
		SRmaterial* mat = model.GetMaterial(i);
		mat->PutMaxSvm(0.0);
	}

	int i;
	SRnode* node;

	int numnode = model.GetNumNodes();

	//default output: Cgx (Calculix) .frd format

	if (!analysis.cgxFrdFile.OutOpenNoFail())
	{
		LOGPRINT("unable to open cgx frd file for postprocessing");
		LOGPRINT("this model may be open already in Cgx postprocessor");
		REPPRINT("unable to open cgx frd file for postprocessing");
		REPPRINT("this model may be open already in Cgx postprocessor");
		ERROREXIT;
	}


	//mesh definition:
	int numNodeOut = MeshToFrd();

	//stresses:
	PostProcessElementStresses();

	//stress at each node is now stored in nodalStress (node num, component num)
	FRDPRINT("  100CL  101%12.5le%12d                     0    1           1", 1.0, numNodeOut);
	FRDPRINT(" -4  STRESS      6    1");
	FRDPRINT(" -5  SXX         1    4    1    1");
	FRDPRINT(" -5  SYY         1    4    2    2");
	FRDPRINT(" -5  SZZ         1    4    3    3");
	FRDPRINT(" -5  SXY         1    4    1    2");
	FRDPRINT(" -5  SYZ         1    4    2    3");
	FRDPRINT(" -5  SZX         1    4    3    1");
	double stressConv = analysis.getstressUnitConversion();
	for (i = 0; i < numnode; i++)
	{
		node = model.GetNode(i);
		if (node->isOrphan())
			continue;
		FRDPRINTNORET(" -1%10d", node->GetUserid());
		FRDPRINTNORET("%12.5le", stressConv*nodalStress.Get(i, xxComponent));
		FRDPRINTNORET("%12.5le", stressConv*nodalStress.Get(i, yyComponent));
		FRDPRINTNORET("%12.5le", stressConv*nodalStress.Get(i, zzComponent));
		FRDPRINTNORET("%12.5le", stressConv*nodalStress.Get(i, xyComponent));
		FRDPRINTNORET("%12.5le", stressConv*nodalStress.Get(i, yzComponent));
		FRDPRINTNORET("%12.5le", stressConv*nodalStress.Get(i, xzComponent));
		FRDPRINTRET;
	}
	FRDPRINT(" -3");

	//displacements
	FRDPRINT("  100CL  101%12.5le%12d                     0    1           1", 1.0, numNodeOut);
	FRDPRINT(" -4  DISP        4    1");
	FRDPRINT(" -5  D1          1    2    1    0");
	FRDPRINT(" -5  D2          1    2    2    0");
	FRDPRINT(" -5  D3          1    2    2    0");
	FRDPRINT(" -5  ALL         1    2    0    0    1ALL");
	SRvec3 disp;
	for (i = 0; i < numnode; i++)
	{
		node = model.GetNode(i);
		if (node->isOrphan())
			continue;
		node->GetDisp(disp);
		disp.Scale(analysis.getlengthUnitConversion());
		FRDPRINT(" -1%10d%12.5le%12.5le%12.5le", node->GetUserid(), disp.d[0], disp.d[1], disp.d[2]);
	}
	FRDPRINT(" -3");

	FRDPRINT("9999"); //end of data

	analysis.cgxFrdFile.Close();

	if (analysis.getoutputf06())
		OutputF06();

	nodalStress.Free();
	nodeDisps.Free();
}


void SRpostProcess::PostProcessElementStresses()
{
	//output stress results for each element
	for (int i = 0; i < model.GetNumMaterials(); i++)
	{
		SRmaterial* mat = model.GetMaterial(i);
		mat->SetVolPercentYielded(0.0);
	}

	int e;
	SRelement* elem;

	int nn = model.GetNumNodes();
	nodalStress.Allocate(nn, 6);
	nodalStressCount.Allocate(nn);

	for (e = 0; e < model.GetNumElements(); e++)
	{
		elem = model.GetElement(e);
		CalculateElementStresses(elem);
	}

	int nel = model.GetNumElements();
	bool doMaxClipping = (nel - analysis.errorChecker.GetnumSacrificialElements() > 0);
	for (e = 0; e < nel; e++)
	{
		elem = model.GetElement(e);
		if (elem->isSacrificial())
			fillSacricialElementNodalStress(elem, doMaxClipping);
	}

	nodalStressCount.Free();
}

void SRpostProcess::CalculateElementStresses(SRelement* elem, bool checkMaxOnly)
{
	//calculate stress results for an element

	//input:
		//elem = pointer to element
		//checkMaxOnly = true for check max only, false to store stress in nodalStress 
	//note:
		//nodalStress stores stress at each node in model
		//if global stress smoothing has not been performed, stress contribution
		//from each element can be different so averaging is performed

	static int numNodesTotalBrick = 27; //nodes, centroid, face centroids
	static int numNodesTotalTet = 15; //nodes, centroid, face centroids
	static int numNodesTotalWedge = 21; //nodes, centroid, face centroids
	int i, nn, ne;
	double r, s, t;
	double stress[6];

	nn = elem->GetNumNodes();
	ne = elem->GetNumLocalEdges();
	SRnode* node;
	SRedge* edge;
	int id;
	double elSvmMax = 0.0;
	elem->SetSvmMax(0.0);
	double svm;
	elem->SetBasisData();
	bool hasAllowable = !checkMaxOnly && analysis.getallMatsHaveAllowable() && elem->GetMaterial()->isAllowableAssigned();
	double Allowable = 0.0;
	if(hasAllowable)
		Allowable = elem->GetMaterial()->GetAllowableStress();
	int numNodesYielded = 0;
	for (i = 0; i < nn; i++)
	{
		id = elem->GetNodeId(i);
		node = model.GetNode(id);
		elem->NodeNaturalCoords(i, r, s, t);
		elem->GetStress(r, s, t, stress);
		if (!elem->isSacrificial())
		{
			svm = StressMaxCheck(elem, node->Position(), stress);
			if (svm > elSvmMax)
				elSvmMax = svm;
			if (hasAllowable && !checkMaxOnly)
			{
				if (svm > Allowable)
					numNodesYielded++;
			}
			if (!checkMaxOnly)
			{
				if (nodalStressCount.d[id] == 0)
				{
					for (int c = 0; c < 6; c++)
						nodalStress.Put(id, c, stress[c]);
					nodalStressCount.d[id]++;
				}
			}
		}
	}
	for(i = 0; i < ne; i++)
	{
		edge = elem->GetEdge(i);
		id = edge->GetMidNodeId();
		node = model.GetNode(id);
		elem->NodeNaturalCoords(i + nn, r, s, t);
		elem->GetStress(r, s, t, stress);
		if (!elem->isSacrificial())
		{
			svm = StressMaxCheck(elem, node->Position(), stress);
			if (svm > elSvmMax)
				elSvmMax = svm;
			if (hasAllowable && !checkMaxOnly)
			{
				if (svm > Allowable)
					numNodesYielded++;
			}
			if (!checkMaxOnly)
			{
				if (nodalStressCount.d[id] == 0)
				{
					for (int c = 0; c < 6; c++)
						nodalStress.Put(id, c, stress[c]);
					nodalStressCount.d[id]++;
				}
			}
		}
	}

	//also sample face centroids and element centroids:
	SRvec3 pos;
	double rf, sf;
	for (i = 0; i < elem->GetNumLocalFaces(); i++)
	{
		SRface* face = elem->GetFace(i);
		model.map.FaceCentroid(face, rf, sf);
		model.map.ElementNaturalCoordsFromFace(elem, i, rf, sf, r, s, t);
		elem->GetStress(r, s, t, stress);
		if (!elem->isSacrificial())
		{
			face->Position(rf, sf, pos);
			svm = StressMaxCheck(elem, pos, stress);
			if (hasAllowable && !checkMaxOnly)
			{
				if (svm > Allowable)
					numNodesYielded++;
			}

			if (svm > elSvmMax)
				elSvmMax = svm;
		}
	}
	model.map.ElementCentroid(elem, r, s, t);
	elem->GetStress(r, s, t, stress);
	if (!elem->isSacrificial())
	{
		elem->Position(r, s, t, pos);
		svm = StressMaxCheck(elem, pos, stress);
		if (hasAllowable && !checkMaxOnly)
		{
			if (svm > Allowable)
				numNodesYielded++;
		}
		if (svm > elSvmMax)
			elSvmMax = svm;
	}

	elem->SetSvmMax(elSvmMax);
	StressMaxCheckvsAllowable(elem, elSvmMax);

	if (!checkMaxOnly)
	{
		int nnTotal = numNodesTotalTet;
		if (elem->GetType() == brick)
			nnTotal = numNodesTotalBrick;
		else if (elem->GetType() == wedge)
			nnTotal = numNodesTotalWedge;
		double percentYielded = ((double)numNodesYielded) / ((double)nnTotal);
		elem->GetMaterial()->AddToVolPercentYielded(percentYielded);
	}
}

void SRpostProcess::fillSacricialElementNodalStress(SRelement* elem, bool doMaxClipping)
{
	//for sacrifical elements only, fill up stress in nodes not owned by nonsacrifical elements.
	//but clip back to max in the material the element is assigned

	//input:
		//elem = pointer to element
	//note:
		//fills contribution of this element to nodalStress
	int i, nn, ne;
	double r, s, t;
	double stress[6];
	nn = elem->GetNumNodes();
	ne = elem->GetNumLocalEdges();
	SRedge* edge;
	int id;
	elem->SetBasisData();
	SRmaterial* elmat = elem->GetMaterial();
	double maxSvm = elmat->GetMaxSvm();
	double targetRatio = 0.9;
	//sanity check maxSvm should be <= max in model
	if (maxSvm > analysis.GetStressMax())
		maxSvm = analysis.GetStressMax();
	//keep it down to targetRatio so it doesn't show up as hot spot in fringes:
	maxSvm *= targetRatio;
	double ratc[6];
	for (i = 0; i < nn; i++)
	{
		id = elem->GetNodeId(i);
		if (doMaxClipping && nodalStressCount.Get(id) != 0)
			continue;
		elem->NodeNaturalCoords(i, r, s, t);
		elem->GetStress(r, s, t, stress);
		double vm = model.math.GetSvm(stress);
		double ratvm = 1.0;
		for (int j = 0; j < 6; j++)
			ratc[j] = 1.0;
		if (doMaxClipping)
		{
			if (vm > maxSvm && vm > TINY)
				ratvm = maxSvm / vm;
			for (int c = 0; c < 6; c++)
			{
				double maxstress = targetRatio*fabs(analysis.GetStressMaxComp(c));
				double fstress = fabs(stress[c]);
				if (fstress > TINY)
				{
					ratc[c] = maxstress / fstress;
					if (ratvm < ratc[c])
						ratc[c] = ratvm;
				}
			}
		}

		nodalStressCount.Put(id, 1);
	}
	for (i = 0; i < ne; i++)
	{
		edge = elem->GetEdge(i);
		if (edge->GetPorder() < 2)
			continue;
		id = edge->GetMidNodeId();
		if (doMaxClipping && nodalStressCount.Get(id) != 0)
			continue;

		elem->NodeNaturalCoords(i + nn, r, s, t);
		elem->GetStress(r, s, t, stress);
		double vm = model.math.GetSvm(stress);
		double ratvm = 1.0;
		for (int j = 0; j < 6; j++)
			ratc[j] = 1.0;
		if (doMaxClipping)
		{
			if (vm > maxSvm && vm > TINY)
				ratvm = maxSvm / vm;
			for (int c = 0; c < 6; c++)
			{
				double maxstress = targetRatio*fabs(analysis.GetStressMaxComp(c));
				double fstress = fabs(stress[c]);
				if (fstress > TINY)
				{
					ratc[c] = maxstress / fstress;
					if (ratvm < ratc[c])
						ratc[c] = ratvm;
				}
			}
			for (int c = 0; c < 6; c++)
			{
				stress[c] *= ratc[c];
				nodalStress.Put(id, c, stress[c]);
			}

			for (int c = 0; c < 6; c++)
			{
				stress[c] *= ratc[c];
				nodalStress.Put(id, c, stress[c]);
			}

			nodalStressCount.Put(id, 1);
		}
		nodalStressCount.Put(id, 1);
	}
}

double SRpostProcess::StressMaxCheck(SRelement* elem, SRvec3& pos, double stress[])
{
	//check if stress exceeds max in model
	//input:
		//pos = xyz position
		//stress = stress tensor stored as vector

	double svm = model.math.GetSvm(stress);
	if (svm> analysis.GetStressMax())
		analysis.SetStressMax(elem, pos, svm);
	analysis.SetStressMaxComp(stress);
	analysis.UpdateCustomCriterion(stress);
	return svm;
}

void SRpostProcess::StressMaxCheckvsAllowable(SRelement* elem, double svmMax)
{
	//check the max stress in an element vs the allowable stress for the material assigned to the element
	//print warning if it is exceeded
	SRmaterial* mat = elem->GetMaterial();
	if (!elem->isSacrificial())
	{
		if (svmMax > mat->GetMaxSvm())
			mat->PutMaxSvm(svmMax);
	}
	if (mat->GetHighStressWarned())
		return;
	double allowstress = mat->GetAllowableStress();
	if (allowstress < TINY)
		return;
	if (svmMax > allowstress)
	{
		LOGPRINT(" warning. allowable stress exceeded for material: %s", mat->GetName());
		mat->SetHighStressWarned(true);
		analysis.setanyMatYielded(true);
	}
}

void SRpostProcess::CalculateMaxStress()
{
	//check for max stress during p loop
	//note:
		//this uses checkMaxOnly with output set to true
		//so stresses are sampled but not output.
		//CalculateElementStresses samples each element node (corner and midside)

	int i;
	SRelement* elem;
	analysis.zeroStressMax();

	bool checkmaxonly = true;

	for (i = 0; i < model.GetNumElements(); i++)
	{
		elem = model.GetElement(i);
		if (!elem->isSacrificial())
			CalculateElementStresses(elem, checkmaxonly);
	}
}

int SRpostProcess::MeshToFrd()
{
	//write the mesh definition to the Cgx frd file
	//mesh definition, nodes and elements:
	//local edge number mapping between SR and cgx:
	int brickSRtoCgx[12] = { 0, 1, 2, 3, 8, 9, 10, 11, 4, 5, 6, 7 };
	int wedgeSRtoCgx[9] = { 0, 1, 2, 6, 7, 8, 3, 4, 5 };
	int tetSRtoCgx[6] = { 0, 1, 2, 3, 4, 5 };

	FRDPRINT("    1C"); // ''1C'' defines a new calc
	FRDPRINT("    1UDATE   26.march.2000"); // ''1U'' stores user job - information, can be any string, ttd put real date, add more lines, e.g. PGM stressrefine, model file path 
	//nodes:
	int i, j;
	SRnode* node;
	int numnode = model.GetNumNodes();
	int numNodeOut = 0;
	for (i = 0; i < numnode; i++)
	{
		node = model.GetNode(i);
		if (node->isOrphan())
			continue;
		numNodeOut++;
	}
	FRDPRINT("    2C                  %12d                                     1", numNodeOut);
	for (i = 0; i < numnode; i++)
	{
		node = model.GetNode(i);
		if (node->isOrphan())
			continue;
		//http://bconverged.com/calculix/doc/cgx/html/node167.html
		//long Format:(1X,'-1',I10,3E12.5)
		FRDPRINT(" -1%10d%12.5le%12.5le%12.5le", node->GetUserid(), node->GetXyz(0), node->GetXyz(1), node->GetXyz(2));
	}
	FRDPRINT(" -3");
	//elements. eltype: 10 node tet 6, wedge 5, brick 4
	int numel = model.GetNumElements();
	FRDPRINT("    3C                  %12d                                     1", numel);
	int *srToCgx;
	for (i = 0; i < numel; i++)
	{
		SRelement* elem = model.GetElement(i);
		int type;
		bool needContinue = false;
		if (elem->GetType() == brick)
		{
			needContinue = true;
			type = 4;
			srToCgx = brickSRtoCgx;
		}
		else if (elem->GetType() == wedge)
		{
			needContinue = true;
			type = 5;

			srToCgx = wedgeSRtoCgx;
		}
		else
		{
			type = 6;
			srToCgx = tetSRtoCgx;
		}
		FRDPRINT(" -1%10d%5d%5d%5d", elem->GetUserid(), type, 0, 1);
		FRDPRINTNORET(" -2");
		int nn = elem->GetNumNodes();
		for (j = 0; j < nn; j++)
			FRDPRINTNORET("%10d", elem->GetNode(j)->GetUserid());
		for (j = 0; j < elem->GetNumLocalEdges(); j++)
		{
			int edgeNum = srToCgx[j];
			FRDPRINTNORET("%10d", elem->GetEdge(edgeNum)->GetMidNodeUserId());
			nn++;
			if (nn == 10 && needContinue)
			{
				FRDPRINTRET;
				FRDPRINTNORET(" -2");
			}
		}
		FRDPRINTRET;
	}
	FRDPRINT(" -3");
	return numNodeOut;
}

double SRpostProcess::GlobalStrainSmooth()
{
	//global strain smoothing
	//since strain continuity is not assured across material interfaces,
	//smooth strains in groups of elements with same material

	int n = analysis.GetNumFunctions();
	model.AllocateSmoothFunctionEquations(n);
	analysis.AllocateSkipFun(n);

	int numel = model.GetNumElements();
	elSmooth.Allocate(numel);
	double strainMax = 0.0;
	SRelement* elem;
	if (analysis.getAllMatsEquivElast())
	{
		for (int j = 0; j < numel; j++)
		{
			elem = model.GetElement(j);
			elSmooth.Put(j, 1);
		}
		strainMax = ElementStrainSmoothByMaterial();
		double *strainVecs[6];
		for (int c = 0; c < 6; c++)
			strainVecs[c] = getSmoothedStrainVec(c);
		for (int j = 0; j < numel; j++)
		{
			elem = model.GetElement(j);
			elem->DownloadSmoothedStrains(strainVecs);
			elem->checkMaterialStrainMax(strainMax);
		}
		return strainMax;
	}

	for (int i = 0; i < model.GetNumMaterials(); i++)
	{
		double strainMaxMat = 0.0;
		elSmooth.Zero();
		bool anyElToSmooth = false;
		for (int j = 0; j < numel; j++)
		{
			elem = model.GetElement(j);
			if (elem->GetMaterialId() == i)
			{
				elSmooth.Put(j, 1);
				anyElToSmooth = true;
			}
		}

		if (anyElToSmooth)
		{
			strainMaxMat = ElementStrainSmoothByMaterial();
			if (strainMaxMat > strainMax)
				strainMax = strainMaxMat;

			double *strainVecs[6];
			for (int c = 0; c < 6; c++)
				strainVecs[c] = getSmoothedStrainVec(c);
			for (int j = 0; j < numel; j++)
			{
				elem = model.GetElement(j);
				if (elSmooth.Get(j) == 1)
				{
					elem->DownloadSmoothedStrains(strainVecs);
					elem->checkMaterialStrainMax(strainMaxMat);
				}
			}
		}
	}

	return strainMax;
}
double SRpostProcess::ElementStrainSmoothByMaterial()
{
	//interpolate Strain with same basis functions as were used in last solution
	//for displacement. Determine coefficients by least-squares fitting to "raw
	//strains"
	//return:
		//maximum "raw" strain in model
	//notes:
		//this is called once for each material. The elements that have the currently active material
		//are stored in model.elSmooth.

	//sacrificial elements are smoothed together with other elements of same material to avoid discontinuity at 
	//boundary. sacrificial elements are frozen at p2 to minimize effects of singularities
#ifdef NOSOLVER
	return 0.0;
#else
	int e, i, j, nint, gp, gi, gj, stiffLoc, comp, symloc, elmax;
	double r, s, t, w, bi, bj, detJ, strain[6], ae, strainMax = 0.0, stress[6], svm, svmmax = 0.0;
	double eT, eTmax = 0.0;
	SRelement* elem;

	analysis.NumberEquationsSmooth();

	int n = analysis.GetNumSmoothEquations();
	if (n == 0)
		return 0.0;

	SRsolver* sol = &(analysis.solver);
	SRpardiso* par = NULL;
	par = sol->parDisoPtr();
	par->smoothBookkeep();
	double* strainRhsvec[6];

	for (i = 0; i < 6; i++)
	{
		smoothedStrains[i].Allocate(n);
		strainRhsvec[i] = smoothedStrains[i].GetVector();
	}
	int nfunel;
	int len, eqi, eqj;

	for (e = 0; e < model.GetNumElements(); e++)
	{
		double elSvmMax = 0.0;
		double* stiff = analysis.GetElementStiffnessVector();
		elem = model.GetElement(e);
		elem->SetBasisData();
		elem->SetSvmMax(0.0);
		nfunel = elem->GetNumFunctions();
		elem->FillStiffDiag(nfunel);
		nint = model.math.FillGaussPoints(elem);
		if (elSmooth.Get(e) == 0)
			continue;
		elem->AllocateRawStrains(nint);
		len = nfunel*(nfunel + 1) / 2;
		for (i = 0; i < len; i++)
			stiff[i] = 0.0;
		double etx, ety, etz;
		for (gp = 0; gp < nint; gp++)
		{
			model.math.GetGP3d(gp, r, s, t, w);
			detJ = elem->FillMapping(r, s, t);
			double* basisvec = elem->FillBasisFuncs(r, s, t, both);

			eT = elem->CalculateRawStrain(r, s, t, strain, etx, ety, etz);

			if (eT > eTmax)
				eTmax = eT;
			elem->StraintoStress(r, s, t, strain, stress);
			svm = model.math.GetSvm(stress);
			if (svm > elSvmMax)
				elSvmMax = svm;
			if (svm > elem->GetSvmMax())
				elem->SetSvmMax(svm);
			if (!elem->isSacrificial())
			{
				for (comp = 0; comp < 6; comp++)
				{
					ae = fabs(strain[comp]);
					if (ae > strainMax)
						strainMax = ae;
				}
			}

			elem->PutRawStrain(gp, strain);
			w *= detJ;
			for (i = 0; i < nfunel; i++)
			{
				gi = elem->GetFunctionNumber(i);
				eqi = model.GetSmoothFunctionEquation(gi);
				bi = basisvec[i] * w;
				for (comp = 0; comp < 6; comp++)
					strainRhsvec[comp][eqi] += (strain[comp] * bi);
				for (j = i; j < nfunel; j++)
				{
					bj = basisvec[j];
					symloc = elem->GetStiffnessLocationAboveDiag(i, j);
					stiff[symloc] += bi*bj;
				}
			}
		}

		if (elSvmMax > svmmax)
		{
			svmmax = elSvmMax;
			elmax = e;
		}

		par->smoothAssemble(e, stiff);
	}

	bool smoothing = true;
	par->solve(smoothing);
	par->clear();

	return strainMax;
#endif
}
