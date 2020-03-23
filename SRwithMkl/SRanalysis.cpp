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
// SRanalysis.cpp: implementation of the SRanalysis class.
//
//////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include "SRanalysis.h"
#include "SRinput.h"
#include "globalWrappers.h"
#ifndef NOSOLVER
#include "mkl.h"
#endif

extern SRmodel model;

//////////////////////////////////////////////////////////////////////

SRanalysis::SRanalysis()
{
	Initialize();
};

void SRanalysis::Initialize()
{
	anybricks = false;
	anytets = false;
	anywedges = false;
	detectSacrificialElements = true;
	anyEnforcedDisplacement = false;
	maxElemUid = 0;
	adaptIt = 0;
	elUidAtMaxError = -1;
	allConstraintsAsPenalty = false; //change to true to use penalty constraints for all constraints
	allMatsHaveAllowable = true;
	maxAllowableAnyActiveMat = 0.0;
	strainMax = 0.0;
	numFlattenedElHighStress = 0;
	allMatsEquivElast = true;
	anyMatYielded = false;
	prevStressMax = 0.0;
	maxsp1 = 0.0;
	minsp2 = BIG;
	maxCustom = 0.0;
	nodeidAtMaxStressPos = -1;
	maxStressElid = -1;
	needSoftSprings = false;
	stressUnitConversion = 1.0;
	stressUnitstr = "Pa (N/m^2)";
	lengthUnitConversion = 1.0;
	lengthUnitstr = "m";
	useUnits = true;
	checkForHotSpotatMax = true;
	anyHotSpotElFound = false;
	outputf06 = false;
	doEnergySmooth = true;
	pOrderUniform = false;
	ErrorTolerance = ERRTOL;
	maxPJump = 10;
	maxPorder = 5;
	maxPorderFinalAdapt = MAXP;
	maxPorderLowStress = 6;
	maxPinModel = 2;
	adaptLoopMax = MAXADAPTLOOPS;
	numSacrElem = 0;
	initializeErrorCheck();
}


static const char *stressComp[6] = { "xx", "yy", "zz", "xy", "xz", "yz"};

void SRanalysis::Run()
{
	//run analysis on this model

	//1. initial operations outside p-loop
	//2. p-loop:
		//function numbering
		//process constraints
		//equation numbering
		//calculate element stiffnesse
		//process forces
		//process enforced displacements
		//solve:
			//assemble and decomp global matrix
			//backsolve force to calculate displacements
		//adapt:
			//calculate raw stresses and perform smoothing
			//calculate errors:
				//1. smoothed vs raw stresses
				//2. traction jumps across shared faces
				//3. traction jumps vs applied loads
			//use errors to determine next p level
		//loop until error within tolerance, all edges at max p, or max iterations reached
	//3. post-process
		//output mesh, displacement results file, and stress results file. output summary quantities

	//note: this outputs to three places in the same folder the input (.msh) file is in:
		//report.txt is a summary that gets displayed by the UI. This output is done by "REPPRINT"
		//filname.log, where filename is the name of the input file, is a detailed log. It will be displayed by the ui if anything has gone wrong
		//This is done by "LOGPRINT"
		//filename.out is more detailed results output. It is done by "LOGPRINT"
		//SCREENPRINT is to the command window which you'll see if you run this standalone. If run through the UI that window is suppressed and the SCREENPRINT
		//output is redirected to the status bar.

	double errForOutput = 0.0;
	if (pOrderUniform)
		passData.Allocate(8);
	else
		passData.Allocate(adaptLoopMax);

	int nits;
	int lastIt = 0;
	SRvec3 dmax;
	int nodeUidAtMaxDisp = 0;

	int its;

	bool anyLcsEnfd = PreProcessPenaltyConstraints();

	if (pOrderUniform)
	{
		maxPorder = maxPorderFinalAdapt;
		nits = maxPorder - 1;
	}
	else
	{
		nits = adaptLoopMax;
		if (maxPorder == 2)
			nits = 1;
	}

	bool adapted = true;

	nodeUidAtMaxDisp = 0;

	for (its = 0; its < nits; its++)
	{
		adaptIt = its;

		LOGPRINT("\nAdaptive Solution. Iteration: %d\n", its + 1);
		LOGPRINT("Maximum Polynomial Order in Model: %d\n", maxPinModel);
		model.setVolume(0.0);

		numberGlobalFunctions();

		ProcessConstraints();

		numEquations = numberEquations();

		allocateSolutionVector();

		model.allocateSmallElementData(anyLcsEnfd);

		checkElementMapping();

		LOGPRINT("Calculating Element Stiffness\n");

		CheckElementsFitInMemory();

		SCREENPRINT("Adaptive Iteration: %d. Calculating element stiffnesses\n", its + 1);
		CalculateElementStiffnesses(anyLcsEnfd);

#ifdef NOSOLVER
		SCREENPRINT("Cannot proceed further without solver\n");
		LOGPRINT("Cannot proceed further without solver\n");
		exit(0);
#endif

		if (needSoftSprings)
			AddSoftSprings();

		ProcessForces();

		EnforcedDisplacementAssemble();

		//solution:
		SCREENPRINT("Adaptive Iteration: %d. Solving Equations\n", its + 1);
		LOGPRINT("Solving %d Equations\n", numEquations);
		LOGPRINT("\nNumber of Equations: %d\n", numEquations);

		solver.DoSolution();

		NodalMaxDisp(nodeUidAtMaxDisp, dmax);
		LOGPRINT("\nMaximum nodal displacement (at node %d)\n", nodeUidAtMaxDisp);
		LOGPRINT("ux: %lg\nuy: %lg\nuz: %lg\n\n", dmax.d[0], dmax.d[1], dmax.d[2]);

		int pbeforeAdapt = maxPinModel;

		SCREENPRINT("Adaptive Iteration: %d. Stress Smoothing and Calculating Stress Error\n", its + 1);
		LOGPRINT("Stress Smoothing and Calculating Stress Error\n");
		if (its != (nits - 1))
			adapted = Adapt(its);
		else
			adapted = Adapt(its, CHECKERRORONLY);//final pass. CheckErroronly is true

		errForOutput = CalculateMaxErrorForOutput();
		LOGPRINT("Estimated Error in Stress Calculation: %6.2lg percent", errForOutput);
		LOGPRINT("Maximum von Mises Stress in Model: %lg\n", stressUnitConversion*stressMax);

		SRPassData& pd = passData.Get(its);
		pd.err = errForOutput;
		pd.maxp = pbeforeAdapt;
		pd.neq = numEquations;
		pd.maxsvm = stressMax;
		pd.maxsp1 = maxsp1;
		pd.maxCustom = maxCustom;

		prevStressMax = stressMax;

		model.FreeSmallElementData();

		lastIt++;

		if (!adapted)
			break;
	}

	bool slopekinkWarnNeeded = false;
	if (maxStressElid != -1)
	{
		SRelement* elemAtMax = model.GetElement(maxStressElid);
		slopekinkWarnNeeded = elemAtMax->checkSlopeKink();
	}

	//clean up memory and disk space that is only needed during adaptive solution:
	CleanUp(true);//partial is true

	model.allocateSmallElementData(numEquations);
	LOGPRINT("\nAdaptive loop complete");
	LOGPRINT("\nAdaptive loop complete");
	LOGPRINT("Post Processing");
	post.PostProcess();

	bool flattenedWarningNeeded = checkFlattenedElementHighStress();

	model.freeNodeFaces();

	SRstring line;
	REPPRINT("Results of StressRefine full model solution\n");

	if (useUnits)
	{
		REPPRINT("Units of Stress: %s\n", stressUnitstr.getStr());
		REPPRINT("Units of Length: %s\n", lengthUnitstr.getStr());
	}
	else
	{
		REPPRINT("\nUnits of Stress are the same as units for Young's Modulus in Nastran MAT1\n");
	}

	if (detectSacrificialElements)
	{
		if (SingStressCheck())
			REPPRINT("Singular\n");
	}
	if (checkForSmallMaxStress())
		REPPRINT("SmallStressDetected\n");
	if (flattenedWarningNeeded)
		REPPRINT("Flattened High StressElement\n");
	if ((errForOutput) > 10.5)//10.5 instead of 10 is for roundoff
		REPPRINT("HighError: %lg\n", errForOutput);
	if (slopekinkWarnNeeded)
		REPPRINT("SlopeKinkAtMax\n");

	REPPRINT("MaxVMIts");
	for (int i = 0; i < lastIt; i++)
		REPPRINT(",%lg", passData.Get(i).maxsvm);
	REPPRINT("\n");
	REPPRINT("MaxSP1Its");
	for (int i = 0; i < lastIt; i++)
		REPPRINT(",%lg", passData.Get(i).maxsp1);
	REPPRINT("\n");
	REPPRINT("MaxCustomIts");
	for (int i = 0; i < lastIt; i++)
		REPPRINT(",%lg", passData.Get(i).maxCustom);
	REPPRINT("\n");

	REPPRINT("Maximum von Mises Stress in Model: %12.3lg\n  ", stressUnitConversion*stressMax);
	REPPRINT("Maximum Principal Stress in Model: %12.\n  ", stressUnitConversion*maxsp1);
	REPPRINT("Minimum Principal Stress in Model: %12.3lg\n  ", stressUnitConversion*minsp2);

	REPPRINT("Estimated Error in Stress Calculation: %6.2lg percent\n", errForOutput);
	double maxPercentYielded = 0.0;
	int nmat = model.GetNumMaterials();
	int numactivemat = 0;
	for (int i = 0; i < nmat; i++)
	{
		SRmaterial* mat = model.GetMaterial(i);
		if (mat->isActive())
			numactivemat++;
	}
	if (numactivemat > 1 && !allMatsEquivElast)
	{
		for (int i = 0; i < nmat; i++)
		{
			SRmaterial* mat = model.GetMaterial(i);
			if (mat->isActive())
			{
				REPPRINT("For Elements with Material: %s\n", mat->GetName());
				double vm = mat->GetMaxSvm();
				REPPRINT("    Max von Mises Stress: %12.3lg\n", stressUnitConversion*vm);
				if (mat->isAllowableAssigned())
				{
					double percentYielded;
					percentYielded = mat->GetVolPercentYielded();
					double fs = 1000.0;
					double allow = mat->GetAllowableStress();
					if (vm > SMALL)
						fs = allow / vm;
					if (fs < 1.0)
					{
						REPPRINT("    Allowable Stress: %lg Factor of Safety: %5.2lg -- Yielded\n", stressUnitConversion*allow, fs);
						if (percentYielded > 1.0)
							percentYielded = 0.9999;
						REPPRINT("    Percentage of Material volume yielded: %5.2lf\n", (percentYielded*100.0));
						if (percentYielded > maxPercentYielded)
							maxPercentYielded = percentYielded;
					}
					else
						REPPRINT("    Allowable Stress: %lg Factor of Safety: %lg\n", stressUnitConversion*allow, fs);
				}
			}
		}
	}

	REPPRINT("Maximum displacement in Model: %lg\n  ", lengthUnitConversion*dmax.Magnitude());
	REPPRINT("\nMaximum Polynomial Order in Model: %d\n", maxPinModel);
	REPPRINT("Number of in Equations in Model: %d\n\n", numEquations);
	if (nits != 1)
	{
		REPPRINT("Adaptive Analysis     Max Polynomial Order    Estimated Stress Error    Max Von Mises Stress\n");
		REPPRINT("Iteration                                                               in Model\n");

		for (int its = 0; its < lastIt; its++)
		{
			SRPassData& pd = passData.Get(its);
			REPPRINT("         %d                      %d                  %-6.2lg                    %12.3lg\n", its + 1, pd.maxp, pd.err, stressUnitConversion*pd.maxsvm);
		}
	}

	model.repFile.Close();

	LOGPRINT("\n\nFinal result for Model: %s\n", fileNameTail.getStr());
	LOGPRINT("Total Volume of Model : %lg\n", model.getVolume());
	PStats();
	LOGPRINT("Estimated Error in Stress Calculation: %6.2lg percent\n", errForOutput);
	LOGPRINT("Maximum Polynomial Order in Model: %d\n", maxPinModel);
	LOGPRINT("Maximum nodal displacement (at node %d)\n", nodeUidAtMaxDisp);
	LOGPRINT("ux: %lg\nuy: %lg\nuz: %lg\n", dmax.d[0], dmax.d[1], dmax.d[2]);
	LOGPRINT("Maximum von Mises Stress in Model: %lg  at position %lg %lg %lg\n", stressMax, maxStressPos.d[0], maxStressPos.d[1], maxStressPos.d[2]);
	for (int i = 0; i < 6; i++)
		LOGPRINT("Maximum Stress Component %s in Model: %lg\n", stressComp[i], stressMaxComp[i]);

	if (model.getMaxFlattened() > TINY)
	{
		LOGPRINT("max element flattening for bad mapping: %lg\n", model.getMaxFlattened());
		LOGPRINT("element: %d\n", model.getMaxFlattenedElUid());
	}

	CleanUp();
}

bool SRanalysis::Adapt(int pIteration, bool checkErrorOnly)
{
	//adapt polynomial orders of all elements after error checking;
	//input:
		//pIteration = p-loop iteration number
		//checkErrorOnly = true to check error only but not increase p
	//notes:
		//adaptivity steps:
		//calculate raw strains and perform smoothing (in errorChecker.SetUp)
		//calculate errors
			//1. traction jumps across shared faces
			//2. traction jumps vs applied loads
		//use errors to determine next p level (in model.FindNextP)

	int i, e;
	SRedge *edge;
	SRelement *elem;
	int p;
	bool pup = false, anypup = false;
	errorMax = 0.0;

	int maxp0 = maxPinModel;

	bool finalAdapt = (pIteration == adaptLoopMax - 2);

	if (finalAdapt)
		maxPorder = maxPorderFinalAdapt;
	else if (maxPorder > maxPorderFinalAdapt)
		maxPorder = maxPorderFinalAdapt;


	for (int e = 0; e < model.GetNumElements(); e++)
	{
		SRelement* elem = model.GetElement(e);
		elem->DownloadDisplacement();
	}
	post.GlobalStrainSmooth();
	post.CalculateMaxStress();
	setupErrorCheck(finalAdapt, stressMax, ErrorTolerance, maxPorder, maxPJump, maxPorderLowStress);


	if (pIteration == 0 && detectSacrificialElements)
	{
		numSacrElem = CheckAutoSacrificialElements();
		if (numSacrElem > 0)
		{
			//this has changed sacrificial element status, so redo max calculation:
			post.CalculateMaxStress();
		}
	}

	if (!checkErrorOnly)
	{
		LOGPRINT("\n");
		LOGPRINT("Maximum von Mises Stress in Model: %lg\n", stressMax);
		LOGPRINT("Maximum Stress Components in Model\n");
		PrintStresses(stressMaxComp);
	}

	//pOrderUniform is off by default, can be turned on by a flag in the settings file
	if(pOrderUniform)
	{
		//AdaptUniform will increase p of all edges in model by 1 if any error exceeds tolerance
		anypup = AdaptUniform(pIteration, checkErrorOnly);
		if (checkErrorOnly)
		{
			CleanUpErrorCheck();
			return false;
		}
		for (e = 0; e < model.GetNumElements(); e++)
		{
			elem = model.GetElement(e);
			elem->SetPChanged(true);
		}

		CleanUpErrorCheck();
		return anypup;
	}
	int maxpSacr = 0;
	for (e = 0; e < model.GetNumElements(); e++)
	{
		elem = model.GetElement(e);
		elem->SetPChanged(false);
		pup = FindNextP(elem->GetId(), &p);
		if (pup)
		{
			if (!elem->isSacrificial())
			{
				anypup = true;
				if (p > maxPinModel)
					maxPinModel = p;
			}
			else if (p > maxpSacr)
				maxpSacr = p;
			elem->PutNewPorder(p);
			elem->SetPChanged(true);
		}
		else
			elem->PutNewPorder(0);
	}

	if (checkErrorOnly)
	{
		CleanUpErrorCheck();
		maxPinModel = maxp0;
		return false;
	}
	else
	{
		if (anypup && maxpSacr > maxPinModel)
			maxPinModel = maxpSacr;
	}

	if (!anypup)
	{
		CleanUpErrorCheck();
		return anypup;
	}

	for (e = 0; e <model.GetNumElements(); e++)
	{
		elem = model.GetElement(e);
		p = elem->GetNewPorder();
		//update P-order of all the element's edges:
		for (i = 0; i < elem->GetNumLocalEdges(); i++)
		{
			edge = elem->GetEdge(i);
			if (p > edge->GetPorder())
				edge->putPorder(p);
		}
	}

	for (e = 0; e < model.GetNumElements(); e++)
	{
		elem = model.GetElement(e);
		if (elem->GetPChanged())
			continue;
		for (i = 0; i < elem->GetNumLocalEdges(); i++)
		{
			edge = elem->GetEdge(i);
			if (edge->PChanged())
			{
				elem->SetPChanged(true);
				break;
			}
		}
	}

	CleanUpErrorCheck();
	return anypup;
}

bool SRanalysis::AdaptUniform(int pIteration, bool checkErrorOnly)
{
	//uniform p adaptivity
	//input:
		//pIteration = p-loop iteration number
		//checkErrorOnly = true to check error only but not increase p
	//notes:
		//finds max error in model. if less than tolerance, return "any-p-increased" = false
		//else:
			//increase p of all edges in model by 1 unless they are already at max,
	//return:
		//"any-p-increased" = true

	int e;
	SRelement* elem;
	SRedge* edge;

	for (e = 0; e < model.GetNumElements(); e++)
	{
		elem = model.GetElement(e);
		if (elem->isSacrificial())
			continue;
		
		double error = FindElementError(elem->GetId());
		setErrorMax(error, elem->GetUserid());
	}

	errorMax = CalculateMaxErrorForOutput();
	errorMax *= 0.01;//converting from percentage to actual

	if (checkErrorOnly)
	{
		CleanUpErrorCheck();
		return false;
	}

	if (errorMax < ErrorTolerance)
		return false;

	SRintVector edgeupped;
	int id, i, p;
	edgeupped.Allocate(model.GetNumEdges());
	for (e = 0; e < model.GetNumElements(); e++)
	{
		elem = model.GetElement(e);
		if (elem->isSacrificial())
			continue;
		for (i = 0; i <elem->GetNumLocalEdges(); i++)
		{
			id = elem->GetLocalEdgeGlobalId(i);
			if (edgeupped.Get(id))
				continue;
			edgeupped.Put(id, 1);
			edge = model.GetEdge(id);
			p = edge->GetPorder();
			p++;
			if (p <= maxPorder)
			{
				edge->putPorder(p);
				if (p > maxPinModel)
					maxPinModel = p;
			}
		}
	}
	edgeupped.Free();
	return true;
}


void SRanalysis::CalculateElementStiffnesses(bool anyLcsEnfd)
{
	//calculate elemental stiffness matrix for all elements in model
	//all elements must fit in memory

	double* enfdForceVec = NULL;
	if(anyLcsEnfd)
		enfdForceVec = model.getEnforcedVec();

	//scratch space needed by elements:
	model.allocateElementData();

	int n = model.GetNumElements();
	LOGPRINT("Calculating %d elements", n);

	int dprogTarget = (int)(0.1* (double)n);
	int progTarget = dprogTarget;
	int progout = 10;
	for (int i = 0; i < n; i++)
	{
		if (i > progTarget)
		{
			progTarget += dprogTarget;
			progout += 10;
		}

		SRelement* elem = model.GetElement(i);
		double* elStiff = CalculateStiffnessMatrix(i, true, enfdForceVec);

		if (elStiff == NULL)
		{
			LOGPRINT("Not enough memory for Element Stiffness matrices\nMemory required: %lg\n", MaxElementMem);
			ERROREXIT;
		}
	}

	LOGPRINT("\n");

	for (int i = 0; i < n; i++)
	{
		SRelement* elem = model.GetElement(i);
		model.updateVolume(elem->GetVolume());
	}

	model.FreeElementData();

	if (anyLcsEnfd)
	{
		for (int j = 0; j < numEquations; j++)
			solution.PlusAssign(j, enfdForceVec[j]);
	}
}

void SRanalysis::CreateModel()
{
	//create model by processing input file
	input.ReadModel();

	for (int i = 0; i < model.GetNumMaterials(); i++)
	{
		SRmaterial* mat = model.GetMaterial(i);
		if (mat->isActive())
		{
			if(!mat->isAllowableAssigned())
				allMatsHaveAllowable = false;
			else
			{
				if (mat->GetAllowableStress() > maxAllowableAnyActiveMat)
					maxAllowableAnyActiveMat = mat->GetAllowableStress();
			}
		}
	}
}

void SRanalysis::AllocateSmoothFunctionEquations(int n)
{
	model.AllocateSmoothFunctionEquations(n);
};


void SRanalysis::SetStressMax(SRelement* elem, SRvec3& pos, double svm)
{
	//set the value of max stress in model, at corresponding element id and position number
	//input:
		//pos = position of stress that is greater than current max
		//svm = von Mises stress value
	maxStressPos.Copy(pos);
	stressMax = svm;
	maxStressElid = elem->GetId();
}

void SRanalysis::UpdateCustomCriterion(double stress[])
{
	double custom;

	//Replace this with whatever stress criterion is desired
	{
		custom = model.math.GetSvm(stress);
	}


	if (custom > maxCustom)
		maxCustom = custom;
}

void SRanalysis::SetStressMaxComp(double *stressComp)
{
	//check for any components of a stress tensor exceed max stress in model
	//input:
		//stressComp = stress tensor stored as vector
	//note:
		//updates class variable stressMaxComp

	for (int i = 0; i < 6; i++)
	{
		if (fabs(stressComp[i]) > fabs(stressMaxComp[i]))
			stressMaxComp[i] = stressComp[i];
	}
	double sp1, sp2;
	model.math.GetPrinStress(stressComp, sp1, sp2);
	if (sp1 > maxsp1)
		maxsp1 = sp1;
	if (sp2 < minsp2)
		minsp2 = sp2;
}




double SRanalysis::GetDisplacementCoeff(int gfun, int dof)
{
	//look up displacement coefficient for a global function and dof
	//input:
		//gfun = global function number
		//dof = degree of freedom
	//return:
		//displacement coefficient
	return model.GetDisplacementCoeff(gfun, dof);
}

void SRanalysis::NumberEquationsSmooth()
{
	//number global equations for smoothing.
	//notes:
		//there is only one dof.
		//skip functions not belonging to elements in current material group
		//(determined by post.elSmooth flag)
	int i, e, gfun;
	SRelement* elem;
	int nfun = GetNumFunctions();
	for(i = 0; i < nfun; i++)
		skipFun.Put(i, 1);
	for (e = 0; e <model.GetNumElements(); e++)
	{
		elem = model.GetElement(e);
		if (post.getElSmooth(e))
		{
			for(i = 0; i < elem->GetNumFunctions(); i++)
			{
				gfun = elem->GetFunctionNumber(i);
				skipFun.Put(gfun, 0);
			}
		}
	}

	int eq = 0;
	for(i = 0; i < nfun; i++)
	{
		if(skipFun.Get(i) == 0)
		{
			model.PutSmoothFunctionEquation(i, eq);
			eq++;
		}
		else
			model.PutSmoothFunctionEquation(i, -1);
	}
	numSmoothEquations = eq;
}


void SRanalysis::ProcessConstraints()
{
	//process constraints
	//notes:
		//number functionEquations(fun,dof) as negative for a gcs constrained global dof.
		//the number is -1 if there is no enforced displacement, else
		//it is the number in the enforcedDisplacements matrix where the enforced displacement is stored
		//class variables functionEquations and enforcedDisplacements are updated
		//lcs constraints are skipped because they are handled via penalty

	int i, dof, gfun;
	SRconstraint* con;
	SRnode* node;

	model.AllocateDofVectors();

	for (i = 0; i < model.GetNumConstraints(); i++)
	{
		con = model.GetConstraint(i);

		//skip lcs constraint, they are handled via penalty:
		if (!con->isGcs())
			continue;
		int eid = con->GetEntityId();
		if (con->GetType() == nodalCon)
		{
			node = model.GetNode(eid);
			gfun = node->GetGlobalFunctionNumber();
			for (dof = 0; dof < 3; dof++)
			{
				if (con->IsConstrainedDof(dof))
				{
					model.PutFunctionEquation(gfun, dof, -1);
					if (con->hasEnforcedDisp())
					{
						double enfd = con->getDisp(0, dof);
						model.PutEnforcedDisplacement(gfun, dof, enfd);
					}
				}
			}
		}
		else if (con->GetType() == faceCon)
		{
			con->ProcessFaceConstraint();
		}
	}

	for (i = 0; i < model.GetNumElements(); i++)
	{
		SRelement* elem = model.GetElement(i);
		for (int f = 0; f < elem->GetNumFunctions(); f++)
		{
			int gfun = elem->GetFunctionNumber(f);
			for (dof = 0; dof < 3; dof++)
			{
				if (GetFunctionEquation(gfun, dof) < 0)
				{
					elem->SetConstrained();
					break;
				}
			}
			if (elem->isConstrained())
				break;
		}
	}

	int nfun = GetNumFunctions();
	funUncon.Allocate(nfun);
	funUncon.Set(1);
	//funUncon is used for solver bookkeeping. it is true only
	//if the function is not owned by any elements that have constraints, even
	//if they are not on this function
	for (i = 0; i < model.GetNumElements(); i++)
	{
		SRelement* elem = model.GetElement(i);
		if (elem->isConstrained())
		{
			for (int f = 0; f < elem->GetNumFunctions(); f++)
			{
				int gfun = elem->GetFunctionNumber(f);
				funUncon.Put(gfun, 0);
			}
		}
	}
}


void SRanalysis::EnforcedDisplacementAssemble()
{
	//assemble contribution of Enforced Displacements to global force vector
	//notes:
		//global force vector is stored in class variable model.solution

	if (!IsAnyEnforcedDisplacement())
		return;


	double* elstiff = NULL;
	for (int el = 0; el < model.GetNumElements(); el++)
	{
		bool stiffhasbeenread = false;
		SRelement* elem = model.GetElement(el);
		for (int collfun = 0; collfun < elem->GetNumFunctions(); collfun++)
		{
			int colgfun = elem->GetFunctionNumber(collfun);
			for (int coldof = 0; coldof < 3; coldof++)
			{
				int eq = model.GetFunctionEquation(colgfun, coldof);
				//equation numbers less than zero are a flag that this
				//dof is constrained:
				if (eq < 0)
				{
					double enfdisp = model.GetEnforcedDisplacement(colgfun, coldof);
					if (!stiffhasbeenread)
					{
						elstiff = ReadElementStiffness(elem);
						stiffhasbeenread = true;
					}
					//this column of element stiffness contributes to the global force vector:
					int colleq = collfun * 3 + coldof;
					for (int rowlfun = 0; rowlfun < elem->GetNumFunctions(); rowlfun++)
					{
						int rowgfun = elem->GetFunctionNumber(rowlfun);
						for (int rowdof = 0; rowdof < 3; rowdof++)
						{
							int roweq = model.GetFunctionEquation(rowgfun, rowdof);
							if (roweq >= 0)
							{
								int rowleq = rowlfun * 3 + rowdof;
								int stiffloc = elem->GetStiffnessLocation(rowleq, colleq);
								double enfdForceVal = -enfdisp*elstiff[stiffloc];
								AddToSolutionVector(roweq, enfdForceVal);
							}
						}
					}
				}
			}
		}
	}
}

void SRanalysis::ProcessForces()
{
	//Process all forces in model: node, edge, face, volume, and thermal
	//and assemble into global force vector
	//notes:
		//global force vector is stored in class variable model.solution

	int nf, eid, dof, fun, eq;
	SRforce* force;
	double* globalForce = GetSolutionVector();
	SRvec3 resF;
	nf = model.GetNumForces();
	for (int f = 0; f < nf; f++)
	{
		force = model.GetForce(f);
		eid = force->GetEntityId();
		if (force->GetType() == nodalForce)
		{
			SRnode* node = model.GetNode(eid);
			if (node->isOrphan())
				continue;
			if (!node->isMidSide())
				fun = node->GetGlobalFunctionNumber();
			else
			{
				int eid = node->GetMidSideEdgeOwner();
				SRedge* edge = model.GetEdge(eid);
				fun = edge->GetGlobalFunctionNumber(2);
			}
			double fv[3];
			for (dof = 0; dof < 3; dof++)
				fv[dof] = force->GetForceVal(0, dof);
			bool needRotate = false;
			SRmat33 R;
			if (force->GetCoordId() != -1)
			{
				needRotate = true;
				SRcoord* coord = model.GetCoord(force->GetCoordId());
				coord->GetRotationMatrix(false, node->Position(), R);
			}

			if (needRotate)
			{
				SRvec3 t = fv;
				for (int i = 0; i < 3; i++)
				{
					fv[i] = 0.0;
					for (int j = 0; j < 3; j++)
						fv[i] += (R.rows[i].d[j] * t.d[j]);
				}
			}

			for (dof = 0; dof < 3; dof++)
			{
				eq = model.GetFunctionEquation(fun, dof);
				if (eq >= 0)
				{
					globalForce[eq] += (fv[dof]);
					resF.d[dof] += fv[dof];
				}
			}
		}
		else if (force->GetType() == edgeForce)
		{
			SRedge* edge = model.GetEdge(eid);
			edge->ProcessForce(force, resF);
		}
		else if (force->GetType() == faceForce)
		{
			SRface* face = model.GetFace(eid);
			face->ProcessForce(force, resF);
		}
	}

	SRthermalForce*tf = model.GetThermalForce();
	if (tf != NULL)
		tf->Process();

	ProcessVolumeForces(resF);
	LOGPRINT("\nResultant load on model:\n");
	LOGPRINT("   Fx: %lg\n", resF.d[0]);
	LOGPRINT("   Fy: %lg\n", resF.d[1]);
	LOGPRINT("   Fz: %lg\n", resF.d[2]);
}


void SRanalysis::ProcessVolumeForces(SRvec3& ResF)
{
	//fill up contribution of volume forces to global force vector
	//input:
		//ResF = resultant force vector on model
	//output:
		//ResF = updated
	//note:
		//global force vector is stored in class variable model.solution
		//volume force applies to entire model is it's element list elList is empty,
		//otherwise it only applies to the elements on the list.

	SRelement* elem;
	SRdoubleVector basvec;
	basvec.Allocate(model.GetmaxNumElementFunctions());
	double* globalForce;
	globalForce = GetSolutionVector();
	SRvolumeForce* vf;
	int nvf = model.getNumVolumeForces();
	for (int v = 0; v < nvf; v++)
	{
		vf = model.getVolumeForce(v);
		int ellistLen = vf->getElListLen();
		if (ellistLen == 0)
		{
			for (int i = 0; i < model.GetNumElements(); i++)
			{
				elem = model.GetElement(i);
				ProcessElemVolumeForce(vf, elem, ResF, basvec, globalForce);
			}
		}
		else
		{
			for (int i = 0; i < ellistLen; i++)
			{
				int eid = vf->getElList(i);
				elem = model.GetElement(eid);
				ProcessElemVolumeForce(vf, elem, ResF, basvec, globalForce);
			}
		}
	}
}

void SRanalysis::ProcessElemVolumeForce(SRvolumeForce* vf, SRelement*elem, SRvec3& ResF, SRdoubleVector& basvec, double* globalForce)
{
	//process contribution of a single element to a volume force
	//input:
	//vf = pointer to volume force
	//elem = pointer to element force
	//scratch:
	//basvec = storage for element basis functions
	//output:
	//ResF has been modified with element's contribution to resultant force on model
	//globalForce has been modified with element's contribution to global force vector

	int nint, gp, fun, dof, eq, nfun, gfun;
	double r, s, t, w, detj;
	double force[3], bw;
	SRvec3 p;

	if (vf->GetType() == gravity)
		vf->GetForceValue(elem, p, force);
	nint = model.math.FillGaussPoints(elem);
	nfun = elem->GetNumFunctions();
	for (gp = 0; gp < nint; gp++)
	{
		model.math.GetGP3d(gp, r, s, t, w);
		detj = elem->FillMapping(r, s, t, true);
		w *= detj;
		model.basis.ElementBasisFuncs(r, s, t, elem, basvec.GetVector());
		if (vf->GetType() == centrifugal)
		{
			elem->Position(r, s, t, p);
			vf->GetForceValue(elem, p, force);
		}

		//ResF += force*w:
		p.Assign(force);
		p.Scale(w);
		ResF.PlusAssign(p);
		SRvec3 pos;
		elem->Position(r, s, t, pos);

		for (fun = 0; fun < nfun; fun++)
		{
			gfun = elem->GetFunctionNumber(fun);
			bw = basvec.Get(fun) * w;
			for (dof = 0; dof < 3; dof++)
			{
				eq = model.GetFunctionEquation(gfun, dof);
				if (eq >= 0)
					globalForce[eq] += (bw*force[dof]);
			}
		}
	}
}


void SRanalysis::CleanUp(bool partial)
{
	//miscelaneous memory and disk cleanup at end of run
	model.CleanUp(partial);
	if (partial)
		return;
	solution.Free();
	skipFun.Free();
	solver.Cleanup();
	input.Cleanup();
}

void SRanalysis::NodalMaxDisp(int& nodeUidAtMaxDisp, SRvec3& umax)
{
	//find max nodal displacement in model
	//output:
		//nodeUidAtMaxDisp = node Uid at which max disp occurred
	//output:
		//umax =3 dof displacement vector with highest magnitude in model
	double umagmax = 0.0;
	umax.Zero();
	SRvec3 u;
	for (int i = 0; i < model.GetNumNodes(); i++)
	{
		SRnode* node = model.GetNode(i);
		node->GetDisp(u);
		double umag = u.Length();
		if (umag > umagmax)
		{
			umagmax = umag;
			umax.Copy(u);
			nodeUidAtMaxDisp = node->GetUserid();
		}
	}
}

void SRanalysis::PrintStresses(double *stress)
{
	//print stress tensor to model the model output file
	//input:
		//stress = stress tensor stored as vector
	for (int i = 0; i < 6; i++)
		LOGPRINT("%s: %lg\n", stressComp[i], stress[i]);
}

void SRanalysis::CheckElementsFitInMemory()
{
	//check if all element stiffness matrices will fit in memory
	//note:
		//sets class variable elementsInMemory
	MaxElementMem = 0.0;
	int numel = model.GetNumElements();
	for (int e = 0; e < numel; e++)
	{
		SRelement *elem = model.GetElement(e);
		int nfun = elem->GetNumFunctions();
		int neq = 3 * nfun;
		int len = neq*(neq + 1) / 2;
		MaxElementMem += (double) len;
	}
	MaxElementMem *= sizeof(double);
}

double* SRanalysis::ReadElementStiffness(SRelement* elem)
{
	//look up element stiffness matrix in memory
	//input:
		//elem = pointer to element
	//return:
		//start of element stiffness matrix for this element
	double *stiff = NULL;
	stiff = elem->GetStiffnessMatrix();
	return stiff;
}


bool SRanalysis::checkFlattenedElementHighStress()
{
	//check if an element whose mapping was flattened is adjacent to max- possible "error pollution"
	bool adjacentElemFlattened = false;

	double scaledFlattenTol = 0.1; // this is like flattening an edge that spans 45 degrees of a circle

	//rescale edge straightenfraction using edge chord length, then update elements flattenfraction
	for (int e = 0; e < model.GetNumEdges(); e++)
	{
		SRedge* edge = model.GetEdge(e);
		edge->ScaleStraightenVsEdgeLength();
	}
	for (int e = 0; e < model.GetNumElements(); e++)
	{
		SRelement* elem = model.GetElement(e);
		elem->rescaleFlattenFraction();
	}

	if (maxStressElid == -1)
		return false;
	SRelement* elemAtMax = model.GetElement(maxStressElid);
	double maxFlattenHighStress = 0.0;
	SRintVector badAdjElems;
	for (int n = 0; n < elemAtMax->GetNumNodesTotal(); n++)
	{
		int nid = elemAtMax->getNodeOrMidNodeId(n);
		int nf = model.GetNumNodeFaces(nid);
		for (int f = 0; f < nf; f++)
		{
			int fid = model.GetNodeFace(nid, f);
			SRface* face = model.GetFace(fid);
			for (int e = 0; e < 2; e++)
			{
				int elid = face->GetElementOwner(e);
				if (elid == -1 || elid == maxStressElid)
					continue;
				SRelement* adjElem = model.GetElement(elid);
				bool bf = false;
				for (int f = 0; f < adjElem->GetNumLocalFaces(); f++)
				{
					SRface* face = adjElem->GetFace(f);
					if (face->IsBoundaryFace())
					{
						bf = true;
						break;
					}
				}
				if (!bf)
					continue;

				if (adjElem->getFlattenFraction() > scaledFlattenTol)
				{
					badAdjElems.PushBack(elid);
					adjacentElemFlattened = true;
					if (adjElem->getFlattenFraction() > maxFlattenHighStress)
					{
						maxFlattenHighStress = adjElem->getFlattenFraction();
					}
				}
			}
		}
	}
	if (adjacentElemFlattened)
	{
		LOGPRINT(" max flattening adjacent to max stress element: %lg", maxFlattenHighStress);
		//ttd: comment out after dbg
		//plot faces of elements that were flattend more than scaledFlattenTol, but only if they have bdry faces:
		//post.PlotElems(badAdjElems.GetNum(), badAdjElems.d);
		return true;
	}
	return false;
}


double* SRanalysis::GetElementStiffnessVector()
{
	SRElementData* eld = model.GetElementData();
	return eld->getElementStiffness().GetVector();
}

void SRanalysis::SetEdgesToPorder(int p)
{
	//set all global edges of model to polynomial order p
	SRedge* edge;
	for (int i = 0; i < model.GetNumEdges(); i++)
	{
		edge = model.GetEdge(i);
		edge->putPorder(p);
	}
}

void SRanalysis::GetFileNameFromExtension(char* ext, SRstring& name)
{
	model.logFile.filename.Left('.', name);
	name += ".";
	name += ext;
}

SRnode* SRanalysis::GetNodeFromUid(int i)
{
	//look up a pointer to a global node using user id "i"
	int id = input.NodeFind(i);
	if (id == -1)
		return NULL;
	else
		return model.GetNode(id);
}

bool SRanalysis::SingStressCheck()
{
	//return true if distance from location of max stress to any sacrficial element is small compared to
	//size of element at max
	double distTol = BIG;
	if (maxStressElid == -1)
		return false;
	SRelement* elem = model.GetElement(maxStressElid);
	distTol = elem->GetSize();
	double sacrDistMin = BIG;
	for (int e = 0; e < model.GetNumElements(); e++)
	{
		SRelement* elem = model.GetElement(e);
		if (!elem->isSacrificial() || elem->isFlattened())
			continue;
		SRvec3 elcent;
		model.map.ElementCentroid(elem, elcent);
		double dist = maxStressPos.Distance(elcent);
		if (dist < sacrDistMin)
		{
			sacrDistMin = dist;
		}
	}
	if (sacrDistMin < distTol)
	{
		return true;
	}
	else
		return false;
}

void SRanalysis::AddSoftSprings()
{
	//add soft springs to each node for model stabilization.
	//approx model average stiffness: K = AE/L. A ~ L*L
	//L = model size. so K = LE. apply ~1e-08 times this to each element
	//face, so 1/3 that to each node
	double kmult = 1.e-08;
	double L = model.GetSize();
	double E = model.GetMaterial(0)->getE();
	double nodestiff = ONETHIRD*kmult*L*E;
	SRintVector nodeDone;
	int nnodes = model.GetNumNodes();
	nodeDone.Allocate(nnodes);
	int nel = model.GetNumElements();
	for (int elnum = 0; elnum < nel; elnum++)
	{
		SRelement* elem = model.GetElement(elnum);
		if (!elem->GetPChanged())
			continue;
		double *elStiff = ReadElementStiffness(elem);
		int nn = elem->GetNumCorners();
		for (int n = 0; n < nn; n++)
		{
			int nid = elem->GetNodeId(n);
			if (nodeDone.d[nid] == 1)
				continue;
			nodeDone.d[nid] = 1;
			int lfun = n;
			for (int dof = 0; dof < 3; dof++)
			{
				int eq = lfun * 3 + dof;
				int rowcolloc = elem->getStiffDiag(eq);
				elStiff[rowcolloc] += nodestiff;
			}
		}
	}
}

void SRanalysis::setAdaptLoopMax(int nits)
{
	//set max number of iterations allowed in adaptivity loop to "nits"
	adaptLoopMax = nits;
	if (nits > 3)
		maxPJump = 2;
}

double SRanalysis::CalculateMaxErrorForOutput()
{
	//calculate error estimate for output to report file
	//error criterion is min of:
	//error in the element at max stress in the model or
	//relative change in max stress in model between this and previous
	//adaptive pass

	SRelement* elemAtMax = NULL;
	double svmmax = 0.0;
	for (int e = 0; e < model.GetNumElements(); e++)
	{
		SRelement* elem = model.GetElement(e);
		if (!elem->isSacrificial())
		{
			if (elem->GetSvmMax() > svmmax)
			{
				svmmax = elem->GetSvmMax();
				elemAtMax = elem;
			}
		}
	}
	if (elemAtMax != NULL)
	{
		double err = elemAtMax->GetError();
		double err2;
		if (fabs(prevStressMax) > TINY)
		{
			err2 = fabs(stressMax - prevStressMax);
			if (fabs(stressMax) > TINY)
				err2 /= fabs(stressMax);
			err = MATHMIN(err, err2);
		}
		return 100.0*err;
	}
	else
		return 0.0;
}

void SRanalysis::setErrorMax(double error, int eluid)
{
	//set the maximum error in model
	//input:
		//error = error value in an element
		// eluid = element userid
	//note:
		//if error is greater than current max, sets class variables errorMax, errorSmoothRawAtMax, errorFaceJumpAtMax, and elUidAtMaxError
	if (error > errorMax)
	{
		errorMax = error;
		elUidAtMaxError = eluid;
	};
}

void SRanalysis::allocateSolutionVector()
{
	solution.Allocate(numEquations);
	model.setSolutionVector(solution.d);
}

void SRanalysis::PStats()
{
	//print p-order statistics to model output file
	LOGPRINT("P-order    percentage of model");
	int nedge = model.GetNumEdges();
	for (int p = 2; p <= maxPinModel; p++)
	{
		int edgesAtP = 0;
		for (int i = 0; i < nedge; i++)
		{
			if (model.GetEdge(i)->GetPorder() == p)
				edgesAtP++;
		}
		double percentage = 100.0* ((double)edgesAtP) / ((double)nedge);
		LOGPRINT("      %d    %4.2lf", p, percentage);
	}
}

void SRanalysis::zeroStressMax()
{
	stressMax = 0.0;
	for (int i = 0; i < 6; i++)
		stressMaxComp[i] = 0.0;
	maxsp1 = 0.0;
	minsp2 = BIG;
	maxCustom = 0.0;
};

int SRanalysis::GetNumFunctions()
{
	return model.GetNumFunctions();
}


