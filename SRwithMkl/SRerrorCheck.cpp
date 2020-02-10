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
// SRerrorCheck.cpp: implementation of the SRerrorCheck class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "SRmodel.h"
#include "SRanalysis.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#endif

extern SRmodel model;
extern SRanalysis analysis;

SRerrorCheck::SRerrorCheck()
{
	smallMaxStressDetected = false;
	lowStressTol = LOWSTRESSTOL; 
	lowStressTolFinalAdapt = LOWSTRESSTOLFINAL;
	checkReentrant = true;
	numSacrificialElements = 0;
}

void SRerrorCheck::Initialize()
{
	lowStressTol = LOWSTRESSTOL;
	lowStressTolFinalAdapt = LOWSTRESSTOLFINAL;
}

double SRerrorCheck::FindError(SRelement *elem)
{
	//look up current error in an element
	//input:
		//elem = pointer to element
		//basis = scratch space for basis functions
	//return:
		//error
	//note:
		//both smooth vs raw stresses at the elements gauss points,
		//and traction jump errors at the element's faces are checked
		//this routine also updates the max error in model, unless the element is sacrificial
	//the error calculated is also stored as a class variable in the element

	double errorSmoothRaw = 0.0, errorFaceJumps = 0.0;

	errorSmoothRaw = FindSmoothedVsRawError(elem);

    errorFaceJumps = FindFaceTractionJumpError(elem);

	double error = MATHMAX(errorSmoothRaw, errorFaceJumps);

	if (!elem->isSacrificial())
	{
		//04/15 errorSmoothRaw can be too conservative. use face jump error for reporting, less conservative
		error = errorFaceJumps;
		analysis.setErrorMax(error, elem->GetUserid(), errorSmoothRaw, errorFaceJumps);
	}

	elem->PutError(error);
    return error;
}

bool SRerrorCheck::FindRequiredPOrder(SRelement *elem, int &p)
{
	//Estimate required P-order for an element based on current error
	//input:
		//elem = pointer to element
	//output:
		//p = P-order
	//return:
		//true if update required else false
	//note:
		//p is not allowed to exceed 2 for sacricial elements
		//if the element is sacricial
		//if the element has low stress compared to the allowable for the material or max in model,
		//the p order is not allowed to exceed the max p order for low stress elements for the model

	int currentP;
	double error, errRatio, expon, xp;
	error = FindError(elem);
	double errTol = analysis.GetErrorTolerance();
	if (error < errTol)
		return false;
	currentP = elem->GetMaxPorder();

	//from Zienkiewicz-Zhu, '89:
	expon = 1.0 / currentP;
	errRatio = error / errTol;
	xp = currentP*pow(errRatio, expon);
	p = (int)xp;
	double dp = xp - (double)p;
	if (dp > 0.5)
		p++;
	if (p <= currentP)
		p = currentP + 1;
	if (p - currentP > analysis.getMaxPJump())
		p = currentP + analysis.getMaxPJump();
	if (elem->isSacrificial())
	{
		if (p > 2)
			p = 2;
	}

	double stressTol = GetStressMaxForErrorCheck(elem);
	double relTol = LOWSTRESSTOL;
	int pLow = analysis.GetMaxPorderLowStress();

	if (elem->GetMaxStress() < relTol*stressTol)
	{
		if (p > pLow)
			p = pLow;
		if (currentP == pLow)
			return false;
	}

	if (p > analysis.GetMaxPorder())
		p = analysis.GetMaxPorder();

	if (currentP == analysis.GetMaxPorder())
		return false;
	else
		return true;
}

double SRerrorCheck::FindSmoothedVsRawError(SRelement *elem)
{
	//determine error in element by comparing globally-smoothed vs raw strains at each gauss point
	//return:
		//max error in element of any strain component at any gauss point

	int i, j, nint;
	double w, r, s, t, err, errmax = 0.0;
	double strain[6], stress[6], svm;
	nint = model.math.FillGaussPoints(elem);
	double strainScale = GetStrainMaxForErrorCheck(elem);

	elem->SetBasisData();

	double rawsvmmax = 0.0;
	double rawstrain[6];
	for (i = 0; i < nint; i++)
	{
		model.math.GetGP3d(i,r,s,t,w);
		elem->GetSmoothedStrain(r, s, t, strain);
		for(j = 0; j < 6; j++)
		{
			rawstrain[j] = elem->GetRawStrain(i, j);
			err = fabs(strain[j] - rawstrain[j]) / strainScale;
			if (err > errmax)
				errmax = err;
		}
		elem->StraintoStress(r, s, t, rawstrain, stress);
		svm = model.math.GetSvm(stress);
		if (svm > rawsvmmax)
			rawsvmmax = svm;
	}
	return errmax;
}

double SRerrorCheck::FindFaceTractionJumpError(SRelement *elem)
{
	//find the maximum traction jump on any face of an element
	//input:
		//elem = pointer to element
	//note:
		//before this routine is called, traction jumps must be stored for each global face
		//by calling fillFaceTractionJumps

	double err,errmax = 0.0;
    SRface* face;
	double stressScale = GetStressMaxForErrorCheck(elem);
	for (int lface = 0; lface < elem->GetNumLocalFaces(); lface++)
	{
		face = elem->GetFace(lface);
        err = face->GetTractionJump();
		if (stressScale > TINY)
			err /= stressScale;
		if(err > errmax)
			errmax = err;
	}
	return errmax;
}

void SRerrorCheck::CleanUp()
{
	//miscellaneous clean up after error checking

	int i;
	SRelement* elem;
	for(i = 0; i < model.GetNumElements(); i++)
	{
		elem = model.GetElement(i);
		elem->FreeRawStrains();
	}
}

bool SRerrorCheck::AutoSacrificialElements()
{
	//mark elements as sacrificial that have:
	//1. nodal forces or constraints
	//2. edge constraints
	//3. edges at reentrant corners with adjacent boundary faces not fully constrained
	//4. edges that share two boundary faces, one is constrained, the other not constrained or incompatibibly constrained
	//   (except for cartesian symmetry constraints)

	SRintVector singNode;
	int nnode = model.GetNumNodes();
	singNode.Allocate(nnode);
	SRintVector singEdge;
	singEdge.Allocate(model.GetNumEdges());
	bool anySingNode = false;
	bool anySingEdge = false;

	int i, id, n;
	for (i = 0; i < model.GetNumNodes(); i++)
	{
		SRnode* node = model.GetNode(i);
		if (node->isOrphan())
			continue;
	}
	n = model.GetNumForces();
	for (i = 0; i < n; i++)
	{
		SRforce* force = model.GetForce(i);
		if (force->GetType() == nodalForce)
		{
			id = force->GetEntityId();
			singNode.Put(id, 1);
			anySingNode = true;
		}
		else if (force->GetType() == edgeForce)
		{
			id = force->GetEntityId();
			singEdge.Put(id, 1);
			anySingEdge = true;
		}
	}
	n = model.GetNumConstraints();
	for (i = 0; i < n; i++)
	{
		SRconstraint* con = model.GetConstraint(i);
		if (con->GetType() == nodalCon)
		{
			id = con->GetEntityId();
			singNode.Put(id, 1);
			anySingNode = true;
		}
	}

	//check for reentrant corners faces on surface; singular unless
	//all 3 dofs of both faces contrained
	int nej = model.GetNumEdges();
	SRintVector edgeface1;
	SRintVector edgeface2;
	SRintVector edgefacelej1;
	SRintVector edgefacelej2;
	edgeface1.Allocate(nej);
	edgeface2.Allocate(nej);
	edgefacelej1.Allocate(nej);
	edgefacelej2.Allocate(nej);
	for (i = 0; i < nej; i++)
	{
		edgeface1.Put(i, -1);
		edgeface2.Put(i, -1);
	}

	SRface* face;
	n = model.GetNumFaces();
	for (i = 0; i < n; i++)
	{
		face = model.GetFace(i);
		if (!face->IsBoundaryFace())
			continue;
		for (int lej = 0; lej < face->GetNumLocalEdges(); lej++)
		{
			int gej = face->GetLocalEdgeGlobalId(lej);
			if (edgeface1.Get(gej) == -1)
			{
				edgeface1.Put(gej, i);
				edgefacelej1.Put(gej, lej);
			}
			else
			{
				edgeface2.Put(gej, i);
				edgefacelej2.Put(gej, lej);
			}
		}
	}

	SRface* face2;
	double dotmin = 1.0;
	bool AnyExoneratedByLowStressJump = false;

	for (i = 0; i < nej; i++)
	{
		SRedge* edge = model.GetEdge(i);

		int f1 = edgeface1.Get(i);
		int f2 = edgeface2.Get(i);
		if (f1 == -1 || f2 == -1)
			continue;
		face = model.GetFace(f1);
		face2 = model.GetFace(f2);
		int lej = edgefacelej1.Get(i);
		int lej2 = edgefacelej2.Get(i);

		bool compatibleConstraints = true;
		bool nonSymCon = false;
		bool nonSymCon2 = false;
		bool bothAllCon = true;
		int numCondof = 0;
		int numCondof2 = 0;
		int connedDof = -1;
		int connedDof2 = -1;
		//corner between two faces on boundary
		SRconstraint* con = face->GetConstraint();
		SRconstraint* con2 = face2->GetConstraint();
		//ok to check for symmetry constraint only if faces are normal to each other:
		bool facesAreNormal = checkFacesNormal(face, lej, face2, lej2);

		if (con != NULL)
		{
			for (int dof = 0; dof < 3; dof++)
			{
				if (!con->IsConstrainedDof(dof))
					bothAllCon = false;
				else
				{
					numCondof++;
					if (connedDof == -1)
						connedDof = dof;
					else
					{
						//not a symmetry constraint if more then one constrained dof
						nonSymCon = true;
					}
				}
			}
			if (!nonSymCon && facesAreNormal)
				nonSymCon = !checkForSymCon(face, con, connedDof);
		}
		else
		{
			bothAllCon = false;
		}
		if (con2 != NULL)
		{
			if (con == NULL)
				compatibleConstraints = false;
			for (int dof = 0; dof < 3; dof++)
			{
				if (!con2->IsConstrainedDof(dof))
				{
					bothAllCon = false;
					if (con != NULL && con->IsConstrainedDof(dof))
						compatibleConstraints = false;
				}
				else
				{
					numCondof2++;
					if (con != NULL && !con->IsConstrainedDof(dof))
						compatibleConstraints = false;
					if (connedDof2 == -1)
						connedDof2 = dof;
					else
					{
						//not a symmetry constraint if more then one constrained dof
						nonSymCon2 = true;
					}
				}
			}
			if (!nonSymCon2 && facesAreNormal)
				nonSymCon2 = !checkForSymCon(face2, con2, connedDof2);
		}
		else
		{
			if (con != NULL)
				compatibleConstraints = false;

			bothAllCon = false;
		}
		if (!compatibleConstraints && (nonSymCon || nonSymCon2))
		{
			//last check: nonsingular if at least one face is flat,
			//can slide in plane, and other is constrained or enforced normal to it:
			if (facesAreNormal && (numCondof < 2) && (numCondof2 < 2))
			{
				if (face->isFlat())
					compatibleConstraints = checkForTangentialCon(face, con, face2, con2, connedDof, connedDof2);
				else if (face2->isFlat())
					compatibleConstraints = checkForTangentialCon(face2, con2, face, con, connedDof2, connedDof);
			}
		}

		if (!compatibleConstraints && (nonSymCon || nonSymCon2))
		{
			//constrained boundary face next to incompatibibly constrained boundary face.
			singEdge.Put(i, 1);
			anySingEdge = true;
		}
		else if (!bothAllCon && checkReentrant)
		{
			//both faces not fully constrained,
			//check for reentrant corner:

			SRvec3 norm, norm2;
			SRvec3 pos1, pos2, v21;
			double r, s;
			face->NaturalCoordinatesNearMidedge(lej, r, s);
			face->Position(r, s, pos1);
			face->OutwardNormal(r, s, norm);
			face2->NaturalCoordinatesNearMidedge(lej2, r, s);
			face2->Position(r, s, pos2);
			face2->OutwardNormal(r, s, norm2);
			pos1.Subtract(pos2, v21);
			v21.Normalize();
			double normdot = norm.Dot(norm2);
			if (normdot < 0.866)//30 degrees
			{
				//kink in the normals to the faces at the corners
				if ((norm2.Dot(v21) > 0.0) && (norm.Dot(v21) < 0.0))
				{
					//reentrant if vector from point on face2 to point on face 1
					//points in same direction as normal to face2 and in opposite direction to normal to face1
					//LOGPRINT(" sing edge face %d norm dot %lg\n", face->GetId(), dot);
					anySingEdge = true;
					if (normdot < dotmin)
					{
						dotmin = normdot;
						SRelement* elem;
						int eid = face->GetElementOwner(0);
						elem = model.GetElement(eid);
						eid = face2->GetElementOwner(0);
						elem = model.GetElement(eid);
					}
					if (normdot < -0.1)
					{
						//this is suspicious. could be crack-like feature, but could also be spurious, e.g. outward normal sign
						//fooled by slivery element.
						singEdge.Put(i, 1);
					}
					else
						singEdge.Put(i, 1);
				}
			}
		}
	}

	if (!anySingEdge && !anySingNode)
		return false;

	//mark nodes owned by singular edges as singular:
	for (i = 0; i < model.GetNumEdges(); i++)
	{
		int pf = singEdge.Get(i);
		if (pf != 0)
		{
			SRedge* edge = model.GetEdge(i);
			for (int j = 0; j < 2; j++)
			{
				id = edge->GetNodeId(j);
				singNode.Put(id, pf);
			}
		}
	}

	//mark all elements that own a singular node as sacrificial:
	for (i = 0; i < model.GetNumElements(); i++)
	{
		SRelement* elem = model.GetElement(i);
		for (int j = 0; j < elem->GetNumNodes(); j++)
		{
			id = elem->GetNodeId(j);
			int pf = singNode.Get(id);
			if (pf > 0)
			{
				elem->SetSacrificial();
				break;
			}
		}
	}

	//any faceFromNodal forces that slipped through conversion to multifaceforces, and did not
	//get smoothed. set their edges to p2:
	for (int i = 0; i < model.GetNumForces(); i++)
	{
		SRforce* force = model.GetForce(i);
		if (force->isFaceFromNodal())
		{
			int fid = force->GetEntityId();
			SRface* face = model.GetFace(fid);
			for (int l = 0; l < face->GetNumLocalEdges(); l++)
			{
				SRedge* edge = face->GetEdge(l);
				edge->putPorder(2);
				edge->setSacrificial();
			}
		}
	}

	//adaptivity won't work unless their are at least two
	//non-sacrificial elements or face jumps can be checked
	//at interfaces with sacrificial elements:
	int numnonsacr = 0;
	int numsacr = 0;
	for (i = 0; i < model.GetNumElements(); i++)
	{
		SRelement* elem = model.GetElement(i);
		if (!elem->isSacrificial())
			numnonsacr++;
		else
			numsacr++;
	}

	numSacrificialElements = numsacr;

	if (numsacr > 0)
	{
		OUTPRINT("Elements detected adjacent to singularities");
		OUTPRINT("These elements will be ignored for error checking and for computing max stress in model");
		OUTPRINT("The stress is theoretically infinite at these singularities.\nPlease review your model geometry, loads, and constraints to verify that this was your intent");
	}
	if (numnonsacr < 2)
	{
		OUTPRINT("Warning: there are less than two nonsingular elements in model. \n");
		OUTPRINT("          There may not be sufficient information for error checking\n");
		OUTPRINT("          and adaptivity\n\n");
	}
	else
		OUTPRINT("\n");

	for (i = 0; i < model.GetNumElements(); i++)
	{
		SRelement* elem = model.GetElement(i);
		if (elem->isSacrificial())
		{
			for (int l = 0; l < elem->GetNumLocalEdges(); l++)
				elem->GetEdge(l)->setSacrificial();
		}
	}

	return (numsacr != 0);
}

double SRerrorCheck::fillFaceTractionJumps()
{
//store traction jump for each global face in model
	//return:
		//max thermal strain on any face
	//note:
		//fills variable tractionJump for each global face in model
		//max traction jump is calculated for any of the p2 (3x3) gauss points for the face

	int i, j, k, nint;
	double rf, sf, w, tract1, tract2, r, s, t;
	double tract[3];
    double eT, eTMax = 0.0;
	SRelement* elem1 = NULL;
	SRelement* elem2 = NULL;
	SRface* face;
	SRvec3 norm;
    SRdoubleMatrix faceTract1, faceTract2;

	bool p2Gps = true; //set this to false to sample gauss points based on current p order of face
	int nintMax = 0;
	int nintMaxP2 = 9;
	int facePoverride = -1;
	for (i = 0; i < model.GetNumFaces(); i++)
	{
		face = model.GetFace(i);
		nint = model.math.CountGaussPoints(face);
		if (nint > nintMax)
			nintMax = nint;
		if (p2Gps && nintMax > nintMaxP2)
		{
			nintMax = nintMaxP2;
			facePoverride = 2;
		}
	}
    faceTract1.Allocate(nintMax, 3);
    faceTract2.Allocate(nintMax, 3);

    for(i = 0; i < model.GetNumFaces(); i++)
    {
        face = model.GetFace(i);
        face->SetTractionJump(0.0);

        int id1 = face->GetElementOwner(0);
        int id2 = face->GetElementOwner(1);
		int lface1 = face->GetElementLocalFace(0);
		int lface2 = face->GetElementLocalFace(1);

		if (id2 != -1)
		{
			elem2 = model.GetElement(id2);
			//skip sacrficicial element so it does not drive error of adjacent non-sacrificial element:
			if (elem2->isSacrificial())
				continue;
		}

        //skip boundary face if constrained:
        if(id2 == -1 && (face->GetConstraintId() != -1) )
            continue;

		nint = model.math.FillGaussPoints(face, facePoverride);

        double tractJump = 0.0;

	    elem1 = model.GetElement(id1);

		//skip sacrficicial element so it does not drive error of adjacent non-sacrificial element:
		if (elem1->isSacrificial())
			continue;


		for (j = 0; j < nint; j++)
		{
			model.math.GetGP2d(j, rf, sf, w);
			eT = getElementFaceTraction(lface1, face, elem1, rf, sf, r, s, t, tract);
            if(eT > eTMax)
                eTMax = eT;
			for(k = 0; k < 3; k++)
                faceTract1.Put(j, k, tract[k]);
		}

        if(id2 == -1)
        {
            //face owned by only one element. check traction vs applied traction (or vs zero if no loading):
			double fv[3], rf, sf, w;
			if (!face->hasForce())
            {
				for(k = 0; k < 3; k++)
                    fv[k] = 0.0;
                        
            }
			for(j = 0; j < nint; j++)
			{
				if(face->hasForce())
				{
					model.math.GetGP2d(j, rf, sf, w);
					face->GetSummedForceValue(rf, sf, fv);
				}
				for(k = 0; k < 3; k++)
                {
					tract1 = faceTract1.Get(j, k);
					tract2 = fv[k];
                    double jumpK = fabs(tract1 - tract2);
                    if(jumpK > tractJump)
                        tractJump = jumpK;
				}
			}
        }
        else
        {
			//compare traction from first element that owns the face with second:
	        elem2 = model.GetElement(id2);

			if (elem2->isSacrificial())
				continue;

		    for(j = 0; j < nint; j++)
		    {
				model.math.GetGP2d(j, rf, sf, w);
				eT = getElementFaceTraction(lface2, face, elem2, rf, sf, r, s, t, tract);
                if(eT > eTMax)
                    eTMax = eT;
			    for(k = 0; k < 3; k++)
                    faceTract2.Put(j, k, tract[k]);
		    }

		    for(j = 0; j < nint; j++)
		    {
				for(k = 0; k < 3; k++)
                {
					tract1 = faceTract1.Get(j, k);
					tract2 = faceTract2.Get(j, k);
                    double jumpK = fabs(tract1 - tract2);
                    if(jumpK > tractJump)
                        tractJump = jumpK;
				}
            }
        }

		face->SetTractionJump(tractJump);
    }

    return eTMax;
}


double SRerrorCheck::getElementFaceTraction(int lface, SRface *face, SRelement* elem, double rf, double sf, double& r, double& s, double& t, double* tract)
{
	//calculate traction at a face at natural coordinates rf, sf
	//input:
		//lface = local face number
		//face = pointer to global face
		//elem = pointer to global element of which face is a local face
		//rf, sf = natural coordinates on face
	//output:
		//tract = traction vector
		//r, s, t = corresponding element natural coordinates to rf, sf
	//note:
		//a warning is printed if the position at rf,sf of face is not the same
		//as position at r,s,t of element. That would indicae a bug in model.map.ElementNaturalCoordsFromFace

	double eT, stress[6], strain[6], etx, ety, etz;
    SRvec3 norm;

	SRvec3 facePos;
	face->Position(rf, sf, facePos);
	model.map.ElementNaturalCoordsFromFace(elem, lface, rf, sf, r, s, t);

	double detj = elem->FillMapping(r, s, t);
	elem->FillBasisFuncs(r, s, t, derivonly);
	eT = elem->CalculateRawStrain(r, s, t, strain, etx, ety, etz);
	elem->StraintoStress(r, s, t, strain, stress);

	for (int i = 0; i < 6; i++)
	{
		double ae = fabs(strain[i]);
		if (ae > strainMax)
		{
			strainMax = ae;
			elem->checkMaterialStrainMax(ae);
		}
		ae = fabs(stress[i]);
		if (ae > stressMax)
		{
			stressMax = ae;
			elem->checkMaterialStressMax(ae);
		}
	}

	//note: see definition of OutwardNormal. It uses the 1st element that
	//owns a face to determine outward normal. This assures that consistent normal
	//is used for tractions on both sides
	face->OutwardNormal(rf, sf, norm);
	model.math.GetTraction(norm.d, stress, tract);
	return eT;
}

double SRerrorCheck::GetStressMaxForErrorCheck(SRelement* elem)
{
	//calculate max stress to use for comparing stress error against in an element
	//input:
		//elem = pointer to element
	//return:
		//max stress
	//note:
		//if small max stress in model was detected, return stressMax
		//else
		//if all active materials have an allowable, and this element's material allowable is less than max in model,
		//max stress in this element's material is returned,
		//else max in model is returned.
	return stressMax;
}

double SRerrorCheck::GetStrainMaxForErrorCheck(SRelement* elem)
{
	//calculate max strain to use for comparing strain error against in an element
	//input:
		//elem = pointer to element
	//return:
		//max strain
	//note:
		//the default is max strain in model
		//if the material has been assigned an allowable stress, the material allowable strain is set to corresponding strain value.
		//and that strain value is used if it is less than max in model

	SRmaterial* elMat = elem->GetMaterial();
	if (!model.getAllMatsHaveAllowable())
		return strainMax;
	double allowStress = elMat->GetAllowableStress();
	if (allowStress > TINY && allowStress < model.getMaxAllowableAnyActiveMat())
		return elMat->GetMaxStrain();
	else
		return strainMax;
}


bool SRerrorCheck::checkForTangentialCon(SRface* face, SRconstraint* con, SRface* face2, SRconstraint* con2, int connedDof, int connedDof2)
{
	//check if constraints are non-singular because: 
		//1st face is flat
		// only direction orthogonal to it is constrained
		// other face is constrained in a direction tangential to 1st face
	//input:
		//face = the 1st face
		//con = its constraint
		//connedDof = its constrained Dof
		//face2 = the 2nd face
		//con2 = its constraint
		//connedDof2 = its constrained Dof
	//true if situation is non sinfular

	if (!face->isFlat())
		return false;
	if (con == NULL)
		return false;

	if (con2 == NULL)
		return checkForOrthogonalCon(face, con, connedDof);

	//check if constrained dof is normal to face:
	double rc, sc;
	model.map.FaceCentroid(face, rc, sc);
	SRvec3 norm;
	face->OutwardNormal(rc, sc, norm);
	SRvec3 dofDirGcs;
	SRvec3 dofDirGcs2;
	//in gcs, unit vector in direction connedDof is corresponding row of identity matrix:
	model.math.getIdentityMatrix().vecCopy(connedDof, dofDirGcs, true);//copy rows flag is true
	if (con2 != NULL)
		model.math.getIdentityMatrix().vecCopy(connedDof2, dofDirGcs2, true);//copy rows flag is true
	SRvec3 dofDir;
	if (!con->isGcs() || !con2->isGcs())
	{
		int nsamp = 1;
		SRvec3 pos[5];
		model.map.FaceCentroid(face, pos[0]);
		//also sample at nodes:
		nsamp += face->GetNumNodes();
		for (int i = 1; i < nsamp; i++)
			pos[i].Copy(face->GetNode(i - 1)->Position());

		SRmat33 R;
		for (int i = 0; i < nsamp; i++)
		{
			//is connedDof from this face normal to it:
			dofDir.Copy(dofDirGcs);
			if (!con->isGcs())
			{
				//rotate unit vector in direction connedDof to line up with lcs:
				con->GetCoord()->GetRotationMatrix(true, pos[i], R);//true for global to local, dummyArg for position because position not needed for cartesian system
				dofDir.Rotate(R);
			}
			double dot = dofDir.Dot(norm);
			if (fabs(dot) < (1.0 - RELSMALL))
				return false;
			//is connedDof from other face normal to this face:
			dofDir.Copy(dofDirGcs2);
			if (!con2->isGcs())
			{
				//rotate unit vector in direction connedDof2 to line up with lcs:
				con2->GetCoord()->GetRotationMatrix(true, pos[i], R);//true for global to local, dummyArg for position because position not needed for cartesian system
				dofDir.Rotate(R);
			}
			dot = dofDir.Dot(norm);
			if (fabs(dot) > RELSMALL)
				return false;
		}
	}
	else
	{
		//is connedDof from this face normal to it:
		double dot = dofDirGcs.Dot(norm);
		if (fabs(dot) < (1.0 - RELSMALL))
			return false;
		//is connedDof from other face normal to this face:
		dot = dofDirGcs2.Dot(norm);
		if (fabs(dot) > RELSMALL)
			return false;
	}
	return true;
}

bool SRerrorCheck::checkForSymCon(SRface* face, SRconstraint* con, int connedDof)
{
	//check if a constraint is a symmetry constraint
	//input:
		//face = global face the constraint is on
		//con = constraint
		//connedDof = the single dof is constrained
	//return:
		//true if con is a symmetry constraint
	//note:
		//criteria for symmetry constraint:
		//gcs or cartesian lcs
		//no enforced displacement, face is flat
		//the single dof that is constrained, connedDof, points in directionnormal to face
	if (con->hasEnforcedDisp())
		return false;
	return checkForOrthogonalCon(face, con, connedDof);
}

bool SRerrorCheck::checkForOrthogonalCon(SRface* face, SRconstraint* con, int connedDof)
{
	//check if a constraint is a orthogonal to a flat face
	//input:
		//face = global face the constraint is on
		//con = constraint
		//connedDof = the single dof is constrained
	//return:
		//true if con is a symmetry constraint
	//note:
		//criteria for orthogonal constraint:
			//gcs or cartesian lcs
			//face is flat
			//the single dof that is constrained, connedDof, points in directionnormal to face
	if (!face->isFlat())
		return false;
	//check if constrained dof is normal to face:
	double rc, sc;
	model.map.FaceCentroid(face, rc, sc);
	SRvec3 norm;
	face->OutwardNormal(rc, sc, norm);
	SRvec3 dofDirGcs;
	//in gcs, unit vector in direction connedDof is corresponding row of identity matrix:
	model.math.getIdentityMatrix().vecCopy(connedDof, dofDirGcs, true);//copy rows flag is true
	SRvec3 dofDir;
	if (!con->isGcs())
	{
		int nsamp = 1;
		SRvec3 pos[5];
		model.map.FaceCentroid(face, pos[0]);
		if (con->GetCoord()->GetType() != cartesian)
		{
			//also sample at nodes:
			nsamp += face->GetNumNodes();
			for (int i = 1; i < nsamp; i++)
				pos[i].Copy(face->GetNode(i - 1)->Position());
		}

		SRmat33 R;
		for (int i = 0; i < nsamp; i++)
		{
			dofDir.Copy(dofDirGcs);
			//rotate unit vector in direction connedDof to line up with lcs:
			con->GetCoord()->GetRotationMatrix(true, pos[i], R);//true for global to local, dummyArg for position because position not needed for cartesian system
			dofDir.Rotate(R);
			double dot = dofDir.Dot(norm);
			if (fabs(dot) < (1.0 - RELSMALL))
				return false;
		}
	}
	else
	{
		double dot = dofDirGcs.Dot(norm);
		if (fabs(dot) < (1.0 - RELSMALL))
			return false;
	}
	return true;
}

bool SRerrorCheck::checkFacesNormal(SRface *face, int lej, SRface* face2, int lej2)
{
	//check if two faces are normal to each other at the center of a shared edge
	//input:
		//face = 1st face
		//lej = local edge on 1st face corresponding to the shared edge
		//face2 = 2nd face
		//lej2 = local edge on 2nd face corresponding to the shared edge
	//return:
		//true if the faces are normal else false
	double r, s;
	SRvec3 norm, norm2;
	face->NaturalCoordinatesNearMidedge(lej, r, s);
	face->OutwardNormal(r, s, norm);
	face2->NaturalCoordinatesNearMidedge(lej2, r, s);
	face2->OutwardNormal(r, s, norm2);
	double dot = norm2.Dot(norm);
	if (fabs(dot) < RELSMALL)
		return true;
	else
		return false;
}

double SRerrorCheck::getElSvmAtMidNode(SRelement* elem, int mid)
{
	//calculate von Mises stress at a midnode belong to one of an
	//element's edges
	//input:
		//elem = pointer to element
		//mid = global node id of mid-node
	//return:
		//von Mises stress at midnode
	double r, s, t, etx, ety, etz;
	double strain[6], stress[6];
	for (int lej = 0; lej < elem->GetNumLocalEdges(); lej++)
	{
		if (elem->GetEdge(lej)->GetMidNodeId() == mid)
		{
			model.map.ElementNaturalCoordsAtMidedge(elem, lej, r, s, t);
			elem->CalculateRawStrain(r, s, t, strain, etx, ety, etz);
			elem->StraintoStress(r, s, t, strain, stress);
			return model.math.GetSvm(stress);
		}
	}
	return 0.0;
}

void SRerrorCheck::SetUp(bool finalAdapt)
{
	//miscellaneous set up before performing error checking:
	//notes:
		//calculate raw strains and perform smoothing
		//store traction jumps at all faces

	SRvec3 pos;
	double eTmax = 0.0;

	if (finalAdapt)
		lowStressTol = lowStressTolFinalAdapt;

	smallMaxStressDetected = false;
	for (int e = 0; e < model.GetNumElements(); e++)
	{
		SRelement* elem = model.GetElement(e);
		elem->DownloadDisplacement();
	}

	strainMax = analysis.post.GlobalStrainSmooth();

	analysis.post.CalculateMaxStress();

	stressMax = analysis.GetStressMax();

	eTmax = fillFaceTractionJumps();

	double minMatScale = BIG;
	for (int m = 0; m < model.GetNumMaterials(); m++)
	{
		double scale = model.GetMaterial(m)->MatScale();
		if (scale < minMatScale)
			minMatScale = scale;
	}

	//catch free thermal expansion case:
	if (strainMax < RELSMALL*eTmax)
	{
		strainMax = eTmax;
		stressMax = minMatScale*strainMax;
	}

	//very low stress (rigid body) case:
	//if 0.2% is yield strain, this is 0.000001 times that:
	double smallStrainTol = 2.e-08;

	if (stressMax < minMatScale*smallStrainTol)
	{
		smallMaxStressDetected = true;
		if (strainMax < smallStrainTol)
			strainMax = smallStrainTol;
		stressMax = minMatScale*smallStrainTol;
	}
}


