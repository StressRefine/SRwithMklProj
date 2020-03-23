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
// SRinput.cpp: implementation of the SRinput class.
//
//////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <search.h>
#include "SRanalysis.h"
#include "SRmachDep.h"
#include "SRinput.h"
#include "globalWrappers.h"

extern SRanalysis analysis;
extern SRmodel model;

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#endif


void SRinput::ReadModel()
{
	//read the model from input file in analysis.inputFile
	//also read any run settings specific to this model

	SRstring filename, line, basename, tail, tok;

	//assign output, results, and combined results names; delete existing
	//files in case future opens will be for "append":
	analysis.outdir.Right(slashChar, tail);
	basename = analysis.outdir;
	basename += slashStr;
	basename += tail;

	filename = basename;
	filename += ".frd";
	analysis.cgxFrdFile.setFileName(filename);

	filename = basename;
	filename += "_SR.f06";
	analysis.f06outFile.setFileName(filename);

	LOGPRINT("\nReading Settings...");
	readSettings();

	if (!analysis.inputFile.Open(SRinputMode))
	{
		const char *tmp = filename.LastChar(slashChar, true);
		LOGPRINT(" file not found: %s", tmp);
		ERROREXIT;
		return;
	}

	LOGPRINT("\nReading mesh...");

	/* read model: */
	nnode = nelem = ncon = nforce = nvolforce = nfacePressure = nfaceTraction = 0;
	nFaceCon = 0;
	nmat = ncoord = 0;
	numFaceFromNodalLoad = 0;
	anyFaceHasMultipleLoads = false;
	int numbricks = 0;
	int numwedges = 0;
	int numtets = 0;

	//create default GCS coordinate system (needed for lcs constraints)
	SRcoord* coord = model.addCoord();
	coord->Create(0.0, 0.0, 0.0);
	coord->setName("SR_DEFAULT_GCS");

	analysis.inputFile.GetLine(line);
	analysis.inputFile.ToTop();

	//count entities: nodes, elements, forces, local output nodes:
	while (1)
	{
		if (!analysis.inputFile.GetLine(line))
			break;

		if (line.isCommentOrBlank())
			continue;

		if (line == "nodes")
			CountEntities(nnode);
		else if (line == "elements")
			nelem = CountElements(numbricks, numwedges, numtets);
		else if (line == "forces")
			CountEntities(nforce);
		else if (line == "facePressures")
			CountEntities(nfacePressure);
		else if (line == "faceTractions")
			CountEntities(nfaceTraction);
		else if (line == "volumeForce")
			CountEntities(nvolforce);
		else if (line == "constraints")
			CountEntities(ncon);
		else if (line == "faceConstraints")
			CountEntities(nFaceCon);
		else if (line.CompareUseLength("coord"))
		{
			//input coordinates and materials on 1st pass because
			//they will be referred to by name:
			InputCoordinates();
		}
		else if (line == "materials")
			InputMaterials();
	}

	SetNumNodes(nnode);

	if (nnode == 0 || nelem == 0)
	{
		LOGPRINT("Mesh definition incomplete. Missing nodes or elements");
		ERROREXIT;
	}

	if (!userPSettings &&  nelem > 40000)
	{
		//auto econsolve:
		analysis.maxPorder = analysis.maxPorderFinalAdapt = 5;
		analysis.adaptLoopMax = 2;
	}

	LOGPRINT("Number of elements: %d", nelem);
	LOGPRINT("Number of nodes: %d", nnode);
	LOGPRINT("\n");

	SetNumElements(numbricks, numwedges, numtets);
	int neltmp = model.GetNumElements();

	model.allocateVolumeForces(nvolforce);
	//worst case of model with all bricks. Not wasteful because
	//these are just pointers
	int nedge = 12 * nelem;
	model.allocateEdges(nedge);
	nnode += nedge;

	nforce += nfacePressure;
	nforce += nfaceTraction;
	model.allocateForces(nforce);

	model.allocateNodes(nnode);
	analysis.inputFile.ToTop();

	while (1)
	{
		if (!analysis.inputFile.GetLine(line))
		{
			LOGPRINT("mesh input file incomplete");
			ERROREXIT;
		}

		if (line.isCommentOrBlank())
			continue;

		if (line == "nodes")
		{
			InputNodes();
			break;
		}
	}

	while (1)
	{
		if (!analysis.inputFile.GetLine(line))
			break;
		if (line.isCommentOrBlank())
			continue;

		if (line == "elements")
		{
			InputElements();
			break;
		}
	}
	allocateConstraints(ncon);

	int nel = model.GetNumElements();

	analysis.SetEdgesToPorder(2);

	createGlobalFaces();

	mapSetup();


	while (1)
	{
		if (!analysis.inputFile.GetLine(line))
			break;

		curLine = line;

		if (line.isCommentOrBlank())
			continue;
		if (line == "forces")
			InputNodalForces();
		else if (line == "facePressures")
			InputFacePressures();
		else if (line == "faceTractions")
			InputFaceTractions();
		else if (line == "volumeForce")
			InputVolumeForces();
		else if (line == "thermal")
			InputThermal();
		else if (line == "constraints")
			InputNodalConstraints();
		else if (line == "faceConstraints")
			InputFaceConstraints(nFaceCon);
	}

	analysis.inputFile.Close();

	model.freeNodeEdges();

	fillMultiFaceForceGroups();
	EnergySmoothNodalForces();
}

void SRinput::InputNodes()
{
	//input nodes

	//Input Spec:
		//userId x y z

	nodeUidOffset = 0;
	SRstring line, tok;
	SRnode* node;
	int nnode = 0, uid;
	double x, y, z, xmin, xmax, ymin, ymax, zmin, zmax;
	xmin = BIG;
	xmax = -BIG;
	ymin = BIG;
	ymax = -BIG;
	zmin = BIG;
	zmax= -BIG;
	while(1)
	{
		if(!analysis.inputFile.GetLine(line))
			break;
        if(line.isCommentOrBlank())
            continue;
		if(line == "end")
			break;
		line.TokRead(uid);
		if (nnode == 0)
			nodeUidOffset = uid;
		else if (nodeUidOffset != -1)
		{
			if (uid - nodeUidOffset != nnode)
				nodeUidOffset = -1;
		}
		line.TokRead(x);
		line.TokRead(y);
		line.TokRead(z);
		tok = line.Token();

		if (x < xmin)
			xmin = x;
		if(x > xmax)
			xmax = x;
		if(y < ymin)
			ymin = y;
		if(y > ymax)
			ymax = y;
		if(z < zmin)
			zmin = z;
		if(z > zmax)
			zmax = z;
		createNode(uid, x, y, z);

		if (uid > model.getMaxNodeUid())
			model.setMaxNodeUid(uid);
		nnode++;
	}
	double dx, dy, dz;
	dx = xmax - xmin;
	dy = ymax - ymin;
	dz = zmax - zmin;
	model.SetSize(sqrt(dx*dx + dy*dy + dz*dz));

	FillAndSortNodeUids();

	model.allocateNodeEdges(nnode);
}

void SRinput::InputElements()
{
	//Input Spec:
		//userId Material-name nodeuserId,i=1 to # of nodes
			//# of nodes=
			//10 for a tet, 15 for a wedge, 20 for a brick
			//for linear mesh # of nodes=
			//4 for a tet, 6 for a wedge, 8 for a brick
	//note:
		//this also creates the model edges vector

	SRstring line, tok;

	int uid, nodeuid, mid, nodev[20], nt, id;
	elemUidOffset = 0;
	while (1)
	{
		if(!analysis.inputFile.GetLine(line))
			break;
        if(line.isCommentOrBlank())
            continue;
		if(line == "end")
			break;
		line.TokRead(uid);
		if (nelem == 0)
			elemUidOffset = uid;
		else if (elemUidOffset != -1)
		{
			if (uid - elemUidOffset != nelem)
				elemUidOffset = -1;
		}

		if (uid > analysis.maxElemUid)
			analysis.maxElemUid = uid;

		tok = line.Token();
		mid = GetMaterialId(tok);
		SRmaterial* mat = model.GetMaterial(mid);
		mat->updateNumElements();
		nt = 0;
		while(1)
		{
			if(!line.TokRead(nodeuid))
				break;
			id = NodeFind(nodeuid);
			nodev[nt] = id;
			nt++;
		}
		if (nt == 4 || nt == 6 || nt == 8)
			ERROREXIT; //linear mesh not supported

		createElementIso(uid, nt, nodev, mat->getE(), mat->getNu());
	}

	bool anyorphan = false;
	SRnode* node;
	SRelement* elem;
	for (int i = 0; i < model.GetNumNodes(); i++)
	{
		node = model.GetNode(i);
		if (node->isOrphan())
		{
			anyorphan = true;
			break;
		}
	}
	if (anyorphan)
		LOGPRINT("Ignoring Nodes that are not Referenced by elements");

	if (elemUidOffset == -1)
	{
		//fill elem uid vector and sort in ascending order of uid for faster
		//elem-finding:
		elemUids.Allocate(nelem);
		SRuidData *euid;
		for (int i = 0; i < nelem; i++)
		{
			elem = model.GetElement(i);
			euid = elemUids.GetPointer(i);
			euid->id = i;
			euid->uid = elem->GetUserid();
		}
		SortElems();
	}
}

void SRinput::InputMaterials()
{
	//input materials

	//Input Spec:
		//name type
			//type = isotropic, orthotropic or general
			//if iso: rho, alpha then E,nu on second line
			//if ortho: rho, alphax, alphay, alphaz, then c11,c12,c13,c22,c23,c33,c44,c55,c66
			//if general: rho, alphax, alphay, alphaz, then full cij matrix, 36 constants,6 per line
				//(must be symmetric)

	SRstring line;
	SRmaterial* mat;
	SRstring type;

	SRfile* f;
	f = &analysis.inputFile;

	while(1)
	{
		if(!f->GetLine(line))
			break;
        if(line.isCommentOrBlank())
            continue;
		if(line == "end")
			break;
		if(model.GetNumMaterials() == MAXNMAT)
		{
			LOGPRINT("Maximum Number of Materials allowed is %d\n", MAXNMAT);
			ERROREXIT;
		}
		mat = model.addMat();
		mat->setId(model.GetNumMaterials() - 1);
		SRstring matname;
		matname = line.Token();
		mat->setName(matname);
		double E, nu, rho, tref, alphax, alphay, alphaz, allowableStress;
		SRcij orthoCij;
		SRgenAnisoCij gcij;
		type = line.Token();
		f->GetLine(line);
		line.TokRead(rho);
		if(type == "iso")
		{ 
			line.TokRead(alphax);
			line.TokRead(tref, CHECKFORTRAILINGCOMMENT);
			line.TokRead(allowableStress, CHECKFORTRAILINGCOMMENT);
			f->GetLine(line);
			line.TokRead(E);
			line.TokRead(nu);
			mat->assignTempRho(tref, alphax, alphax, alphax, rho, allowableStress);
			mat->assignIso(E, nu);
		}
		else if(type == "ortho")
		{
			line.TokRead(alphax);
			line.TokRead(alphay);
			line.TokRead(alphaz);
			line.TokRead(tref);
			line.TokRead(allowableStress);
			f->GetLine(line);
			line.TokRead(orthoCij.c11);
			line.TokRead(orthoCij.c12);
			line.TokRead(orthoCij.c13);
			line.TokRead(orthoCij.c22);
			line.TokRead(orthoCij.c23);
			line.TokRead(orthoCij.c33);
			line.TokRead(orthoCij.c44);
			line.TokRead(orthoCij.c55);
			line.TokRead(orthoCij.c66);
			mat->assignTempRho(tref, alphax, alphay, alphaz, rho, allowableStress);
			mat->assignOrtho(orthoCij);
		}
		else if(type == "gen")
		{
			line.TokRead(alphax);
			line.TokRead(alphay);
			line.TokRead(alphaz);
			line.TokRead(tref);
			line.TokRead(allowableStress);
			
			for(int i = 0; i < 6; i++)
			{
				f->GetLine(line);
				double cij;
				for (int j = 0; j < 6; j++)
				{
					line.TokRead(cij);
					gcij.setC(i, j, cij);
				}
			}
			if (!gcij.symCheck())
			{
				LOGPRINT("improper definition of general anisotropic material %s", matname.getStr());
				LOGPRINT("material is not symmetric");
				REPPRINT("improper definition of general anisotropic material %s", matname.getStr());
				REPPRINT("material is not symmetric");
				ERROREXIT;
			}
			mat->assignTempRho(tref, alphax, alphay, alphaz, rho, allowableStress);
			mat->assignGenAniso(gcij);
		}
		else
		{
			LOGPRINT("improper type in definition of material");
			ERROREXIT;
		}
	}

	equivMatTest();
}

void SRinput::InputCoordinates()
{
	//input local coordinate systems
	//Input Spec: name type "NotGcsAligned"
	//type = cartesian,spherical,cylindrical
	//x0,y0,z0 (origin)
	//if NotGcsAligned:
	//p1, p3 are points along local e1 and e3 axes

	SRstring line, type, tok;
	SRcoord* coord;
	int ncoord = 0;
	double x0, y0, z0;
	bool gcsaligned;

	while (1)
	{
		if (!analysis.inputFile.GetLine(line))
			break;
		if (line.isCommentOrBlank())
			continue;
		if (line == "end")
			break;
		if (ncoord > MAXNCOORD)
		{
			LOGPRINT("Maximum Number of Coordinate Systems allowed is %d\n", MAXNCOORD);
			LOGPRINT("Maximum Number of Coordinate Systems allowed is %d\n", MAXNCOORD);
			ERROREXIT;
		}
		coord = model.addCoord();
		coord->setName(line.Token());
		type = line.Token();
		tok = line.Token();
		if (tok == "NotGcsAligned")
			gcsaligned = false;
		else
			gcsaligned = true;
		analysis.inputFile.GetLine(line);
		line.TokRead(x0);
		line.TokRead(y0);
		line.TokRead(z0);
		if (type == "cartesian")
			coord->setType(cartesian);
		else if (type == "spherical")
			coord->setType(spherical);
		else if (type == "cylindrical")
			coord->setType(cylindrical);
		else
			ERROREXIT;

		if (gcsaligned)
			coord->Create(x0, y0, z0);
		else
		{
			SRvec3 p13, p3;
			analysis.inputFile.GetLine(line);
			line.TokRead(p13.d[0]);
			line.TokRead(p13.d[1]);
			line.TokRead(p13.d[2]);
			line.TokRead(p3.d[0]);
			line.TokRead(p3.d[1]);
			line.TokRead(p3.d[2]);
			coord->Create(x0, y0, z0, p13, p3);
		}
		ncoord++;
	}
}

void SRinput::InputNodalForces()
{
	//input nodal forces

	//input spec (all on one line):
	//node-Id
	//pressure "pressure" if pressure load (only valid for nodal loads that will become face loads) OR
	//or "coord" coord name
	//or "gcs"
	//if "pressure", 
	// pressure value
	// else
	//force vals for 3 dofs; 0 if not forced
	//notes:
	// edge and face forces will be automatically created if all nodes are loaded compatibly
	// and if doNodalToFace

	SRstring tok, line;

	SRforce* force;
	int nodeuid;
	SRforce forceTmp;

	int nnode = model.GetNumNodes();
	currentNodeForceStore.Allocate(nnode);

	currentForces.Allocate(nnode);
	if (nodeForceStore.GetNum() == 0)
		nodeForceStore.Allocate(nnode);

	bool rotateNodalLcs = true;

	while (1)
	{
		if (!analysis.inputFile.GetLine(line))
			break;
		if (line.isCommentOrBlank())
			continue;
		if (line == "end")
			break;
		line.TokRead(nodeuid);
		int nId = NodeFind(nodeuid);
		if (nId == -1)
			continue;//node doesn't exist in model

		//read in the force, then check for duplicate forces on this node:
		forceTmp.Clear();
		forceTmp.SetEntityId(nId);
		tok = line.Token();
		if (tok == "pressure")
		{
			forceTmp.setPressure(true);
			forceTmp.AllocateForceVals(1, 1);
			double p;
			line.TokRead(p);
			forceTmp.catForceVal(0, 0, p);
		}
		else if (tok == "coord")
		{
			tok = line.Token();
			int cId = GetCoordId(tok);
			forceTmp.setCoordId(cId);
		}
		else if (!tok.Compare("gcs"))
		{
			LOGPRINT("incorrect force type on node %d. must be 'pressure', 'coord Name' or 'gcs'", nodeuid);
			ERROREXIT;
		}

		if (!forceTmp.isPressure())
		{
			//for non pressure forces there are 3 force values. if the user omits the 2nd or 3rd they are set to 0
			forceTmp.AllocateForceVals(1, 3);
			double ft;
			for (int dof = 0; dof < 3; dof++)
			{
				if (!line.TokRead(ft))
					ft = 0.0;
				forceTmp.catForceVal(0, dof, ft);
			}
			if (forceTmp.GetCoordId() != -1 && rotateNodalLcs)
			{
				SRvec3 fl, fg;
				SRmat33 R;
				for (int dof = 0; dof < 3; dof++)
					fl.d[dof] = forceTmp.GetForceVal(0, dof);
				SRcoord* coord = model.GetCoord(forceTmp.GetCoordId());
				SRnode* node = model.GetNode(nId);
				coord->GetRotationMatrix(false, node->getPos(), R);

				for (int i = 0; i < 3; i++)
				{
					for (int j = 0; j < 3; j++)
						fg.d[i] += (R.rows[i].d[j] * fl.d[j]);
				}
				for (int dof = 0; dof < 3; dof++)
					forceTmp.setForceVal(0, dof, fg.d[dof]);
				forceTmp.setCoordId(-1);
			}
		}

		int nnodalForce = currentNodeForceStore.d[nId].forces.GetNum();

		//check for duplicate force, same node, could happen e.g. at a corner
		//if two edges sharing a face loaded.
		bool dupe = false;
		for (int iprev = 0; iprev < nnodalForce; iprev++)
		{
			int prevId = currentNodeForceStore.d[nId].GetForceId(iprev);
			force = currentForces.GetPointer(prevId);
			if (force->GetCoordId() == forceTmp.GetCoordId() && force->isPressure() == forceTmp.isPressure())
			{
				//duplicate and compatible force
				//add the force vectors
				dupe = true;
				force->AddNodalForce(forceTmp);
				break;
			}
		}
		if (!dupe)
		{
			//new force
			int forceid = currentForces.GetNum();
			currentNodeForceStore.d[nId].forces.PushBack(forceid);
			force = currentForces.Add();
			force->Copy(forceTmp);
		}
	}

	//AND currentForces and model.forces and update nodeForceStore:
	for (int i = 0; i < nnode; i++)
	{
		int nnodalForce = nodeForceStore.d[i].forces.GetNum();
		int curnumnodalForce = currentNodeForceStore.d[i].forces.GetNum();
		for (int curforcenum = 0; curforcenum < curnumnodalForce; curforcenum++)
		{
			int curforceid = currentNodeForceStore.d[i].GetForceId(curforcenum);
			SRforce* curforce = currentForces.d[curforceid];
			SRforce* force;
			bool dupe = false;
			for (int forcenum = 0; forcenum < nnodalForce; forcenum++)
			{
				force = nodeForceStore.d[i].GetForce(forcenum);
				if (force->GetCoordId() == curforce->GetCoordId() && force->isPressure() == curforce->isPressure())
				{
					force->AddNodalForce(*curforce, true); //summing sets is true
					dupe = true;
				}
			}
			if (!dupe)
			{
				//new force:
				int forceid = model.GetNumForces();
				force = model.addForce();
				nodeForceStore.d[i].forces.PushBack(forceid);
				force->Copy(*curforce);
			}
		}
	}

	//catch case that only midside nodes were loaded. This can occur for constant force on a face,
	//equiv nodal load at corners is 0 so femap doesn't put it out
	int numNewCornerForces = 0;
	for (int i = 0; i < nnode; i++)
	{
		int nnodalForce = nodeForceStore.d[i].forces.GetNum();
		SRnode* node = model.GetNode(i);
		if (nnodalForce == 0)
			continue;
		int eid = node->GetMidSideEdgeOwner();
		if (eid == -1)
			continue;
		SRedge* edge = model.GetEdge(eid);
		for (int n = 0; n < 2; n++)
		{
			int nid = edge->GetNodeId(n);
			if (nodeForceStore.d[nid].forces.GetNum() == 0)
				numNewCornerForces++;
		}
	}
	if (numNewCornerForces != 0)
	{
		int nforce = model.GetNumAllocatedForces();
		nforce += numNewCornerForces;
		nforce += nfacePressure;
		nforce += nfaceTraction;
		model.allocateForces(nforce);
		for (int i = 0; i < nnode; i++)
		{
			int nnodalForce = nodeForceStore.d[i].forces.GetNum();
			SRnode* node = model.GetNode(i);
			if (nnodalForce == 0)
				continue;
			int eid = node->GetMidSideEdgeOwner();
			if (eid == -1)
				continue;
			SRedge* edge = model.GetEdge(eid);
			for (int n = 0; n < 2; n++)
			{
				int nid = edge->GetNodeId(n);
				if (nodeForceStore.d[nid].forces.GetNum() == 0)
				{
					for (int f = 0; f < nnodalForce; f++)
					{
						force = nodeForceStore.d[i].GetForce(f);
						int forceid = model.GetNumForces();
						SRforce* cornerforce = model.addForce();
						cornerforce->Copy(*force);
						cornerforce->SetEntityId(nid);
						cornerforce->zeroForceVals();
						nodeForceStore.d[nid].forces.PushBack(forceid);
					}
				}
			}
		}
	}
	doNodalToFaceForces();
}

void SRinput::InputFacePressures()
{
	//input pressures on faces

	//input spec (all on one line):
	//eluid, n1, n2, n3, n4, p1, p2, p3, p4
	//eluid = uid of element that owns the face
	//n1, n2, n3, n4 = nodes at corner of face, n4 = 1 for tri;
	//p1, p2, p3, p4 = pressures at corner of face, p4 omitted for tri face

	SRstring tok, line;

	SRforce* force;
	int gno[4];
	int i;
	int nfacep = 0;
	while (1)
	{
		if (!analysis.inputFile.GetLine(line))
			break;
		if (line.isCommentOrBlank())
			continue;
		if (line == "end")
			break;
		int forceId = model.GetNumForces();
		force = model.addForce();
		force->setType(faceForce);
		force->setPressure(true);
		int eluid;
		line.TokRead(eluid);
		int elId = elemFind(eluid);
		SRelement* elem = model.GetElement(elId);
		int nv[4], nuid[4];
		for (i = 0; i < 4; i++)
		{
			line.TokRead(nuid[i]);
			if (nuid[i] != -1)
				nv[i] = NodeFind(nuid[i]);
			else
				nv[i] = -1;
		}
		int fId = -1;
		fId = model.elemFaceFind(elem, nv, gno);
		if (fId == -1)
		{
			LOGPRINT(" improperly defined pressure load on face. element user id: %d", eluid);
			LOGPRINT(" face node user ids: %d %d %d %d", nuid[0], nuid[1], nuid[2], nuid[3]);
			LOGPRINT(" face not found in the mesh");
			ERROREXIT;
		}

		force->SetEntityId(fId);
		SRface *face = model.GetFace(fId);
		face->pushBackForceids(forceId);
		force->setElemId(elId);
		double pv[8];
		int nn = 4;
		if (nv[3] == -1)
			nn = 3;
		for (i = 0; i < nn; i++)
			line.TokRead(pv[i]);

		force->AllocateForceVals(2 * nn, 1);
		for (i = 0; i < nn; i++)
		{
			int idg = gno[i];
			force->setForceVal(idg, 0, pv[i]);
		}
		//interpolate on face to get pv of midnodes:
		if (nn == 3)
		{
			//tri face
			pv[3] = 0.5*(force->GetForceVal(0, 0) + force->GetForceVal(1, 0));
			pv[4] = 0.5*(force->GetForceVal(1, 0) + force->GetForceVal(2, 0));
			pv[5] = 0.5*(force->GetForceVal(0, 0) + force->GetForceVal(2, 0));
		}
		else
		{
			//quad face
			pv[4] = 0.5*(force->GetForceVal(0, 0) + force->GetForceVal(1, 0));
			pv[5] = 0.5*(force->GetForceVal(1, 0) + force->GetForceVal(2, 0));
			pv[6] = 0.5*(force->GetForceVal(3, 0) + force->GetForceVal(2, 0));
			pv[7] = 0.5*(force->GetForceVal(0, 0) + force->GetForceVal(3, 0));
		}

		for (int l = 0; l < nn; l++)
		{
			int m = l + nn;
			//linearly interpolate from corners to get midside values:
			force->setForceVal(m, 0, pv[m]);
		}
	}
	nfacep++;
}


void SRinput::InputFaceTractions()
{
	//input tractions on faces

	//input spec (all on one line):
	//eluid, n1, n2, n3, n4
	//eluid = uid of element that owns the face
	//n1, n2, n3, n4 = uids of nodes at corner of face, n4 = -1 for tri
	//then 3 continuation lines of t1, t2, t3, t4
	// (tractions at corner of face, t4 omitted for tri face); 1 for each dof

	SRstring tok, line;

	SRforce* force = NULL;
	int gno[4];
	int i;

	while (1)
	{
		if (!analysis.inputFile.GetLine(line))
			break;
		if (line.isCommentOrBlank())
			continue;
		if (line == "end")
			break;
		int eluid;
		line.TokRead(eluid);
		int elId = elemFind(eluid);
		SRelement* elem = model.GetElement(elId);
		int nv[4], nuid[4];
		for (i = 0; i < 4; i++)
		{
			line.TokRead(nuid[i]);
			nv[i] = NodeFind(nuid[i]);
		}

		bool found = false;
		int fId = -1;
		fId = model.elemFaceFind(elem, nv, gno);
		if (fId == -1)
		{
			LOGPRINT(" improperly defined traction load on face. element user id: %d", eluid);
			LOGPRINT(" face node user ids: %d %d %d %d", nuid[0], nuid[1], nuid[2], nuid[3]);
			LOGPRINT(" face not found in the mesh");
			ERROREXIT;
		}

		//see if a force with fId exists:
		for (i = 0; i < model.GetNumForces(); i++)
		{
			force = model.GetForce(i);
			if (force->GetEntityId() == fId)
			{
				if (force->isPressure())
				{
					LOGPRINT(" duplicate pressure and traction load on face. element user id: %d", eluid);
					LOGPRINT(" face node user ids: %d %d %d %d", nuid[0], nuid[1], nuid[2], nuid[3]);
					ERROREXIT;
				}

				found = true;
				break;
			}
		}

		if (!found)
		{
			int forceId = model.GetNumForces();
			force = model.addForce();
			force->SetEntityId(fId);
			force->setType(faceForce);
			SRface* face = model.GetFace(fId);
			face->pushBackForceids(forceId);
		}
		force->setElemId(elId);

		force->setType(faceForce);
		for (int dof = 0; dof < 3; dof++)
		{

			if (!analysis.inputFile.GetLine(line) || !line.CompareUseLength("#"))
			{
				if (dof < 2)
				{
					LOGPRINT(" not enough degrees of freedom specified for traction load on face. element user id: %d", eluid);
					LOGPRINT(" face node user ids: %d %d %d %d", nuid[0], nuid[1], nuid[2], nuid[3]);
					ERROREXIT;
				}
			}
			tok = line.Token(); //skip the # character
			double tv[8];
			int nn = 4;
			if (nv[3] == -1)
				nn = 3;
			line.TokRead(tv[0]);
			for (i = 1; i < nn; i++)
				tv[i] = tv[0];
			for (i = 1; i < nn; i++)
			{
				double tt;
				if (!line.TokRead(tt))
					break;
				tv[i] = tt;
			}
			if (dof == 0)
				force->AllocateForceVals(2 * nn, 3);
			for (i = 0; i < nn; i++)
			{
				int idg = gno[i];
				force->setForceVal(idg, dof, tv[i]);
			}
			//interpolate on face to get tv of midnodes:
			if (nn == 3)
			{
				//tri face
				tv[3] = 0.5*(force->GetForceVal(0, dof) + force->GetForceVal(1, dof));
				tv[4] = 0.5*(force->GetForceVal(1, dof) + force->GetForceVal(2, dof));
				tv[5] = 0.5*(force->GetForceVal(0, dof) + force->GetForceVal(2, dof));
			}
			else
			{
				//quad face
				tv[4] = 0.5*(force->GetForceVal(0, dof) + force->GetForceVal(1, dof));
				tv[5] = 0.5*(force->GetForceVal(1, dof) + force->GetForceVal(2, dof));
				tv[6] = 0.5*(force->GetForceVal(3, dof) + force->GetForceVal(2, dof));
				tv[7] = 0.5*(force->GetForceVal(0, dof) + force->GetForceVal(3, dof));
			}

			for (int l = 0; l < nn; l++)
			{
				int m = l + nn;
				//linearly interpolate from corners to get midside values:
				force->setForceVal(m, dof, tv[m]);
			}
		} // for (int dof = 0; dof < 3; dof++)

	} // while (1)
}

void SRinput::InputVolumeForces()
{
	//input volume force (gravity or centrifugal)
	//input spec:
		//type (gravity,centrifugal)
		//g1,g2,g3 for gravity or
		//omega, axis, origin for centrifugal

	SRstring tok, line;

	SRvolumeForce* force;
	double omega;
	bool nextLineWasRead = false;
	while(1)
	{
		if (!nextLineWasRead)
		{
			if (!analysis.inputFile.GetLine(line))
				break;
		}
        if(line.isCommentOrBlank())
            continue;
		if (line == "end")
			break;
		force = model.addVolumeForce();
		if (line == "gravity")
		{
			force->setType(gravity);
			tok = line.Token();
			SRvec3 g;
			line.TokRead(g.d[0]);
			line.TokRead(g.d[1]);
			line.TokRead(g.d[2]);
			force->setG(g);
		}
		else if (line == "centrifugal")
		{
			force->setType(centrifugal);
			tok = line.Token();
			line.TokRead(omega);
			force->setOmega2(omega*omega);
			SRvec3 v;
			line.TokRead(v.d[0]);
			line.TokRead(v.d[1]);
			line.TokRead(v.d[2]);
			force->setAxis(v);
			line.TokRead(v.d[0]);
			line.TokRead(v.d[1]);
			line.TokRead(v.d[2]);
			force->setOrigin(v);
			double alpha = 0.0;
			line.TokRead(alpha);
			force->setAlpha(alpha);
		}
		else
		{
			LOGPRINT("improper volume force type. should be 'gravity' or 'centrifugal'");
			ERROREXIT;
		}
		if (!analysis.inputFile.GetLine(line))
			break;
		//optional continuation lines of element list
		if (line.CompareUseLength("#"))
		{
			nextLineWasRead = false;
			tok = line.Token(); //skip the # character
			tok = line.Token(); //skip the elements keyword
			int numel;
			line.TokRead(numel);
			force->allocateElList(numel);
			int nread = 0;
			int eluid, eid;
			while (nread < numel)
			{
				if (!analysis.inputFile.GetLine(line))
					break;
				tok = line.Token(); // skip the # character
				while (1)
				{
					if (!line.TokRead(eluid))
						break;
					eid = elemFind(eluid);
					force->putElList(nread, eid);
					nread++;
					if (nread == numel)
						break;
				}
			}
		}
		else
			nextLineWasRead = true;
	}
}

void SRinput::InputThermal()
{
	//input thermal loads
	//input spec:
		//"constant" temp
		//else "variable", then multiple lines of
		//nodeid, temperature, 1 for each node in model, one per line,
		//then "end"

	SRstring tok, line;

	model.allocateThermalForce();
	SRthermalForce* therm = model.GetThermalForce();
	double temp;
	analysis.inputFile.GetLine(line);
	int uid;
	if (line == "constant")
	{
		tok = line.Token();
		therm->setConstantTemp(true);
		line.TokRead(temp);
		therm->setTemp(temp);
	}
	else if (line == "variable")
	{
		therm->setConstantTemp(false);
		int n = model.GetNumNodes();
		therm->allocateNodalTemp(n);

		while(1)
		{
			analysis.inputFile.GetLine(line);
			if (line == "end")
				break;
			line.TokRead(uid);
			int i = NodeFind(uid);
			line.TokRead(temp);
			if (i != -1)
				therm->putNodalTemp(i, temp);
		}
	}
	else
		ERROREXIT;
}

void SRinput::InputNodalConstraints()
{
	//input nodal constraints

	//input spec: the following are all on one line:
	//node user id
	//enforced disp vals for each constrained dof, 0 if not enforced, "-" if not constrained
	//"coord" coord-name
	//notes:
		//edge and face constraints will be automatically created
		//if all nodes of the edge or face are constrained compatibly

	SRstring tok, line;
	SRconstraint* constraint;
	int dof;
	int nodeuid;

	SRconstraint conTmp;

	int nnode = model.GetNumNodes();

	currentNodeConStore.Allocate(nnode);

	currentConstraints.Allocate(nnode);

	if (nodeConStore.GetNum() == 0)
		nodeConStore.Allocate(nnode);


	while (1)
	{
		if (!analysis.inputFile.GetLine(line))
			break;
		if (line == "end")
			break;
		line.TokRead(nodeuid);
		int nId = NodeFind(nodeuid);
		if (nId == -1)
			continue;//node doesn't exist in model

		//read in the constraint, then check for duplicate constraints on this node:
		conTmp.Clear();
		conTmp.SetType(nodalCon);
		conTmp.SetEntityId(nId);
		double enfd;
		for (dof = 0; dof < 3; dof++)
		{
			conTmp.SetConstrainedDof(dof,false);
			tok = line.Token();
			if (tok.getLength() == 1 && tok == "-")
				continue;
			conTmp.SetConstrainedDof(dof, true);
			enfd = tok.RealRead();
			if (fabs(enfd) > TINY)
				analysis.anyEnforcedDisplacement = true;
			if (!conTmp.hasEnforcedDisp())
				conTmp.allocateEnforcedDisplacementData(1);
			conTmp.PutEnforcedDisplacementData(0, dof, enfd);
		}
		tok = line.Token();
		if (tok == "coord")
		{
			tok = line.Token();
			int id = GetCoordId(tok);
			conTmp.SetCoordId(id);
		}

		int nnodalCon = currentNodeConStore.d[nId].cons.GetNum();

		//check for duplicate constraint, same node, could happen e.g. at a corner
		//if two edges sharing a face constrained.
		bool dupe = false;
		for (int iprev = 0; iprev < nnodalCon; iprev++)
		{
			int prevId = currentNodeConStore.d[nId].GetConId(iprev);
			constraint = currentConstraints.GetPointer(prevId);
			if (constraint->GetCoordId() == conTmp.GetCoordId())
			{
				//duplicate and compatible constraint
				//add the constrainedDof flags and the enforced displacement vectors
				dupe = true;
				constraint->AddNodalConstraint(conTmp);
				break;
			}
		}
		if (!dupe)
		{
			//new constraint:
			int conid = currentConstraints.GetNum();
			currentNodeConStore.d[nId].cons.PushBack(conid);
			constraint = currentConstraints.Add();
			constraint->Copy(conTmp);
		}
	}

	//AND currentConstraints and model.constraints and update nodeConStore:
	for (int i = 0; i < nnode; i++)
	{
		int nnodalCon = nodeConStore.d[i].cons.GetNum();
		int curnumnodalCon = currentNodeConStore.d[i].cons.GetNum();
		for (int curconnum = 0; curconnum < curnumnodalCon; curconnum++)
		{
			int curconid = currentNodeConStore.d[i].GetConId(curconnum);
			SRconstraint* curcon = currentConstraints.GetPointer(curconid);
			SRconstraint* con;
			bool dupe = false;
			for (int connum = 0; connum < nnodalCon; connum++)
			{
				con = nodeConStore.d[i].GetCon(connum);
				if (con->GetCoordId() == curcon->GetCoordId())
				{
					con->AddNodalConstraint(*curcon, true); //summing sets is true
					dupe = true;
				}
			}
			if (!dupe)
			{
				//new constraint:
				int conid = model.GetNumConstraints();
				bool constrainedDof[3];
				SRvec3 enfd;
				for (int dof = 0; dof < 3; dof++)
					constrainedDof[dof] = curcon->IsConstrainedDof(dof);
				curcon->getDisp(0, enfd);
				inputNodalConstraint(curcon->GetEntityId(), constrainedDof, enfd.d);
				nodeConStore.d[i].cons.PushBack(conid);
				SRnode* node = model.GetNode(curcon->GetEntityId());
				node->SetConstraintId(conid);
			}
		}
	}

	currentNodeConStore.Free();
	currentConstraints.Free();

	doNodalToFaceConstraints();
}



void SRinput::InputFaceConstraints(int nfc)
{
	//input face constraints
	//input spec:
	//1st line: elem uid of the face owner, then nodeuids of the face corners (4th nodeuid =-1 for triangular face), then 3 constrained dof flags (0 or 1)
	//then optional "coord" coord-name. constraint is gcs is coord flag omitted
	//Then for each constrained dof, continuation lines of enforced disp vals for each corner (or just one value if constant)
	//note:
	//the continuation lines must start with "#"

	SRstring tok, line;

	if (nfc == 0)
		return;

	int ncon = model.GetNumConstraints();
	ncon += nfc;
	allocateConstraints(ncon);

	int gno[4];
	int i;
	bool condof[3];

	while (1)
	{
		if (!analysis.inputFile.GetLine(line))
			break;
		if (line == "end")
			break;
		int euid;
		line.TokRead(euid);
		int eid = elemFind(euid);
		if (eid == -1)
			ERROREXIT;
		SRelement *elem = model.GetElement(eid);
		int nv[4];
		nv[3] = -1;
		for (i = 0; i < 4; i++)
		{
			if (!line.TokRead(nv[i]))
				break;
			nv[i] = NodeFind(nv[i]);
		}
		int nnodes = nv[3] == -1 ? 3 : 4;
		if (nv[3] == -1)
			nnodes = 3;
		for (i = 0; i < 3; i++)
		{
			condof[i] = false;
			int ct;
			if (!line.TokRead(ct))
				break;
			if (ct == 1)
				condof[i] = true;
		}

		int coordId = -1;
		tok = line.Token();
		if (tok.CompareUseLength("coord"))
		{
			tok = line.Token();
			coordId = GetCoordId(tok);
		}

		int conId = model.GetNumConstraints();
		int fId = addFaceConstraint(eid, nv, condof, coordId);
		if (fId == -1)
			ERROREXIT;
		SRface* face = model.GetFace(fId);
		face->setConstraintId(conId);
		int nnodesTotal = face->GetNumNodesTotal();
		double enfd[8];
		double enfdm[8][3];
		bool allEndZero = true;
		bool enfdZero[3];
		for (int dof = 0; dof < 3; dof++)
		{
			enfdZero[dof] = true;
			if (!condof[dof])
				continue;
			if (!analysis.inputFile.GetLine(line))
				ERROREXIT;
			tok = line.Token();//skip continuation character ("#")
			int nread = 0;
			for (i = 0; i < nnodesTotal; i++)
			{
				if (!line.TokRead(enfd[i]))
					break;
				if (fabs(enfd[i]) > TINY)
				{
					enfdZero[dof] = false;
					allEndZero = false;
				}
				nread++;
			}
			if (nread == 1)
			{
				for (i = 1; i < nnodes; i++)
					enfd[i] = enfd[0];

			}
			else if (nread > nnodesTotal)
				ERROREXIT;
			if (enfdZero[dof])
				continue;

			if (nread < nnodesTotal)
			{
				//interpolate on face to get enfd of midnodes:
				if (nnodes == 3)
				{
					//tri face
					enfd[3] = 0.5*(enfdm[0][dof] + enfdm[1][dof]);
					enfd[4] = 0.5*(enfdm[1][dof] + enfdm[2][dof]);
					enfd[5] = 0.5*(enfdm[0][dof] + enfdm[2][dof]);
				}
				else
				{
					//quad face
					enfd[4] = 0.5*(enfdm[0][dof] + enfdm[1][dof]);
					enfd[5] = 0.5*(enfdm[1][dof] + enfdm[2][dof]);
					enfd[6] = 0.5*(enfdm[3][dof] + enfdm[2][dof]);
					enfd[7] = 0.5*(enfdm[0][dof] + enfdm[3][dof]);
				}
			}
			for (i = 0; i < nnodesTotal; i++)
				enfdm[i][dof] = enfd[i];
		}
		if (!allEndZero)
		{
			for (i = 0; i < nnodesTotal; i++)
			{
				SRvec3 enfdn;
				for (int dof = 0; dof < 3; dof++)
					enfdn.d[dof] = enfdm[i][dof];
				inputFaceNodeEnfd(conId, i, enfdn.d);
			}
		}
	}
}


void SRinput::readSettings()
{
	//check for run settings in SRsettings file
	//and set corresponding model settings flags
	SRstring filename, line, basename, tail, tok;
	bool econSolve = false;
	bool highAcc = false;

	SRfile& settingsFile = analysis.settingsFile;
	bool settingsOpened = settingsFile.Open(SRinputMode);
	if (!settingsOpened)
		return;

	LOGPRINT("Stress Refine\n");
	LOGPRINT("\nCustom Settings for this run:");
	bool anyCustom = false;
	userPSettings = false;
	customPmax = 0;
	numSettingsStrings = 0;

	bool lineWasRead = false;

	LOGPRINT("Stress Refine\n");
	LOGPRINT("\nCustom Settings for this run:");

	while (settingsOpened)
	{
		if (!lineWasRead)
		{
			if (!settingsFile.GetLine(line))
				break;
		}
		lineWasRead = false;
		if (line.isCommentOrBlank(SKIPCONTINATION))
			continue;
		else if (line.CompareUseLength("NOUNITS"))
			analysis.useUnits = false;
		else if (line.CompareUseLength("stress_conversion"))
		{
			line.Token();
			line.TokRead(analysis.stressUnitConversion);
			analysis.stressUnitstr = line.Token();
		}
		else if (line.CompareUseLength("length_conversion"))
		{
			line.Token();
			line.TokRead(analysis.lengthUnitConversion);
			analysis.lengthUnitstr = line.Token();
		}
		else if (line.CompareUseLength("needSoftSprings"))
		{
			analysis.needSoftSprings = true;
			LOGPRINT("use soft springs for stabilization");
		}
		else if (line == "PMAX")
		{
			tok = line.Token();
			int maxp;
			line.TokRead(maxp);
			analysis.setMaxPorder(maxp);
			LOGPRINT(" max p order %d\n", analysis.maxPorderFinalAdapt);
			anyCustom = true;
		}
		else if (line.CompareUseLength("MAXITERATIONS"))
		{
			tok = line.Token();
			int nits;
			line.TokRead(nits);
			analysis.setAdaptLoopMax(nits);
			LOGPRINT(" maximum adaptive loop iterations: %d\n", nits);
		}
		else if (line.CompareUseLength("maxPJump"))
		{
			tok = line.Token();
			line.TokRead(analysis.maxPJump);
		}
		else if (line.CompareUseLength("uniformP"))
		{
			analysis.pOrderUniform = true;
			LOGPRINT(" Uniform P-adaptivity\n");
			anyCustom = true;
		}
		else if (line.CompareUseLength("ERRORTOL"))
		{
			tok = line.Token();
			double errTol;
			line.TokRead(errTol);
			if (errTol < 0 || errTol > 25)
			{
				LOGPRINT(" error Tolerance must be in range 0 to 25 percent");
				ERROREXIT;
			}
			analysis.setErrorTolerance(errTol / 100.0);
			LOGPRINT(" Error Tolerance: %lg\n", analysis.GetErrorTolerance());
			anyCustom = true;
		}
		else if (line.CompareUseLength("NOSACRIFICIAL"))
		{
			LOGPRINT("Sacrificial Element detection disabled\n");
			anyCustom = true;
			analysis.SetDetectSacr(false);
		}
		else if (line.CompareUseLength("LOWSTRESSTOL"))
		{
			double tol;
			tok = line.Token();
			line.TokRead(tol);
			SetLowStressTolerance(tol);
			LOGPRINT("Low Stress Tolerance: %lg\n", tol);
			anyCustom = true;
		}
		else if (line.CompareUseLength("ALLCONSTRAINTSPENALTY"))
		{
			analysis.allConstraintsAsPenalty = true;
			LOGPRINT("calculate all constraints using penalty method\n");
			anyCustom = true;
		}
		else if (line.CompareUseLength("econSolve"))
		{
			econSolve = true;
			LOGPRINT(" economy solution setting");
			anyCustom = true;
		}
		else if (line.CompareUseLength("highAccuracy"))
		{
			highAcc = true;
			LOGPRINT(" high accuracy setting");
			anyCustom = true;
		}
		else if (line.CompareUseLength("noEnergySmooth"))
		{
			analysis.doEnergySmooth = false;
			LOGPRINT(" no energy smoothing of nodal loads");
			anyCustom = true;
		}
		else if (line.CompareUseLength("outputF06"))
		{
			analysis.outputf06 = true;
			LOGPRINT(" Output Ascii Nastran Stresses (.f06)");
			anyCustom = true;
		}
		else
			LOGPRINT("unrecognized setting: %s", line.getStr());
	}
	if (!anyCustom)
		LOGPRINT("   --None");

	settingsFile.Close();

	if (econSolve)
	{
		analysis.maxPorder = analysis.maxPorderFinalAdapt = 5;
		analysis.adaptLoopMax = 2;
	}
	if (highAcc)
	{
		analysis.ErrorTolerance = 0.03;
		analysis.maxPorder = analysis.maxPorderFinalAdapt = 8;
		SetLowStressToleranceFinalAdapt(0.5);
		SetLowStressTolerance(0.5);
	}
}

void SRinput::doNodalToFaceForces()
{
	//convert nodal forces to face forces if all the nodes of a face are compatibly loaded
	int nfaceForce = nodalToFaceForces(countOnlyFlag);
	if (nfaceForce != 0)
	{
		int nforce0 = model.GetNumForces();
		int nforceTotal = nforce0 + nfaceForce;
		model.allocateForces(nforceTotal);
		nodalToFaceForces(addTheForcesFlag);//false to add the face forces
	}

	if (nfaceForce == 0)
		return;

	//free the nodal forces which becaem inactive:
	for (int f = 0; f < model.GetNumForces(); f++)
	{
		SRforce* force = model.GetForce(f);
		if (force->GetType() == inactiveForce)
			model.freeForce(f);
	}
	model.packForces();
	//face tractions and face pressures may not have been read yet. reallocate for them:
	model.allocateForces(model.GetNumForces() + nfacePressure + nfaceTraction);

	//correct force ids of face forces
	//(may have changed because of removal of inactive nodal forces):
	for (int i = 0; i < model.GetNumForces(); i++)
	{
		SRforce* f = model.GetForce(i);
		if (f == NULL)
			continue;
		int eid = f->GetEntityId();
		if (f->GetType() == faceForce)
			model.GetFace(eid)->pushBackForceids(i);
	}

	//error check. there should be no "orphan" nodal forces with pressure loads, or orphan midside nodes
	for (int f = 0; f < model.GetNumForces(); f++)
	{
		SRforce* force = model.GetForce(f);
		bool okforce = true;
		if ((force->GetType() == nodalForce))
		{
			if (force->isPressure())
			{
				okforce = false;
				break;
			}
			int nid = force->GetEntityId();
			SRnode* node = model.GetNode(nid);
			if (node->isMidSide())
			{
				okforce = false;
				break;
			}
		}
		if (!okforce)
		{
			SRnode* node = model.GetNode(force->GetEntityId());
			LOGPRINT(" incorrect nodal force definition for node %d", node->GetUserid());
			LOGPRINT(" pressure loads may not be applied to individual nodes, only to all nodes of faces");
			ERROREXIT;
		}
	}

	nodeForceStore.Free();

}

void SRinput::doNodalToFaceConstraints()
{
	//convert nodal constraints to face constraints if all the nodes of a face are compatibly constrained
	if (model.GetNumFaces() == 0)
		return;
	//create face constraints from nodal constraints:
	int nfaceCon = nodalToFaceConstraints(countOnlyFlag);
	if (nfaceCon == 0)
		return;
	//reallocate model.constraints with room for the face constraints, then add the face constraints.
	int ncon0 = model.GetNumConstraints();
	int nconTotal = ncon0 + nfaceCon;
	allocateConstraints(nconTotal);
	nodalToFaceConstraints(addTheConstraintsFlag);

	//fatal error if any nodes remain with more than one active constraint;
	for (int i = 0; i < nodeConStore.GetNum(); i++)
	{
		int numactive = 0;
		int n = nodeConStore.d[i].cons.GetNum();
		for (int j = 0; j < n; j++)
		{
			SRconstraint* con = nodeConStore.Get(i).GetCon(j);
			if (con->GetType() == nodalCon)
				numactive++;
		}
		if (numactive > 1)
			ERROREXIT;
	}

	//get rid of any remaining constraints that aren't sensible:
	//nodal constraint on orphan node, or constrained midside node whose edge is not constrained.
	for (int i = 0; i < model.GetNumConstraints(); i++)
	{
		SRconstraint* con = model.GetConstraint(i);
		if (con->GetType() == nodalCon)
		{
			int nodeId = con->GetEntityId();
			SRnode* node = model.GetNode(nodeId);
			if (node->isOrphan())
				con->SetType(inactiveCon);
			else if (node->isMidSide())
			{
				int eid = model.GetNode(nodeId)->GetMidSideEdgeOwner();
				SRedge* edge = model.GetEdge(eid);
				if (edge->GetNode(0)->getConstraintId() == -1 || edge->GetNode(1)->getConstraintId() == -1)
					con->SetType(inactiveCon);
			}
		}
	}

	//delete the inactive nodal constraints:
	for (int i = 0; i < model.GetNumConstraints(); i++)
	{
		SRconstraint* con = model.GetConstraint(i);
		if (con->GetType() == inactiveCon)
		{
			int nodeId = con->GetEntityId();
			SRnode* node = model.GetNode(nodeId);
			node->SetConstraintId(-1);
			model.freeConstraints(i);
		}
	}

	model.packConstraints();

	//correct constraint ids of entities involved in constraints
	//(may have changed because of removal of inactive nodal constraints):
	for (int i = 0; i < model.GetNumConstraints(); i++)
	{
		SRconstraint* con = model.GetConstraint(i);
		if (con == NULL)
			continue;
		int eid = con->GetEntityId();
		if (con->GetType() == nodalCon)
			model.GetNode(eid)->SetConstraintId(i);
		if (con->GetType() == faceCon)
			model.GetFace(eid)->setConstraintId(i);
		if (analysis.allConstraintsAsPenalty && con->isGcs())
			con->SetCoordId(0); //0 is reserved for default gcs coord system
	}

	nodeConStore.Free();

}



