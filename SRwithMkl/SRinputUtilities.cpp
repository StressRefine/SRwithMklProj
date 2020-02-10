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
// SRinputUtilities.cpp: implementation of utility routines for the SRinput class.
//
//////////////////////////////////////////////////////////////////////



#include "stdafx.h"
#include <stdlib.h>
#include <search.h>
#include "SRmodel.h"
#include "SRanalysis.h"
#include "SRmachDep.h"
#include "SRinput.h"

extern SRmodel model;
extern SRanalysis analysis;

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#endif

int SRuidDataCompareFunc(const void *v1, const void *v2)
{
	//compare function for sorting and searching of node userids
	//input:
		//v1 = pointer to uid data for node 1
		//v2 = pointer to uid data for node 2
	//return:
		//-1 if node1 uid is less than node2 uid
		//0  if node1 uid is equal to node2 uid
		//+1 if node1 uid is greater than node2 uid

	SRuidData* node1 = (SRuidData*)v1;
	SRuidData* node2 = (SRuidData*) v2;
	if (node1->uid < node2->uid)
		return -1;
	else if (node1->uid == node2->uid)
		return 0;
	else
		return 1;
}

void SRinput::Cleanup()
{
	nodeUids.Free();
	elemUids.Free();
}


int SRinput::NodeFind(int uid)
{
	//find node with user Id uid
	//input:
		//uid = user id to match
	//return:
		//number of the node that matches uid, -1 if not found

	if (uid == -1)
		return -1;
	else if (nodeUidOffset != -1)
	{
		int nid = uid - nodeUidOffset;
		//make sure nid is in bounds:
		if (nid < 0 || nid >= model.GetNumNodes())
			return -1;
		return nid;
	}
	else
	{
		//binary search:
		SRuidData* nuid;
		SRuidData uidt;
		//"search key" has to be same data type expected by compare function see SRuidDataCompareFunc:
		uidt.uid = uid;
		nuid = (SRuidData *)bsearch(&uidt, nodeUids.GetVector(), nodeUids.GetNum(), sizeof(SRuidData), SRuidDataCompareFunc);
		if (nuid == NULL)
			return -1;
		else
			return nuid->id;
	}
}


int SRinput::elemFind(int uid)
{
	//find elem with user Id uid
	//input:
		//uid = user id to match
	//return:
		//number of the node that matches uid, -1 if not found

	//binary search:
	if (elemUidOffset != -1)
		return uid - elemUidOffset;
	else
	{
		SRuidData* euid;
		SRuidData uidt;
		//"search key" has to be same data type expected by compare function see SRuidDataCompareFunc:
		uidt.uid = uid;
		euid = (SRuidData *)bsearch(&uidt, elemUids.GetVector(), elemUids.GetNum(), sizeof(SRuidData), SRuidDataCompareFunc);
		if (euid == NULL)
			return -1;
		else
			return euid->id;
	}
}

void SRinput::CountEntities(int &num)
{	
	//count the entities currently being read in input file by counting 
	//lines (except blank or comment) until end is encountered

	SRstring line;
	while (1)
	{
		if(!analysis.inputFile.GetLine(line))
			break;
		if (line.isCommentOrBlank(SKIPCONTINATION))
            continue;
		if(line=="end")
			break;
		num++;
	}
}

void SRinput::CountElements(int &num)
{
	//count the elements currently being read in input file by counting 
	//lines (except blank or comment) until end is encountered
	//note:
		//this is the same as CountEntities but also checks for linear mesh

	SRstring line, tok;
	while (1)
	{
		if (!analysis.inputFile.GetLine(line))
			break;
		if (num == 0)
		{
			tok = line.Token();//skip uid field
			tok = line.Token();//skip material name field
			int nt = 0, nuid;
			while (1)
			{
				if (!line.TokRead(nuid))
					break;
				nt++;
			}
			if (nt == 4 || nt == 6 || nt == 8)
				ERROREXIT; //linear mesh not allowed
		}
		if (line.isCommentOrBlank(SKIPCONTINATION))
			continue;
		if (line == "end")
			break;
		num++;
	}
}

int SRinput::GetMaterialId(SRstring &name)
{
	//look up the material id with "name"
	//return:
		//material number, -1 if not found

	SRmaterial* mat;
	for (int i = 0; i < model.GetNumMaterials(); i++)
	{
		mat = model.GetMaterial(i);
		if (name.Compare(mat->GetName()))
			return i;
	}
	return -1;
}

int SRinput::GetCoordId(SRstring &name)
{
	//look up the coordinate system id with "name"
	//return:
		//coordinate system number, -1 if not found

	SRcoord*     coord;
	for (int i = 0; i < model.GetNumCoords(); i++)
	{
		coord = model.GetCoord(i);
		if (name == coord->GetName())
			return i;
	}
	SRASSERT(false);
	return -1;
}
void SRinput::SortNodes()
{
	//sort "nodeUids" array in ascending order of uid
	int n = nodeUids.GetNum();
	qsort(nodeUids.GetVector(), n, sizeof(SRuidData), SRuidDataCompareFunc);
}

void SRinput::SortElems()
{
	//sort "elemUids" array in ascending order of uid
	int n = elemUids.GetNum();
	qsort(elemUids.GetVector(), n, sizeof(SRuidData), SRuidDataCompareFunc);
}

int SRinput::nodalToFaceConstraints(bool countOnly)
{
	//convert nodal to face constraints if all nodes of the face are constrained
	//input:
		//countOnly = true to just count the new face constraints, not add them
	//notes:
		//adds new face constraints to model.constraints
		//use of nodeconstore: duplicate nodal constraints are blended together if they are compatible (bot gcs or both lcs with same coord sys)
		//there will only be duplicate constraints in nodeconstore if they are not compatible.
		//So it's necessary to use nodeconstore to loop over all the duplicate constraints of each node

	int nface = model.GetNumFaces();

	int nfaceCon = 0;
	for (int f = 0; f < nface; f++)
	{
		SRface* face = model.GetFace(f);
		if (face->GetElementOwner(1) != -1)//not boundary face
			continue;

		int nn = face->GetNumNodes();

		bool allNodesConstrained = true;
		bool anyNodeConstrained = false;
		for (int n = 0; n < nn; n++)
		{
			int nid = face->GetNodeId(n);
			int ncon = nodeConStore.Get(nid).cons.GetNum();
			if (ncon == 0)
			{
				allNodesConstrained = false;
				break;
			}
		}

		if (allNodesConstrained)
		{
			//check midside nodes:
			for (int l = 0; l < nn; l++)
			{

				int mid = face->getLocalEdge(l)->GetEdge()->GetMidNodeId();
				if (nodeConStore.Get(mid).cons.GetNum() == 0)
				{
					allNodesConstrained = false;
					break;
				}
			}
		}

		if (!allNodesConstrained)
			continue;

		//convert to face constraint only if dofs constrained are compatible and coordIds are compatible:
		bool constrainedDof0[3];
		int nnodes = nn;
		int nnodesTotal = 2 * nnodes;
		int nid = face->GetNodeId(0);
		int ncon0 = nodeConStore.Get(nid).cons.GetNum();
		SRconstraint* newFaceCon = NULL;

		for (int icon0 = 0; icon0 < ncon0; icon0++)
		{
			bool compatible = true;
			SRconstraint* nodeCon0 = nodeConStore.Get(nid).GetCon(icon0);
			int coord0 = nodeCon0->GetCoordId();
			for (int dof = 0; dof < 3; dof++)
				constrainedDof0[dof] = nodeCon0->IsConstrainedDof(dof);
			bool allNodesConstrainedDofSame[3];
			bool anyNodeEnforcedDof[3];
			for (int dof = 0; dof < 3; dof++)
			{
				allNodesConstrainedDofSame[dof] = true;
				if (fabs(nodeCon0->getDisp(0, dof)) > TINY)
					anyNodeEnforcedDof[dof] = true;
				else
					anyNodeEnforcedDof[dof] = false;
			}
			for (int n = 1; n < nnodesTotal; n++)
			{
				int nid = face->GetNodeOrMidNodeId(n);
				int ncon = nodeConStore.Get(nid).cons.GetNum();
				for (int icon = 0; icon < ncon; icon++)
				{
					compatible = true;
					SRconstraint* nodeCon = nodeConStore.Get(nid).GetCon(icon);
					int coord = nodeCon->GetCoordId();
					if (coord != coord0)
					{
						compatible = false;
						continue;
					}
					for (int dof = 0; dof < 3; dof++)
					{
						if (!constrainedDof0[dof] || (constrainedDof0[dof] && !nodeCon->IsConstrainedDof(dof)))
							allNodesConstrainedDofSame[dof] = false;
						if (fabs(nodeCon->getDisp(0, dof)) > TINY)
							anyNodeEnforcedDof[dof] = true;
					}
					if (!compatible)
						continue;
					bool anydofsame = false;
					for (int dof = 0; dof < 3; dof++)
					{
						if (allNodesConstrainedDofSame[dof])
						{
							anydofsame = true;
							break;
						}
					}

					if (!anydofsame)
						compatible = false;
					if (compatible)
						break;
				}  //endif: for (icon = 0; icon < node->GetNumConstraints(); icon++)
				if (!compatible)
					break;

			} // endif: for (int n = 1; n < nnodesTotal; n++)

			if (compatible)
			{
				//can't have more than one constraint on same face:
				if (newFaceCon != NULL)
					ERROREXIT;
				nfaceCon++;
				if (countOnly)
					break;
				face->setConstraintId(model.GetNumConstraints());
				newFaceCon = model.addConstraint();
				newFaceCon->SetEntityId(f);
				for (int dof = 0; dof < 3; dof++)
				{
					if (allNodesConstrainedDofSame[dof])
						newFaceCon->SetConstrainedDof(dof, 1);
				}
				newFaceCon->SetType(faceCon);
				newFaceCon->SetCoordId(coord0);
				bool anyEnforcedDisp = false;
				for (int dof = 0; dof < 3; dof++)
				{
					if (anyNodeEnforcedDof[dof])
						anyEnforcedDisp = true;
				}

				if (anyEnforcedDisp)
				{
					newFaceCon->allocateEnforcedDisplacementData(nnodesTotal);
					for (int dof = 0; dof < 3; dof++)
					{
						if (!newFaceCon->IsConstrainedDof(dof))
							continue;
						for (int n = 0; n < nnodesTotal; n++)
						{
							if (anyNodeEnforcedDof[dof])
							{
								double enfd = nodeCon0->getDisp(0, dof);
								newFaceCon->PutEnforcedDisplacementData(n, dof, enfd);
							}
						}
					}
				}
			} //endif: if (compatible)
		} // endif: for (int icon0 = 0; icon0 < node->GetNumConstraints(); icon0++)
	} //endif: for (int f = 0; f < nface; f++)

	if (!countOnly)
	{
		//turn off constrained dofs of nodes involved in face constraints:
		for (int f = 0; f < nface; f++)
		{
			SRface* face = model.GetFace(f);
			if (face->GetElementOwner(1) != -1)//not boundary face
				continue;
			SRconstraint* faceCon = face->GetConstraint();
			if (faceCon == NULL)
				continue;

			int nnodesTotal = 2 * face->GetNumNodes();;

			for (int n = 0; n < nnodesTotal; n++)
			{
				int nid = face->GetNodeOrMidNodeId(n);
				int ncon = nodeConStore.Get(nid).cons.GetNum();
				for (int icon = 0; icon < ncon; icon++)
				{
					SRconstraint* nodeCon = nodeConStore.Get(nid).GetCon(icon);
					for (int dof = 0; dof < 3; dof++)
					{
						if (faceCon->IsConstrainedDof(dof))
							nodeCon->SetConstrainedDof(dof, 0);
					}
				}
			}
		}

		for (int n = 0; n < model.GetNumConstraints(); n++)
		{
			SRconstraint* nodeCon = model.GetConstraint(n);
			if (nodeCon->GetType() == nodalCon)
			{
				bool anydofconstrained = false;
				for (int dof = 0; dof < 3; dof++)
				{
					if (nodeCon->IsConstrainedDof(dof))
					{
						anydofconstrained = true;
						break;
					}
				}
				if (!anydofconstrained)
					nodeCon->SetType(inactiveCon);
			}
		}
	}

	return nfaceCon;
}

int SRinput::nodalToFaceForces(bool countOnly)
{
	//convert nodal forces to face forces if all the nodes of a face are loaded compatibly
	//and it is a boundary face

	//input:
		//countonly = true to just count the number of new face forces else false
	//note:
		//adds new face forces to model.forces

	int nface = model.GetNumFaces();

	SRintVector compatibleForceNum;
	compatibleForceNum.Allocate(8);
	compatibleForceNum.Set(-1);

	int numFaceForce = 0;
	for (int f = 0; f < nface; f++)
	{
		SRface* face = model.GetFace(f);
		for (int i = 0; i < face->GetNumLocalEdges(); i++)
			int mid = face->getLocalEdge(i)->GetEdge()->GetMidNodeId();
		if (face->GetElementOwner(1) != -1) //not boundary face
			continue;

		int nn = face->GetNumNodes();

		bool allNodesLoaded = true;
		for (int n = 0; n < nn; n++)
		{
			int nid = face->GetNodeId(n);
			int nforce = nodeForceStore.Get(nid).forces.GetNum();
			if (nforce == 0)
			{
				allNodesLoaded = false;
				break;
			}
		}
		if (allNodesLoaded)
		{
			//check also for midside nodes
			for (int l = 0; l < nn; l++)
			{

				int mid = face->getLocalEdge(l)->GetEdge()->GetMidNodeId();
				int nforce = nodeForceStore.Get(mid).forces.GetNum();
				if (nforce == 0)
				{
					allNodesLoaded = false;
					break;
				}
			}

		}
		if (!allNodesLoaded)
			continue;

		//convert to face force only if coordIds and pressure flags are compatible:
		int nnodes = nn;
		int nnodesTotal = 2 * nnodes;
		int nid0 = face->GetNodeId(0);
		int nforce0 = nodeForceStore.Get(nid0).forces.GetNum();
		SRforce* newFaceForce = NULL;
		for (int iforce0 = 0; iforce0 < nforce0; iforce0++)
		{
			compatibleForceNum.Set(-1);
			compatibleForceNum.Put(0, iforce0);
			SRforce* nodeForce0 = nodeForceStore.Get(nid0).GetForce(iforce0);
			int coord0 = nodeForce0->GetCoordId();
			bool pressure0 = nodeForce0->isPressure();
			double fv0[3];
			int ndof = 3;
			if (pressure0)
				ndof = 1;
			for (int dof = 0; dof < ndof; dof++)
				fv0[dof] = nodeForce0->GetForceVal(0, dof);

			SRforce* nodeForce;

			for (int n = 1; n < nnodesTotal; n++)
			{
				int nid = face->GetNodeOrMidNodeId(n);
				int nforce = nodeForceStore.Get(nid).forces.GetNum();
				for (int iforce = 0; iforce < nforce; iforce++)
				{
					SRforce* nodeForce = nodeForceStore.Get(nid).GetForce(iforce);
					int coord = nodeForce->GetCoordId();
					if ((coord == coord0) && (nodeForce->isPressure() == pressure0))
					{
						compatibleForceNum.Put(n, iforce);
						break;
					}
				}

			}
			bool compatible = true;
			int nn = nnodesTotal;
			for (int n = 1; n < nn; n++)
			{
				if (compatibleForceNum.Get(n) == -1)
				{
					compatible = false;
					break;
				}
			}
			if (compatible)
			{
				//can't have more than one force on same face:
				if (newFaceForce != NULL)
					ERROREXIT;
				numFaceForce++;
				if (countOnly)
					break;
				face->pushBackForceids(model.GetNumForces());
				if (face->getNumForceIds() > 1)
					anyFaceHasMultipleLoads = true;
				newFaceForce = model.addForce();
				newFaceForce->setType(faceForce);
				newFaceForce->setFaceFromNodal(true);
				numFaceFromNodalLoad++;
				newFaceForce->setCoordId(coord0);
				newFaceForce->setPressure(pressure0);
				newFaceForce->SetEntityId(f);
				newFaceForce->AllocateForceVals(nnodesTotal, ndof);
				SRforce nodeForceTmp;
				nodeForceTmp.AllocateForceVals(1, 3);
				for (int n = 0; n < nnodesTotal; n++)
				{
					int iforce = compatibleForceNum.Get(n);
					int nid;
					if (n == 0)
						nodeForce = nodeForce0;
					else
					{
						nid = face->GetNodeOrMidNodeId(n);
						nodeForce = nodeForceStore.Get(nid).GetForce(iforce);
					}
					for (int dof = 0; dof < ndof; dof++)
					{
						double ft = nodeForce->GetForceVal(0, dof);
						newFaceForce->setForceVal(n, dof, ft);
					}
					//set this node force inactive to flag it should be removed below:
					nodeForce->setType(inactiveForce);
				}
			} //endif: if (compatible)
		} //for (int iforce0 = 0; iforce0 < nforce0; iforce0++)
	} //for (int f = 0; f < nface; f++)

	if (countOnly)
		return numFaceForce;

	//for biaxial loading at corner of two shared faces, an adjacent face might pick up the loaded dof at the edge only
	//filter this out:
	for (int f = 0; f < model.GetNumFaces(); f++)
	{
		SRface* face = model.GetFace(f);
		if (face->getNumForceIds() == 0)
			continue;
		SRforce* faceForce = model.GetForce(face->getForceId(0));
		int nnodesTotal = face->GetNumNodesTotal();
		int ndof = 3;
		if (faceForce->isPressure())
			ndof = 1;

		for (int dof = 0; dof < ndof; dof++)
		{
			int nej = face->GetNumLocalEdges();
			bool edgeNode[8];
			for (int lej = 0; lej < nej; lej++)
			{
				for (int n = 0; n < nnodesTotal; n++)
					edgeNode[n] = false;
				bool entireEdgeLoaded = true;
				for (int ejnode = 0; ejnode < 3; ejnode++)
				{
					int ln = model.map.GetFaceEdgeLocalNode(nej, lej, ejnode);
					edgeNode[ln] = true;
					double ft = faceForce->GetForceVal(ln, dof);
					if (fabs(ft) < TINY)
					{
						entireEdgeLoaded = false;
						break;
					}
				}
				if (!entireEdgeLoaded)
					continue;
				bool anyOtherNodeLoaded = false;
				for (int n = 0; n < nnodesTotal; n++)
				{
					if (edgeNode[n])
						continue;
					double ft = faceForce->GetForceVal(n, dof);
					if (fabs(ft) > TINY)
					{
						anyOtherNodeLoaded = true;
						break;
					}
				}
				if (!anyOtherNodeLoaded)
				{
					//only one edge is loaded for this dof. filter out by setting all forcevals for this dof to 0:
					for (int n = 0; n < nnodesTotal; n++)
						faceForce->setForceVal(n, dof, 0.0);
				}
			}

			//same check for only one node loaded, would happen at corner where 3 faces meet
			int numLoadedNodes = 0;
			for (int n = 0; n < nnodesTotal; n++)
			{
				double ft = faceForce->GetForceVal(n, dof);
				if (fabs(ft) > TINY)
					numLoadedNodes++;
			}
			if (numLoadedNodes == 1)
			{
				for (int n = 0; n < nnodesTotal; n++)
					faceForce->setForceVal(n, dof, 0.0);
			}
		} // end: for (int dof = 0; dof < 3; dof++)
	} // end: for (int f = 0; f < model.GetNumFaces(); f++)

	for (int f = 0; f < model.GetNumFaces(); f++)
	{
		SRface* face = model.GetFace(f);
		if (face->getNumForceIds() == 0)
			continue;
		//empty forceIds, they will be reassigned after inactive nodal forces removed
		face->freeForceIds();
	}

	return numFaceForce;
}

int SRinput::elemFaceFind(SRelement* elem, int nv[], int gno[])
{
	//find the face on an element that matches nodes in nv array
	//input:
		//elem = pointer to element
		//nv = corner nodes of faces. nv[3] = -1 for tri face
	//output:
		//gno = global node orders of the nodes on the face corresponding to nv
	//return
		//id of the global faces in the model with corner nodes nv; -1 if no match
	for (int l = 0; l < elem->GetNumLocalFaces(); l++)
	{
		SRface* face = elem->GetFace(l);
		if (SRfaceUtil::FaceMatch(nv, face->getNodeIds(), gno))
			return elem->GetLocalFaceGlobalId(l);
	}
	return -1;
}

void SRinput::FillAndSortNodeUids()
{
	//fill node uid vector and sort in ascending order of uid for faster
	//node-finding
	int nnode = model.GetNumNodes();
	if (nodeUidOffset == -1)
	{
		nodeUids.Allocate(nnode);
		SRuidData *nuid;
		for (int i = 0; i < nnode; i++)
		{
			SRnode* node = model.GetNode(i);
			nuid = nodeUids.GetPointer(i);
			nuid->id = i;
			nuid->uid = node->GetUserid();
		}
		SortNodes();
	}
}

SRconstraint* SRnodeConInputStore::GetCon(int i)
{
	int id = cons.Get(i);
	return model.GetConstraint(id);
}

SRforce* SRnodeForceInputStore::GetForce(int i)
{
	int id = forces.Get(i);
	return model.GetForce(id);
}

void SRinput::equivMatTest()
{
	//check all materials for equivalence with other materials
	int nm = model.GetNumMaterials();
	for (int i = 0; i < nm; i++)
	{
		SRmaterial* mati = model.GetMaterial(i);
		for (int j = i + 1; j < nm; j++)
		{
			SRmaterial* matj = model.GetMaterial(j);
			if (mati->diffElast(matj))
			{
				analysis.setAllMatsEquivElast(false);
				return;
			}
		}

	}
}

void SRinput::fillMultiFaceForceGroups()
{
	//fill multifaceforcegroups, which are groups of face forces that are detected to be on the same surface
	//note:
		//criterion for "on the same surface" is the faces sharing an adjacent edge do not have a kink in normals
		//exceeding tolerance of 10 degrees
	if (numFaceFromNodalLoad == 0)
		return;

	if (anyFaceHasMultipleLoads)
		return;

	nodeNeedsNodalForce.Allocate(model.GetNumNodes());
	nodeNeedsNodalForce.Set(-1);

	int nface = model.GetNumFaces();
	model.allocateFaceForceGroups(nface);
	for (int f = 0; f < nface; f++)
	{
		SRface* face = model.GetFace(f);
		if (!face->hasForce())
			continue;
		if (face->getMultifaceForceGroupId() != -1)
			continue;
		//create new surface for this face:
		int ffid = model.getNumFaceForceGroups();
		SRFaceForceGroup* ff = model.addFaceForceGroup();
		face->setMultifaceForceGroupId(ffid);
		ff->addFace(face);
		checkBoundaryEdgeNeighbors(ffid, face);
		int fid = face->getForceId(0);
		SRforce* force = model.GetForce(fid);
		if (force->isFaceFromNodal())
			ff->fromNodalForces = true;
		force->setFaceFromNodal(false);
	}

	SRintVector nodeonff;
	int nnode = model.GetNumNodes();
	nodeonff.Allocate(nnode);
	for (int ffn = 0; ffn < model.getNumFaceForceGroups(); ffn++)
	{
		nodeonff.Zero();
		SRFaceForceGroup* ff = model.getFaceForceGroup(ffn);
		for (int f = 0; f < ff->faceIds.GetNum(); f++)
		{
			SRface* face = model.GetFace(ff->faceIds.Get(f));
			for (int n = 0; n < face->GetNumNodesTotal(); n++)
			{
				int nid = face->GetNodeOrMidNodeId(n);
				nodeonff.Put(nid, 1);
			}
		}
		for (int n = 0; n < nnode; n++)
		{
			if (nodeonff.Get(n) == 1)
			{
				ff->nodeIds.PushBack(n);
				nodeNeedsNodalForce.Put(n, 1);
			}
		}
		ff->nodeIds.Sort();
	}
}

void SRinput::checkBoundaryEdgeNeighbors(int parentMultifaceForceGroupId, SRface* parentFace)
{
	//loop over this face's local edges. see if face neighbors are on same surface as parentMultifaceForceGroupId (using normal kink test)

	for (int lej = 0; lej < parentFace->GetNumLocalEdges(); lej++)
	{
		SRedge* edge = parentFace->GetEdge(lej);
		for (int ef = 0; ef < edge->getNumBoundaryfaceData(); ef++)
		{
			int fid2 = edge->getBoundaryFaceData(ef).faceId;
			if (fid2 == parentFace->GetId())
				continue;
			int lej2 = edge->getBoundaryFaceData(ef).localEdgeId;
			SRface* face2 = model.GetFace(fid2);
			if (!face2->hasForce())
				continue;
			if (!edge->checkKinkOK(parentFace, lej, face2, lej2))
				continue;
			if (face2->getMultifaceForceGroupId() != -1)
				continue;
			face2->setMultifaceForceGroupId(parentMultifaceForceGroupId);
			SRFaceForceGroup* ff = model.getFaceForceGroup(parentMultifaceForceGroupId);
			ff->addFace(face2);
			checkBoundaryEdgeNeighbors(parentMultifaceForceGroupId, face2);
		}
	}
}

void SRinput::EnergySmoothNodalForces()
{
	//convert nodel force values on faces belonging to a group to smooth tractions, using consistent energy approach

	if (!analysis.doEnergySmooth)
		return;

	int nfg = model.getNumFaceForceGroups();

	for (int g = 0; g < nfg; g++)
	{
		SRFaceForceGroup* ffg = model.getFaceForceGroup(g);
		if (!ffg->fromNodalForces)
			continue;
		SRpardiso parSolver;
		ffg->faceFunLoc.Allocate(8,ffg->faceIds.GetNum());
		double N[8];
		double rf, sf, w;
		int nfun = ffg->nodeIds.GetNum();
		ffg->ForceDof.Allocate(3, ffg->nodeIds.GetNum());
		for (int f = 0; f < ffg->faceIds.GetNum(); f++)
		{
			int fid = ffg->faceIds.Get(f);
			SRface* face = model.GetFace(fid);
			int forceid = face->getForceId(0);
			SRforce* faf = model.GetForce(forceid);
			int nfunface = face->GetNumNodesTotal();
			int fafmatlen = nfunface*(nfunface + 1) / 2;
			faf->allocateStiffMat(fafmatlen);
			faf->allocateStiffDiag(nfunface);
			int nz = 0;
			for (int r = 0; r < nfunface; r++)
			{
				faf->putStiffDiag(r, nz);
				nz += (nfunface - r);
			}
			for (int i = 0; i < nfunface; i++)
			{
				int nid = face->GetNodeOrMidNodeId(i);
				int loc = ffg->nodeIds.Find(nid);
				if (loc == -1)
					ERROREXIT;
				ffg->faceFunLoc.Put(i, f, loc);
			}
			int nint = model.math.FillGaussPoints(face);
			for (int gp = 0; gp < nint; gp++)
			{
				model.math.GetGP2d(gp, rf, sf, w);
				w *= (face->Jacobian(rf, sf));
				model.map.FaceShapeFunctions(face, rf, sf, N);
				for (int i = 0; i < nfunface; i++)
				{
					double niw = N[i] * w;
					int iloc = ffg->faceFunLoc.Get(i, f);
					for (int dof = 0; dof < 3; dof++)
						ffg->ForceDof.Put(dof, iloc, faf->GetForceVal(i, dof));
					for (int j = i; j < nfunface; j++)
					{
						double NN = N[j] * niw;
						int ij = parSolver.getFFGElemStiffLocation(faf, i, j);
						faf->catStiffMat(ij, NN);
					}
				}
			}
		}
		parSolver.solveFFG(ffg);
		SRdoubleVector Fv;
		Fv.Allocate(nfun);
		double *F = Fv.d;

		for (int dof = 0; dof < 3; dof++)
		{
			for (int j = 0; j < nfun; j++)
				F[j] = ffg->ForceDof.Get(dof, j);
			//F was overwritten with tractions:
			for (int f = 0; f < ffg->faceIds.GetNum(); f++)
			{
				int fid = ffg->faceIds.Get(f);
				SRface* face = model.GetFace(fid);
				int forceid = face->getForceId(0);
				SRforce* faf = model.GetForce(forceid);
				faf->setCoordId(-1);
				int nfunface = face->GetNumNodesTotal();
				for (int i = 0; i < nfunface; i++)
				{
					int iloc = ffg->faceFunLoc.Get(i, f);
					faf->setForceVal(i, dof, F[iloc]);
				}
				faf->freeStiffMat();
				faf->freeStiffDiag();
			}
		}
		ffg->faceFunLoc.Free();
		ffg->smoothStiff.Free();
	}

	bool echoSmoothTractions = false;
	if (!echoSmoothTractions)
		return;
	LOGPRINT("\nsmoothed face tractions");
	LOGPRINT(" node uid   tx      ty     tz");
	for (int g = 0; g < nfg; g++)
	{
		SRFaceForceGroup* ffg = model.getFaceForceGroup(g);
		if (!ffg->fromNodalForces)
			continue;
		for (int f = 0; f < ffg->faceIds.GetNum(); f++)
		{
			int fid = ffg->faceIds.d[f];
			SRface* face = model.GetFace(fid);
			int forceid = face->getForceId(0);
			SRforce* faf = model.GetForce(forceid);
			LOGPRINT("local face: %d", f);
			for (int n = 0; n < faf->GetNumForces(); n++)
			{
				int nuid = face->GetNodeOrMidnode(n)->GetUserid();
				LOGPRINTNORET("%d", nuid);
				for (int dof = 0; dof < 3; dof++)
				{
					double ft = faf->GetForceVal(n, dof);
					LOGPRINTNORET(" %lg", ft);
				}
				LOGPRINTRET;
			}
		}
	}
}




