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
// SRinput.h: interface for the SRinput class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(SRINPUT_INCLUDED)
#define SRINPUT_INCLUDED

#include "SRmodel.h"
#include "SRfile.h"
#include "SRstring.h"

#define countOnlyFlag true
#define addTheForcesFlag false
#define addTheConstraintsFlag false

int SRuidDataCompareFunc(const void *v1, const void *v2);

struct SRuidData
{
	int id;
	int uid;
};

class SRnodeConInputStore
{
public:
	SRintVector cons;
	SRconstraint* GetCon(int i);
	int GetConId(int i) {return cons.Get(i); };
};

class SRnodeForceInputStore
{
public:
	SRintVector forces;
	SRforce* GetForce(int i);
	int GetForceId(int i) { return forces.Get(i); };
};


class SRFaceForceGroup;

class SRinput  
{
public:
	void Cleanup();
	void SortNodes();
	void SortElems();
	void InputNodalConstraints();
	void InputFaceConstraints(int nfc);
	void InputThermal();
	int GetCoordId(SRstring& name);
	int GetMaterialId(SRstring& name);
	void InputCoordinates();
	void InputVolumeForces();
	void InputNodalForces();
	void InputFacePressures();
	void InputMaterials();
	void CountEntities(int& num);
	void CountElements(int& num);
	void InputElements();
	void InputNodes();
	int NodeFind(int uid);
	int elemFind(int uid);
	void ReadModel();
	int nodalToFaceConstraints(bool countOnly);
	int nodalToFaceForces(bool countOnly);
	int elemFaceFind(SRelement* elem, int nv[], int gno[]);
	void InputFaceTractions();
	void FillAndSortNodeUids();
	void equivMatTest();
	void InputMultiFaceConstraints(int nmultifaceConGroups);
	void InputMultiFaceForces(int nmultifaceForceGroups);
	void EnergySmoothNodalForces();
	void fillMultiFaceForceGroups();
	void checkBoundaryEdgeNeighbors(int parentMultifaceForceGroupId, SRface*parentFace);
	void readSettings();
	void doNodalToFaceForces();
	void doNodalToFaceConstraints();

	int CoordUidOffset;
	int MatUidOffset;
	int elPropUidOffset;
	int lastCoordUid;
	int lastCoordId;
	int lastMatUid;
	int lastMatId;
	int lastElPropUid;
	int lastElPropId;
	SRvector <SRuidData> coordUids;
	SRvector <SRuidData> matUids;
	SRvector <SRuidData> elpropUids;

	SRstring curLine;
	SRintVector nodeNeedsNodalForce;
	SRvector <SRnodeForceInputStore> currentNodeForceStore;
	SRpointerVector <SRforce> currentForces;
	SRvector <SRnodeConInputStore> currentNodeConStore;
	SRpointerVector <SRconstraint> currentConstraints;


	double val;
	int nodeUidOffset;
	int elemUidOffset;
	int numFaceFromNodalLoad;
	bool anyFaceHasMultipleLoads;
	bool userPSettings;
	int numSettingsStrings;
	int customPmax;
	int nnode;
	int nelem;
	int ncon;
	int nforce;
	int nvolforce;
	int nfacePressure;
	int nfaceTraction;
	int nFaceCon;
	int nmat;
	int ncoord;
	bool anyCoordsReferenceGrids;
	bool anyNodalForce;

	SRvector <SRuidData> nodeUids;
	SRvector <SRuidData> elemUids;
	SRvector <SRnodeConInputStore> nodeConStore;
	SRvector <SRnodeForceInputStore> nodeForceStore;
};

#endif // !defined(SRINPUT_INCLUDED)
