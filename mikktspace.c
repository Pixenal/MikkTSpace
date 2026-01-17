/** \file mikktspace/mikktspace.c
 *  \ingroup mikktspace
 */
/**
 *  Copyright (C) 2011 by Morten S. Mikkelsen
 *
 *  This software is provided 'as-is', without any express or implied
 *  warranty.  In no event will the authors be held liable for any damages
 *  arising from the use of this software.
 *
 *  Permission is granted to anyone to use this software for any purpose,
 *  including commercial applications, and to alter it and redistribute it
 *  freely, subject to the following restrictions:
 *
 *  1. The origin of this software must not be misrepresented; you must not
 *     claim that you wrote the original software. If you use this software
 *     in a product, an acknowledgment in the product documentation would be
 *     appreciated but is not required.
 *  2. Altered source versions must be plainly marked as such, and must not be
 *     misrepresented as being the original software.
 *  3. This notice may not be removed or altered from any source distribution.
 */

#include <assert.h>
#include <stdlib.h>

#include "mikktspace.h"

static int32_t Build4RuleGroups(STriInfo pTriInfos[], SGroup pGroups[], int32_t piGroupTrianglesBuffer[], const int32_t piTriListIn[], const int32_t iNrTrianglesIn);

int32_t MakeIndex(const int32_t iFace, const int32_t iVert)
{
	assert(iVert>=0 && iVert<4 && iFace>=0);
	return (iFace<<2) | (iVert&0x3);
}

void IndexToData(int32_t * piFace, int32_t * piVert, const int32_t iIndexIn)
{
	piVert[0] = iIndexIn&0x3;
	piFace[0] = iIndexIn>>2;
}

STSpace AvgTSpace(const STSpace * pTS0, const STSpace * pTS1)
{
	STSpace ts_res;

	// this if is important. Due to floating point precision
	// averaging when ts0==ts1 will cause a slight difference
	// which results in tangent space splits later on
	if (pTS0->fMagS==pTS1->fMagS && pTS0->fMagT==pTS1->fMagT &&
	   veq(pTS0->vOs,pTS1->vOs)	&& veq(pTS0->vOt, pTS1->vOt))
	{
		ts_res.fMagS = pTS0->fMagS;
		ts_res.fMagT = pTS0->fMagT;
		ts_res.vOs = pTS0->vOs;
		ts_res.vOt = pTS0->vOt;
	}
	else
	{
		ts_res.fMagS = 0.5f*(pTS0->fMagS+pTS1->fMagS);
		ts_res.fMagT = 0.5f*(pTS0->fMagT+pTS1->fMagT);
		ts_res.vOs = vadd(pTS0->vOs,pTS1->vOs);
		ts_res.vOt = vadd(pTS0->vOt,pTS1->vOt);
		if ( VNotZero(ts_res.vOs) ) ts_res.vOs = Normalize(ts_res.vOs);
		if ( VNotZero(ts_res.vOt) ) ts_res.vOt = Normalize(ts_res.vOt);
	}

	return ts_res;
}

// based on the 4 rules, identify groups based on connectivity
bool makeGroups(SMikkTState *pState) {
	pState->iNrMaxGroups = pState->iNrTrianglesIn*3;
	pState->pGroups = (SGroup *) pState->pCtx->alloc.fpMalloc(sizeof(SGroup)*pState->iNrMaxGroups);
	pState->piGroupTrianglesBuffer = (int32_t *) pState->pCtx->alloc.fpMalloc(sizeof(int32_t)*pState->iNrTrianglesIn*3);
	if (pState->pGroups==NULL || pState->piGroupTrianglesBuffer==NULL) {
		if (pState->pGroups!=NULL) {
			pState->pCtx->alloc.fpFree(pState->pGroups);
		}
		if (pState->piGroupTrianglesBuffer!=NULL) {
			pState->pCtx->alloc.fpFree(pState->piGroupTrianglesBuffer);
		}
		pState->pCtx->alloc.fpFree(pState->piTriListIn);
		pState->pCtx->alloc.fpFree(pState->pTriInfos);
		return false;
	}
	pState->iNrActiveGroups = Build4RuleGroups(
		pState->pTriInfos,
		pState->pGroups,
		pState->piGroupTrianglesBuffer,
		pState->piTriListIn,
		pState->iNrTrianglesIn
	);
	return true;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// it is IMPORTANT that this function is called to evaluate the hash since
// inlining could potentially reorder instructions and generate different
// results for the same effective input value fVal.
int32_t FindGridCell(const float fMin, const float fMax, const float fVal)
{
	const float fIndex = MIKKT_CELLS * ((fVal-fMin)/(fMax-fMin));
	const int32_t iIndex = (int32_t)fIndex;
	return iIndex < MIKKT_CELLS ? (iIndex >= 0 ? iIndex : 0) : (MIKKT_CELLS - 1);
}

static bool AssignRecur(const int32_t piTriListIn[], STriInfo psTriInfos[], const int32_t iMyTriIndex, SGroup * pGroup);
static void AddTriToGroup(SGroup * pGroup, const int32_t iTriIndex);

static int32_t Build4RuleGroups(STriInfo pTriInfos[], SGroup pGroups[], int32_t piGroupTrianglesBuffer[], const int32_t piTriListIn[], const int32_t iNrTrianglesIn)
{
	const int32_t iNrMaxGroups = iNrTrianglesIn*3;
	int32_t iNrActiveGroups = 0;
	int32_t iOffset = 0, f=0, i=0;
	(void)iNrMaxGroups;  /* quiet warnings in non debug mode */
	for (f=0; f<iNrTrianglesIn; f++)
	{
		for (i=0; i<3; i++)
		{
			// if not assigned to a group
			if ((pTriInfos[f].iFlag&GROUP_WITH_ANY)==0 && pTriInfos[f].AssignedGroup[i]==NULL)
			{
				bool bOrPre;
				int32_t neigh_indexL, neigh_indexR;
				const int32_t vert_index = piTriListIn[f*3+i];
				assert(iNrActiveGroups<iNrMaxGroups);
				pTriInfos[f].AssignedGroup[i] = &pGroups[iNrActiveGroups];
				pTriInfos[f].AssignedGroup[i]->iVertexRepresentitive = vert_index;
				pTriInfos[f].AssignedGroup[i]->bOrientPreservering = (pTriInfos[f].iFlag&ORIENT_PRESERVING)!=0;
				pTriInfos[f].AssignedGroup[i]->iNrFaces = 0;
				pTriInfos[f].AssignedGroup[i]->pFaceIndices = &piGroupTrianglesBuffer[iOffset];
				++iNrActiveGroups;

				AddTriToGroup(pTriInfos[f].AssignedGroup[i], f);
				bOrPre = (pTriInfos[f].iFlag&ORIENT_PRESERVING)!=0 ? true : false;
				neigh_indexL = pTriInfos[f].FaceNeighbors[i];
				neigh_indexR = pTriInfos[f].FaceNeighbors[i>0?(i-1):2];
				if (neigh_indexL>=0) // neighbor
				{
					const bool bAnswer =
						AssignRecur(piTriListIn, pTriInfos, neigh_indexL,
									pTriInfos[f].AssignedGroup[i] );
					
					const bool bOrPre2 = (pTriInfos[neigh_indexL].iFlag&ORIENT_PRESERVING)!=0 ? true : false;
					const bool bDiff = bOrPre!=bOrPre2 ? true : false;
					assert(bAnswer || bDiff);
					(void)bAnswer, (void)bDiff;  /* quiet warnings in non debug mode */
				}
				if (neigh_indexR>=0) // neighbor
				{
					const bool bAnswer =
						AssignRecur(piTriListIn, pTriInfos, neigh_indexR,
									pTriInfos[f].AssignedGroup[i] );

					const bool bOrPre2 = (pTriInfos[neigh_indexR].iFlag&ORIENT_PRESERVING)!=0 ? true : false;
					const bool bDiff = bOrPre!=bOrPre2 ? true : false;
					assert(bAnswer || bDiff);
					(void)bAnswer, (void)bDiff;  /* quiet warnings in non debug mode */
				}

				// update offset
				iOffset += pTriInfos[f].AssignedGroup[i]->iNrFaces;
				// since the groups are disjoint a triangle can never
				// belong to more than 3 groups. Subsequently something
				// is completely screwed if this assertion ever hits.
				assert(iOffset <= iNrMaxGroups);
			}
		}
	}

	return iNrActiveGroups;
}

static void AddTriToGroup(SGroup * pGroup, const int32_t iTriIndex)
{
	pGroup->pFaceIndices[pGroup->iNrFaces] = iTriIndex;
	++pGroup->iNrFaces;
}

static bool AssignRecur(const int32_t piTriListIn[], STriInfo psTriInfos[],
				 const int32_t iMyTriIndex, SGroup * pGroup)
{
	STriInfo * pMyTriInfo = &psTriInfos[iMyTriIndex];

	// track down vertex
	const int32_t iVertRep = pGroup->iVertexRepresentitive;
	const int32_t * pVerts = &piTriListIn[3*iMyTriIndex+0];
	int32_t i=-1;
	if (pVerts[0]==iVertRep) i=0;
	else if (pVerts[1]==iVertRep) i=1;
	else if (pVerts[2]==iVertRep) i=2;
	assert(i>=0 && i<3);

	// early out
	if (pMyTriInfo->AssignedGroup[i] == pGroup) return true;
	else if (pMyTriInfo->AssignedGroup[i]!=NULL) return false;
	if ((pMyTriInfo->iFlag&GROUP_WITH_ANY)!=0)
	{
		// first to group with a group-with-anything triangle
		// determines it's orientation.
		// This is the only existing order dependency in the code!!
		if ( pMyTriInfo->AssignedGroup[0] == NULL &&
			pMyTriInfo->AssignedGroup[1] == NULL &&
			pMyTriInfo->AssignedGroup[2] == NULL )
		{
			pMyTriInfo->iFlag &= (~ORIENT_PRESERVING);
			pMyTriInfo->iFlag |= (pGroup->bOrientPreservering ? ORIENT_PRESERVING : 0);
		}
	}
	{
		const bool bOrient = (pMyTriInfo->iFlag&ORIENT_PRESERVING)!=0 ? true : false;
		if (bOrient != pGroup->bOrientPreservering) return false;
	}

	AddTriToGroup(pGroup, iMyTriIndex);
	pMyTriInfo->AssignedGroup[i] = pGroup;

	{
		const int32_t neigh_indexL = pMyTriInfo->FaceNeighbors[i];
		const int32_t neigh_indexR = pMyTriInfo->FaceNeighbors[i>0?(i-1):2];
		if (neigh_indexL>=0)
			AssignRecur(piTriListIn, psTriInfos, neigh_indexL, pGroup);
		if (neigh_indexR>=0)
			AssignRecur(piTriListIn, psTriInfos, neigh_indexR, pGroup);
	}
	return true;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////

bool CompareSubGroups(const SSubGroup * pg1, const SSubGroup * pg2)
{
	bool bStillSame=true;
	int32_t i=0;
	if (pg1->iNrFaces!=pg2->iNrFaces) return false;
	while (i<pg1->iNrFaces && bStillSame)
	{
		bStillSame = pg1->pTriMembers[i]==pg2->pTriMembers[i] ? true : false;
		if (bStillSame) ++i;
	}
	return bStillSame;
}

void QuickSort(int32_t* pSortBuffer, int32_t iLeft, int32_t iRight, uint32_t uSeed)
{
	int32_t iL, iR, n, index, iMid, iTmp;

	// Random
	uint32_t t=uSeed&31;
	t=(uSeed<<t)|(uSeed>>(32-t));
	uSeed=uSeed+t+3;
	// Random end

	iL=iLeft; iR=iRight;
	n = (iR-iL)+1;
	assert(n>=0);
	index = (int32_t) (uSeed%n);

	iMid=pSortBuffer[index + iL];

	do
	{
		while (pSortBuffer[iL] < iMid)
			++iL;
		while (pSortBuffer[iR] > iMid)
			--iR;

		if (iL <= iR)
		{
			iTmp = pSortBuffer[iL];
			pSortBuffer[iL] = pSortBuffer[iR];
			pSortBuffer[iR] = iTmp;
			++iL; --iR;
		}
	}
	while (iL <= iR);

	if (iLeft < iR)
		QuickSort(pSortBuffer, iLeft, iR, uSeed);
	if (iL < iRight)
		QuickSort(pSortBuffer, iL, iRight, uSeed);
}

/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

static void QuickSortEdges(SEdge * pSortBuffer, int32_t iLeft, int32_t iRight, const int32_t channel, uint32_t uSeed);
static void GetEdge(int32_t * i0_out, int32_t * i1_out, int32_t * edgenum_out, const int32_t indices[], const int32_t i0_in, const int32_t i1_in);

void BuildNeighborsFast(STriInfo pTriInfos[], SEdge * pEdges, const int32_t piTriListIn[], const int32_t iNrTrianglesIn)
{
	// build array of edges
	uint32_t uSeed = INTERNAL_RND_SORT_SEED;				// could replace with a random seed?
	int32_t iEntries=0, iCurStartIndex=-1, f=0, i=0;
	for (f=0; f<iNrTrianglesIn; f++)
		for (i=0; i<3; i++)
		{
			const int32_t i0 = piTriListIn[f*3+i];
			const int32_t i1 = piTriListIn[f*3+(i<2?(i+1):0)];
			pEdges[f*3+i].i0 = i0 < i1 ? i0 : i1;			// put minimum index in i0
			pEdges[f*3+i].i1 = !(i0 < i1) ? i0 : i1;		// put maximum index in i1
			pEdges[f*3+i].f = f;							// record face number
		}

	// sort over all edges by i0, this is the pricy one.
	QuickSortEdges(pEdges, 0, iNrTrianglesIn*3-1, 0, uSeed);	// sort channel 0 which is i0

	// sub sort over i1, should be fast.
	// could replace this with a 64 bit int sort over (i0,i1)
	// with i0 as msb in the quicksort call above.
	iEntries = iNrTrianglesIn*3;
	iCurStartIndex = 0;
	for (i=1; i<iEntries; i++)
	{
		if (pEdges[iCurStartIndex].i0 != pEdges[i].i0)
		{
			const int32_t iL = iCurStartIndex;
			const int32_t iR = i-1;
			//const int32_t iElems = i-iL;
			iCurStartIndex = i;
			QuickSortEdges(pEdges, iL, iR, 1, uSeed);	// sort channel 1 which is i1
		}
	}

	// sub sort over f, which should be fast.
	// this step is to remain compliant with BuildNeighborsSlow() when
	// more than 2 triangles use the same edge (such as a butterfly topology).
	iCurStartIndex = 0;
	for (i=1; i<iEntries; i++)
	{
		if (pEdges[iCurStartIndex].i0 != pEdges[i].i0 || pEdges[iCurStartIndex].i1 != pEdges[i].i1)
		{
			const int32_t iL = iCurStartIndex;
			const int32_t iR = i-1;
			//const int32_t iElems = i-iL;
			iCurStartIndex = i;
			QuickSortEdges(pEdges, iL, iR, 2, uSeed);	// sort channel 2 which is f
		}
	}

	// pair up, adjacent triangles
	for (i=0; i<iEntries; i++)
	{
		const int32_t i0=pEdges[i].i0;
		const int32_t i1=pEdges[i].i1;
		const int32_t f = pEdges[i].f;
		bool bUnassigned_A;

		int32_t i0_A, i1_A;
		int32_t edgenum_A, edgenum_B=0;	// 0,1 or 2
		GetEdge(&i0_A, &i1_A, &edgenum_A, &piTriListIn[f*3], i0, i1);	// resolve index ordering and edge_num
		bUnassigned_A = pTriInfos[f].FaceNeighbors[edgenum_A] == -1 ? true : false;

		if (bUnassigned_A)
		{
			// get true index ordering
			int32_t j=i+1, t;
			bool bNotFound = true;
			while (j<iEntries && i0==pEdges[j].i0 && i1==pEdges[j].i1 && bNotFound)
			{
				bool bUnassigned_B;
				int32_t i0_B, i1_B;
				t = pEdges[j].f;
				// flip i0_B and i1_B
				GetEdge(&i1_B, &i0_B, &edgenum_B, &piTriListIn[t*3], pEdges[j].i0, pEdges[j].i1);	// resolve index ordering and edge_num
				//assert(!(i0_A==i1_B && i1_A==i0_B));
				bUnassigned_B =  pTriInfos[t].FaceNeighbors[edgenum_B]==-1 ? true : false;
				if (i0_A==i0_B && i1_A==i1_B && bUnassigned_B)
					bNotFound = false;
				else
					++j;
			}

			if (!bNotFound)
			{
				int32_t t = pEdges[j].f;
				pTriInfos[f].FaceNeighbors[edgenum_A] = t;
				//assert(pTriInfos[t].FaceNeighbors[edgenum_B]==-1);
				pTriInfos[t].FaceNeighbors[edgenum_B] = f;
			}
		}
	}
}

void BuildNeighborsSlow(STriInfo pTriInfos[], const int32_t piTriListIn[], const int32_t iNrTrianglesIn)
{
	int32_t f=0, i=0;
	for (f=0; f<iNrTrianglesIn; f++)
	{
		for (i=0; i<3; i++)
		{
			// if unassigned
			if (pTriInfos[f].FaceNeighbors[i] == -1)
			{
				const int32_t i0_A = piTriListIn[f*3+i];
				const int32_t i1_A = piTriListIn[f*3+(i<2?(i+1):0)];

				// search for a neighbor
				bool bFound = false;
				int32_t t=0, j=0;
				while (!bFound && t<iNrTrianglesIn)
				{
					if (t!=f)
					{
						j=0;
						while (!bFound && j<3)
						{
							// in rev order
							const int32_t i1_B = piTriListIn[t*3+j];
							const int32_t i0_B = piTriListIn[t*3+(j<2?(j+1):0)];
							//assert(!(i0_A==i1_B && i1_A==i0_B));
							if (i0_A==i0_B && i1_A==i1_B)
								bFound = true;
							else
								++j;
						}
					}
					
					if (!bFound) ++t;
				}

				// assign neighbors
				if (bFound)
				{
					pTriInfos[f].FaceNeighbors[i] = t;
					//assert(pTriInfos[t].FaceNeighbors[j]==-1);
					pTriInfos[t].FaceNeighbors[j] = f;
				}
			}
		}
	}
}

static void QuickSortEdges(SEdge * pSortBuffer, int32_t iLeft, int32_t iRight, const int32_t channel, uint32_t uSeed)
{
	uint32_t t;
	int32_t iL, iR, n, index, iMid;

	// early out
	SEdge sTmp;
	const int32_t iElems = iRight-iLeft+1;
	if (iElems<2) return;
	else if (iElems==2)
	{
		if (pSortBuffer[iLeft].array[channel] > pSortBuffer[iRight].array[channel])
		{
			sTmp = pSortBuffer[iLeft];
			pSortBuffer[iLeft] = pSortBuffer[iRight];
			pSortBuffer[iRight] = sTmp;
		}
		return;
	}

	// Random
	t=uSeed&31;
	t=(uSeed<<t)|(uSeed>>(32-t));
	uSeed=uSeed+t+3;
	// Random end

	iL = iLeft;
	iR = iRight;
	n = (iR-iL)+1;
	assert(n>=0);
	index = (int32_t) (uSeed%n);

	iMid=pSortBuffer[index + iL].array[channel];

	do
	{
		while (pSortBuffer[iL].array[channel] < iMid)
			++iL;
		while (pSortBuffer[iR].array[channel] > iMid)
			--iR;

		if (iL <= iR)
		{
			sTmp = pSortBuffer[iL];
			pSortBuffer[iL] = pSortBuffer[iR];
			pSortBuffer[iR] = sTmp;
			++iL; --iR;
		}
	}
	while (iL <= iR);

	if (iLeft < iR)
		QuickSortEdges(pSortBuffer, iLeft, iR, channel, uSeed);
	if (iL < iRight)
		QuickSortEdges(pSortBuffer, iL, iRight, channel, uSeed);
}

// resolve ordering and edge number
static void GetEdge(int32_t * i0_out, int32_t * i1_out, int32_t * edgenum_out, const int32_t indices[], const int32_t i0_in, const int32_t i1_in)
{
	*edgenum_out = -1;
	
	// test if first index is on the edge
	if (indices[0]==i0_in || indices[0]==i1_in)
	{
		// test if second index is on the edge
		if (indices[1]==i0_in || indices[1]==i1_in)
		{
			edgenum_out[0]=0;	// first edge
			i0_out[0]=indices[0];
			i1_out[0]=indices[1];
		}
		else
		{
			edgenum_out[0]=2;	// third edge
			i0_out[0]=indices[2];
			i1_out[0]=indices[0];
		}
	}
	else
	{
		// only second and third index is on the edge
		edgenum_out[0]=1;	// second edge
		i0_out[0]=indices[1];
		i1_out[0]=indices[2];
	}
}


/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// Degenerate triangles ////////////////////////////////////

void DegenPrologue(STriInfo pTriInfos[], int32_t piTriList_out[], const int32_t iNrTrianglesIn, const int32_t iTotTris)
{
	int32_t iNextGoodTriangleSearchIndex=-1;
	bool bStillFindingGoodOnes;

	// locate quads with only one good triangle
	int32_t t=0;
	while (t<(iTotTris-1))
	{
		const int32_t iFO_a = pTriInfos[t].iOrgFaceNumber;
		const int32_t iFO_b = pTriInfos[t+1].iOrgFaceNumber;
		if (iFO_a==iFO_b)	// this is a quad
		{
			const bool bIsDeg_a = (pTriInfos[t].iFlag&MARK_DEGENERATE)!=0 ? true : false;
			const bool bIsDeg_b = (pTriInfos[t+1].iFlag&MARK_DEGENERATE)!=0 ? true : false;
			if ((bIsDeg_a^bIsDeg_b)!=0)
			{
				pTriInfos[t].iFlag |= QUAD_ONE_DEGEN_TRI;
				pTriInfos[t+1].iFlag |= QUAD_ONE_DEGEN_TRI;
			}
			t += 2;
		}
		else
			++t;
	}

	// reorder list so all degen triangles are moved to the back
	// without reordering the good triangles
	iNextGoodTriangleSearchIndex = 1;
	t=0;
	bStillFindingGoodOnes = true;
	while (t<iNrTrianglesIn && bStillFindingGoodOnes)
	{
		const bool bIsGood = (pTriInfos[t].iFlag&MARK_DEGENERATE)==0 ? true : false;
		if (bIsGood)
		{
			if (iNextGoodTriangleSearchIndex < (t+2))
				iNextGoodTriangleSearchIndex = t+2;
		}
		else
		{
			int32_t t0, t1;
			// search for the first good triangle.
			bool bJustADegenerate = true;
			while (bJustADegenerate && iNextGoodTriangleSearchIndex<iTotTris)
			{
				const bool bIsGood = (pTriInfos[iNextGoodTriangleSearchIndex].iFlag&MARK_DEGENERATE)==0 ? true : false;
				if (bIsGood) bJustADegenerate=false;
				else ++iNextGoodTriangleSearchIndex;
			}

			t0 = t;
			t1 = iNextGoodTriangleSearchIndex;
			++iNextGoodTriangleSearchIndex;
			assert(iNextGoodTriangleSearchIndex > (t+1));

			// swap triangle t0 and t1
			if (!bJustADegenerate)
			{
				int32_t i=0;
				for (i=0; i<3; i++)
				{
					const int32_t index = piTriList_out[t0*3+i];
					piTriList_out[t0*3+i] = piTriList_out[t1*3+i];
					piTriList_out[t1*3+i] = index;
				}
				{
					const STriInfo tri_info = pTriInfos[t0];
					pTriInfos[t0] = pTriInfos[t1];
					pTriInfos[t1] = tri_info;
				}
			}
			else
				bStillFindingGoodOnes = false;	// this is not supposed to happen
		}

		if (bStillFindingGoodOnes) ++t;
	}

	assert(bStillFindingGoodOnes);	// code will still work.
	assert(iNrTrianglesIn == t);
}
