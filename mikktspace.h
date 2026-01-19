/** \file mikktspace/mikktspace.h
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

#pragma once

#ifdef WIN32
#define MIKKT_FORCE_INLINE __forceinline
#define MIKKT_NOINLINE __declspec(noinline)
#else
#define MIKKT_FORCE_INLINE __attribute__((always_inline)) static inline
#define MIKKT_NOINLINE __attribute__ ((noinline))
#endif

#ifdef __cplusplus
extern "C" {
#endif

#ifndef M_PI
#define M_PI	3.1415926535897932384626433832795
#endif

#define INTERNAL_RND_SORT_SEED		39871946

#define MARK_DEGENERATE				1
#define QUAD_ONE_DEGEN_TRI			2
#define GROUP_WITH_ANY				4
#define ORIENT_PRESERVING			8

#define MIKKT_CELLS 2048

/* Author: Morten S. Mikkelsen
 * Version: 1.0
 *
 * The files mikktspace.h and mikktspace.c are designed to be
 * stand-alone files and it is important that they are kept this way.
 * Not having dependencies on structures/classes/libraries specific
 * to the program, in which they are used, allows them to be copied
 * and used as is into any tool, program or plugin.
 * The code is designed to consistently generate the same
 * tangent spaces, for a given mesh, in any tool in which it is used.
 * This is done by performing an internal welding step and subsequently an order-independent evaluation
 * of tangent space for meshes consisting of triangles and quads.
 * This means faces can be received in any order and the same is true for
 * the order of vertices of each face. The generated result will not be affected
 * by such reordering. Additionally, whether degenerate (vertices or texture coordinates)
 * primitives are present or not will not affect the generated results either.
 * Once tangent space calculation is done the vertices of degenerate primitives will simply
 * inherit tangent space from neighboring non degenerate primitives.
 * The analysis behind this implementation can be found in my master's thesis
 * which is available for download --> http://image.diku.dk/projects/media/morten.mikkelsen.08.pdf
 * Note that though the tangent spaces at the vertices are generated in an order-independent way,
 * by this implementation, the interpolated tangent space is still affected by which diagonal is
 * chosen to split each quad. A sensible solution is to have your tools pipeline always
 * split quads by the shortest diagonal. This choice is order-independent and works with mirroring.
 * If these have the same length then compare the diagonals defined by the texture coordinates.
 * XNormal which is a tool for baking normal maps allows you to write your own tangent space plugin
 * and also quad triangulator plugin.
 */

#include <assert.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include <float.h>

#include <pixenals_thread_utils.h>

typedef struct SMikkTSpaceContext SMikkTSpaceContext;

typedef struct {
	// Returns the number of faces (triangles/quads) on the mesh to be processed.
	int (*m_getNumFaces)(const SMikkTSpaceContext * pContext);

	// Returns the number of vertices on face number iFace
	// iFace is a number in the range {0, 1, ..., getNumFaces()-1}
	int (*m_getNumVerticesOfFace)(const SMikkTSpaceContext * pContext, const int iFace);

	// returns the position/normal/texcoord of the referenced face of vertex number iVert.
	// iVert is in the range {0,1,2} for triangles and {0,1,2,3} for quads.
	void (*m_getPosition)(const SMikkTSpaceContext * pContext, float fvPosOut[], const int iFace, const int iVert);
	void (*m_getNormal)(const SMikkTSpaceContext * pContext, float fvNormOut[], const int iFace, const int iVert);
	void (*m_getTexCoord)(const SMikkTSpaceContext * pContext, float fvTexcOut[], const int iFace, const int iVert);

	// either (or both) of the two setTSpace callbacks can be set.
	// The call-back m_setTSpaceBasic() is sufficient for basic normal mapping.

	// This function is used to return the tangent and fSign to the application.
	// fvTangent is a unit length vector.
	// For normal maps it is sufficient to use the following simplified version of the bitangent which is generated at pixel/vertex level.
	// bitangent = fSign * cross(vN, tangent);
	// Note that the results are returned unindexed. It is possible to generate a new index list
	// But averaging/overwriting tangent spaces by using an already existing index list WILL produce INCRORRECT results.
	// DO NOT! use an already existing index list.
	void (*m_setTSpaceBasic)(const SMikkTSpaceContext * pContext, const float fvTangent[], const float fSign, const int iFace, const int iVert);

	// This function is used to return tangent space results to the application.
	// fvTangent and fvBiTangent are unit length vectors and fMagS and fMagT are their
	// true magnitudes which can be used for relief mapping effects.
	// fvBiTangent is the "real" bitangent and thus may not be perpendicular to fvTangent.
	// However, both are perpendicular to the vertex normal.
	// For normal maps it is sufficient to use the following simplified version of the bitangent which is generated at pixel/vertex level.
	// fSign = bIsOrientationPreserving ? 1.0f : (-1.0f);
	// bitangent = fSign * cross(vN, tangent);
	// Note that the results are returned unindexed. It is possible to generate a new index list
	// But averaging/overwriting tangent spaces by using an already existing index list WILL produce INCRORRECT results.
	// DO NOT! use an already existing index list.
	void (*m_setTSpace)(const SMikkTSpaceContext * pContext, const float fvTangent[], const float fvBiTangent[], const float fMagS, const float fMagT,
						const bool bIsOrientationPreserving, const int iFace, const int iVert);
} SMikkTSpaceInterface;

typedef struct MikKTFPtrs {
	void *(*fpMalloc)(size_t);
	void *(*fpCalloc)(size_t, size_t);
	void (*fpFree)(void *);
	void *(*fpRealloc)(void *, size_t);
} MikkTFPtrs;

struct SMikkTSpaceContext {
	SMikkTSpaceInterface * m_pInterface;	// initialized with callback functions
	void * m_pUserData;						// pointer to client side mesh data etc. (passed as the first parameter with every interface call)
	PixalcFPtrs alloc;
};

typedef struct {
	float x, y, z;
} SVec3;

typedef struct {
	int32_t iNrFaces;
	int32_t * pFaceIndices;
	int32_t iVertexRepresentitive;
	bool bOrientPreservering;
} SGroup;

typedef struct {
	int32_t FaceNeighbors[3];
	SGroup * AssignedGroup[3];
	
	// normalized first order face derivatives
	SVec3 vOs, vOt;
	float fMagS, fMagT;	// original magnitudes

	// determines if the current and the next triangle are a quad.
	int32_t iOrgFaceNumber;
	int32_t iFlag, iTSpacesOffs;
	unsigned char vert_num[4];
} STriInfo;

typedef struct {
	SVec3 vOs;
	float fMagS;
	SVec3 vOt;
	float fMagT;
	int32_t iCounter;	// this is to average back into quads.
	int32_t face;
	uint32_t corner : 31;
	uint32_t bOrient : 1;
} STSpace;

typedef struct SMikkTState {
	const SMikkTSpaceContext *pCtx;
	int32_t *piTriListIn;
	int32_t *piGroupTrianglesBuffer;
	STriInfo *pTriInfos;
	SGroup *pGroups;
	int32_t iNrTrianglesIn;
	int32_t f;
	int32_t t;
	int32_t i;
	int32_t iTotTris;
	int32_t iDegenTriangles;
	int32_t iNrMaxGroups;
	int32_t iNrActiveGroups;
	int32_t index;
	const int32_t iNrFaces;
	bool bRes;
	const float fThresCos;
} SMikkTState;

typedef struct {
	int32_t iNrFaces;
	int32_t * pTriMembers;
} SSubGroup;

typedef union {
	struct
	{
		int32_t i0, i1, f;
	};
	int32_t array[3];
} SEdge;

typedef struct {
	float vert[3];
	int32_t index;
} STmpVert;



static inline
bool veq( const SVec3 v1, const SVec3 v2 ) {
	return (v1.x == v2.x) && (v1.y == v2.y) && (v1.z == v2.z);
}

static inline
SVec3 vadd( const SVec3 v1, const SVec3 v2 ) {
	SVec3 vRes;
	vRes.x = v1.x + v2.x;
	vRes.y = v1.y + v2.y;
	vRes.z = v1.z + v2.z;
	return vRes;
}

static inline
SVec3 vsub( const SVec3 v1, const SVec3 v2 ) {
	SVec3 vRes;
	vRes.x = v1.x - v2.x;
	vRes.y = v1.y - v2.y;
	vRes.z = v1.z - v2.z;
	return vRes;
}

static inline
SVec3 vscale(const float fS, const SVec3 v) {
	SVec3 vRes;
	vRes.x = fS * v.x;
	vRes.y = fS * v.y;
	vRes.z = fS * v.z;
	return vRes;
}

static inline
float LengthSquared( const SVec3 v ) {
	return v.x*v.x + v.y*v.y + v.z*v.z;
}

static inline
float Length( const SVec3 v ) {
	return sqrtf(LengthSquared(v));
}

static inline
SVec3 Normalize( const SVec3 v ) {
	return vscale(1 / Length(v), v);
}

static inline
float vdot( const SVec3 v1, const SVec3 v2) {
	return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}

static inline
bool NotZero(const float fX) {
	// could possibly use FLT_EPSILON instead
	return fabsf(fX) > FLT_MIN;
}

static inline
bool VNotZero(const SVec3 v) {
	// might change this to an epsilon based test
	return NotZero(v.x) || NotZero(v.y) || NotZero(v.z);
}

int32_t MakeIndex(const int32_t iFace, const int32_t iVert);
void IndexToData(int32_t * piFace, int32_t * piVert, const int32_t iIndexIn);
STSpace AvgTSpace(const STSpace * pTS0, const STSpace * pTS1);
void DegenPrologue(STriInfo pTriInfos[], int32_t piTriList_out[], const int32_t iNrTrianglesIn, const int32_t iTotTris);
bool makeGroups(SMikkTState *pState);
MIKKT_NOINLINE int32_t FindGridCell(const float fMin, const float fMax, const float fVal);
void BuildNeighborsFast(STriInfo pTriInfos[], SEdge * pEdges, const int32_t piTriListIn[], const int32_t iNrTrianglesIn);
void BuildNeighborsSlow(STriInfo pTriInfos[], const int32_t piTriListIn[], const int32_t iNrTrianglesIn);
bool CompareSubGroups(const SSubGroup * pg1, const SSubGroup * pg2);
void QuickSort(int32_t* pSortBuffer, int32_t iLeft, int32_t iRight, uint32_t uSeed);


MIKKT_FORCE_INLINE
SVec3 GetPosition(const SMikkTSpaceContext * pContext, const int32_t index)
{
	int32_t iF, iI;
	SVec3 res; float pos[3];
	IndexToData(&iF, &iI, index);
	pContext->m_pInterface->m_getPosition(pContext, pos, iF, iI);
	res.x=pos[0]; res.y=pos[1]; res.z=pos[2];
	return res;
}

MIKKT_FORCE_INLINE
SVec3 GetNormal(const SMikkTSpaceContext * pContext, const int32_t index)
{
	int32_t iF, iI;
	SVec3 res; float norm[3];
	IndexToData(&iF, &iI, index);
	pContext->m_pInterface->m_getNormal(pContext, norm, iF, iI);
	res.x=norm[0]; res.y=norm[1]; res.z=norm[2];
	return res;
}

MIKKT_FORCE_INLINE
SVec3 GetTexCoord(const SMikkTSpaceContext * pContext, const int32_t index)
{
	int32_t iF, iI;
	SVec3 res; float texc[2];
	IndexToData(&iF, &iI, index);
	pContext->m_pInterface->m_getTexCoord(pContext, texc, iF, iI);
	res.x=texc[0]; res.y=texc[1]; res.z=1.0f;
	return res;
}

// returns the texture area times 2
MIKKT_FORCE_INLINE
float CalcTexArea(const SMikkTSpaceContext * pContext, const int32_t indices[])
{
	const SVec3 t1 = GetTexCoord(pContext, indices[0]);
	const SVec3 t2 = GetTexCoord(pContext, indices[1]);
	const SVec3 t3 = GetTexCoord(pContext, indices[2]);

	const float t21x = t2.x-t1.x;
	const float t21y = t2.y-t1.y;
	const float t31x = t3.x-t1.x;
	const float t31y = t3.y-t1.y;

	const float fSignedAreaSTx2 = t21x*t31y - t21y*t31x;

	return fSignedAreaSTx2<0 ? (-fSignedAreaSTx2) : fSignedAreaSTx2;
}

static inline
void MergeVertsFast(int32_t piTriList_in_and_out[], STmpVert pTmpVert[], const SMikkTSpaceContext * pContext, const int32_t iL_in, const int32_t iR_in)
{
	// make bbox
	int32_t c=0, l=0, channel=0;
	float fvMin[3], fvMax[3];
	float dx=0, dy=0, dz=0, fSep=0;
	for (c=0; c<3; c++)
	{	fvMin[c]=pTmpVert[iL_in].vert[c]; fvMax[c]=fvMin[c];	}
	for (l=(iL_in+1); l<=iR_in; l++) {
		for (c=0; c<3; c++) {
			if (fvMin[c]>pTmpVert[l].vert[c]) fvMin[c]=pTmpVert[l].vert[c];
			if (fvMax[c]<pTmpVert[l].vert[c]) fvMax[c]=pTmpVert[l].vert[c];
		}
	}

	dx = fvMax[0]-fvMin[0];
	dy = fvMax[1]-fvMin[1];
	dz = fvMax[2]-fvMin[2];

	channel = 0;
	if (dy>dx && dy>dz) channel=1;
	else if (dz>dx) channel=2;

	fSep = 0.5f*(fvMax[channel]+fvMin[channel]);

	// stop if all vertices are NaNs
	if (!isfinite(fSep))
		return;

	// terminate recursion when the separation/average value
	// is no longer strictly between fMin and fMax values.
	if (fSep>=fvMax[channel] || fSep<=fvMin[channel])
	{
		// complete the weld
		for (l=iL_in; l<=iR_in; l++)
		{
			int32_t i = pTmpVert[l].index;
			const int32_t index = piTriList_in_and_out[i];
			const SVec3 vP = GetPosition(pContext, index);
			const SVec3 vN = GetNormal(pContext, index);
			const SVec3 vT = GetTexCoord(pContext, index);

			bool bNotFound = true;
			int32_t l2=iL_in, i2rec=-1;
			while (l2<l && bNotFound)
			{
				const int32_t i2 = pTmpVert[l2].index;
				const int32_t index2 = piTriList_in_and_out[i2];
				const SVec3 vP2 = GetPosition(pContext, index2);
				const SVec3 vN2 = GetNormal(pContext, index2);
				const SVec3 vT2 = GetTexCoord(pContext, index2);
				i2rec=i2;

				//if (vP==vP2 && vN==vN2 && vT==vT2)
				if (vP.x==vP2.x && vP.y==vP2.y && vP.z==vP2.z &&
					vN.x==vN2.x && vN.y==vN2.y && vN.z==vN2.z &&
					vT.x==vT2.x && vT.y==vT2.y && vT.z==vT2.z)
					bNotFound = false;
				else
					++l2;
			}
			
			// merge if previously found
			if (!bNotFound)
				piTriList_in_and_out[i] = piTriList_in_and_out[i2rec];
		}
	}
	else
	{
		int32_t iL=iL_in, iR=iR_in;
		assert((iR_in-iL_in)>0);	// at least 2 entries

		// separate (by fSep) all points between iL_in and iR_in in pTmpVert[]
		while (iL < iR)
		{
			bool bReadyLeftSwap = false, bReadyRightSwap = false;
			while ((!bReadyLeftSwap) && iL<iR)
			{
				assert(iL>=iL_in && iL<=iR_in);
				bReadyLeftSwap = !(pTmpVert[iL].vert[channel]<fSep);
				if (!bReadyLeftSwap) ++iL;
			}
			while ((!bReadyRightSwap) && iL<iR)
			{
				assert(iR>=iL_in && iR<=iR_in);
				bReadyRightSwap = pTmpVert[iR].vert[channel]<fSep;
				if (!bReadyRightSwap) --iR;
			}
			assert( (iL<iR) || !(bReadyLeftSwap && bReadyRightSwap) );

			if (bReadyLeftSwap && bReadyRightSwap)
			{
				const STmpVert sTmp = pTmpVert[iL];
				assert(iL<iR);
				pTmpVert[iL] = pTmpVert[iR];
				pTmpVert[iR] = sTmp;
				++iL; --iR;
			}
		}

		assert(iL==(iR+1) || (iL==iR));
		if (iL==iR)
		{
			const bool bReadyRightSwap = pTmpVert[iR].vert[channel]<fSep;
			if (bReadyRightSwap) ++iL;
			else --iR;
		}

		// only need to weld when there is more than 1 instance of the (x,y,z)
		if (iL_in < iR)
			MergeVertsFast(piTriList_in_and_out, pTmpVert, pContext, iL_in, iR);	// weld all left of fSep
		if (iL < iR_in)
			MergeVertsFast(piTriList_in_and_out, pTmpVert, pContext, iL, iR_in);	// weld all right of (or equal to) fSep
	}
}

MIKKT_FORCE_INLINE
void MergeVertsSlow(int32_t piTriList_in_and_out[], const SMikkTSpaceContext * pContext, const int32_t pTable[], const int32_t iEntries)
{
	// this can be optimized further using a tree structure or more hashing.
	int32_t e=0;
	for (e=0; e<iEntries; e++)
	{
		int32_t i = pTable[e];
		const int32_t index = piTriList_in_and_out[i];
		const SVec3 vP = GetPosition(pContext, index);
		const SVec3 vN = GetNormal(pContext, index);
		const SVec3 vT = GetTexCoord(pContext, index);

		bool bNotFound = true;
		int32_t e2=0, i2rec=-1;
		while (e2<e && bNotFound)
		{
			const int32_t i2 = pTable[e2];
			const int32_t index2 = piTriList_in_and_out[i2];
			const SVec3 vP2 = GetPosition(pContext, index2);
			const SVec3 vN2 = GetNormal(pContext, index2);
			const SVec3 vT2 = GetTexCoord(pContext, index2);
			i2rec = i2;

			if (veq(vP,vP2) && veq(vN,vN2) && veq(vT,vT2))
				bNotFound = false;
			else
				++e2;
		}
		
		// merge if previously found
		if (!bNotFound)
			piTriList_in_and_out[i] = piTriList_in_and_out[i2rec];
	}
}

MIKKT_FORCE_INLINE
STSpace EvalTspace(int32_t face_indices[], const int32_t iFaces, const int32_t piTriListIn[], const STriInfo pTriInfos[],
                          const SMikkTSpaceContext * pContext, const int32_t iVertexRepresentitive)
{
	STSpace res;
	float fAngleSum = 0;
	int32_t face=0;
	res.vOs.x=0.0f; res.vOs.y=0.0f; res.vOs.z=0.0f;
	res.vOt.x=0.0f; res.vOt.y=0.0f; res.vOt.z=0.0f;
	res.fMagS = 0; res.fMagT = 0;

	for (face=0; face<iFaces; face++)
	{
		const int32_t f = face_indices[face];

		// only valid triangles get to add their contribution
		if ( (pTriInfos[f].iFlag&GROUP_WITH_ANY)==0 )
		{
			SVec3 n, vOs, vOt, p0, p1, p2, v1, v2;
			float fCos, fAngle, fMagS, fMagT;
			int32_t i=-1, index=-1, i0=-1, i1=-1, i2=-1;
			if (piTriListIn[3*f+0]==iVertexRepresentitive) i=0;
			else if (piTriListIn[3*f+1]==iVertexRepresentitive) i=1;
			else if (piTriListIn[3*f+2]==iVertexRepresentitive) i=2;
			assert(i>=0 && i<3);

			// project
			index = piTriListIn[3*f+i];
			n = GetNormal(pContext, index);
			vOs = vsub(pTriInfos[f].vOs, vscale(vdot(n,pTriInfos[f].vOs), n));
			vOt = vsub(pTriInfos[f].vOt, vscale(vdot(n,pTriInfos[f].vOt), n));
			if ( VNotZero(vOs) ) vOs = Normalize(vOs);
			if ( VNotZero(vOt) ) vOt = Normalize(vOt);

			i2 = piTriListIn[3*f + (i<2?(i+1):0)];
			i1 = piTriListIn[3*f + i];
			i0 = piTriListIn[3*f + (i>0?(i-1):2)];

			p0 = GetPosition(pContext, i0);
			p1 = GetPosition(pContext, i1);
			p2 = GetPosition(pContext, i2);
			v1 = vsub(p0,p1);
			v2 = vsub(p2,p1);

			// project
			v1 = vsub(v1, vscale(vdot(n,v1),n)); if ( VNotZero(v1) ) v1 = Normalize(v1);
			v2 = vsub(v2, vscale(vdot(n,v2),n)); if ( VNotZero(v2) ) v2 = Normalize(v2);

			// weight contribution by the angle
			// between the two edge vectors
			fCos = vdot(v1,v2); fCos=fCos>1?1:(fCos<(-1) ? (-1) : fCos);
			fAngle = (float) acos(fCos);
			fMagS = pTriInfos[f].fMagS;
			fMagT = pTriInfos[f].fMagT;

			res.vOs=vadd(res.vOs, vscale(fAngle,vOs));
			res.vOt=vadd(res.vOt,vscale(fAngle,vOt));
			res.fMagS+=(fAngle*fMagS);
			res.fMagT+=(fAngle*fMagT);
			fAngleSum += fAngle;
		}
	}

	// normalize
	if ( VNotZero(res.vOs) ) res.vOs = Normalize(res.vOs);
	if ( VNotZero(res.vOt) ) res.vOt = Normalize(res.vOt);
	if (fAngleSum>0)
	{
		res.fMagS /= fAngleSum;
		res.fMagT /= fAngleSum;
	}

	return res;
}

MIKKT_FORCE_INLINE
bool GenerateTSpaces(PixalcLinAlloc *pTSpaceAlloc, PixtyRange range, const STriInfo pTriInfos[], const SGroup pGroups[],
                             const int32_t piTriListIn[], const float fThresCos,
                             const SMikkTSpaceContext * pContext)
{
	typedef struct TSpaceLookup {
		uint32_t idx : 31;
		uint32_t valid : 1;
	} TSpaceLookup;
	typedef struct TSpaceTable {
		TSpaceLookup *pArr;
		int32_t size;
		int32_t count;
	} TSpaceTable;
	STSpace * pSubGroupTspace = NULL;
	SSubGroup * pUniSubGroups = NULL;
	int32_t * pTmpMembers = NULL;
	int32_t iMaxNrFaces=0, iUniqueTspaces=0, g=0, i=0;
	for (g=range.start; g<range.end; g++)
		if (iMaxNrFaces < pGroups[g].iNrFaces)
			iMaxNrFaces = pGroups[g].iNrFaces;

	if (iMaxNrFaces == 0) return true;

	// make initial allocations
	pSubGroupTspace =  pContext->alloc.fpMalloc(sizeof(STSpace)*iMaxNrFaces);
	pUniSubGroups =  pContext->alloc.fpMalloc(sizeof(SSubGroup)*iMaxNrFaces);
	pTmpMembers =  pContext->alloc.fpMalloc(sizeof(int32_t)*iMaxNrFaces);
	if (pSubGroupTspace==NULL || pUniSubGroups==NULL || pTmpMembers==NULL)
	{
		if (pSubGroupTspace!=NULL) pContext->alloc.fpFree(pSubGroupTspace);
		if (pUniSubGroups!=NULL) pContext->alloc.fpFree(pUniSubGroups);
		if (pTmpMembers!=NULL) pContext->alloc.fpFree(pTmpMembers);
		return false;
	}
	
	TSpaceTable table = {0};

	iUniqueTspaces = 0;
	for (g=range.start; g<range.end; g++)
	{
		const SGroup * pGroup = &pGroups[g];
		int32_t iUniqueSubGroups = 0, s=0;

		for (i=0; i<pGroup->iNrFaces; i++)	// triangles
		{
			const int32_t f = pGroup->pFaceIndices[i];	// triangle number
			int32_t index=-1, iVertIndex=-1, iOF_1=-1, iMembers=0, j=0, l=0;
			SSubGroup tmp_group;
			bool bFound;
			SVec3 n, vOs, vOt;
			if (pTriInfos[f].AssignedGroup[0]==pGroup) index=0;
			else if (pTriInfos[f].AssignedGroup[1]==pGroup) index=1;
			else if (pTriInfos[f].AssignedGroup[2]==pGroup) index=2;
			assert(index>=0 && index<3);

			iVertIndex = piTriListIn[f*3+index];
			assert(iVertIndex==pGroup->iVertexRepresentitive);

			// is normalized already
			n = GetNormal(pContext, iVertIndex);
			
			// project
			vOs = vsub(pTriInfos[f].vOs, vscale(vdot(n,pTriInfos[f].vOs), n));
			vOt = vsub(pTriInfos[f].vOt, vscale(vdot(n,pTriInfos[f].vOt), n));
			if ( VNotZero(vOs) ) vOs = Normalize(vOs);
			if ( VNotZero(vOt) ) vOt = Normalize(vOt);

			// original face number
			iOF_1 = pTriInfos[f].iOrgFaceNumber;
			
			iMembers = 0;
			for (j=0; j<pGroup->iNrFaces; j++)
			{
				const int32_t t = pGroup->pFaceIndices[j];	// triangle number
				const int32_t iOF_2 = pTriInfos[t].iOrgFaceNumber;

				// project
				SVec3 vOs2 = vsub(pTriInfos[t].vOs, vscale(vdot(n,pTriInfos[t].vOs), n));
				SVec3 vOt2 = vsub(pTriInfos[t].vOt, vscale(vdot(n,pTriInfos[t].vOt), n));
				if ( VNotZero(vOs2) ) vOs2 = Normalize(vOs2);
				if ( VNotZero(vOt2) ) vOt2 = Normalize(vOt2);

				{
					const bool bAny = ( (pTriInfos[f].iFlag | pTriInfos[t].iFlag) & GROUP_WITH_ANY )!=0 ? true : false;
					// make sure triangles which belong to the same quad are joined.
					const bool bSameOrgFace = iOF_1==iOF_2 ? true : false;

					const float fCosS = vdot(vOs,vOs2);
					const float fCosT = vdot(vOt,vOt2);

					assert(f!=t || bSameOrgFace);	// sanity check
					if (bAny || bSameOrgFace || (fCosS>fThresCos && fCosT>fThresCos))
						pTmpMembers[iMembers++] = t;
				}
			}

			// sort pTmpMembers
			tmp_group.iNrFaces = iMembers;
			tmp_group.pTriMembers = pTmpMembers;
			if (iMembers>1)
			{
				uint32_t uSeed = INTERNAL_RND_SORT_SEED;	// could replace with a random seed?
				QuickSort(pTmpMembers, 0, iMembers-1, uSeed);
			}

			// look for an existing match
			bFound = false;
			l=0;
			while (l<iUniqueSubGroups && !bFound)
			{
				bFound = CompareSubGroups(&tmp_group, &pUniSubGroups[l]);
				if (!bFound) ++l;
			}
			
			// assign tangent space index
			assert(bFound || l==iUniqueSubGroups);
			//piTempTangIndices[f*3+index] = iUniqueTspaces+l;

			// if no match was found we allocate a new subgroup
			if (!bFound)
			{
				// insert new subgroup
				int32_t * pIndices = (int32_t *) pContext->alloc.fpMalloc(sizeof(int32_t)*iMembers);
				if (pIndices==NULL)
				{
					// clean up and return false
					int32_t s=0;
					for (s=0; s<iUniqueSubGroups; s++)
						pContext->alloc.fpFree(pUniSubGroups[s].pTriMembers);
					pContext->alloc.fpFree(pUniSubGroups);
					pContext->alloc.fpFree(pTmpMembers);
					pContext->alloc.fpFree(pSubGroupTspace);
					return false;
				}
				pUniSubGroups[iUniqueSubGroups].iNrFaces = iMembers;
				pUniSubGroups[iUniqueSubGroups].pTriMembers = pIndices;
				memcpy(pIndices, tmp_group.pTriMembers, iMembers*sizeof(int32_t));
				pSubGroupTspace[iUniqueSubGroups] =
					EvalTspace(tmp_group.pTriMembers, iMembers, piTriListIn, pTriInfos, pContext, pGroup->iVertexRepresentitive);
				++iUniqueSubGroups;
			}

			assert(((pTriInfos[f].iFlag&ORIENT_PRESERVING)!=0) == pGroup->bOrientPreservering);
			// output tspace
			{
				const int32_t iOffs = pTriInfos[f].iTSpacesOffs;//idx of first corner in face
				const int32_t iVert = pTriInfos[f].vert_num[index];
				int32_t idx = iOffs + iVert;
				TSpaceLookup *pLookup = NULL;
				if (idx < table.count) {
					pLookup = table.pArr + idx;
				}
				else {
					//this macro resizes exponentially to fit the size specified
					int32_t oldSize = table.size;
					PIXALC_DYN_ARR_RESIZE(TSpaceLookup, &pContext->alloc, &table, idx + 1);
					if (table.size != oldSize) {
						PIX_ERR_ASSERT("", table.size > oldSize);
						memset(table.pArr + oldSize, 0, sizeof(TSpaceLookup) * (table.size - oldSize));
					}
					table.count = idx + 1;
					pLookup = table.pArr + idx;
				}
				STSpace *pTS_out = NULL;
				if (pLookup->valid) {
					pTS_out = pixalcLinAllocIdx(pTSpaceAlloc, pLookup->idx);
					assert(pTS_out->iCounter < 2);
					*pTS_out = AvgTSpace(pTS_out, &pSubGroupTspace[l]);
					pTS_out->iCounter = 2;
				}
				else {
					pLookup->idx = pixalcLinAlloc(pTSpaceAlloc, (void **)&pTS_out, 1);
					pLookup->valid = true;
					assert(pTS_out->iCounter==0);
					*pTS_out = pSubGroupTspace[l];
					pTS_out->iCounter = 1;
				}
				pTS_out->face = iOF_1;
				pTS_out->corner = iVert;
				pTS_out->bOrient = pGroup->bOrientPreservering;
			}
		}

		// clean up and offset iUniqueTspaces
		for (s=0; s<iUniqueSubGroups; s++)
			pContext->alloc.fpFree(pUniSubGroups[s].pTriMembers);
		iUniqueTspaces += iUniqueSubGroups;
	}

	// clean up
	pContext->alloc.fpFree(pUniSubGroups);
	pContext->alloc.fpFree(pTmpMembers);
	pContext->alloc.fpFree(pSubGroupTspace);
	if (table.pArr) {
		pContext->alloc.fpFree(table.pArr);
	}
	return true;
}

typedef struct GenTSpaceArgsShared {
	const SMikkTSpaceContext *pContext;
	const STriInfo *pTriInfos;
	const SGroup *pGroups;
	const int32_t *piTriListIn;
	const int32_t iNrActiveGroups;
	const float fThresCos;
} GenTSpaceArgsShared;

typedef struct GenTSpaceArgs {
	const GenTSpaceArgsShared *pShared;
	PixalcLinAlloc alloc;
	PixtyRange range;
} GenTSpaceArgs;

static inline
PixErr genTSpacesForJob(void *pArgsVoid) {
	PixErr err = PIX_ERR_SUCCESS;
	GenTSpaceArgs *pArgs = pArgsVoid;
	bool result = GenerateTSpaces(
		&pArgs->alloc,
		pArgs->range,
		pArgs->pShared->pTriInfos,
		pArgs->pShared->pGroups,
		pArgs->pShared->piTriListIn,
		pArgs->pShared->fThresCos,
		pArgs->pShared->pContext
	);
	PIX_ERR_RETURN_IFNOT_COND(err, result, "");
	return err;
}

// count triangles on supported faces
MIKKT_FORCE_INLINE
void countTris(SMikkTState *pState) {
	for (int32_t i = 0; i < pState->iNrFaces; i++) {
		const int32_t verts = pState->pCtx->m_pInterface->m_getNumVerticesOfFace(pState->pCtx, i);
		if (verts == 3) {
			++pState->iNrTrianglesIn;
		}
		else if (verts == 4) {
			pState->iNrTrianglesIn += 2;
		}
	}
}

MIKKT_FORCE_INLINE
int32_t GenerateInitialVerticesIndexList(SMikkTState *pState) {
	int32_t iTSpacesOffs = 0, f=0, t=0;
	int32_t iDstTriIndex = 0;
	for (f=0; f<pState->pCtx->m_pInterface->m_getNumFaces(pState->pCtx); f++) {
		const int32_t verts = pState->pCtx->m_pInterface->m_getNumVerticesOfFace(pState->pCtx, f);
		if (verts!=3 && verts!=4) {
			continue;
		}
		pState->pTriInfos[iDstTriIndex].iOrgFaceNumber = f;
		pState->pTriInfos[iDstTriIndex].iTSpacesOffs = iTSpacesOffs;

		if (verts==3) {
			unsigned char * pVerts = pState->pTriInfos[iDstTriIndex].vert_num;
			pVerts[0]=0; pVerts[1]=1; pVerts[2]=2;
			pState->piTriListIn[iDstTriIndex*3+0] = MakeIndex(f, 0);
			pState->piTriListIn[iDstTriIndex*3+1] = MakeIndex(f, 1);
			pState->piTriListIn[iDstTriIndex*3+2] = MakeIndex(f, 2);
			++iDstTriIndex;	// next
		}
		else {
			{
				pState->pTriInfos[iDstTriIndex+1].iOrgFaceNumber = f;
				pState->pTriInfos[iDstTriIndex+1].iTSpacesOffs = iTSpacesOffs;
			}
			{
				// need an order independent way to evaluate
				// tspace on quads. This is done by splitting
				// along the shortest diagonal.
				const int32_t i0 = MakeIndex(f, 0);
				const int32_t i1 = MakeIndex(f, 1);
				const int32_t i2 = MakeIndex(f, 2);
				const int32_t i3 = MakeIndex(f, 3);
				const SVec3 T0 = GetTexCoord(pState->pCtx, i0);
				const SVec3 T1 = GetTexCoord(pState->pCtx, i1);
				const SVec3 T2 = GetTexCoord(pState->pCtx, i2);
				const SVec3 T3 = GetTexCoord(pState->pCtx, i3);
				const float distSQ_02 = LengthSquared(vsub(T2,T0));
				const float distSQ_13 = LengthSquared(vsub(T3,T1));
				bool bQuadDiagIs_02;
				if (distSQ_02<distSQ_13) {
					bQuadDiagIs_02 = true;
				}
				else if (distSQ_13<distSQ_02) {
					bQuadDiagIs_02 = false;
				}
				else {
					const SVec3 P0 = GetPosition(pState->pCtx, i0);
					const SVec3 P1 = GetPosition(pState->pCtx, i1);
					const SVec3 P2 = GetPosition(pState->pCtx, i2);
					const SVec3 P3 = GetPosition(pState->pCtx, i3);
					const float distSQ_02 = LengthSquared(vsub(P2,P0));
					const float distSQ_13 = LengthSquared(vsub(P3,P1));

					bQuadDiagIs_02 = distSQ_13<distSQ_02 ? false : true;
				}

				if (bQuadDiagIs_02) {
					{
						unsigned char * pVerts_A = pState->pTriInfos[iDstTriIndex].vert_num;
						pVerts_A[0]=0; pVerts_A[1]=1; pVerts_A[2]=2;
					}
					pState->piTriListIn[iDstTriIndex*3+0] = i0;
					pState->piTriListIn[iDstTriIndex*3+1] = i1;
					pState->piTriListIn[iDstTriIndex*3+2] = i2;
					++iDstTriIndex;	// next
					{
						unsigned char * pVerts_B = pState->pTriInfos[iDstTriIndex].vert_num;
						pVerts_B[0]=0; pVerts_B[1]=2; pVerts_B[2]=3;
					}
					pState->piTriListIn[iDstTriIndex*3+0] = i0;
					pState->piTriListIn[iDstTriIndex*3+1] = i2;
					pState->piTriListIn[iDstTriIndex*3+2] = i3;
					++iDstTriIndex;	// next
				}
				else {
					{
						unsigned char * pVerts_A = pState->pTriInfos[iDstTriIndex].vert_num;
						pVerts_A[0]=0; pVerts_A[1]=1; pVerts_A[2]=3;
					}
					pState->piTriListIn[iDstTriIndex*3+0] = i0;
					pState->piTriListIn[iDstTriIndex*3+1] = i1;
					pState->piTriListIn[iDstTriIndex*3+2] = i3;
					++iDstTriIndex;	// next
					{
						unsigned char * pVerts_B = pState->pTriInfos[iDstTriIndex].vert_num;
						pVerts_B[0]=1; pVerts_B[1]=2; pVerts_B[2]=3;
					}
					pState->piTriListIn[iDstTriIndex*3+0] = i1;
					pState->piTriListIn[iDstTriIndex*3+1] = i2;
					pState->piTriListIn[iDstTriIndex*3+2] = i3;
					++iDstTriIndex;	// next
				}
			}
		}
		iTSpacesOffs += verts;
		assert(iDstTriIndex<=pState->iNrTrianglesIn);
	}
	for (t=0; t<pState->iNrTrianglesIn; t++) {
		pState->pTriInfos[t].iFlag = 0;
	}
	// return total amount of tspaces
	return iTSpacesOffs;
}

MIKKT_FORCE_INLINE
void GenerateSharedVerticesIndexListSlow(int32_t piTriList_in_and_out[], const SMikkTSpaceContext * pContext, const int32_t iNrTrianglesIn) {
	int32_t iNumUniqueVerts = 0, t=0, i=0;
	for (t=0; t<iNrTrianglesIn; t++)
	{
		for (i=0; i<3; i++)
		{
			const int32_t offs = t*3 + i;
			const int32_t index = piTriList_in_and_out[offs];

			const SVec3 vP = GetPosition(pContext, index);
			const SVec3 vN = GetNormal(pContext, index);
			const SVec3 vT = GetTexCoord(pContext, index);

			bool bFound = false;
			int32_t t2=0, index2rec=-1;
			while (!bFound && t2<=t)
			{
				int32_t j=0;
				while (!bFound && j<3)
				{
					const int32_t index2 = piTriList_in_and_out[t2*3 + j];
					const SVec3 vP2 = GetPosition(pContext, index2);
					const SVec3 vN2 = GetNormal(pContext, index2);
					const SVec3 vT2 = GetTexCoord(pContext, index2);
					
					if (veq(vP,vP2) && veq(vN,vN2) && veq(vT,vT2))
						bFound = true;
					else
						++j;
				}
				if (!bFound) ++t2;
			}

			assert(bFound);
			// if we found our own
			if (index2rec == index) { ++iNumUniqueVerts; }

			piTriList_in_and_out[offs] = index2rec;
		}
	}
}


MIKKT_FORCE_INLINE
void GenerateSharedVerticesIndexList(SMikkTState *pState)
{
	// Generate bounding box
	int32_t * piHashTable=NULL, * piHashCount=NULL, * piHashOffsets=NULL, * piHashCount2=NULL;
	STmpVert * pTmpVert = NULL;
	int32_t i=0, iChannel=0, k=0, e=0;
	int32_t iMaxCount=0;
	SVec3 vMin = GetPosition(pState->pCtx, 0), vMax = vMin, vDim;
	float fMin, fMax;
	for (i=1; i<(pState->iNrTrianglesIn*3); i++)
	{
		const int32_t index = pState->piTriListIn[i];

		const SVec3 vP = GetPosition(pState->pCtx, index);
		if (vMin.x > vP.x) vMin.x = vP.x;
		else if (vMax.x < vP.x) vMax.x = vP.x;
		if (vMin.y > vP.y) vMin.y = vP.y;
		else if (vMax.y < vP.y) vMax.y = vP.y;
		if (vMin.z > vP.z) vMin.z = vP.z;
		else if (vMax.z < vP.z) vMax.z = vP.z;
	}

	vDim = vsub(vMax,vMin);
	iChannel = 0;
	fMin = vMin.x; fMax=vMax.x;
	if (vDim.y>vDim.x && vDim.y>vDim.z)
	{
		iChannel=1;
		fMin = vMin.y;
		fMax = vMax.y;
	}
	else if (vDim.z>vDim.x)
	{
		iChannel=2;
		fMin = vMin.z;
		fMax = vMax.z;
	}

	// make allocations
	piHashTable = (int32_t *)   pState->pCtx->alloc.fpMalloc(sizeof(int32_t)*pState->iNrTrianglesIn*3);
	piHashCount = (int32_t *)   pState->pCtx->alloc.fpMalloc(sizeof(int32_t)*MIKKT_CELLS);
	piHashOffsets = (int32_t *) pState->pCtx->alloc.fpMalloc(sizeof(int32_t)*MIKKT_CELLS);
	piHashCount2 = (int32_t *)  pState->pCtx->alloc.fpMalloc(sizeof(int32_t)*MIKKT_CELLS);

	if (piHashTable==NULL || piHashCount==NULL || piHashOffsets==NULL || piHashCount2==NULL)
	{
		if (piHashTable!=NULL) pState->pCtx->alloc.fpFree(piHashTable);
		if (piHashCount!=NULL) pState->pCtx->alloc.fpFree(piHashCount);
		if (piHashOffsets!=NULL) pState->pCtx->alloc.fpFree(piHashOffsets);
		if (piHashCount2!=NULL) pState->pCtx->alloc.fpFree(piHashCount2);
		GenerateSharedVerticesIndexListSlow(pState->piTriListIn, pState->pCtx, pState->iNrTrianglesIn);
		return;
	}
	memset(piHashCount, 0, sizeof(int32_t)*MIKKT_CELLS);
	memset(piHashCount2, 0, sizeof(int32_t)*MIKKT_CELLS);

	// count amount of elements in each cell unit
	for (i=0; i<(pState->iNrTrianglesIn*3); i++)
	{
		const int32_t index = pState->piTriListIn[i];
		const SVec3 vP = GetPosition(pState->pCtx, index);
		const float fVal = iChannel==0 ? vP.x : (iChannel==1 ? vP.y : vP.z);
		const int32_t iCell = FindGridCell(fMin, fMax, fVal);
		++piHashCount[iCell];
	}

	// evaluate start index of each cell.
	piHashOffsets[0]=0;
	for (k=1; k<MIKKT_CELLS; k++)
		piHashOffsets[k]=piHashOffsets[k-1]+piHashCount[k-1];

	// insert vertices
	for (i=0; i<(pState->iNrTrianglesIn*3); i++)
	{
		const int32_t index = pState->piTriListIn[i];
		const SVec3 vP = GetPosition(pState->pCtx, index);
		const float fVal = iChannel==0 ? vP.x : (iChannel==1 ? vP.y : vP.z);
		const int32_t iCell = FindGridCell(fMin, fMax, fVal);
		int32_t * pTable = NULL;

		assert(piHashCount2[iCell]<piHashCount[iCell]);
		pTable = &piHashTable[piHashOffsets[iCell]];
		pTable[piHashCount2[iCell]] = i;	// vertex i has been inserted.
		++piHashCount2[iCell];
	}
	for (k=0; k<MIKKT_CELLS; k++)
		assert(piHashCount2[k] == piHashCount[k]);	// verify the count
	pState->pCtx->alloc.fpFree(piHashCount2);

	// find maximum amount of entries in any hash entry
	iMaxCount = piHashCount[0];
	for (k=1; k<MIKKT_CELLS; k++)
		if (iMaxCount<piHashCount[k])
			iMaxCount=piHashCount[k];
	pTmpVert = (STmpVert *) pState->pCtx->alloc.fpMalloc(sizeof(STmpVert)*iMaxCount);
	

	// complete the merge
	for (k=0; k<MIKKT_CELLS; k++)
	{
		// extract table of cell k and amount of entries in it
		int32_t * pTable = &piHashTable[piHashOffsets[k]];
		const int32_t iEntries = piHashCount[k];
		if (iEntries < 2) continue;

		if (pTmpVert!=NULL)
		{
			for (e=0; e<iEntries; e++)
			{
				int32_t i = pTable[e];
				const SVec3 vP = GetPosition(pState->pCtx, pState->piTriListIn[i]);
				pTmpVert[e].vert[0] = vP.x; pTmpVert[e].vert[1] = vP.y;
				pTmpVert[e].vert[2] = vP.z; pTmpVert[e].index = i;
			}
			MergeVertsFast(pState->piTriListIn, pTmpVert, pState->pCtx, 0, iEntries-1);
		}
		else
			MergeVertsSlow(pState->piTriListIn, pState->pCtx, pTable, iEntries);
	}

	if (pTmpVert!=NULL) { pState->pCtx->alloc.fpFree(pTmpVert); }
	pState->pCtx->alloc.fpFree(piHashTable);
	pState->pCtx->alloc.fpFree(piHashCount);
	pState->pCtx->alloc.fpFree(piHashOffsets);
}


// Mark all degenerate triangles
MIKKT_FORCE_INLINE
void markDegenTris(SMikkTState *pState) {
	pState->iTotTris = pState->iNrTrianglesIn;
	pState->iDegenTriangles = 0;
	for (int32_t i = 0; i < pState->iTotTris; ++i) {
		const int32_t i0 = pState->piTriListIn[i * 3 + 0];
		const int32_t i1 = pState->piTriListIn[i * 3 + 1];
		const int32_t i2 = pState->piTriListIn[i * 3 + 2];
		const SVec3 p0 = GetPosition(pState->pCtx, i0);
		const SVec3 p1 = GetPosition(pState->pCtx, i1);
		const SVec3 p2 = GetPosition(pState->pCtx, i2);
		if (veq(p0,p1) || veq(p0,p2) || veq(p1,p2)) {
			// degenerate
			pState->pTriInfos[i].iFlag |= MARK_DEGENERATE;
			++pState->iDegenTriangles;
		}
	}
	pState->iNrTrianglesIn = pState->iTotTris - pState->iDegenTriangles;
}

MIKKT_FORCE_INLINE
void InitTriInfo(STriInfo pTriInfos[], const int32_t piTriListIn[], const SMikkTSpaceContext * pContext, const int32_t iNrTrianglesIn)
{
	int32_t f=0, i=0, t=0;
	// pTriInfos[f].iFlag is cleared in GenerateInitialVerticesIndexList() which is called before this function.

	// generate neighbor info list
	for (f=0; f<iNrTrianglesIn; f++)
		for (i=0; i<3; i++)
		{
			pTriInfos[f].FaceNeighbors[i] = -1;
			pTriInfos[f].AssignedGroup[i] = NULL;

			pTriInfos[f].vOs.x=0.0f; pTriInfos[f].vOs.y=0.0f; pTriInfos[f].vOs.z=0.0f;
			pTriInfos[f].vOt.x=0.0f; pTriInfos[f].vOt.y=0.0f; pTriInfos[f].vOt.z=0.0f;
			pTriInfos[f].fMagS = 0;
			pTriInfos[f].fMagT = 0;

			// assumed bad
			pTriInfos[f].iFlag |= GROUP_WITH_ANY;
		}

	// evaluate first order derivatives
	for (f=0; f<iNrTrianglesIn; f++)
	{
		// initial values
		const SVec3 v1 = GetPosition(pContext, piTriListIn[f*3+0]);
		const SVec3 v2 = GetPosition(pContext, piTriListIn[f*3+1]);
		const SVec3 v3 = GetPosition(pContext, piTriListIn[f*3+2]);
		const SVec3 t1 = GetTexCoord(pContext, piTriListIn[f*3+0]);
		const SVec3 t2 = GetTexCoord(pContext, piTriListIn[f*3+1]);
		const SVec3 t3 = GetTexCoord(pContext, piTriListIn[f*3+2]);

		const float t21x = t2.x-t1.x;
		const float t21y = t2.y-t1.y;
		const float t31x = t3.x-t1.x;
		const float t31y = t3.y-t1.y;
		const SVec3 d1 = vsub(v2,v1);
		const SVec3 d2 = vsub(v3,v1);

		const float fSignedAreaSTx2 = t21x*t31y - t21y*t31x;
		//assert(fSignedAreaSTx2!=0);
		SVec3 vOs = vsub(vscale(t31y,d1), vscale(t21y,d2));	// eq 18
		SVec3 vOt = vadd(vscale(-t31x,d1), vscale(t21x,d2)); // eq 19

		pTriInfos[f].iFlag |= (fSignedAreaSTx2>0 ? ORIENT_PRESERVING : 0);

		if ( NotZero(fSignedAreaSTx2) )
		{
			const float fAbsArea = fabsf(fSignedAreaSTx2);
			const float fLenOs = Length(vOs);
			const float fLenOt = Length(vOt);
			const float fS = (pTriInfos[f].iFlag&ORIENT_PRESERVING)==0 ? (-1.0f) : 1.0f;
			if ( NotZero(fLenOs) ) pTriInfos[f].vOs = vscale(fS/fLenOs, vOs);
			if ( NotZero(fLenOt) ) pTriInfos[f].vOt = vscale(fS/fLenOt, vOt);

			// evaluate magnitudes prior to normalization of vOs and vOt
			pTriInfos[f].fMagS = fLenOs / fAbsArea;
			pTriInfos[f].fMagT = fLenOt / fAbsArea;

			// if this is a good triangle
			if ( NotZero(pTriInfos[f].fMagS) && NotZero(pTriInfos[f].fMagT))
				pTriInfos[f].iFlag &= (~GROUP_WITH_ANY);
		}
	}

	// force otherwise healthy quads to a fixed orientation
	while (t<(iNrTrianglesIn-1))
	{
		const int32_t iFO_a = pTriInfos[t].iOrgFaceNumber;
		const int32_t iFO_b = pTriInfos[t+1].iOrgFaceNumber;
		if (iFO_a==iFO_b)	// this is a quad
		{
			const bool bIsDeg_a = (pTriInfos[t].iFlag&MARK_DEGENERATE)!=0 ? true : false;
			const bool bIsDeg_b = (pTriInfos[t+1].iFlag&MARK_DEGENERATE)!=0 ? true : false;
			
			// bad triangles should already have been removed by
			// DegenPrologue(), but just in case check bIsDeg_a and bIsDeg_a are false
			if ((bIsDeg_a||bIsDeg_b)==false)
			{
				const bool bOrientA = (pTriInfos[t].iFlag&ORIENT_PRESERVING)!=0 ? true : false;
				const bool bOrientB = (pTriInfos[t+1].iFlag&ORIENT_PRESERVING)!=0 ? true : false;
				// if this happens the quad has extremely bad mapping!!
				if (bOrientA!=bOrientB)
				{
					bool bChooseOrientFirstTri = false;
					if ((pTriInfos[t+1].iFlag&GROUP_WITH_ANY)!=0) bChooseOrientFirstTri = true;
					else if ( CalcTexArea(pContext, &piTriListIn[t*3+0]) >= CalcTexArea(pContext, &piTriListIn[(t+1)*3+0]) )
						bChooseOrientFirstTri = true;

					// force match
					{
						const int32_t t0 = bChooseOrientFirstTri ? t : (t+1);
						const int32_t t1 = bChooseOrientFirstTri ? (t+1) : t;
						pTriInfos[t1].iFlag &= (~ORIENT_PRESERVING);	// clear first
						pTriInfos[t1].iFlag |= (pTriInfos[t0].iFlag&ORIENT_PRESERVING);	// copy bit
					}
				}
			}
			t += 2;
		}
		else
			++t;
	}
	
	// match up edge pairs
	{
		SEdge * pEdges = (SEdge *) pContext->alloc.fpMalloc(sizeof(SEdge)*iNrTrianglesIn*3);
		if (pEdges==NULL)
			BuildNeighborsSlow(pTriInfos, piTriListIn, iNrTrianglesIn);
		else
		{
			BuildNeighborsFast(pTriInfos, pEdges, piTriListIn, iNrTrianglesIn);
	
			pContext->alloc.fpFree(pEdges);
		}
	}
}

MIKKT_FORCE_INLINE
void DegenEpilogue(STSpace psTspace[], STriInfo pTriInfos[], int32_t piTriListIn[], const SMikkTSpaceContext * pContext, const int32_t iNrTrianglesIn, const int32_t iTotTris)
{
	int32_t t=0, i=0;
	// deal with degenerate triangles
	// punishment for degenerate triangles is O(N^2)
	for (t=iNrTrianglesIn; t<iTotTris; t++)
	{
		// degenerate triangles on a quad with one good triangle are skipped
		// here but processed in the next loop
		const bool bSkip = (pTriInfos[t].iFlag&QUAD_ONE_DEGEN_TRI)!=0 ? true : false;

		if (!bSkip)
		{
			for (i=0; i<3; i++)
			{
				const int32_t index1 = piTriListIn[t*3+i];
				// search through the good triangles
				bool bNotFound = true;
				int32_t j=0;
				while (bNotFound && j<(3*iNrTrianglesIn))
				{
					const int32_t index2 = piTriListIn[j];
					if (index1==index2) bNotFound=false;
					else ++j;
				}

				if (!bNotFound)
				{
					const int32_t iTri = j/3;
					const int32_t iVert = j%3;
					const int32_t iSrcVert=pTriInfos[iTri].vert_num[iVert];
					const int32_t iSrcOffs=pTriInfos[iTri].iTSpacesOffs;
					const int32_t iDstVert=pTriInfos[t].vert_num[i];
					const int32_t iDstOffs=pTriInfos[t].iTSpacesOffs;
					
					// copy tspace
					psTspace[iDstOffs+iDstVert] = psTspace[iSrcOffs+iSrcVert];
				}
			}
		}
	}

	// deal with degenerate quads with one good triangle
	for (t=0; t<iNrTrianglesIn; t++)
	{
		// this triangle belongs to a quad where the
		// other triangle is degenerate
		if ( (pTriInfos[t].iFlag&QUAD_ONE_DEGEN_TRI)!=0 )
		{
			SVec3 vDstP;
			int32_t iOrgF=-1, i=0;
			bool bNotFound;
			unsigned char * pV = pTriInfos[t].vert_num;
			int32_t iFlag = (1<<pV[0]) | (1<<pV[1]) | (1<<pV[2]);
			int32_t iMissingIndex = 0;
			if ((iFlag&2)==0) iMissingIndex=1;
			else if ((iFlag&4)==0) iMissingIndex=2;
			else if ((iFlag&8)==0) iMissingIndex=3;

			iOrgF = pTriInfos[t].iOrgFaceNumber;
			vDstP = GetPosition(pContext, MakeIndex(iOrgF, iMissingIndex));
			bNotFound = true;
			i=0;
			while (bNotFound && i<3)
			{
				const int32_t iVert = pV[i];
				const SVec3 vSrcP = GetPosition(pContext, MakeIndex(iOrgF, iVert));
				if (veq(vSrcP, vDstP)==true)
				{
					const int32_t iOffs = pTriInfos[t].iTSpacesOffs;
					psTspace[iOffs+iMissingIndex] = psTspace[iOffs+iVert];
					bNotFound=false;
				}
				else
					++i;
			}
			assert(!bNotFound);
		}
	}
}

// set data
MIKKT_FORCE_INLINE
void setData(const SMikkTSpaceContext *pCtx, GenTSpaceArgs *pJobArgs) {
	PixalcLinAllocIter iter = {0};
	pixalcLinAllocIterInit(&pJobArgs->alloc, (PixtyRange){0, INT32_MAX}, &iter);
	for (; !pixalcLinAllocIterAtEnd(&iter); pixalcLinAllocIterInc(&iter)) {
		const STSpace * pTSpace = pixalcLinAllocGetItem(&iter);
		float tang[] = {pTSpace->vOs.x, pTSpace->vOs.y, pTSpace->vOs.z};
		if (pCtx->m_pInterface->m_setTSpace!=NULL) {
			float bitang[] = {pTSpace->vOt.x, pTSpace->vOt.y, pTSpace->vOt.z};
			pCtx->m_pInterface->m_setTSpace(
				pCtx,
				tang,
				bitang,
				pTSpace->fMagS,
				pTSpace->fMagT,
				pTSpace->bOrient,
				pTSpace->face, pTSpace->corner
			);
		}
		if (pCtx->m_pInterface->m_setTSpaceBasic!=NULL) {
			pCtx->m_pInterface->m_setTSpaceBasic(
				pCtx,
				tang,
				pTSpace->bOrient ? 1.0f : -1.0f,
				pTSpace->face, pTSpace->corner
			);
		}
	}
}

static inline
int32_t getTSpaceCountForGroup(PixtyRange range, const SGroup *pGroups) {
	int32_t count = 0;
	for (int32_t i = range.start; i < range.end; ++i) {
		count += pGroups[i].iNrFaces * 3; //iNrFaces == number of tris
	}
	return count;
}

//thread safe!
MIKKT_FORCE_INLINE
PixErr genTangSpace(const SMikkTSpaceContext * pCtx, void *pThreadPool, const float fAngularThreshold) {
	PixErr err = PIX_ERR_SUCCESS;
	PIX_ERR_THROW_IFNOT_COND(
		err,
		pCtx->alloc.fpMalloc &&
		pCtx->alloc.fpCalloc &&
		pCtx->alloc.fpRealloc &&
		pCtx->alloc.fpFree,
		"",
		0
	);
	PIX_ERR_THROW_IFNOT_COND(
		err,
		pCtx->m_pInterface->m_getNumFaces &&
		pCtx->m_pInterface->m_getNumVerticesOfFace &&
		pCtx->m_pInterface->m_getPosition &&
		pCtx->m_pInterface->m_getNormal &&
		pCtx->m_pInterface->m_getTexCoord &&
		(pCtx->m_pInterface->m_setTSpace || pCtx->m_pInterface->m_setTSpaceBasic),
		"",
		0
	);
	// count nr_triangles
	SMikkTState state = {
		.pCtx = pCtx,
		.iNrFaces = pCtx->m_pInterface->m_getNumFaces(pCtx),
		.fThresCos = (float) cos((fAngularThreshold*(float)M_PI)/180.0f)
	};
	countTris(&state);
	PIX_ERR_THROW_IFNOT_COND(err, state.iNrTrianglesIn > 0, "", 0);

	// allocate memory for an index list
	state.piTriListIn = (int32_t *) pCtx->alloc.fpMalloc(sizeof(int32_t) * 3 * state.iNrTrianglesIn);
	state.pTriInfos = (STriInfo *) pCtx->alloc.fpMalloc(sizeof(STriInfo) * state.iNrTrianglesIn);
	PIX_ERR_THROW_IFNOT_COND(err, state.piTriListIn && state.pTriInfos,"", 0);
	
	// make an initial triangle --> face index list
	int32_t iNrTSPaces = GenerateInitialVerticesIndexList(&state);
	// make a welded index list of identical positions and attributes (pos, norm, texc)
	GenerateSharedVerticesIndexList(&state);
	markDegenTris(&state);
	// mark all triangle pairs that belong to a quad with only one
	// good triangle. These need special treatment in DegenEpilogue().
	// Additionally, move all good triangles to the start of
	// pTriInfos[] and piTriListIn[] without changing order and
	// put the degenerate triangles last.
	DegenPrologue(state.pTriInfos, state.piTriListIn, state.iNrTrianglesIn, state.iTotTris);
	// evaluate triangle level attributes and neighbor list
	InitTriInfo(state.pTriInfos, state.piTriListIn, pCtx, state.iNrTrianglesIn);
	
	if (!makeGroups(&state)) {
		PIX_ERR_THROW(err, "", 0);
	}

	// make tspaces, each group is split up into subgroups if necessary
	// based on fAngularThreshold. Finally a tangent space is made for
	// every resulting subgroup
	
	GenTSpaceArgs jobArgs[PIX_THREAD_MAX_THREADS] = {0};
	int32_t jobCount = 4;
	if (jobCount > PIX_THREAD_MAX_THREADS) {
		jobCount = PIX_THREAD_MAX_THREADS;
	}
	int32_t perJob = state.iNrActiveGroups / jobCount;
	if (!perJob) {
		jobCount = state.iNrActiveGroups;
		perJob = 1;
	}
	
	GenTSpaceArgsShared shared = {
		.fThresCos = state.fThresCos,
		.iNrActiveGroups = state.iNrActiveGroups,
		.pContext = pCtx,
		.pGroups = state.pGroups,
		.piTriListIn = state.piTriListIn,
		.pTriInfos = state.pTriInfos
	};
	/*
	GenTSpaceArgs jobArgs = {
		.pShared = &shared,
		.range = {0, state.iNrActiveGroups},
	};
	pixalcLinAllocInit(&pCtx->alloc, &jobArgs.alloc, sizeof(STSpace), 256, false);
	err = genTSpacesForJob(&jobArgs);
	PIX_ERR_THROW_IFNOT(err, "", 0);
	*/
	
	void *argPtrs[PIX_THREAD_MAX_THREADS] = {0};
	for (int32_t i = 0; i < jobCount; ++i) {
		GenTSpaceArgs *pEntry = jobArgs + i;
		argPtrs[i] = pEntry;
		pEntry->pShared = &shared;
		pEntry->range.start = perJob * i;
		pEntry->range.end = i == jobCount - 1 ?
			state.iNrActiveGroups : pEntry->range.start + perJob;
		pixalcLinAllocInit(&pCtx->alloc, &pEntry->alloc, sizeof(STSpace), 128, false);
	}
	void *jobHandles[PIX_THREAD_MAX_THREADS] = {0};
	err = pixthJobStackPushJobs(
		pThreadPool,
		jobCount,
		jobHandles,
		genTSpacesForJob,
		argPtrs
	);
	PIX_ERR_THROW_IFNOT(err, "", 0);
	bool done = false;
	err = pixthWaitForJobsIntern(pThreadPool, jobCount, jobHandles, true, &done);
	PIX_ERR_THROW_IFNOT(err, "", 0);
	
	
	// degenerate quads with one good triangle will be fixed by copying a space from
	// the good triangle to the coinciding vertex.
	// all other degenerate triangles will just copy a space from any good triangle
	// with the same welded index in piTriListIn[].
	/*
	TODO this should be able to go into genTSpace job
	DegenEpilogue(
		state.psTspace,
		state.pTriInfos,
		state.piTriListIn,
		pCtx,
		state.iNrTrianglesIn,
		state.iTotTris
	);
	*/

	for (int32_t i = 0; i < jobCount; ++i) {
		setData(pCtx, jobArgs + i);
		pixalcLinAllocDestroy(&jobArgs[i].alloc);
	}

	PIX_ERR_CATCH(0, err, ;);
	if (state.piTriListIn) {
		pCtx->alloc.fpFree(state.piTriListIn);
	}
	if (state.pTriInfos) {
		pCtx->alloc.fpFree(state.pTriInfos);
	}
	if (state.pGroups) {
		pCtx->alloc.fpFree(state.pGroups);
	}
	if (state.piGroupTrianglesBuffer) {
		pCtx->alloc.fpFree(state.piGroupTrianglesBuffer);
	}
	return err;
}

//thread safe!
// Default (recommended) fAngularThreshold is 180 degrees (which means threshold disabled)
MIKKT_FORCE_INLINE
bool genTangSpaceDefault(const SMikkTSpaceContext * pContext, void *pThreadPool)
{
	return genTangSpace(pContext, pThreadPool, 180.0f);
}


// To avoid visual errors (distortions/unwanted hard edges in lighting), when using sampled normal maps, the
// normal map sampler must use the exact inverse of the pixel shader transformation.
// The most efficient transformation we can possibly do in the pixel shader is
// achieved by using, directly, the "unnormalized" interpolated tangent, bitangent and vertex normal: vT, vB and vN.
// pixel shader (fast transform out)
// vNout = normalize( vNt.x * vT + vNt.y * vB + vNt.z * vN );
// where vNt is the tangent space normal. The normal map sampler must likewise use the
// interpolated and "unnormalized" tangent, bitangent and vertex normal to be compliant with the pixel shader.
// sampler does (exact inverse of pixel shader):
// float3 row0 = cross(vB, vN);
// float3 row1 = cross(vN, vT);
// float3 row2 = cross(vT, vB);
// float fSign = dot(vT, row0)<0 ? -1 : 1;
// vNt = normalize( fSign * float3(dot(vNout,row0), dot(vNout,row1), dot(vNout,row2)) );
// where vNout is the sampled normal in some chosen 3D space.
//
// Should you choose to reconstruct the bitangent in the pixel shader instead
// of the vertex shader, as explained earlier, then be sure to do this in the normal map sampler also.
// Finally, beware of quad triangulations. If the normal map sampler doesn't use the same triangulation of
// quads as your renderer then problems will occur since the interpolated tangent spaces will differ
// eventhough the vertex level tangent spaces match. This can be solved either by triangulating before
// sampling/exporting or by using the order-independent choice of diagonal for splitting quads suggested earlier.
// However, this must be used both by the sampler and your tools/rendering pipeline.

#ifdef __cplusplus
}
#endif
