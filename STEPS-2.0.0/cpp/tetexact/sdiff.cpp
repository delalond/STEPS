////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2007-2013ÊOkinawa Institute of Science and Technology, Japan.
// Copyright (C) 2003-2006ÊUniversity of Antwerp, Belgium.
//
// See the file AUTHORS for details.
//
// This file is part of STEPS.
//
// STEPSÊis free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// STEPSÊis distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.ÊIf not, see <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////


// Standard library & STL headers.
#include <vector>

// STEPS headers.
#include "../common.h"
#include "../math/constants.hpp"
#include "../solver/surfdiffdef.hpp"
#include "../solver/patchdef.hpp"
#include "sdiff.hpp"
#include "tri.hpp"
#include "kproc.hpp"
#include "tetexact.hpp"

#include <iostream>

////////////////////////////////////////////////////////////////////////////////

NAMESPACE_ALIAS(steps::tetexact, stex);
NAMESPACE_ALIAS(steps::solver, ssolver);
NAMESPACE_ALIAS(steps::math, smath);

////////////////////////////////////////////////////////////////////////////////

stex::SDiff::SDiff(ssolver::SurfDiffdef * sdef, stex::Tri * tri)
: KProc()
, pSDiffdef(sdef)
, pTri(tri)
, pUpdVec()
, pScaledDcst(0.0)
, pDcst(0.0)
, pCDFSelector()

{
	assert(pSDiffdef != 0);
	assert(pTri != 0);
	stex::Tri * next[3] =
	{
		pTri->nextTri(0),
		pTri->nextTri(1),
		pTri->nextTri(2),
	};

    ligGIdx = pSDiffdef->lig();
    ssolver::Patchdef * pdef = pTri->patchdef();
    lidxTri = pdef->specG2L(ligGIdx);


    /*
    uint gidx = pSurfDiffdef->lig();
    for (uint i = 0; i < 4; ++i)
    {
        if (next[i] == 0)
        {
        	pNeighbCompLidx[i] = -1;
        	continue;
        }
        else
        {
        	pNeighbCompLidx[i] = next[i]->compdef()->specG2L(ligGIdx);
        }
    }
	*/

    // Precalculate part of the scaled diffusion constant.
	uint ldidx = pTri->patchdef()->surfdiffG2L(pSDiffdef->gidx());
	double dcst = pTri->patchdef()->dcst(ldidx);
    pDcst = dcst;

    double d[3] = { 0.0, 0.0, 0.0};
    for (uint i = 0; i < 3; ++i)
    {
        // Compute the scaled diffusion constant.
    	// Need to here check if the direction is a diffusion boundary
        double dist = pTri->dist(i);
        if ((dist > 0.0) && (next[i] != 0))
        {
        	if (next[i]->patchdef() == pdef)	d[i] = (pTri->length(i) * dcst) / (pTri->area() * dist);
        	else d[i] = 0;
        }
    }
    // Compute scaled "diffusion constant".
    for (uint i = 0; i < 3; ++i)
    {
		pScaledDcst += d[i];
	}

    // Should not be negative!
    assert(pScaledDcst >= 0);

    // Setup the selector distribution.
    if (pScaledDcst == 0.0)
    {
        pCDFSelector[0] = 0.0;
        pCDFSelector[1] = 0.0;
    }
    else
    {
        pCDFSelector[0] = d[0] / pScaledDcst;
        pCDFSelector[1] = pCDFSelector[0] + (d[1] / pScaledDcst);
    }
}

////////////////////////////////////////////////////////////////////////////////

stex::SDiff::~SDiff(void)
{
}

////////////////////////////////////////////////////////////////////////////////

void stex::SDiff::checkpoint(std::fstream & cp_file)
{
    cp_file.write((char*)&rExtent, sizeof(uint));
    cp_file.write((char*)&pFlags, sizeof(uint));

    cp_file.write((char*)&pScaledDcst, sizeof(double));
    cp_file.write((char*)&pDcst, sizeof(double));
    cp_file.write((char*)pCDFSelector, sizeof(double) * 2);

    cp_file.write((char*)&(crData.recorded), sizeof(bool));
    cp_file.write((char*)&(crData.pow), sizeof(int));
    cp_file.write((char*)&(crData.pos), sizeof(unsigned));
    cp_file.write((char*)&(crData.rate), sizeof(double));
}

////////////////////////////////////////////////////////////////////////////////

void stex::SDiff::restore(std::fstream & cp_file)
{
    cp_file.read((char*)&rExtent, sizeof(uint));
    cp_file.read((char*)&pFlags, sizeof(uint));

    cp_file.read((char*)&pScaledDcst, sizeof(double));
    cp_file.read((char*)&pDcst, sizeof(double));
    cp_file.read((char*)pCDFSelector, sizeof(double) * 2);


    cp_file.read((char*)&(crData.recorded), sizeof(bool));
    cp_file.read((char*)&(crData.pow), sizeof(int));
    cp_file.read((char*)&(crData.pos), sizeof(unsigned));
    cp_file.read((char*)&(crData.rate), sizeof(double));
}

////////////////////////////////////////////////////////////////////////////////

void stex::SDiff::setupDeps(void)
{
    // We will check all KProcs of the following simulation elements:
    //   * the 'source' triangle
    //   * any neighbouring tetrahedrons
    //
    // But also in the possible 'destination' triangles (leading to
    // four different dependency lists, each containing a copy of the
    // dependencies in the 'source' tet):
    //   * any neighbouring triangles
    //   * any neighbouring tetrahedrons of these neighbouring tris
    //


    // Search for dependencies in the 'source' triangle.
    std::set<stex::KProc*> local;

    KProcPVecCI kprocend = pTri->kprocEnd();
    for (KProcPVecCI k = pTri->kprocBegin(); k != kprocend; ++k)
    {
        // Check locally.
        if ((*k)->depSpecTri(ligGIdx, pTri) == true) {
            //local.insert((*k)->getSSARef());
            local.insert(*k);
        }
    }

    // Check the neighbouring tetrahedrons.
    stex::WmVol * itet = pTri->iTet();
    if (itet != 0)
    {
    	kprocend = itet->kprocEnd();
    	for (KProcPVecCI k = itet->kprocBegin(); k != kprocend; ++k)
    	{
    		if ((*k)->depSpecTri(ligGIdx, pTri) == true) {
    			//local.insert((*k)->getSSARef());
    			local.insert(*k);
    		}
    	}
    }

    stex::WmVol * otet = pTri->oTet();
    if (otet != 0)
    {
    	kprocend = otet->kprocEnd();
    	for (KProcPVecCI k = otet->kprocBegin(); k != kprocend; ++k)
    	{
    		if ((*k)->depSpecTri(ligGIdx, pTri) == true) {
    			//local.insert((*k)->getSSARef());
    			local.insert(*k);
    		}
    	}
    }


    // Search for dependencies in neighbouring triangles.
    for (uint i = 0; i < 3; ++i)
    {
        // Fetch next triangle, if it exists.
        stex::Tri * next = pTri->nextTri(i);
        if (next == 0) continue;

        // Copy local dependencies.
        std::set<stex::KProc*> local2(local.begin(), local.end());
        //std::set<stex::KProc*> local2_objs(local_objs.begin(), local_objs.end());

        // Find the ones 'locally' in the next tri.
        kprocend = next->kprocEnd();
        for (KProcPVecCI k = next->kprocBegin(); k != kprocend; ++k)
        {
            if ((*k)->depSpecTri(ligGIdx, next) == true) {
                //local.insert((*k)->getSSARef());
                local2.insert(*k);
            }
        }

        /*
        // Find deps in neighbouring tetrahedrons in the next tri.
        for (uint j = 0; j < 2; ++j)
        {
            // Fetch next tetrahedron, if it exists.
            stex::Tet * next2 = next->nextTet(j);
            if (next2 == 0) continue;

            // Find deps.
            kprocend = next2->kprocEnd();
            for (KProcPVecCI k = next2->kprocBegin(); k != kprocend; ++k)
            {
                if ((*k)->depSpecTri(ligGIdx, next) == true) {
                    //local.insert((*k)->getSSARef());
                    local2.insert(*k);
                }
            }
        }
        */
        // Fetch inner tetrahedron, if it exists.
        stex::WmVol * itet = next->iTet();
        if (itet != 0)
        {
        	// Find deps.
        	kprocend = itet->kprocEnd();
        	for (KProcPVecCI k = itet->kprocBegin(); k != kprocend; ++k)
        	{
        		if ((*k)->depSpecTri(ligGIdx, next) == true) {
        			//local.insert((*k)->getSSARef());
        			local2.insert(*k);
        		}
        	}
        }

        // Fetch outer tetrahedron, if it exists.
        stex::WmVol * otet = next->oTet();
        if (otet != 0)
        {
        	// Find deps.
        	kprocend = otet->kprocEnd();
        	for (KProcPVecCI k = otet->kprocBegin(); k != kprocend; ++k)
        	{
        		if ((*k)->depSpecTri(ligGIdx, next) == true) {
        			//local.insert((*k)->getSSARef());
        			local2.insert(*k);
        		}
        	}
        }

        // Copy the set to the update vector.
        pUpdVec[i].assign(local2.begin(), local2.end());
        //pUpdObjVec[i].assign(local2_objs.begin(), local2_objs.end());

    }

}

////////////////////////////////////////////////////////////////////////////////

bool stex::SDiff::depSpecTet(uint gidx, stex::WmVol * tet)
{
    return false;
}

////////////////////////////////////////////////////////////////////////////////

bool stex::SDiff::depSpecTri(uint gidx, stex::Tri * tri)
{
    if (pTri != tri) return false;
    if (gidx != ligGIdx) return false;
    return true;
}

////////////////////////////////////////////////////////////////////////////////

void stex::SDiff::reset(void)
{
    resetExtent();


    uint ldidx = pTri->patchdef()->surfdiffG2L(pSDiffdef->gidx());
	double dcst = pTri->patchdef()->dcst(ldidx);
    setDcst(dcst);

    setActive(true);

    crData.recorded = false;
    crData.pow = 0;
    crData.pos = 0;
    crData.rate = 0.0;

}

////////////////////////////////////////////////////////////////////////////////

void stex::SDiff::setDcst(double dcst)
{
	assert(dcst >= 0.0);
	pDcst = dcst;

	stex::Tri * next[3] =
	{
		pTri->nextTri(0),
		pTri->nextTri(1),
		pTri->nextTri(2),
	};

    double d[3] = { 0.0, 0.0, 0.0};

    for (uint i = 0; i < 3; ++i)
    {
        // Compute the scaled diffusion constant.
    	// Need to here check if the direction is a diffusion boundary
        double dist = pTri->dist(i);
        if ((dist > 0.0) && (next[i] != 0))
        {
    		if (next[i]->patchdef() == pTri->patchdef())	d[i] = (pTri->length(i) * dcst) / (pTri->area() * dist);
    		else d[i] = 0;
        }
    }

    // Compute scaled "diffusion constant".
    pScaledDcst = 0.0;

    for (uint i = 0; i < 3; ++i)
    {
    	pScaledDcst += d[i];
	}
    // Should not be negative!
    assert(pScaledDcst >= 0);

    // Setup the selector distribution.
    if (pScaledDcst == 0.0)
    {
        pCDFSelector[0] = 0.0;
        pCDFSelector[1] = 0.0;
    }
    else
    {
        pCDFSelector[0] = d[0] / pScaledDcst;
        pCDFSelector[1] = pCDFSelector[0] + (d[1] / pScaledDcst);
    }
}

////////////////////////////////////////////////////////////////////////////////

double stex::SDiff::rate(steps::tetexact::Tetexact * solver)
{
    if (inactive()) return 0.0;

    // Compute the rate.
    double rate = (pScaledDcst) * static_cast<double>(pTri->pools()[lidxTri]);
    assert(std::isnan(rate) == false);

    // Return.
    return rate;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<stex::KProc*> const & stex::SDiff::apply(steps::rng::RNG * rng, double dt, double simtime)
{
    //uint lidxTet = this->lidxTet;
    // Pre-fetch some general info.


    // Apply local change.
    uint * local = pTri->pools() + lidxTri;
	bool clamped = pTri->clamped(lidxTri);

    if (clamped == false)
    {
        assert(*local > 0);
    }

    // Apply change in next voxel: select a direction.
    double sel = rng->getUnfEE();

    if (sel < pCDFSelector[0])
    {
        // Direction 1.
        stex::Tri * nexttri = pTri->nextTri(0);
        // If there is no next tet 0, pCDFSelector[0] should be zero
        // So we can assert that nextet 0 does indeed exist
        assert (nexttri != 0);

        if (nexttri->clamped(lidxTri) == false)
        {
        	nexttri->incCount(lidxTri,1);
        }
        if (clamped == false) {pTri->incCount(lidxTri, -1); }

        rExtent++;

        return pUpdVec[0];
    }
    else if (sel < pCDFSelector[1])
    {
        // Direction 2.
        stex::Tri * nexttri = pTri->nextTri(1);
        // If there is no next tet 1, pCDFSelector[1] should be zero
        // So we can assert that nextet 1 does indeed exist
        assert (nexttri != 0);

        if (nexttri->clamped(lidxTri) == false)
        {
        	nexttri->incCount(lidxTri,1);
        }
        if (clamped == false) {pTri->incCount(lidxTri, -1); }

        rExtent++;

        return pUpdVec[1];
    }
    else
    {
        // Direction 3.
        stex::Tri * nexttri = pTri->nextTri(2);
        assert (nexttri != 0);

        if (nexttri->clamped(lidxTri) == false)
        {
        	nexttri->incCount(lidxTri,1);
        }
        if (clamped == false) {pTri->incCount(lidxTri, -1); }

        rExtent++;

        return pUpdVec[2];

    }

    // This should never happen!
    assert(0);
    std::cerr << "Cannot find a suitable direction for diffusion!\n";
    throw;
    return pUpdVec[0];
}

////////////////////////////////////////////////////////////////////////////////

uint stex::SDiff::updVecSize(void) const
{
	uint maxsize = pUpdVec[0].size();
	for (uint i=1; i <= 2; ++i)
	{
		if (pUpdVec[i].size() > maxsize) maxsize = pUpdVec[i].size();
	}
	return maxsize;
}

////////////////////////////////////////////////////////////////////////////////

// END
