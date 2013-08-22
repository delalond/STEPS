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


/*
 *  Last Changed Rev:  $Rev: 489 $
 *  Last Changed Date: $Date: 2013-04-19 19:37:42 +0900 (Fri, 19 Apr 2013) $
 *  Last Changed By:   $Author: iain $
 */

// Standard library & STL headers.
#include <cassert>
#include <cmath>
#include <algorithm>
#include <functional>
#include <iostream>

// STEPS headers.
#include "../common.h"
#include "../solver/compdef.hpp"
#include "../solver/diffdef.hpp"
#include "../solver/reacdef.hpp"
#include "diff.hpp"
#include "reac.hpp"
#include "tet.hpp"
#include "tri.hpp"
#include "kproc.hpp"
#include "tetexact.hpp"
#include "wmvol.hpp"

////////////////////////////////////////////////////////////////////////////////

NAMESPACE_ALIAS(steps::tetexact, stex);
NAMESPACE_ALIAS(steps::solver, ssolver);

////////////////////////////////////////////////////////////////////////////////

stex::Tet::Tet
(
    uint idx, solver::Compdef * cdef, double vol,
    double a0, double a1, double a2, double a3,
    double d0, double d1, double d2, double d3,
    int tet0, int tet1, int tet2, int tet3
)
: WmVol(idx, cdef, vol)
, pTets()
//, pTris()
, pNextTet()
, pAreas()
, pDist()
{
	assert (a0 > 0.0 && a1 > 0.0 && a2 > 0.0 && a3 > 0.0);
    assert (d0 >= 0.0 && d1 >= 0.0 && d2 >= 0.0 && d3 >= 0.0);

    pNextTris.resize(4);

    // At this point we don't have neighbouring tet pointers,
    // but we can store their indices
    for (uint i=0; i <= 3; ++i)
    {
    	pNextTet[i] = 0;
    	pNextTris[i] = 0;
    }
    pTets[0] = tet0;
    pTets[1] = tet1;
    pTets[2] = tet2;
    pTets[3] = tet3;

    pAreas[0] = a0;
    pAreas[1] = a1;
    pAreas[2] = a2;
    pAreas[3] = a3;

    pDist[0] = d0;
    pDist[1] = d1;
    pDist[2] = d2;
    pDist[3] = d3;

    std::fill_n(pDiffBndDirection, 4, false);
    kprocs().resize(compdef()->countDiffs() + compdef()->countReacs());

}

////////////////////////////////////////////////////////////////////////////////

stex::Tet::~Tet(void)
{

}

////////////////////////////////////////////////////////////////////////////////

void stex::Tet::checkpoint(std::fstream & cp_file)
{
    cp_file.write((char*)pDiffBndDirection, sizeof(bool) * 4);
    WmVol::checkpoint(cp_file);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tet::restore(std::fstream & cp_file)
{
    cp_file.read((char*)pDiffBndDirection, sizeof(bool) * 4);
    WmVol::restore(cp_file);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tet::setNextTet(uint i, stex::Tet * t)
{
	/*
    if (t->compdef() != compdef())
    {
        pNextTet[i] = 0;
    }
    else
    {
        pNextTet[i] = t;
        if (pNextTri[i] != 0) std::cout << "WARNING: writing over nextTri index " << i;
        pNextTri[i] = 0;
    }
    */

	// Now adding all tets, even those from other compartments, due to the diffusion boundaries
    pNextTet[i] = t;

    //if (pNextTris[i] != 0) std::cout << "WARNING: writing over nextTri index " << i;
    pNextTris[i] = 0;

}

////////////////////////////////////////////////////////////////////////////////

void stex::Tet::setDiffBndDirection(uint i)
{
	assert(i < 4);

	pDiffBndDirection[i] = true;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tet::setNextTri(stex::Tri *t)
{
	assert(false);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tet::setNextTri(uint i, stex::Tri * t)
{
	assert (pNextTris.size() == 4);
	assert (i <= 3);

    // This is too common now to include this message- for any internal patch this happens
	//if (pNextTet[i] != 0) std::cout << "WARNING: writing over nextTet index " << i;

    pNextTet[i] = 0;
    pNextTris[i]= t;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tet::setupKProcs(stex::Tetexact * tex)
{
    uint j = 0;

    // Create reaction kproc's.
    uint nreacs = compdef()->countReacs();
    for (uint i = 0; i < nreacs; ++i)
    {
        ssolver::Reacdef * rdef = compdef()->reacdef(i);
        stex::Reac * r = new stex::Reac(rdef, this);
        kprocs()[j++] = r;
        tex->addKProc(r);
    }

    // Create diffusion kproc's.
    // NOTE: The order is important here- diffs should come after reacs,
    // because diffs will not be stored in WmVols and the Comp will call the
    // parent method often.
    uint ndiffs = compdef()->countDiffs();
    for (uint i = 0; i < ndiffs; ++i)
    {
        ssolver::Diffdef * ddef = compdef()->diffdef(i);
        stex::Diff * d = new stex::Diff(ddef, this);
        kprocs()[j++] = d;
        tex->addKProc(d);
    }
}

////////////////////////////////////////////////////////////////////////////////

stex::Diff * stex::Tet::diff(uint lidx) const
{
    assert(lidx < compdef()->countDiffs());
    return dynamic_cast<stex::Diff*>(pKProcs[compdef()->countReacs() + lidx]);
}

////////////////////////////////////////////////////////////////////////////////

// END
