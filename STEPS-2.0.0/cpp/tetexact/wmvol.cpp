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
#include <cassert>
#include <cmath>
#include <algorithm>
#include <functional>
#include <iostream>

// STEPS headers.
#include "../common.h"
#include "../math/constants.hpp"
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

stex::WmVol::WmVol
(
    uint idx, solver::Compdef * cdef, double vol
)
: pIdx(idx)
, pCompdef(cdef)
, pVol(vol)
, pPoolCount(0)
, pPoolFlags(0)
, pKProcs()
, pNextTris()
{
    assert(pCompdef != 0);
	assert (pVol > 0.0);

    // Based on compartment definition, build other structures.
    uint nspecs = compdef()->countSpecs();
    pPoolCount = new uint[nspecs];
    pPoolFlags = new uint[nspecs];
    std::fill_n(pPoolCount, nspecs, 0);
    std::fill_n(pPoolFlags, nspecs, 0);
    pKProcs.resize(compdef()->countReacs());

}

////////////////////////////////////////////////////////////////////////////////

stex::WmVol::~WmVol(void)
{
    // Delete species pool information.
    delete[] pPoolCount;
    delete[] pPoolFlags;

    // Delete reaction rules.
    KProcPVecCI e = pKProcs.end();
    for (KProcPVecCI i = pKProcs.begin(); i != e; ++i) delete *i;
}

////////////////////////////////////////////////////////////////////////////////

void stex::WmVol::checkpoint(std::fstream & cp_file)
{
    uint nspecs = compdef()->countSpecs();
    cp_file.write((char*)pPoolCount, sizeof(uint) * nspecs);
    cp_file.write((char*)pPoolFlags, sizeof(uint) * nspecs);
}

////////////////////////////////////////////////////////////////////////////////

void stex::WmVol::restore(std::fstream & cp_file)
{
    uint nspecs = compdef()->countSpecs();
    cp_file.read((char*)pPoolCount, sizeof(uint) * nspecs);
    cp_file.read((char*)pPoolFlags, sizeof(uint) * nspecs);
}

////////////////////////////////////////////////////////////////////////////////

void stex::WmVol::setNextTri(stex::Tri * t)
{
	uint index = pNextTris.size();
	pNextTris.push_back(t);
}

////////////////////////////////////////////////////////////////////////////////

void stex::WmVol::setupKProcs(stex::Tetexact * tex)
{

    uint j = 0;

    // Note: ignoring diffusion KProcs

    // Create reaction kproc's.
    uint nreacs = compdef()->countReacs();
    for (uint i = 0; i < nreacs; ++i)
    {
        ssolver::Reacdef * rdef = compdef()->reacdef(i);
        stex::Reac * r = new stex::Reac(rdef, this);
        pKProcs[j++] = r;
        tex->addKProc(r);
    }

}

////////////////////////////////////////////////////////////////////////////////

void stex::WmVol::reset(void)
{
    uint nspecs = compdef()->countSpecs();
    std::fill_n(pPoolCount, nspecs, 0);
    std::fill_n(pPoolFlags, nspecs, 0);
    std::for_each(pKProcs.begin(), pKProcs.end(),
    		std::mem_fun(&stex::KProc::reset));
}

////////////////////////////////////////////////////////////////////////////////

double stex::WmVol::conc(uint gidx) const
{
	uint lspidx = compdef()->specG2L(gidx);
	double n = pPoolCount[lspidx];
	return (n/(1.0e3*pVol*steps::math::AVOGADRO));
}

////////////////////////////////////////////////////////////////////////////////

void stex::WmVol::setCount(uint lidx, uint count)
{
	assert (lidx < compdef()->countSpecs());
	uint oldcount = pPoolCount[lidx];
	pPoolCount[lidx] = count;

	/*
	// 16/01/10 IH: Counts now not stored in compartment object.
	// Now update the count in this tet's comp
	int diff = count - oldcount;
	double newcount = (compdef()->pools()[lidx]) + static_cast<double>(diff);
	// Compdef method will do the checking on the double argument
	// (should be positive or zero!)
	compdef()->setCount(lidx, newcount);
	*/
}

////////////////////////////////////////////////////////////////////////////////

void stex::WmVol::incCount(uint lidx, int inc)
{
	assert (lidx < compdef()->countSpecs());
	pPoolCount[lidx] += inc;
    assert(pPoolCount[lidx] >= 0);


	/* 16/01/10 IH: Counts now not stored in compartment object.
	compdef()->setCount(lidx, (compdef()->pools()[lidx] + static_cast<double>(inc)));
	assert(compdef()->pools()[lidx] >= 0.0);
	*/
}

////////////////////////////////////////////////////////////////////////////////

void stex::WmVol::setClamped(uint lidx, bool clamp)
{
    if (clamp == true) pPoolFlags[lidx] |= CLAMPED;
    else pPoolFlags[lidx] &= ~CLAMPED;
}

////////////////////////////////////////////////////////////////////////////////

stex::Reac * stex::WmVol::reac(uint lidx) const
{
    assert(lidx < compdef()->countReacs());
    return dynamic_cast<stex::Reac*>(pKProcs[lidx]);
}

////////////////////////////////////////////////////////////////////////////////

// END
