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


// STL headers.
#include <string>
#include <cassert>

// STEPS headers.
#include "../common.h"
#include "types.hpp"
#include "../error.hpp"
#include "statedef.hpp"
#include "specdef.hpp"
#include "surfdiffdef.hpp"
#include "../model/spec.hpp"

////////////////////////////////////////////////////////////////////////////////

NAMESPACE_ALIAS(steps::solver, ssolver);
NAMESPACE_ALIAS(steps::model, smod);

////////////////////////////////////////////////////////////////////////////////

ssolver::SurfDiffdef::SurfDiffdef(Statedef * sd, uint idx, steps::model::Diff * d)
: pStatedef(sd)
, pIdx(idx)
, pName()
, pDcst()
, pLig()
, pSetupdone(false)
, pSpec_DEP(0)
{
    assert(pStatedef != 0);
    assert(d != 0);

    pName = d->getID();
    pDcst = d->getDcst();
    pLig = d->getLig()->getID();

    uint nspecs = pStatedef->countSpecs();
    if (nspecs == 0) return;
    pSpec_DEP = new int[nspecs];
    std::fill_n(pSpec_DEP, nspecs, DEP_NONE);

}

////////////////////////////////////////////////////////////////////////////////

ssolver::SurfDiffdef::~SurfDiffdef(void)
{
	if (pStatedef->countSpecs() > 0) delete[] pSpec_DEP;
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::SurfDiffdef::checkpoint(std::fstream & cp_file)
{
    cp_file.write((char*)&pDcst, sizeof(double));
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::SurfDiffdef::restore(std::fstream & cp_file)
{
    cp_file.read((char*)&pDcst, sizeof(double));
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::SurfDiffdef::setup(void)
{
	assert (pSetupdone == false);

	pSpec_DEP[lig()] = DEP_STOICH;

	pSetupdone = true;
}

////////////////////////////////////////////////////////////////////////////////

std::string const ssolver::SurfDiffdef::name(void) const
{
	return pName;
}

////////////////////////////////////////////////////////////////////////////////

double ssolver::SurfDiffdef::dcst(void) const
{
	return pDcst;
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::SurfDiffdef::setDcst(double d)
{
	assert (d >= 0.0);
	pDcst = d;
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::SurfDiffdef::lig(void) const
{
	assert (pStatedef != 0);
	return pStatedef->getSpecIdx(pLig);
}

////////////////////////////////////////////////////////////////////////////////
/*
void ssolver::Diffdef::setLig(uint gidx)
{
	assert (gidx < pStatedef->countSpecs());
	ssolver::Specdef * spec = pStatedef->specdef(gidx);
	pLig = spec->name();
}
*/
////////////////////////////////////////////////////////////////////////////////

int ssolver::SurfDiffdef::dep(uint gidx) const
{
    assert(pSetupdone == true);
    assert(gidx < pStatedef->countSpecs());
    return pSpec_DEP[gidx];
}

////////////////////////////////////////////////////////////////////////////////

bool ssolver::SurfDiffdef::reqspec(uint gidx) const
{
    assert(pSetupdone == true);
    assert(gidx < pStatedef->countSpecs());
    if (pSpec_DEP[gidx] != DEP_NONE) return true;
    return false;
}

////////////////////////////////////////////////////////////////////////////////

// END
