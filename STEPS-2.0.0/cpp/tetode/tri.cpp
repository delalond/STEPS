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
#include "../solver/patchdef.hpp"

#include "tet.hpp"
#include "tri.hpp"
#include "../math/constants.hpp"

////////////////////////////////////////////////////////////////////////////////

NAMESPACE_ALIAS(steps::tetode, stode);
NAMESPACE_ALIAS(steps::solver, ssolver);

////////////////////////////////////////////////////////////////////////////////

stode::Tri::Tri(uint idx, steps::solver::Patchdef * patchdef, double area,
			double l0, double l1, double l2, double d0, double d1, double d2,
			int tetinner, int tetouter, int tri0, int tri1, int tri2)
: pIdx(idx)
, pPatchdef(patchdef)
, pArea(area)
, pLengths()
, pDist()
, pInnerTet(0)
, pOuterTet(0)
, pTets()
, pNextTri()
{
	assert(pPatchdef != 0);
	assert(pArea > 0.0);

	assert (l0 > 0.0 && l1 > 0.0 && l2 > 0.0);
    assert (d0 >= 0.0 && d1 >= 0.0 && d2 >= 0.0);

	pTets[0] = tetinner;
	pTets[1] = tetouter;

	pTris[0] = tri0;
	pTris[1] = tri1;
	pTris[2] = tri2;

	pNextTri[0] = 0;
	pNextTri[1] = 0;
	pNextTri[2] = 0;

	pLengths[0] = l0;
	pLengths[1] = l1;
	pLengths[2] = l2;

    pDist[0] = d0;
    pDist[1] = d1;
    pDist[2] = d2;

}

////////////////////////////////////////////////////////////////////////////////

stode::Tri::~Tri(void)
{

}

////////////////////////////////////////////////////////////////////////////////

void stode::Tri::setInnerTet(stode::Tet * t)
{
	pInnerTet = t;
}

////////////////////////////////////////////////////////////////////////////////

void stode::Tri::setOuterTet(stode::Tet * t)
{
	pOuterTet = t;
}

////////////////////////////////////////////////////////////////////////////////

void stode::Tri::setNextTri(uint i, stode::Tri * t)
{
	assert (i <= 2);

    pNextTri[i]= t;
}

////////////////////////////////////////////////////////////////////////////////

void stode::Tri::checkpoint(std::fstream & cp_file)
{
}

////////////////////////////////////////////////////////////////////////////////

void stode::Tri::restore(std::fstream & cp_file)
{
}

////////////////////////////////////////////////////////////////////////////////
//END
