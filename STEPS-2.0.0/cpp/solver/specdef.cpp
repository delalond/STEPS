////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2007-2013�Okinawa Institute of Science and Technology, Japan.
// Copyright (C) 2003-2006�University of Antwerp, Belgium.
//
// See the file AUTHORS for details.
//
// This file is part of STEPS.
//
// STEPS�is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// STEPS�is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.�If not, see <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////

/*
 *  Last Changed Rev:  $Rev: 489 $
 *  Last Changed Date: $Date: 2013-04-19 19:37:42 +0900 (Fri, 19 Apr 2013) $
 *  Last Changed By:   $Author: iain $
 */

// STL headers.
#include <string>
#include <cassert>

// STEPS headers.
#include "../common.h"
#include "types.hpp"
#include "../error.hpp"
#include "statedef.hpp"
#include "specdef.hpp"
#include "../model/spec.hpp"

////////////////////////////////////////////////////////////////////////////////

NAMESPACE_ALIAS(steps::solver, ssolver);

////////////////////////////////////////////////////////////////////////////////

ssolver::Specdef::Specdef(Statedef * sd, uint idx, steps::model::Spec * s)
: pStatedef(sd)
, pIdx(idx)
, pName()
, pSetupdone(false)
{
	assert(pStatedef != 0);
	assert(s != 0);
	pName = s->getID();

}

////////////////////////////////////////////////////////////////////////////////

ssolver::Specdef::~Specdef(void)
{

}

////////////////////////////////////////////////////////////////////////////////

void ssolver::Specdef::checkpoint(std::fstream & cp_file)
{
    // reserve
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::Specdef::restore(std::fstream & cp_file)
{
    // reserve
}

////////////////////////////////////////////////////////////////////////////////

std::string const ssolver::Specdef::name(void) const
{
	return pName;
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::Specdef::setup(void)
{

}

////////////////////////////////////////////////////////////////////////////////


// END

