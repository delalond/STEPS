/*
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
#include <steps/common.h>
#include <steps/solver/types.hpp>
#include <steps/error.hpp>
#include <steps/solver/statedef.hpp>
#include <steps/solver/chanstatedef.hpp>
#include <steps/model/chanstate.hpp>

////////////////////////////////////////////////////////////////////////////////

NAMESPACE_ALIAS(steps::solver, ssolver);

////////////////////////////////////////////////////////////////////////////////

ssolver::ChanStatedef::ChanStatedef(Statedef * sd, uint idx, steps::model::ChanState * cs)
: pStatedef(sd)
, pIdx(idx)
, pName()
, pChanState(0)
, pSetupdone(false)
{
	assert(pStatedef != 0);
	assert(cs != 0);
	pName = cs->getID();
	pChanState = cs;
														////// anything else????
}

////////////////////////////////////////////////////////////////////////////////

ssolver::ChanStatedef::~ChanStatedef(void)
{

}

////////////////////////////////////////////////////////////////////////////////

void ssolver::ChanStatedef::setup(void)
{

}

////////////////////////////////////////////////////////////////////////////////


// END
*/
