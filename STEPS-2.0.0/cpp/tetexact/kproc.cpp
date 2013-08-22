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
#include <vector>
#include <cassert>

// STEPS headers.
#include "../common.h"
#include "kproc.hpp"

////////////////////////////////////////////////////////////////////////////////

NAMESPACE_ALIAS(steps::tetexact, stex);

////////////////////////////////////////////////////////////////////////////////

stex::KProc::KProc(void)
: rExtent(0)
, pFlags(0)
, pSchedIDX(0)
, crData()
{
}

////////////////////////////////////////////////////////////////////////////////

stex::KProc::~KProc(void)
{
}

////////////////////////////////////////////////////////////////////////////////

void stex::KProc::setActive(bool active)
{
    if (active == true) pFlags &= ~INACTIVATED;
    else pFlags |= INACTIVATED;
}

////////////////////////////////////////////////////////////////////////////////

uint stex::KProc::getExtent(void) const
{
	return rExtent;
}

////////////////////////////////////////////////////////////////////////////////

void stex::KProc::resetExtent(void)
{
	rExtent = 0;
}
////////////////////////////////////////////////////////////////////////////////

void stex::KProc::resetCcst(void) const
{
	// This should never get called on base object
	assert (false);
}

////////////////////////////////////////////////////////////////////////////////

double stex::KProc::c(void) const
{
    // Should never get called on base object
	assert (false);
}

////////////////////////////////////////////////////////////////////////////////

double stex::KProc::h(void)
{
    // Should never get called on base object
	assert (false);
}

////////////////////////////////////////////////////////////////////////////////

// END
