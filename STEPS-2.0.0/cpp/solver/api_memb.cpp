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
#include <sstream>

// STEPS headers.
#include "../common.h"
#include "../error.hpp"
#include "api.hpp"
#include "statedef.hpp"

////////////////////////////////////////////////////////////////////////////////

USING(std, string);
USING_NAMESPACE(steps::solver);

////////////////////////////////////////////////////////////////////////////////

void API::setMembPotential(string const & m, double v)
{
	// the following may raise exceptions if string is unused
	uint midx = pStatedef->getMembIdx(m);

	_setMembPotential(midx, v);
}

////////////////////////////////////////////////////////////////////////////////

void API::setMembCapac(std::string const & m, double cm)
{
	// the following may raise exceptions if string is unused
	uint midx = pStatedef->getMembIdx(m);

	_setMembCapac(midx, cm);
}

////////////////////////////////////////////////////////////////////////////////

void API::setMembVolRes(std::string const & m, double ro)
{
	// the following may raise exceptions if string is unused
	uint midx = pStatedef->getMembIdx(m);

	_setMembVolRes(midx, ro);
}

////////////////////////////////////////////////////////////////////////////////

void API::setMembRes(std::string const & m, double ro, double vrev)
{
	// the following may raise exceptions if string is unused
	uint midx = pStatedef->getMembIdx(m);

	_setMembRes(midx, ro, vrev);
}

////////////////////////////////////////////////////////////////////////////////

void API::_setMembPotential(uint midx, double v)
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void API::_setMembCapac(uint midx, double cm)
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void API::_setMembVolRes(uint midx, double ro)
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void API::_setMembRes(uint midx, double ro, double vrev)
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

// END




