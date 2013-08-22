/*
 *
 * NOT NEEDED? BECAUSE CHANNEL STATESARE JUST LIKE SPECIES OBJECTS AT THIS LEVEL, THOUGH CHANDEF OBJECTS ARE AVAILABLE TO GROUP CHANNEL STATES
 *
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

#ifndef STEPS_SOLVER_CHANSTATEDEF_HPP
#define STEPS_SOLVER_CHANSTATEDEF_HPP 1

// STL headers.
#include <string>

// STEPS headers.
#include "../common.h"
#include "statedef.hpp"
#include "../model/chanstate.hpp"

////////////////////////////////////////////////////////////////////////////////

START_NAMESPACE(steps)
START_NAMESPACE(solver)

// Forwards declarations

////////////////////////////////////////////////////////////////////////////////

/// Defined Channel State
class ChanStatedef
{

public:
    /// Constructor
    ///
    /// \param sd State of the solver.
    /// \param idx Global index of the channel.
    /// \param cs Pointer to the assocaited ChanState object.
	ChanStatedef(Statedef * sd, uint idx, steps::model::ChanState * cs);

    /// Destructor
	~ChanStatedef(void);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: SPECIES
    ////////////////////////////////////////////////////////////////////////

	/// Return the global index of this species.
	inline uint gidx(void) const
	{ return pIdx; }

    /// Return the name of the species.
	inline std::string const name(void) const
	{ return pName; }

    ////////////////////////////////////////////////////////////////////////
    // SOLVER METHODS: SETUP
    ////////////////////////////////////////////////////////////////////////
    /// Setup the object.
    ///
	/// This method is included for consistency with other def objects,
	/// but currently does nothing.
	void setup(void);

    ////////////////////////////////////////////////////////////////////////

private:

    ////////////////////////////////////////////////////////////////////////

	Statedef                          * pStatedef;
	uint                                pIdx;
	std::string                         pName;
	bool								pSetupdone;
	// ChanState                         * pChanState;

    ////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

END_NAMESPACE(solver)
END_NAMESPACE(steps)

#endif
// STEPS_SOLVER_CHANSTATEDEF_HPP

// END

*/
