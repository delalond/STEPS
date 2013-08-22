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

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// STL headers.
#include <cassert>
#include <sstream>
#include <string>

// STEPS headers.
#include "../common.h"
#include "../error.hpp"
#include "chan.hpp"
#include "chanstate.hpp"
#include "spec.hpp"

////////////////////////////////////////////////////////////////////////////////

USING_NAMESPACE(std);
USING_NAMESPACE(steps::model);

////////////////////////////////////////////////////////////////////////////////

ChanState::ChanState(string const & id, Model * model, Chan * chan)
: Spec(id, model)
, pChan(chan)
{
    if (pChan == 0)
    {
        ostringstream os;
        os << "No channel provided to ChanState initializer function";
        throw steps::ArgErr(os.str());
    }
    if (model != chan->getModel())
    {
    	ostringstream os;
    	os << "Channel is unknown in this model.";
    	throw steps::ArgErr(os.str());
    }

    pChan->_handleChanStateAdd(this);
}

////////////////////////////////////////////////////////////////////////////////

ChanState::~ChanState(void)
{
    if (pChan == 0) return;
	_handleSelfDelete();
}

////////////////////////////////////////////////////////////////////////////////

void ChanState::_handleSelfDelete(void)
{
	// Base method
	Spec::_handleSelfDelete();

    pChan->_handleChanStateDel(this);
    pChan = 0;
}

////////////////////////////////////////////////////////////////////////////////

void ChanState::setID(string const & id)
{
    assert(pChan != 0);
    if (id == getID()) return;
    // The following might raise an exception, e.g. if the new ID is not
    // valid or not unique. If this happens, we don't catch but simply let
    // it pass by into the Python layer.
    pChan->_handleChanStateIDChange(getID(), id);
    // This line will only be executed if the previous call didn't raise
    // an exception.
    Spec::setID(id);
}

////////////////////////////////////////////////////////////////////////////////

// END
