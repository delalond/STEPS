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

#ifndef STEPS_TETEXACT_GHKCURR_HPP
#define STEPS_TETEXACT_GHKCURR_HPP 1

////////////////////////////////////////////////////////////////////////////////

// Standard library & STL headers.
#include <string>
#include <vector>

// STEPS headers.
#include "../common.h"
#include "../math/constants.hpp"
#include "../solver/ghkcurrdef.hpp"
#include "kproc.hpp"
//#include "tetexact.hpp"


////////////////////////////////////////////////////////////////////////////////

START_NAMESPACE(steps)
START_NAMESPACE(tetexact)

////////////////////////////////////////////////////////////////////////////////

// Forward declarations.
class Tri;
class Tetexact;

////////////////////////////////////////////////////////////////////////////////

class GHKcurr
: public steps::tetexact::KProc
{

public:

    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    GHKcurr(steps::solver::GHKcurrdef * ghkdef, steps::tetexact::Tri * tri);
    ~GHKcurr(void);

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream & cp_file);

    /// restore data
    void restore(std::fstream & cp_file);

    ////////////////////////////////////////////////////////////////////////
    // VIRTUAL INTERFACE METHODS
    ////////////////////////////////////////////////////////////////////////

    void setupDeps(void);
    bool depSpecTet(uint gidx, steps::tetexact::WmVol * tet);
    bool depSpecTri(uint gidx, steps::tetexact::Tri * tri);
    void reset(void);


    double rate(steps::tetexact::Tetexact * solver);

    // double rate(double v, double T);
    std::vector<KProc*> const & apply(steps::rng::RNG * rng, double dt, double simtime);

    inline bool efflux(void) const
    { return pEffFlux; }

    void setEffFlux(bool efx)
    { pEffFlux = efx; }

    uint updVecSize(void) const
    { return pUpdVec.size(); }

    ////////////////////////////////////////////////////////////////////////

private:

    ////////////////////////////////////////////////////////////////////////

    steps::solver::GHKcurrdef         * pGHKcurrdef;
    steps::tetexact::Tri              * pTri;
    std::vector<KProc*>                 pUpdVec;

    // Flag if flux is outward, positive flux (true) or inward, negative flux (false)
    bool								pEffFlux;

    ////////////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

END_NAMESPACE(tetexact)
END_NAMESPACE(steps)

////////////////////////////////////////////////////////////////////////////////

#endif

// STEPS_TETEXACT_GHKCURR_HPP
