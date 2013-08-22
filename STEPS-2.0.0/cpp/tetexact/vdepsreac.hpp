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
 *  Last Changed Rev:  $Rev: 338 $
 *  Last Changed Date: $Date: 2010-05-18 12:48:52 +0900 (Tue, 18 May 2010) $
 *  Last Changed By:   $Author: iain $
 */

#ifndef STEPS_TETEXACT_VDEPSREAC_HPP
#define STEPS_TETEXACT_VDEPSREAC_HPP 1

////////////////////////////////////////////////////////////////////////////////

// Standard library & STL headers.
#include <map>
#include <string>
#include <vector>

// STEPS headers.
#include "../common.h"
#include "../math/constants.hpp"
#include "../solver/vdepsreacdef.hpp"
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

class VDepSReac
: public steps::tetexact::KProc
{

public:

    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    VDepSReac(steps::solver::VDepSReacdef * vdsrdef, steps::tetexact::Tri * tri);
    ~VDepSReac(void);

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

    double rate(steps::tetexact::Tetexact * solver = 0);
    std::vector<KProc*> const & apply(steps::rng::RNG * rng, double dt, double simtime);

    uint updVecSize(void) const
    { return pUpdVec.size(); }

    ////////////////////////////////////////////////////////////////////////

private:

    ////////////////////////////////////////////////////////////////////////

    steps::solver::VDepSReacdef       * pVDepSReacdef;
    steps::tetexact::Tri              * pTri;
    std::vector<KProc*>                   pUpdVec;

    // The information about the size of the comaprtment or patch, and the
    // dimensions. Important for scaling the constant.
    // As volumes and areas currently don't change this can be stored as
    // a constant.
    double 								pScaleFactor;

    ////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

END_NAMESPACE(tetexact)
END_NAMESPACE(steps)

////////////////////////////////////////////////////////////////////////////////

#endif

// STEPS_TETEXACT_VDEPSREAC_HPP
