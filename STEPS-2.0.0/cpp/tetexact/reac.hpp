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

#ifndef STEPS_TETEXACT_REAC_HPP
#define STEPS_TETEXACT_REAC_HPP 1


// Standard library & STL headers.
#include <map>
#include <string>
#include <vector>
#include <fstream>

// STEPS headers.
#include "../common.h"
#include "../math/constants.hpp"
#include "../solver/reacdef.hpp"
#include "kproc.hpp"
//#include "tetexact.hpp"


////////////////////////////////////////////////////////////////////////////////

START_NAMESPACE(steps)
START_NAMESPACE(tetexact)

////////////////////////////////////////////////////////////////////////////////

// Forward declarations.
class WmVol;
class Tri;
class Tetexact;

////////////////////////////////////////////////////////////////////////////////

class Reac
: public steps::tetexact::KProc
{

public:

    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    Reac(steps::solver::Reacdef * rdef, steps::tetexact::WmVol * tet);
    ~Reac(void);

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream & cp_file);

    /// restore data
    void restore(std::fstream & cp_file);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS
    ////////////////////////////////////////////////////////////////////////

    double c(void) const
    { return pCcst; }
    void resetCcst(void);

    inline double kcst(void) const
    { return pKcst; }
    void setKcst(double k);

    double h(void)
    { return (rate()/pCcst); }

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

    steps::solver::Reacdef                              * pReacdef;
    steps::tetexact::WmVol                              * pTet;
    std::vector<KProc*>                                   pUpdVec;
    /// Properly scaled reaction constant.
    double                                                pCcst;
    // Also store the K constant for convenience
    double                                                pKcst;

    ////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

END_NAMESPACE(tetexact)
END_NAMESPACE(steps)

////////////////////////////////////////////////////////////////////////////////

#endif

// STEPS_TETEXACT_REAC_HPP