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

#ifndef STEPS_TETEXACT_DIFF_HPP
#define STEPS_TETEXACT_DIFF_HPP 1


// Standard library & STL headers.
#include <map>
#include <string>
#include <vector>
#include <fstream>

// STEPS headers.
#include "../common.h"
#include "../math/constants.hpp"
#include "../solver/diffdef.hpp"
#include "kproc.hpp"
#include "tetexact.hpp"

////////////////////////////////////////////////////////////////////////////////

START_NAMESPACE(steps)
START_NAMESPACE(tetexact)

////////////////////////////////////////////////////////////////////////////////

// Forward declarations.
class Tet;
class Tri;
class WmVol;

////////////////////////////////////////////////////////////////////////////////

class Diff
: public steps::tetexact::KProc
{

public:

    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    Diff(steps::solver::Diffdef * ddef, steps::tetexact::Tet * tet);
    ~Diff(void);

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

    inline steps::solver::Diffdef * def(void) const
    { return pDiffdef; }

    inline double dcst(void) const
    { return pDcst; }
    void setDcst(double d);

    void setupDeps(void);
    bool depSpecTet(uint gidx, steps::tetexact::WmVol * tet);
    bool depSpecTri(uint gidx, steps::tetexact::Tri * tri);
    void reset(void);
    double rate(steps::tetexact::Tetexact * solver = 0);
    std::vector<KProc*> const & apply(steps::rng::RNG * rng, double dt, double simtime);

    uint updVecSize(void) const;

    ////////////////////////////////////////////////////////////////////////

    void setDiffBndActive(uint i, bool active);

    bool getDiffBndActive(uint i) const;

    ////////////////////////////////////////////////////////////////////////

    //inline steps::solver::Reacdef * defr(void) const
    //{ return pReacdef; }

    ////////////////////////////////////////////////////////////////////////

private:

    ////////////////////////////////////////////////////////////////////////

    uint                                ligGIdx;
    uint                                lidxTet;
    steps::solver::Diffdef            * pDiffdef;
    steps::tetexact::Tet              * pTet;
    std::vector<KProc*>                 pUpdVec[4];

    // Storing the species local index for each neighbouring tet: Needed
    // because neighbours may belong to different compartments
    // and therefore have different spec indices
    int 							  	pNeighbCompLidx[4];

    /// Properly scaled diffusivity constant.
    double                              pScaledDcst;
    // Compartmental dcst. Stored for convenience
    double                              pDcst;
    /// Used in selecting which directory the molecule should go.
    double                              pCDFSelector[3];

    // A flag to see if the species can move between compartments
    bool 								pDiffBndActive[4];

    // Flags to store if a direction is a diffusion boundary direction
    bool 							    pDiffBndDirection[4];

    ////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

END_NAMESPACE(tetexact)
END_NAMESPACE(steps)

////////////////////////////////////////////////////////////////////////////////

#endif // STEPS_TETEXACT_DIFF_HPP
