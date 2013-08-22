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

#ifndef STEPS_TETEXACT_SDIFF_HPP
#define STEPS_TETEXACT_SDIFF_HPP 1


// Standard library & STL headers.
#include <map>
#include <string>
#include <vector>
#include <fstream>

// STEPS headers.
#include "../common.h"
#include "../math/constants.hpp"
#include "../solver/surfdiffdef.hpp"
#include "kproc.hpp"
#include "tetexact.hpp"

////////////////////////////////////////////////////////////////////////////////

START_NAMESPACE(steps)
START_NAMESPACE(tetexact)

////////////////////////////////////////////////////////////////////////////////

// Forward declarations.
class Tri;

////////////////////////////////////////////////////////////////////////////////

class SDiff
: public steps::tetexact::KProc
{

public:

    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    SDiff(steps::solver::SurfDiffdef * sdef, steps::tetexact::Tri * tri);
    ~SDiff(void);

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

    inline steps::solver::SurfDiffdef * sdef(void) const
    { return pSDiffdef; }

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

    //inline steps::solver::Reacdef * defr(void) const
    //{ return pReacdef; }

    ////////////////////////////////////////////////////////////////////////

private:

    ////////////////////////////////////////////////////////////////////////

    uint                                ligGIdx;
    uint                                lidxTri;
    steps::solver::SurfDiffdef        * pSDiffdef;
    steps::tetexact::Tri              * pTri;
    std::vector<KProc*>                 pUpdVec[3];

    /*
    // Storing the species local index for each neighbouring tri: Needed
    // because neighbours may belong to different patches if we ever
    // implement diffusion boundary for surfaces
    // and therefore have different spec indices
    int 							  	pNeighbPatchLidx[3];

     */


    /// Properly scaled diffusivity constant.
    double                              pScaledDcst;
    // Compartmental dcst. Stored for convenience
    double                              pDcst;

    /// Used in selecting which directory the molecule should go.
    double                              pCDFSelector[2];

    ////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

END_NAMESPACE(tetexact)
END_NAMESPACE(steps)

////////////////////////////////////////////////////////////////////////////////

#endif // STEPS_TETEXACT_SDIFF_HPP
