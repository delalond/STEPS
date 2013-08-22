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

#ifndef STEPS_TETODE_PATCH_HPP
#define STEPS_TETODE_PATCH_HPP 1


// STL headers.
#include <cassert>
#include <vector>
#include <fstream>
#include <map>

// STEPS headers.
#include "../common.h"
#include "../solver/patchdef.hpp"
#include "tri.hpp"
#include "../solver/types.hpp"

////////////////////////////////////////////////////////////////////////////////

START_NAMESPACE(steps)
START_NAMESPACE(tetode)

////////////////////////////////////////////////////////////////////////////////

NAMESPACE_ALIAS(steps::tetode, stode);

////////////////////////////////////////////////////////////////////////////////

// Forward declarations.
class Patch;

// Auxiliary declarations.
typedef Patch *                         PatchP;
typedef std::vector<PatchP>             PatchPVec;
typedef PatchPVec::iterator             PatchPVecI;
typedef PatchPVec::const_iterator       PatchPVecCI;

////////////////////////////////////////////////////////////////////////////////

class Patch
{
public:

    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    Patch(steps::solver::Patchdef * patchdef);
    ~Patch(void);

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream & cp_file);

    /// restore data
    void restore(std::fstream & cp_file);

    /// Checks whether Tri::patchdef() corresponds to this object's
    /// PatchDef. There is no check whether the Tri object has already
    /// been added to this Patch object before (i.e. no duplicate
    /// checking).
    ///
    void addTri(stode::Tri * tri);

    ////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS
    ////////////////////////////////////////////////////////////////////////

    inline double area(void) const
    { return pArea; }

    // Return the local index of a tri given by global index
    uint getTri_GtoL(uint gidx);

    // Return the tri of a given local index
    Tri * getTri(uint lidx);

    inline steps::solver::Patchdef * def(void) const
    { return pPatchdef; }

    inline uint countTris(void) const
    { return pTris.size(); }


    inline TriPVecCI bgnTri(void) const
    { return pTris.begin(); }
    inline TriPVecCI endTri(void) const
    { return pTris.end(); }

    ////////////////////////////////////////////////////////////////////////

private:

    ////////////////////////////////////////////////////////////////////////

    steps::solver::Patchdef           * pPatchdef;

    TriPVec                             pTris;

    double 								pArea;

    // A map storing global index to local
    std::map<uint, uint>					pTris_GtoL;

};

////////////////////////////////////////////////////////////////////////////////

END_NAMESPACE(tetode)
END_NAMESPACE(steps)

#endif
// STEPS_TETODE_PATCH_HPP

// END
