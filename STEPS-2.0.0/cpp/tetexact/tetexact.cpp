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
#include <cmath>
#include <vector>
#include <cassert>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <limits>
#include <queue>
#include <fstream>
#include <iomanip>

// STEPS headers.
#include "../common.h"
#include "tetexact.hpp"
#include "kproc.hpp"
#include "tet.hpp"
#include "tri.hpp"
#include "reac.hpp"
#include "sreac.hpp"
#include "diff.hpp"
#include "comp.hpp"
#include "patch.hpp"
#include "wmvol.hpp"
#include "ghkcurr.hpp"
#include "vdeptrans.hpp"
#include "vdepsreac.hpp"
#include "diffboundary.hpp"
#include "../math/constants.hpp"
#include "../error.hpp"
#include "../solver/statedef.hpp"
#include "../solver/compdef.hpp"
#include "../solver/patchdef.hpp"
#include "../solver/reacdef.hpp"
#include "../solver/sreacdef.hpp"
#include "../solver/diffdef.hpp"
#include "../solver/chandef.hpp"
#include "../solver/ghkcurrdef.hpp"
#include "../solver/ohmiccurrdef.hpp"
#include "../solver/vdeptransdef.hpp"
#include "../solver/vdepsreacdef.hpp"
#include "../solver/types.hpp"
#include "../geom/tetmesh.hpp"
#include "../geom/tet.hpp"
#include "../geom/tri.hpp"

#include "../solver/efield/efield.hpp"

////////////////////////////////////////////////////////////////////////////////


NAMESPACE_ALIAS(steps::tetexact, stex);
NAMESPACE_ALIAS(steps::solver, ssolver);
NAMESPACE_ALIAS(steps::math, smath);

////////////////////////////////////////////////////////////////////////////////

/// Unary function that calls the array delete[] operator on pointers. Easy
/// to use with STL/Boost (see steps::tools::DeletePointer).
///
struct DeleteArray
{
    template <typename Type> void operator() (Type * pointer) const
    {
        delete[] pointer;
    }
};

////////////////////////////////////////////////////////////////////////////////

void stex::schedIDXSet_To_Vec(stex::SchedIDXSet const & s, stex::SchedIDXVec & v)
{
    v.resize(s.size());
    std::copy(s.begin(), s.end(), v.begin());
}

////////////////////////////////////////////////////////////////////////////////

stex::Tetexact::Tetexact(steps::model::Model * m, steps::wm::Geom * g, steps::rng::RNG * r,
						 bool calcMembPot)
: API(m, g, r)
, pMesh(0)
, pKProcs()
, pComps()
, pCompMap()
, pPatches()
, pDiffBoundaries()
, pTets()
, pTris()
, pWmVols()
, pA0(0.0)
//, pBuilt(false)
, pEFflag(calcMembPot)
, pEField(0)
, pTemp(0.0)
, pEFDT(1.0e-5)
, pEFNVerts(0)
, pEFVerts(0)
, pEFNTris(0)
, pEFTris(0)
, pEFTris_vec(0)
, pEFNTets(0)
, pEFTets(0)
, pEFVert_GtoL()
, pEFTri_GtoL()
, pEFTet_GtoL()
, pEFTri_LtoG()
{
    if (rng() == 0)
    {
        std::ostringstream os;
        os << "No RNG provided to solver initializer function";
        throw steps::ArgErr(os.str());
    }

	// All initialization code now in _setup() to allow EField solver to be
	// derived and create EField local objects within the constructor
    _setup();
}

////////////////////////////////////////////////////////////////////////////////

stex::Tetexact::~Tetexact(void)
{
    CompPVecCI comp_e = pComps.end();
    for (CompPVecCI c = pComps.begin(); c != comp_e; ++c) delete *c;
    PatchPVecCI patch_e = pPatches.end();
    for (PatchPVecCI p = pPatches.begin(); p != patch_e; ++p) delete *p;

    DiffBoundaryPVecCI db_e = pDiffBoundaries.end();
    for (DiffBoundaryPVecCI db = pDiffBoundaries.begin(); db!=db_e; ++db) delete *db;

    WmVolPVecCI wvol_e = pWmVols.end();
    for (WmVolPVecCI wvol = pWmVols.begin(); wvol != wvol_e; ++wvol)
    {
        if ((*wvol) != 0) delete (*wvol);
    }

    TetPVecCI tet_e = pTets.end();
    for (TetPVecCI t = pTets.begin(); t != tet_e; ++t)
    {
        if ((*t) != 0) delete (*t);
    }

    TriPVecCI tri_e = pTris.end();
    for (TriPVecCI t = pTris.end(); t != tri_e; ++t)
    {
        if ((*t) != 0) delete (*t);
    }

    uint ngroups = nGroups.size();
    for (uint i = 0; i < ngroups; i++) {
        free(nGroups[i]->indices);
        delete nGroups[i];
    }

    ngroups = pGroups.size();
    for (uint i = 0; i < ngroups; i++) {
        free(pGroups[i]->indices);
        delete pGroups[i];
    }

    if (efflag())
    {
        delete[] pEFVerts;
        delete[] pEFTris;
        delete[] pEFTets;
        delete pEField;
        delete[] pEFVert_GtoL;
        delete[] pEFTri_GtoL;
        delete[] pEFTet_GtoL;
        delete[] pEFTri_LtoG;
    }
}

///////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::checkpoint(std::string const & file_name)
{
    std::cout << "Checkpoint to " << file_name  << "...";
	std::fstream cp_file;

    cp_file.open(file_name.c_str(),
                std::fstream::out | std::fstream::binary | std::fstream::trunc);

	statedef()->checkpoint(cp_file);

    CompPVecCI comp_e = pComps.end();
    for (CompPVecCI c = pComps.begin(); c != comp_e; ++c) (*c)->checkpoint(cp_file);

    PatchPVecCI patch_e = pPatches.end();
    for (PatchPVecCI p = pPatches.begin(); p != patch_e; ++p) (*p)->checkpoint(cp_file);

    DiffBoundaryPVecCI db_e = pDiffBoundaries.end();
    for (DiffBoundaryPVecCI db = pDiffBoundaries.begin(); db != db_e; ++db) {
        (*db)->checkpoint(cp_file);
    }

    WmVolPVecCI wmv_e = pWmVols.end();
    for (WmVolPVecCI wmv = pWmVols.begin(); wmv != wmv_e; ++wmv)
    {
        if ((*wmv) != 0) {
            (*wmv)->checkpoint(cp_file);
        }
    }


    TetPVecCI tet_e = pTets.end();
    for (TetPVecCI t = pTets.begin(); t != tet_e; ++t)
    {
        if ((*t) != 0) {
        (*t)->checkpoint(cp_file);
        }
    }

    TriPVecCI tri_e = pTris.end();
    for (TriPVecCI t = pTris.begin(); t != tri_e; ++t)
    {
        if ((*t) != 0) {
            (*t)->checkpoint(cp_file);
        }
    }

    KProcPVecCI e = pKProcs.end();
    for (KProcPVecCI i = pKProcs.begin(); i != e; ++i) (*i)->checkpoint(cp_file);

    if (pEFflag) {
        cp_file.write((char*)&pTemp, sizeof(double));
        cp_file.write((char*)&pEFDT, sizeof(double));
        pEField->checkpoint(cp_file);
    }

    cp_file.write((char*)&nEntries, sizeof(uint));

    // checkpoint CR SSA

    cp_file.write((char*)&pSum, sizeof(double));
    cp_file.write((char*)&nSum, sizeof(double));
    cp_file.write((char*)&pA0, sizeof(double));

    uint n_ngroups = nGroups.size();
    uint n_pgroups = pGroups.size();

    cp_file.write((char*)&n_ngroups, sizeof(uint));
    cp_file.write((char*)&n_pgroups, sizeof(uint));

    for (uint i = 0; i < n_ngroups; i++) {
        CRGroup* group = nGroups[i];
        cp_file.write((char*)&(group->capacity), sizeof(unsigned));
        cp_file.write((char*)&(group->size), sizeof(unsigned));
        cp_file.write((char*)&(group->max), sizeof(double));
        cp_file.write((char*)&(group->sum), sizeof(double));

        for (uint j = 0; j < group->size; j++) {
            uint idx = group->indices[j]->schedIDX();
            cp_file.write((char*)&idx, sizeof(uint));
        }
    }

    for (uint i = 0; i < n_pgroups; i++) {
        CRGroup* group = pGroups[i];
        cp_file.write((char*)&(group->capacity), sizeof(unsigned));
        cp_file.write((char*)&(group->size), sizeof(unsigned));
        cp_file.write((char*)&(group->max), sizeof(double));
        cp_file.write((char*)&(group->sum), sizeof(double));

        for (uint j = 0; j < group->size; j++) {
            uint idx = group->indices[j]->schedIDX();
            cp_file.write((char*)&idx, sizeof(uint));
        }
    }

    cp_file.close();
    std::cout << "complete.\n";
}

///////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::restore(std::string const & file_name)
{
    std::cout << "Restore from " << file_name << "...";
	std::fstream cp_file;

    cp_file.open(file_name.c_str(),
                std::fstream::in | std::fstream::binary);

    cp_file.seekg(0);

	statedef()->restore(cp_file);

    CompPVecCI comp_e = pComps.end();
    for (CompPVecCI c = pComps.begin(); c != comp_e; ++c) (*c)->restore(cp_file);

    PatchPVecCI patch_e = pPatches.end();
    for (PatchPVecCI p = pPatches.begin(); p != patch_e; ++p) (*p)->restore(cp_file);

    DiffBoundaryPVecCI db_e = pDiffBoundaries.end();
    for (DiffBoundaryPVecCI db = pDiffBoundaries.begin(); db != db_e; ++db) {
        (*db)->restore(cp_file);
    }

    WmVolPVecCI wmv_e = pWmVols.end();
    for (WmVolPVecCI wmv = pWmVols.begin(); wmv != wmv_e; ++wmv)
    {
        if ((*wmv) != 0) {
            (*wmv)->restore(cp_file);
        }
    }

    TetPVecCI tet_e = pTets.end();
    for (TetPVecCI t = pTets.begin(); t != tet_e; ++t)
    {
        if ((*t) != 0) {
            (*t)->restore(cp_file);
        }
    }
    TriPVecCI tri_e = pTris.end();
    for (TriPVecCI t = pTris.begin(); t != tri_e; ++t)
    {
        if ((*t) != 0) {
            (*t)->restore(cp_file);
        }
    }

    KProcPVecCI e = pKProcs.end();
    for (KProcPVecCI i = pKProcs.begin(); i != e; ++i) (*i)->restore(cp_file);



    if (pEFflag) {
        cp_file.read((char*)&pTemp, sizeof(double));
        cp_file.read((char*)&pEFDT, sizeof(double));
        pEField->restore(cp_file);
    }

    uint stored_entries = 0;
    cp_file.read((char*)&stored_entries, sizeof(uint));

    if (stored_entries != nEntries) {
        std::ostringstream os;
		os << "Unknown Restore Error!";
		throw steps::ArgErr(os.str());
    }

    // restore CR SSA
    cp_file.read((char*)&pSum, sizeof(double));
    cp_file.read((char*)&nSum, sizeof(double));
    cp_file.read((char*)&pA0, sizeof(double));

    uint n_ngroups;
    uint n_pgroups;

    cp_file.read((char*)&n_ngroups, sizeof(uint));
    cp_file.read((char*)&n_pgroups, sizeof(uint));

    nGroups.resize(n_ngroups);
    pGroups.resize(n_pgroups);

    for (uint i = 0; i < n_ngroups; i++) {
        unsigned capacity;
        unsigned size;
        double max;
        double sum;

        cp_file.read((char*)&capacity, sizeof(unsigned));
        cp_file.read((char*)&size, sizeof(unsigned));
        cp_file.read((char*)&max, sizeof(double));
        cp_file.read((char*)&sum, sizeof(double));

        nGroups[i] = new CRGroup(0, capacity);
        nGroups[i]->size = size;
        nGroups[i]->max = max;
        nGroups[i]->sum = sum;

        for (uint j = 0; j < size; j++) {
            uint idx;
            cp_file.read((char*)&idx, sizeof(uint));
            nGroups[i]->indices[j] = pKProcs[idx];
        }
    }

    for (uint i = 0; i < n_pgroups; i++) {
        unsigned capacity;
        unsigned size;
        double max;
        double sum;

        cp_file.read((char*)&capacity, sizeof(unsigned));
        cp_file.read((char*)&size, sizeof(unsigned));
        cp_file.read((char*)&max, sizeof(double));
        cp_file.read((char*)&sum, sizeof(double));

        pGroups[i] = new CRGroup(0, capacity);
        pGroups[i]->size = size;
        pGroups[i]->max = max;
        pGroups[i]->sum = sum;

        for (uint j = 0; j < size; j++) {
            uint idx;
            cp_file.read((char*)&idx, sizeof(uint));
            pGroups[i]->indices[j] = pKProcs[idx];
        }
    }

    cp_file.close();

    std::cout << "complete.\n";
}

////////////////////////////////////////////////////////////////////////////////


std::string stex::Tetexact::getSolverName(void) const
{
	return "tetexact";
}

////////////////////////////////////////////////////////////////////////////////

std::string stex::Tetexact::getSolverDesc(void) const
{
	return "SSA Composition and Rejection Exact Method in tetrahedral mesh";
}

////////////////////////////////////////////////////////////////////////////////

std::string stex::Tetexact::getSolverAuthors(void) const
{
	return "Stefan Wils, Iain Hepburn, Weiliang Chen";
}

////////////////////////////////////////////////////////////////////////////////

std::string stex::Tetexact::getSolverEmail(void) const
{
	return "steps.dev@gmail.com";
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setup(void)
{
	// Perform upcast.
	if  (! (pMesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom())))
	{
		std::ostringstream os;
		os << "Geometry description to steps::solver::Tetexact solver";
		os << " constructor is not a valid steps::tetmesh::Tetmesh object.";
		throw steps::ArgErr(os.str());
	}

	// First initialise the pTets, pTris vector, because
	// want tets and tris to maintain indexing from Geometry
	uint ntets = mesh()->countTets();
	uint ntris = mesh()->countTris();
	uint ncomps = mesh()->_countComps();

	pTets.assign(ntets, NULL);
	pTris.assign(ntris, NULL);
	pWmVols.assign(ncomps, NULL);

	/*
     // Now check that every Tetrahedron has been added to a compartment.
     // This is a requirement for the forseeable future
     uint ntmcomps = mesh()->_countComps();
     uint tetcount = 0;

     for (uint c = 0; c < ntmcomps; ++c)
     {
     steps::wm::Comp * wmcomp = mesh()->_getComp(c);
     steps::tetmesh::TmComp * tmcomp = dynamic_cast<steps::tetmesh::TmComp*>(wmcomp);
     tetcount += tmcomp->countTets();
     }
     if (tetcount != ntets)
     {
     std::ostringstream os;
     os << "Cannot create solver object with this ";
     os << "steps.geom.Tetmesh geometry description object";
     os << " All tetrahedrons in Tetmesh must be assigned to a";
     os << " steps.geom.TmComp compartment object.";
     throw steps::ArgErr(os.str());
     }
     */

	// Now create the actual compartments.
	ssolver::CompDefPVecCI c_end = statedef()->endComp();
    for (ssolver::CompDefPVecCI c = statedef()->bgnComp(); c != c_end; ++c)
    {
        uint compdef_gidx = (*c)->gidx();
        uint comp_idx = _addComp(*c);
        assert(compdef_gidx == comp_idx);
    }
    // Create the actual patches.
    ssolver::PatchDefPVecCI p_end = statedef()->endPatch();
    for (ssolver::PatchDefPVecCI p = statedef()->bgnPatch(); p != p_end; ++p)
    {
        uint patchdef_gidx = (*p)->gidx();
        uint patch_idx = _addPatch(*p);
        assert(patchdef_gidx == patch_idx);
    }

    // Create the diffusion boundaries
    ssolver::DiffBoundarydefPVecCI db_end = statedef()->endDiffBoundary();
    for (ssolver::DiffBoundaryDefPVecCI db = statedef()->bgnDiffBoundary(); db != db_end; ++db)
    {
    	uint diffboundary_gidx = (*db)->gidx();
    	uint diffb_idx = _addDiffBoundary(*db);
    	assert(diffboundary_gidx == diffb_idx);
    }

    /*

    got to here. Now have to figure out how to connect patch triangles to the WmComps if
    one of the compartments is well-mixed (or is that done automatically??)
     */

    uint npatches = pPatches.size();
    assert (mesh()->_countPatches() == npatches);
    for (uint p = 0; p < npatches; ++p)
    {
    	// Now add the tris for this patch
    	// We have checked the indexing- p is the global index
    	steps::wm::Patch * wmpatch = mesh()->_getPatch(p);
        // Perform upcast
        if (steps::tetmesh::TmPatch * tmpatch = dynamic_cast<steps::tetmesh::TmPatch*>(wmpatch))
        {
        	steps::tetexact::Patch * localpatch = pPatches[p];
			std::vector<uint> const & triindcs = tmpatch->_getAllTriIndices();
			std::vector<uint>::const_iterator t_end = triindcs.end();
			for (std::vector<uint>::const_iterator t = triindcs.begin();
				 t != t_end; ++t)
			{
				steps::tetmesh::Tri * tri = new steps::tetmesh::Tri(mesh(), (*t));
				assert (tri->getPatch() == tmpatch);
				double area = tri->getArea();

				double l0 = tri->getBar0Length();
				double l1 = tri->getBar1Length();
				double l2 = tri->getBar2Length();

				std::vector<int> tris = tri->getTriIdxs(tmpatch);

				int tri0 = tris[0];
				int tri1 = tris[1];
				int tri2 = tris[2];

				double d0 = tri->getTriDist(0, tri0);
				double d1 = tri->getTriDist(1, tri1);
				double d2 = tri->getTriDist(2, tri2);


				int tetinner = tri->getTet0Idx();
				int tetouter = tri->getTet1Idx();

				_addTri((*t), localpatch, area, l0, l1, l2, d0, d1, d2, tetinner, tetouter, tri0, tri1, tri2);
				delete tri;
			}
        }
        else
        {
    		std::ostringstream os;
    		os << "Well-mixed patches not supported in steps::solver::Tetexact solver";
    		throw steps::ArgErr(os.str());
        }
    }

    ncomps = pComps.size();
    assert (mesh()->_countComps() == ncomps);

    for (uint c = 0; c < ncomps; ++c)
    {
        // Now add the tets for this comp
     	// We have checked the indexing- c is the global index
        steps::wm::Comp * wmcomp = mesh()->_getComp(c);
        // Perform upcast
        if (steps::tetmesh::TmComp * tmcomp = dynamic_cast<steps::tetmesh::TmComp*>(wmcomp))
        {
         	steps::tetexact::Comp * localcomp = pComps[c];
           	std::vector<uint> const & tetindcs = tmcomp->_getAllTetIndices();
           	std::vector<uint>::const_iterator t_end = tetindcs.end();
           	for (std::vector<uint>::const_iterator t = tetindcs.begin();
   				t != t_end; ++t)
           	{
           		steps::tetmesh::Tet * tet = new steps::tetmesh::Tet(mesh(), (*t));
           		assert (tet->getComp() == tmcomp);
           		double vol = tet->getVol();
           		double a0 = tet->getTri0Area();
           		double a1 = tet->getTri1Area();
           		double a2 = tet->getTri2Area();
           		double a3 = tet->getTri3Area();
           		double d0 = tet->getTet0Dist();
           		double d1 = tet->getTet1Dist();
           		double d2 = tet->getTet2Dist();
           		double d3 = tet->getTet3Dist();
           		// At this point fetch the indices of neighbouring tets too
           		int tet0 = tet->getTet0Idx();
           		int tet1 = tet->getTet1Idx();
           		int tet2 = tet->getTet2Idx();
           		int tet3 = tet->getTet3Idx();
           		_addTet((*t), localcomp, vol, a0, a1, a2, a3, d0, d1, d2, d3,
							tet0, tet1, tet2, tet3);
           		delete tet;
           	}
        }
        else
        {
         	// This means that this compartment is a well-mixed compartment
           	// It will behave like a tetrahedral-based compartment, but
           	// contain only one 'voxel' that is connected to all surface
           	// triangles and has the same volume as the whole compartment
           	steps::tetexact::Comp * localcomp = pComps[c];
           	uint cidx = c;
           	_addWmVol(cidx, localcomp, localcomp->def()->vol());

           	// Now find all the triangles that reference this well-mixed volume
           	// and set the inner or outer tetrahedron index accordingly.

           	uint nopatches = wmcomp->_countOPatches();
           	for (uint i = 0; i < nopatches; ++i)
           	{
           		steps::wm::Patch * op = wmcomp->_getOPatch(i);
           		// 	Comp may have no outer patch
           		if (op != 0)
           		{
           			if (steps::tetmesh::TmPatch * comp_opatch = dynamic_cast<steps::tetmesh::TmPatch*>(op))
           			{
           				std::vector<uint> const & triindcs = comp_opatch->_getAllTriIndices();
           				std::vector<uint>::const_iterator t_end = triindcs.end();
           				for (std::vector<uint>::const_iterator t = triindcs.begin();
							t != t_end; ++t)
           				{
           					// Well-mixed volume should have been added above
           					assert (pWmVols[c] != 0);
           					pTris[(*t)]->setInnerTet(pWmVols[c]);
							// Add triangle to WmVols' table of neighbouring triangles.
							pWmVols[c]->setNextTri(pTris[(*t)]);
           				}
           			}
           			else
           			{
           				// This should never happen, but ok throw an exception
           				std::ostringstream os;
           				os << "Well-mixed patches not supported in steps::solver::Tetexact solver";
           				throw steps::ProgErr(os.str());
           			}
           		}
           	}

           	uint nipatches = wmcomp->_countIPatches();
           	for (uint i = 0; i < nipatches; ++i)
           	{
           		steps::wm::Patch * ip = wmcomp->_getIPatch(i);
				// Comp may not have an inner patch
				if (ip != 0)
				{
					if (steps::tetmesh::TmPatch * comp_ipatch = dynamic_cast<steps::tetmesh::TmPatch*>(ip))
					{
						std::vector<uint> const & triindcs = comp_ipatch->_getAllTriIndices();
						std::vector<uint>::const_iterator t_end = triindcs.end();
						for (std::vector<uint>::const_iterator t = triindcs.begin();
							t != t_end; ++t)
						{
							// Well-mixed volume should have been added above
							assert (pWmVols[c] != 0);
							pTris[(*t)]->setOuterTet(pWmVols[c]);
							// Add triangle to WmVols' table of neighbouring triangles.
							pWmVols[c]->setNextTri(pTris[(*t)]);
						}
					}
					else
					{
						// This should never happen, but ok throw an exception
						std::ostringstream os;
						os << "Well-mixed patches not supported in steps::solver::Tetexact solver";
						throw steps::ProgErr(os.str());
					}
				}
           	}

        }
    }


    // All tets and tris that belong to some comp or patch have been created
    // locally- now we can connect them locally
    // NOTE: currently if a tetrahedron's neighbour belongs to a different
    // comp they do not talk to each other (see stex::Tet::setNextTet())
    //

    assert (ntets == pTets.size());
    // pTets member size of all tets in geometry, but may not be filled with
    // local tets if they have not been added to a compartment
    /*
	 std::vector<steps::tetexact::Tet *>::const_iterator t_end = pTets.end();
	 for (std::vector<steps::tetexact::Tet *>::const_iterator t = pTets.begin();
	 t < nlocaltets; ++t)
	 */
    for (uint t = 0; t < ntets; ++t)
    {
    	if (pTets[t] == 0) continue;
    	int tet0 = pTets[t]->tet(0);
    	int tet1 = pTets[t]->tet(1);
    	int tet2 = pTets[t]->tet(2);
    	int tet3 = pTets[t]->tet(3);
    	// DEBUG 19/03/09 : again, incorrectly didn't allow for tet index == 0
    	if (tet0 >= 0)
    	{
			if (pTets[tet0] != 0) pTets[t]->setNextTet(0, pTets[tet0]);
    	}
    	if (tet1 >= 0)
    	{
			if (pTets[tet1] != 0) pTets[t]->setNextTet(1, pTets[tet1]);
    	}
    	if (tet2 >= 0)
    	{
    		if (pTets[tet2] != 0) pTets[t]->setNextTet(2, pTets[tet2]);
    	}
    	if (tet3 >= 0)
    	{
    		if (pTets[tet3] != 0) pTets[t]->setNextTet(3, pTets[tet3]);
    	}
    	// Not setting Tet triangles at this point- only want to set
    	// for surface triangles
    }
    assert (ntris == pTris.size());

    for (uint t = 0; t < ntris; ++t)
    {
    	// Looping over all possible tris, but only some have been added to a patch
    	if (pTris[t] == 0) continue;

    	int tri0 = pTris[t]->tri(0);
    	int tri1 = pTris[t]->tri(1);
    	int tri2 = pTris[t]->tri(2);

    	if (tri0 >= 0)
    	{
			if (pTris[tri0] != 0) pTris[t]->setNextTri(0, pTris[tri0]);
    	}
    	if (tri1 >= 0)
    	{
			if (pTris[tri1] != 0) pTris[t]->setNextTri(1, pTris[tri1]);
    	}
    	if (tri2 >= 0)
    	{
    		if (pTris[tri2] != 0) pTris[t]->setNextTri(2, pTris[tri2]);
    	}

    	// By convention, triangles in a patch should have an inner tetrahedron defined
    	// (neighbouring tets 'flipped' if necessary in Tetmesh)
    	// but not necessarily an outer tet
    	// 17/3/10- actually this is not the case any more with well-mixed compartments
    	//
    	int tetinner = pTris[t]->tet(0);
    	int tetouter = pTris[t]->tet(1);

    	/* INNER MAY NOT EXIST FOR EFIELD VERSION NOW INNER TETS MAY BE WELL_MIXED COMPARTMENTS
    	// DEBUG 18/03/09:
    	// Now correct check, previously didn't allow for tet index == 0
    	assert(tetinner >= 0);
    	assert(pTets[tetinner] != 0 );
		*/

    	// Now inside and outside tetrahedrons may be normal tetrahedrons, which
    	// means compartment is a mesh compartment, or wmvols describing a
    	// well-mixed compartment with multiple triangle connections.


    	if (tetinner >= 0)
    	{
    		// NEW FOR THIS VERSION: Tris store index of inner and outer tet (outer may not exist if on
    		// surface) but tets may not belong to a compartment, even inner tets now
    		// since they may be well-mixed compartments
    		//
    		if (pTets[tetinner] != 0)
    		{
    			// A triangle may already have an inner tet defined as a well-mixed
    			// volume, but that should not be the case here:
    			assert (pTris[t]->iTet() == 0);

    			pTris[t]->setInnerTet(pTets[tetinner]);
    			// Now add this triangle to inner tet's list of neighbours
    			for (uint i=0; i <= 4; ++i)
    			{
    				// include assert for debugging purposes and remove
    				// once this is tested
    				assert (i < 4);														//////////
    				// check if there is already a neighbouring tet or tri
    				// In theory if there is a tri to add, the tet should
    				// have less than 4 neighbouring tets added because
    				// a neighbouring tet(s) is in a different compartment

    				// THIS IS NOT THE CASE ANYMORE: tets in different compartments can be neighbours
    				// so as to allow for diffusion boundaries

    				// 	Also check tris because in some cases a surface tet
    				// may have more than 1 neighbouring tri
    				// NOTE: The order here will end up being different to the
    				// neighbour order at the Tetmesh level

    				// Now with diffusion boundaries, meaning tets can have neighbours that
    				// are in different comps, we must check the compartment
    				steps::tetexact::Tet * tet_in = pTets[tetinner];
    				if (tet_in->nextTet(i) != 0 && tet_in->compdef() == tet_in->nextTet(i)->compdef()) continue;

    				if (tet_in->nextTri(i) != 0) continue;
    				tet_in->setNextTri(i, pTris[t]);
    				break;
    			}
    		}
    	}

    	// DEBUG 18/03/09:
    	// Now correct check, previously didn't allow for tet index == 0
    	if (tetouter >= 0)
    	{
    		if (pTets[tetouter] != 0)
    		{
    			// A triangle may already have an inner tet defined as a well-mixed
    			// volume, but that should not be the case here:
    			assert (pTris[t]->oTet() == 0);

    			pTris[t]->setOuterTet(pTets[tetouter]);
    			// Add this triangle to outer tet's list of neighbours
    			for (uint i=0; i <= 4; ++i)
    			{
    				assert (i < 4);

    				// See above in that tets now store tets from different comps
    				steps::tetexact::Tet * tet_out = pTets[tetouter];

    				if (tet_out->nextTet(i) != 0 && tet_out->compdef() == tet_out->nextTet(i)->compdef()) continue;

    				if (tet_out->nextTri(i) != 0) continue;
    				tet_out->setNextTri(i, pTris[t]);
    				break;
				}
    		}
    	}
    }

    // Now loop over the diffusion boundaries:
    // 1) get all the triangles and get the two tetrahedrons
    // 2) figure out which direction is the direction for a tetrahedron
    // 3) add the tetrahedron and the direction to local object

    // This is here because we need all tets to have been assigned correctly
    // to compartments. Check every one and set the compA and compB for the db
    uint ndiffbnds = pDiffBoundaries.size();
    assert(ndiffbnds ==	mesh()->_countDiffBoundaries());

    for (uint db = 0; db < ndiffbnds; ++db)
    {
    	steps::tetexact::DiffBoundary * localdiffb = pDiffBoundaries[db];
        std::vector<uint> dbtrisvec = localdiffb->def()->tris();

        uint compAidx = localdiffb->def()->compa();
        uint compBidx = localdiffb->def()->compb();
        steps::solver::Compdef * compAdef = statedef()->compdef(compAidx);
        steps::solver::Compdef * compBdef = statedef()->compdef(compBidx);

        std::vector<uint>::const_iterator dbtris_end = dbtrisvec.end();
        for (std::vector<uint>::const_iterator dbtris = dbtrisvec.begin(); dbtris != dbtris_end; ++dbtris)
        {
    		steps::tetmesh::Tri * tri = new steps::tetmesh::Tri(mesh(), (*dbtris));

    		uint tetAidx = tri->getTet0Idx();
    		uint tetBidx = tri->getTet1Idx();
    		assert (tetAidx >= 0 and tetBidx >= 0);

    		delete tri;

    		steps::tetexact::Tet * tetA = _tet(tetAidx);
    		steps::tetexact::Tet * tetB = _tet(tetBidx);
    		assert(tetA != 0 and tetB != 0);

    		steps::solver::Compdef * tetA_cdef = tetA->compdef();
    		steps::solver::Compdef * tetB_cdef = tetB->compdef();
    		assert (tetA_cdef != 0);
    		assert (tetB_cdef != 0);

    		if (tetA_cdef != compAdef)
    		{
    			assert (tetB_cdef == compAdef);
    			assert (tetA_cdef == compBdef);
    		}
    		else
    		{
    			assert (tetB_cdef == compBdef);
    			assert (tetA_cdef == compAdef);
    		}

    		// Ok, checks over, lets get down to business
    		int direction_idx_a = -1;
    		int direction_idx_b = -1;

    		steps::tetmesh::Tet * tetA_mesh = new steps::tetmesh::Tet(mesh(), (tetAidx));
    		steps::tetmesh::Tet * tetB_mesh = new steps::tetmesh::Tet(mesh(), (tetBidx));

    		for (uint i = 0; i < 4; ++i)
    		{
    			if (tetA_mesh->getTriIdx(i) == (*dbtris))
    			{
    				assert(direction_idx_a == -1);
    				direction_idx_a = i;
    			}
    			if (tetB_mesh->getTriIdx(i) == (*dbtris))
    			{
    				assert(direction_idx_b == -1);
    				direction_idx_b = i;
    			}
    		}
    		assert (direction_idx_a != -1);
    		assert (direction_idx_b != -1);

    		// Set the tetrahedron and direction to the Diff Boundary object
    		localdiffb->setTetDirection(tetAidx, direction_idx_a);
    		localdiffb->setTetDirection(tetBidx, direction_idx_b);

    		delete tetA_mesh;
    		delete tetB_mesh;
        }

        localdiffb->setComps(_comp(compAidx), _comp(compBidx));

        // Before the kprocs are set up ( in _setup) the tetrahedrons need to know the diffusion
        // boundary direction, so let's do it here  - the diff bounday has had all
        // tetrahedrons added

        // Might as well copy the vectors because we need to index through
        std::vector<uint> tets = localdiffb->getTets();
        std::vector<uint> tets_direction = localdiffb->getTetDirection();

        ntets = tets.size();
        assert (ntets <= pTets.size());
        assert (tets_direction.size() == ntets);

        for (uint t = 0; t < ntets; ++t)
        {
        	steps::tetexact::Tet * localtet = _tet(tets[t]);
        	localtet->setDiffBndDirection(tets_direction[t]);
        }

    }


	TetPVecCI tet_end = pTets.end();
	for (TetPVecCI t = pTets.begin(); t != tet_end; ++t)
	{
		// DEBUG: vector holds all possible tetrahedrons,
		// but they have not necessarily been added to a compartment.
		if ((*t) == 0) continue;

		(*t)->setupKProcs(this);
	}

	WmVolPVecCI wmv_end = pWmVols.end();
	for (WmVolPVecCI wmv = pWmVols.begin(); wmv != wmv_end; ++wmv)
	{
		// Vector allows for all compartments to be well-mixed, so
		// hold null-pointer for mesh compartments
		if ((*wmv) == 0) continue;

		(*wmv)->setupKProcs(this);
	}

	TriPVecCI tri_end = pTris.end();
	for (TriPVecCI t = pTris.begin(); t != tri_end; ++t)
	{
		// DEBUG: vector holds all possible triangles, but
		// only patch triangles are filled
		if ((*t) == 0) continue;

		(*t)->setupKProcs(this, efflag());
	}
	// Resolve all dependencies
	for (TetPVecCI t = pTets.begin(); t != tet_end; ++t)
	{
		// DEBUG: vector holds all possible tetrahedrons,
		// but they have not necessarily been added to a compartment.
		if ((*t) == 0) continue;

		KProcPVecCI kprocend = (*t)->kprocEnd();
		for (KProcPVecCI k = (*t)->kprocBegin(); k != kprocend; ++k)
		{
		    (*k)->setupDeps();
		}
	}


	for (WmVolPVecCI wmv = pWmVols.begin(); wmv != wmv_end; ++wmv)
	{
		// Vector allows for all compartments to be well-mixed, so
		// hold null-pointer for mesh compartments
		if ((*wmv) == 0) continue;
		KProcPVecCI kprocend = (*wmv)->kprocEnd();
		for (KProcPVecCI k = (*wmv)->kprocBegin(); k != kprocend; ++k)
		{
		    (*k)->setupDeps();
		}
	}

	for (TriPVecCI t = pTris.begin(); t != tri_end; ++t)
	{
		// DEBUG: vector holds all possible triangles, but
		// only patch triangles are filled
		if ((*t) == 0) continue;

	    KProcPVecCI kprocend = (*t)->kprocEnd();
	    for (KProcPVecCI k = (*t)->kprocBegin(); k != kprocend; ++k)
	    {
	        (*k)->setupDeps();
	    }
	}

	// Create EField structures if EField is to be calculated
	if (efflag() == true) _setupEField();

	nEntries = pKProcs.size();
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setupEField(void)
{

	//// Note to self: for now roughly following flow from original code in sim/controller.py.
	//// code for setting up a mesh was in func_tetmesh constructor and called functions.

	assert(efflag() == true);

	// Give temperature a default value of 20c
	pTemp = 293.15;

	uint nmembs = mesh()->_countMembs();

	if  (nmembs != 1)
	{
		std::ostringstream os;
		os << "Membrane potential solver currently supports only one ";
		os << "membrane description object.";
		throw steps::ArgErr(os.str());
	}

	steps::tetmesh::Memb * memb = mesh()->_getMemb(0);
	assert(memb != 0);

	// TODO: Decide what checks are needed for the membrane and implement them here

	pEFNTets = memb->countVolTets();
	pEFNTris = memb->countTris();
	pEFNVerts = memb->countVerts();

    pEFTets = new uint[neftets() * 4];
    assert(pEFTets != 0);

    // All the triangles we will count here are the real membrane triangles,
    // virtual triangles will not require a capacitance.
    pEFTris = new uint[neftris() * 3];
    assert(pEFTris != 0);

    pEFVerts = new double[nefverts() * 3];
    assert(pEFVerts != 0);

    uint nverts = mesh()->countVertices();
    uint ntris = mesh()->countTris();
    uint ntets= mesh()->countTets();

    pEFVert_GtoL = new int[nverts];
    for (uint i=0; i < nverts; ++i) pEFVert_GtoL[i] = -1;
    pEFTri_GtoL = new int[ntris];
    for (uint i=0; i< ntris; ++i) pEFTri_GtoL[i] = -1;
    pEFTet_GtoL = new int[ntets];
    for (uint i=0; i < ntets; ++i) pEFTet_GtoL[i] = -1;

    pEFTri_LtoG = new uint[neftris()];

    // Copy the data to local structures.

    std::vector<uint> membverts = memb->_getAllVertIndices();
    assert(membverts.size() == nefverts());
    for (uint efv = 0; efv < nefverts(); ++efv)
    {
    	uint vertidx = membverts[efv];
    	double* verttemp = mesh()->_getVertex(vertidx);
    	uint efv2 = efv*3;

    	// CONVERTING TO MICRONS HERE. EFIELD OBJECT WILL NOT PERFORM THIS CONVERSION
    	pEFVerts[efv2] = verttemp[0] * 1.0e6;
    	pEFVerts[efv2+1] = verttemp[1] * 1.0e6;
    	pEFVerts[efv2+2] = verttemp[2] * 1.0e6;

    	pEFVert_GtoL[vertidx] = efv;
    }

    std::vector<uint> membtets = memb->_getAllVolTetIndices();
    assert(membtets.size() == neftets());
    for (uint eft=0; eft < neftets(); ++eft)
    {
    	uint tetidx = membtets[eft];
    	uint* tettemp = mesh()->_getTet(tetidx);
    	uint eft2 = eft*4;

    	// Convert to indices used by EField object
    	int tv0 =  pEFVert_GtoL[tettemp[0]];
    	int tv1 = pEFVert_GtoL[tettemp[1]];
    	int tv2 = pEFVert_GtoL[tettemp[2]];
    	int tv3 = pEFVert_GtoL[tettemp[3]];
    	if  (tv0 ==-1 || tv1 == -1 || tv2 == -1 || tv3 == -1)
    	{
    		std::ostringstream os;
    		os << "Failed to create EField structures.";
    		throw steps::ProgErr(os.str());
    	}

    	pEFTets[eft2] = tv0;
    	pEFTets[eft2+1] = tv1;
    	pEFTets[eft2+2] = tv2;
    	pEFTets[eft2+3] = tv3;

    	pEFTet_GtoL[tetidx] = eft;
    }

    std::vector<uint> membtris = memb->_getAllTriIndices();
    assert(membtris.size() == neftris());

    pEFTris_vec.resize(neftris());

    for (uint eft = 0; eft < neftris(); ++eft)
    {
    	uint triidx = membtris[eft];
    	uint* tritemp = mesh()->_getTri(triidx);
    	uint eft2 = eft*3;

    	// Convert to indices used by EField object
    	int tv0 =  pEFVert_GtoL[tritemp[0]];
    	int tv1 = pEFVert_GtoL[tritemp[1]];
    	int tv2 = pEFVert_GtoL[tritemp[2]];
    	if  (tv0 ==-1 || tv1 == -1 || tv2 == -1)
    	{
    		std::ostringstream os;
    		os << "Failed to create EField structures.";
    		throw steps::ProgErr(os.str());
    	}

    	pEFTris[eft2] = tv0;
    	pEFTris[eft2+1] = tv1;
    	pEFTris[eft2+2] = tv2;

    	pEFTri_GtoL[triidx] = eft;
    	pEFTri_LtoG[eft] = triidx;

    	// This is added now for quicker iteration during run()
    	// Extremely important for larger meshes, orders of magnitude times faster
    	pEFTris_vec[eft] = pTris[triidx];
    }

    pEField = new steps::solver::efield::EField(nefverts(), pEFVerts, neftris(), pEFTris, neftets(), pEFTets, memb->_getOpt_method(), memb->_getOpt_file_name());
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::saveMembOpt(std::string const & opt_file_name)
{
	if  (efflag() != true)
	{
		std::ostringstream os;
		os << "saveMembOpt method only available if running EField ";
		throw steps::ArgErr(os.str());
	}

	pEField->saveOptimal(opt_file_name);

}

////////////////////////////////////////////////////////////////////////////////

uint stex::Tetexact::_addComp(steps::solver::Compdef * cdef)
{
    stex::Comp * comp = new Comp(cdef);
    assert(comp != 0);
    uint compidx = pComps.size();
    pComps.push_back(comp);
    pCompMap[cdef] = comp;
    return compidx;
}

////////////////////////////////////////////////////////////////////////////////

uint stex::Tetexact::_addPatch(steps::solver::Patchdef * pdef)
{
    /* Comp * icomp = 0;
	 Comp * ocomp = 0;
	 if (pdef->icompdef()) icomp = pCompMap[pdef->icompdef()];
	 if (pdef->ocompdef()) ocomp = pCompMap[pdef->ocompdef()];
	 */
    stex::Patch * patch = new Patch(pdef);
    assert(patch != 0);
    uint patchidx = pPatches.size();
    pPatches.push_back(patch);
    return patchidx;
}

////////////////////////////////////////////////////////////////////////////////

uint stex::Tetexact::_addDiffBoundary(steps::solver::DiffBoundarydef * dbdef)
{
	stex::DiffBoundary * diffb = new DiffBoundary(dbdef);
	assert(diffb != 0);
	uint dbidx = pDiffBoundaries.size();
	pDiffBoundaries.push_back(diffb);
	return dbidx;
}


////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_addTet(uint tetidx,
							 steps::tetexact::Comp * comp, double vol,
							 double a1, double a2, double a3, double a4,
							 double d1, double d2, double d3, double d4,
							 int tet0, int tet1, int tet2, int tet3)
{
	steps::solver::Compdef * compdef  = comp->def();
    stex::Tet * localtet = new stex::Tet(tetidx, compdef, vol, a1, a2, a3, a4, d1, d2, d3, d4,
									     tet0, tet1, tet2, tet3);
    assert(localtet != 0);
    assert(tetidx < pTets.size());
    assert(pTets[tetidx] == 0);
    pTets[tetidx] = localtet;
    comp->addTet(localtet);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_addWmVol(uint cidx, steps::tetexact::Comp * comp, double vol)
{
	steps::solver::Compdef * compdef  = comp->def();
	stex::WmVol * localtet = new stex::WmVol(cidx, compdef, vol);
	assert(localtet != 0);
	assert(cidx < pWmVols.size());
	pWmVols[cidx] = localtet;
	comp->addTet(localtet);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_addTri(uint triidx, steps::tetexact::Patch * patch, double area,
							 double l0, double l1, double l2, double d0, double d1, double d2,  int tinner, int touter, int tri0, int tri1, int tri2)
{
    steps::solver::Patchdef * patchdef = patch->def();
    stex::Tri * tri = new stex::Tri(triidx, patchdef, area, l0, l1, l2, d0, d1, d2,  tinner, touter, tri0, tri1, tri2);
    assert(tri != 0);
    assert (triidx < pTris.size());
    assert (pTris[triidx] == 0);
    pTris[triidx] = tri;
    patch->addTri(tri);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::reset(void)
{
	std::for_each(pComps.begin(), pComps.end(), std::mem_fun(&Comp::reset));
	std::for_each(pPatches.begin(), pPatches.end(), std::mem_fun(&Patch::reset));

    TetPVecCI tet_end = pTets.end();
    for (TetPVecCI tet = pTets.begin(); tet != tet_end; ++tet)
    {
    	if (*tet == 0) continue;
    	(*tet)->reset();
    }

    WmVolPVecCI wmvol_end = pWmVols.end();
    for (WmVolPVecCI wmvol = pWmVols.begin(); wmvol != wmvol_end; ++wmvol)
    {
    	if ((*wmvol) == 0) continue;
    	(*wmvol)->reset();
    }

	TriPVecCI tri_end = pTris.end();
	for (TriPVecCI t = pTris.begin(); t != tri_end; ++t)
	{
		if ((*t) == 0) continue;
		(*t)->reset();
	}

    uint ngroups = nGroups.size();
    for (uint i = 0; i < ngroups; i++) {
        free(nGroups[i]->indices);
        delete nGroups[i];
    }
    nGroups.clear();

    ngroups = pGroups.size();
    for (uint i = 0; i < ngroups; i++) {
        free(pGroups[i]->indices);
        delete pGroups[i];
    }
    pGroups.clear();

    pSum = 0.0;
    nSum = 0.0;
    pA0 = 0.0;

	statedef()->resetTime();
	statedef()->resetNSteps();
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::run(double endtime)
{
	if (efflag() == false)
	{
		if (endtime < statedef()->time())
		{
			std::ostringstream os;
			os << "Endtime is before current simulation time";
			throw steps::ArgErr(os.str());
		}
		while (statedef()->time() < endtime)
		{
			stex::KProc * kp = _getNext();
			if (kp == 0) break;
			double a0 = getA0();
			if (a0 == 0.0) break;
			double dt = rng()->getExp(a0);
			if ((statedef()->time() + dt) > endtime) break;
			_executeStep(kp, dt);
		}
		statedef()->setTime(endtime);
	}
	else if (efflag() == true)
	{
		// Run the simulation, including the EField calculation.
		// This loop will assume that the SSA dt is sufficiently small so
		// that a number of SSA events execute between every EField calculation.
		// A warning message will be printed if the SSA dt is larger than the EField dt.
		// The EField dt is actually a MAXIMUM dt- the actual time for the EField
		// calculation will be exact with respect to the last event time in the
		// SSA before reaching the EField dt.
		while (statedef()->time() < endtime)
		{
			// The zero propensity
			double a0 = getA0();
			// We need a bool to check if the SSA contains no possible events. In
			// this rare case (continue to) execute the EField calculation to the endtime.
			bool ssa_on = true;
			double ssa_dt = 0.0;
			if (a0 != 0.0) ssa_dt = rng()->getExp(a0);
			else (ssa_on = false);
			// Set the actual efield dt. This value will take a maximum pEFDT.
			double ef_dt = 0.0;

			while (ssa_on && (ef_dt + ssa_dt) < pEFDT )
			{
				stex::KProc * kp = _getNext();
				if (kp == 0) break;
				_executeStep(kp, ssa_dt);
				ef_dt += ssa_dt;

				a0 = getA0();
				if (a0 != 0.0) ssa_dt = rng()->getExp(a0);
				else (ssa_on = false);

			}
			assert(ef_dt < pEFDT);

			// It's possible that ef_dt is zero here: ssa_dt is large, or has become large.
			// In that case print a warning but continue, running the EField simulation for EFDT
			if (ssa_on == false || ef_dt == 0.0)
			{
	    		std::ostringstream os;
	    		//os << "\nWARNING: SSA tau is larger than EField dt.";
	    		//std::cout << os << std::endl;
			}

			if (ef_dt == 0.0)
			{
				// This means that tau is larger than EField dt. We have no choice but to
				// increase the state time by pEFDT.
				ef_dt = pEFDT;
				statedef()->incTime(pEFDT);
			}

			// Now to perform the EField calculation. This means finding ohmic and GHK
			// currents from triangles during the ef_dt and applying these to the EField
			// object.

			/*
			 * Older method is much, much slower for larger meshes. e.g. for ~200,000
			 * tets this method is 2 or 3 orders of magnitude slower than pEFTris_vec method
			 *
			for (uint ltidx = 0; ltidx < neftris(); ++ltidx)
			{
				uint gtidx = pEFTri_LtoG[ltidx];

				Tri * tritemp = tris()[gtidx];
				// We can use local indices for faster EField method
				double v = pEField->getTriV(ltidx);
				// The triangle will compute all GHK and Ohmic currents
				double cur = tritemp->computeI(v, ef_dt, statedef()->time());
				pEField->setTriI(ltidx, cur);
			}
			*/

			TriPVecCI eftri_end = pEFTris_vec.end();
			uint tlidx = 0;
			double sttime = statedef()->time();
			for (TriPVecCI eft = pEFTris_vec.begin(); eft != eftri_end; ++eft)
			{
				double v = pEField->getTriV(tlidx);
				double cur = (*eft)->computeI(v, ef_dt, sttime);
				pEField->setTriI(tlidx, cur);
				tlidx++;
			}

			pEField->advance(ef_dt);
			// TODO: Replace this with something that only resets voltage-dependent things
			_update();
		}
	}

	else assert(false);

}

////////////////////////////////////////////////////////////////////////

void stex::Tetexact::advance(double adv)
{
	if (adv < 0.0)
	{
		std::ostringstream os;
		os << "Time to advance cannot be negative";
	    throw steps::ArgErr(os.str());
	}

	double endtime = statedef()->time() + adv;
	run(endtime);
}

////////////////////////////////////////////////////////////////////////
/*
void stex::Tetexact::advanceSteps(uint nsteps)
{
	while (nsteps != 0) {
        stex::KProc * kp = _getNext();
        if (kp == 0) return;
        double a0 = getA0();
        if (a0 == 0.0) return;
        double dt = rng()->getExp(a0);
        _executeStep(kp, dt);
        nsteps--;
    }
}
*/
////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::step(void)
{
	if (efflag() == true)
	{
		std::ostringstream os;
		os << "Method not available with EField calculation.";
		throw steps::ArgErr(os.str());
	}

	stex::KProc * kp = _getNext();
	if (kp == 0) return;
	double a0 = getA0();
	if (a0 == 0.0) return;
	double dt = rng()->getExp(a0);
	_executeStep(kp, dt);
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::getTime(void) const
{
	return statedef()->time();
}

////////////////////////////////////////////////////////////////////////

uint stex::Tetexact::getNSteps(void) const
{
    return statedef()->nsteps();
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::setTime(double time)
{
	statedef()->setTime(time);
}

////////////////////////////////////////////////////////////////////////

void stex::Tetexact::setNSteps(uint nsteps)
{
    statedef()->setNSteps(nsteps);
}

////////////////////////////////////////////////////////////////////////

void stex::Tetexact::setTemp(double t)
{
	if (efflag() == false)
	{
		std::ostringstream os;
		os << "\nWARNING: Temperature set in simulation without membrane ";
		os << "potential calculation will be ignored.\n";
		std::cout << os << std::endl;
	}
	assert(t >= 0.0);
	pTemp = t;
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getCompVol(uint cidx) const
{
	assert(cidx < statedef()->countComps());
	assert (statedef()->countComps() == pComps.size());
	stex::Comp * comp = _comp(cidx);
	assert(comp != 0);
	return comp->vol();
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getCompCount(uint cidx, uint sidx) const
{
	assert(cidx < statedef()->countComps());
	assert(sidx < statedef()->countSpecs());
	assert (statedef()->countComps() == pComps.size());
	stex::Comp * comp = _comp(cidx);
	assert(comp != 0);
	uint slidx = comp->def()->specG2L(sidx);
	if (slidx == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "Species undefined in compartment.\n";
		throw steps::ArgErr(os.str());
	}

	uint count = 0;
	WmVolPVecCI t_end = comp->endTet();
	for (WmVolPVecCI t = comp->bgnTet(); t != t_end; ++t)
	{
		count += (*t)->pools()[slidx];
	}

	return count;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setCompCount(uint cidx, uint sidx, double n)
{
	assert(cidx < statedef()->countComps());
	assert(sidx < statedef()->countSpecs());
	assert (n >= 0.0);
	assert (statedef()->countComps() == pComps.size());
	stex::Comp * comp = _comp(cidx);
	uint slidx = comp->def()->specG2L(sidx);
	if (slidx == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "Species undefined in compartment.\n";
		throw steps::ArgErr(os.str());
	}
	if (n > std::numeric_limits<unsigned int>::max( ))
	{
		std::ostringstream os;
		os << "Can't set count greater than maximum unsigned integer (";
		os << std::numeric_limits<unsigned int>::max( ) << ").\n";
		throw steps::ArgErr(os.str());
	}

	double totalvol = comp->def()->vol();
	WmVolPVecCI t_end = comp->endTet();

	double n_int = std::floor(n);
	double n_frc = n - n_int;
	uint c = static_cast<uint>(n_int);
	if (n_frc > 0.0)
	{
		double rand01 = rng()->getUnfIE();
		if (rand01 < n_frc) c++;
	}

    uint nremoved = 0;
    for (WmVolPVecCI t = comp->bgnTet(); t != t_end; ++t)
    {
        WmVol * tet = *t;
        if (n == 0.0) {
            tet->setCount(slidx, 0);
            continue;
        }

		// New method (allowing ceiling) means we have to set the counts
		// to zero for any tets after all molecules have been injected
		if (nremoved == c)
		{
			tet->setCount(slidx, 0);
			continue;
		}

        double fract = static_cast<double>(c) * (tet->vol() / totalvol);
        uint n3 = static_cast<uint>(std::floor(fract));

		// BUGFIX 29/09/2010 IH. By not allowing the ceiling here
		// concentration gradients could appear in the injection.
        double n3_frac = fract - static_cast<double>(n3);
        if (n3_frac > 0.0)
        {
            double rand01 = rng()->getUnfIE();
            if (rand01 < n3_frac) n3++;
        }

        ///
        // BUGFIX 18/11/09 IH. By reducing c here we were not giving all
        // tets an equal share. Tets with low index would have a
        // greater share than those with higher index.
        //c -= n3;
        nremoved += n3;

        if (nremoved >= c)
        {
            n3 -= (nremoved-c);
            nremoved = c;
        }

        tet->setCount(slidx, n3);
    }
    assert(nremoved <= c);
    c -= nremoved;
	while (c != 0)
	{
        WmVol * tet = comp->pickTetByVol(rng()->getUnfIE());

        assert (tet != 0);
        tet->setCount(slidx, (tet->pools()[slidx] + 1.0));
        c--;
	}
	for (WmVolPVecCI t = comp->bgnTet(); t != t_end; ++t)
	{
		_updateSpec(*t, slidx);
	}
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getCompAmount(uint cidx, uint sidx) const
{
	// the following method does all the necessary argument checking
	double count = _getCompCount(cidx, sidx);
	return (count / smath::AVOGADRO);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setCompAmount(uint cidx, uint sidx, double a)
{
	// convert amount in mols to number of molecules
	double a2 = a * steps::math::AVOGADRO;
	// the following method does all the necessary argument checking
	_setCompCount(cidx, sidx, a2);
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getCompConc(uint cidx, uint sidx) const
{
	// the following method does all the necessary argument checking
	double count = _getCompCount(cidx, sidx);
	ssolver::Compdef * comp = statedef()->compdef(cidx);
	assert(comp != 0);
	double vol = comp->vol();
	return count/ (1.0e3 * vol * steps::math::AVOGADRO);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setCompConc(uint cidx, uint sidx, double c)
{
	assert(c >= 0.0);
	assert (cidx < statedef()->countComps());
	ssolver::Compdef * comp = statedef()->compdef(cidx);
	assert(comp != 0);
	double count = c * (1.0e3 * comp->vol() * steps::math::AVOGADRO);
	// the following method does all the necessary argument checking
	_setCompCount(cidx, sidx, count);
}

////////////////////////////////////////////////////////////////////////////////

bool stex::Tetexact::_getCompClamped(uint cidx, uint sidx) const
{
	assert(cidx < statedef()->countComps());
	assert(sidx < statedef()->countSpecs());
	assert(statedef()->countComps() == pComps.size());
	stex::Comp * comp = _comp(cidx);
	assert(comp != 0);
	uint lsidx = comp->def()->specG2L(sidx);
    if (lsidx == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "Species undefined in compartment.\n";
		throw steps::ArgErr(os.str());
	}

    WmVolPVecCI t_end = comp->endTet();
    for (WmVolPVecCI t = comp->bgnTet(); t != t_end; ++t)
    {
    	if ((*t)->clamped(lsidx) == false) return false;
    }
    return true;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setCompClamped(uint cidx, uint sidx, bool b)
{
	assert(cidx < statedef()->countComps());
	assert(sidx < statedef()->countSpecs());
	assert(statedef()->countComps() == pComps.size());
	stex::Comp * comp = _comp(cidx);
	assert(comp != 0);
	uint lsidx = comp->def()->specG2L(sidx);
	if (lsidx == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "Species undefined in compartment.\n";
		throw steps::ArgErr(os.str());
	}

	// Set the flag in def object, though this may not be necessary
	comp->def()->setClamped(lsidx, b);

	WmVolPVecCI t_end = comp->endTet();
	for (WmVolPVecCI t = comp->bgnTet(); t != t_end; ++t)
	{
		(*t)->setClamped(lsidx, b);
	}
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getCompReacK(uint cidx, uint ridx) const
{
	assert(cidx < statedef()->countComps());
	assert(ridx < statedef()->countReacs());
	assert(statedef()->countComps() == pComps.size());
	stex::Comp * comp = _comp(cidx);
	assert(comp != 0);
	uint lridx = comp->def()->reacG2L(ridx);
	if (lridx == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "Reaction undefined in compartment.\n";
		throw steps::ArgErr(os.str());
	}

	// We're just returning the default value for this comp, individual
	// tets may have different Kcsts set individually
	return (comp->def()->kcst(lridx));
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setCompReacK(uint cidx, uint ridx, double kf)
{
	assert(cidx < statedef()->countComps());
	assert(ridx < statedef()->countReacs());
	assert(statedef()->countComps() == pComps.size());
	assert(kf >= 0.0);
	stex::Comp * comp = _comp(cidx);
	assert(comp != 0);
	uint lridx = comp->def()->reacG2L(ridx);
	if (lridx == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "Reaction undefined in compartment.\n";
		throw steps::ArgErr(os.str());
	}

	// First set the default value for the comp
	comp->def()->setKcst(lridx, kf);

	// Now update all tetrahedra in this comp
	WmVolPVecCI t_end = comp->endTet();
	for (WmVolPVecCI t = comp->bgnTet(); t != t_end; ++t)
	{
		(*t)->reac(lridx)->setKcst(kf);
	}

	// Rates have changed
	_update();
}

////////////////////////////////////////////////////////////////////////////////

bool stex::Tetexact::_getCompReacActive(uint cidx, uint ridx) const
{
	assert(cidx < statedef()->countComps());
	assert(ridx < statedef()->countReacs());
	assert(statedef()->countComps() == pComps.size());
	stex::Comp * comp = _comp(cidx);
	assert(comp != 0);
	uint lridx = comp->def()->reacG2L(ridx);
	if (lridx == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "Reaction undefined in compartment.\n";
		throw steps::ArgErr(os.str());
	}

	WmVolPVecCI t_end = comp->endTet();
	for (WmVolPVecCI t = comp->bgnTet(); t != t_end; ++t)
	{
		if ((*t)->reac(lridx)->inactive() == true) return false;
	}
	return true;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setCompReacActive(uint cidx, uint ridx, bool a)
{
	assert(cidx < statedef()->countComps());
	assert(ridx < statedef()->countReacs());
	assert(statedef()->countComps() == pComps.size());
	stex::Comp * comp = _comp(cidx);
	assert(comp != 0);
	uint lridx = comp->def()->reacG2L(ridx);
	if (lridx == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "Reaction undefined in compartment.\n";
		throw steps::ArgErr(os.str());
	}

	// Set the default value for the comp, though this is not entirely
	// necessary
	comp->def()->setActive(lridx, a);

	WmVolPVecCI t_end = comp->endTet();
	for (WmVolPVecCI t = comp->bgnTet(); t != t_end; ++t)
	{
		(*t)->reac(lridx)->setActive(a);
	}
    // It's cheaper to just recompute everything.
    _update();
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getCompDiffD(uint cidx, uint didx) const
{
	assert (cidx < statedef()->countComps());
	assert (didx < statedef()->countDiffs());
	assert (statedef()->countComps() == pComps.size());
	stex::Comp * comp = _comp(cidx);
	assert (comp != 0);
	uint ldidx = comp->def()->diffG2L(didx);
	if (ldidx == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "Diffusion rule undefined in compartment.\n";
		throw steps::ArgErr(os.str());
	}

	// We're just returning the default value for this comp, individual
	// tets may have different Dcsts set individually
	return (comp->def()->dcst(ldidx));

}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setCompDiffD(uint cidx, uint didx, double dk)
{
	assert (cidx < statedef()->countComps());
	assert (didx < statedef()->countDiffs());
	assert (statedef()->countComps() == pComps.size());
	assert(dk >= 0.0);
	stex::Comp * comp = _comp(cidx);
	assert (comp != 0);
	uint ldidx = comp->def()->diffG2L(didx);
	if (ldidx == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "Diffusion rule undefined in compartment.\n";
		throw steps::ArgErr(os.str());
	}

	// First set the default value for the comp
	comp->def()->setDcst(ldidx, dk);

	// Now update all tets in this comp
	WmVolPVecCI t_end = comp->endTet();
	for (WmVolPVecCI t = comp->bgnTet(); t != t_end; ++t)
	{
		if (stex::Tet * tet = dynamic_cast<stex::Tet *>(*t))
		{
			tet->diff(ldidx)->setDcst(dk);
		}
		else
		{
			std::ostringstream os;
			os << "Cannot change diffusion constant in well-mixed compartment.";
			throw steps::ArgErr(os.str());
		}
	}

	// Rates have changed
	_update();
}

////////////////////////////////////////////////////////////////////////////////

bool stex::Tetexact::_getCompDiffActive(uint cidx, uint didx) const
{
	assert (cidx < statedef()->countComps());
	assert (didx < statedef()->countDiffs());
	assert (statedef()->countComps() == pComps.size());
	stex::Comp * comp = _comp(cidx);
	assert (comp != 0);
	uint ldidx = comp->def()->diffG2L(didx);
	if (ldidx == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "Diffusion rule undefined in compartment.\n";
		throw steps::ArgErr(os.str());
	}

	WmVolPVecCI t_end = comp->endTet();
	for (WmVolPVecCI t = comp->bgnTet(); t != t_end; ++t)
	{
		//stex::WmVol * wmcomp = (*t);
		if (stex::Tet * tet = dynamic_cast<stex::Tet *>((*t)))
		{
			if (tet->diff(ldidx)->inactive() == true) return false;
		}
		else
		{
			std::ostringstream os;
			os << "Diffusion activation not defined in well-mixed compartment.\n";
			throw steps::ArgErr(os.str());
		}
	}
	return true;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setCompDiffActive(uint cidx, uint didx, bool act)
{
	assert (cidx < statedef()->countComps());
	assert (didx < statedef()->countDiffs());
	assert (statedef()->countComps() == pComps.size());
	stex::Comp * comp = _comp(cidx);
	assert (comp != 0);
	uint ldidx = comp->def()->diffG2L(didx);
	if (ldidx == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "Diffusion rule undefined in compartment.\n";
		throw steps::ArgErr(os.str());
	}

	WmVolPVecCI t_end = comp->endTet();
	for (WmVolPVecCI t = comp->bgnTet(); t != t_end; ++t)
	{
		//stex::WmVol * wmcomp = (*t);
		if (stex::Tet * tet = dynamic_cast<stex::Tet *>((*t)))
		{
			tet->diff(ldidx)->setActive(act);
		}
		else
		{
			std::ostringstream os;
			os << "Cannot change diffusion constant in well-mixed compartment.\n";
			throw steps::ArgErr(os.str());
		}
	}
    // It's cheaper to just recompute everything.
    _update();
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getPatchArea(uint pidx) const
{
	assert (pidx < statedef()->countPatches());
	assert (statedef()->countPatches() == pPatches.size());
	stex::Patch * patch = _patch(pidx);
	assert (patch != 0);
	return patch->area();
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getPatchCount(uint pidx, uint sidx) const
{
	assert (pidx < statedef()->countPatches());
	assert (sidx < statedef()->countSpecs());
	assert (statedef()->countPatches() == pPatches.size());
	stex::Patch * patch = _patch(pidx);
	assert (patch != 0);
	uint slidx = patch->def()->specG2L(sidx);
	if (slidx == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "Species undefined in patch.\n";
		throw steps::ArgErr(os.str());
	}

	uint count = 0;
	TriPVecCI t_end = patch->endTri();
	for (TriPVecCI t = patch->bgnTri(); t != t_end; ++t)
	{
		count += (*t)->pools()[slidx];
	}

	return count;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setPatchCount(uint pidx, uint sidx, double n)
{
	assert (pidx < statedef()->countPatches());
	assert (sidx < statedef()->countSpecs());
	assert (statedef()->countPatches() == pPatches.size());
	assert (n >= 0.0);
	stex::Patch * patch = _patch(pidx);
	assert (patch != 0);
	uint slidx = patch->def()->specG2L(sidx);
	if (slidx == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "Species undefined in patch.\n";
		throw steps::ArgErr(os.str());
	}
	if (n > std::numeric_limits<unsigned int>::max( ))
	{
		std::ostringstream os;
		os << "Can't set count greater than maximum unsigned integer (";
		os << std::numeric_limits<unsigned int>::max( ) << ").\n";
		throw steps::ArgErr(os.str());
	}

	double totalarea = patch->def()->area();
	TriPVecCI t_end = patch->endTri();

	double n_int = std::floor(n);
	double n_frc = n - n_int;
	uint c = static_cast<uint>(n_int);
	if (n_frc > 0.0)
	{
		double rand01 = rng()->getUnfIE();
		if (rand01 < n_frc) c++;
	}

	uint nremoved = 0;
	for (TriPVecCI t = patch->bgnTri(); t != t_end; ++t)
	{
		Tri * tri = *t;

        if (n == 0.0) {
            tri->setCount(slidx, 0);
            continue;
        }
		// New method (allowing ceiling) means we have to set the counts
		// to zero for any triangles after all molecules have been injected
		if (nremoved == c)
		{
			tri->setCount(slidx, 0);
			continue;
		}

		double fract = static_cast<double>(c) * (tri->area()/totalarea);
        uint n3 = static_cast<uint>(std::floor(fract));

		// BUGFIX 29/09/2010 IH. By not allowing the ceiling here
		// concentration gradients could appear in the injection.
        double n3_frac = fract - static_cast<double>(n3);
        if (n3_frac > 0.0)
        {
            double rand01 = rng()->getUnfIE();
            if (rand01 < n3_frac) n3++;
        }

        // BUGFIX 18/11/09 IH. By reducing c here we were not giving all
        // triangles an equal share. Triangles with low index would have a
        // greater share than those with higher index.
        // c -= n3;
        nremoved += n3;

        if (nremoved >= c)
        {
            n3 -= (nremoved-c);
            nremoved = c;
        }

        tri->setCount(slidx, n3);

	}
	assert(nremoved <= c);
	c -= nremoved;

	while(c != 0)
	{
		Tri * tri = patch->pickTriByArea(rng()->getUnfIE());
		assert (tri != 0);
		tri->setCount(slidx, (tri->pools()[slidx] + 1.0));
		c--;
	}

	for (TriPVecCI t = patch->bgnTri(); t != t_end; ++t)
	{
		_updateSpec(*t, slidx);
	}

	// Rates have changed
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getPatchAmount(uint pidx, uint sidx) const
{
	// the following method does all the necessary argument checking
	double count = _getPatchCount(pidx, sidx);
	return (count / steps::math::AVOGADRO);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setPatchAmount(uint pidx, uint sidx, double a)
{
	assert(a >= 0.0);
	// convert amount in mols to number of molecules
	double a2 = a * steps::math::AVOGADRO;
	// the following method does all the necessary argument checking
	_setPatchCount(pidx, sidx, a2);
}

////////////////////////////////////////////////////////////////////////////////

bool stex::Tetexact::_getPatchClamped(uint pidx, uint sidx) const
{
	assert (pidx < statedef()->countPatches());
	assert (sidx < statedef()->countSpecs());
	assert(statedef()->countPatches() == pPatches.size());
	stex::Patch * patch = _patch(pidx);
	assert(patch != 0);
	uint lsidx = patch->def()->specG2L(sidx);
    if (lsidx == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "Species undefined in patch.\n";
		throw steps::ArgErr(os.str());
	}

    TriPVecCI t_end = patch->endTri();
    for (TriPVecCI t = patch->bgnTri(); t != t_end; ++t)
    {
        if ((*t)->clamped(lsidx) == false) return false;
    }
    return true;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setPatchClamped(uint pidx, uint sidx, bool buf)
{
	assert (pidx < statedef()->countPatches());
	assert (sidx < statedef()->countSpecs());
	assert(statedef()->countPatches() == pPatches.size());
	stex::Patch * patch = _patch(pidx);
	assert(patch != 0);
	uint lsidx = patch->def()->specG2L(sidx);
    if (lsidx == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "Species undefined in patch.\n";
		throw steps::ArgErr(os.str());
	}

    // Set the flag in def object for consistency, though this is not
    // entirely necessary
    patch->def()->setClamped(lsidx, buf);

    TriPVecCI t_end = patch->endTri();
    for (TriPVecCI t = patch->bgnTri(); t != t_end; ++t)
    {
        (*t)->setClamped(lsidx, buf);
    }

}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getPatchSReacK(uint pidx, uint ridx) const
{
	assert (pidx < statedef()->countPatches());
	assert (ridx < statedef()->countSReacs());
	assert(statedef()->countPatches() == pPatches.size());
	stex::Patch * patch = _patch(pidx);
	assert(patch != 0);
	uint lsridx = patch->def()->sreacG2L(ridx);
    if (lsridx == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "Surface reaction undefined in patch.\n";
		throw steps::ArgErr(os.str());
	}

    // We're just returning the default value for this patch, individual
    // triangles may have different Kcsts set
    return (patch->def()->kcst(lsridx));
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setPatchSReacK(uint pidx, uint ridx, double kf)
{
	assert (pidx < statedef()->countPatches());
	assert (ridx < statedef()->countSReacs());
	assert(statedef()->countPatches() == pPatches.size());
	assert(kf >= 0.0);
	stex::Patch * patch = _patch(pidx);
	assert(patch != 0);
	uint lsridx = patch->def()->sreacG2L(ridx);
    if (lsridx == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "Surface reaction undefined in patch.\n";
		throw steps::ArgErr(os.str());
	}

    // First set the default values for this patch
    patch->def()->setKcst(lsridx, kf);

    // Now update all triangles in this patch
    TriPVecCI t_end = patch->endTri();
    for (TriPVecCI t = patch->bgnTri(); t != t_end; ++t)
    {
        (*t)->sreac(lsridx)->setKcst(kf);
    }

	// Rates have changed
	_update();
}

////////////////////////////////////////////////////////////////////////////////

bool stex::Tetexact::_getPatchSReacActive(uint pidx, uint ridx) const
{
	assert (pidx < statedef()->countPatches());
	assert (ridx < statedef()->countSReacs());
	assert(statedef()->countPatches() == pPatches.size());
	stex::Patch * patch = _patch(pidx);
	assert(patch != 0);
	uint lsridx = patch->def()->sreacG2L(ridx);
    if (lsridx == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "Surface reaction undefined in patch.\n";
		throw steps::ArgErr(os.str());
	}

    TriPVecCI t_end = patch->endTri();
    for (TriPVecCI t = patch->bgnTri(); t != t_end; ++t)
    {
        if ((*t)->sreac(lsridx)->inactive() == true) return false;
    }
    return true;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setDiffBoundaryDiffusionActive(uint dbidx, uint sidx, bool act)
{
	assert (dbidx < statedef()->countDiffBoundaries());
	assert (sidx < statedef()->countSpecs());

	// Need to do two things:
	// 1) check if the species is defined in both compartments conencted
	// by the diffusion boundary
	// 2) loop over all tetrahedrons around the diff boundary and then the
	// diffusion rules and activate diffusion if the diffusion rule
	// relates to this species

	stex::DiffBoundary * diffb = _diffboundary(dbidx);
	stex::Comp * compA = diffb->compA();
	stex::Comp * compB = diffb->compB();

	/*
	ssolver::Diffdef * diffdef = statedef()->diffdef(didx);
	uint specgidx = diffdef->lig();
	*/
	uint lsidxA = compA->def()->specG2L(sidx);
	uint lsidxB = compB->def()->specG2L(sidx);


	if (lsidxA == ssolver::LIDX_UNDEFINED or lsidxB == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "Species undefined in compartments connected by diffusion boundary.\n";
		throw steps::ArgErr(os.str());
	}

	std::vector<uint> bdtets = diffb->getTets();
	std::vector<uint> bdtetsdir = diffb->getTetDirection();

	// Have to use indices rather than iterator because need access to the
	// tet direction
	uint ntets = bdtets.size();

	for (uint bdt = 0; bdt != ntets; ++bdt)
	{
		stex::Tet * tet = _tet(bdtets[bdt]);
		uint direction = bdtetsdir[bdt];
		assert(direction >= 0 and direction < 4);

		// Each diff kproc then has access to the species through it's defined parent
		uint ndiffs = tet->compdef()->countDiffs();
		for (uint d = 0; d != ndiffs; ++d)
		{
			stex::Diff * diff = tet->diff(d);
			// sidx is the global species index; so is the lig() return from diffdef
			uint specgidx = diff->def()->lig();
			if (specgidx == sidx)
			{
				diff->setDiffBndActive(direction, act);
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////////

bool stex::Tetexact::_getDiffBoundaryDiffusionActive(uint dbidx, uint sidx) const
{
	assert (dbidx < statedef()->countDiffBoundaries());
	assert (sidx < statedef()->countSpecs());

	stex::DiffBoundary * diffb = _diffboundary(dbidx);
	stex::Comp * compA = diffb->compA();
	stex::Comp * compB = diffb->compB();

	uint lsidxA = compA->def()->specG2L(sidx);
	uint lsidxB = compB->def()->specG2L(sidx);

	if (lsidxA == ssolver::LIDX_UNDEFINED or lsidxB == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "Species undefined in compartments connected by diffusion boundary.\n";
		throw steps::ArgErr(os.str());
	}

	std::vector<uint> bdtets = diffb->getTets();
	std::vector<uint> bdtetsdir = diffb->getTetDirection();

	// Have to use indices rather than iterator because need access to the
	// tet direction

	uint ntets = bdtets.size();

	for (uint bdt = 0; bdt != ntets; ++bdt)
	{
		stex::Tet * tet = _tet(bdtets[bdt]);
		uint direction = bdtetsdir[bdt];
		assert(direction >= 0 and direction < 4);

		// Each diff kproc then has access to the species through it's defined parent
		uint ndiffs = tet->compdef()->countDiffs();
		for (uint d = 0; d != ndiffs; ++d)
		{
			stex::Diff * diff = tet->diff(d);
			// sidx is the global species index; so is the lig() return from diffdef
			uint specgidx = diff->def()->lig();
			if (specgidx == sidx)
			{
				// Just need to check the first one
				if (diff->getDiffBndActive(direction)) return true;
				else return false;
			}
		}
	}
	return false;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setPatchSReacActive(uint pidx, uint ridx, bool a)
{
	assert (pidx < statedef()->countPatches());
	assert (ridx < statedef()->countSReacs());
	assert(statedef()->countPatches() == pPatches.size());
	stex::Patch * patch = _patch(pidx);
	assert(patch != 0);
	uint lsridx = patch->def()->sreacG2L(ridx);
    if (lsridx == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "Surface reaction undefined in patch.\n";
		throw steps::ArgErr(os.str());
	}

    // First set the flags in def object for consistency, though this is
    // not entirely necessary for this solver
    patch->def()->setActive(lsridx, a);

    TriPVecCI t_end = patch->endTri();
    for (TriPVecCI t = patch->bgnTri(); t != t_end; ++t)
    {
        (*t)->sreac(lsridx)->setActive(a);
    }
    // It's cheaper to just recompute everything.
    _update();
}

////////////////////////////////////////////////////////////////////////////////

bool stex::Tetexact::_getPatchVDepSReacActive(uint pidx, uint vsridx) const
{
	assert (pidx < statedef()->countPatches());
	assert (vsridx < statedef()->countVDepSReacs());
	assert(statedef()->countPatches() == pPatches.size());
	stex::Patch * patch = _patch(pidx);
	assert(patch != 0);
	uint lvsridx = patch->def()->vdepsreacG2L(vsridx);
    if (lvsridx == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "Voltage-dependent surface reaction undefined in patch.\n";
		throw steps::ArgErr(os.str());
	}

    TriPVecCI t_end = patch->endTri();
    for (TriPVecCI t = patch->bgnTri(); t != t_end; ++t)
    {
        if ((*t)->vdepsreac(lvsridx)->inactive() == true) return false;
    }
    return true;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setPatchVDepSReacActive(uint pidx, uint vsridx, bool a)
{
	assert (pidx < statedef()->countPatches());
	assert (vsridx < statedef()->countVDepSReacs());
	assert(statedef()->countPatches() == pPatches.size());
	stex::Patch * patch = _patch(pidx);
	assert(patch != 0);
	uint lvsridx = patch->def()->vdepsreacG2L(vsridx);
    if (lvsridx == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "Voltage-dependent surface reaction undefined in patch.\n";
		throw steps::ArgErr(os.str());
	}

    // Not necessary and not possible to set the flags in def object

    TriPVecCI t_end = patch->endTri();
    for (TriPVecCI t = patch->bgnTri(); t != t_end; ++t)
    {
        (*t)->vdepsreac(lvsridx)->setActive(a);
    }
    // It's cheaper to just recompute everything.
    _update();
}
////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::addKProc(steps::tetexact::KProc * kp)
{
	assert (kp != 0);

	SchedIDX nidx = pKProcs.size();
	pKProcs.push_back(kp);
	kp->setSchedIDX(nidx);
}

////////////////////////////////////////////////////////////////////////////////
/*
void stex::Tetexact::_build(void)
{
	assert (pBuilt == false);

    pBuilt = true;
}
*/
////////////////////////////////////////////////////////////////////////////////

steps::tetexact::KProc * stex::Tetexact::_getNext(void) const
{
#ifdef SSA_DEBUG
    std::cout << "SSA: Search for next event\n";
#endif
    assert(pA0 >= 0.0);
    // Quick check to see whether nothing is there.
    if (pA0 == 0.0) return NULL;

    double selector = pA0 * rng()->getUnfII();

#ifdef SSA_DEBUG
    std::cout << "selector: " << selector << " in total " << pA0 << "\n";
#endif

    double partial_sum = 0.0;

    uint n_neg_groups = nGroups.size();
    uint n_pos_groups = pGroups.size();

    for (uint i = 0; i < n_neg_groups; i++) {
        CRGroup* group = nGroups[i];
        if (group->size == 0) continue;

        if (selector > partial_sum + group->sum) {
            partial_sum += group->sum;
#ifdef SSA_DEBUG
            std::cout << "increase partial sum to " << partial_sum;
            std::cout << " by neg group " << i << "\n";
#endif
            continue;
        }

        double g_max = group->max;
        double random_rate = g_max * rng()->getUnfII();;
        uint group_size = group->size;
        uint random_pos = rng()->get() % group_size;
        KProc* random_kp = group->indices[random_pos];

#ifdef SSA_DEBUG
        std::cout << "search for event\n";
        std::cout << "group max: " << g_max << "\n";
        std::cout << "random rate: " << random_rate << "\n";
        std::cout << "random pos " << random_pos << "\n";
        std::cout << "random kp rate " << random_kp->crData.rate << "\n";
#endif

        while (random_kp->crData.rate <= random_rate) {
            random_rate = g_max * rng()->getUnfII();
            random_pos = rng()->get() % group_size;
            random_kp = group->indices[random_pos];
#ifdef SSA_DEBUG
            std::cout << "renew search\n";
            std::cout << "random rate: " << random_rate << "\n";
            std::cout << "random pos " << random_pos << "\n";
            std::cout << "random kp rate " << random_kp->crData.rate << "\n";
#endif
        }

#ifdef SSA_DEBUG
        std::cout << "selected kp index: " << random_kp->schedIDX() << "\n";
        std::cout << "--------------------------------------------------------\n";
#endif
        return random_kp;
    }

    for (uint i = 0; i < n_pos_groups; i++) {
        CRGroup* group = pGroups[i];
        if (group->size == 0) continue;

        if (selector > partial_sum + group->sum) {
            partial_sum += group->sum;
#ifdef SSA_DEBUG
            std::cout << "increase partial sum to " << partial_sum;
            std::cout << " by neg group " << i << "\n";
#endif
            continue;
        }

        double g_max = group->max;
        double random_rate = g_max * rng()->getUnfII();;
        uint group_size = group->size;
        uint random_pos = rng()->get() % group_size;
        KProc* random_kp = group->indices[random_pos];

#ifdef SSA_DEBUG
        std::cout << "search for event\n";
        std::cout << "group max: " << g_max << "\n";
        std::cout << "random rate: " << random_rate << "\n";
        std::cout << "random pos " << random_pos << "\n";
        std::cout << "random kp rate " << random_kp->crData.rate << "\n";
#endif

        while (random_kp->crData.rate <= random_rate) {
            random_rate = g_max * rng()->getUnfII();
            random_pos = rng()->get() % group_size;
            random_kp = group->indices[random_pos];
#ifdef SSA_DEBUG
            std::cout << "renew search\n";
            std::cout << "random rate: " << random_rate << "\n";
            std::cout << "random pos " << random_pos << "\n";
            std::cout << "random kp rate " << random_kp->crData.rate << "\n";
#endif
        }

#ifdef SSA_DEBUG
        std::cout << "selected kp index: " << random_kp->schedIDX() << "\n";
        std::cout << "--------------------------------------------------------\n";
#endif
        return random_kp;
    }


    std::cerr << "Cannot find any suitable entry.\n";
    std::cerr << "A0: " << std::setprecision (15) << pA0 << "\n";
    std::cerr << "Selector: " << std::setprecision (15) << selector << "\n";
    std::cerr << "Current Partial Sum: " << std::setprecision (15) << partial_sum << "\n";

    std::cerr << "Distribution of group sums\n";
    std::cerr << "Negative groups\n";

    for (uint i = 0; i < n_neg_groups; i++) {
        std::cerr << i << ": " << std::setprecision (15) << nGroups[i]->sum << "\n";
    }
    std::cerr << "Positive groups\n";
    for (uint i = 0; i < n_pos_groups; i++) {
        std::cerr << i << ": " << std::setprecision (15) << pGroups[i]->sum << "\n";
    }

    throw;
}

////////////////////////////////////////////////////////////////////////////////
/*
void stex::Tetexact::_reset(void)
{

}
*/
////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_executeStep(steps::tetexact::KProc * kp, double dt)
{
    std::vector<KProc*> const & upd = kp->apply(rng(), dt, statedef()->time());
    _update(upd);
    statedef()->incTime(dt);
    statedef()->incNSteps(1);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_updateSpec(steps::tetexact::WmVol * tet, uint spec_lidx)
{
    std::set<KProc*> updset;

    // Loop over tet.
    KProcPVecCI kproc_end = tet->kprocEnd();
    for (KProcPVecCI k = tet->kprocBegin(); k != kproc_end; ++k)
    {
        updset.insert(*k);
    }

    std::vector<stex::Tri *>::const_iterator tri_end = tet->nexttriEnd();
    for (std::vector<stex::Tri *>::const_iterator tri = tet->nexttriBegin();
			 tri != tri_end; ++tri)
    {
    	if ((*tri) == 0) continue;
        kproc_end = (*tri)->kprocEnd();
        for (KProcPVecCI k = (*tri)->kprocBegin(); k != kproc_end; ++k)
        {
            updset.insert(*k);
        }
    }

    // Send the list of kprocs that need to be updated to the schedule.
    if (updset.empty()) return;
    std::vector<KProc*> updvec(updset.begin(), updset.end());

    _update(updvec);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_updateSpec(steps::tetexact::Tri * tri, uint spec_lidx)
{
    std::set<stex::KProc*> updset;

    KProcPVecCI kproc_end = tri->kprocEnd();
    for (KProcPVecCI k = tri->kprocBegin(); k != kproc_end; ++k)
    {
        updset.insert(*k);
    }

    // Send the list of kprocs that need to be updated to the schedule.
    if (updset.empty()) return;

    std::vector<stex::KProc*> updvec(updset.begin(), updset.end());

    _update(updvec);

}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getCompReacH(uint cidx, uint ridx) const
{
	assert(cidx < statedef()->countComps());
	assert(ridx < statedef()->countReacs());
	ssolver::Compdef * comp = statedef()->compdef(cidx);
	assert(comp != 0);
	uint lridx = comp->reacG2L(ridx);
	if (lridx == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "Reaction undefined in compartment.\n";
		throw steps::ArgErr(os.str());
	}

	// The 'local' Comp object has same index as solver::Compdef object
	stex::Comp * lcomp = pComps[cidx];
	assert (lcomp->def() == comp);

	WmVolPVecCI t_bgn = lcomp->bgnTet();
	WmVolPVecCI t_end = lcomp->endTet();
	if (t_bgn == t_end) return 0.0;

	double h = 0.0;
	for (WmVolPVecCI t = t_bgn; t != t_end; ++t)
	{
		stex::Reac * reac = (*t)->reac(lridx);
		h += reac->h();
	}

	return h;
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getCompReacC(uint cidx, uint ridx) const
{
	assert(cidx < statedef()->countComps());
	assert(ridx < statedef()->countReacs());
	ssolver::Compdef * comp = statedef()->compdef(cidx);
	assert(comp != 0);
	uint lridx = comp->reacG2L(ridx);
	if (lridx == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "Reaction undefined in compartment.\n";
		throw steps::ArgErr(os.str());
	}

	// The 'local' Comp object has same index as solver::Compdef object
	stex::Comp * lcomp = pComps[cidx];
	assert (lcomp->def() == comp);

	WmVolPVecCI t_bgn = lcomp->bgnTet();
	WmVolPVecCI t_end = lcomp->endTet();
	if (t_bgn == t_end) return 0.0;
	double c2 = 0.0;
	double v = 0.0;
	for (WmVolPVecCI t = t_bgn; t != t_end; ++t)
	{
		double v2 = (*t)->vol();
		stex::Reac * reac = (*t)->reac(lridx);
		c2 += (reac->c() * v2);
		v += v2;
	}
	assert(v > 0.0);
	return c2/v;
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getCompReacA(uint cidx, uint ridx) const
{
	assert(cidx < statedef()->countComps());
	assert(ridx < statedef()->countReacs());
	ssolver::Compdef * comp = statedef()->compdef(cidx);
	assert(comp != 0);
	uint lridx = comp->reacG2L(ridx);
	if (lridx == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "Reaction undefined in compartment.\n";
		throw steps::ArgErr(os.str());
	}

	// The 'local' Comp object has same index as solver::Compdef object
	stex::Comp * lcomp = pComps[cidx];
	assert (lcomp->def() == comp);

	WmVolPVecCI t_bgn = lcomp->bgnTet();
	WmVolPVecCI t_end = lcomp->endTet();
	if (t_bgn == t_end) return 0.0;

	double a = 0.0;
	for (WmVolPVecCI t = t_bgn; t != t_end; ++t)
	{
		stex::Reac * reac = (*t)->reac(lridx);
		a += reac->rate();
	}

	return a;
}

////////////////////////////////////////////////////////////////////////////////

uint stex::Tetexact::_getCompReacExtent(uint cidx, uint ridx) const
{
	assert(cidx < statedef()->countComps());
	assert(ridx < statedef()->countReacs());
	ssolver::Compdef * comp = statedef()->compdef(cidx);
	assert(comp != 0);
	uint lridx = comp->reacG2L(ridx);
	if (lridx == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "Reaction undefined in compartment.\n";
		throw steps::ArgErr(os.str());
	}

	// The 'local' Comp object has same index as solver::Compdef object
	stex::Comp * lcomp = pComps[cidx];
	assert (lcomp->def() == comp);

	WmVolPVecCI t_bgn = lcomp->bgnTet();
	WmVolPVecCI t_end = lcomp->endTet();
	if (t_bgn == t_end) return 0.0;

	uint x = 0.0;
	for (WmVolPVecCI t = t_bgn; t != t_end; ++t)
	{
		stex::Reac * reac = (*t)->reac(lridx);
		x += reac->getExtent();
	}

	return x;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_resetCompReacExtent(uint cidx, uint ridx)
{
	assert(cidx < statedef()->countComps());
	assert(ridx < statedef()->countReacs());
	ssolver::Compdef * comp = statedef()->compdef(cidx);
	assert(comp != 0);
	uint lridx = comp->reacG2L(ridx);
	if (lridx == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "Reaction undefined in compartment.\n";
		throw steps::ArgErr(os.str());
	}

	// The 'local' Comp object has same index as solver::Compdef object
	stex::Comp * lcomp = pComps[cidx];
	assert (lcomp->def() == comp);

	WmVolPVecCI t_bgn = lcomp->bgnTet();
	WmVolPVecCI t_end = lcomp->endTet();
	if (t_bgn == t_end) return;

	for (WmVolPVecCI t = t_bgn; t != t_end; ++t)
	{
		stex::Reac * reac = (*t)->reac(lridx);
		reac->resetExtent();
	}
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getPatchSReacH(uint pidx, uint ridx) const
{
	assert (pidx < statedef()->countPatches());
	assert (ridx < statedef()->countSReacs());
	ssolver::Patchdef * patch = statedef()->patchdef(pidx);
	uint lsridx = patch->sreacG2L(ridx);
	if (lsridx == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "Surface reaction undefined in patch.\n";
		throw steps::ArgErr(os.str());
	}

	// The 'local' Patch object has same index as solver::Patchdef object
	stex::Patch * lpatch = pPatches[pidx];
	assert (lpatch->def() == patch);

	TriPVecCI t_bgn = lpatch->bgnTri();
	TriPVecCI t_end = lpatch->endTri();
	if (t_bgn == t_end) return 0.0;

	double h = 0.0;
	for (TriPVecCI t = t_bgn; t != t_end; ++t)
	{
		stex::SReac * sreac = (*t)->sreac(lsridx);
		h += sreac->h();
	}

	return h;
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getPatchSReacC(uint pidx, uint ridx) const
{
	assert (pidx < statedef()->countPatches());
	assert (ridx < statedef()->countSReacs());
	ssolver::Patchdef * patch = statedef()->patchdef(pidx);
	uint lsridx = patch->sreacG2L(ridx);
	if (lsridx == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "Surface reaction undefined in patch.\n";
		throw steps::ArgErr(os.str());
	}

	// The 'local' Patch object has same index as solver::Patchdef object
	stex::Patch * lpatch = pPatches[pidx];
	assert (lpatch->def() == patch);

	TriPVecCI t_bgn = lpatch->bgnTri();
	TriPVecCI t_end = lpatch->endTri();
	if (t_bgn == t_end) return 0.0;

	double c = 0.0;
	double a = 0.0;
	for (TriPVecCI t = t_bgn; t != t_end; ++t)
	{
		double a2 = (*t)->area();
		stex::SReac * sreac = (*t)->sreac(lsridx);
		c += (sreac->c() * a2);
		a += a2;
	}
	assert (a > 0.0);
	return c/a;
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getPatchSReacA(uint pidx, uint ridx) const
{
	assert (pidx < statedef()->countPatches());
	assert (ridx < statedef()->countSReacs());
	ssolver::Patchdef * patch = statedef()->patchdef(pidx);
	uint lsridx = patch->sreacG2L(ridx);
	if (lsridx == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "Surface reaction undefined in patch.\n";
		throw steps::ArgErr(os.str());
	}

	// The 'local' Patch object has same index as solver::Patchdef object
	stex::Patch * lpatch = pPatches[pidx];
	assert (lpatch->def() == patch);

	TriPVecCI t_bgn = lpatch->bgnTri();
	TriPVecCI t_end = lpatch->endTri();
	if (t_bgn == t_end) return 0.0;

	double a = 0.0;
	for (TriPVecCI t = t_bgn; t != t_end; ++t)
	{
		stex::SReac * sreac = (*t)->sreac(lsridx);
		a += sreac->rate();
	}
	return a;
}

////////////////////////////////////////////////////////////////////////////////

uint stex::Tetexact::_getPatchSReacExtent(uint pidx, uint ridx) const
{
	assert (pidx < statedef()->countPatches());
	assert (ridx < statedef()->countSReacs());
	ssolver::Patchdef * patch = statedef()->patchdef(pidx);
	uint lsridx = patch->sreacG2L(ridx);
	if (lsridx == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "Surface reaction undefined in patch.\n";
		throw steps::ArgErr(os.str());
	}

	// The 'local' Patch object has same index as solver::Patchdef object
	stex::Patch * lpatch = pPatches[pidx];
	assert (lpatch->def() == patch);

	TriPVecCI t_bgn = lpatch->bgnTri();
	TriPVecCI t_end = lpatch->endTri();
	if (t_bgn == t_end) return 0;

	double x = 0.0;
	for (TriPVecCI t = t_bgn; t != t_end; ++t)
	{
		stex::SReac * sreac = (*t)->sreac(lsridx);
		x += sreac->getExtent();
	}
	return x;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_resetPatchSReacExtent(uint pidx, uint ridx)
{
	assert (pidx < statedef()->countPatches());
	assert (ridx < statedef()->countSReacs());
	ssolver::Patchdef * patch = statedef()->patchdef(pidx);
	uint lsridx = patch->sreacG2L(ridx);
	if (lsridx == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "Surface reaction undefined in patch.\n";
		throw steps::ArgErr(os.str());
	}

	// The 'local' Patch object has same index as solver::Patchdef object
	stex::Patch * lpatch = pPatches[pidx];
	assert (lpatch->def() == patch);

	TriPVecCI t_bgn = lpatch->bgnTri();
	TriPVecCI t_end = lpatch->endTri();
	if (t_bgn == t_end) return;

	for (TriPVecCI t = t_bgn; t != t_end; ++t)
	{
		stex::SReac * sreac = (*t)->sreac(lsridx);
		sreac->resetExtent();
	}
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getTetVol(uint tidx) const
{
	assert (tidx < pTets.size());
	if (pTets[tidx] == 0)
	{
    	std::ostringstream os;
    	os << "Tetrahedron " << tidx << " has not been assigned to a compartment.";
    	throw steps::ArgErr(os.str());
	}
	return pTets[tidx]->vol();
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setTetVol(uint tidx, double vol)
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

bool stex::Tetexact::_getTetSpecDefined(uint tidx, uint sidx) const
{
	assert (tidx < pTets.size());
	assert (sidx < statedef()->countSpecs());

	if (pTets[tidx] == 0) return false;

	stex::Tet * tet = pTets[tidx];
	uint lsidx = tet->compdef()->specG2L(sidx);
	if (lsidx == ssolver::LIDX_UNDEFINED) return false;
	else return true;
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getTetCount(uint tidx, uint sidx) const
{
	assert (tidx < pTets.size());
	assert (sidx < statedef()->countSpecs());

	if (pTets[tidx] == 0)
	{
    	std::ostringstream os;
    	os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
    	throw steps::ArgErr(os.str());
	}

	stex::Tet * tet = pTets[tidx];
	uint lsidx = tet->compdef()->specG2L(sidx);
	if (lsidx == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "Species undefined in tetrahedron.\n";
		throw steps::ArgErr(os.str());
	}

	return tet->pools()[lsidx];
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setTetCount(uint tidx, uint sidx, double n)
{
	assert (tidx < pTets.size());
	assert (sidx < statedef()->countSpecs());
	assert (n >= 0.0);

	if (pTets[tidx] == 0)
	{
    	std::ostringstream os;
    	os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
    	throw steps::ArgErr(os.str());
	}
	if (n > std::numeric_limits<unsigned int>::max( ))
	{
		std::ostringstream os;
		os << "Can't set count greater than maximum unsigned integer (";
		os << std::numeric_limits<unsigned int>::max( ) << ").\n";
		throw steps::ArgErr(os.str());
	}

	stex::Tet * tet = pTets[tidx];

	uint lsidx = tet->compdef()->specG2L(sidx);
	if (lsidx == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "Species undefined in tetrahedron.\n";
		throw steps::ArgErr(os.str());
	}

	double n_int = std::floor(n);
	double n_frc = n - n_int;
	uint c = static_cast<uint>(n_int);
	if (n_frc > 0.0)
	{
		double rand01 = rng()->getUnfIE();
		if (rand01 < n_frc) c++;
	}

	// Tet object updates def level Comp object counts
	tet->setCount(lsidx, c);
	_updateSpec(tet, lsidx);
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getTetAmount(uint tidx, uint sidx) const
{
	// following method does all necessary argument checking
	double count = _getTetCount(tidx, sidx);
	return count/steps::math::AVOGADRO;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setTetAmount(uint tidx, uint sidx, double m)
{
	// convert amount in mols to number of molecules
	double m2 = m * steps::math::AVOGADRO;
	// the following method does all the necessary argument checking
	_setTetCount(tidx, sidx, m2);
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getTetConc(uint tidx, uint sidx) const
{
	// following method does all necessary argument checking
	double count = _getTetCount(tidx, sidx);
	stex::Tet * tet = pTets[tidx];
	double vol = tet->vol();
	return (count/(1.0e3 * vol * steps::math::AVOGADRO));
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setTetConc(uint tidx, uint sidx, double c)
{
	assert (c >= 0.0);
	assert (tidx < pTets.size());

	if (pTets[tidx] == 0)
	{
    	std::ostringstream os;
    	os << "Tetrahedron " << tidx << " has not been assigned to a compartment.";
    	throw steps::ArgErr(os.str());
	}

	stex::Tet * tet = pTets[tidx];
	double count = c * (1.0e3 * tet->vol() * steps::math::AVOGADRO);
	// the following method does all the necessary argument checking
	_setTetCount(tidx, sidx, count);
}

////////////////////////////////////////////////////////////////////////////////

bool stex::Tetexact::_getTetClamped(uint tidx, uint sidx) const
{
	assert (tidx < pTets.size());
	assert (sidx < statedef()->countSpecs());

	if (pTets[tidx] == 0)
	{
    	std::ostringstream os;
    	os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
    	throw steps::ArgErr(os.str());
	}

	stex::Tet * tet = pTets[tidx];

	uint lsidx = tet->compdef()->specG2L(sidx);
	if (lsidx == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "Species undefined in tetrahedron.\n";
		throw steps::ArgErr(os.str());
	}

	return tet->clamped(lsidx);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setTetClamped(uint tidx, uint sidx, bool buf)
{
	assert (tidx < pTets.size());
	assert (sidx < statedef()->countSpecs());

	if (pTets[tidx] == 0)
	{
    	std::ostringstream os;
    	os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
    	throw steps::ArgErr(os.str());
	}

	stex::Tet * tet = pTets[tidx];

	uint lsidx = tet->compdef()->specG2L(sidx);
	if (lsidx == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "Species undefined in tetrahedron.\n";
		throw steps::ArgErr(os.str());
	}

	tet->setClamped(lsidx, buf);
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getTetReacK(uint tidx, uint ridx) const
{
	assert (tidx < pTets.size());
	assert (ridx < statedef()->countReacs());

	if (pTets[tidx] == 0)
	{
    	std::ostringstream os;
    	os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
    	throw steps::ArgErr(os.str());
	}

	stex::Tet * tet = pTets[tidx];

	uint lridx = tet->compdef()->reacG2L(ridx);
	if (lridx == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "Reaction undefined in tetrahedron.\n";
		throw steps::ArgErr(os.str());
	}

	return (tet->reac(lridx)->kcst());
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setTetReacK(uint tidx, uint ridx, double kf)
{
	assert (tidx < pTets.size());
	assert (ridx < statedef()->countReacs());
	assert (kf >= 0.0);

	if (pTets[tidx] == 0)
	{
    	std::ostringstream os;
    	os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
    	throw steps::ArgErr(os.str());
	}

	stex::Tet * tet = pTets[tidx];

	uint lridx = tet->compdef()->reacG2L(ridx);
	if (lridx == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "\nReaction undefined in tetrahedron.";
		throw steps::ArgErr(os.str());
	}

	tet->reac(lridx)->setKcst(kf);

	_updateElement(tet->reac(lridx));
    _updateSum();
}

////////////////////////////////////////////////////////////////////////////////

bool stex::Tetexact::_getTetReacActive(uint tidx, uint ridx) const
{
	assert (tidx < pTets.size());
	assert (ridx < statedef()->countReacs());

	if (pTets[tidx] == 0)
	{
    	std::ostringstream os;
    	os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
    	throw steps::ArgErr(os.str());
	}

	stex::Tet * tet = pTets[tidx];

	uint lridx = tet->compdef()->reacG2L(ridx);
	if (lridx == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "Reaction undefined in tetrahedron.\n";
		throw steps::ArgErr(os.str());
	}

	if (tet->reac(lridx)->inactive() == true) return false;
	return true;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setTetReacActive(uint tidx, uint ridx, bool act)
{
	assert (tidx < pTets.size());
	assert (ridx < statedef()->countReacs());

	if (pTets[tidx] == 0)
	{
    	std::ostringstream os;
    	os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
    	throw steps::ArgErr(os.str());
	}

	stex::Tet * tet = pTets[tidx];

	uint lridx = tet->compdef()->reacG2L(ridx);
	if (lridx == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "Reaction undefined in tetrahedron.\n";
		throw steps::ArgErr(os.str());
	}

	tet->reac(lridx)->setActive(act);

	_updateElement(tet->reac(lridx));
    _updateSum();
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getTetDiffD(uint tidx, uint didx) const
{
	assert (tidx < pTets.size());
	assert (didx < statedef()->countDiffs());

	if (pTets[tidx] == 0)
	{
    	std::ostringstream os;
    	os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
    	throw steps::ArgErr(os.str());
	}

	stex::Tet * tet = pTets[tidx];

	uint ldidx = tet->compdef()->diffG2L(didx);
	if (ldidx == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "Diffusion rule undefined in tetrahedron.\n";
		throw steps::ArgErr(os.str());
	}

	return (tet->diff(ldidx)->dcst());

}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setTetDiffD(uint tidx, uint didx, double dk)
{
	assert (tidx < pTets.size());
	assert (didx < statedef()->countDiffs());

	if (pTets[tidx] == 0)
	{
    	std::ostringstream os;
    	os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
    	throw steps::ArgErr(os.str());
	}

	stex::Tet * tet = pTets[tidx];

	uint ldidx = tet->compdef()->diffG2L(didx);
	if (ldidx == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "Diffusion rule undefined in tetrahedron.\n";
		throw steps::ArgErr(os.str());
	}

	tet->diff(ldidx)->setDcst(dk);

	_updateElement(tet->diff(ldidx));
    _updateSum();
}

////////////////////////////////////////////////////////////////////////////////

bool stex::Tetexact::_getTetDiffActive(uint tidx, uint didx) const
{
	assert (tidx < pTets.size());
	assert (didx < statedef()->countDiffs());

	if (pTets[tidx] == 0)
	{
    	std::ostringstream os;
    	os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
    	throw steps::ArgErr(os.str());
	}

	stex::Tet * tet = pTets[tidx];

	uint ldidx = tet->compdef()->diffG2L(didx);
	if (ldidx == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "Diffusion rule undefined in tetrahedron.\n";
		throw steps::ArgErr(os.str());
	}

	if (tet->diff(ldidx)->inactive() == true) return false;
	return true;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setTetDiffActive(uint tidx, uint didx, bool act)
{
	assert (tidx < pTets.size());
	assert (didx < statedef()->countDiffs());

	if (pTets[tidx] == 0)
	{
    	std::ostringstream os;
    	os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
    	throw steps::ArgErr(os.str());
	}

	stex::Tet * tet = pTets[tidx];

	uint ldidx = tet->compdef()->diffG2L(didx);
	if (ldidx == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "Diffusion rule undefined in tetrahedron.\n";
		throw steps::ArgErr(os.str());
	}

	tet->diff(ldidx)->setActive(act);

	_updateElement(tet->diff(ldidx));
    _updateSum();
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getTetReacH(uint tidx, uint ridx) const
{
	assert (tidx < pTets.size());
	assert (ridx < statedef()->countReacs());

	if (pTets[tidx] == 0)
	{
    	std::ostringstream os;
    	os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
    	throw steps::ArgErr(os.str());
	}

	stex::Tet * tet = pTets[tidx];

	uint lridx = tet->compdef()->reacG2L(ridx);
	if (lridx == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "Reaction undefined in tetrahedron.\n";
		throw steps::ArgErr(os.str());
	}

	return tet->reac(lridx)->h();
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getTetReacC(uint tidx, uint ridx) const
{
	assert (tidx < pTets.size());
	assert (ridx < statedef()->countReacs());

	if (pTets[tidx] == 0)
	{
    	std::ostringstream os;
    	os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
    	throw steps::ArgErr(os.str());
	}

	stex::Tet * tet = pTets[tidx];

	uint lridx = tet->compdef()->reacG2L(ridx);
	if (lridx == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "Reaction undefined in tetrahedron.\n";
		throw steps::ArgErr(os.str());
	}

	return tet->reac(lridx)->c();
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getTetReacA(uint tidx, uint ridx) const
{
	assert (tidx < pTets.size());
	assert (ridx < statedef()->countReacs());

	if (pTets[tidx] == 0)
	{
    	std::ostringstream os;
    	os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
    	throw steps::ArgErr(os.str());
	}

	stex::Tet * tet = pTets[tidx];

	uint lridx = tet->compdef()->reacG2L(ridx);
	if (lridx == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "Reaction undefined in tetrahedron.\n";
		throw steps::ArgErr(os.str());
	}

	return tet->reac(lridx)->rate();
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getTetDiffA(uint tidx, uint didx) const
{
	assert (tidx < pTets.size());
	assert (didx < statedef()->countDiffs());

	if (pTets[tidx] == 0)
	{
    	std::ostringstream os;
    	os << "Tetrahedron " << tidx << " has not been assigned to a compartment.\n";
    	throw steps::ArgErr(os.str());
	}

	stex::Tet * tet = pTets[tidx];

	uint ldidx = tet->compdef()->diffG2L(didx);
	if (ldidx == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "Diffusion rule undefined in tetrahedron.\n";
		throw steps::ArgErr(os.str());
	}

	return tet->diff(ldidx)->rate();
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getTriArea(uint tidx) const
{
	assert (tidx < pTris.size());

	if (pTris[tidx] == 0)
	{
    	std::ostringstream os;
    	os << "Triangle " << tidx << " has not been assigned to a patch.";
    	throw steps::ArgErr(os.str());
	}

	return pTris[tidx]->area();
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setTriArea(uint tidx, double area)
{
	throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

bool stex::Tetexact::_getTriSpecDefined(uint tidx, uint sidx) const
{
	assert (tidx < pTris.size());
	assert (sidx < statedef()->countSpecs());

	if (pTris[tidx] == 0) return false;

	stex::Tri * tri = pTris[tidx];
	uint lsidx = tri->patchdef()->specG2L(sidx);
	if (lsidx == ssolver::LIDX_UNDEFINED) return false;
	else return true;
}

////////////////////////////////////////////////////////////////////////////////


double stex::Tetexact::_getTriCount(uint tidx, uint sidx) const
{
	assert (tidx < pTris.size());
	assert (sidx < statedef()->countSpecs());

	if (pTris[tidx] == 0)
	{
    	std::ostringstream os;
    	os << "Triangle " << tidx << " has not been assigned to a patch.\n";
    	throw steps::ArgErr(os.str());
	}

	stex::Tri * tri = pTris[tidx];
	uint lsidx = tri->patchdef()->specG2L(sidx);
	if (lsidx == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "Species undefined in triangle.\n";
		throw steps::ArgErr(os.str());
	}

	return tri->pools()[lsidx];
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setTriCount(uint tidx, uint sidx, double n)
{
	assert (tidx < pTris.size());
	assert (sidx < statedef()->countSpecs());
	assert (n >= 0.0);

	if (pTris[tidx] == 0)
	{
    	std::ostringstream os;
    	os << "Triangle " << tidx << " has not been assigned to a patch.\n";
    	throw steps::ArgErr(os.str());
	}
	if (n > std::numeric_limits<unsigned int>::max( ))
	{
		std::ostringstream os;
		os << "Can't set count greater than maximum unsigned integer (";
		os << std::numeric_limits<unsigned int>::max( ) << ").\n";
		throw steps::ArgErr(os.str());
	}

	stex::Tri * tri = pTris[tidx];
	uint lsidx = tri->patchdef()->specG2L(sidx);
	if (lsidx == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "Species undefined in triangle.\n";
		throw steps::ArgErr(os.str());
	}

	double n_int = std::floor(n);
	double n_frc = n - n_int;
	uint c = static_cast<uint>(n_int);
	if (n_frc > 0.0)
	{
		double rand01 = rng()->getUnfIE();
		if (rand01 < n_frc) c++;
	}

	// Tri object updates counts in def level Comp object
	tri->setCount(lsidx, c);
	_updateSpec(tri, lsidx);
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getTriAmount(uint tidx, uint sidx) const
{
	// following method does all necessary argument checking
	double count = _getTriCount(tidx, sidx);
	return count/steps::math::AVOGADRO;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setTriAmount(uint tidx, uint sidx, double m)
{
	// convert amount in mols to number of molecules
	double m2 = m * steps::math::AVOGADRO;
	// the following method does all the necessary argument checking
	_setTriCount(tidx, sidx, m2);
}

////////////////////////////////////////////////////////////////////////////////


bool stex::Tetexact::_getTriClamped(uint tidx, uint sidx) const
{
	assert (tidx < pTris.size());
	assert (sidx < statedef()->countSpecs());

	if (pTris[tidx] == 0)
	{
    	std::ostringstream os;
    	os << "Triangle " << tidx << " has not been assigned to a patch.\n";
    	throw steps::ArgErr(os.str());
	}

	stex::Tri * tri = pTris[tidx];

	uint lsidx = tri->patchdef()->specG2L(sidx);
	if (lsidx == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "Species undefined in triangle.\n";
		throw steps::ArgErr(os.str());
	}

	return tri->clamped(lsidx);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setTriClamped(uint tidx, uint sidx, bool buf)
{
	assert (tidx < pTris.size());
	assert (sidx < statedef()->countSpecs());

	if (pTris[tidx] == 0)
	{
    	std::ostringstream os;
    	os << "Triangle " << tidx << " has not been assigned to a patch.\n";
    	throw steps::ArgErr(os.str());
	}

	stex::Tri * tri = pTris[tidx];

	uint lsidx = tri->patchdef()->specG2L(sidx);
	if (lsidx == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "Species undefined in triangle.\n";
		throw steps::ArgErr(os.str());
	}

	tri->setClamped(lsidx, buf);
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getTriSReacK(uint tidx, uint ridx) const
{
	assert (tidx < pTris.size());
	assert (ridx < statedef()->countSReacs());

	if (pTris[tidx] == 0)
	{
    	std::ostringstream os;
    	os << "Triangle " << tidx << " has not been assigned to a patch.\n";
    	throw steps::ArgErr(os.str());
	}

	stex::Tri * tri = pTris[tidx];

	uint lsridx = tri->patchdef()->sreacG2L(ridx);
	if (lsridx == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "Surface reaction undefined in triangle.\n";
		throw steps::ArgErr(os.str());
	}

	return (tri->sreac(lsridx)->kcst());
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setTriSReacK(uint tidx, uint ridx, double kf)
{
	assert (tidx < pTris.size());
	assert (ridx < statedef()->countSReacs());

	if (pTris[tidx] == 0)
	{
    	std::ostringstream os;
    	os << "Triangle " << tidx << " has not been assigned to a patch.\n";
    	throw steps::ArgErr(os.str());
	}

	stex::Tri * tri = pTris[tidx];

	uint lsridx = tri->patchdef()->sreacG2L(ridx);
	if (lsridx == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "Surface reaction undefined in triangle.\n";
		throw steps::ArgErr(os.str());
	}

	tri->sreac(lsridx)->setKcst(kf);

	_updateElement(tri->sreac(lsridx));
    _updateSum();

}

////////////////////////////////////////////////////////////////////////////////

bool stex::Tetexact::_getTriSReacActive(uint tidx, uint ridx) const
{
	assert (tidx < pTris.size());
	assert (ridx < statedef()->countSReacs());

	if (pTris[tidx] == 0)
	{
    	std::ostringstream os;
    	os << "Triangle " << tidx << " has not been assigned to a patch.\n";
    	throw steps::ArgErr(os.str());
	}

	stex::Tri * tri = pTris[tidx];

	uint lsridx = tri->patchdef()->sreacG2L(ridx);
	if (lsridx == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "Surface reaction undefined in triangle.\n";
		throw steps::ArgErr(os.str());
	}

	if (tri->sreac(lsridx)->inactive() == true) return false;
	return true;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setTriSReacActive(uint tidx, uint ridx, bool act)
{
	assert (tidx < pTris.size());
	assert (ridx < statedef()->countSReacs());

	if (pTris[tidx] == 0)
	{
    	std::ostringstream os;
    	os << "Triangle " << tidx << " has not been assigned to a patch.\n";
    	throw steps::ArgErr(os.str());
	}

	stex::Tri * tri = pTris[tidx];

	uint lsridx = tri->patchdef()->sreacG2L(ridx);
	if (lsridx == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "Surface reaction undefined in triangle.\n";
		throw steps::ArgErr(os.str());
	}

	tri->sreac(lsridx)->setActive(act);

	_updateElement(tri->sreac(lsridx));
    _updateSum();
}

////////////////////////////////////////////////////////////////////////////////

bool stex::Tetexact::_getTriVDepSReacActive(uint tidx, uint vsridx) const
{
	assert (tidx < pTris.size());
	assert (vsridx < statedef()->countVDepSReacs());

	if (pTris[tidx] == 0)
	{
    	std::ostringstream os;
    	os << "Triangle " << tidx << " has not been assigned to a patch.\n";
    	throw steps::ArgErr(os.str());
	}

	stex::Tri * tri = pTris[tidx];

	uint lvsridx = tri->patchdef()->vdepsreacG2L(vsridx);
	if (lvsridx == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "Voltage-dependent surface reaction undefined in triangle.\n";
		throw steps::ArgErr(os.str());
	}

	if (tri->vdepsreac(lvsridx)->inactive() == true) return false;
	return true;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setTriVDepSReacActive(uint tidx, uint vsridx, bool act)
{
	assert (tidx < pTris.size());
	assert (vsridx < statedef()->countVDepSReacs());

	if (pTris[tidx] == 0)
	{
    	std::ostringstream os;
    	os << "Triangle " << tidx << " has not been assigned to a patch.\n";
    	throw steps::ArgErr(os.str());
	}

	stex::Tri * tri = pTris[tidx];

	uint lvsridx = tri->patchdef()->vdepsreacG2L(vsridx);
	if (lvsridx == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "Voltage-dependent surface reaction undefined in triangle.\n";
		throw steps::ArgErr(os.str());
	}

	tri->vdepsreac(lvsridx)->setActive(act);

	_updateElement(tri->vdepsreac(lvsridx));
    _updateSum();
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getTriSReacH(uint tidx, uint ridx) const
{
	assert (tidx < pTris.size());
	assert (ridx < statedef()->countSReacs());

	if (pTris[tidx] == 0)
	{
    	std::ostringstream os;
    	os << "Triangle " << tidx << " has not been assigned to a patch.\n";
    	throw steps::ArgErr(os.str());
	}

	stex::Tri * tri = pTris[tidx];

	uint lsridx = tri->patchdef()->sreacG2L(ridx);
	if (lsridx == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "Surface reaction undefined in triangle.\n";
		throw steps::ArgErr(os.str());
	}

	return tri->sreac(lsridx)->h();
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getTriSReacC(uint tidx, uint ridx) const
{
	assert (tidx < pTris.size());
	assert (ridx < statedef()->countSReacs());

	if (pTris[tidx] == 0)
	{
    	std::ostringstream os;
    	os << "Triangle " << tidx << " has not been assigned to a patch.\n";
    	throw steps::ArgErr(os.str());
	}

	stex::Tri * tri = pTris[tidx];

	uint lsridx = tri->patchdef()->sreacG2L(ridx);
	if (lsridx == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "Surface reaction undefined in triangle.\n";
		throw steps::ArgErr(os.str());
	}

	return tri->sreac(lsridx)->c();
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getTriSReacA(uint tidx, uint ridx) const
{
	assert (tidx < pTris.size());
	assert (ridx < statedef()->countSReacs());

	if (pTris[tidx] == 0)
	{
    	std::ostringstream os;
    	os << "Triangle " << tidx << " has not been assigned to a patch.\n";
    	throw steps::ArgErr(os.str());
	}

	stex::Tri * tri = pTris[tidx];

	uint lsridx = tri->patchdef()->sreacG2L(ridx);
	if (lsridx == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "Surface reaction undefined in triangle.\n";
		throw steps::ArgErr(os.str());
	}

	return tri->sreac(lsridx)->rate();
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::setEfieldDT(double efdt)
{
	if (efflag() != true)
	{
    	std::ostringstream os;
    	os << "Method not available: EField calculation not included in simulation.";
    	throw steps::ArgErr(os.str());
	}
	if (efdt <= 0.0)
	{
    	std::ostringstream os;
    	os << "EField dt must be graeter than zero.";
    	throw steps::ArgErr(os.str());
	}
	pEFDT = efdt;
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getTetV(uint tidx) const
{
	if (efflag() != true)
	{
    	std::ostringstream os;
    	os << "Method not available: EField calculation not included in simulation.";
    	throw steps::ArgErr(os.str());
	}
	int loctidx = pEFTet_GtoL[tidx];
	if (loctidx == -1)
	{
    	std::ostringstream os;
    	os << "Tetrahedron index " << tidx << " not assigned to a conduction volume.";
    	throw steps::ArgErr(os.str());
	}

    // EField object should convert value to base s.i. units
    return pEField->getTetV(loctidx);

}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setTetV(uint tidx, double v)
{
	if (efflag() != true)
	{
    	std::ostringstream os;
    	os << "Method not available: EField calculation not included in simulation.";
    	throw steps::ArgErr(os.str());
	}
	int loctidx = pEFTet_GtoL[tidx];
	if (loctidx == -1)
	{
    	std::ostringstream os;
    	os << "Tetrahedron index " << tidx << " not assigned to a conduction volume.";
    	throw steps::ArgErr(os.str());
	}

	// EField object should convert to millivolts
    pEField->setTetV(loctidx, v);
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getTetVClamped(uint tidx) const
{
	if (efflag() != true)
	{
    	std::ostringstream os;
    	os << "Method not available: EField calculation not included in simulation.";
    	throw steps::ArgErr(os.str());
	}
	int loctidx = pEFTet_GtoL[tidx];
	if (loctidx == -1)
	{
    	std::ostringstream os;
    	os << "Tetrahedron index " << tidx << " not assigned to a conduction volume.";
    	throw steps::ArgErr(os.str());
	}

    return pEField->getTetVClamped(loctidx);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setTetVClamped(uint tidx, bool cl)
{
	if (efflag() != true)
	{
    	std::ostringstream os;
    	os << "Method not available: EField calculation not included in simulation.";
    	throw steps::ArgErr(os.str());
	}
	int loctidx = pEFTet_GtoL[tidx];
	if (loctidx == -1)
	{
    	std::ostringstream os;
    	os << "Tetrahedron index " << tidx << " not assigned to a conduction volume.";
    	throw steps::ArgErr(os.str());
	}

    pEField->setTetVClamped(loctidx, cl);
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getTriV(uint tidx) const
{
	if (efflag() != true)
	{
    	std::ostringstream os;
    	os << "Method not available: EField calculation not included in simulation.";
    	throw steps::ArgErr(os.str());
	}
	int loctidx = pEFTri_GtoL[tidx];
	if (loctidx == -1)
	{
    	std::ostringstream os;
    	os << "Triangle index " << tidx << " not assigned to a membrane.";
    	throw steps::ArgErr(os.str());
	}
    // EField object should convert value to base s.i. units
    return pEField->getTriV(loctidx);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setTriV(uint tidx, double v)
{
	if (efflag() != true)
	{
    	std::ostringstream os;
    	os << "Method not available: EField calculation not included in simulation.";
    	throw steps::ArgErr(os.str());
	}
	int loctidx = pEFTri_GtoL[tidx];
	if (loctidx == -1)
	{
    	std::ostringstream os;
    	os << "Triangle index " << tidx << " not assigned to a membrane.";
    	throw steps::ArgErr(os.str());
	}

    // EField object should convert to millivolts
    pEField->setTriV(loctidx, v);
}

////////////////////////////////////////////////////////////////////////////////

bool stex::Tetexact::_getTriVClamped(uint tidx) const
{
	if (efflag() != true)
	{
    	std::ostringstream os;
    	os << "Method not available: EField calculation not included in simulation.";
    	throw steps::ArgErr(os.str());
	}
	int loctidx = pEFTri_GtoL[tidx];
	if (loctidx == -1)
	{
    	std::ostringstream os;
    	os << "Triangle index " << tidx << " not assigned to a membrane.";
    	throw steps::ArgErr(os.str());
	}

    return pEField->getTriVClamped(loctidx);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setTriVClamped(uint tidx, bool cl)
{
	if (efflag() != true)
	{
    	std::ostringstream os;
    	os << "Method not available: EField calculation not included in simulation.";
    	throw steps::ArgErr(os.str());
	}
	int loctidx = pEFTri_GtoL[tidx];
	if (loctidx == -1)
	{
    	std::ostringstream os;
    	os << "Triangle index " << tidx << " not assigned to a membrane.";
    	throw steps::ArgErr(os.str());
	}

    pEField->setTriVClamped(loctidx, cl);
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getTriOhmicI(uint tidx) const
{
	if (efflag() != true)
	{
    	std::ostringstream os;
    	os << "Method not available: EField calculation not included in simulation.";
    	throw steps::ArgErr(os.str());
	}

	stex::Tri * tri = pTris[tidx];

	int loctidx = pEFTri_GtoL[tidx];
	if (loctidx == -1)
	{
    	std::ostringstream os;
    	os << "Triangle index " << tidx << " not assigned to a membrane.";
    	throw steps::ArgErr(os.str());
	}

	return tri->getOhmicI(pEField->getTriV(loctidx), efdt());
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getTriOhmicI(uint tidx, uint ocidx) const
{
	assert (tidx < pTris.size());
	assert (ocidx < statedef()->countOhmicCurrs());

	int loctidx = pEFTri_GtoL[tidx];
	if (loctidx == -1)
	{
    	std::ostringstream os;
    	os << "Triangle index " << tidx << " not assigned to a membrane.";
    	throw steps::ArgErr(os.str());
	}

	stex::Tri * tri = pTris[tidx];

	uint locidx = tri->patchdef()->ohmiccurrG2L(ocidx);
	if (locidx == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "Ohmic current undefined in triangle.\n";
		throw steps::ArgErr(os.str());
	}

	return tri->getOhmicI(locidx, pEField->getTriV(loctidx), efdt());
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getTriGHKI(uint tidx) const
{
	if (efflag() != true)
	{
    	std::ostringstream os;
    	os << "Method not available: EField calculation not included in simulation.";
    	throw steps::ArgErr(os.str());
	}

	stex::Tri * tri = pTris[tidx];

	return tri->getGHKI(efdt());
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getTriGHKI(uint tidx, uint ghkidx) const
{
	if (efflag() != true)
	{
    	std::ostringstream os;
    	os << "Method not available: EField calculation not included in simulation.";
    	throw steps::ArgErr(os.str());
	}

	stex::Tri * tri = pTris[tidx];

	uint locidx = tri->patchdef()->ghkcurrG2L(ghkidx);
	if (locidx == ssolver::LIDX_UNDEFINED)
	{
		std::ostringstream os;
		os << "GHK current undefined in triangle.\n";
		throw steps::ArgErr(os.str());
	}

	return tri->getGHKI(locidx, efdt());
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getTriI(uint tidx) const
{
	if (efflag() != true)
	{
    	std::ostringstream os;
    	os << "Method not available: EField calculation not included in simulation.";
    	throw steps::ArgErr(os.str());
	}

	int loctidx = pEFTri_GtoL[tidx];
	if (loctidx == -1)
	{
    	std::ostringstream os;
    	os << "Triangle index " << tidx << " not assigned to a membrane.";
    	throw steps::ArgErr(os.str());
	}

    // EField object should convert to required units
    return pEField->getTriI(loctidx);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setVertI(uint vidx, double cur)
{
	if (efflag() != true)
	{
    	std::ostringstream os;
    	os << "Method not available: EField calculation not included in simulation.";
    	throw steps::ArgErr(os.str());
	}
	int locvidx = pEFVert_GtoL[vidx];
	if (locvidx == -1)
	{
    	std::ostringstream os;
    	os << "Vertex index " << vidx << " not assigned to a conduction volume or membrane.";
    	throw steps::ArgErr(os.str());
	}

    // EField object should convert to required units
    pEField->setVertI(locvidx, cur);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setVertIClamp(uint vidx, double cur)
{
	if (efflag() != true)
	{
    	std::ostringstream os;
    	os << "Method not available: EField calculation not included in simulation.";
    	throw steps::ArgErr(os.str());
	}
	int locvidx = pEFVert_GtoL[vidx];
	if (locvidx == -1)
	{
    	std::ostringstream os;
    	os << "Vertex index " << vidx << " not assigned to a conduction volume or membrane.";
    	throw steps::ArgErr(os.str());
	}

    // EField object should convert to required units
    pEField->setVertIClamp(locvidx, cur);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setTriI(uint tidx, double cur)
{
	if (efflag() != true)
	{
    	std::ostringstream os;
    	os << "Method not available: EField calculation not included in simulation.";
    	throw steps::ArgErr(os.str());
	}
	int loctidx = pEFTri_GtoL[tidx];
	if (loctidx == -1)
	{
    	std::ostringstream os;
    	os << "Triangle index " << tidx << " not assigned to a membrane.";
    	throw steps::ArgErr(os.str());
	}

    // EField object should convert to required units
    pEField->setTriI(loctidx, cur);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setTriIClamp(uint tidx, double cur)
{
	if (efflag() != true)
	{
    	std::ostringstream os;
    	os << "Method not available: EField calculation not included in simulation.";
    	throw steps::ArgErr(os.str());
	}
	int loctidx = pEFTri_GtoL[tidx];
	if (loctidx == -1)
	{
    	std::ostringstream os;
    	os << "Triangle index " << tidx << " not assigned to a membrane.";
    	throw steps::ArgErr(os.str());
	}

    // EField object should convert to required units
    pEField->setTriIClamp(loctidx, cur);
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getVertV(uint vidx) const
{
	if (efflag() != true)
	{
    	std::ostringstream os;
    	os << "Method not available: EField calculation not included in simulation.";
    	throw steps::ArgErr(os.str());
	}
	int locvidx = pEFVert_GtoL[vidx];
	if (locvidx == -1)
	{
    	std::ostringstream os;
    	os << "Vertex index " << vidx << " not assigned to a conduction volume or membrane.";
    	throw steps::ArgErr(os.str());
	}
    // EField object should convert value to base s.i. units
    return pEField->getVertV(locvidx);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setVertV(uint vidx, double v)
{
	if (efflag() != true)
	{
    	std::ostringstream os;
    	os << "Method not available: EField calculation not included in simulation.";
    	throw steps::ArgErr(os.str());
	}
	int locvidx = pEFVert_GtoL[vidx];
	if (locvidx == -1)
	{
    	std::ostringstream os;
    	os << "Vertex index " << vidx << " not assigned to a conduction volume or membrane.";
    	throw steps::ArgErr(os.str());
	}
    // EField object should convert to millivolts
    pEField->setVertV(locvidx, v);
}

////////////////////////////////////////////////////////////////////////////////

bool stex::Tetexact::_getVertVClamped(uint vidx) const
{
	if (efflag() != true)
	{
    	std::ostringstream os;
    	os << "Method not available: EField calculation not included in simulation.";
    	throw steps::ArgErr(os.str());
	}
	int locvidx = pEFVert_GtoL[vidx];
	if (locvidx == -1)
	{
    	std::ostringstream os;
    	os << "Vertex index " << vidx << " not assigned to a conduction volume or membrane.";
    	throw steps::ArgErr(os.str());
	}

    return pEField->getVertVClamped(locvidx);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setVertVClamped(uint vidx, bool cl)
{
	if (efflag() != true)
	{
    	std::ostringstream os;
    	os << "Method not available: EField calculation not included in simulation.";
    	throw steps::ArgErr(os.str());
	}
	int locvidx = pEFVert_GtoL[vidx];
	if (locvidx == -1)
	{
    	std::ostringstream os;
    	os << "Vertex index " << vidx << " not assigned to a conduction volume or membrane.";
    	throw steps::ArgErr(os.str());
	}
    // EField object should convert to millivolts
    pEField->setVertVClamped(locvidx, cl);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setMembRes(uint midx, double ro, double vrev)
{
	if (efflag() != true)
	{
    	std::ostringstream os;
    	os << "Method not available: EField calculation not included in simulation.";
    	throw steps::ArgErr(os.str());
	}
	if (ro <= 0.0)
	{
    	std::ostringstream os;
    	os << "Resistivity must be greater than zero.";
    	throw steps::ArgErr(os.str());
	}
    // EField object should convert to required units
	assert (midx == 0);
    pEField->setSurfaceResistivity(midx, ro, vrev);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setMembPotential(uint midx, double v)
{
	if (efflag() != true)
	{
    	std::ostringstream os;
    	os << "Method not available: EField calculation not included in simulation.";
    	throw steps::ArgErr(os.str());
	}
    // EField object should convert to millivolts
	assert (midx == 0);
    pEField->setMembPotential(midx, v);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setMembCapac(uint midx, double cm)
{
	if (efflag() != true)
	{
    	std::ostringstream os;
    	os << "Method not available: EField calculation not included in simulation.";
    	throw steps::ArgErr(os.str());
	}
	if (cm < 0.0)
	{
    	std::ostringstream os;
    	os << "Capacitance must be greater than or equal to zero.";
    	throw steps::ArgErr(os.str());
	}


    // EField object should convert to required units
	assert (midx == 0);
    pEField->setMembCapac(midx, cm);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setMembVolRes(uint midx, double ro)
{
	if (efflag() != true)
	{
    	std::ostringstream os;
    	os << "Method not available: EField calculation not included in simulation.";
    	throw steps::ArgErr(os.str());
	}
	if (ro < 0.0)
	{
    	std::ostringstream os;
    	os << "Resistivity must be greater than or equal to zero.";
    	throw steps::ArgErr(os.str());
	}
    // EField object should convert to required units
	assert (midx == 0);
    pEField->setMembVolRes(midx, ro);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_updateElement(KProc* kp)
{

    #ifdef SSA_DEBUG
    std::cout << "SSA: Update KProc element " << kp->schedIDX() << "\n";
    #endif

    double new_rate = kp->rate(this);

    #ifdef SSA_DEBUG
    std::cout << "new rate: " << new_rate << "\n";
    #endif

    CRKProcData & data = kp->crData;
    double old_rate = data.rate;

    data.rate = new_rate;

    #ifdef SSA_DEBUG
    std::cout << "data recorded: " << (data.recorded ? "Yes" : "No") << "\n";
    std::cout << "pow: " << data.pow << "\n";
    std::cout << "pos: " << data.pos << "\n";
    std::cout << "rate: " << data.rate << "\n";
    #endif

    if (old_rate == new_rate)  return;


    // new rate in positive groups
    if (new_rate >= 0.5) {
        #ifdef SSA_DEBUG
        std::cout << "new rate in positive group\n";
        #endif

        // pow is the same
        int old_pow = data.pow;
        int new_pow;
        double temp = frexp(new_rate, &new_pow);

        #ifdef SSA_DEBUG
        std::cout << "power changes from " << old_pow << " to " << new_pow << "\n";
        #endif

        if (old_pow == new_pow && data.recorded) {

            CRGroup* old_group = _getGroup(old_pow);

            #ifdef SSA_DEBUG
            std::cout << "in the same power group, change group sum and positive sum\n";
            std::cout << "old group sum: " << old_group->sum << " old pSum " << pSum << "\n";
            #endif

            old_group->sum += (new_rate - old_rate);

            #ifdef SSA_DEBUG
            std::cout << "new group sum: " << old_group->sum << " new pSum " << pSum << "\n";
            #endif
        }
        // pow is not the same
        else {
            data.pow = new_pow;


            if (data.recorded) {
                #ifdef SSA_DEBUG
                std::cout << "remove data from group with power " << old_pow << "\n";
                #endif

                // remove old
                CRGroup* old_group = _getGroup(old_pow);
                (old_group->size) --;

                if (old_group->size == 0) old_group->sum = 0.0;
                else {
                    old_group->sum -= old_rate;

                    KProc* last = old_group->indices[old_group->size];
                    old_group->indices[data.pos] = last;
                    last->crData.pos = data.pos;
                }
            }

            // add new
            #ifdef SSA_DEBUG
            std::cout << "add data to group with power " << new_pow << "\n";
            #endif
            if (pGroups.size() <= new_pow) _extendPGroups(new_pow + 1);

            CRGroup* new_group = pGroups[new_pow];

            assert(new_group != NULL);
            if (new_group->size == new_group->capacity) _extendGroup(new_group);
            uint pos = new_group->size;
            new_group->indices[pos] = kp;
            new_group->size++;
            new_group->sum += new_rate;
            data.pos = pos;

        }
        data.recorded = true;

    }
    // new rate in negative group
    else if (new_rate < 0.5 && new_rate > 1e-20) {
        int old_pow = data.pow;
        int new_pow;
        double temp = frexp(new_rate, &new_pow);

        #ifdef SSA_DEBUG
        std::cout << "power changes from " << old_pow << " to " << new_pow << "\n";
        #endif

        if (old_pow == new_pow && data.recorded) {

            CRGroup* old_group = _getGroup(old_pow);
            #ifdef SSA_DEBUG
            std::cout << "in the same power group, change group sum and negative sum\n";
            std::cout << "old group sum: " << old_group->sum << " old nSum " << nSum << "\n";
            #endif

            old_group->sum += (new_rate - old_rate);

            #ifdef SSA_DEBUG
            std::cout << "new group sum: " << old_group->sum << " new nSum " << nSum << "\n";
            #endif
        }
        // pow is not the same
        else {
            data.pow = new_pow;

            if (data.recorded) {
                #ifdef SSA_DEBUG
                std::cout << "remove data from group with power " << old_pow << "\n";
                #endif

                CRGroup* old_group = _getGroup(old_pow);
                (old_group->size) --;

                if (old_group->size == 0) old_group->sum = 0.0;
                else {
                    old_group->sum -= old_rate;

                    KProc* last = old_group->indices[old_group->size];
                    old_group->indices[data.pos] = last;
                    last->crData.pos = data.pos;
                }
            }

            // add new
            #ifdef SSA_DEBUG
            std::cout << "add data to group with power " << new_pow << "\n";
            #endif

            if (nGroups.size() <= -new_pow) _extendNGroups(-new_pow + 1);

            CRGroup* new_group = nGroups[-new_pow];

            if (new_group->size == new_group->capacity) _extendGroup(new_group);
            uint pos = new_group->size;
            new_group->indices[pos] = kp;
            new_group->size++;
            new_group->sum += new_rate;
            data.pos = pos;

        }
        data.recorded = true;
    }

    else {

        if (data.recorded) {
            #ifdef SSA_DEBUG
            std::cout << "remove data from group with power " << data.pow << "\n";
            #endif
            CRGroup* old_group = _getGroup(data.pow);

            // remove old
            old_group->size --;

            if (old_group->size == 0) old_group->sum = 0.0;
            else {
                old_group->sum -= old_rate;

                KProc* last = old_group->indices[old_group->size];
                old_group->indices[data.pos] = last;
                last->crData.pos = data.pos;
            }
        }
        data.recorded = false;
    }

    #ifdef SSA_DEBUG
    std::cout << "--------------------------------------------------------\n";
    #endif
}
////////////////////////////////////////////////////////////////////////////////

// END
