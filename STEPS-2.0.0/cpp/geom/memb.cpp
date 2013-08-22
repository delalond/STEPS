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
#include <cassert>
#include <algorithm>
#include <vector>
#include <sstream>
#include <string>
#include <queue>
#include <set>

// STEPS headers.
#include "../common.h"
#include "../error.hpp"
#include "memb.hpp"
#include "tri.hpp"
#include "tmpatch.hpp"
#include "tmcomp.hpp"

NAMESPACE_ALIAS(steps::tetmesh, stetmesh);

////////////////////////////////////////////////////////////////////////////////

stetmesh::Memb::Memb(std::string const & id, Tetmesh * container,
					 std::vector<stetmesh::TmPatch *> const & patches,
					 bool verify, uint opt_method, std::string const & opt_file_name)
: pID(id)
, pTetmesh(container)
, pTri_indices()
, pTrivirt_indices()
, pTet_indices()
, pVert_indices()
, pTrisN(0)
, pTetsN(0)
, pTriVirtsN(0)
, pVertsN(0)
, pOpen(false)
, pOpt_method(opt_method)
, pOpt_file_name(opt_file_name)
{
    if (pTetmesh == 0)
    {
  	    std::ostringstream os;
  	    os << "No mesh provided to Membrane initializer function.\n";
  	    throw steps::ArgErr(os.str());
    }
    if (patches.size() == 0)
    {
  	    std::ostringstream os;
  	    os << "No Patches provided to Membrane initializer function.\n";
  	    throw steps::ArgErr(os.str());
    }

    if (pOpt_method != 1 and pOpt_method != 2)
    {
    	std::ostringstream os;
		os << "Unknown optimization method. Choices are 1 or 2.\n";
		throw steps::ArgErr(os.str());
    }

    // BUGFIX, IH: This should be a set because the patches can share an inner compartment,
    // which would otherwise result in multiple occurances of tets
    std::set<uint> tets_set = std::set<uint>();

    std::vector<TmPatch *>::const_iterator p_end = patches.end();
    for (std::vector<TmPatch *>::const_iterator p = patches.begin(); p != p_end; ++p)
    {
    	std::vector<uint> tris_temp = (*p)->getAllTriIndices();
    	pTri_indices.insert(pTri_indices.end(), tris_temp.begin(), tris_temp.end());
    	pTrisN += tris_temp.size();

    	if (stetmesh::TmComp * comp_temp = dynamic_cast<stetmesh::TmComp *>((*p)->getIComp()))
    	{
    		std::vector<uint> tets_temp = comp_temp->getAllTetIndices();
    		tets_set.insert(tets_temp.begin(), tets_temp.end());
			//pTet_indices.insert(pTet_indices.end(), tets_temp.begin(), tets_temp.end());
    		//pTetsN += tets_temp.size();
    	}
    	else
    	{
      	    std::ostringstream os;
      	    os << "Conduction volume (inner compartment(s) to membrane patch(es) must be ";
      	    os << "tetrahedral and not well-mixed.\n";
      	    throw steps::ArgErr(os.str());
    	}
    }

    pTetsN = tets_set.size();

    pTet_indices.resize(pTetsN);

    std::copy(tets_set.begin(), tets_set.end(), pTet_indices.begin());

    assert (pTrisN == pTri_indices.size());
    assert(pTetsN == pTet_indices.size());


    /* Previous code with user-given triangle indices
    // The maximum triangle index in tetrahedral mesh
    uint maxidx = (pTetmesh->countTris() -1);

    for (uint i=0; i <tris.size(); ++i)
    {
    	// perform some checks on this triangle
    	bool included = false;
    	for(uint j = 0; j <pTrisN; ++j)
    	{
    		// check if tri has already occured in this list (duplicate)
    		if (tris[i] == pTri_indices[j])
    		{
    			included = true;
    			break;
    		}
    	}
    	if (included == true) continue;

    	if (tris[i] > maxidx)
    	{
    		std::ostringstream os;
    		os << "Invalid index supplied for triangle #" << i << " in list.";
    		throw steps::ArgErr(os.str());
    	}

		// add the triangle to this patch
		pTri_indices.push_back(tris[i]);
		++pTrisN;
	}
	*/

    if (verify)
    {
    	std::ostringstream os;
    	os << "Verifying membrane. This can take some time...\n";
    	// Now perform some checks on the triangle and print warning if surface is open, or has more than 3 neighbours
    	std::vector<uint>::const_iterator tri_end = pTri_indices.end();
    	for (std::vector<uint>::const_iterator tri = pTri_indices.begin(); tri != tri_end; ++tri)
    	{
    		uint neighbours = 0;

    		// Fetch all triangle neighbours- we don't know how many neighbours there are
    		std::set<uint> tritrineighbs = pTetmesh->getTriTriNeighbs(*(tri));

    		std::vector<uint>::const_iterator t_end = pTri_indices.end();
    		for (std::vector<uint>::const_iterator t = pTri_indices.begin(); t != t_end; ++t)
    		{
    			// continue of this is the same triangle
    			if ((*t) == (*tri)) continue;
    			if (tritrineighbs.find(*t) != tritrineighbs.end()) neighbours +=1;
    		}

    		// If 3 neighbours weren't found this is not a closed surface
    		// MOST SURFACE TRIANGLES WILL HAVE 3 NEIGHBOURS, BUT
    		// IF A TRIANGLE IS ON A SHARP EDGE IT CAN HAVE 5 OR EVEN 7 NEIGHBOURS
    		if ((!open()) && neighbours < 3)
    		{
    			os << "WARNING: Triangular surface provided to Membrane initializer function";
    			os << " is thought to be an open surface.\n";
    			pOpen = true;
    			// throw steps::ArgErr(os.str());
    		}
    		// What can we say if tri has more than 3 neighbours? Probably a bad surface, though with particularly
    		// 'pointed' regions neighbours can be 5, 7, or even higher even number and still
    		// be a closed surface. If one of these triangles further happens to be on an edge
    		// may have 4, 6, or higher even number, though these numbers would usually
    		// suggest a bad surface. However, don't complain yet- perform the multiple surface
    		// test later to catch that.
    		if (neighbours > 3)
    		{
    			os << "WARNING: Triangle " << (*tri) << " has " << neighbours << " neighbours.\n";
    		}

    		/*
			// If neighbours are even, this is a bad mesh.
			else if (neighbours > 3 && neighbours%2 == 0 )
			{
				std::ostringstream os;
				os << "Triangular surface provided to Membrane initializer function";
				os << " is thought to be an open surface.";
				throw steps::ArgErr(os.str());
			}
    		 */
    	}


    	// Now to perform a basic multiple surfaces test

    	// Create the queue of triangles- each triangle's neighbours will be
    	// added to the queue only once
    	std::queue<uint> tri_queue;
		// Let's start with the zeroth element
		tri_queue.push(pTri_indices[0]);

		// Create the set of triangle neighbours. If this set contains all triangle
		// indices at the end of the loop, this is a single closed surface
		std::set<uint> tri_set;

		while (! tri_queue.empty())
		{
			uint triidx = tri_queue.front();
			std::set<uint> trineighbs = pTetmesh->getTriTriNeighbs(triidx);
			std::vector<uint>::const_iterator t_end = pTri_indices.end();
			for (std::vector<uint>::const_iterator t = pTri_indices.begin(); t != t_end; ++t)
			{
				// continue if this is the same triangle
				if ((*t) == triidx) continue;
				if (trineighbs.find(*t) != trineighbs.end())
				{
					// A triangle's neighbour. Now need to check if it is already in the queue
					if (tri_set.find((*t)) == tri_set.end())
					{
						tri_queue.push(*t);
						tri_set.insert(*t);
					}
				}
			}
			tri_queue.pop();
		}

		if (tri_set.size() != pTri_indices.size())
		{
			std::ostringstream os;
			os << "Triangular surface provided to Membrane initializer function";
			os << " is a multiple surface.\n";
			throw steps::ArgErr(os.str());
		}


		/* Previous code when conduction volume was explicitly given by the user
		// Now perform checks on the conduction volume.
		// The maximum tetrahedron index in tetrahedral mesh
		uint maxtetidx = (pTetmesh->countTets() -1);

		for (uint i=0; i < vol_tets.size(); ++i)
		{
			// perform some checks on this triangle
			bool included = false;
			for(uint j = 0; j <pTetsN; ++j)
			{
				// check if tet has already occured in this list (duplicate)
				if (vol_tets[i] == pTet_indices[j])
				{
					included = true;
					break;
				}
			}
			if (included == true) continue;

			if (vol_tets[i] > maxtetidx)
			{
				std::ostringstream os;
				os << "Invalid index supplied for tetrahedron #" << i << " in list.";
				throw steps::ArgErr(os.str());
			}

			// add the triangle to this patch
			pTet_indices.push_back(vol_tets[i]);
			++pTetsN;
		}
		 */


		// Loop over all triangles and test that all triangles have one tetrahedron
		// neighbour in the volume.
		tri_end = pTri_indices.end();
		for (std::vector<uint>::const_iterator tri = pTri_indices.begin(); tri != tri_end; ++tri)
		{
			uint tetneighbours = 0;
			// Fetch the two triangle tet neighbours
			int * tritetneighbs = pTetmesh->_getTriTetNeighb(*(tri));

			std::vector<uint>::const_iterator t_end = pTet_indices.end();
			for (std::vector<uint>::const_iterator t = pTet_indices.begin(); t!= t_end; ++t)
			{
				if (tritetneighbs[0] == (*t) || tritetneighbs[1] == (*t)) tetneighbours+=1;
			}

			if (tetneighbours < 1)
			{
				std::ostringstream os;
				os << "Conduction volume provided to Membrane initializer function";
				os << " is not connected to membrane triangle # " << (*tri) << "\n";
				throw steps::ArgErr(os.str());
			}
			if (tetneighbours > 1)
			{
				std::ostringstream os;
				os << "Conduction volume provided to Membrane initializer function";
				os << " is connected twice to membrane triangle # " << (*tri) << "\n";
				throw steps::ArgErr(os.str());
			}
		}

		// Now perform the 'multiple volume' test

		// Create the queue of tetrahedrons- each tetrahedrons' neighbours
		// will be added to the queue only once
		std::queue<uint> tet_queue;
		// Start with the zeroth element
		tet_queue.push(pTet_indices[0]);

		// Create the set of tetrahedron's neighbours. If this set contains all
		// tetrahedrons at the end of the loop, this is a single volume.
		std::set<uint> tet_set;

		while(! tet_queue.empty())
		{
			uint tetidx = tet_queue.front();
			int * tetneighbs = pTetmesh->_getTetTetNeighb(tetidx);
			std::vector<uint>::const_iterator t_end = pTet_indices.end();
			for (std::vector<uint>::const_iterator t = pTet_indices.begin(); t != t_end; ++t)
			{
				// continue if this is the same tet
				if ((*t) == tetidx) continue;
				for (uint i = 0; i < 4; ++i)
				{
					if (tetneighbs[i] == (*t))
					{
						if (tet_set.find((*t)) == tet_set.end())
						{
							tet_queue.push(*t);
							tet_set.insert(*t);
						}

						// No need to test other neighbours
						break;
					}
				}
			}
			tet_queue.pop();
		}

		if (tet_set.size() != pTet_indices.size())
		{
			std::ostringstream os;
			os << "Conduction volume provided to Membrane initializer function";
			os << " is a multiple volume.\n";
			throw steps::ArgErr(os.str());
		}


		// Now to set up the 'virtual membrane triangles'. That is, with an open
		// surface, we need to know which triangles are in the open region.
		if (open())
		{
			std::vector<uint> trivirttemp = std::vector<uint>();
			// Loop over tetrahedrons and add triangles from surface tets that
			// are not already in pTri_indices
			std::vector<uint>::const_iterator tet_end = pTet_indices.end();
			for (std::vector<uint>::const_iterator tet = pTet_indices.begin(); tet != tet_end; ++tet)
			{
				// Fetch neighbours
				std::vector<int> tettetneighbs = pTetmesh->getTetTetNeighb(*tet);
				// Loop over all tets again and find how many neighbours this tet has.
				// If it has less than 4 it is a candidate for a surface tet whose
				// triangle (in the correct direction) makes up part of the
				// 'virtual surface'
				std::vector<uint>::const_iterator tettemp_end = pTet_indices.end();
				bool gotneighb[4] = {false, false, false, false};
				for (std::vector<uint>::const_iterator tettemp = pTet_indices.begin(); tettemp!= tettemp_end; ++tettemp)
				{
					if ((*tettemp) == (*tet)) continue;
					for (uint i = 0; i < 4; ++i)
					{
						if (tettetneighbs[i] == (*tettemp))
						{
							assert(gotneighb[i] == false);
							gotneighb[i] = true;
							break;
						}
					}
				}
				// Now need to check if the triangle in this direction is in
				// pTri_indices
				std::vector<uint> tettrineighb = pTetmesh->getTetTriNeighb(*tet);
				for (uint j = 0; j < 4; ++j)
				{
					if (gotneighb[j] == false)
					{
						// Now need to check if the triangle in this direction is in
						// pTri_indices
						tri_end = pTri_indices.end();
						bool membtri = false;
						for (std::vector<uint>::const_iterator tri = pTri_indices.begin(); tri != tri_end; ++tri)
						{
							if (tettrineighb[j] == *tri)
							{
								membtri = true;
								break;
							}
						}
						if (membtri == false) trivirttemp.push_back(pTetmesh->getTetTriNeighb(*tet)[j]);
					}
				}
			}

			// Ok, let's now go on a little walk and add the correct surface triangles
			// Create the queue of triangles- each triangle's neighbours will be
			// added to the queue only once
			std::queue<uint> tri_queue;
			// Let's start with the zeroth element
			tri_queue.push(pTri_indices[0]);

			// Create the set of triangle neighbours. This set will contain all the
			// surface triangles at the end- real and virtual
			std::set<uint> tri_set;

			// The virtual triangles, in a set
			std::set<uint> trivirt_set;

			while (! tri_queue.empty())
			{
				uint triidx = tri_queue.front();
				std::set<uint> trineighbs = pTetmesh->getTriTriNeighbs(triidx);

				std::vector<uint>::const_iterator t_end = pTri_indices.end();
				for (std::vector<uint>::const_iterator t = pTri_indices.begin(); t != t_end; ++t)
				{
					// continue if this is the same triangle
					if ((*t) == triidx) continue;
					if (trineighbs.find(*t) != trineighbs.end())
					{
						// A triangle's neighbour. Now need to check if it is already in the queue
						if (tri_set.find((*t)) == tri_set.end())
						{
							tri_queue.push(*t);
							tri_set.insert(*t);
						}
					}
				}
				std::vector<uint>::const_iterator tv_end = trivirttemp.end();
				for (std::vector<uint>::const_iterator tv = trivirttemp.begin(); tv!= tv_end; ++tv)
				{
					if ((*tv) == triidx) continue;
					if (trineighbs.find(*tv) != trineighbs.end())
					{
						if (tri_set.find((*tv)) == tri_set.end())
						{
							tri_queue.push(*tv);
							tri_set.insert(*tv);
							trivirt_set.insert(*tv);
						}
					}
				}
				tri_queue.pop();
			}

			// So trivirt_set holds all the virtual triangles. Just need to copy to
			// vector
			std::set<uint>::const_iterator tvset_end = trivirt_set.end();
			for (std::set<uint>::const_iterator tvset = trivirt_set.begin(); tvset != tvset_end; ++tvset)
			{
				pTrivirt_indices.push_back(*tvset);
			}
			pTriVirtsN = pTrivirt_indices.size();
		}

	}
	// Create empty set of vertex indices
	std::set<uint> vertset = std::set<uint>();

	std::vector<uint>::const_iterator tet_end = pTet_indices.end();
    for (std::vector<uint>::const_iterator tet = pTet_indices.begin(); tet != tet_end; ++tet)
    {
    	uint * tettemp = getContainer()->_getTet(*tet);
    	for (uint i = 0; i < 4; ++i)
    	{
    		vertset.insert(tettemp[i]);
    	}
    }

    std::set<uint>::const_iterator vset_end = vertset.end();
    for (std::set<uint>::const_iterator vset = vertset.begin(); vset != vset_end; ++vset)
    {
    	pVert_indices.push_back(*vset);
    }
    pVertsN = pVert_indices.size();

    pTetmesh->_handleMembAdd(this);
}

////////////////////////////////////////////////////////////////////////////////

stetmesh::Memb::~Memb(void)
{
	if (pTetmesh == 0) return;
	_handleSelfDelete();
}

////////////////////////////////////////////////////////////////////////////////

void stetmesh::Memb::_handleSelfDelete(void)
{
	pTetmesh->_handleMembDel(this);

	pTri_indices.clear();
	pTet_indices.clear();
	pTrivirt_indices.clear();
	pVert_indices.clear();

	pTrisN = 0;
	pTetsN = 0;
	pTriVirtsN = 0;
	pVertsN = 0;

	pTetmesh = 0;

}

////////////////////////////////////////////////////////////////////////////////

std::vector<bool> stetmesh::Memb::isTriInside(std::vector<uint> tri) const
{
	uint notris = tri.size();
	std::vector<bool> inside(notris);
	for (uint i=0; i < notris; ++i)
	{
		bool triinside = false;
		for (uint j=0; j< pTri_indices.size() ; ++j)
		{
			if (tri[i] == pTri_indices[j])
			{
				triinside = true;
				break;
			}
		}
		if (triinside == true) inside[i] = true;
		else inside[i] = false;
	}
	return inside;
}

////////////////////////////////////////////////////////////////////////////////

// END
