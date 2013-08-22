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

// STL headers.
#include <cassert>
#include <cmath>
#include <sstream>


// STEPS headers.
#include "../common.h"
#include "../math/triangle.hpp"
#include "tri.hpp"
#include "../error.hpp"

NAMESPACE_ALIAS(steps::tetmesh, stetmesh);

////////////////////////////////////////////////////////////////////////////////

stetmesh::Tri::Tri(Tetmesh * mesh, uint tidx)
: pTetmesh(mesh)
, pTidx(tidx)
, pVerts()
, pBaryc()
{
	if (pTetmesh == 0)
    {
		std::ostringstream os;
        os << "No mesh provided to Tri initializer function.\n";
        throw steps::ArgErr(os.str());
    }

    uint * tri_temp = pTetmesh->_getTri(tidx);
    pVerts[0] = tri_temp[0];
    pVerts[1] = tri_temp[1];
    pVerts[2] = tri_temp[2];


    pBaryc = new double[3];

}

////////////////////////////////////////////////////////////////////////////////

stetmesh::Tri::~Tri(void)
{
	delete pBaryc;
}

////////////////////////////////////////////////////////////////////////////////

double stetmesh::Tri::getArea(void) const
{
	assert(pTetmesh != 0);
	return (pTetmesh->getTriArea(pTidx));
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> stetmesh::Tri::getBarycenter(void) const
{
	double * v0 = pTetmesh->_getVertex(pVerts[0]);
	double * v1 = pTetmesh->_getVertex(pVerts[1]);
	double * v2 = pTetmesh->_getVertex(pVerts[2]);
	/*double v0[3], v1[3], v2[3];	// Defunct code. Pointers directly fetched
	for (uint i=0; i < 3; ++i)
	{
		v0[i] = v0vec[i];
		v1[i] = v1vec[i];
		v2[i] = v2vec[i];
	}*/
	double baryc[3];
	steps::math::triBarycenter(v0, v1, v2, baryc);
	std::vector<double> barycentre(3);
	barycentre[0] = baryc[0];
	barycentre[1] = baryc[1];
	barycentre[2] = baryc[2];
	return barycentre;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> stetmesh::Tri::getNorm(void) const
{
	assert (pTetmesh != 0);
	return (pTetmesh->getTriNorm(pTidx));
}

////////////////////////////////////////////////////////////////////////////////

stetmesh::TmPatch * stetmesh::Tri::getPatch(void) const
{
	assert (pTetmesh != 0);
	return (pTetmesh->getTriPatch(pTidx));
}

////////////////////////////////////////////////////////////////////////////////

stetmesh::Tet stetmesh::Tri::getTet(uint i) const
{
	assert(i <= 1);
	int tetidx = pTetmesh->_getTriTetNeighb(pTidx)[i];
	assert(tetidx != -1);
	return (Tet(pTetmesh, tetidx));
}

////////////////////////////////////////////////////////////////////////////////

int stetmesh::Tri::getTetIdx(uint i) const
{
	assert(i <= 1);
	return (pTetmesh->_getTriTetNeighb(pTidx)[i]);
}

////////////////////////////////////////////////////////////////////////////////

uint stetmesh::Tri::getVertexIdx(uint i) const
{
	assert (i <= 2);
	return pVerts[i];
}

////////////////////////////////////////////////////////////////////////////////

stetmesh::Tet stetmesh::Tri::getTet0(void) const
{
	return getTet(0);
}

////////////////////////////////////////////////////////////////////////////////

stetmesh::Tet stetmesh::Tri::getTet1(void) const
{
	return getTet(1);
}

////////////////////////////////////////////////////////////////////////////////

stetmesh::Tet stetmesh::Tri::getInnerTet(void) const
{
	return getTet(0);
}

////////////////////////////////////////////////////////////////////////////////

stetmesh::Tet stetmesh::Tri::getOuterTet(void) const
{
	return getTet(1);
}

////////////////////////////////////////////////////////////////////////////////

double * stetmesh::Tri::_getNorm(void) const
{
	assert (pTetmesh != 0);
	return (pTetmesh->_getTriNorm(pTidx));
}

////////////////////////////////////////////////////////////////////////////////

double * stetmesh::Tri::_getBarycenter(void) const
{
	double * v0 = pTetmesh->_getVertex(pVerts[0]);
	double * v1 = pTetmesh->_getVertex(pVerts[1]);
	double * v2 = pTetmesh->_getVertex(pVerts[2]);

	steps::math::triBarycenter(v0, v1, v2, pBaryc);

	return pBaryc;
}

////////////////////////////////////////////////////////////////////////////////

double stetmesh::Tri::getTriDist(uint i, int triidx) const
{
	assert(i <= 3);

	if (triidx != -1)
	{
		Tri * tritemp = new Tri(pTetmesh, triidx);
		assert(tritemp->getPatch() == getPatch());
		double * bary1 = _getBarycenter();
		double * bary2 = tritemp->_getBarycenter();
		double xdist = bary1[0] - bary2[0];
		double ydist = bary1[1] - bary2[1];
		double zdist = bary1[2] - bary2[2];
		delete tritemp;
		return (sqrt((xdist*xdist) + (ydist*ydist) + (zdist*zdist)));
	}
	else return 0.0;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<int> stetmesh::Tri::getTriIdxs(stetmesh::TmPatch* tmpatch ) const
{
	assert(getPatch() != 0);
	assert(tmpatch != 0);

	// Tetmesh::getTriTriNeighb conveniently now puts all the neighbouring triangles in the correct order
	return  pTetmesh->getTriTriNeighb(pTidx, tmpatch);



	/*
	std::vector<int> triidxs = std::vector<int>(3);
	triidxs[0]=-1;
	triidxs[1]=-1;
	triidxs[2]=-1;

    std::set<uint> tritrineighbs = pTetmesh->getTriTriNeighb(pTidx, tmpatch);

    std::vector<uint> tribars = pTetmesh->getTriBars(getIdx());

    uint neighbcount=0;
    std::set<uint>::const_iterator tn_end = tritrineighbs.end();
    for (std::set<uint>::const_iterator tn = tritrineighbs.begin(); tn!=tn_end; ++tn)
    {
    	stetmesh::Tri * tritemp = new stetmesh::Tri(pTetmesh, *tn);
    	if (tritemp->getPatch() == getPatch())
    	{
    		if (neighbcount > 2)
    		{
    			std::ostringstream os;
    	        os << "Error in Patch initialisation for '" << getPatch()->getID();
    	        os <<"'. Patch triangle idx " << getIdx() << " found to have more than 3 neighbours. \n";
    	        throw steps::ArgErr(os.str());
    		}
    		else
    		{
				// Find the right direction by comparing bars
    			std::vector<uint>  neighbtribars = pTetmesh->getTriBars(*tn);

    			bool found = false;
    			for (uint i = 0; i < 3; ++i)
    			{
    				if (found == true) continue;
    	    		for (uint j=0; j<3; ++j)
    				{
    	    			if (neighbtribars[j] == tribars[i])
    	    			{
    	    				assert(triidxs[i] == -1);
    	    				triidxs[i] = *tn;
    	    				found = true;
    	    				continue;
    	    			}
    				}
    			}
    			assert(found == true);

    			neighbcount++;
    		}
    	}
    }

	return triidxs;
	*/
}

////////////////////////////////////////////////////////////////////////////////

/*
int stetmesh::Tri::getTriIdx(uint i) const
{
	assert(i < 3);

	return getTriIdxs()[i];
}
*/

////////////////////////////////////////////////////////////////////////////////

double stetmesh::Tri::getBarLength(uint i) const
{
	assert(i <= 2);
	uint baridx = pTetmesh->getTriBars(getIdx())[i];

	std::vector<uint> barverts =  pTetmesh->getBar(baridx);

	std::vector<double> vert0 = pTetmesh->getVertex(barverts[0]);
	std::vector<double> vert1 = pTetmesh->getVertex(barverts[1]);

	double xdist  = vert0[0]-vert1[0];
	double ydist  = vert0[1]-vert1[1];
	double zdist  = vert0[2]-vert1[2];

	return (sqrt((xdist*xdist) + (ydist*ydist) + (zdist*zdist)));
}

////////////////////////////////////////////////////////////////////////////////

// END
