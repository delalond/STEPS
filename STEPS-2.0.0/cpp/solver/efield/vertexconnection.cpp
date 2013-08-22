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
#include <iostream>
#include <sstream>
#include <cassert>

// STEPS headers.
#include "../../common.h"
#include "vertexconnection.hpp"
#include "vertexelement.hpp"

////////////////////////////////////////////////////////////////////////////////

NAMESPACE_ALIAS(steps::solver::efield, sefield);
USING_NAMESPACE(std);

////////////////////////////////////////////////////////////////////////////////

sefield::VertexConnection::VertexConnection(sefield::VertexElement * v1, sefield::VertexElement * v2)
: pVert1(v1)
, pVert2(v2)
, pGeomCC(0.0)
{
    assert(v1 != 0);
    assert(v2 != 0);
    pVert1->addConnection(this);
    pVert2->addConnection(this);
}

////////////////////////////////////////////////////////////////////////////////

sefield::VertexConnection::~VertexConnection(void)
{
}

////////////////////////////////////////////////////////////////////////////////

void sefield::VertexConnection::checkpoint(std::fstream & cp_file)
{
    cp_file.write((char*)&pGeomCC, sizeof(double));
}

////////////////////////////////////////////////////////////////////////////////

void sefield::VertexConnection::restore(std::fstream & cp_file)
{
    cp_file.read((char*)&pGeomCC, sizeof(double));
}

////////////////////////////////////////////////////////////////////////////////

sefield::VertexElement * sefield::VertexConnection::getOther(sefield::VertexElement * element)
{
    VertexElement * ret;
    if (pVert1 == element)
    {
        ret = pVert2;
    }
    else if (pVert2 == element)
    {
        ret = pVert1;
    }
    else
    {
        assert(0);
        ret = 0;
    }
    return ret;
}

////////////////////////////////////////////////////////////////////////////////

//bool VertexConnection::hasInternalEnd(void)
//{
//    return (vea->isInternal() || veb->isInternal());
//}

////////////////////////////////////////////////////////////////////////////////

//bool VertexConnection::isEdge(void)
//{
//    return (vea->isEdge() && veb->isEdge());
//}

////////////////////////////////////////////////////////////////////////////////

// END
