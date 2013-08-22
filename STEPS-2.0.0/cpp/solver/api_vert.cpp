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
#include <string>
#include <sstream>

// STEPS headers.
#include "../common.h"
#include "../error.hpp"
#include "api.hpp"
#include "statedef.hpp"
#include "../geom/tetmesh.hpp"

////////////////////////////////////////////////////////////////////////////////

USING(std, string);
USING_NAMESPACE(steps::solver);

////////////////////////////////////////////////////////////////////////////////

double API::getVertV(uint vidx) const
{
	if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
	{
		if (vidx >= mesh->countVertices())
		{
			std::ostringstream os;
			os << "Vertex index out of range.";
			throw steps::ArgErr(os.str());
		}
		return _getVertV(vidx);
	}

	else
	{
		std::ostringstream os;
		os << "Method not available for this solver.";
		throw steps::NotImplErr();
	}
}

////////////////////////////////////////////////////////////////////////////////

void API::setVertV(uint vidx, double v)
{
	if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
	{
		if (vidx >= mesh->countVertices())
		{
			std::ostringstream os;
			os << "Vertex index out of range.";
			throw steps::ArgErr(os.str());
		}
		_setVertV(vidx, v);
	}

	else
	{
		std::ostringstream os;
		os << "Method not available for this solver.";
		throw steps::NotImplErr();
	}
}

////////////////////////////////////////////////////////////////////////////////

bool API::getVertVClamped(uint vidx) const
{
	if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
	{
		if (vidx >= mesh->countVertices())
		{
			std::ostringstream os;
			os << "Vertex index out of range.";
			throw steps::ArgErr(os.str());
		}
		_getVertVClamped(vidx);
	}

	else
	{
		std::ostringstream os;
		os << "Method not available for this solver.";
		throw steps::NotImplErr();
	}
}

////////////////////////////////////////////////////////////////////////////////

void API::setVertVClamped(uint vidx, bool cl)
{
	if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
	{
		if (vidx >= mesh->countVertices())
		{
			std::ostringstream os;
			os << "Vertex index out of range.";
			throw steps::ArgErr(os.str());
		}
		_setVertVClamped(vidx, cl);
	}

	else
	{
		std::ostringstream os;
		os << "Method not available for this solver.";
		throw steps::NotImplErr();
	}
}

////////////////////////////////////////////////////////////////////////////////

void API::setVertIClamp(uint vidx, double i)
{
	if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
	{
		if (vidx >= mesh->countVertices())
		{
			std::ostringstream os;
			os << "Vertex index out of range.";
			throw steps::ArgErr(os.str());
		}
		_setVertIClamp(vidx, i);
	}

	else
	{
		std::ostringstream os;
		os << "Method not available for this solver.";
		throw steps::NotImplErr();
	}
}

////////////////////////////////////////////////////////////////////////////////

double API::_getVertV(uint vidx) const
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void API::_setVertV(uint vidx, double v)
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

bool API::_getVertVClamped(uint vidx) const
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void API::_setVertVClamped(uint vidx, bool cl)
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

void API::_setVertIClamp(uint vidx, double i)
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

// END
