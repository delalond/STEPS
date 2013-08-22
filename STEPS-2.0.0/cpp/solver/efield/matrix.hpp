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

#ifndef STEPS_SOLVER_EFIELD_MATRIX_HPP
#define STEPS_SOLVER_EFIELD_MATRIX_HPP 1

#include <fstream>
#include <iostream>

// STEPS headers.
#include "../../common.h"

////////////////////////////////////////////////////////////////////////////////

START_NAMESPACE(steps)
START_NAMESPACE(solver)
START_NAMESPACE(efield)

////////////////////////////////////////////////////////////////////////////////

/// \todo Clean up (especially get rid of the pointer-to-pointer storage);
/// could be a useful addition for Boost STEPS.
///
/// \author Robert Cannon
///
class Matrix
{

public:

    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    /// Constructor that creates an uninitialized n0 * n0 matrix.
    ///
	Matrix(uint n0);

	/// Constructor that creates an nn * nn matrix and initializes it
	/// by copying the contents of da.
	///
	Matrix(uint nn, double ** da);

	/// Destructor.
	///
	~Matrix(void);

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream & cp_file);

    /// restore data
    void restore(std::fstream & cp_file);

	////////////////////////////////////////////////////////////////////////
	// MATRIX OPERATIONS
	////////////////////////////////////////////////////////////////////////

	/// Makes a deep copy of the matrix.
	///
	Matrix * copy(void);

	/// Computes left-hand vector product.
	///
	/// The resulting array needs to deallocated by the caller!
	///
	double * lvprod(double * v);

	/// Returns the transpose of this matrix.
	///
	Matrix * transpose(void);

	/// Computes the determinant of this matrix.
	///
	double det(void);

	/// Returns the inverse of this matrix.
	///
	Matrix * inverse(void);

	/// Compute the LU decomposition.
	///
	void LU(void);

	///
	double * lubksb(double*);

	/// NOT IMPLEMENTED????
	double * rvprod(double*);

	////////////////////////////////////////////////////////////////////////

private:

    double       ** pA;
    double        * pWS;
    uint            pN;
    int           * pPerm;
    int             pSign;

};

////////////////////////////////////////////////////////////////////////////////

END_NAMESPACE(efield)
END_NAMESPACE(solver)
END_NAMESPACE(steps)

////////////////////////////////////////////////////////////////////////////////

#endif
// STEPS_SOLVER_EFIELD_MATRIX_HPP

// END
