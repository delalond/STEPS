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

#ifndef STEPS_SOLVER_EFIELD_BANDDIAGONALMATRIX_HPP
#define STEPS_SOLVER_EFIELD_BANDDIAGONALMATRIX_HPP 1

// STL headers.
#include <string>
#include <cmath>
#include <fstream>
#include <iostream>

// STEPS headers.
#include "../../common.h"

////////////////////////////////////////////////////////////////////////////////

START_NAMESPACE(steps)
START_NAMESPACE(solver)
START_NAMESPACE(efield)

////////////////////////////////////////////////////////////////////////////////

/// Reference: Numerical recipes, 2nd ed, section 2.4 (p. 50).
///
/// Author: Robert Cannon
///
class BandDiagonalMatrix
{

public:

    /// Constructor.
    ///
    /*
    BandDiagonalMatrix(int, int, double**, double**);
     */
    BandDiagonalMatrix(int, int, double*, double*);


    /// Destructor.
    ///
	~BandDiagonalMatrix(void);

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream & cp_file);

    /// restore data
    void restore(std::fstream & cp_file);

	/// Perform LU decomposition.
	///
	void lu(void);


	void lubksb(double*, double*);

private:

    /*
    double** a;
     */
    double* a;

    /*
    double** al;
     */
    double* al;


    double* ws;
    int* perm;
    int n;
    int halfbw;

};

////////////////////////////////////////////////////////////////////////////////

END_NAMESPACE(efield)
END_NAMESPACE(solver)
END_NAMESPACE(steps)

////////////////////////////////////////////////////////////////////////////////

#endif
// STEPS_SOLVER_EFIELD_BANDDIAGONALMATRIX_HPP

// END
