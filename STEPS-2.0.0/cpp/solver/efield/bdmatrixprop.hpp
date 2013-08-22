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

#ifndef STEPS_SOLVER_EFIELD_BDMATRIXPROP_HPP
#define STEPS_SOLVER_EFIELD_BDMATRIXPROP_HPP 1

// STL headers.
#include <ctime>
#include <fstream>
#include <iostream>

// STEPS headers.
#include "../../common.h"
#include "bdmatrix.hpp"
#include "tetmesh.hpp"

////////////////////////////////////////////////////////////////////////////////

START_NAMESPACE(steps)
START_NAMESPACE(solver)
START_NAMESPACE(efield)

////////////////////////////////////////////////////////////////////////////////

/// The main current problem with this class is that it has some implicit
/// assumptions about what happens at each time step. Important task:
/// document or otherwise 'fix' them to be more in line with STEPS.
///
class BandedMatrixProp
{

public:

    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    /// Constructor; creates a banded diagonal matrix for a certain mesh.
    ///
	BandedMatrixProp(TetMesh * msh);

	/// Destructor.
	///
	~BandedMatrixProp(void);

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream & cp_file);

    /// restore data
    void restore(std::fstream & cp_file);

    /// Set the initial potential (in millivolt).
	///
    //void setInitialPotential(double v);

    void setSurfaceConductance(double, double);

    std::string makeOutputLine(double, double*);

    //void reset(void);

    // New method to replace setInitialPotential and reset()
    void setPotential(double v);

	////////////////////////////////////////////////////////////////////////
	// METHODS
	////////////////////////////////////////////////////////////////////////

    /// This is the crucial method -- advance the system with a certain
	/// timestep dt.
	///
    void advance(double dt);	// converted to ms in Efield object

    ////////////////////////////////////////////////////////////////////////
    // METHODS: OBJECT ACCESS
    ////////////////////////////////////////////////////////////////////////

    /// Return the potential at vertex i in mV
    ///
    double getV(int i) const
    { return pV[i]; }

    /// Set the potential at vertex i in mV.
    ///
	void setV(int i, double d)
	{ pV[i] = d;}

	////////////////////////////////////////////////////////////////////////

	/// Return whether vertex i is clamped.
	///
	bool getClamped(int i) const
	{ return pVertexClamp[i]; }

	/// (De-)activate the voltage clamp on some vertex i.
	///
	void setClamped(int i, bool b)
    { pVertexClamp[i] = b; }

	////////////////////////////////////////////////////////////////////////

	/// Return the amount of current being injected in a vertex in picoamp
	///
	bool getInj(int i) const
    { return -pVertexInj[i]; }

	/// Set the amount of current being injected in a vertex in picoamp
	/// NOTE: Have to change the sign here to match the terms in the
	/// matrix solution
	void setInj(int i, double d)
	{ pVertexInj[i] = -d; }

	/// Set the amount of current that is constantly injected through a
	/// vertex element until explicitly cancelled or changed.
	/// NOTE: Have to change the sign here to match the terms in the
	/// matrix solution
    void setVertIClamp(int i, double c)
	{
		pVertCurClamp[i] = -c;
	}

	////////////////////////////////////////////////////////////////////////

	/// Returns the current that will be injected in a triangle in picoamp
	///
    double getTriI(int i) const
    { return -pTriCur[i]; }

    /// Set an amount of current (expressed in ...) which is 'injected'
    /// through a triangular surface element over one EField dt. This current is divided
    /// over the 3 vertices that comprise the triangle and added to the
    /// current injected in vertex, on top of the current set with
    /// BandedMatrixProp::setInj. Usually these triangle currents
    /// result from open channels and receptors computed in the
    /// biochemical part of STEPS.
	/// NOTE: Have to change the sign here to match the terms in the
	/// matrix solution
    void setTriI(int i, double d)
	{
		pTriCur[i] = -d;
	}

	/// Set the amount of current that is constantly injected through a
	/// triangular element until explicitly cancelled or changed.
	/// NOTE: Have to change the sign here to match the terms in the
	/// matrix solution
    void setTriIClamp(int i, double c)
	{
		pTriCurClamp[i] = -c;
	}

	////////////////////////////////////////////////////////////////////////

private:

    ////////////////////////////////////////////////////////////////////////
    // INTERNAL METHODS
    ////////////////////////////////////////////////////////////////////////

    /// Construct the matrix.
    ///
    void populateMRHS(double);

    /// What the hell does this do???? And when?? (See also: pTotGSurf,
    /// pInjTot in this class.)
    ///
    void chargeCheck(double);

    ////////////////////////////////////////////////////////////////////////
    // DATA FIELDS
    ////////////////////////////////////////////////////////////////////////

    /// Pointer to the mesh.
    ///
    TetMesh *                   pMesh;

    /// Number of vertices in the mesh, stored locally.
    ///
    uint                        pNVerts;

    ////////////////////////////////////////////////////////////////////////

    /// The initial potential set to all vertex points.
    ///
    // double                      pV0;

    /// The local potential over all mesh vertex points.
    ///
    double *                    pV;

    double *                    pGExt;

    /// Extracellular potential (?)
    ///
    double                      pVExt;

    /// Vertex-based current injections.
    ///
    double *                    pVertexInj;

    /// An array storing, for each mesh vertex, whether the potential at
    /// that vertex point should be clamped (true) or not (false).
    ///
    bool *                      pVertexClamp;

    /// Currents over triangles.
    ///
    double *                    pTriCur;

    /// Current clamps over triangles
    ///
    double * 					pTriCurClamp;

    ///
    double *                    pVertCur;

    /// Current clamp over vertices. This will complement any triangle clamp

    double *					pVertCurClamp;

    ////////////////////////////////////////////////////////////////////////
    // CONSISTENCY CHECKING
    ////////////////////////////////////////////////////////////////////////

    /// Total surface conductance... what does it mean, from where and in
    /// which context is it set and why doesn't it appear to be used
    /// anywhere outside of the chargeCheck method?
    //double                      pTotGSurf;

    ///
    double                      pInjTot;

    ////////////////////////////////////////////////////////////////////////
    // SOLUTION WORKSPACE
    ////////////////////////////////////////////////////////////////////////

    // Replace with flat arrays as they are far quicker
    /* *******************************
    double **                   pRawBDM;
    double **                   pMWK;
    */

    double *                    pRawBDM;

    double *                    pMWK;

    double *                    pDV;
    double *                    pRHS;

    int                         pHalfBW;
    int                         pBW;
    int                         maxdi;
    BandDiagonalMatrix *        pBDM;

    ////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

END_NAMESPACE(efield)
END_NAMESPACE(solver)
END_NAMESPACE(steps)

////////////////////////////////////////////////////////////////////////////////

#endif
// STEPS_SIM_EFIELD_BDMATRIXPROP_HPP

// END
