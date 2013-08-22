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
#include <algorithm>
#include <iostream>
#include <string>
#include <sstream>
#include <ctime>

// STEPS headers.
#include "../../common.h"
#include "bdmatrix.hpp"
#include "bdmatrixprop.hpp"
#include "tetmesh.hpp"
#include "vertexconnection.hpp"
#include "vertexelement.hpp"

////////////////////////////////////////////////////////////////////////////////

USING_NAMESPACE(std);
NAMESPACE_ALIAS(steps::solver::efield, sefield);

////////////////////////////////////////////////////////////////////////////////

sefield::BandedMatrixProp::BandedMatrixProp(TetMesh * msh)
: pMesh(msh)
, pNVerts(pMesh->countVertices())
//, pV0(0.0)
, pGExt(0)
, pVExt(0.0)
, pVertexInj(0)
, pVertexClamp(0)
, pTriCur(0)
, pVertCur(0)
, pTriCurClamp(0)
, pVertCurClamp(0)
//, pTotGSurf(0.0)
, pInjTot(0.0)
, pRawBDM(0)
, pMWK(0)
, pDV(0)
, pRHS(0)
, pHalfBW(0)
, pBW(0)
, pBDM(0)
{
	pMesh->reindexElements(); // Is this necessary?

	pV = new double[pNVerts];
	fill_n(pV, pNVerts, 0.0);

	pGExt = new double[pNVerts];
	fill_n(pGExt, pNVerts, 0.0);

	pVertexInj = new double[pNVerts];
	fill_n(pVertexInj, pNVerts, 0.0);

    pVertexClamp = new bool[pNVerts];
    fill_n(pVertexClamp, pNVerts, false);

    int ntri = pMesh->getNTri();
	pTriCur = new double[ntri];
	fill_n(pTriCur, ntri, 0.0);
	pTriCurClamp = new double[ntri];
	fill_n(pTriCurClamp, ntri, 0.0);

	pVertCur = new double[pNVerts];
	fill_n(pVertCur, pNVerts, 0.0);

	pVertCurClamp = new double[pNVerts];
	fill_n(pVertCurClamp, pNVerts, 0.0);

	// Find out how big the band of the band-diagonal matrix should be.
	maxdi = 0;
	for (int iv = 0; iv < pNVerts; ++iv)
	{
		VertexElement* ve = pMesh->getVertex(iv);
		int ind = ve->getIDX();

		for (int i = 0; i < ve->getNCon(); ++i)
		{
			int inbr = ve->nbrIdx(i);
			int di = ind - inbr;
			if (di < 0)
			{
				di = -di;
			}
			if (di > maxdi)
			{
				maxdi = di;
			}
		}
	}
	/*
	// Report on standard output (remove later).
	stringstream ss;
	ss << "banded matrix, mesh size=";
	ss << pNVerts << "  band half width=" << maxdi;
	cout << ss.str() << endl;
	*/

	pHalfBW = maxdi;
	pBW = 2 * maxdi + 1;



    /* *******************************
	pRawBDM = new double*[pNVerts];
	for (int i = 0; i < pNVerts; ++i)
	{
		pRawBDM[i] = new double[pBW];
		fill_n(pRawBDM[i], pBW, 0.0);
	}
    */
    pRawBDM = new double[pBW*pNVerts];
    fill_n(pRawBDM, pBW*pNVerts, 0.0);

    /***************************
	pMWK = new double*[pNVerts];
	for (int i = 0; i < pNVerts; ++i)
	{
		pMWK[i] = new double[maxdi + 1];
		fill_n(pMWK[i], maxdi + 1, 0.0);
	}
     */
    pMWK = new double[pNVerts*(maxdi + 1)];
    fill_n(pMWK, pNVerts*(maxdi + 1), 0.0);



	pRHS = new double[pNVerts];
	fill_n(pRHS, pNVerts, 0.0);

	pDV = new double[pNVerts];
	fill_n(pDV, pNVerts, 0.0);

	pBDM = new BandDiagonalMatrix(pNVerts, pBW, pRawBDM, pMWK);
}

////////////////////////////////////////////////////////////////////////////////

sefield::BandedMatrixProp::~BandedMatrixProp(void)
{
    delete[] pV;
    delete[] pGExt;
    delete[] pVertexInj;
    delete[] pVertexClamp;
    // bug(swils) 22-Oct-2008: forgotten
    delete[] pTriCur;
    delete[] pVertCur;

    /* *********************
    for (int i = 0; i < pNVerts; ++i)
    {
        delete[] pRawBDM[i];
    }
    */
    delete[] pRawBDM;

    /************************
    for (int i = 0; i < pNVerts; ++i)
    {
        delete[] pMWK[i];
    }
    delete[] pMWK;
     */
    delete[] pMWK;



    delete[] pDV;
    delete[] pRHS;

    // bug(swils) 06-Aug-2008: forgotten
    delete pBDM;
}

////////////////////////////////////////////////////////////////////////////////

void sefield::BandedMatrixProp::checkpoint(std::fstream & cp_file)
{
    cp_file.write((char*)pV, sizeof(double) * pNVerts);
    cp_file.write((char*)pGExt, sizeof(double) * pNVerts);
    cp_file.write((char*)&pVExt, sizeof(double));
    cp_file.write((char*)pVertexInj, sizeof(double) * pNVerts);
    cp_file.write((char*)pVertexClamp, sizeof(bool) * pNVerts);
    int ntri = pMesh->getNTri();
    cp_file.write((char*)&ntri, sizeof(int));
    cp_file.write((char*)pTriCur, sizeof(double) * ntri);
    cp_file.write((char*)pTriCurClamp, sizeof(double) * ntri);
    cp_file.write((char*)pVertCur, sizeof(double) * pNVerts);
    cp_file.write((char*)pVertCurClamp, sizeof(double) * pNVerts);
    cp_file.write((char*)&pInjTot, sizeof(double));

    cp_file.write((char*)&pHalfBW, sizeof(int));
    cp_file.write((char*)&pBW, sizeof(int));
    cp_file.write((char*)&maxdi, sizeof(int));

    cp_file.write((char*)pRawBDM, sizeof(double) * pNVerts * pBW);
    cp_file.write((char*)pMWK, sizeof(double) * pNVerts * (maxdi + 1));
    cp_file.write((char*)pDV, sizeof(double) * pNVerts);
    cp_file.write((char*)pRHS, sizeof(double) * pNVerts);

    pBDM->checkpoint(cp_file);
}

////////////////////////////////////////////////////////////////////////////////

void sefield::BandedMatrixProp::restore(std::fstream & cp_file)
{
    cp_file.read((char*)pV, sizeof(double) * pNVerts);
    cp_file.read((char*)pGExt, sizeof(double) * pNVerts);
    cp_file.read((char*)&pVExt, sizeof(double));
    cp_file.read((char*)pVertexInj, sizeof(double) * pNVerts);
    cp_file.read((char*)pVertexClamp, sizeof(bool) * pNVerts);
    int ntri = 0;
    cp_file.read((char*)&ntri, sizeof(int));
    cp_file.read((char*)pTriCur, sizeof(double) * ntri);
    cp_file.read((char*)pTriCurClamp, sizeof(double) * ntri);
    cp_file.read((char*)pVertCur, sizeof(double) * pNVerts);
    cp_file.read((char*)pVertCurClamp, sizeof(double) * pNVerts);
    cp_file.read((char*)&pInjTot, sizeof(double));

    cp_file.read((char*)&pHalfBW, sizeof(int));
    cp_file.read((char*)&pBW, sizeof(int));
    cp_file.read((char*)&maxdi, sizeof(int));

    cp_file.read((char*)pRawBDM, sizeof(double) * pNVerts * pBW);
    cp_file.read((char*)pMWK, sizeof(double) * pNVerts * (maxdi + 1));
    cp_file.read((char*)pDV, sizeof(double) * pNVerts);
    cp_file.read((char*)pRHS, sizeof(double) * pNVerts);

    pBDM->restore(cp_file);
}

////////////////////////////////////////////////////////////////////////////////

void sefield::BandedMatrixProp::advance(double dt)
{
    populateMRHS(dt);

    // Performance: the vast majority of the work is in bdm->lu()
    pBDM->lu();
    pBDM->lubksb(pRHS, pDV);

    for (uint i = 0; i < pNVerts; ++i)
    {
        // bugfix(wils) 06-Aug-2008: the voltage-clamp variable was being
        // stored for each vertex (in variable pVertexClamp), but not used
        // during the update step.
        if (pVertexClamp[i] == false)
        {
            pV[i] += pDV[i];
        }
    }

    // Has strange behaviours, uses unexplainable variables such as this
    // pInjTot. Look at this later...
    //chargeCheck(dt);
}

////////////////////////////////////////////////////////////////////////////////
/*
void BandedMatrixProp::setInitialPotential(double d)
{
	// NOTE: getting rid of this- replace with setPotential, to behave a little
	// like reset(), and no need to store pV0;
	pV0 = d;
}
*/
////////////////////////////////////////////////////////////////////////////////

void sefield::BandedMatrixProp::setSurfaceConductance(double gsa, double vr)
{
	pVExt = vr;
	// pTotGSurf = 0.;
	double totsurf = 0.;
	for (int i = 0; i < pNVerts; ++i)
	{
		VertexElement* ve = pMesh->getVertex(i);
		double ge = gsa * ve->getSurfaceArea();
		pGExt[ve->getIDX()] = ge;
		// pTotGSurf += ge;

		totsurf += ve->getSurfaceArea();
	}

	/*
	stringstream ss;
	ss << "total surface area " << totsurf;
	cout << ss.str() << endl;
	*/
}

////////////////////////////////////////////////////////////////////////////////

string sefield::BandedMatrixProp::makeOutputLine(double time, double* v)
{
	stringstream s;
	s << time << " " << pV[0] << " " << pV[pNVerts-1];
	return s.str();
}

////////////////////////////////////////////////////////////////////////////////
/*
void BandedMatrixProp::reset(void)
{
    for (int i = 0; i < pNVerts; ++i)
    {
        pV[i] = pV0;
    }
}
*/
////////////////////////////////////////////////////////////////////////////////

void sefield::BandedMatrixProp::setPotential(double v)
{
    for (int i = 0; i < pNVerts; ++i)
    {
        pV[i] = v;
    }
}

////////////////////////////////////////////////////////////////////////////////

void sefield::BandedMatrixProp::populateMRHS(double dt)
{
	// Currents into vertices are the per-vertex injections plus
	// the per-triangle currents divided over their neighbors.

	// NOTE: at the moment pVertexInj is not used at all, and
	// this is a little inefficient at the moment- i.e. looping
	// over all vertices when typically only a few will have
	// a current clamp, if any
	for (uint i = 0; i < pNVerts; ++i)
	{
		pVertCur[i] = pVertexInj[i];
		pVertCur[i] += pVertCurClamp[i];

	}
	uint ntris = pMesh->getNTri();
	for (uint i = 0; i < ntris; ++i)
	{
		uint * triv = pMesh->getTriangle(i);
		double c = pTriCur[i] / 3.0;
		double cc = pTriCurClamp[i] / 3.0;

		pVertCur[triv[0]] += (c+cc);
		pVertCur[triv[1]] += (c+cc);
		pVertCur[triv[2]] += (c+cc);
	}


	// NOTE: Vertex currents at this point are given in pA

	for (uint i = 0; i < pNVerts; ++i)
	{
		pRHS[i] = pVertCur[i] * dt;
	}
	// NOTE: time is in units of ms

	int nw = 2 * pHalfBW + 1;
	for (int i = 0; i < pNVerts; ++i)
	{
		for (int j = 0; j < nw; ++j)
		{
            /* *************************
			pRawBDM[i][j] = 0.0;
            */
            pRawBDM[(i*nw)+j] = 0.0;
		}
	}
	for (int ivert = 0; ivert < pNVerts; ++ivert)
	{
		VertexElement * ve = pMesh->getVertex(ivert);

		int ind = ve->getIDX();

		// NOTE: In following units are:
		// Time: ms
		// External conductance: ??
		// Potentials: mV
		// Capacitance: pF

		// the diagonal term (zero for all the internal points)
        /* ***************************************
		pRawBDM[ind][pHalfBW] += ve->getCapacitance() + dt * pGExt[ind];
         */
        pRawBDM[(ind*nw)+pHalfBW] += ve->getCapacitance() + dt * pGExt[ind];

		pRHS[ind] += dt * pGExt[ind] * (pVExt - pV[ind]);

		// Now, loop through all the neighbours adding on contributions
		// to matrix and pRHS.
		for (int inbr = 0; inbr < ve->getNCon(); ++inbr)
		{
			int k = ve->nbrIdx(inbr);
			double cc = ve->getCC(inbr);
			// right hand side
			pRHS[ind] += dt * cc * (pV[k] - pV[ind]);

			/* if (pV[k] != pV[ind])
			{
			    //std::cerr << ind << "->" << k << ": ";
			    //std::cerr << "{cc=" << cc << "} * ";
			    //std::cerr << "{dt=" << dt << "} * ";
			    //std::cerr << "{pV=" << pV[ind] << "} * ";
			    //std::cerr << "{pVn=" << pV[k] << "} = ";
			    //std::cerr << dt * cc * (pV[k] - pV[ind]) << "\n";
			}
			*/

			// conductance terms for the doagonal
            /* ***************************************
			pRawBDM[ind][pHalfBW] += dt * cc;
             */
            pRawBDM[(ind*nw)+pHalfBW] += dt * cc;

			// other ends of the conductances
            /* ***************************************
			pRawBDM[ind][k - ind + pHalfBW] -= dt * cc;
             */
            pRawBDM[(ind*nw)+ (k-ind+pHalfBW)]-= dt * cc;

		}
	}

	// Iain 31/8/2011 Now reset the currents
	fill_n(pVertexInj, pNVerts, 0.0);
	fill_n(pTriCur, ntris, 0.0);
	fill_n(pVertCur, pNVerts, 0.0);
}

////////////////////////////////////////////////////////////////////////////////

void sefield::BandedMatrixProp::chargeCheck(double dt)
{
    //if (pTotGSurf == 0.0)
    //{
        // Can do a charge conservation check -
        // TODO: smarter check that adds currents
        double cvt = 0.;
        for (int i = 0; i < pNVerts; ++i)
        {
            VertexElement * ve = pMesh->getVertex(i);
            pInjTot += pVertexInj[i] * dt;
            // NOTE: IH 5/10 removed pV0 in below:
            // cvt += (pV[i] - pV0) * ve->getCapacitance();
            cvt += (pV[i]) * ve->getCapacitance();
        }
        double dabs = abs((cvt - pInjTot) / (cvt + pInjTot));
        if (dabs > 1.e-7)
        {
            stringstream s;

            s << "NEAR FATAL sum error too large ";
            s << dabs << " " << cvt << " " << pInjTot;
            cout << s.str() << endl;
        }
    //}
}

////////////////////////////////////////////////////////////////////////////////

// END

