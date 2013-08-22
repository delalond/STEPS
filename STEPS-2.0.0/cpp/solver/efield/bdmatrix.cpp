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
/////////////////////////////////////////////////////////////////////////////

// STL headers.
#include <iostream>
#include <cmath>

// STEPS headers.
#include "../../common.h"
#include "bdmatrix.hpp"

////////////////////////////////////////////////////////////////////////////////

USING_NAMESPACE(std);
NAMESPACE_ALIAS(steps::solver::efield, sefield);

////////////////////////////////////////////////////////////////////////////////

sefield::BandDiagonalMatrix::BandDiagonalMatrix
(
    int nrow,
    int ncol,
    /*
    double** ain,
    */
    double* ain,

    /*
	double** alin
     */
    double* alin

)
{
	a = ain;
	al = alin;

	n = nrow;
	halfbw = (ncol - 1) / 2;

    std::cout << "\nHalf bandwidth: " << halfbw;
	perm = new int[nrow];
	ws = new double[nrow];
}

////////////////////////////////////////////////////////////////////////////////

sefield::BandDiagonalMatrix::~BandDiagonalMatrix(void)
{
    delete[] perm;
    delete[] ws;
}

////////////////////////////////////////////////////////////////////////////////

void sefield::BandDiagonalMatrix::checkpoint(std::fstream & cp_file)
{
    cp_file.write((char*)perm, sizeof(int) * n);
    cp_file.write((char*)ws, sizeof(double) * n);
}

////////////////////////////////////////////////////////////////////////////////

void sefield::BandDiagonalMatrix::restore(std::fstream & cp_file)
{
    cp_file.read((char*)perm, sizeof(int) * n);
    cp_file.read((char*)ws, sizeof(double) * n);
}

////////////////////////////////////////////////////////////////////////////////

void sefield::BandDiagonalMatrix::lu(void)
{
	double TINY = 1.0e-20;

	int w = 2 * halfbw + 1;
    // w is equivalent to pBW in the prop. 'a' dimensions are nrow ('n') x 'w'

	int p = halfbw;

	for (int i = 0; i < halfbw; ++i)
	{
		for (int j = halfbw - i; j < w; ++j)
		{
            /*
			a[i][j - p] = a[i][j];
             */
            a[(i*w) +(j - p)] = a[(i*w) +j];
		}
		p = p - 1;
		for (int j = w - p - 1; j < w; ++j)
		{
            /*
			a[i][j] = 0.;
             */
            a[(i*w) +j] = 0.0;
		}
	}

	double d = 1.;

	p = halfbw;

    double * dpointerk = &a[0];

	for (int k = 0; k < n; ++k)
	{
        /*
		double dum = a[k][0];
		*/
        double dum = a[k*w];

        int ipiv = k;
		if (p < n)
		{
			p=p+1;
		}
		for (int j=k+1; j < p; ++j)
		{
			//Find the pivot element.
            /*
			if (abs(a[j][0]) > abs(dum))
             */
            if (abs(a[j*w]) > abs(dum))
			{
                /*
				dum=a[j][0];
                 */
                dum=a[j*w];

				ipiv = j;
			}
		}

		perm[k] = ipiv;

		if (dum == 0.)
		{
            /*
			a[k][0] = TINY;
             */
            a[k*w] = TINY;

		}

		if (ipiv != k)
		{
			// swap rowS
			d=-d;
			for (int j = 0; j < w; ++j)
			{
                /*
				dum=a[k][j];
				a[k][j] = a[ipiv][j];
				a[ipiv][j] = dum;
                 */
                dum=a[(k*w)+j];
				a[(k*w)+j] = a[(ipiv*w)+j];
				a[(ipiv*w)+j] = dum;

			}
		}


        double * dpointeri = &a[0];

        dpointeri = dpointeri+((k+1)*w);

		// now for the eliminiation
		for (int i = k+1; i < p; ++i)
		{
            /*
			dum=a[i][0] / a[k][0];
             */
            //dum=a[i*w] / a[k*w];

            dum = (*dpointeri)/(*dpointerk);

            /*
			al[k][i-k-1] = dum;
             */
            al[k*(halfbw+1)+i-k-1] = dum;


			for (int j = 1; j < w; ++j)
			{
                /*
                double negt =  a[k][j] * dum;
                double numb = a[i][j] - negt;
                a[i][j-1] = numb;
                 */
                //int idx1 = (k*w)+j;
                //int idx2 = (i*w)+j;
                //double negt =  a[idx1] * dum;
                //double numb = a[idx2] - negt;
                //a[idx2-1] = numb;

                double * dptempk = dpointerk+j;
                double negt = *(dptempk) * dum;
                double numb = *(dpointeri+1) - negt;
                *(dpointeri)= numb;

                dpointeri++;

                // original
				//a[i][j-1] = a[i][j] - dum * a[k][j];
			}
            // We have to increment one more time for a total of w
			// increments because j starts at 1 in above loop
            dpointeri++;
            /*
			a[i][w-1] = 0.;
             */
            a[(i*w)+(w-1)]=0.0;

		}
    dpointerk = dpointerk+w;
	}

}

////////////////////////////////////////////////////////////////////////////////

void sefield::BandDiagonalMatrix::lubksb(double * bin, double * b)
{
	for (int i = 0; i < n; ++i)
	{
		b[i] = bin[i];
	}

	int w = 2 * halfbw + 1;
	int p = halfbw;

	for (int k = 0; k < n; ++k)
	{
		int ipiv = perm[k];
		if (ipiv != k)
		{
			double dum = b[k];
			b[k] = b[ipiv];
			b[ipiv] = dum;
		}
		if (p < n)
		{
			p=p+1;
		}
		for (int i= k+1; i < p; ++i)
		{
			/*
            b[i] -= al[k][i-k-1] * b[k];
             */
            b[i] -= al[k*(halfbw+1)+i-k-1]*b[k];
		}
	}

	p = 1;
	// Backsubstitution.
	for (int i = n-1; i >= 0; --i)
	{
		double dum = b[i];
		for (int k=1; k < p; ++k)
		{
            /*
			dum = dum - a[i][k] * b[k+i];
             */
            dum = dum - a[(i*w)+k] * b[k+i];
		}

        /******************
		b[i] = dum / a[i][0];
         */
        b[i] = dum/a[i*w];

		if (p < w)
		{
			p=p+1;
		}
	}
}

////////////////////////////////////////////////////////////////////////////////

// END
