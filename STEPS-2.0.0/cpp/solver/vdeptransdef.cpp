////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2007-2013�Okinawa Institute of Science and Technology, Japan.
// Copyright (C) 2003-2006�University of Antwerp, Belgium.
//
// See the file AUTHORS for details.
//
// This file is part of STEPS.
//
// STEPS�is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// STEPS�is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.�If not, see <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////

// STL headers.
#include <string>
#include <cassert>
#include <iostream>
#include <sstream>

// STEPS headers.
#include "../common.h"
#include "types.hpp"
#include "../error.hpp"
#include "statedef.hpp"
#include "vdeptransdef.hpp"
#include "../model/chanstate.hpp"
#include "../model/vdeptrans.hpp"


NAMESPACE_ALIAS(steps::solver, ssolver);
NAMESPACE_ALIAS(steps::model, smod);

////////////////////////////////////////////////////////////////////////////////

ssolver::VDepTransdef::VDepTransdef(Statedef * sd, uint idx, smod::VDepTrans * vdt)
: pStatedef(sd)
, pIdx(idx)
, pName()
, pSetupdone(false)
, pVMin(0.0)
, pVMax(0.0)
, pDV(0.0)
, pSrc()
, pDst()
, pVRateTab(0)
, pSpec_DEP(0)
, pSpec_SRCCHAN(GIDX_UNDEFINED)
, pSpec_DSTCHAN(GIDX_UNDEFINED)
{
	assert (pStatedef != 0);
	assert (vdt != 0);

    pName = vdt->getID();

    pSrc = vdt->getSrc()->getID();
    pDst = vdt->getDst()->getID();

    // Copy rate information from model object
    pVMin = vdt->_getVMin();
    pVMax = vdt->_getVMax();
    pDV = vdt->_getDV();
    uint tablesize = vdt->_getTablesize();
    assert(tablesize == static_cast<uint>(std::floor((pVMax - pVMin) / pDV)) + 1);

    pVRateTab = new double[tablesize];
    // Just temporarily store the pointer:
    //double * rates = vdt->_getRate();

    for (uint i = 0; i < tablesize; ++i)
    {
    	pVRateTab[i] = vdt->_getRate()[i];
    }

    uint nspecs = pStatedef->countSpecs();
    if (nspecs == 0) return; // Would be weird, but okay.
    pSpec_DEP = new int[nspecs];
    std::fill_n(pSpec_DEP, nspecs, DEP_NONE);
}

////////////////////////////////////////////////////////////////////////////////

ssolver::VDepTransdef::~VDepTransdef(void)
{
	delete[] pVRateTab;

	if (pStatedef->countSpecs() > 0)
	{
		delete[] pSpec_DEP;
	}

}

////////////////////////////////////////////////////////////////////////////////

void ssolver::VDepTransdef::checkpoint(std::fstream & cp_file)
{
    cp_file.write((char*)&pVMin, sizeof(double));
    cp_file.write((char*)&pVMax, sizeof(double));
    cp_file.write((char*)&pDV, sizeof(double));
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::VDepTransdef::restore(std::fstream & cp_file)
{
    cp_file.read((char*)&pVMin, sizeof(double));
    cp_file.read((char*)&pVMax, sizeof(double));
    cp_file.read((char*)&pDV, sizeof(double));
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::VDepTransdef::setup(void)
{
	assert(pSetupdone == false);

	uint sidx = pStatedef->getSpecIdx(pSrc);
	uint didx = pStatedef->getSpecIdx(pDst);

	pSpec_SRCCHAN = sidx;
	pSpec_DSTCHAN = didx;
	pSpec_DEP[sidx] |= DEP_STOICH;
	pSpec_DEP[didx] |= DEP_STOICH;

	pSetupdone = true;
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::VDepTransdef::srcchanstate(void) const
{
	assert(pSetupdone == true);
	return pSpec_SRCCHAN;
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::VDepTransdef::dstchanstate(void) const
{
	assert(pSetupdone == true);
	return pSpec_DSTCHAN;
}

////////////////////////////////////////////////////////////////////////////////

double ssolver::VDepTransdef::getVDepRate(double v) const
{
	assert(pSetupdone == true);
	assert(pVRateTab != 0);
	if (v > pVMax)
    {
        std::ostringstream os;
        os << "Voltage is higher than maximum for VDepTrans, " << name() << ": ";
        os << v << " > " << pVMax;
        throw steps::ProgErr(os.str());
    }
	if (v < pVMin)
    {
        std::ostringstream os;
        os << "Voltage is lower than maximum for VDepTrans, " << name() << ": ";
        os << v << " < " << pVMin;
        throw steps::ProgErr(os.str());
    }

	double v2 = ((v - pVMin) / pDV);
	double lv = floor(v2);
    uint lvidx = static_cast<uint>(lv);
    uint uvidx = static_cast<uint>(ceil(v2));
    double r = v2-lv;

    return (((1.0 - r) * pVRateTab[lvidx]) + (r * pVRateTab[uvidx]));

}

////////////////////////////////////////////////////////////////////////////////

int ssolver::VDepTransdef::dep(uint gidx) const
{
	assert(pSetupdone == true);
	assert(gidx < pStatedef->countSpecs());
	return pSpec_DEP[gidx];
}

////////////////////////////////////////////////////////////////////////////////

bool ssolver::VDepTransdef::req(uint gidx) const
{
	assert(pSetupdone == true);
	assert(gidx < pStatedef->countSpecs());
	if (pSpec_DEP[gidx] != DEP_NONE) return true;
	return false;
}

////////////////////////////////////////////////////////////////////////////////

// END
