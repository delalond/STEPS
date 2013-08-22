# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# STEPS - STochastic Engine for Pathway Simulation
# Copyright (C) 2007-2013 Okinawa Institute of Science and Technology, Japan.
# Copyright (C) 2003-2006 University of Antwerp, Belgium.
#
# See the file AUTHORS for details.
#
# This file is part of STEPS.
#
# STEPS is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# STEPS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#  Last Changed Rev:  $Rev: 489 $
#  Last Changed Date: $Date: 2013-04-19 19:37:42 +0900 (Fri, 19 Apr 2013) $
#  Last Changed By:   $Author: iain $

try:
    from setuptools import setup, Extension
    
except ImportError:
    from distutils.core import setup
    from distutils.extension import Extension


def name():
    return 'STEPS'
    
def version():
    return '2.0.0'
    
def author():
    return 'STEPS Development Team'
    
def email():
    return 'steps.dev@gmail.com'
    
def url():
    return 'http://steps.sourceforge.net'
    
def desc():
    return 'STochastic Engine for Pathway Simulation'

def download():
    return 'http://sourceforge.net/projects/steps/files/src'

def platforms():
    return ['Mac OS X', 'Windows XP', 'Windows Vista', 'Linux', 'Unix']
    
def license():
    return 'GNU General Public License Version 3.0'
    
def packages():
    return ['steps', 'steps/utilities']
  
def steps_ext():
    return dict(
        name='_steps_swig',
        
        sources=['cpp/error.cpp',
                 
                 'cpp/tetode/comp.cpp', 'cpp/tetode/patch.cpp', 'cpp/tetode/tet.cpp', 
                 'cpp/tetode/tri.cpp', 'cpp/tetode/tetode.cpp',
                 
                 'third_party/cvode-2.6.0/src/cvode/cvode_band.c',
                 'third_party/cvode-2.6.0/src/cvode/cvode_bandpre.c',
                 'third_party/cvode-2.6.0/src/cvode/cvode_bbdpre.c',
                 'third_party/cvode-2.6.0/src/cvode/cvode_dense.c',
                 'third_party/cvode-2.6.0/src/cvode/cvode_diag.c',
                 'third_party/cvode-2.6.0/src/cvode/cvode_direct.c',
                 'third_party/cvode-2.6.0/src/cvode/cvode_io.c', #'third_party/cvode-2.6.0/src/cvode/cvode_lapack.c',
                 'third_party/cvode-2.6.0/src/cvode/cvode_spbcgs.c',
                 'third_party/cvode-2.6.0/src/cvode/cvode_spgmr.c',
                 'third_party/cvode-2.6.0/src/cvode/cvode_spils.c',
                 'third_party/cvode-2.6.0/src/cvode/cvode_sptfqmr.c',
                 'third_party/cvode-2.6.0/src/cvode/cvode.c',
                 
                 'third_party/cvode-2.6.0/src/nvec_ser/fnvector_serial.c',
                 'third_party/cvode-2.6.0/src/nvec_ser/nvector_serial.c',
                 
                 'third_party/cvode-2.6.0/src/sundials/sundials_band.c',
                 'third_party/cvode-2.6.0/src/sundials/sundials_dense.c',
                 'third_party/cvode-2.6.0/src/sundials/sundials_direct.c',
                 'third_party/cvode-2.6.0/src/sundials/sundials_iterative.c',
                 'third_party/cvode-2.6.0/src/sundials/sundials_math.c',
                 'third_party/cvode-2.6.0/src/sundials/sundials_nvector.c',
                 'third_party/cvode-2.6.0/src/sundials/sundials_spbcgs.c',
                 'third_party/cvode-2.6.0/src/sundials/sundials_spgmr.c',
                 'third_party/cvode-2.6.0/src/sundials/sundials_sptfqmr.c',
                 
                 'cpp/solver/api_comp.cpp','cpp/solver/api_main.cpp', 'cpp/solver/api_memb.cpp',
                 'cpp/solver/api_patch.cpp','cpp/solver/api_tet.cpp', 'cpp/solver/api_vert.cpp',
                 'cpp/solver/api_tri.cpp', 'cpp/solver/api_diffboundary.cpp',
                 'cpp/solver/compdef.cpp',
                 'cpp/solver/diffdef.cpp','cpp/solver/surfdiffdef.cpp','cpp/solver/patchdef.cpp',
                 'cpp/solver/reacdef.cpp','cpp/solver/specdef.cpp',
                 'cpp/solver/sreacdef.cpp','cpp/solver/statedef.cpp',
                 'cpp/solver/chandef.cpp', 'cpp/solver/ghkcurrdef.cpp','cpp/solver/diffboundarydef.cpp',
                 'cpp/solver/ohmiccurrdef.cpp', 'cpp/solver/vdeptransdef.cpp', 'cpp/solver/vdepsreacdef.cpp',
                 
                 'cpp/solver/efield/bdmatrix.cpp', 'cpp/solver/efield/bdmatrixprop.cpp',
                 'cpp/solver/efield/efield.cpp', 'cpp/solver/efield/matrix.cpp',
                 'cpp/solver/efield/tetcoupler.cpp', 'cpp/solver/efield/tetmesh.cpp',
                 'cpp/solver/efield/vertexconnection.cpp', 'cpp/solver/efield/vertexelement.cpp',
                 
                 'cpp/tetexact/comp.cpp','cpp/tetexact/diff.cpp', 'cpp/tetexact/sdiff.cpp',
                 'cpp/tetexact/kproc.cpp','cpp/tetexact/patch.cpp',
                 'cpp/tetexact/reac.cpp','cpp/tetexact/sreac.cpp',
                 'cpp/tetexact/tet.cpp','cpp/tetexact/tetexact.cpp',
                 'cpp/tetexact/tri.cpp',
                 'cpp/tetexact/ghkcurr.cpp',
                 'cpp/tetexact/vdeptrans.cpp', 'cpp/tetexact/vdepsreac.cpp',
                 'cpp/tetexact/diffboundary.cpp', 
                 'cpp/tetexact/wmvol.cpp',
                 
                 'cpp/wmdirect/comp.cpp','cpp/wmdirect/kproc.cpp',
                 'cpp/wmdirect/patch.cpp','cpp/wmdirect/reac.cpp',
                 'cpp/wmdirect/sreac.cpp','cpp/wmdirect/wmdirect.cpp',
                 
                 'cpp/wmrk4/wmrk4.cpp',
                                  
                 'cpp/model/model.cpp', 'cpp/model/diff.cpp', 'cpp/model/chan.cpp',
                 'cpp/model/reac.cpp','cpp/model/spec.cpp','cpp/model/sreac.cpp',
                 'cpp/model/surfsys.cpp','cpp/model/volsys.cpp',
                 'cpp/model/chanstate.cpp','cpp/model/ohmiccurr.cpp', 
                 'cpp/model/ghkcurr.cpp', 'cpp/model/vdeptrans.cpp', 'cpp/model/vdepsreac.cpp',
                 
                 'cpp/math/tetrahedron.cpp', 'cpp/math/tools.cpp',
                 'cpp/math/linsolve.cpp','cpp/math/triangle.cpp','cpp/math/ghk.cpp',
                 
                 'cpp/geom/comp.cpp','cpp/geom/geom.cpp','cpp/geom/patch.cpp',
                 'cpp/geom/tetmesh.cpp','cpp/geom/tet.cpp',
                 'cpp/geom/tmcomp.cpp','cpp/geom/tmpatch.cpp','cpp/geom/tri.cpp',
                 'cpp/geom/memb.cpp',  'cpp/geom/diffboundary.cpp',
                 
                 'cpp/rng/rng.cpp', 'cpp/rng/mt19937.cpp',
                 
                 'swig/steps_wrap.cpp'],
                #define_macros=[('SSA_DEBUG', 'None')],
            undef_macros=['NDEBUG']
        )
        
def ext_modules():
    modules = [steps_ext()]
    return modules

ExtModule = lambda extension:  Extension(**extension)

setup(name = name(),
      version = version(),
      author = author(),
      author_email = email(),
      url = url(),
      description = desc(),
      download_url = download(),
      platforms = platforms(),
      license = license(),
      
      packages = packages(), 
      
      ext_package = 'steps',
      
      ext_modules  = [ExtModule(ext) for ext in ext_modules()]
      )


