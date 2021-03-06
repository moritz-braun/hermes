// This file is part of Hermes3D
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// Hermes3D is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation; either version 2 of the License,
// or (at your option) any later version.
//
// Hermes3D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes3D; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#ifndef _HERMES_3D_
#define _HERMES_3D_

// hermes_common
#include "../../hermes_common/array.h"
#include "../../hermes_common/bitarray.h"
#include "../../hermes_common/map.h"
#include "../../hermes_common/maphs.h"
#include "../../hermes_common/mapord.h"
#include "../../hermes_common/trace.h"
#include "../../hermes_common/utils.h"
#include "../../hermes_common/bctypes.h"

// solvers
#include "../../hermes_common/solver/solver.h"
#include "../../hermes_common/solver/umfpack_solver.h"
#include "../../hermes_common/solver/superlu.h"
#include "../../hermes_common/solver/pardiso.h"
#include "../../hermes_common/solver/petsc.h"
#include "../../hermes_common/solver/epetra.h"
#include "../../hermes_common/solver/amesos.h"
#include "../../hermes_common/solver/aztecoo.h"
#include "../../hermes_common/solver/nox.h"
#include "../../hermes_common/solver/mumps.h"

// preconditioners
#include "../../hermes_common/solver/precond.h"
#include "../../hermes_common/solver/precond_ifpack.h"
#include "../../hermes_common/solver/precond_ml.h"

#include "h3d_common.h"

// mesh
#include "mesh.h"

// mesh loaders
#include "meshloader.h"
#include "loader/mesh3d.h"
#include "loader/hdf5.h"
#include "loader/exodusii.h"

// spaces
#include "space/space.h"
#include "space/h1.h"
#include "space/hcurl.h"

#include "order.h"
// quadrature
#include "quad.h"
#include "quadstd.h"
#include "forms.h"

#include "refmap.h"
#include "integrals/h1.h"
#include "integrals/hcurl.h"

#include "refdomain.h"

#include "shapefn.h"

// shapesets
#include "shapeset/shapeset.h"
#include "shapeset/common.h"
#include "shapeset/h1lobattotetra.h"
#include "shapeset/h1lobattohex.h"
#include "shapeset/hcurllobattohex.h"

#include "norm.h"

// output
#include "output.h"
#include "output/gmsh.h"
#include "output/vtk.h"
#include "output/graph.h"

#include "asmlist.h"
#include "solution.h"
#include "filter.h"
#include "weakform.h"
#include "discrete_problem.h"

// adapt
#include "adapt/adapt.h"
#include "adapt/h1proj.h"
#include "adapt/h1projipol.h"

#include "ogprojection.h"

// global H3D iface
void set_verbose(bool verb = true);

#endif
