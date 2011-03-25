// This file is part of Hermes2D.
//
// Hermes2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D.  If not, see <http://www.gnu.org/licenses/>.

#ifndef __H2D_LINEAR_ELASTICITY_WEAK_FORMS_H
#define __H2D_LINEAR_ELASTICITY_WEAK_FORMS_H

#include "../integrals/integrals_h1.h"

/* Default weak form for linear elasticity (Lame equations)
   with Dirichlet and/or zero Neumann BC

   Nonzero Neumann and Newton boundary conditions can be enabled 
   by creating a descendant and adding surface forms to it. 
*/

class MatrixFormVolLinearElasticity_x_x : public WeakForm::MatrixFormVol
{
public:
  MatrixFormVolLinearElasticity_x_x(unsigned int i, unsigned int j, double lambda, double mu)
        : WeakForm::MatrixFormVol(i, j, HERMES_SYM), lambda(lambda), mu(mu) { }
  MatrixFormVolLinearElasticity_x_x(unsigned int i, unsigned int j, std::string area, double lambda, double mu)
        : WeakForm::MatrixFormVol(i, j, HERMES_SYM, area), lambda(lambda), mu(mu) { }

  template<typename Real, typename Scalar>
  Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                     Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
    return (lambda + 2*mu) * int_dudx_dvdx<Real, Scalar>(n, wt, u, v) +
                        mu * int_dudy_dvdy<Real, Scalar>(n, wt, u, v);
  }

  scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
               Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
    return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
  }

  Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
          Geom<Ord> *e, ExtData<Ord> *ext) {
    return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
  }

private:
    double lambda, mu;
};

class MatrixFormVolLinearElasticity_x_y : public WeakForm::MatrixFormVol
{
public:
  MatrixFormVolLinearElasticity_x_y(unsigned int i, unsigned int j, double lambda, double mu)
          : WeakForm::MatrixFormVol(i, j, HERMES_SYM), lambda(lambda), mu(mu) {}
  MatrixFormVolLinearElasticity_x_y(unsigned int i, unsigned int j, std::string area, double lambda, double mu)
        : WeakForm::MatrixFormVol(i, j, HERMES_SYM, area), lambda(lambda), mu(mu) { }

  template<typename Real, typename Scalar>
  Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                     Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
  {
    return lambda * int_dudy_dvdx<Real, Scalar>(n, wt, u, v) +
               mu * int_dudx_dvdy<Real, Scalar>(n, wt, u, v);
  }

  scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
               Func<double> *v, Geom<double> *e, ExtData<scalar> *ext)
  {
    return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
  }

  Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
          Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
  {
     return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
  }

private:
  double lambda, mu;
};

class MatrixFormVolLinearElasticity_y_y : public WeakForm::MatrixFormVol
{
public:
  MatrixFormVolLinearElasticity_y_y(unsigned int i, unsigned int j, double lambda, double mu)
          : WeakForm::MatrixFormVol(i, j, HERMES_SYM), lambda(lambda), mu(mu) { }
  MatrixFormVolLinearElasticity_y_y(unsigned int i, unsigned int j, std::string area, double lambda, double mu)
        : WeakForm::MatrixFormVol(i, j, HERMES_SYM, area), lambda(lambda), mu(mu) { }

  template<typename Real, typename Scalar>
  Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                     Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
  {
    return              mu * int_dudx_dvdx<Real, Scalar>(n, wt, u, v) +
           (lambda + 2*mu) * int_dudy_dvdy<Real, Scalar>(n, wt, u, v);
  }

  scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
               Func<double> *v, Geom<double> *e, ExtData<scalar> *ext)
  {
    return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
  }

  Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
          Geom<Ord> *e, ExtData<Ord> *ext)
  {
     return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
  }

private:
  double lambda, mu;
};

#endif
