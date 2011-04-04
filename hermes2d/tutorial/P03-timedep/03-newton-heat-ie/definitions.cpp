#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"

class CustomWeakFormHeatTransferNonlinear : public WeakForm
{
public:
  CustomWeakFormHeatTransferNonlinear(double alpha, double tau, Solution* sln_prev_time) : WeakForm(1) {
    add_matrix_form(new CustomMatrixFormVolHeatTransferNonlinear(0, 0, alpha, tau));

    VectorFormVolHeatTransfer* vector_form = new VectorFormVolHeatTransfer(0, alpha, tau);
    vector_form->ext.push_back(sln_prev_time);
    add_vector_form(vector_form);
  };

private:
  class CustomMatrixFormVolHeatTransferNonlinear : public WeakForm::MatrixFormVol
  {
  public:
    // This weak form is custom since it contains a nonlinearity in the diffusion term.
    CustomMatrixFormVolHeatTransferNonlinear(int i, int j, double alpha, double tau) 
          : WeakForm::MatrixFormVol(i, j, HERMES_NONSYM), alpha(alpha), tau(tau) { }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
      Scalar result = 0;
      Func<Scalar>* u_prev_newton = u_ext[0];
      for (int i = 0; i < n; i++)
        result += wt[i] * (u->val[i] * v->val[i] / tau + dlam_du<Real>(u_prev_newton->val[i]) * u->val[i] *
                           (u_prev_newton->dx[i] * v->dx[i] + u_prev_newton->dy[i] * v->dy[i])
                           + lam<Real>(u_prev_newton->val[i]) * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]));
      return result;
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
      return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, ExtData<Ord> *ext) const {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }

    // Thermal conductivity (temperature-dependent)
    // Note: for any u, this function has to be positive.
    template<typename Real>
    Real lam(Real u) { 
      return 1 + pow(u, alpha); 
    }

    // Derivative of the thermal conductivity with respect to 'u'.
    template<typename Real>
    Real dlam_du(Real u) { 
      return alpha*pow(u, alpha - 1); 
    }
    
    // Members.
    double alpha;
    double tau;
  };

  // This form (residual) is custom since it contains a nonlinear term.
  class VectorFormVolHeatTransfer : public WeakForm::VectorFormVol
  {
  public:
    VectorFormVolHeatTransfer(int i, double alpha, double tau) : WeakForm::VectorFormVol(i), alpha(alpha), tau(tau) { }

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
      Scalar result = 0;
      Func<Scalar>* u_prev_newton = u_ext[0];
      Func<Scalar>* u_prev_time = ext->fn[0];
      for (int i = 0; i < n; i++)
        result += wt[i] * ((u_prev_newton->val[i] - u_prev_time->val[i]) * v->val[i] / tau +
                          lam<Real>(u_prev_newton->val[i]) * (u_prev_newton->dx[i] * v->dx[i] + u_prev_newton->dy[i] * v->dy[i])
		           - heat_src<Real>(e->x[i], e->y[i]) * v->val[i]);
      return result;
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
      return vector_form<scalar, scalar>(n, wt, u_ext, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }

    // Heat sources (can be a general function of 'x' and 'y').
    template<typename Real>
    Real heat_src(Real x, Real y) {
      return 1.0;
    }

    // Thermal conductivity (temperature-dependent)
    // Note: for any u, this function has to be positive.
    template<typename Real>
    Real lam(Real u)  { 
      return 1 + pow(u, alpha); 
    }

    // Members.
    double alpha;
    double tau;
  };
};

/* Essential boundary condition */

class EssentialBCNonConst : public EssentialBoundaryCondition {
public:
  EssentialBCNonConst(std::string marker) : EssentialBoundaryCondition(Hermes::vector<std::string>())
  {
    markers.push_back(marker);
  }

  ~EssentialBCNonConst() {};

  inline EssentialBCValueType get_value_type() const { return EssentialBoundaryCondition::BC_FUNCTION; }

  virtual scalar value(double x, double y) const
  {
    return (x+10)*(y+10)/100.;
  }
};

/* Initial condition */

class CustomInitialCondition : public ExactSolutionScalar
{
public:
  CustomInitialCondition(Mesh* mesh) : ExactSolutionScalar(mesh) {};

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const {
    dx = (y+10)/10.;
    dy = (x+10)/10.;
  };

  virtual scalar value (double x, double y) const {
    return (x+10)*(y+10)/100.;
  }

  virtual Ord ord(Ord x, Ord y) const {
    return (x+10)*(y+10)/100.;
  }
};
