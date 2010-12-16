#define H2D_REPORT_INFO
#include "hermes2d.h"
#include <stdio.h>
#include <stdio.h>
#include <cmath>
using namespace std;
#include "python_api.h"
// This is a solver for the 2D Schroedinger equation

int P_INIT = 6;                                   
// Uniform polynomial degree of mesh elements.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_UMFPACK, SOLVER_PETSC,
                                                  // SOLVER_MUMPS, and more are coming.

// Boundary condition types.
// Note: "essential" means that solution value is prescribed.
BCType bc_types(int marker)
{
  return BC_ESSENTIAL;
}

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
{
  return 0;
}

double pot(double x, double y) {
  return x*x+y*y;
}

// Weak forms.

Python p;// python interpreter object.

#include "forms.cpp"
int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &mesh);

  // Perform initial mesh refinements (optional).
  for (int i=0;i<3;i++) mesh.refine_all_elements();

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, bc_types, essential_bc_values, P_INIT);
  int ndof = get_num_dofs(&space);
  info("ndof = %d", ndof);

  // Initialize the weak formulation for the left hand side i.e. H 
  WeakForm wfH;
  wfH.add_matrix_form(callback(bilinear_form_H));
  //wfH.add_vector_form(callback(linear_form));
  // Initialize the linear problem.
  LinearProblem lpH(&wfH, &space);

  // Select matrix solver.
  Matrix *Hmat,*Vmat,*Umat; Vector* eivec; CommonSolver* solver;
  init_matrix_solver(matrix_solver, ndof, Hmat, eivec, solver);

  // Assemble stiffness matrix
  lpH.assemble(Hmat, eivec);


  // Initialize the weak formulation for the potential matrix 
  WeakForm wfPot;
  wfPot.add_matrix_form(bilinear_form_V, bilinear_form_V_ord, H2D_SYM);
  // Initialize the linear problem.
  LinearProblem lpPot(&wfPot, &space);
  init_matrix_solver(matrix_solver, ndof, Vmat, eivec, solver);
  // Assemble potential matrix
  lpPot.assemble(Vmat,eivec);

  // Initialize the weak formulation for the right hand side i.e. U 
  WeakForm wfU;
  wfU.add_matrix_form(callback(bilinear_form_U));
  // Initialize the linear problem.
  LinearProblem lpU(&wfU, &space);
  init_matrix_solver(matrix_solver, ndof, Umat, eivec, solver);
  // Assemble overlap matrix 
  lpU.assemble(Umat,eivec);
  CSRMatrix  mat1(Hmat);
  CSRMatrix  mat2(Umat);
  CSRMatrix  mat3(Vmat);
  p.exec("print 'Python initialized'");
  p.push("hmat",c2py_CSRMatrix(&mat1));
  p.push("umat",c2py_CSRMatrix(&mat2));
  p.push("vmat",c2py_CSRMatrix(&mat3));
  p.exec("h_csr=hmat.to_scipy_csr()");
  p.exec("u_csr=umat.to_scipy_csr()");
  p.exec("v_csr=vmat.to_scipy_csr()");
  p.exec("H=h_csr.tocoo()");
  p.exec("U=u_csr.tocoo()");
  p.exec("V=v_csr.tocoo()");
  p.exec("from solvers.eigen import solve_eig_numpy, solve_eig_pysparse");
  p.exec("eigs = solve_eig_pysparse(H+V,U,6.0,n_eigs=1)");
  p.exec("E, v = eigs[0]");
  double E = py2c_double(p.pull("E"));
  int n;
  double *evec;
  numpy2c_double_inplace(p.pull("v"), &evec, &n);
  printf("E=%.16f", E);
  for (int i=0;i<ndof;i++) eivec->set(i,abs(evec[i]));// ensure that the plotted eigenvector is positive
  Solution* sln = new Solution(&space, eivec);
  ScalarView view("Solution", new WinGeom(0, 0, 1024, 768));
  // Visualize the solution.
  view.set_3d_mode();
  view.show(sln);
  // Wait for the view to be closed.
  View::wait();  
  return 0; 
};

