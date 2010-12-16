#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#include "python_api.h"
#include "hermes2d.h"
#include <stdio.h>
#include <cmath>
using namespace std;

// This is a solver for the Schroedinger equation of
// the Hydrogen molecular Ion in cylindrical coordinates
// the nuclei are at (0,1) and (0,-1) in atomic units.
// to make the potential less singular 
// the wave function is written as a product of a function satisfying
// the cusp condition at the nuclei and a smoother function to be determined via the method of finite elements
// on the domain [0,16]x[0,16] in cylindrical coordinates rho and z.
// The boundary conditions are 
// essential for rho=rhomax and z=zmax and
// natural  for rho=0 and z=0.

int P_INIT =6;  // Uniform polynomial degree of mesh elements.
int zp = 1;                                       // zparity +1/-1
const int INIT_REF_NUM = 4; // DON'T CHANGE THIS!!
const int FINAL_REF_NUM = 1; // final global refinement 
//int final_ref_num=2; // now a variable !
const double E_LITERATURE=-2.2052684289898; 
// value from 
// J. Chem. Phys. 43, 3004 (1965); doi:10.1063/1.1697265.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_UMFPACK, SOLVER_PETSC,
                                                  // SOLVER_MUMPS, and more are coming.

// Boundary condition types.
// Note: "essential" means that solution value is prescribed.
BCType bc_types(int marker)
{
  if (marker == 2 or marker == 3)
    return BC_ESSENTIAL;
  if (marker == 1 and zp == -1) 
    return BC_ESSENTIAL;
  if (marker == 1 and zp == 1)
    return BC_NATURAL;
  if (marker == 4 )
    return BC_NATURAL;
}

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
{
  return 0;
}
const int Nnuc=2; // number of nuclei, here 2
double C[Nnuc]= { 1.0186573603637741, 1.0186573603637741  }; // coefficients in cusp factor 
double Rnuc[Nnuc][3]={  {0.0,0.0,-1.0}, {0.0,0.0,1.0} }; // coordinates of the two nuclei
// the following functions are needed for the cusp factor  version of h2plus
double ri(int i,double rho,double z){
  return sqrt(rho*rho+(z-Rnuc[i][2])*(z-Rnuc[i][2]));
}// this function gives the distances from the point (rho,z)  to the i-th nucleus

double f(double rho,double z){
  return 1.0+C[0]*exp(-2.0*ri(0,rho,z))+C[1]*exp(-2.0*ri(1,rho,z));
}// this is the factor that satisfies the cusp conditions at both nuclei
double laplacef(double rho,double z){
  return 4.0*C[0]*exp(-2.0*ri(0,rho,z))+4.0*C[1]*exp(-2.0*ri(1,rho,z))\
    -4.0/ri(0,rho,z)*C[0]*exp(-2.0*ri(0,rho,z))-4.0/ri(1,rho,z)*C[1]*exp(-2.0*ri(1,rho,z));
}// this is laplace of f

double pot(double rho,double z){
  return -2.0/ri(0,rho,z)-2.0/ri(1,rho,z)-laplacef(rho,z)/f(rho,z);
}// this is the potential with the singularities at the nuclei removed and the 1/r replace by a kink 

double  wfun (double rho,double z){
  return f(rho,z)*f(rho,z);
}// this is the weight function f**2

Python p;// python interpreter object;

#include "forms.cpp"
int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;

  // If supplied take argv[1] as the name of the mesh file.
  char* meshname; 
  if (argc > 1){
    meshname=argv[1];
  }
  else{
    meshname="domain.mesh";
  }
  
  mloader.load(meshname, &mesh);
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();
  // Perform initial mesh refinements.
  // mesh.refine_element(0,0);
  // This refines element [0,8]x[0,8] 
  for (int i=0;i<2;i++) mesh.refine_all_elements();
  //The following code finds the 
  // vertex id at the nuclear coordinate (0.0,1.0).
  int NUM_ELEMS=mesh.get_max_element_id();
  int vertex_id_nucleus=-1;
  for (int id=0;id<NUM_ELEMS;id++) {
    Element* elem=mesh.get_element(id);
    for(int i=0;i<4;i++){
      Node* node= elem->vn[i];
      if ( node->x == 0.0 and node->y == 1.0) 
	vertex_id_nucleus=node->id;
    }
  }
  if (vertex_id_nucleus == -1)
    {
      verbose("nucleus is not a vertex! Slow convergence would result");
      exit(0);
    }
  printf("vertex_id_nucleus=%d\n",vertex_id_nucleus);
  mesh.refine_towards_vertex(vertex_id_nucleus,10);
  for (int i=0;i<FINAL_REF_NUM;i++) mesh.refine_all_elements();
  MeshView mview("Final Grid", new WinGeom(0, 0, 1200,900));
  mview.show(&mesh);
  View::wait();
  exit(1);
  // Create an H1 space with default shapeset.
  H1Space space(&mesh, bc_types, essential_bc_values, P_INIT);
  int ndof = get_num_dofs(&space);
  NUM_ELEMS=mesh.get_max_element_id();
  info("ndof = %d", ndof);
  info("nelems = %d", NUM_ELEMS);
  // Initialize the weak formulation for the left hand side i.e. H 
  WeakForm wfH;;
  wfH.add_matrix_form(bilinear_form_H, bilinear_form_ord, H2D_SYM);
  // Initialize the linear problem.
  LinearProblem lpH(&wfH, &space);

  // Select matrix solver.
  Matrix *Hmat,*Vmat,*Umat;
  Vector* eivec; CommonSolver* solver;
  init_matrix_solver(matrix_solver, ndof, Hmat, eivec, solver);

  // Assemble stiffness matrix
  lpH.assemble(Hmat, eivec);

  // Initialize the weak formulation for the potential matrix 
  WeakForm wfPot;
  wfPot.add_matrix_form(bilinear_form_V, bilinear_form_ord, H2D_SYM);
  // Initialize the linear problem.
  LinearProblem lpPot(&wfPot, &space);
  init_matrix_solver(matrix_solver, ndof, Vmat, eivec, solver);
  // Assemble potential matrix
  lpPot.assemble(Vmat,eivec);
  // Initialize the weak formulation for the right hand side i.e. U 
  WeakForm wfU;
  wfU.add_matrix_form(bilinear_form_U, bilinear_form_ord, H2D_SYM);
  // Initialize the linear problem.
  LinearProblem lpU(&wfU, &space);
  init_matrix_solver(matrix_solver, ndof, Umat, eivec, solver);
  // Assemble overlap matrix 
  lpU.assemble(Umat,eivec);
  cpu_time.tick();
  verbose("Total running time for assembling matrices : %g s", cpu_time.accumulated());
  cpu_time.reset();
  CSRMatrix  mat1(Hmat);
  CSRMatrix  mat2(Umat);
  CSRMatrix  mat3(Vmat);
  cpu_time.tick();
  verbose("Time for converting matrices to CSR  : %g s", cpu_time.accumulated());
  cpu_time.reset();
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
  p.exec("eigs = solve_eig_pysparse(H+V,U,-2.1,n_eigs=1)");
  p.exec("E, v = eigs[0]");
  double E = py2c_double(p.pull("E"));
  int n;
  double *evec;
  numpy2c_double_inplace(p.pull("v"), &evec, &n);
  printf("E=%.16f   Delta E=%.16e\n", E,E-E_LITERATURE);
  cpu_time.tick();
  verbose("Total running time for solving generalized eigenvalue problem: %g s", cpu_time.accumulated());
  // Convert double array to vector and coefficient vector into a Solution.
  for (int i=0;i<ndof;i++) eivec->set(i,abs(evec[i]));// ensure that the plotted eigenvector is positive
  Solution* sln = new Solution(&space, eivec);
  ScalarView view("Solution", new WinGeom(0, 0, 1024, 768));
  /*
  // this code prints the values of the solution on a grid 
  // of points along the z axis.
  double xp;
  for (int i=0;i<4000;i++){
    xp=i*0.001;
    printf("%24.15e    %24.15e\n",xp,sln->get_pt_value(0.0,xp,H2D_FN_VAL_0));
  }
  */
  // Visualize the solution.
  //view.set_3d_mode();
  view.show(sln);
  // Wait for the view to be closed.
  View::wait(H2DV_WAIT_KEYPRESS);
  //  }  
  return 0; 
};

