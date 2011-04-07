#define HERMES_REPORT_ALL
#include "hermes2d.h"
#include <stdio.h>

using Teuchos::RCP;
using Teuchos::rcp;
using Hermes::EigenSolver;

//  This example solves the eigenproblem for the Laplace operator in 
//  a square with zero boundary conditions. Python and Pysparse must
//  be installed. 
//
//  PDE: (-Laplace +x*x+y*y) u = lambda_k u,
//  where lambda_0, lambda_1, ... are the eigenvalues.
//
//  Domain: Square (-4, 4)^2.
//
//  BC:  Homogeneous Dirichlet.
//
//  The following parameters can be changed:

const int NUMBER_OF_EIGENVALUES = 1;             // Desired number of eigenvalues.
const int P_INIT = 4;                             // Uniform polynomial degree of mesh elements.
const int INIT_REF_NUM = 3;                       // Number of initial mesh refinements.
const int REF_TOWARDS_ORIGIN=1;                   // Number of mesh refinements towards origin.
const double TARGET_VALUE = 2.0;                  // PySparse parameter: Eigenvalues in the vicinity of 
                                                  // this number will be computed. 
const double TOL = 1e-10;                         // Pysparse parameter: Error tolerance.
const int MAX_ITER = 1000;                        // PySparse parameter: Maximum number of iterations.

// Boundary markers.
const std::string BDY = "Zero Dirichlet";

// Weak forms.
#include "definitions.cpp"

int main(int argc, char* argv[])
{
  if (argc <2)
    error("Not enough parameters, provide a test number!");
  int test_type=atoi(argv[1]);
  info("Desired number of eigenvalues: %d.", NUMBER_OF_EIGENVALUES);

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &mesh);

  // Perform initial mesh refinements (optional).
  for (int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();
  //The following code finds the 
  // vertex id at the origin (0.0,0.0).
  int NUM_ELEMS=mesh.get_max_element_id();
  info("NUM_ELEMS = %d",NUM_ELEMS);
  int elem_ids_to_be_refined[4];
  int vertex_id_origin=-1;
  int j=0;
  for (int id=0;id<NUM_ELEMS;id++) {
    Element* elem=mesh.get_element(id);
    //info("active=%d",elem->active);
    if (elem->active) {
      for(int i=0;i<4;i++){
	Node* node= elem->vn[i];
	if ( node->x == 0.0 && node->y == 0.0){ 
	  vertex_id_origin=node->id;
	  info("id= %d",id);
	  elem_ids_to_be_refined[j]=id;
	  j++;
	}
      }
    }
  }
  if (vertex_id_origin == -1) error("This should not happen!.");
  info("Vertex_id_origin = %d",vertex_id_origin);
  if (test_type == 0) {
      info("refining towards vertex using builtin method");
      mesh.refine_towards_vertex(vertex_id_origin,REF_TOWARDS_ORIGIN);
  }else if (test_type ==1){
    info("refining towards vertex via 4 calls to refine_element_id  method");
    for (int i=0;i<4;i++) {
      info("%d %d\n",i,elem_ids_to_be_refined[i]);
      mesh.refine_element_id(elem_ids_to_be_refined[i],0);
    }
  } else {
    error("illegal value of test_type was specified");
  }
  // Initialize boundary conditions. 
  DefaultEssentialBCConst bc_essential(BDY, 0.0);
  EssentialBCs bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, &bcs, P_INIT);
  int ndof = space.get_num_dofs();
  info("ndof: %d.", ndof);

  // Initialize the weak formulation.
  WeakFormEigenLeft wf_left;
  WeakFormEigenRight wf_right;

  // Initialize matrices.
  RCP<SparseMatrix> matrix_left = rcp(new CSCMatrix());
  RCP<SparseMatrix> matrix_right = rcp(new CSCMatrix());

  // Assemble the matrices.
  bool is_linear = true;
  DiscreteProblem dp_left(&wf_left, &space, is_linear);
  dp_left.assemble(matrix_left.get());
  DiscreteProblem dp_right(&wf_right, &space, is_linear);
  dp_right.assemble(matrix_right.get());

  EigenSolver es(matrix_left, matrix_right);
  info("Calling Pysparse...");
  es.solve(NUMBER_OF_EIGENVALUES, TARGET_VALUE, TOL, MAX_ITER);
  info("Pysparse finished.");
  es.print_eigenvalues();

  // Initializing solution vector, solution and ScalarView.
  double* coeff_vec;
  Solution sln;
  // ScalarView view("Solution", new WinGeom(0, 0, 440, 350));

  // Reading solution vectors and visualizing.
  double* eigenval = new double[NUMBER_OF_EIGENVALUES];
  int neig = es.get_n_eigs();
  if (neig != NUMBER_OF_EIGENVALUES) error("Mismatched number of eigenvectors in the eigensolver output file.");  
  for (int ieig = 0; ieig < neig; ieig++) {
    eigenval[ieig] = es.get_eigenvalue(ieig);
    info("eigenvalue [%d]=%24.15f",ieig,eigenval[ieig]);
    int n;
    es.get_eigenvector(ieig, &coeff_vec, &n);
    // Convert coefficient vector into a Solution.
    Solution::vector_to_solution(coeff_vec, &space, &sln);

    // Visualize the solution.
    char title[100];
    sprintf(title, "Solution %d, val = %g", ieig, eigenval[ieig]);
    //view.set_title(title);
    //view.show(&sln);

    // Wait for keypress.
    //View::wait(HERMES_WAIT_KEYPRESS);
  }

  return 0; 
};

