#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"
#include "weakform_library/laplace.h"

class CustomWeakFormPoisson : public WeakForm
{
public:
  CustomWeakFormPoisson(double const_f) : WeakForm(1)
  {
    add_matrix_form(new DefaultMatrixFormStiffness(0, 0));
    add_vector_form(new DefaultVectorFormVolConst(0, const_f));
  };
};

