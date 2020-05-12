// Deal.ii
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/convergence_table.h>
#include <deal.II/base/numbers.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/affine_constraints.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

// STL
#include <array>
#include <fstream>
#include <iostream>
#include <cmath>


namespace Helmholtz
{
  using namespace dealii;

  ///////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////

  template <int dim>
  class Solution : public Function<dim>
  {
  public:
    virtual double value(const Point<dim> & p,
                         const unsigned int component = 0) const override;

    virtual Tensor<1, dim>
    gradient(const Point<dim> & p,
             const unsigned int component = 0) const override;
  };


  template <int dim>
  double Solution<dim>::value(const Point<dim> &p, const unsigned int) const
  {
    double return_value = 1;

    for (unsigned int i = 0; i<dim; ++i)
        	return_value *= cos(2*numbers::PI*p(i));

    return return_value;
  }


  template <>
  Tensor<1, 2> Solution<2>::gradient(const Point<2> &p,
                                         const unsigned int) const
  {
    Tensor<1, 2> return_value;

	return_value[0] = - 2 * numbers::PI * sin(2*numbers::PI*p(0)) * cos(2*numbers::PI*p(1));
	return_value[1] = - 2 * numbers::PI * cos(2*numbers::PI*p(0)) * sin(2*numbers::PI*p(1));

    return return_value;
  }


  template <>
    Tensor<1, 3> Solution<3>::gradient(const Point<3> &p,
                                           const unsigned int) const
    {
      Tensor<1, 3> return_value;

      return_value[0] = - 2 * numbers::PI * sin(2*numbers::PI*p(0)) * cos(2*numbers::PI*p(1)) * cos(2*numbers::PI*p(2));
      return_value[1] = - 2 * numbers::PI * cos(2*numbers::PI*p(0)) * sin(2*numbers::PI*p(1)) * cos(2*numbers::PI*p(2));
      return_value[2] = - 2 * numbers::PI * cos(2*numbers::PI*p(0)) * cos(2*numbers::PI*p(1)) * sin(2*numbers::PI*p(2));

      return return_value;
    }


  ///////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////


  template <int dim>
  class RightHandSide : public Function<dim>
  {
  public:
    virtual double value(const Point<dim> & p,
                         const unsigned int component = 0) const override;
  };


  template <int dim>
  double RightHandSide<dim>::value(const Point<dim> &p,
                                   const unsigned int) const
  {
    double return_value = (1 + 4 * dim * std::pow(numbers::PI,2));

    for (unsigned int i = 0; i<dim; ++i)
    	return_value *= cos(2*numbers::PI*p(i));

    return return_value;
  }


  ///////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////


  template <int dim>
  class HelmholtzSolver
  {
  public:
    HelmholtzSolver(const FiniteElement<dim> &fe);

    ~HelmholtzSolver();

    void run();

  private:
    void setup_system();
    void assemble_system();
    void solve();
    void process_solution(const unsigned int cycle);
    void output_vtk(const unsigned int cycle);

    Triangulation<dim> triangulation;
    DoFHandler<dim>    dof_handler;

    SmartPointer<const FiniteElement<dim>> fe;

    AffineConstraints<double> constraints;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double> solution;
    Vector<double> system_rhs;

    ConvergenceTable convergence_table;
  };


  template <int dim>
  HelmholtzSolver<dim>::HelmholtzSolver(const FiniteElement<dim> &fe)
    : dof_handler(triangulation)
    , fe(&fe)
  {}


  template <int dim>
  HelmholtzSolver<dim>::~HelmholtzSolver()
  {
    dof_handler.clear();
  }


  template <int dim>
  void HelmholtzSolver<dim>::setup_system()
  {
    dof_handler.distribute_dofs(*fe);
    DoFRenumbering::Cuthill_McKee(dof_handler);

    constraints.clear();
    DoFTools::make_hanging_node_constraints(dof_handler,
                                            constraints);

    Solution<dim>	exact_solution;

    VectorTools::interpolate_boundary_values(dof_handler,
                                               0,
                                               exact_solution,
                                               constraints);

    constraints.close();

    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp);
    constraints.condense(dsp);
    sparsity_pattern.copy_from(dsp);

    system_matrix.reinit(sparsity_pattern);

    solution.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());
  }


  template <int dim>
  void HelmholtzSolver<dim>::assemble_system()
  {
    QGauss<dim>     quadrature_formula(fe->degree + 1);

    const unsigned int n_q_points      = quadrature_formula.size();

    const unsigned int dofs_per_cell = fe->dofs_per_cell;

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>     cell_rhs(dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    FEValues<dim> fe_values(*fe,
                            quadrature_formula,
                            update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);

    RightHandSide<dim>  right_hand_side;
    std::vector<double> rhs_values(n_q_points);

    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        cell_matrix = 0.;
        cell_rhs    = 0.;

        fe_values.reinit(cell);

        right_hand_side.value_list(fe_values.get_quadrature_points(),
                                   rhs_values);

        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
              for (unsigned int j = 0; j < dofs_per_cell; ++j)
                cell_matrix(i, j) +=
                  ((fe_values.shape_grad(i, q_point) *     // grad phi_i(x_q)
                      fe_values.shape_grad(j, q_point)     // grad phi_j(x_q)
                    +                                      //
                    fe_values.shape_value(i, q_point) *    // phi_i(x_q)
                      fe_values.shape_value(j, q_point)) * // phi_j(x_q)
                   fe_values.JxW(q_point));                // dx


              cell_rhs(i) += (fe_values.shape_value(i, q_point) * // phi_i(x_q)
                              rhs_values[q_point] *               // f(x_q)
                              fe_values.JxW(q_point));            // dx
            }

        cell->get_dof_indices(local_dof_indices);
        constraints.distribute_local_to_global(
                cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);
      }
  }


  template <int dim>
  void HelmholtzSolver<dim>::solve()
  {
    SolverControl solver_control(1000, 1e-12);
    SolverCG<>    cg(solver_control);

    PreconditionSSOR<> preconditioner;
    preconditioner.initialize(system_matrix, 1.2);

    cg.solve(system_matrix, solution, system_rhs, preconditioner);

    std::cout << "   " << solver_control.last_step()
                << " CG iterations needed to obtain convergence." << std::endl;

    constraints.distribute(solution);
  }


  template <int dim>
  void HelmholtzSolver<dim>::process_solution(const unsigned int cycle)
  {
	Solution<dim> exact_solution;

    Vector<float> difference_per_cell(triangulation.n_active_cells());
    VectorTools::integrate_difference(dof_handler,
                                      solution,
                                      exact_solution,
                                      difference_per_cell,
                                      QGauss<dim>(fe->degree + 1),
                                      VectorTools::L2_norm);
    const double L2_error =
      VectorTools::compute_global_error(triangulation,
                                        difference_per_cell,
                                        VectorTools::L2_norm);


    VectorTools::integrate_difference(dof_handler,
                                      solution,
                                      Solution<dim>(),
                                      difference_per_cell,
                                      QGauss<dim>(fe->degree + 1),
                                      VectorTools::H1_seminorm);
    const double H1_error =
      VectorTools::compute_global_error(triangulation,
                                        difference_per_cell,
                                        VectorTools::H1_seminorm);


    const unsigned int n_active_cells = triangulation.n_active_cells();
    const unsigned int n_dofs         = dof_handler.n_dofs();

    std::cout << "Cycle " << cycle << ':' << std::endl
              << "   Number of active cells:       " << n_active_cells
              << std::endl
              << "   Number of degrees of freedom: " << n_dofs << std::endl;

    convergence_table.add_value("cycle", cycle);
    convergence_table.add_value("cells", n_active_cells);
    convergence_table.add_value("dofs", n_dofs);
    convergence_table.add_value("L2", L2_error);
    convergence_table.add_value("H1", H1_error);
  }


  template <int dim>
  void HelmholtzSolver<dim>::output_vtk(const unsigned int cycle)
  {
	std::string vtk_filename
		= "solution-" + Utilities::int_to_string(dim,1) + "D";

	switch (fe->degree)
	  {
		case 1:
		  vtk_filename += "_q1";
		  break;
		case 2:
		  vtk_filename += "_q2";
		  break;
		default:
		  Assert(false, ExcNotImplemented());
	  }

	vtk_filename += "_cycle-" + Utilities::int_to_string(cycle,1);
	vtk_filename += ".vtk";

	std::ofstream output(vtk_filename);

	DataOut<dim> data_out;
	data_out.attach_dof_handler(dof_handler);
	data_out.add_data_vector(solution, "solution");

	data_out.build_patches(fe->degree);
	data_out.write_vtk(output);
}


  template <int dim>
  void HelmholtzSolver<dim>::run()
  {

    const unsigned int n_cycles = 5;

    for (unsigned int cycle = 0; cycle < n_cycles; ++cycle)
      {
        if (cycle == 0)
          {
            GridGenerator::hyper_cube(triangulation, 0, 1);
            triangulation.refine_global(1);
          }
        else
        	triangulation.refine_global();

        setup_system();

        assemble_system();
        solve();

        process_solution(cycle);

        output_vtk(cycle);
      }


    convergence_table.set_precision("L2", 3);
    convergence_table.set_precision("H1", 3);

    convergence_table.set_scientific("L2", true);
    convergence_table.set_scientific("H1", true);

    convergence_table.set_tex_caption("cells", "\\# cells");
    convergence_table.set_tex_caption("dofs", "\\# dofs");
    convergence_table.set_tex_caption("L2", "$L^2$-error");
    convergence_table.set_tex_caption("H1", "$H^1$-error");

    convergence_table.set_tex_format("cells", "r");
    convergence_table.set_tex_format("dofs", "r");

    std::cout << std::endl;
    convergence_table.write_text(std::cout);

    // Write table into a LaTeX file.
    std::string error_filename
		= "error-" + Utilities::int_to_string(dim,1) + "D";

    switch (fe->degree)
      {
        case 1:
          error_filename += "_q1";
          break;
        case 2:
          error_filename += "_q2";
          break;
        default:
          Assert(false, ExcNotImplemented());
      }

    error_filename += ".tex";
    std::ofstream error_table_file(error_filename);

    convergence_table.write_tex(error_table_file);

    /////////////////

	convergence_table.add_column_to_supercolumn("cycle", "n cells");
	convergence_table.add_column_to_supercolumn("cells", "n cells");

	std::vector<std::string> new_order;
	new_order.emplace_back("n cells");
	new_order.emplace_back("H1");
	new_order.emplace_back("L2");
	convergence_table.set_column_order(new_order);

	convergence_table.evaluate_convergence_rates(
	  "L2", ConvergenceTable::reduction_rate);
	convergence_table.evaluate_convergence_rates(
	  "L2", ConvergenceTable::reduction_rate_log2);
	convergence_table.evaluate_convergence_rates(
	  "H1", ConvergenceTable::reduction_rate);
	convergence_table.evaluate_convergence_rates(
	  "H1", ConvergenceTable::reduction_rate_log2);

	std::cout << std::endl;
	convergence_table.write_text(std::cout);

	std::string conv_filename
		= "convergence-" + Utilities::int_to_string(dim,1) + "D";

	switch (fe->degree)
	  {
		case 1:
		  conv_filename += "_q1";
		  break;
		case 2:
		  conv_filename += "_q2";
		  break;
		default:
		  Assert(false, ExcNotImplemented());
	  }
	conv_filename += ".tex";

	std::ofstream table_file(conv_filename);
	convergence_table.write_tex(table_file);
  }

} // namespace Helmholtz


  ///////////////////////////////////////////
  ///////////////////////////////////////////
  ///////////////////////////////////////////


int main()
{
  const unsigned int dim = 3;

  try
    {
      using namespace dealii;

      {
        std::cout << "Solving with Q1 elements in " << dim << "D" << std::endl
                  << std::endl
                  << "=============================="
                  << std::endl
                  << std::endl;

        FE_Q<dim>             fe(1);
        Helmholtz::HelmholtzSolver<dim> helmholtz_problem_2d(fe);

        helmholtz_problem_2d.run();

        std::cout << std::endl;
      }

      {
        std::cout << "Solving with Q2 elements in " << dim << "D" << std::endl
                  << "==============================" << std::endl
                  << std::endl;

        FE_Q<dim>             fe(2);
        Helmholtz::HelmholtzSolver<dim> helmholtz_problem_2d(fe);

        helmholtz_problem_2d.run();

        std::cout << std::endl;
      }
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
