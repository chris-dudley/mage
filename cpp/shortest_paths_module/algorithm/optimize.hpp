#pragma once

#include <vector>

#include <mg_graph.hpp>
#include "Highs.h"
#include "libs/Eigen/Core"

#include "path.hpp"
#include "polyedge.hpp"

namespace shortest_paths {

struct OptimizationResult {
  std::vector<double> input_ratios;
  std::vector<double> input_amounts;
  double expected_output;
  bool success;
};

struct OptimizationOptions {
  OptimizationOptions();

  double ftol;
  double ptol;
  double xtol;
  uint max_iter;
};

OptimizationOptions::OptimizationOptions() : ftol(1e-6), ptol(1e-6), xtol(1e-4), max_iter(50) {}

// TODO: Placeholder interface
template <typename TSize = std::uint64_t>
OptimizationResult OptimizeFlows(const mg_graph::GraphView<TSize> &graph, const std::vector<Path<TSize>> &paths) {
  OptimizationResult result{std::vector<double>(paths.size(), 0.0), std::vector<double>(paths.size(), 0.0), 0.0, false};

  // TODO: Actually implement
  result.success = true;

  return result;
}

template <typename TSize = std::uint64_t>
class FrankWolfe {
  FrankWolfe(const uint size, const EdgeNetwork<TSize> &ptr_edge_network, Eigen::VectorXd initial_guess,
             OptimizationOptions = OptimizationOptions());

  const EdgeNetwork<TSize> &ptr_edge_network;
  OptimizationOptions options;

  HighsModel model;
  Highs highs;
  HighsStatus status;
  const HighsLp &ptr_lp;

  Eigen::VectorXd solution;
  Eigen::VectorXd solution_old;

  double fx;
  double fx_old;

  bool ftol_flag;
  bool ptol_flag;
  bool xtol_flag;

  /// @brief Updates the weights given a gradient vector
  /// @param gradient Gradient vector
  void update_costs(Eigen::VectorXd &gradient);

  /// @brief Finds a vector on the unit simplex (probability
  /// vectors) which maximizes dot product between it and the
  /// gradient. Recalling that:
  ///
  ///   g \cdot v = ||g|| ||v|| \cos(\theta)
  ///
  /// This means that v is parallel and in the same direction as
  /// the gradient and v is as large as possible. This can be
  /// thought of as the furthest we can move along the direction of
  /// the gradient before leaving the unit simplex.
  ///
  /// @param gradient Gradient vector
  Eigen::VectorXd minimize_dot(Eigen::VectorXd &gradient);

  /// @brief Determines the step size to take between two candidate
  /// solutions using a back-stepping proceedure.
  ///
  /// @param current The current solution
  /// @param candidate The candidate solution
  double determine_step(Eigen::VectorXd &candidate);

  bool check_exit_flags();

  /// @brief Runs the Frank-Wolfe algorithm
  void run();

  /// @brief Returns the run result
  OptimizationResult result();
};

template <typename TSize>
FrankWolfe<TSize>::FrankWolfe(const uint size, const EdgeNetwork<TSize> &ptr_edge_network,
                              Eigen::VectorXd initial_guess, OptimizationOptions optimization_options)
    : ptr_edge_network(ptr_edge_network),
      options(optimization_options),
      model(),
      highs(),
      ptr_lp(highs.getLp()),
      ftol_flag(false),
      ptol_flag(false),
      xtol_flag(false) {
  // Initialize values
  this->solution = initial_guess;
  this->solution_old = initial_guess;
  this->fx = ptr_edge_network->evaluate_allocations(this->solution);
  this->fx_old = this->fx;

  // Create and populate a HighsModel instance for the LP
  // Starts by creating the following problem, which can
  // then be updated to implement a Frank-Wolfe method
  //
  // Min    f  =  sum_i 0 * x[i]
  // s.t.
  //           sum(x) = 1

  std::vector<double> costs(size, 0);
  std::vector<double> col_lower(size, 0);
  std::vector<double> col_upper(size, 1);

  std::vector<HighsInt> starts(size + 1, 0);
  for (uint i = 0; i <= size; i++) starts[i] = i;

  std::vector<HighsInt> index(size, 0);
  std::vector<double> values(size, 1.0);

  this->model.lp_.num_col_ = size;
  this->model.lp_.num_row_ = 1;

  this->model.lp_.sense_ = ObjSense::kMaximize;
  this->model.lp_.col_cost_ = costs;
  this->model.lp_.col_lower_ = col_lower;
  this->model.lp_.col_upper_ = col_upper;

  this->model.lp_.row_lower_ = {1.0};
  this->model.lp_.row_upper_ = {1.0};

  this->model.lp_.a_matrix_.format_ = MatrixFormat::kColwise;
  this->model.lp_.a_matrix_.start_ = starts;
  this->model.lp_.a_matrix_.index_ = index;
  this->model.lp_.a_matrix_.value_ = values;

  this->status = highs.passModel(model);
  assert(this->status == HighsStatus::kOk && "God hates us, and the model was specified incorrectly.");
}

template <typename TSize>
void FrankWolfe<TSize>::update_costs(Eigen::VectorXd &gradient) {
  assert(gradient.size() == this->ptr_lp.num_col_ && "Number of variables does not match the size of the gradient");
  for (HighsInt i = 0; i < gradient.size(); i++) {
    highs.changeColsCost(i, i, &gradient[i]);
  }
  assert(this->status == HighsStatus::kOk && "Failed to update column costs");
}

template <typename TSize>
Eigen::VectorXd FrankWolfe<TSize>::minimize_dot(Eigen::VectorXd &gradient) {
  this->update_costs(gradient);

  this->highs.run();
  assert(this->status == HighsStatus::kOk && "Solver broke during run step");

  const HighsSolution &ptr_h_sol = highs.getSolution();

  Eigen::VectorXd solution(this->ptr_lp.num_col_);

  for (uint i = 0; i < solution.size(); i++) {
    solution[i] = ptr_h_sol.col_value[i];
  }

  return solution;
}

template <typename TSize>
double FrankWolfe<TSize>::determine_step(Eigen::VectorXd &candidate) {
  Eigen::VectorXd direction = candidate - this->solution;

  double gamma = 1;
  for (uint i = 0; i < 10; i++) {
    double temp = this->ptr_edge_network->evaluate_allocations(this->solution + gamma * direction);
    if (temp < this->fx) break;
    gamma /= 2.0;
  }

  return gamma;
}

template <typename TSize>
bool FrankWolfe<TSize>::check_exit_flags() {
  // Check that there was a meaningful change in the output
  if (abs(this->fx - this->fx_old) < options.ftol and options.ftol > 0) this->ftol_flag = true;

  // Check that there was a meaningful percent change in the output
  if (abs((this->fx - this->fx_old) / this->fx) < options.ptol and options.ptol > 0) this->ptol_flag = true;

  // Check that there was a meaningful change in the allocation
  if (abs((this->solution - this->solution_old).norm()) < options.xtol and options.xtol > 0) this->xtol_flag = true;
}

template <typename TSize>
void FrankWolfe<TSize>::run() {
  for (uint i = 0; i < options.max_iter; i++) {
    this->gradient = ptr_edge_network.evaluate_gradient(this->solution);
    this->candidate = this->minimize_dot(this->gradient);
    Eigen::VectorXd direction = this->candidate - this->solution;
    double gamma = this->determine_step(this->candidate);

    double temp = this->evaluate_allocations(this->solution + gamma * direction);
    if (this->fx < temp) {
      this->solution_old = this->solution;
      this->fx_old = this->fx;

      this->solution = this->solution + gamma * direction;
      this->fx = temp;
    }

    this->check_exit_flags();
    if (ftol_flag | ptol_flag | xtol_flag) break;
  }
}

}  // namespace shortest_paths
