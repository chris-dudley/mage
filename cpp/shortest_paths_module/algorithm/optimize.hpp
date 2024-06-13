#pragma once

#include <vector>
#include <optional>
#include <functional>

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
  double ftol;
  double ptol;
  double xtol;
  uint max_iter;
};

template <typename TSize = std::uint64_t>
class FrankWolfe {
public:
  using TargetFunc = std::function<double(const Eigen::VectorXd&)>;

  FrankWolfe(const size_t num_variables, TargetFunc target_function, Eigen::VectorXd initial_guess,
             const OptimizationOptions& optimization_options);

  /// @brief Updates the weights given a gradient vector
  /// @param gradient Gradient vector
  inline void update_costs(Eigen::VectorXd &gradient);

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
  inline Eigen::VectorXd minimize_dot(Eigen::VectorXd &gradient);

  /// @brief Determines the gradient by shifting each value of `allocation` individually up and down by `eps`
  /// and evaulating the target function, then taking the average.
  /// Used if no separate gradient function is defined.
  inline Eigen::VectorXd evaluate_gradient(const Eigen::VectorXd& allocation, const double eps = 1.0e-6) const;

  /// @brief Determines the step size to take between two candidate
  /// solutions using a back-stepping proceedure.
  ///
  /// @param current The current solution
  /// @param candidate The candidate solution
  inline double determine_step(Eigen::VectorXd &candidate);

  inline void check_exit_flags();

  bool is_valid() const noexcept {
    return this->status == HighsStatus::kOk;
  }

  /// @brief Runs the Frank-Wolfe algorithm
  bool run();

  /// @brief Returns the run result
  const Eigen::VectorXd& result() const;

  /// @brief Returns the expected outout for the computed solution input.
  double expected_output() const;

private:
  TargetFunc target_func;
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

};

template <typename TSize>
FrankWolfe<TSize>::FrankWolfe(
  const size_t num_variables, TargetFunc target_function,
  Eigen::VectorXd initial_guess, const OptimizationOptions& optimization_options
)
    : target_func(target_function),
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
  this->fx = target_function(this->solution);
  this->fx_old = this->fx;

  // Create and populate a HighsModel instance for the LP
  // Starts by creating the following problem, which can
  // then be updated to implement a Frank-Wolfe method
  //
  // Min    f  =  sum_i 0 * x[i]
  // s.t.
  //           sum(x) = 1

  std::vector<double> costs(num_variables, 0);
  std::vector<double> col_lower(num_variables, 0);
  std::vector<double> col_upper(num_variables, 1);

  std::vector<HighsInt> starts(num_variables + 1, 0);
  for (uint i = 0; i <= num_variables; i++) starts[i] = i;

  std::vector<HighsInt> index(num_variables, 0);
  std::vector<double> values(num_variables, 1.0);

  this->model.lp_.num_col_ = num_variables;
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
  if (this->status != HighsStatus::kOk) {
    throw std::runtime_error("Highs model was not correctly specified");
  }
}

template <typename TSize>
inline void FrankWolfe<TSize>::update_costs(Eigen::VectorXd &gradient) {
  assert(gradient.size() == this->ptr_lp.num_col_ && "Number of variables does not match the size of the gradient");
  this->status = highs.changeColsCost(0, gradient.size()-1, gradient.data());
}

template <typename TSize>
inline Eigen::VectorXd FrankWolfe<TSize>::minimize_dot(Eigen::VectorXd &gradient) {
  this->update_costs(gradient);

  this->status = this->highs.run();
  if (this->status != HighsStatus::kOk) {
    throw std::runtime_error("Highs solver broke during the run step");
  }

  const HighsSolution &ptr_h_sol = highs.getSolution();

  Eigen::VectorXd solution(this->ptr_lp.num_col_);

  for (uint i = 0; i < solution.size(); i++) {
    solution[i] = ptr_h_sol.col_value[i];
  }

  return solution;
}

template <typename TSize>
inline Eigen::VectorXd FrankWolfe<TSize>::evaluate_gradient(const Eigen::VectorXd& allocation, const double eps) const {
  Eigen::VectorXd gradient(allocation.size());
  Eigen::VectorXd shift_up = allocation;
  Eigen::VectorXd shift_down = allocation;

  for (auto i = 0; i < gradient.size(); i++) {
    // Shift for gradient
    shift_up[i] += eps;
    shift_down[i] -= eps;

    gradient[i] = (target_func(shift_up) - target_func(shift_down)) / (2 * eps);

    // Shift back
    shift_up[i] = allocation[i];
    shift_down[i] = allocation[i];
  }

  return gradient;
}

template <typename TSize>
inline double FrankWolfe<TSize>::determine_step(Eigen::VectorXd &candidate) {
  Eigen::VectorXd direction = candidate - this->solution;

  double gamma = 1;
  for (uint i = 0; i < 10; i++) {
    double temp = this->target_func(this->solution + gamma * direction);
    if (temp < this->fx) break;
    gamma /= 2.0;
  }

  return gamma;
}

template <typename TSize>
inline void FrankWolfe<TSize>::check_exit_flags() {
  // Check that there was a meaningful change in the output
  if (abs(fx - fx_old) < options.ftol and options.ftol > 0) ftol_flag = true;

  // Check that there was a meaningful percent change in the output
  if (options.ptol > 0 && fx == 0.0 && fx != fx_old) {
    ptol_flag = ((fx_old - fx) / fx_old) < options.ptol;
  } else if (options.ptol > 0 && fx != 0.0) {
    ptol_flag = ((fx - fx_old) / fx) < options.ptol;
  }

  // Check that there was a meaningful change in the allocation
  if (abs((solution - solution_old).norm()) < options.xtol and options.xtol > 0) xtol_flag = true;
}

template <typename TSize>
bool FrankWolfe<TSize>::run() {
  for (uint i = 0; i < options.max_iter; i++) {
    auto gradient = this->evaluate_gradient(this->solution);
    if (this->status != HighsStatus::kOk)
      break;
    auto candidate = this->minimize_dot(gradient);
    if (this->status != HighsStatus::kOk)
      break;

    Eigen::VectorXd direction = candidate - this->solution;
    double gamma = this->determine_step(candidate);

    double temp = this->target_func(this->solution + gamma * direction);
    if (this->fx < temp) {
      this->solution_old = this->solution;
      this->fx_old = this->fx;

      this->solution = this->solution + gamma * direction;
      this->fx = temp;
    }

    this->check_exit_flags();
    if (ftol_flag | ptol_flag | xtol_flag) break;
  }

  return is_valid();
}

template <typename TSize>
const Eigen::VectorXd& FrankWolfe<TSize>::result() const {
  return solution;
}

template <typename TSize>
double FrankWolfe<TSize>::expected_output() const {
  return fx;
}

template <typename TSize = std::uint64_t>
EdgeNetwork<TSize> PathsToEdgeNetwork(
  const std::vector<Path<TSize>>& paths,
  const std::vector<std::vector<double>>& x_coords, std::vector<std::vector<double>>& y_coords
) {
  std::vector<std::vector<PolyEdge<TSize>>> edges;

  for (const auto& path : paths) {
    std::vector<PolyEdge<TSize>> poly_path;
    for (auto edge_id: path.edges) {
      poly_path.emplace_back(edge_id, x_coords.at(edge_id), y_coords.at(edge_id));
    }
    edges.emplace_back(std::move(poly_path));
  }

  return EdgeNetwork(std::move(edges));
}

template <typename TSize = std::uint64_t>
OptimizationResult OptimizeFlows(
  const mg_graph::GraphView<TSize> &graph, const std::vector<Path<TSize>> &paths,
  const std::vector<std::vector<double>>& x_coords, std::vector<std::vector<double>>& y_coords,
  const double input_flow,
  const std::optional<std::vector<double>>& initial_allocations = std::nullopt,
  const bool initial_allocation_ratios = false,
  const double ftol = 1.0e-6, const double ptol = 1.0e-6, const double xtol = 1.0e-4, const uint max_iter = 50
) {
  OptimizationResult result{
    std::vector<double>(paths.size(), 0.0),
    std::vector<double>(paths.size(), 0.0),
    0.0,
    false
  };
  OptimizationOptions options{ftol, ptol, xtol, max_iter};
  if (paths.size() == 0) {
    // Nothing to optimize over.
    result.success = true;
    return result;
  }

  // Set up initial guesses
  Eigen::VectorXd allocations(paths.size());
  if (initial_allocations) {
    if (initial_allocations->size() != paths.size()) {
      throw std::domain_error(fmt::format("Initial allocations ({}) must match number of paths ({})", initial_allocations->size(), paths.size()));
    }
    Eigen::Map<const Eigen::VectorXd> init_alloc_vec(initial_allocations->data(), initial_allocations->size());

    if (initial_allocation_ratios) {
      // Initial allocations are value amounts, make sure they add up to the input flow
      double alloc_diff = abs(init_alloc_vec.sum() - 1.0); 
      if (alloc_diff > xtol) {
        throw std::domain_error(fmt::format("Sum of initial ratio allocations ({}) must be 1.0", init_alloc_vec.sum()));
      }
      allocations = init_alloc_vec;
    } else {
      // Initial allocations are value amounts, make sure they add up to the input flow
      double alloc_diff = abs(init_alloc_vec.sum() - input_flow); 
      if (alloc_diff > xtol) {
        throw std::domain_error(fmt::format("Sum of initial allocations ({}) must match input flow ({})", init_alloc_vec.sum(), input_flow));
      }
      allocations = init_alloc_vec / input_flow;
    }
  } else {
    // Default to even spread
    allocations.setConstant(1.0 / paths.size());
  }

  EdgeNetwork<TSize> edge_network = PathsToEdgeNetwork(paths, x_coords, y_coords);
  auto target_func = [&edge_network, input_flow](const Eigen::VectorXd& allocs) -> double {
    return edge_network.evaluate_allocations(allocs * input_flow);
  };

  FrankWolfe<TSize> optimizer(paths.size(), target_func, allocations, options);

  result.success = optimizer.run();
  if (result.success) {
    const auto& solution = optimizer.result();
    if (solution.size() > paths.size()) {
      throw std::runtime_error(fmt::format(
        "Optimizer returned more solutions ({}) than we have paths ({})!",
        solution.size(),
        paths.size()
      ));
    }
    for (auto i = 0; i < solution.size(); i++) {
      result.input_ratios[i] = solution[i];
      result.input_amounts[i] = solution[i] * input_flow;
      result.expected_output = optimizer.expected_output();
    }
  }
  
  return result;
}

}  // namespace shortest_paths
