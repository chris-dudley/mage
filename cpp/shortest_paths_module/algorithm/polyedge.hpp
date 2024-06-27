#pragma once

#include "fmt/core.h"
#include "fmt/ranges.h"
#include "libs/Eigen/Core"
#include "libs/unsupported/Eigen/NNLS"

#include <cmath>
#include <limits>
#include <ranges>
#include <stdexcept>
#include <unordered_map>
#include <utility>

namespace shortest_paths {

template <typename TSize = std::uint64_t>
class PolyEdge {
 public:
  static constexpr const TSize INVALID_ID = std::numeric_limits<TSize>::max();

  using SelfType = PolyEdge<TSize>;

  // Default constructor with an invalid edge id and single coefficient of 1.0
  PolyEdge() : id_(INVALID_ID), coefficients_{{1.0}} {}
  PolyEdge(TSize id) : id_(id), coefficients_{{1.0}} {}
  // Basic constructors from parts
  PolyEdge(TSize id, const Eigen::VectorXd &coefficients) : id_(id), coefficients_(coefficients) {}
  PolyEdge(TSize id, Eigen::VectorXd &&coefficients) : id_(id), coefficients_(std::move(coefficients)) {}

  // Copy and move constructors
  PolyEdge(const SelfType &other) : id_(other.id_), coefficients_(other.coefficients_) {}
  PolyEdge(SelfType &&other) : id_(std::move(other.id_)), coefficients_(std::move(other.coefficients_)) {}

  /// @brief Constructs a PolyEdge from an edge ID and set of points upon which to fit the polynomial.
  /// @param id ID of the edge.
  /// @param x Vector of x coordinates to which to fit the polynomial.
  /// @param y Vector of y coordinates to which to fit the polynomial.
  /// @param fit_2point_quadratic Whether to attempt to fit 2 points to f(x) = a + bx^2 instead of f(x) = a + bx
  PolyEdge(TSize id, const Eigen::VectorXd &x, const Eigen::VectorXd &y, bool fit_2point_quadratic = false)
      : id_(id), coefficients_() {
    if (x.size() != y.size()) {
      throw std::domain_error(fmt::format("Number of x and y values must match ({} != {})", x.size(), y.size()));
    }
    calculate_coefficients(x, y, fit_2point_quadratic);
  }

  /// @brief Constructs a PolyEdge from an edge ID and two contiguous containers of x and y values.
  /// @param id ID of the edge.
  /// @param x Contiguous container of x coordinates to which to fit the polynomial.
  /// @param y Contiguous container of y coordinates to which to fit the polynomial.
  /// @param fit_2point_quadratic Whether to attempt to fit 2 points to f(x) = a + bx^2 instead of f(x) = a + bx
  template <std::ranges::contiguous_range XContainer, std::ranges::contiguous_range YContainer>
  PolyEdge(TSize id, const XContainer &x, const YContainer &y, bool fit_2point_quadratic = false)
      : id_(id), coefficients_() {
    static_assert(std::same_as<double, std::ranges::range_value_t<XContainer>>,
                  "X values container must contain same floating point type as edge");
    static_assert(std::same_as<double, std::ranges::range_value_t<YContainer>>,
                  "Y values container must contain same floating point type as edge");
    const auto x_size = std::ranges::size(x);
    const auto y_size = std::ranges::size(y);

    if (x_size != y_size) {
      throw std::domain_error(fmt::format("Number of x and y values must match ({} != {})", x_size, y_size));
    }

    Eigen::Map<const Eigen::VectorXd> x_map(std::ranges::cdata(x), x_size);
    Eigen::Map<const Eigen::VectorXd> y_map(std::ranges::cdata(y), y_size);

    calculate_coefficients(x_map, y_map, fit_2point_quadratic);
  }

  // Copy and move assignment
  SelfType &operator=(const SelfType &other) {
    id_ = other.id_;
    coefficients_ = other.coefficients_;
    return *this;
  }
  SelfType &operator=(SelfType &&other) noexcept {
    id_ = other.id_;
    coefficients_ = std::move(other.coefficients_);
    return *this;
  }

  /// @brief Swaps the values of this edge with the other edge.
  /// @param other The edge with which to swap values.
  void swap(SelfType &other) {
    std::swap(id_, other.id_);
    std::swap(coefficients_, other.coefficients_);
  }

  // Field accessors
  TSize id() const { return id_; }
  Eigen::VectorXd coefficients() const { return coefficients_; }

  /// @brief Update the coefficients of the fitted polynomial to the given coefficients.
  /// @param coefficients The new polynomial coefficients.
  /// @return This edge.
  SelfType &set_coefficients(const Eigen::VectorXd &coefficients) {
    coefficients_ = coefficients;
    return *this;
  }

  /// @brief Calculates the coefficients using an NNLS solver. This ensures that the polynomial
  ///     is strictly monotonically increasing.
  /// @param x The vector of x-values for points along the curve.
  /// @param y The vector of y-values for points along the curve.
  /// @param fit_2point_quadratic Whether to attempt to fit 2 points to f(x) = a + bx^2 instead of f(x) = a + bx
  /// @throws std::logic_error if the NNLS solver somehow returns negative coefficients.
  /// @throws std::invalid_argument if the sizes of the vectors don't match or contain non-finite values.
  template <typename DerivedX, typename DerivedY>
  void calculate_coefficients(const Eigen::MatrixBase<DerivedX> &x, const Eigen::MatrixBase<DerivedY> &y,
                              bool fit_2point_quadratic = false) {
    if (x.size() < 2 || y.size() < 2 || x.size() != y.size()) {
      throw std::invalid_argument(
          fmt::format("Invalid vector sizes for edge {}: dim(x)={}, dim(y)={}, require dim(x) = dim(y) > 1", id_,
                      x.size(), y.size()));
    }

    // Initialize the matrix as required by the NNLS solver.
    Eigen::MatrixXd matrix(x.size(), y.size());
    if (x.size() > 2 or not fit_2point_quadratic) {
      for (auto i = 0; i < x.size(); i++) {
        for (auto j = 0; j < y.size(); j++) {
          matrix(i, j) = pow(x(i), j);
        }
      }
    } else {
      // If dim(x) = 2, fit to f(x) = a + bx^2 instead of f(x) = a + bx
      matrix(0, 0) = 1.0;  // pow(x(0), 0)
      matrix(0, 1) = pow(x(0), 2);
      matrix(1, 0) = 1.0;  // pow(x(0), 0)
      matrix(1, 1) = pow(x(1), 2);
    }

    Eigen::NNLS<Eigen::MatrixXd> nnls_solver(matrix);
    if (x.size() == 2 and fit_2point_quadratic) {
      const auto tmp_coef = nnls_solver.solve(y);
      coefficients_ = Eigen::VectorXd{{tmp_coef(0), 0.0, tmp_coef(1)}};
    } else {
      coefficients_ = nnls_solver.solve(y);
    }

    for (auto i = 0; i < coefficients_.size(); i++) {
      if (coefficients_[i] < 0.0) {
        throw std::logic_error(
            fmt::format("Non-negative LLS solver somehow returned negative coefficient[{}]={}", i, coefficients_[i]));
      }
    }
  }

  /// @brief Evaluates the polynomial at a given x value.
  /// @param x The x value at which to evaluate the polynomial.
  /// @return The corresponding y value of the polynomial.
  double evaluate(double x) const {
    double retval = 0.0;
    for (auto i = 0; i < coefficients_.size(); i++) {
      if (i != 0) {
        retval += coefficients_[i] * pow(x, i + 1);
      } else {
        retval += coefficients_[i];
      }
    }
    if (not std::isfinite(retval)) {
      throw std::domain_error(fmt::format("edge {}: f({}) = non-finite value {}", id_, x, retval));
    }
    return retval;
  }

  /// @brief Evaluates the inverse of the polynomial at the given x coordinate.
  /// @param x The x value at which to evaulate the polynomial.
  /// @return The reciprocal of the corresponding y value of the polynomial.
  double evaluate_inv(double x) const { return 1.0 / evaluate(x); }

  /// @brief Evaluate the integral of the polynomial at the given x value.
  /// @param x The x value at which to evaluate the integral of the polynomial.
  /// @return The y value of the integral of the polynomial.
  double evaluate_integral(double x) const {
    double retval = 0;
    for (auto i = 0; i < coefficients_.size(); i++) {
      if (i != 0) {
        retval += coefficients_[i] * pow(x, i + 1) / i;
      } else {
        retval += coefficients_[i] * pow(x, 1);
      }
    }
    if (not std::isfinite(retval)) {
      throw std::domain_error(fmt::format("edge {}: integral of f({}) = non-finite value {}", id_, x, retval));
    }
    return retval;
  }

  double evaluate_avg_price(double x, double prev_swapped = 0.0) const {
    if (x == 0.0) {
      return 0;
    }
    // Average weight (price) is the integral of our polynomial divided by the
    // amount swapped to get there.
    double integral = evaluate_integral(x + prev_swapped) - evaluate_integral(prev_swapped);
    if (not std::isfinite(integral)) {
      throw std::domain_error(
          fmt::format("edge {}: avg_price({}, {}) integral = non-finite value {}", id_, x, prev_swapped, integral));
    }
    double retval = integral / x;
    if (not std::isfinite(retval)) {
      throw std::domain_error(
          fmt::format("edge {}: avg_price({}, {}) = non-finite value {}", id_, x, prev_swapped, retval));
    }
    return integral / x;
  }

  double evaluate_swap(double x, double prev_swapped = 0.0) const {
    if (x == 0.0) {
      return 0;
    }
    double avg_price = evaluate_avg_price(x, prev_swapped);
    // A swap (Out / In) * In execuated at the "average" price paid.
    double retval = x / avg_price;
    if (not std::isfinite(retval)) {
      throw std::domain_error(fmt::format("edge {}: swap({}, {}) = non-finite value {}", id_, x, prev_swapped, retval));
    }
    fmt::println(stderr, "[{}] swap({}, {}) = {}, avg_price = {}", id_, x, prev_swapped, retval, avg_price);
    return retval;
  }

 private:
  TSize id_;
  Eigen::VectorXd coefficients_;
};

template <typename TSize>
void swap(PolyEdge<TSize> &v1, PolyEdge<TSize> &v2) {
  v1.swap(v2);
}

template <typename TSize = std::uint64_t>
class EdgeNetwork {
 public:
  using SelfType = EdgeNetwork<TSize>;
  using EdgeType = PolyEdge<TSize>;
  using Vector = Eigen::VectorXd;

  EdgeNetwork(const std::vector<std::vector<EdgeType>> &paths) : paths_(paths) {}
  EdgeNetwork(std::vector<std::vector<EdgeType>> &&paths) : paths_(std::move(paths)) {}

  EdgeNetwork(const SelfType &) = default;
  EdgeNetwork(SelfType &&) noexcept = default;

  SelfType &operator=(const SelfType &) = default;
  SelfType &operator=(SelfType &&) = default;

  double evaluate_allocations(const Vector &allocation) const {
    // Vector path_outputs = allocation;
    Vector path_outputs = Vector::Zero(allocation.size());
    // Vector temp(allocation.size());

    std::unordered_map<TSize, double> amount_previously_swapped;

    for (size_t i = 0; i < paths_.size(); i++) {
      // temp.setZero();
      double in = allocation(i);
      double out = 0;
      for (auto &edge : paths_[i]) {
        double prev_swapped =
            amount_previously_swapped.contains(edge.id()) ? amount_previously_swapped.at(edge.id()) : 0.0;
        // temp[i] += edge.evaluate_swap(path_outputs[i], prev_swapped);
        // amount_previously_swapped[edge.id()] = path_outputs[i];
        // path_outputs = temp;
        out = edge.evaluate_swap(in, prev_swapped);
        amount_previously_swapped[edge.id()] += in;
        in = out;
      }
      path_outputs(i) = out;
    }

    auto result = path_outputs.sum();
    if (not std::isfinite(result)) {
      throw std::domain_error(fmt::format("edge network: f({}) is non-finite: path_outputs = {}",
                                          fmt::join(allocation, ","), fmt::join(path_outputs, ",")));
    }
    fmt::println(stderr, "f({}) = {} = {}", fmt::join(allocation, ","), fmt::join(path_outputs, "+"), result);
    return result;
  }

  Vector evaluate_gradient(const Vector &allocation, const double eps = 1.0e-6) const {
    Vector gradient(allocation.size());
    Vector shift_up = allocation;
    Vector shift_down = allocation;

    for (auto i = 0; i < gradient.size(); i++) {
      // Shift for gradient
      shift_up[i] += eps;
      shift_down[i] -= eps;

      gradient[i] = (evaluate_allocations(shift_up) - evaluate_allocations(shift_down)) / (2 * eps);

      // Shift back
      shift_up[i] -= eps;
      shift_down[i] += eps;
    }

    return gradient;
  }

  const std::vector<std::vector<EdgeType>> &paths() const { return paths_; }

  const size_t num_paths() const noexcept { return paths_.size(); }

 private:
  std::vector<std::vector<EdgeType>> paths_;
};

}  // namespace shortest_paths
