#pragma once

#include "libs/Eigen/Core"
#include "libs/unsupported/Eigen/NNLS"
#include "fmt/core.h"

#include <limits>
#include <utility>
#include <exception>
#include <ranges>
#include <unordered_map>

namespace shortest_paths {

template <typename TSize = std::uint64_t>
class PolyEdge {
public:
    static constexpr const TSize INVALID_ID = std::numeric_limits<TSize>::max();

    using SelfType = PolyEdge<TSize>;

    // Default constructor with an invalid edge id and single coefficient of 1.0
    PolyEdge(): id_(INVALID_ID), coefficients_{{1.0}} {}
    PolyEdge(TSize id): id_(id), coefficients_{{1.0}} {}
    // Basic constructors from parts
    PolyEdge(TSize id, const Eigen::VectorXd &coefficients): id_(id), coefficients_(coefficients) {}
    PolyEdge(TSize id, Eigen::VectorXd &&coefficients): id_(id), coefficients_(std::move(coefficients)) {}

    // Copy and move constructors
    PolyEdge(const SelfType &other): id_(other.id_), coefficients_(other.coefficients_) {}
    PolyEdge(SelfType &&other): id_(std::move(other.id_)), coefficients_(std::move(other.coefficients_)) { }

    /// @brief Constructs a PolyEdge from an edge ID and set of points upon which to fit the polynomial.
    /// @param id ID of the edge.
    /// @param x Vector of x coordinates to which to fit the polynomial.
    /// @param y Vector of y coordinates to which to fit the polynomial.
    PolyEdge(TSize id, const Eigen::VectorXd& x, const Eigen::VectorXd& y): id_(id), coefficients_() {
        if (x.size() != y.size()) {
            throw std::domain_error(fmt::format("Number of x and y values must match ({} != {})", x.size(), y.size()));
        }
        calculate_coefficients(x, y);
    }

    /// @brief Constructs a PolyEdge from an edge ID and two contiguous containers of x and y values.
    /// @param id ID of the edge.
    /// @param x Contiguous container of x coordinates to which to fit the polynomial.
    /// @param y Contiguous container of y coordinates to which to fit the polynomial.
    template <std::ranges::contiguous_range XContainer, std::ranges::contiguous_range YContainer>
    PolyEdge(TSize id, const XContainer &x, const YContainer &y): id_(id), coefficients_() {
        static_assert (
            std::same_as<double, std::ranges::range_value_t<XContainer>>,
            "X values container must contain same floating point type as edge"
        );
        static_assert (
            std::same_as<double, std::ranges::range_value_t<YContainer>>,
            "Y values container must contain same floating point type as edge"
        );
        const auto x_size = std::ranges::size(x);
        const auto y_size = std::ranges::size(y);

        if (x_size != y_size) {
            throw std::domain_error(fmt::format("Number of x and y values must match ({} != {})", x_size, y_size));
        }

        Eigen::Map<const Eigen::VectorXd> x_map(std::ranges::cdata(x), x_size);
        Eigen::Map<const Eigen::VectorXd> y_map(std::ranges::cdata(y), y_size);

        calculate_coefficients(x_map, y_map);
    }

    // Copy and move assignment
    SelfType& operator=(const SelfType& other) {
        id_ = other.id_;
        coefficients_ = other.coefficients_;
        return *this;
    }
    SelfType& operator=(SelfType &&other) {
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
    TSize id() const {
        return id_;
    }
    Eigen::VectorXd coefficients() const {
        return coefficients_;
    }

    /// @brief Update the coefficients of the fitted polynomial to the given coefficients.
    /// @param coefficients The new polynomial coefficients.
    /// @return This edge.
    SelfType& set_coefficients(const Eigen::VectorXd &coefficients) {
        coefficients_ = coefficients;
        return *this;
    }

    /// @brief Calculates the coefficients using an NNLS solver. This ensures that the polynomial
    ///     is strictly monotonically increasing.
    /// @param x The vector of x-values for points along the curve.
    /// @param y The vector of y-values for points along the curve.
    /// @throws std::logic_error if the NNLS solver somehow returns negative coefficients.
    template <typename DerivedX, typename DerivedY>
    void calculate_coefficients(const Eigen::MatrixBase<DerivedX> &x, const Eigen::MatrixBase<DerivedY> &y) {
        // Initialize the matrix as required by the NNLS solver.
        Eigen::MatrixXd matrix(x.size(), y.size());
        for (auto i = 0; i < x.size(); i++) {
            for (auto j = 0; j < y.size(); j++) {
                matrix(i, j) = pow(x(i), j);
            }
        }

        Eigen::NNLS<Eigen::MatrixXd> nnls_solver(matrix);
        this->coefficients_ = nnls_solver.solve(y);

        for (auto i = 0; i < this->coefficients_.size(); i++) {
            if (coefficients_[i] < 0.0) {
                throw std::logic_error(fmt::format("Non-negative LLS solver somehow returned negative coefficient[{}]={}", i, coefficients_[i]));
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
                retval += coefficients_[i] * pow(x, i+1);
            } else {
                retval += coefficients_[i];
            }
        }
        return retval;
    }

    /// @brief Evaluates the inverse of the polynomial at the given x coordinate.
    /// @param x The x value at which to evaulate the polynomial.
    /// @return The reciprocal of the corresponding y value of the polynomial.
    double evaluate_inv(double x) const {
        return 1.0 / evaluate(x);
    }

    double evaluate_int(double x) const {
        double retval = 0;
        for (auto i = 0; i < coefficients_.size(); i++) {
            if (i != 0) {
                retval += coefficients_[i] * pow(x, i+1) / i;
            } else {
                retval += coefficients_[i] * pow(x, 1);
            }
        }
        return retval;
    }

    double evaluate_avg_price(double x, double prev_swapped = 0.0) {
        // Average weight (price) is the integral of our polynomial divided by the
        // amount swapped to get there.
        double integral = evaluate_int(x + prev_swapped) - evaluate_int(prev_swapped);
        return integral / x;
    }

    double evaluate_swap(double x, double prev_swapped = 0.0) {
        // A swap (Out / In) * In execuated at the "average" price paid.
        return x / evaluate_avg_price(x);
    }

private:
    TSize id_;
    Eigen::VectorXd coefficients_;
};

template <typename TSize>
void swap(PolyEdge<TSize>& v1, PolyEdge<TSize>& v2) {
    v1.swap(v2);
}

template <typename TSize = std::uint64_t>
class EdgeNetwork {
public:
    using SelfType = EdgeNetwork<TSize>;
    using EdgeType = PolyEdge<TSize>;
    using Vector = Eigen::VectorXd;

    EdgeNetwork(const std::vector<std::vector<EdgeType>> &paths): paths_(paths) { }
    EdgeNetwork(std::vector<std::vector<EdgeType>> &&paths): paths_(std::move(paths)) { }

    double evaluate_allocations(const Vector &allocation) {
        Vector path_outputs = allocation;
        Vector temp(allocation.size());

        std::unordered_map<TSize, double> amount_previously_swapped;

        for (size_t i = 0; i < paths_.size(); i++) {
            temp.setZero();
            for (auto &edge : paths_[i]) {
                double prev_swapped = amount_previously_swapped.contains(edge.id()) ? amount_previously_swapped.at(edge.id()) : 0.0;
                temp[i] += edge.evaluate_swap(path_outputs[i], prev_swapped);
                amount_previously_swapped[edge.id()] = path_outputs[i];
            }
            path_outputs = temp;
        }

        return path_outputs.sum();
    }

    Vector evaluate_gradient(const Vector &allocation, const double eps = 1.0e-6) {
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

    const std::vector<std::vector<EdgeType>>& paths() const {
        return paths_;
    }
private:
    std::vector<std::vector<EdgeType>> paths_;
};

} // namespace shortest_paths