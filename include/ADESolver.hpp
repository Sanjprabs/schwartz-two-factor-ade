#ifndef ADE_SOLVER_HPP
#define ADE_SOLVER_HPP

#include "Matrix2D.hpp"
#include "SchwartzModel.hpp"
#include <cstddef>

namespace schwartz {

enum class SweepVariant {
    DoubleSweep,  // Barakat & Clark (1966) symmetric upward + downward average (Section 3.2)
    Saulyev       // Saul'yev (1964) single-sweep alternating direction in time (Section 3.3)
};

enum class ExerciseStyle {
    EuropeanSpot,     // Unconstrained PDE solve for European put on spot (Section 6.2)
    AmericanSpot,     // Variational inequality projection on spot max(K - S, 0) (Chapter 5)
    AmericanFutures   // Exercise obstacle on a futures contract max(K - F(S,delta;tau+theta), 0) (Section 8.1)
};

struct ADEResult {
    double price       = 0.0;
    double delta_S     = 0.0;  // dV / dS = V_x / S (Eq. 52)
    double gamma_S     = 0.0;  // d^2V / dS^2 = (V_xx - V_x) / S^2 (Eq. 52)
    double delta_yield = 0.0;  // dV / d(delta) = V_y (Eq. 52)
    double runtime_sec = 0.0;
};

/**
 * @brief Two-Dimensional Alternating Direction Explicit (ADE) PDE Solver
 * for the transformed Schwartz (1997) Two-Factor Commodity Model.
 */
class ADESolver {
public:
    ADESolver(const SchwartzParameters& params,
              std::size_t NX = 200,
              std::size_t NY = 200,
              std::size_t NT = 600,
              SweepVariant variant = SweepVariant::DoubleSweep,
              ExerciseStyle style = ExerciseStyle::AmericanSpot,
              bool use_seasonal = false,
              double futures_residual_theta = 0.25);

    ADEResult solve() const;

private:
    SchwartzParameters params_;
    std::size_t NX_;
    std::size_t NY_;
    std::size_t NT_;
    SweepVariant variant_;
    ExerciseStyle style_;
    bool use_seasonal_;
    double theta_; // Residual futures maturity T_F - T at option expiry (Section 8.1)

    // Domain half-widths centered at (ln S0, delta0) (Section 4.1, Eq. 34)
    static constexpr double X_HALF_WIDTH = 1.2;
    static constexpr double Y_HALF_WIDTH = 0.4;

    double payoffObstacle(double x, double y, double tau) const noexcept;
    void applySpotBoundaries(Matrix2D& grid, double tau,
                             const std::vector<double>& x_nodes,
                             const std::vector<double>& y_nodes) const noexcept;
    void applyYieldBoundaries(Matrix2D& grid) const noexcept;
};

} // namespace schwartz

#endif // ADE_SOLVER_HPP
