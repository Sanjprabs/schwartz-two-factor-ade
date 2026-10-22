#include "ADESolver.hpp"
#include <cmath>
#include <algorithm>
#include <chrono>
#include <stdexcept>

namespace schwartz {

ADESolver::ADESolver(const SchwartzParameters& params,
                     std::size_t NX,
                     std::size_t NY,
                     std::size_t NT,
                     SweepVariant variant,
                     ExerciseStyle style,
                     bool use_seasonal,
                     double futures_residual_theta)
    : params_(params),
      NX_(NX),
      NY_(NY),
      NT_(NT),
      variant_(variant),
      style_(style),
      use_seasonal_(use_seasonal),
      theta_(futures_residual_theta) {

    if (NX_ < 4 || NY_ < 4 || NT_ < 2) {
        throw std::invalid_argument("Grid dimensions must satisfy NX >= 4, NY >= 4, NT >= 2.");
    }
    if (variant_ == SweepVariant::Saulyev && (NT_ % 2 != 0)) {
        // Section 6.1: Saul'yev requires an even number of time steps so upward/downward sweeps pair
        NT_ += 1;
    }
}

double ADESolver::payoffObstacle(double x, double y, double tau) const noexcept {
    const double S = std::exp(x);
    if (style_ == ExerciseStyle::AmericanFutures) {
        // Section 8.1, Eq. (50): Obstacle on a futures contract with time-to-delivery tau + theta
        const double F = SchwartzAnalytical::futuresPrice(S, y, tau + theta_, params_);
        return std::max(params_.K - F, 0.0);
    }
    return std::max(params_.K - S, 0.0);
}

void ADESolver::applySpotBoundaries(Matrix2D& grid, double tau,
                                    const std::vector<double>& x_nodes,
                                    const std::vector<double>& y_nodes) const noexcept {
    const double S_min = std::exp(x_nodes[0]);
    for (std::size_t j = 0; j <= NY_; ++j) {
        // Lower spot boundary x_min (S -> 0) via Fichera theory (Section 4.3, Eq. 38-39)
        if (style_ == ExerciseStyle::EuropeanSpot) {
            const double F_min = SchwartzAnalytical::futuresPrice(S_min, y_nodes[j], tau, params_);
            grid(0, j) = std::max(params_.K * std::exp(-params_.r * tau)
                                - F_min * std::exp(-params_.r * tau), 0.0);
        } else {
            grid(0, j) = payoffObstacle(x_nodes[0], y_nodes[j], tau);
        }
        // Upper spot boundary x_max (S -> inf): V(x_max, y, tau) = 0
        grid(NX_, j) = 0.0;
    }
}

void ADESolver::applyYieldBoundaries(Matrix2D& grid) const noexcept {
    // Section 4.3, Eq. (40): Linear extrapolation V_yy = 0 at convenience yield boundaries
    for (std::size_t i = 0; i <= NX_; ++i) {
        grid(i, 0)   = 2.0 * grid(i, 1) - grid(i, 2);
        grid(i, NY_) = 2.0 * grid(i, NY_ - 1) - grid(i, NY_ - 2);
    }
}

ADEResult ADESolver::solve() const {
    const auto t_start = std::chrono::high_resolution_clock::now();

    // Construct uniform mesh centered at (ln S0, delta0) (Section 4.1, Eq. 34)
    const double x0 = std::log(params_.S0);
    const double y0 = params_.delta0;
    const double x_min = x0 - X_HALF_WIDTH;
    const double y_min = y0 - Y_HALF_WIDTH;

    const double hx   = (2.0 * X_HALF_WIDTH) / static_cast<double>(NX_);
    const double hy   = (2.0 * Y_HALF_WIDTH) / static_cast<double>(NY_);
    const double dtau = params_.T / static_cast<double>(NT_);

    std::vector<double> x_nodes(NX_ + 1);
    std::vector<double> y_nodes(NY_ + 1);
    for (std::size_t i = 0; i <= NX_; ++i) x_nodes[i] = x_min + static_cast<double>(i) * hx;
    for (std::size_t j = 0; j <= NY_; ++j) y_nodes[j] = y_min + static_cast<double>(j) * hy;

    // Terminal condition V(x, y, tau = 0) = psi(x, y, 0)
    Matrix2D u(NX_ + 1, NY_ + 1, 0.0);
    for (std::size_t i = 0; i <= NX_; ++i) {
        for (std::size_t j = 0; j <= NY_; ++j) {
            u(i, j) = payoffObstacle(x_nodes[i], y_nodes[j], 0.0);
        }
    }

    Matrix2D U(NX_ + 1, NY_ + 1, 0.0);
    Matrix2D V(NX_ + 1, NY_ + 1, 0.0);

    const double inv_2hx = 1.0 / (2.0 * hx);
    const double inv_2hy = 1.0 / (2.0 * hy);
    const double inv_hx2 = 1.0 / (hx * hx);
    const double inv_hy2 = 1.0 / (hy * hy);
    const double inv_2hxhy = 1.0 / (2.0 * hx * hy);

    // Time-marching loop n = 0, ..., NT - 1 (Section 6.1)
    for (std::size_t n = 0; n < NT_; ++n) {
        const double tau_curr = static_cast<double>(n) * dtau;
        const double tau_next = static_cast<double>(n + 1) * dtau;

        // Evaluate constant or seasonal coefficients (Section 4.2 & Section 8.2)
        const double sig1 = use_seasonal_
            ? SeasonalSchedule::sigma1(tau_curr, params_.sigma1) : params_.sigma1;
        const double sig2 = use_seasonal_
            ? SeasonalSchedule::sigma2(tau_curr, params_.sigma2) : params_.sigma2;
        const double alpha_h = use_seasonal_
            ? SeasonalSchedule::alpha_hat(tau_curr, params_.alpha_h) : params_.alpha_h;

        // Transformed diffusion tensor entries (Section 4.2, Eq. 36)
        const double a11 = 0.5 * sig1 * sig1;
        const double a12 = 0.5 * params_.rho * sig1 * sig2;
        const double a22 = 0.5 * sig2 * sig2;

        // Discrete finite difference coefficients (Section 3.4, Eq. 29)
        const double cx   = a11 * inv_hx2;
        const double cy   = a22 * inv_hy2;
        const double cxy  = a12 * inv_2hxhy;
        const double beta = cx + cy;

        const double denom      = 1.0 + dtau * beta;
        const double inv_denom  = 1.0 / denom;
        const double diag_coeff = 1.0 - dtau * (beta + params_.r);

        if (variant_ == SweepVariant::DoubleSweep) {
            U = u;
            V = u;
            applySpotBoundaries(U, tau_next, x_nodes, y_nodes);
            applySpotBoundaries(V, tau_next, x_nodes, y_nodes);
            applyYieldBoundaries(U);
            applyYieldBoundaries(V);

            // Step 2a: Upward sweep U^{n+1} (Eq. 31, i and j increasing)
            for (std::size_t i = 1; i < NX_; ++i) {
                for (std::size_t j = 1; j < NY_; ++j) {
                    const double y  = y_nodes[j];
                    const double b1 = params_.r - y - a11;
                    const double b2 = params_.kappa * (alpha_h - y);

                    const double d_xy = u(i + 1, j + 1) - u(i + 1, j - 1)
                                      - u(i - 1, j + 1) + u(i - 1, j - 1);

                    const double rhs = diag_coeff * u(i, j) + dtau * (
                        cx * (u(i + 1, j) + U(i - 1, j)) +
                        cy * (u(i, j + 1) + U(i, j - 1)) +
                        b1 * inv_2hx * (u(i + 1, j) - U(i - 1, j)) +
                        b2 * inv_2hy * (u(i, j + 1) - U(i, j - 1)) +
                        cxy * d_xy
                    );
                    U(i, j) = rhs * inv_denom;
                }
            }

            // Step 2b: Downward sweep V^{n+1} (Eq. 32, i and j decreasing)
            for (std::size_t i = NX_ - 1; i >= 1; --i) {
                for (std::size_t j = NY_ - 1; j >= 1; --j) {
                    const double y  = y_nodes[j];
                    const double b1 = params_.r - y - a11;
                    const double b2 = params_.kappa * (alpha_h - y);

                    const double d_xy = u(i + 1, j + 1) - u(i + 1, j - 1)
                                      - u(i - 1, j + 1) + u(i - 1, j - 1);

                    const double rhs = diag_coeff * u(i, j) + dtau * (
                        cx * (u(i - 1, j) + V(i + 1, j)) +
                        cy * (u(i, j - 1) + V(i, j + 1)) +
                        b1 * inv_2hx * (V(i + 1, j) - u(i - 1, j)) +
                        b2 * inv_2hy * (V(i, j + 1) - u(i, j - 1)) +
                        cxy * d_xy
                    );
                    V(i, j) = rhs * inv_denom;
                }
            }

            // Step 3: Average U^{n+1} and V^{n+1} (Eq. 25) and refresh boundaries
            for (std::size_t i = 1; i < NX_; ++i) {
                for (std::size_t j = 1; j < NY_; ++j) {
                    u(i, j) = 0.5 * (U(i, j) + V(i, j));
                }
            }
            applySpotBoundaries(u, tau_next, x_nodes, y_nodes);
            applyYieldBoundaries(u);

        } else {
            // Saul'yev Alternating Single-Sweep Scheme (Section 3.3 & Section 6.1)
            const Matrix2D u_old = u;
            applySpotBoundaries(u, tau_next, x_nodes, y_nodes);
            applyYieldBoundaries(u);

            if (n % 2 == 0) {
                // Even step: Upward sweep (Eq. 31)
                for (std::size_t i = 1; i < NX_; ++i) {
                    for (std::size_t j = 1; j < NY_; ++j) {
                        const double y  = y_nodes[j];
                        const double b1 = params_.r - y - a11;
                        const double b2 = params_.kappa * (alpha_h - y);

                        const double d_xy = u_old(i + 1, j + 1) - u_old(i + 1, j - 1)
                                          - u_old(i - 1, j + 1) + u_old(i - 1, j - 1);

                        const double rhs = diag_coeff * u_old(i, j) + dtau * (
                            cx * (u_old(i + 1, j) + u(i - 1, j)) +
                            cy * (u_old(i, j + 1) + u(i, j - 1)) +
                            b1 * inv_2hx * (u_old(i + 1, j) - u(i - 1, j)) +
                            b2 * inv_2hy * (u_old(i, j + 1) - u(i, j - 1)) +
                            cxy * d_xy
                        );
                        u(i, j) = rhs * inv_denom;
                    }
                }
            } else {
                // Odd step: Downward sweep (Eq. 32)
                for (std::size_t i = NX_ - 1; i >= 1; --i) {
                    for (std::size_t j = NY_ - 1; j >= 1; --j) {
                        const double y  = y_nodes[j];
                        const double b1 = params_.r - y - a11;
                        const double b2 = params_.kappa * (alpha_h - y);

                        const double d_xy = u_old(i + 1, j + 1) - u_old(i + 1, j - 1)
                                          - u_old(i - 1, j + 1) + u_old(i - 1, j - 1);

                        const double rhs = diag_coeff * u_old(i, j) + dtau * (
                            cx * (u_old(i - 1, j) + u(i + 1, j)) +
                            cy * (u_old(i, j - 1) + u(i, j + 1)) +
                            b1 * inv_2hx * (u(i + 1, j) - u_old(i - 1, j)) +
                            b2 * inv_2hy * (u(i, j + 1) - u_old(i, j - 1)) +
                            cxy * d_xy
                        );
                        u(i, j) = rhs * inv_denom;
                    }
                }
            }
            applySpotBoundaries(u, tau_next, x_nodes, y_nodes);
            applyYieldBoundaries(u);
        }

        // Step 4: Discrete enforcement of the American exercise constraint (Section 5.3, Eq. 43)
        if (style_ != ExerciseStyle::EuropeanSpot) {
            for (std::size_t i = 0; i <= NX_; ++i) {
                for (std::size_t j = 0; j <= NY_; ++j) {
                    u(i, j) = std::max(u(i, j), payoffObstacle(x_nodes[i], y_nodes[j], tau_next));
                }
            }
        }
    }

    const auto t_end = std::chrono::high_resolution_clock::now();
    const double elapsed = std::chrono::duration<double>(t_end - t_start).count();

    // Extract price and central-difference two-factor Greeks at center node (ic, jc) (Section 8.3, Eq. 52)
    const std::size_t ic = NX_ / 2;
    const std::size_t jc = NY_ / 2;

    const double V_center = u(ic, jc);
    const double V_x  = (u(ic + 1, jc) - u(ic - 1, jc)) / (2.0 * hx);
    const double V_xx = (u(ic + 1, jc) - 2.0 * V_center + u(ic - 1, jc)) / (hx * hx);
    const double V_y  = (u(ic, jc + 1) - u(ic, jc - 1)) / (2.0 * hy);

    ADEResult res;
    res.price       = V_center;
    res.delta_S     = V_x / params_.S0;
    res.gamma_S     = (V_xx - V_x) / (params_.S0 * params_.S0);
    res.delta_yield = V_y;
    res.runtime_sec = elapsed;
    return res;
}

} // namespace schwartz
