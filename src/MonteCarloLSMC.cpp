#include "MonteCarloLSMC.hpp"
#include <vector>
#include <random>
#include <cmath>
#include <algorithm>
#include <chrono>

namespace schwartz {

namespace {

// Evaluates the 6-term two-factor weighted Laguerre & cross-product basis (Section 6.3, Eq. 47)
constexpr std::size_t NUM_BASIS = 6;

inline void evaluateBasis(double S_norm, double delta, double* phi) noexcept {
    const double w  = std::exp(-0.5 * S_norm);
    const double L0 = 1.0;
    const double L1 = 1.0 - S_norm;
    const double L2 = 1.0 - 2.0 * S_norm + 0.5 * S_norm * S_norm;

    phi[0] = 1.0;
    phi[1] = w * L0;
    phi[2] = w * L1;
    phi[3] = w * L2;
    phi[4] = delta;
    phi[5] = delta * delta + S_norm * delta;
}

// Solves A * beta = b via Gaussian elimination with partial pivoting
std::vector<double> solveLinearSystem(std::vector<std::vector<double>> A,
                                      std::vector<double> b) {
    const std::size_t n = b.size();
    for (std::size_t k = 0; k < n; ++k) {
        std::size_t max_row = k;
        double max_val = std::abs(A[k][k]);
        for (std::size_t i = k + 1; i < n; ++i) {
            if (std::abs(A[i][k]) > max_val) {
                max_val = std::abs(A[i][k]);
                max_row = i;
            }
        }
        std::swap(A[k], A[max_row]);
        std::swap(b[k], b[max_row]);

        if (std::abs(A[k][k]) < 1e-12) continue;

        for (std::size_t i = k + 1; i < n; ++i) {
            const double factor = A[i][k] / A[k][k];
            for (std::size_t j = k; j < n; ++j) {
                A[i][j] -= factor * A[k][j];
            }
            b[i] -= factor * b[k];
        }
    }

    std::vector<double> x(n, 0.0);
    for (int i = static_cast<int>(n) - 1; i >= 0; --i) {
        const std::size_t ui = static_cast<std::size_t>(i);
        if (std::abs(A[ui][ui]) < 1e-12) continue;
        double sum = b[ui];
        for (std::size_t j = ui + 1; j < n; ++j) {
            sum -= A[ui][j] * x[j];
        }
        x[ui] = sum / A[ui][ui];
    }
    return x;
}

} // anonymous namespace

MCResult MonteCarloLSMC::priceEuropeanPut(const SchwartzParameters& p,
                                          std::size_t num_paths,
                                          std::size_t num_steps,
                                          std::uint64_t seed) {
    const auto t_start = std::chrono::high_resolution_clock::now();

    std::mt19937_64 rng(seed);
    std::normal_distribution<double> norm(0.0, 1.0);

    const double dt = p.T / static_cast<double>(num_steps);
    const double sqrt_dt = std::sqrt(dt);
    const double rho_comp = std::sqrt(1.0 - p.rho * p.rho);
    const double half_sig1_sq = 0.5 * p.sigma1 * p.sigma1;

    double sum_payoff = 0.0;
    double sum_sq_payoff = 0.0;

    for (std::size_t m = 0; m < num_paths; ++m) {
        double x = std::log(p.S0);
        double delta = p.delta0;

        for (std::size_t step = 0; step < num_steps; ++step) {
            const double z1 = norm(rng);
            const double z2 = p.rho * z1 + rho_comp * norm(rng);

            const double dx = (p.r - delta - half_sig1_sq) * dt + p.sigma1 * sqrt_dt * z1;
            const double d_delta = p.kappa * (p.alpha_h - delta) * dt + p.sigma2 * sqrt_dt * z2;

            x += dx;
            delta += d_delta;
        }

        const double ST = std::exp(x);
        const double payoff = std::max(p.K - ST, 0.0);
        sum_payoff += payoff;
        sum_sq_payoff += payoff * payoff;
    }

    const double disc = std::exp(-p.r * p.T);
    const double mean_payoff = sum_payoff / static_cast<double>(num_paths);
    const double var_payoff = (sum_sq_payoff / static_cast<double>(num_paths))
                            - (mean_payoff * mean_payoff);

    const auto t_end = std::chrono::high_resolution_clock::now();

    MCResult res;
    res.price = disc * mean_payoff;
    res.std_error = disc * std::sqrt(var_payoff / static_cast<double>(num_paths));
    res.runtime_sec = std::chrono::duration<double>(t_end - t_start).count();
    return res;
}

MCResult MonteCarloLSMC::priceAmericanPutLSMC(const SchwartzParameters& p,
                                              std::size_t num_paths,
                                              std::size_t num_exercise_dates,
                                              std::uint64_t seed) {
    const auto t_start = std::chrono::high_resolution_clock::now();

    std::mt19937_64 rng(seed);
    std::normal_distribution<double> norm(0.0, 1.0);

    const double dt = p.T / static_cast<double>(num_exercise_dates);
    const double sqrt_dt = std::sqrt(dt);
    const double disc_step = std::exp(-p.r * dt);
    const double rho_comp = std::sqrt(1.0 - p.rho * p.rho);
    const double half_sig1_sq = 0.5 * p.sigma1 * p.sigma1;

    // Store paths: S_paths[step * num_paths + m], D_paths[step * num_paths + m]
    std::vector<double> S_paths((num_exercise_dates + 1) * num_paths);
    std::vector<double> D_paths((num_exercise_dates + 1) * num_paths);

    for (std::size_t m = 0; m < num_paths; ++m) {
        S_paths[m] = p.S0;
        D_paths[m] = p.delta0;
    }

    for (std::size_t k = 1; k <= num_exercise_dates; ++k) {
        const std::size_t prev_offset = (k - 1) * num_paths;
        const std::size_t curr_offset = k * num_paths;

        for (std::size_t m = 0; m < num_paths; ++m) {
            const double S_prev = S_paths[prev_offset + m];
            const double d_prev = D_paths[prev_offset + m];

            const double z1 = norm(rng);
            const double z2 = p.rho * z1 + rho_comp * norm(rng);

            const double S_next = S_prev * std::exp((p.r - d_prev - half_sig1_sq) * dt
                                                    + p.sigma1 * sqrt_dt * z1);
            const double d_next = d_prev + p.kappa * (p.alpha_h - d_prev) * dt
                                         + p.sigma2 * sqrt_dt * z2;

            S_paths[curr_offset + m] = S_next;
            D_paths[curr_offset + m] = d_next;
        }
    }

    // Terminal cash flows at maturity k = num_exercise_dates
    std::vector<double> cash_flows(num_paths);
    const std::size_t last_offset = num_exercise_dates * num_paths;
    for (std::size_t m = 0; m < num_paths; ++m) {
        cash_flows[m] = std::max(p.K - S_paths[last_offset + m], 0.0);
    }

    // Backward induction through exercise dates k = num_exercise_dates - 1 down to 1
    for (std::size_t k = num_exercise_dates - 1; k >= 1; --k) {
        const std::size_t offset = k * num_paths;

        std::vector<std::vector<double>> XTX(NUM_BASIS, std::vector<double>(NUM_BASIS, 0.0));
        std::vector<double> XTY(NUM_BASIS, 0.0);
        std::vector<std::size_t> itm_indices;
        itm_indices.reserve(num_paths / 2);

        double phi[NUM_BASIS];

        for (std::size_t m = 0; m < num_paths; ++m) {
            cash_flows[m] *= disc_step;
            const double S_curr = S_paths[offset + m];
            if (p.K > S_curr) {
                itm_indices.push_back(m);
                evaluateBasis(S_curr / p.K, D_paths[offset + m], phi);
                const double y_val = cash_flows[m];
                for (std::size_t r = 0; r < NUM_BASIS; ++r) {
                    XTY[r] += phi[r] * y_val;
                    for (std::size_t c = r; c < NUM_BASIS; ++c) {
                        XTX[r][c] += phi[r] * phi[c];
                    }
                }
            }
        }

        for (std::size_t r = 0; r < NUM_BASIS; ++r) {
            for (std::size_t c = 0; c < r; ++c) {
                XTX[r][c] = XTX[c][r];
            }
        }

        if (itm_indices.size() > NUM_BASIS) {
            const std::vector<double> beta = solveLinearSystem(XTX, XTY);

            for (std::size_t m : itm_indices) {
                const double S_curr = S_paths[offset + m];
                const double immediate = p.K - S_curr;
                evaluateBasis(S_curr / p.K, D_paths[offset + m], phi);

                double continuation = 0.0;
                for (std::size_t b = 0; b < NUM_BASIS; ++b) {
                    continuation += beta[b] * phi[b];
                }
                if (immediate > continuation) {
                    cash_flows[m] = immediate;
                }
            }
        }
    }

    // Discount from k = 1 to k = 0 and check immediate exercise at t = 0
    double sum_cf = 0.0;
    double sum_sq_cf = 0.0;
    for (std::size_t m = 0; m < num_paths; ++m) {
        const double val = cash_flows[m] * disc_step;
        sum_cf += val;
        sum_sq_cf += val * val;
    }

    const double mean_val = std::max(sum_cf / static_cast<double>(num_paths),
                                     std::max(p.K - p.S0, 0.0));
    const double var_val = (sum_sq_cf / static_cast<double>(num_paths))
                         - (mean_val * mean_val);

    const auto t_end = std::chrono::high_resolution_clock::now();

    MCResult res;
    res.price = mean_val;
    res.std_error = std::sqrt(std::max(var_val, 0.0) / static_cast<double>(num_paths));
    res.runtime_sec = std::chrono::duration<double>(t_end - t_start).count();
    return res;
}

} // namespace schwartz
