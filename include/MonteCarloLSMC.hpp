#ifndef MONTE_CARLO_LSMC_HPP
#define MONTE_CARLO_LSMC_HPP

#include "SchwartzModel.hpp"
#include <cstddef>
#include <cstdint>

namespace schwartz {

struct MCResult {
    double price       = 0.0;
    double std_error   = 0.0;
    double runtime_sec = 0.0;
};

/**
 * @brief Standard European Monte Carlo and Longstaff-Schwartz (2001) LSMC Engine
 * for benchmarking the 2D Schwartz commodity option PDE solver (Sections 6.2 & 6.3).
 */
class MonteCarloLSMC {
public:
    /**
     * @brief Simulates M paths of the joint (S_t, delta_t) process to price a European Put (Eq. 46).
     */
    static MCResult priceEuropeanPut(const SchwartzParameters& params,
                                     std::size_t num_paths = 250000,
                                     std::size_t num_steps = 100,
                                     std::uint64_t seed = 42);

    /**
     * @brief Prices an American Put using Longstaff-Schwartz (2001) LSMC (Section 6.3, Eq. 47)
     * with weighted Laguerre / two-factor cross-product basis functions.
     */
    static MCResult priceAmericanPutLSMC(const SchwartzParameters& params,
                                         std::size_t num_paths = 250000,
                                         std::size_t num_exercise_dates = 50,
                                         std::uint64_t seed = 42);
};

} // namespace schwartz

#endif // MONTE_CARLO_LSMC_HPP
