#ifndef SCHWARTZ_MODEL_HPP
#define SCHWARTZ_MODEL_HPP

#include <cmath>
#include <utility>
#include <stdexcept>

namespace schwartz {

/**
 * @brief Parameters for the Schwartz (1997) Two-Factor Commodity Model.
 * Default values correspond to the WTI Crude Oil calibration (Chapter 7.1).
 */
struct SchwartzParameters {
    double S0      = 90.0;    // Initial spot price
    double delta0  = 0.10;    // Initial instantaneous convenience yield
    double K       = 100.0;   // Option strike price
    double T       = 0.50;    // Option time to maturity (years)
    double r       = 0.05;    // Risk-free interest rate
    double kappa   = 1.49;    // Convenience yield mean-reversion speed
    double alpha_h = 0.0161;  // Risk-adjusted long-run convenience yield (\hat{\alpha} = \alpha - \lambda/\kappa)
    double sigma1  = 0.158;   // Spot volatility
    double sigma2  = 0.115;   // Convenience yield volatility
    double rho     = 0.189;   // Correlation between spot and convenience yield Brownian motions
};

/**
 * @brief Seasonal parameter specifications (Chapter 8.2, Eq. 51).
 * Evaluates periodic risk-adjusted convenience yield level and volatilities at time-to-maturity tau.
 */
struct SeasonalSchedule {
    static constexpr double PI = 3.14159265358979323846;

    static inline double alpha_hat(double tau, double base_alpha = 0.0161) noexcept {
        return base_alpha + 0.005 * std::sin(2.0 * PI * tau);
    }

    static inline double sigma1(double tau, double base_sigma1 = 0.158) noexcept {
        return base_sigma1 + 0.03 * std::cos(2.0 * PI * tau);
    }

    static inline double sigma2(double tau, double base_sigma2 = 0.115) noexcept {
        return base_sigma2 + 0.03 * std::abs(std::sin(PI * tau));
    }
};

/**
 * @brief Closed-form analytical formulas and hedging utilities under the Schwartz (1997) model.
 */
class SchwartzAnalytical {
public:
    static inline double normalCDF(double x) noexcept {
        return 0.5 * std::erfc(-x / std::sqrt(2.0));
    }

    // Equation (17): B(T) = (1 - exp(-kappa * T)) / kappa
    static inline double B(double tau, double kappa) noexcept {
        if ( std::abs(kappa) < 1e-12 ) return tau;
        return (1.0 - std::exp(-kappa * tau)) / kappa;
    }

    // Equation (18): A(T) closed-form affine factor for futures pricing
    static inline double A(double tau, const SchwartzParameters& p) noexcept {
        const double k  = p.kappa;
        const double k2 = k * k;
        const double k3 = k2 * k;
        const double s1 = p.sigma1;
        const double s2 = p.sigma2;
        const double s2_sq = s2 * s2;

        const double term1 = (p.r - p.alpha_h + 0.5 * s2_sq / k2 - (p.rho * s1 * s2) / k) * tau;
        const double term2 = s2_sq * (1.0 - std::exp(-2.0 * k * tau)) / (4.0 * k3);
        const double term3 = (p.alpha_h * k + p.rho * s1 * s2 - s2_sq / k)
                           * (1.0 - std::exp(-k * tau)) / k2;
        return term1 + term2 + term3;
    }

    // Equation (16): Closed-form Futures price F(S, delta; tau)
    static inline double futuresPrice(double S, double delta, double tau,
                                      const SchwartzParameters& p) noexcept {
        return S * std::exp(-delta * B(tau, p.kappa) + A(tau, p));
    }

    // Equation (20): Integrated variance v^2 of ln(S_T) under risk-neutral measure
    static inline double logSpotVariance(double tau, const SchwartzParameters& p) noexcept {
        const double k = p.kappa;
        const double exp_kT  = (1.0 - std::exp(-k * tau)) / k;
        const double exp_2kT = (1.0 - std::exp(-2.0 * k * tau)) / (2.0 * k);

        const double var_spot  = p.sigma1 * p.sigma1 * tau;
        const double var_yield = (p.sigma2 * p.sigma2 / (k * k)) * (tau - 2.0 * exp_kT + exp_2kT);
        const double cov_term  = (2.0 * p.rho * p.sigma1 * p.sigma2 / k) * (tau - exp_kT);
        return var_spot + var_yield - cov_term;
    }

    // Equation (19): Closed-form European Put option price on spot via Black (1976)
    static inline double europeanPutClosedForm(const SchwartzParameters& p) noexcept {
        const double F0 = futuresPrice(p.S0, p.delta0, p.T, p);
        const double v2 = logSpotVariance(p.T, p);
        const double v  = std::sqrt(v2);

        const double d1 = (std::log(F0 / p.K) + 0.5 * v2) / v;
        const double d2 = d1 - v;

        return std::exp(-p.r * p.T) * (p.K * normalCDF(-d2) - F0 * normalCDF(-d1));
    }

    /**
     * @brief Solves the 2-Futures Hedging Linear System (Chapter 8.3, Eq. 54):
     *   [  F1 / S       F2 / S     ] [ h1 ] = [ Delta_S     ]
     *   [ -B(T1)*F1    -B(T2)*F2   ] [ h2 ]   [ Delta_delta ]
     * Returns pair (h1, h2), where -h_k is the hedge position in contract F_k.
     */
    static inline std::pair<double, double> solveTwoFuturesHedge(
        double delta_S, double delta_yield, double T1, double T2,
        const SchwartzParameters& p) {

        if (std::abs(T1 - T2) < 1e-8) {
            throw std::invalid_argument("Futures maturities T1 and T2 must be distinct.");
        }

        const double F1 = futuresPrice(p.S0, p.delta0, T1, p);
        const double F2 = futuresPrice(p.S0, p.delta0, T2, p);
        const double B1 = B(T1, p.kappa);
        const double B2 = B(T2, p.kappa);

        const double det = (F1 * F2 / p.S0) * (B1 - B2);
        const double h1  = (-B2 * F2 * delta_S - (F2 / p.S0) * delta_yield) / det;
        const double h2  = ( B1 * F1 * delta_S + (F1 / p.S0) * delta_yield) / det;
        return {h1, h2};
    }
};

} // namespace schwartz

#endif // SCHWARTZ_MODEL_HPP
