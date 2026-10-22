#include "SchwartzModel.hpp"
#include "ADESolver.hpp"
#include "MonteCarloLSMC.hpp"
#include <iostream>
#include <cmath>
#include <cassert>

using namespace schwartz;

int main() {
    SchwartzParameters p; // Default WTI Crude Oil parameters (Chapter 7.1)

    // Test 1: Terminal condition of futures price F(S, delta; 0) == S
    const double F_zero = SchwartzAnalytical::futuresPrice(p.S0, p.delta0, 0.0, p);
    assert(std::abs(F_zero - p.S0) < 1e-12);
    std::cout << "[PASS] Test 1: Futures terminal condition F(S, delta, 0) == S\n";

    // Test 2: Closed-form European Put matches paper benchmark (11.601861)
    const double euro_exact = SchwartzAnalytical::europeanPutClosedForm(p);
    assert(std::abs(euro_exact - 11.601861) < 1e-4);
    std::cout << "[PASS] Test 2: Closed-form European Put matches 11.601861 (" << euro_exact << ")\n";

    // Test 3: European ADE solver converges to closed-form price within 5e-3 on 100x100x600 grid
    ADESolver euro_ade(p, 100, 100, 600, SweepVariant::DoubleSweep, ExerciseStyle::EuropeanSpot);
    const double euro_ade_price = euro_ade.solve().price;
    assert(std::abs(euro_ade_price - euro_exact) < 5e-3);
    std::cout << "[PASS] Test 3: European ADE solver matches analytical closed-form (" << euro_ade_price << ")\n";

    // Test 4: American early-exercise premium is strictly positive (V_Amer > V_Euro)
    ADESolver amer_ds(p, 100, 100, 600, SweepVariant::DoubleSweep, ExerciseStyle::AmericanSpot);
    const ADEResult res_ds = amer_ds.solve();
    assert(res_ds.price > euro_ade_price);
    assert(std::abs(res_ds.price - 11.8907) < 0.02);
    std::cout << "[PASS] Test 4: American Put early-exercise premium verified (" << res_ds.price << ")\n";

    // Test 5: Double-Sweep and Saul'yev schemes agree within 5e-5 (Table 3)
    ADESolver amer_saulyev(p, 100, 100, 600, SweepVariant::Saulyev, ExerciseStyle::AmericanSpot);
    const ADEResult res_s = amer_saulyev.solve();
    assert(std::abs(res_ds.price - res_s.price) < 5e-5);
    std::cout << "[PASS] Test 5: Saul'yev vs. Double-Sweep agreement |S - DS| = "
              << std::abs(res_ds.price - res_s.price) << "\n";

    // Test 6: Two-futures hedge linear system (Eq. 54) reconstructs target sensitivities
    const double T1 = 0.25, T2 = 0.75;
    const auto [h1, h2] = SchwartzAnalytical::solveTwoFuturesHedge(
        res_ds.delta_S, res_ds.delta_yield, T1, T2, p);
    const double F1 = SchwartzAnalytical::futuresPrice(p.S0, p.delta0, T1, p);
    const double F2 = SchwartzAnalytical::futuresPrice(p.S0, p.delta0, T2, p);
    const double B1 = SchwartzAnalytical::B(T1, p.kappa);
    const double B2 = SchwartzAnalytical::B(T2, p.kappa);

    const double check_delta_S     = (F1 / p.S0) * h1 + (F2 / p.S0) * h2;
    const double check_delta_yield = -B1 * F1 * h1 - B2 * F2 * h2;
    assert(std::abs(check_delta_S - res_ds.delta_S) < 1e-10);
    assert(std::abs(check_delta_yield - res_ds.delta_yield) < 1e-10);
    std::cout << "[PASS] Test 6: Two-Futures hedging system (Eq. 54) exact inversion verified\n";

    std::cout << "\nAll 6 unit test suites passed successfully!\n";
    return 0;
}
