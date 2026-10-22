#include "SchwartzModel.hpp"
#include "ADESolver.hpp"
#include "MonteCarloLSMC.hpp"
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace schwartz;

int main() {
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "============================================================================\n";
    std::cout << "  SCHWARTZ TWO-FACTOR COMMODITY OPTION PRICER: 2D ADE & LSMC BENCHMARK SUITE\n";
    std::cout << "  Authors: Sanjith Prabhu & Dr. Daniel J. Duffy (arXiv:2610.21145)\n";
    std::cout << "============================================================================\n\n";

    SchwartzParameters base_params; // WTI Crude Oil Calibration (Chapter 7.1)

    // -------------------------------------------------------------------------
    // 1. Section 7.2: European Option Validation
    // -------------------------------------------------------------------------
    std::cout << "[1] SECTION 7.2: EUROPEAN PUT VALIDATION (S0 = 90, delta0 = 0.1, K = 100)\n";
    std::cout << "----------------------------------------------------------------------------\n";
    const double euro_exact = SchwartzAnalytical::europeanPutClosedForm(base_params);
    const ADESolver euro_ade_solver(base_params, 400, 400, 1200,
                                    SweepVariant::DoubleSweep, ExerciseStyle::EuropeanSpot);
    const ADEResult euro_ade = euro_ade_solver.solve();
    const MCResult  euro_mc  = MonteCarloLSMC::priceEuropeanPut(base_params, 250000, 100, 42);

    std::cout << "  Analytical Closed-Form (Eq. 19) : " << euro_exact << "\n";
    std::cout << "  Monte Carlo (250k paths)        : " << euro_mc.price
              << "  (SE: " << euro_mc.std_error << ")\n";
    std::cout << "  Double-Sweep ADE (400x400x1200) : " << euro_ade.price
              << "  (|Error|: " << std::scientific << std::abs(euro_ade.price - euro_exact)
              << std::fixed << ")\n\n";

    // -------------------------------------------------------------------------
    // 2. Section 7.3: Table 1 Double-Sweep ADE Grid Convergence Study
    // -------------------------------------------------------------------------
    std::cout << "[2] SECTION 7.3 (TABLE 1): DOUBLE-SWEEP ADE AMERICAN PUT CONVERGENCE\n";
    std::cout << "----------------------------------------------------------------------------\n";
    const std::vector<std::size_t> NY_list = {50, 100, 200, 400};
    const std::vector<std::size_t> NT_list = {100, 300, 600, 1200};

    std::cout << std::setw(10) << "NY \\ NT";
    for (std::size_t nt : NT_list) std::cout << std::setw(14) << nt;
    std::cout << "\n";

    for (std::size_t ny : NY_list) {
        std::cout << std::setw(10) << ny;
        for (std::size_t nt : NT_list) {
            ADESolver solver(base_params, ny, ny, nt,
                             SweepVariant::DoubleSweep, ExerciseStyle::AmericanSpot);
            std::cout << std::setw(14) << solver.solve().price;
        }
        std::cout << "\n";
    }
    std::cout << "\n";

    // -------------------------------------------------------------------------
    // 3. Section 7.4: Table 2 Market Scenario Analysis
    // -------------------------------------------------------------------------
    std::cout << "[3] SECTION 7.4 (TABLE 2): MARKET SCENARIO ANALYSIS (ADE vs. 250k LSMC)\n";
    std::cout << "----------------------------------------------------------------------------\n";
    struct Scenario {
        double S0;
        double delta0;
        const char* label;
    };
    const std::vector<Scenario> scenarios = {
        {90.0,  0.10, "ITM (Base)"},
        {100.0, 0.10, "ATM"},
        {90.0,  0.30, "ITM, High Yield"}
    };

    std::cout << std::left << std::setw(18) << "Scenario"
              << std::right << std::setw(8) << "S0"
              << std::setw(10) << "delta0"
              << std::setw(14) << "ADE Price"
              << std::setw(14) << "LSMC Price"
              << std::setw(14) << "Abs Diff" << "\n";

    for (const auto& sc : scenarios) {
        SchwartzParameters p = base_params;
        p.S0     = sc.S0;
        p.delta0 = sc.delta0;

        ADESolver ade_sc(p, 400, 400, 1200, SweepVariant::DoubleSweep, ExerciseStyle::AmericanSpot);
        const double v_ade  = ade_sc.solve().price;
        const double v_lsmc = MonteCarloLSMC::priceAmericanPutLSMC(p, 250000, 50, 42).price;

        std::cout << std::left << std::setw(18) << sc.label
                  << std::right << std::setprecision(2) << std::setw(8) << sc.S0
                  << std::setw(10) << sc.delta0
                  << std::setprecision(6) << std::setw(14) << v_ade
                  << std::setw(14) << v_lsmc
                  << std::setw(14) << std::abs(v_ade - v_lsmc) << "\n";
    }
    std::cout << "\n";

    // -------------------------------------------------------------------------
    // 4. Sections 7.5 & 7.6: Double-Sweep vs. Saul'yev & Runtime Benchmark
    // -------------------------------------------------------------------------
    std::cout << "[4] SECTIONS 7.5 & 7.6 (TABLES 3 & 4): DOUBLE-SWEEP vs. SAUL'YEV & RUNTIMES\n";
    std::cout << "----------------------------------------------------------------------------\n";
    struct GridConfig {
        std::size_t NY;
        std::size_t NT;
    };
    const std::vector<GridConfig> bench_grids = {
        {50,  300},
        {100, 600},
        {200, 600},
        {200, 1200}
    };

    std::cout << std::left << std::setw(16) << "Grid (NY x NT)"
              << std::right << std::setw(14) << "Double-Sweep"
              << std::setw(14) << "Saul'yev"
              << std::setw(14) << "|S - DS|"
              << std::setw(12) << "DS Time(s)"
              << std::setw(12) << "S Time(s)" << "\n";

    for (const auto& g : bench_grids) {
        ADESolver ds_solver(base_params, g.NY, g.NY, g.NT,
                            SweepVariant::DoubleSweep, ExerciseStyle::AmericanSpot);
        ADESolver s_solver(base_params, g.NY, g.NY, g.NT,
                           SweepVariant::Saulyev, ExerciseStyle::AmericanSpot);

        // Best of 3 runs as specified in Section 7.6
        ADEResult ds_res = ds_solver.solve();
        ADEResult s_res  = s_solver.solve();
        for (int r = 0; r < 2; ++r) {
            ds_res.runtime_sec = std::min(ds_res.runtime_sec, ds_solver.solve().runtime_sec);
            s_res.runtime_sec  = std::min(s_res.runtime_sec,  s_solver.solve().runtime_sec);
        }

        const std::string grid_str = std::to_string(g.NY) + " x " + std::to_string(g.NT);
        std::cout << std::left << std::setw(16) << grid_str
                  << std::right << std::setw(14) << ds_res.price
                  << std::setw(14) << s_res.price
                  << std::scientific << std::setprecision(2)
                  << std::setw(14) << std::abs(s_res.price - ds_res.price)
                  << std::fixed << std::setprecision(4)
                  << std::setw(12) << ds_res.runtime_sec
                  << std::setw(12) << s_res.runtime_sec << "\n";
        std::cout << std::setprecision(6);
    }
    std::cout << "\n";

    // -------------------------------------------------------------------------
    // 5. Chapter 8: Practical Trader Extensions (Futures Obstacle, Seasonality, Greeks)
    // -------------------------------------------------------------------------
    std::cout << "[5] CHAPTER 8 (TABLE 5): TRADER EXTENSIONS (ATM S0 = K = 100, delta0 = 0.1)\n";
    std::cout << "----------------------------------------------------------------------------\n";
    SchwartzParameters atm_params = base_params;
    atm_params.S0 = 100.0;

    ADESolver atm_spot_solver(atm_params, 200, 200, 600,
                              SweepVariant::DoubleSweep, ExerciseStyle::AmericanSpot, false);
    const ADEResult atm_res = atm_spot_solver.solve();

    ADESolver atm_seasonal_solver(atm_params, 200, 200, 600,
                                  SweepVariant::DoubleSweep, ExerciseStyle::AmericanSpot, true);
    const ADEResult seasonal_res = atm_seasonal_solver.solve();

    ADESolver atm_futures_solver(atm_params, 200, 200, 600,
                                 SweepVariant::DoubleSweep, ExerciseStyle::AmericanFutures, false, 0.25);
    const ADEResult futures_obs_res = atm_futures_solver.solve();

    const auto [h1, h2] = SchwartzAnalytical::solveTwoFuturesHedge(
        atm_res.delta_S, atm_res.delta_yield, 0.25, 0.75, atm_params);

    std::cout << "  ATM American Spot Put Price         : " << atm_res.price << "\n";
    std::cout << "  ATM Spot Delta (Delta_S, Eq. 52)    : " << atm_res.delta_S << "\n";
    std::cout << "  ATM Spot Gamma (Gamma_S, Eq. 52)    : " << atm_res.gamma_S << "\n";
    std::cout << "  ATM Yield Delta (Delta_delta, Eq.52): " << atm_res.delta_yield << "\n";
    std::cout << "  Seasonal Parameter Put Price (Eq.51): " << seasonal_res.price
              << "  (+" << std::setprecision(2)
              << 100.0 * (seasonal_res.price - atm_res.price) / atm_res.price << "%)\n"
              << std::setprecision(6);
    std::cout << "  Futures-Obstacle Put (theta = 0.25) : " << futures_obs_res.price << "\n";
    std::cout << "  2-Futures Hedge Ratios (Eq. 54)     : h1 (T1=0.25) = " << h1
              << ", h2 (T2=0.75) = " << h2 << "\n";
    std::cout << "============================================================================\n";

    return 0;
}
