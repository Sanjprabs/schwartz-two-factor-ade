# Pricing American Commodity Options Under the Schwartz Two-Factor Model Using the Alternating Direction Explicit (ADE) Method

[![CI Build & Tests](https://github.com/Sanjprabs/schwartz-two-factor-ade/actions/workflows/ci.yml/badge.svg)](https://github.com/Sanjprabs/schwartz-two-factor-ade/actions)
[![C++17](https://img.shields.io/badge/C%2B%2B-17-blue.svg)](https://isocpp.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![arXiv](https://img.shields.io/badge/arXiv-2610.21145-b31b1b.svg)](https://arxiv.org/abs/2610.21145)
[![SSRN](https://img.shields.io/badge/SSRN-7604358-navy.svg)](https://ssrn.com//sol3/papers.cfm?abstract_id=7604358)

**Authors:** Sanjith Prabhu & Dr. Daniel J. Duffy  
**Paper:** *Pricing American Commodity Options Under the Schwartz Two-Factor Model Using the Alternating Direction Explicit (ADE) Method* (`arXiv:2610.21145` | `SSRN: 7604358`)

---

## Overview

This repository contains the official C++17 numerical implementation accompanying our working paper on valuing American-style commodity options under the **Schwartz (1997) two-factor model** (stochastic spot price $S_t$ and mean-reverting instantaneous convenience yield $\delta_t$).

Traditional implicit finite difference schemes (such as Crank–Nicolson or ADI) require solving coupled linear systems or Projected Successive Over-Relaxation (PSOR) iterations at every time step to enforce the American early-exercise constraint. We implement both the **Barakat–Clark (1966) double-sweep** and **Saul'yev (1964) single-sweep alternating** variants of the **Alternating Direction Explicit (ADE)** method. By updating grid nodes sequentially and incorporating current-time-level neighbors, ADE achieves unconditional stability for diffusion with strictly explicit $O(N_X N_Y N_T)$ complexity, executing **210x to 390x faster** than Longstaff–Schwartz Least Squares Monte Carlo (LSMC) at equivalent accuracy ($\approx 2 \times 10^{-3}$).

---

## Mathematical Formulation

### 1. Risk-Neutral Joint Dynamics
Under the risk-neutral measure $\mathbb{Q}$, the commodity spot price $S_t$ and instantaneous convenience yield $\delta_t$ evolve according to:

$$dS_t = (r - \delta_t)S_t\,dt + \sigma_1 S_t\,dz_{1,t}$$

$$d\delta_t = \kappa(\hat{\alpha} - \delta_t)\,dt + \sigma_2\,dz_{2,t}, \qquad \hat{\alpha} = \alpha - \frac{\lambda}{\kappa}$$

where $dz_{1,t}\,dz_{2,t} = \rho\,dt$, $\kappa > 0$ is the speed of mean reversion, and $\lambda$ is the market price of convenience yield risk.

### 2. Log-Transformed Governing 2D PDE & Variational Inequality
Applying the logarithmic spot transformation $x = \ln S$, keeping $y = \delta$, and setting time-to-maturity $\tau = T - t$, the state-dependent diffusion coefficients reduce to constants:

$$V_\tau = \mathcal{L}V = \frac{1}{2}\sigma_1^2 V_{xx} + \rho\sigma_1\sigma_2 V_{xy} + \frac{1}{2}\sigma_2^2 V_{yy} + \left(r - y - \frac{1}{2}\sigma_1^2\right)V_x + \kappa(\hat{\alpha} - y)V_y - rV$$

For an American put with payoff obstacle $\psi(x, y, \tau)$, the free-boundary problem satisfies the variational inequality:

$$V_\tau - \mathcal{L}V \ge 0, \qquad V - \psi \ge 0, \qquad (V_\tau - \mathcal{L}V)(V - \psi) = 0, \qquad V(\cdot, 0) = \psi$$

### 3. Fichera Boundary Treatment & Computational Domain
The computational mesh is centered at the initial state $(x_0, y_0) = (\ln S_0, \delta_0)$ on $[x_0 - 1.2, x_0 + 1.2] \times [\delta_0 - 0.4, \delta_0 + 0.4]$ ($\approx \pm 10.7\sigma$ in log-spot and $\pm 6\sigma$ in the stationary Ornstein–Uhlenbeck distribution):
* **Lower Spot Boundary ($S \to 0$):** By Fichera theory, the Fichera function vanishes ($\mathcal{F}\vert{}_{S=0} = 0$), rendering the boundary characteristic; Dirichlet asymptotics are applied ($V \to K e^{-r\tau}$ for European puts; $V = \psi$ for American puts).
* **Upper Spot Boundary ($S \to \infty$):** Homogeneous Dirichlet condition $V(x_{\max}, y, \tau) = 0$.
* **Convenience Yield Boundaries ($y_{\min}, y_{\max}$):** Linear extrapolation ($V_{yy} = 0 \implies u_{i,0} = 2u_{i,1} - u_{i,2}$ and $u_{i,N_Y} = 2u_{i,N_Y-1} - u_{i,N_Y-2}$).

---

## Empirical Results (WTI Crude Oil Calibration)

All benchmarks use the Schwartz (1997) Kalman-filtered WTI crude oil parameters:
$$\kappa = 1.49, \quad \hat{\alpha} = 0.0161, \quad \sigma_1 = 0.158, \quad \sigma_2 = 0.115, \quad \rho = 0.189, \quad r = 0.05, \quad K = 100, \quad T = 0.5$$

### 1. European Put Validation ($S_0 = 90, \delta_0 = 0.1$)
| Pricing Method | Computed Value | Benchmark Error / SE |
| :--- | :---: | :---: |
| **Schwartz (1997) Closed-Form Analytical** | `11.601861` | Exact Reference |
| **Standard Monte Carlo (250,000 paths)** | `11.602859` | $\text{SE} = 0.016473$ |
| **2D Double-Sweep ADE ($400 \times 400 \times 1200$)** | `11.601279` | $\vert{}V_{\text{ADE}} - V_{\text{Exact}}\vert{} = 5.82 \times 10^{-4}$ |

### 2. Double-Sweep ADE Grid Convergence (American Put, $S_0 = 90, \delta_0 = 0.1$)
| $N_Y \backslash N_T$ | 100 | 300 | 600 | 1200 |
| :---: | :---: | :---: | :---: | :---: |
| **50** | 11.8684 | 11.8711 | 11.8715 | 11.8716 |
| **100** | 11.8821 | 11.8898 | 11.8907 | 11.8909 |
| **200** | 11.8598 | 11.8880 | 11.8908 | 11.8916 |
| **400** | 11.7705 | 11.8788 | 11.8894 | **11.8921** |

> **Consistency Note:** Because the ADE truncation error scales with $O(\Delta\tau / h)$ rather than $O(\Delta\tau / h^2)$, temporal resolution $N_T$ must be refined alongside spatial resolution $N_Y$.

### 3. Market Scenario Analysis vs. Longstaff–Schwartz LSMC (250,000 Paths)
| $S_0$ | $\delta_0$ | Market Regime | ADE Price | LSMC Price | Absolute Diff |
| :---: | :---: | :--- | :---: | :---: | :---: |
| **90** | **0.1** | In-the-Money (Base) | `11.892050` | `11.888911` | `0.003139` |
| **100** | **0.1** | At-the-Money (ATM) | `5.020935` | `5.024249` | `0.003314` |
| **90** | **0.3** | ITM, High Convenience Yield | `16.952677` | `16.989201` | `0.036524` |

### 4. Computational Efficiency: Double-Sweep ADE vs. Saul'yev ADE vs. LSMC
| Method | Grid / Path Configuration | $\vert{}V - V_{\text{ref}}\vert{}$ or SE | Runtime (s) | Speedup vs. Equal-Accuracy LSMC |
| :--- | :--- | :---: | :---: | :---: |
| **Double-Sweep ADE** | $N_Y = 50, N_T = 300$ | $2.2 \times 10^{-2}$ | `0.093 s` | — |
| **Saul'yev ADE** | $N_Y = 50, N_T = 300$ | $2.2 \times 10^{-2}$ | `0.050 s` | — |
| **Double-Sweep ADE** | $N_Y = 100, N_T = 600$ | $2.1 \times 10^{-3}$ | `0.749 s` | **~210x faster** |
| **Saul'yev ADE** | $N_Y = 100, N_T = 600$ | $2.1 \times 10^{-3}$ | `0.410 s` | **~390x faster** |
| **Double-Sweep ADE** | $N_Y = 200, N_T = 600$ | $1.9 \times 10^{-3}$ | `2.986 s` | — |
| **Saul'yev ADE** | $N_Y = 200, N_T = 600$ | $1.9 \times 10^{-3}$ | `1.655 s` | — |
| **LSMC (50 dates)** | $M = 10,000$ paths | $\pm 0.0708$ | `0.088 s` | — |
| **LSMC (50 dates)** | $M = 50,000$ paths | $\pm 0.0327$ | `0.632 s` | — |
| **LSMC (50 dates)** | $M = 250,000$ paths | $\pm 0.0144$ | `3.365 s` | Baseline ($1.2 \times 10^7$ paths $\approx 160\text{s}$ needed for $2.1 \times 10^{-3}$) |

---

## Trader Desk Extensions (Chapter 8)
Following feedback from an agricultural commodities trader, the solver includes three desk-oriented extensions:
1. **Exercise Obstacle on a Futures Contract:** Evaluates early exercise against a futures contract $F(S, \delta; \tau + \theta)$ with residual delivery lag $\theta = T_F - T$:
   $$\psi_F(x, \delta, \tau) = \max\left\{K - e^{x - \delta B(\tau + \theta) + A(\tau + \theta)}, 0\right\}$$
2. **Seasonal Parameter Dynamics:** Time-dependent convenience yield equilibrium $\hat{\alpha}(\tau)$ and volatilities $\sigma_1(\tau), \sigma_2(\tau)$ capturing planting and harvest cycles (~2–3% option value impact).
3. **Two-Factor Greeks & Two-Futures Hedging System:** Computes $\Delta_S = V_x / S$, $\Gamma_S = (V_{xx} - V_x)/S^2$, and $\Delta_\delta = V_y$ (yielding $\Delta_S = -0.520$ and $\Delta_\delta = 16.814$ at ATM), and solves the $2 \times 2$ linear hedging system across nearby ($F_1$) and deferred ($F_2$) futures maturities.

---

## Repository Structure

```text
schwartz-two-factor-ade/
├── CMakeLists.txt               # Cross-platform CMake build configuration (-O3)
├── README.md                    # Documentation, mathematical formulation, and results
├── .github/workflows/ci.yml     # Automated GitHub Actions build and unit test pipeline
├── include/
│   ├── Matrix2D.hpp             # Contiguous 2D finite-difference grid container
│   ├── SchwartzModel.hpp        # Model parameters, closed-form futures/European, 2-futures hedge
│   ├── ADESolver.hpp            # Double-Sweep and Saul'yev 2D ADE PDE solver interface
│   └── MonteCarloLSMC.hpp       # European MC & Longstaff-Schwartz LSMC engine interface
├── src/
│   ├── ADESolver.cpp            # ADE directional sweeps, Fichera boundaries, Greeks
│   ├── MonteCarloLSMC.cpp       # Correlated path generator & weighted Laguerre OLS regression
│   └── main.cpp                 # Reproduces all tables and benchmarks from the paper
├── tests/
│   └── unit_tests.cpp           # Automated verification suite
└── paper/                       # Working paper PDF and LaTeX manuscript
```

---

## Build & Run Instructions

### Prerequisites
* A C++17 compatible compiler (`g++ >= 9`, `clang++ >= 10`, or `MSVC >= 2019`)
* `CMake >= 3.14`

### Compiling and Running with CMake
```bash
git clone [https://github.com/Sanjprabs/schwartz-two-factor-ade.git](https://github.com/Sanjprabs/schwartz-two-factor-ade.git)
cd schwartz-two-factor-ade
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build . -j4

# Run unit tests
./schwartz_unit_tests

# Run full paper replication and benchmark suite
./schwartz_ade_suite
```

---

## Citation

```bibtex
@techreport{prabhu2026schwartz,
  title   = {Pricing American Commodity Options Under the Schwartz Two-Factor Model Using the Alternating Direction Explicit (ADE) Method},
  author  = {Prabhu, Sanjith and Duffy, Daniel J.},
  year    = {2026},
  note    = {Working Paper, arXiv:2610.21145 and SSRN: 7604358}
}
```
