# 🚀 PROJECT ENTERPRISE // Orbital Mission Architect & Flight Dynamics Simulator

[![Language: C11](https://img.shields.io/badge/Core_Engine-C11-00599C?style=for-the-badge&logo=c&logoColor=white)](https://en.wikipedia.org/wiki/C11_(C_standard_revision))
[![Frontend: HTML5 / JS / Plotly](https://img.shields.io/badge/Web_UI-Plotly.js_3D-00e5ff?style=for-the-badge&logo=javascript&logoColor=black)](https://plotly.com/javascript/)
[![Standard: WGS84 / ESA-NASA](https://img.shields.io/badge/Astrodynamics-WGS84_%7C_RK4_%7C_Lambert-00e676?style=for-the-badge)](#)
[![Live Demo](https://img.shields.io/badge/Live_Platform-GitHub_Pages-ffab00?style=for-the-badge&logo=github)](https://mattechiane.github.io/Enterprise/)

**Project ENTERPRISE** is a high-precision astrodynamics simulator, global maneuver optimizer, and interactive 3D mission control platform inspired by industry-standard flight dynamics suites (**AGI STK** and **NASA GMAT**).

Built around a fast, zero-dependency **C11 numerical core** coupled with a real-time **3D WebGL/Plotly Mission Control interface**, ENTERPRISE enables end-to-end orbital mission design—from initial parking orbit assessment and optimal inter-orbital transfers to low-thrust electric spirals, terminal proximity rendezvous, and stochastic Monte Carlo dispersion analysis.

🌐 **Live Interactive Web Platform:** [**https://mattechiane.github.io/Enterprise/**](https://mattechiane.github.io/Enterprise/)

---

## ✨ Key Engineering Capabilities

### 1. 🌍 High-Fidelity Orbital Propagation & Perturbations
* **4th-Order Runge-Kutta (RK4) Numerical Integrator:** Propagates full 6-DOF Cartesian state vectors $\mathbf{x} = [\mathbf{r}^T, \mathbf{v}^T, m]^T$ in the Earth-Centered Inertial (ECI) frame.
* **WGS84 Oblate Earth Reference ($J_2$ Zonal Harmonics):** Models nodal regression ($\dot{\Omega}$) and apsidal rotation ($\dot{\omega}$) induced by Earth's equatorial bulge ($R_{\text{eq}} = 6378.137\text{ km}$, $f = 1/298.257$).
* **Exponential Atmospheric Drag & Station-Keeping Budget:** Evaluates orbital decay ($\text{m/day}$) at perigee and computes the annual $\Delta V$ and propellant mass required to counteract atmospheric drag.
* **Bidirectional State Resolver & Physical Guardrails:** Accepts both **Parametric Keplerian** $(a, e, i, \Omega, \omega, \theta)$ and **Cartesian ECI** $(\mathbf{r}, \mathbf{v})$ inputs, enforcing ESA/NASA closed-orbit safety envelopes ($h_p \ge 150\text{ km}$, $0 \le e < 0.95$, apogee within Earth's Hill sphere).

### 2. 🎯 Global Multi-Strategy Maneuver Optimizer & Propulsion Database
Automatically benchmarks multiple orbital transfer architectures in 3D space (accounting for full dihedral plane changes across both inclination $\Delta i$ and right ascension of the ascending node $\Delta \Omega$, as well as orbital eccentricities $e_1, e_2$):
* **3D Optimal Split-Plane Hohmann Transfer (2-Burn):** Optimally partitions the dihedral plane change angle between perigee departure ($\alpha \Delta \theta_{\text{plane}}$) and apogee insertion ($(1-\alpha)\Delta \theta_{\text{plane}}$) to minimize total $\Delta V$.
* **Optimal 3-Burn Bi-Elliptic Transfer:** Evaluates high-apogee intermediate ellipses to exploit low apogee velocities during large plane changes or high radius ratios.
* **Universal Variable Lambert Solver & Trade-Off Curves:** Solves the two-point boundary value problem $(\mathbf{r}_1, \mathbf{r}_2, \Delta t)$ using universal variables and Stumpff functions $C(z), S(z)$, generating interactive **$\Delta V$ & Propellant vs. Time-of-Flight (Porkchop Slice)** trade-off plots.
* **Continuous Low-Thrust Electric Spiral (Edelbaum Theory):** Simulates multi-revolution continuous orbit raising and plane changing for high-$I_{sp}$ electric thrusters.
* **Built-in Propulsion Catalog:**
  * *Solid Composite:* HTPB/AP Star-48B ($I_{sp} = 292\text{ s}$)
  * *Liquid Monopropellant:* Hydrazine $\text{N}_2\text{H}_4$ ($I_{sp} = 235\text{ s}$)
  * *Liquid Storable Bipropellant:* MMH / NTO ($I_{sp} = 316\text{ s}$)
  * *Liquid Methalox:* $\text{LOX} / \text{LCH}_4$ ($I_{sp} = 365\text{ s}$)
  * *Liquid Cryogenic Hydrolox:* $\text{LOX} / \text{LH}_2$ ($I_{sp} = 450\text{ s}$)
  * *Electric Hall-Effect Thruster (HET):* Xenon ($I_{sp} = 1650\text{ s}$)
  * *Electric Gridded Ion Engine:* NEXT-C Xenon ($I_{sp} = 3200\text{ s}$)

### 3. 🛰️ Onboard Subsystems, Ground Radar & Proximity Operations
* **Topocentric Ground Station Access:** Computes real-time azimuth, elevation ($\text{El}$), slant range ($\rho$), and 24-hour contact windows over ground stations (default: *Milano-Bovisa*, $\text{Lat } 45.502^\circ\text{ N}, \text{Lon } 9.156^\circ\text{ E}$, mask $\ge 10^\circ$).
* **Electrical Power Subsystem (EPS) Simulation:** Tracks cylindrical Earth shadow (umbra) transits, Triple-Junction GaAs solar array generation ($W$), and Li-Ion battery State-of-Charge ($\text{SoC } \%$) across sunlight and eclipse cycles.
* **Terminal Rendezvous & Proximity Operations (Clohessy-Wiltshire / Hill Frame):** Analytically propagates relative motion in the Local-Vertical Local-Horizontal (LVLH: R-bar, V-bar, H-bar) frame for two-impulse soft-docking trajectories.
* **Monte Carlo Navigation Uncertainty & TCM Budget:** Performs 100-run Box-Muller stochastic simulations of thrust magnitude and attitude pointing errors, rendering $3\sigma$ (99.7%) arrival dispersion ellipsoids and Mid-Course Trajectory Correction Maneuver ($\Delta V_{\text{TCM}}$) probability distributions.

---

## 🖥️ Interactive Mission Control Dashboard Architecture

The web interface (`index.html`) is structured as a vertical engineering workflow:
1. **3D WGS84 Earth Ellipsoid & Trajectory View:** Interactive 3D inertial globe with real continental coastlines, animated spacecraft marker, and highlighted 3D burn execution points ($\Delta \mathbf{V}_1, \Delta \mathbf{V}_2$).
2. **Mission Engineering Data Sheet:** Complete breakdown of optimal maneuver selection, ECI/Geodetic burn coordinates, RTN thrust unit vectors, propellant mass budget, and initial vs. target regime comparison.
3. **Lambert $\Delta V$ & Propellant Trade-Off Curve:** Interactive cost vs. time-of-flight curve allowing one-click mission retargeting.
4. **2D World Ground Track (WGS84):** Sub-satellite geodetic trace with continuous Earth rotation shift (GMST), burn projections, and ground station horizon mask.
5. **Onboard Subsystems Time-History Plots:** Synchronized telemetry charts for altitude/velocity kinematics, ground radar elevation passes, and solar power / battery SoC.
6. **LVLH Terminal Rendezvous & Monte Carlo Dispersion Analysis:** Relative V-bar/R-bar docking corridor and $3\sigma$ stochastic error analysis.

---

## 🛠️ Project Structure & Build Instructions

```text
Enterprise/
├── CMakeLists.txt          # CMake build configuration (C11 standard)
├── orbit_types.h           # Core data structures, physical constants & propellant catalog
├── astrodynamics.h / .c    # Vector math, Kepler's equation (Newton-Raphson) & frame transforms
├── propagator.h / .c       # RK4 numerical integrator, J2 gravity, atmospheric drag & eclipses
├── mission_analysis.h / .c # Global maneuver optimizer, Lambert solver, Edelbaum spiral, CW & Monte Carlo
├── telemetry.h / .c        # ECI/ECEF/WGS84 geodetic converters & JSON/CSV telemetry exporter
├── main.c                  # CLI mission execution entry point
├── mission_report.json     # Generated mission telemetry payload for the Web UI
└── index.html              # Self-contained 3D Web Mission Control & in-browser solver
```

### Building and Running the C11 Core Engine
Requirements: A C11-compliant compiler (`gcc`, `clang`, or `MSVC`) and `CMake >= 3.20`.

```bash
# 1. Clone the repository
git clone https://github.com/mattechiane/Enterprise.git
cd Enterprise

# 2. Configure and build via CMake
mkdir build && cd build
cmake ..
cmake --build .

# 3. Run the astrodynamics engine to generate mission_report.json
./Enterprise
```

### Launching the Web Platform
* **Cloud Version:** Visit [**https://mattechiane.github.io/Enterprise/**](https://mattechiane.github.io/Enterprise/) directly in any modern browser.
* **Local Development:** Open `index.html` via CLion's built-in web preview or any local HTTP server (e.g., `python -m http.server 8000`) to automatically fetch `mission_report.json` or design custom missions using the **Configure New Mission** modal.

---

## 📜 Author & License
Developed by **Matteo Chianese** — Aerospace Engineering Project (*Project ENTERPRISE*).