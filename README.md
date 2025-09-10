# Runway Scheduling: GA vs ILP

This project compares **Genetic Algorithms (GA)** and **Integer Linear Programming (ILP)** for solving a **runway scheduling optimization problem**.  
The goal is to schedule aircraft arrivals and departures while respecting minimum separation constraints, minimizing delays, and prioritizing important flights.

---

## 📌 Features
- **Synthetic Data Generation**: Creates randomized arrival/departure scenarios.
- **Genetic Algorithm (GA)**:
  - Multi-objective optimization (delay + priority penalty).
  - Uses stochastic delays to simulate real-world uncertainties.
- **Integer Linear Programming (ILP)**:
  - Deterministic optimization of scheduling with separation constraints.
  - Solves using MATLAB’s `intlinprog`.
- **Performance Analysis**:
  - Runtime comparison of GA vs ILP.
  - Objective value comparison of GA vs ILP.
  - Visualization of scaling with increasing aircraft.

---

## 📊 Results Summary
- **Runtime Scaling**:
  - GA scales more smoothly with increasing aircraft.
  - ILP runtime grows much faster, becoming computationally expensive for larger instances.
- **Objective Value**:
  - GA provides better (lower) average objective values.
  - ILP yields higher values as it optimizes deterministically without stochastic handling.

Plots from experiments:

1. **Runtime Scaling**  
   ![Runtime Scaling](figures/runtime_scaling.png)

2. **Objective Comparison**  
   ![Objective Comparison](figures/objective_comparison.png)

---

## 🛠️ Requirements
- MATLAB (R2021a or later recommended)
- Optimization Toolbox
- Global Optimization Toolbox

---

## 🚀 Usage
1. Clone this repository:
   ```bash
   git clone https://github.com/your-username/runway-scheduling-GA-vs-ILP.git
   cd runway-scheduling-GA-vs-ILP
