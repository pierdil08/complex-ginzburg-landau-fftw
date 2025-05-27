# Complex Ginzburg–Landau Equation Solver (FFTW, C)

This project solves the Complex Ginzburg–Landau equation in 2D using the pseudospectral method and a 4th-order Runge–Kutta scheme. It was written in C and executed on Northwestern's high-performance computing cluster.

## 🔬 Equation
The equation takes the form:

∂A/∂t = A + (1 + i*c1)∇²A - (1 - i*c3)|A|²A

## ⚙️ Tools Used
- Language: C
- Libraries: FFTW3
- HPC system: Lagrange @ Northwestern

## 🚀 Features
- FFT-based spectral solver
- RK4 time integration
- Performance tested with varying N

## 📈 Results
Simulation shows the emergence of spiral wave patterns and phase separation.

## 🛠️ How to Run
```bash
make
./cgl 128 1.5 0.25 100000 12345
