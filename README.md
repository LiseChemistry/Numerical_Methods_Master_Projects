# 📚 **Numerical Methods – Master Projects**

[![Fortran](https://img.shields.io/badge/code-Fortran-blue?style=flat-square&logo=fortran)](https://en.wikipedia.org/wiki/Fortran)
[![Python](https://img.shields.io/badge/code-Python-yellow?style=flat-square&logo=python)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-MIT-green?style=flat-square)](LICENSE)

## 📘 Overview

This repository contains numerical methods projects completed during my Master’s program. The projects involve scientific programming in **Fortran** and **Python**, focusing on:

🔹 **Project 1: Interpolation / Extrapolation**  
&nbsp;&nbsp;&nbsp;&nbsp;• Potential of O₂  
&nbsp;&nbsp;&nbsp;&nbsp;• Potential of the C₂H₅ radical  

🔹 **Project 2: Integration & Differentiation**  
&nbsp;&nbsp;&nbsp;&nbsp;• Average value  
&nbsp;&nbsp;&nbsp;&nbsp;• Numerical derivatives of VO₂(R)  
&nbsp;&nbsp;&nbsp;&nbsp;• Integration  

🔹 **Project 3: Model Fitting**  
&nbsp;&nbsp;&nbsp;&nbsp;• Spectroscopic Hamiltonian analysis  

> 🧠 Developed using **Fortran** for high-performance numerical tasks and **Python** for data fitting and visualization.

---

### Project I: Interpolation/Extrapolation (Fortran)  

#### 📌Exercise 1: Potential of O₂ 

We consider the internuclear potential of the O₂ molecule:

```math
V(R) = -D \left\{1 + c_1(R - R^*) + c_2(R - R^*)^2 + c_3(R - R^*)^3 \right\} e^{-\alpha(R - R^*)}
```

Where:

- D = 3.8886 eV  
- R* = 2.2818 a₀
- α = 3.3522498 a₀⁻¹  
- c₁ = 3.6445906 a₀⁻¹  
- c₂ = 3.9281238 a₀⁻²  
- c₃ = 2.0986689 a₀⁻³

This is a model adjusted to experimental spectroscopic data.

We want to interpolate the potential over the interval [1.8, 7]a₀ using these points:

```math
\{ R_q \} = \{1.8, 2.1, 2.4, 2.8, 3.3, 4.0, 5.0, 7.0\}
```
##### 🧩Interpolation Methods:
- Lagrange Polynomial ( L(R) )
- Polynomial in powers of ( 1/R )
- Rational Polynomial (`RatInt.f`)
- Cubic Splines (`Spline.f`)

Plot the interpolation error:

```math
|V(R_i) - V_{\text{Int}}(R_i)|
```
on a **logarithmic scale** for R in [1.8, 10], and compare the behavior of the interpolation and extrapolation methods.

---
 
#### 📌Exercise 2: Potential of the C₂H₅ Radical

We use the model potential:

```math
V(\phi) = e^{-\alpha \cos(6\phi)}, \quad \alpha = 0.2
```

This potential has the same symmetry properties as the exact one.

##### 🔧Tasks:
- Build an interpolation scheme V̄(φ) using `Minv.f`.
- Determine the minimum number of interpolation points {φ_q} such that the standard deviation:

```math
\phi_{q+\frac{1}{2}} = \frac{1}{2}(\phi_{q+1} + \phi_q)
```

```math
\sigma = \sqrt{ \frac{1}{N - 1} \sum_{q=1}^{N-1} \left[ V(\phi_{q+1/2}) - \bar{V}(\phi_{q+1/2}) \right]^2 } < 10^{-5}
```

---

### Project II: Integration – Differentiation (Fortran)
  
#### 📌Exercise 1: Average Value

We want to compute the average distance of the electron to the nucleus in the hydrogen atom’s 3s orbital:

```math
\langle r \rangle_{3s} = \langle \phi_{3s} | r | \phi_{3s} \rangle = 4\pi \int_0^\infty \phi_{3s}^2(\vec{r}) r^3 dr
```

With:

```math
\phi_{3s}(\vec{r}) = \frac{1}{\sqrt{4\pi}} \cdot \frac{2}{81\sqrt{3}} \left(\frac{Z}{a_0}\right)^{3/2} (27 - 18\rho + 2\rho^2) e^{-\rho/3}, \quad \rho = \frac{Zr}{a_0}
```

##### 🧩Methods:
- Gauss-Laguerre Quadrature (`Laguerre_Quad.c`)  
- Monte Carlo Method (with same RNG seed)  
- (Optional) Monte Carlo with Importance Sampling

Compare with the analytical value:

```math
\langle r \rangle_{ns} = \frac{n^2 a_0}{Z} \left(1 + \frac{1}{2}\left(1 - \frac{\ell(\ell + 1)}{n^2} \right)\right)
```
---

#### 📌Exercise 2: Numerical Derivatives of the VO₂ Potential

Reuse the O₂ potential from Project I and compute its first and second derivatives at:

```math
R = R^*
```

Study the effect of:
- Number of points used
- Spacing of points

Compare with the analytical derivatives.

---

#### 📌Exercise 3: Laser Pulse Integration

Consider a laser pulse defined by the function:

```math
A(t) = A_0 \sin(\omega t) e^{-\alpha(t - T)^2}
```

Where:

- ω = 2  
- T = 5π  
- α = 0.05

##### 🔧Tasks:
- Compute the **fluence** Φ:

```math
\Phi = \int_0^{2T} A^2(t) dt
```

##### 🔧Using:
- Trapezoidal method
- Simpson's rule

Increase the number of integration points until the **relative convergence** reaches 10⁻⁸.

---

### Project III: Model Fitting (Python)

#### 📌Exercise 1: Study of a Spectroscopic Hamiltonian

We consider the vibrational spectroscopic Hamiltonian:

```math
E_{\text{sp}}(n) = T_0 + \sum_{\alpha=1}^{6} \omega_\alpha \left(n_\alpha + \frac{1}{2} \right) + \sum_{\alpha \leq \beta} x_{\alpha\beta} \left(n_\alpha + \frac{1}{2} \right)\left(n_\beta + \frac{1}{2} \right) + \sum_{\alpha \leq \beta \leq \gamma} y_{\alpha\beta\gamma} \left(n_\alpha + \frac{1}{2} \right)\left(n_\beta + \frac{1}{2} \right)\left(n_\gamma + \frac{1}{2} \right)
```
Adjust the parameters **T₀**, {**ω<sub>α</sub>**}, {**x<sub>αβ</sub>**}, and {**y<sub>αβγ</sub>**} using the **least squares method**. These parameters can be determined based on the list {**E<sub>exp</sub>(n)**} (see the file `HFCO_exp.dat`), which contains the first **150 calculated vibrational levels** of the HFCO molecule.

##### 🔧Tasks:
- Fit Only ω<sub>α</sub> and x<sub>αβ</sub> and use a regular matrix inversion method to perform the least squares fit:

```python
numpy.linalg.inv
```

- Include y<sub>αβγ</sub> terms. Some of these parameters cannot be determined.
if the state **(n₁ = 1, n₂ = 1, n₃ = 1, n₄ = n₅ = n₆ = 0)** does not appear in the list, then **y₁₂₃** cannot be fitted.

- Test with **Singular Value Decomposition**:

```python
numpy.linalg.svd
```

##### Print:
- The fitted parameters
- The residuals $E_{\text{exp}}(n) - E_{\text{sp}}(n)$ for all levels

---

### License
The package is licensed under the MIT License.
