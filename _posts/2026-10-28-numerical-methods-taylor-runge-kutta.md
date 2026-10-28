---
title: "Numerical Methods : Taylor and Runge-Kutta Methods"
date: 2026-10-28 18:20:00 +0100
categories:
  - Numerical Analysis
tags:
  - Taylor-Methods
  - Runge-Kutta
  - Heun-Method
  - Accuracy
math: false
---

## Introduction

Forward Euler는 계산이 간단하지만 1차 정확도를 가진다. 이번 글에서는 Taylor Method에서 고차 근사가 어떻게 만들어지는지, 함수의 고차 미분을 직접 계산하지 않고도 RK2와 RK4가 높은 정확도를 달성하는 이유를 살펴본다. `Lec-2_3_scheme.pdf`의 5–9쪽을 중심으로 정리하고, Lab Handout Exercise 4의 Heun 방법을 직접 유도한다.

## 1. Higher-Order Taylor Methods

정확한 해가 충분히 매끄럽다면 `tₙ` 주변에서 Taylor 전개를 적용할 수 있다.

![Taylor expansion](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7Dy%28t_n%2Bh%29%3Dy%28t_n%29%2Bhy%5Cprime%28t_n%29%2B%5Cfrac%7Bh%5E2%7D%7B2%7Dy%5Cprime%5Cprime%28t_n%29%2B%5Ccdots%2B%5Cfrac%7Bh%5Em%7D%7Bm%21%7Dy%5E%7B%28m%29%7D%28t_n%29%2BO%28h%5E%7Bm%2B1%7D%29)

미분방정식 `y′=f(t,y)`를 대입하면 첫 도함수는 이미 알고 있다. 2차 미분은 합성함수 미분을 통해 구한다.

![Second derivative from ODE](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7Dy%5Cprime%5Cprime%3Df_t%2Bf_yf)

**m차 Taylor 방법**은 이 전개에서 m차까지 남기고 나머지항을 버리는 방식이다. Forward Euler는 m=1인 특수한 사례이다. 다만 고차 도함수를 반복 계산해야 하므로 실제로는 구현 부담이 크다.

### 1.1 Lecture Example — Third-Order Taylor Method

강의의 `y′ = y − t³ + t + 1`, `y(0)=0.5`에서 `f=y−t³+t+1`이다.

![First total derivative](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7Df%5Cprime%3Df_t%2Bf_y%20f%3Dy-t%5E3-3t%5E2%2Bt%2B2)

![Second total derivative](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7Df%5Cprime%5Cprime%3Dy-t%5E3-3t%5E2-5t%2B2)

따라서 `T₃ = f + (h/2)f′ + (h²/6)f″`이며 `yₙ₊₁ = yₙ + h T₃(tₙ,yₙ)`으로 업데이트한다.

## 2. Runge-Kutta Methods

Runge–Kutta 방법은 Taylor와 유사한 정확도를 얻되 고차 도함수 대신 **여러 지점에서 계산한 기울기의 가중평균**을 사용한다.

![Runge Kutta weighted slopes](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7Dy_%7Bn%2B1%7D%3Dy_n%2Bh%5Csum_%7Bi%3D1%7D%5E%7Bm%7Dw_iK_i)

기울기를 어느 지점에서 평가하고 어떤 비율로 섞을지 선택하여 Taylor 전개 계수와 일치시키는 것이 핵심이다.

### 2.1 Second-Order Runge-Kutta (RK2)

일반적인 2단계 RK 방법에서 `K₁=f(tₙ,yₙ)`, `K₂=f(tₙ+αh,yₙ+βhK₁)`라 두면 2차 정확도를 위해 `w₁+w₂=1`, `w₂α=1/2`, `w₂β=1/2`가 필요하다.

### 2.2 Heun's Method — Lab Handout Exercise 4

`w₁=w₂=1/2`, `α=β=1`로 두면 Heun의 Predictor–Corrector 공식이 된다.

![Heun predictor](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cwidetilde%20y_%7Bn%2B1%7D%3Dy_n%2Bhf%28t_n%2Cy_n%29)

![Heun corrector](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7Dy_%7Bn%2B1%7D%3Dy_n%2B%5Cfrac%20h2%5Cleft%5Bf%28t_n%2Cy_n%29%2Bf%28t_%7Bn%2B1%7D%2C%5Cwidetilde%20y_%7Bn%2B1%7D%29%5Cright%5D)

Lab 문제는 `y′=2y`, `y(0)=1`, `T=1`, `h=1/n`이다. 먼저 예측값을 구하면 `ỹₙ₊₁=(1+2h)yₙ`이고, 교정 단계에 대입하면 다음과 같다.

![Heun update factor](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7Dy_%7Bn%2B1%7D%3D%281%2B2h%2B2h%5E2%29y_n)

![Heun limit](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7Dy_n%3D%5Cleft%281%2B%5Cfrac2n%2B%5Cfrac2%7Bn%5E2%7D%5Cright%29%5En%5Clongrightarrow%20e%5E2)

로그를 취하면 `n ln(1+2/n+2/n²) → 2`이므로 정확해 `y(1)=e²`로 수렴한다. Forward Euler보다 한 단계 안에서 두 기울기를 사용하므로 이 문제의 전역 오차는 2차이다.

### 2.3 Modified Euler (Midpoint RK2)

다른 RK2 선택 `w₁=0`, `w₂=1`, `α=β=1/2`는 구간 중간점의 기울기를 사용하는 방법이다.

![Midpoint RK2](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7Dy_%7Bn%2B1%7D%3Dy_n%2Bh%20f%5Cleft%28t_n%2B%5Cfrac%20h2%2C%5C%2Cy_n%2B%5Cfrac%20h2%20f%28t_n%2Cy_n%29%5Cright%29)

Heun과 Midpoint RK2는 세부 업데이트는 다르지만 둘 다 일반적인 조건 아래 **2차 방법**이다.

## 3. Fourth-Order Runge-Kutta (RK4)

RK4는 구간 시작, 중간 두 번, 끝에서 계산한 기울기 네 개를 조합한다.

![RK4 first stages](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7DK_1%3Df%28t_n%2Cy_n%29%2C%5Cquad%20K_2%3Df%28t_n%2B%5Cfrac%20h2%2Cy_n%2B%5Cfrac%20h2K_1%29)

![RK4 last stages](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7DK_3%3Df%28t_n%2B%5Cfrac%20h2%2Cy_n%2B%5Cfrac%20h2K_2%29%2C%5Cquad%20K_4%3Df%28t_n%2Bh%2Cy_n%2BhK_3%29)

![RK4 update](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7Dy_%7Bn%2B1%7D%3Dy_n%2B%5Cfrac%20h6%28K_1%2B2K_2%2B2K_3%2BK_4%29)

RK4의 한 단계 결함(적절한 의미의 local error)은 `O(h⁵)`이고, 유한 구간의 전역 오차는 `O(h⁴)`가 된다.

## 4. Accuracy Comparison — Python Visualization

강의에서는 Euler, Heun, RK4의 수렴 차수를 실제 오차로 비교한다. 아래 그래프는 **별도의 시각화 예제**로 `y′=2y`, `y(0)=1`, `T=1`을 사용해 각 방법의 오차를 비교한 것이다.

![Accuracy comparison for Euler Heun RK4](/assets/images/ode_schemes_accuracy_comparison.png)

*Figure 1. 추가 Python 시각화. 시간 간격을 줄일 때 Euler는 대략 1차, Heun은 2차, RK4는 4차 속도로 오차가 감소한다.*

**Python Implementation**

```python
import numpy as np
import matplotlib.pyplot as plt

hs = np.array([1/10, 1/20, 1/40, 1/80, 1/160])
exact = np.exp(2)

def stability_factor(z, method):
    if method == "Euler": return 1 + z
    if method == "Heun": return 1 + z + z*z/2
    if method == "RK4": return 1 + z + z*z/2 + z**3/6 + z**4/24

for method in ["Euler", "Heun", "RK4"]:
    errors = [abs(stability_factor(2*h, method)**round(1/h) - exact) for h in hs]
    plt.loglog(hs, errors, 'o-', label=method)

plt.xlabel("Step size h")
plt.ylabel("Absolute error at t=1")
plt.legend()
plt.grid(True, which="both")
plt.show()
```

## 5. Adaptive Step-Size Methods

시간 간격을 무조건 작게 하면 오차는 줄어들 수 있지만 계산량이 늘어난다. **Adaptive Method**는 오차 추정치를 바탕으로 필요한 구간에서는 작은 시간 간격을 쓰고, 변화가 완만한 구간에서는 큰 시간 간격을 사용한다. 강의에서는 Runge–Kutta–Fehlberg (RKF45)를 대표적인 예로 소개한다.


## Computational Experiments — Lab 3

이 절은 `Lab-3.pdf`의 세 실습과 `LAB_3.ipynb`의 실제 구현·저장 출력을 반영한다. 하나의 스칼라 미분방정식에서 벗어나 벡터 상태를 구성하는 방법, 고차 Runge–Kutta의 차수를 실험으로 검증하는 방법, 결합 진동계의 변위를 시각화하는 방법을 연결한다.

### Lab 3 — Exercise 1: Forced Spring-Mass System with RK4

*Source: Lab-3.pdf, Exercise 1; Lec-2_3_scheme.pdf, RK4 and Systems*

강제진동계 `y″ + (κ/m)y = (F₀/m)cos(μt)`를 위치와 속도의 **2차원 1계 시스템**으로 변환한다. `m = 1`, `F₀ = 40`, `κ = 16π²`, `μ = 2π`, `T = 1`을 사용하고, 정확해의 위상 궤적과 RK4의 계산 궤적을 비교한다. 초기 변위는 강의자료의 `A = 1` 선택에서 계산한다.

**Python Implementation**

```python
import numpy as np
import matplotlib.pyplot as plt

# Given parameters
m = 1.0
F0 = 40.0
k = 16 * np.pi**2
mu = 2 * np.pi

omega = np.sqrt(k / m)

# The instructor's choice
A = 1.0
B = F0 / (m * (omega**2 - mu**2))
c0 = A + B

print(f"omega = {omega:.6f}")
print(f"c0 = {c0:.6f}")


# State vector: Y = [y1, y2],
# where y1 = y and y2 = -y'/omega
def forced_spring_rhs(t, Y):
    y1, y2 = Y

    dy1 = -omega * y2
    dy2 = omega * y1 - (F0 / (m * omega)) * np.cos(mu * t)

    return np.array([dy1, dy2], dtype=float)


# General vector-form RK4 solver
def rk4(f, t, Y0):
    Y = np.zeros((len(t), len(Y0)), dtype=float)
    Y[0] = Y0

    for n in range(len(t) - 1):
        h = t[n + 1] - t[n]

        k1 = f(t[n], Y[n])
        k2 = f(t[n] + h / 2, Y[n] + h * k1 / 2)
        k3 = f(t[n] + h / 2, Y[n] + h * k2 / 2)
        k4 = f(t[n] + h, Y[n] + h * k3)

        Y[n + 1] = Y[n] + (h / 6) * (
            k1 + 2*k2 + 2*k3 + k4
        )

    return Y


# Exact solution
def exact_solution(t):
    y1 = A * np.cos(omega * t) + B * np.cos(mu * t)
    y2 = A * np.sin(omega * t) + B * (mu / omega) * np.sin(mu * t)
    return y1, y2


# Solve on [0, 1]
T = 1.0
N = 1000
t = np.linspace(0.0, T, N + 1)

Y0 = np.array([c0, 0.0])
Y_num = rk4(forced_spring_rhs, t, Y0)

y1_exact, y2_exact = exact_solution(t)

# Phase-space trajectory: y1 versus y2
plt.figure(figsize=(8, 5))
plt.plot(Y_num[:, 0], Y_num[:, 1], label="RK4 numerical solution")
plt.plot(y1_exact, y2_exact, "--", label="Exact solution")
plt.xlabel(r"$y_1(t)=y(t)$")
plt.ylabel(r"$y_2(t)=-y'(t)/\omega$")
plt.title("Phase-space trajectory of the forced spring-mass system")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
```

**실행 결과**

```text
omega = 12.566371
c0 = 1.337737
```

![Lab 3 — Exercise 1: Forced Spring-Mass System with RK4 — Figure 1](/assets/images/ode_lab3_ex1_1.png)

*Figure. Lab 3 — Exercise 1: Forced Spring-Mass System with RK4 (1/1).*

### Lab 3 — Exercise 2: Heun vs RK4 Convergence

*Source: Lab-3.pdf, Exercise 2; Lec-2_3_scheme.pdf, Accuracy Comparison*

비선형 진자 `θ″ + sin(θ) = 0`를 위치·각속도의 1계 벡터 시스템으로 바꾼다. `solve_ivp`를 이용한 고정밀 기준해와 비교하여 `N = 100, 200, 400, 800`에서 **유클리드(ℓ₂) 전역오차**를 계산한다. Heun은 차수 약 2, RK4는 차수 약 4를 보여준다. RK4 오차가 극도로 작아지는 구간에서는 부동소수점 오차 때문에 추정 차수가 약간 흔들릴 수 있다.

**Python Implementation**

```python
import numpy as np
from scipy.integrate import solve_ivp

# x = [theta, omega], where omega = theta'
def f(t, x):
    theta, omega = x
    return np.array([omega, -np.sin(theta)])


def heun(N):
    T = 1.0
    h = T / N
    x = np.array([1.0, 0.0])  # theta(0)=1, omega(0)=0

    for n in range(N):
        t = n * h

        k1 = f(t, x)
        k2 = f(t + h, x + h * k1)

        x = x + (h / 2) * (k1 + k2)

    return x


def rk4(N):
    T = 1.0
    h = T / N
    x = np.array([1.0, 0.0])

    for n in range(N):
        t = n * h

        k1 = f(t, x)
        k2 = f(t + h / 2, x + (h / 2) * k1)
        k3 = f(t + h / 2, x + (h / 2) * k2)
        k4 = f(t + h, x + h * k3)

        x = x + (h / 6) * (k1 + 2*k2 + 2*k3 + k4)

    return x


# High-precision reference solution at T=1
reference = solve_ivp(
    f,
    (0.0, 1.0),
    [1.0, 0.0],
    method="DOP853",
    rtol=3e-14,
    atol=1e-15
).y[:, -1]


def print_convergence_table(method, method_name):
    Ns = [100, 200, 400, 800]
    errors = []

    for N in Ns:
        numerical = method(N)
        error = np.linalg.norm(numerical - reference, ord=2)
        errors.append(error)

    print(f"\n{method_name}")
    print("N       h           E_N                 alpha")

    for i, N in enumerate(Ns):
        h = 1.0 / N

        if i == 0:
            alpha_text = "---"
        else:
            alpha = np.log(errors[i - 1] / errors[i]) / np.log(2)
            alpha_text = f"{alpha:.6f}"

        print(f"{N:<7} {h:<11.6g} {errors[i]:.12e}   {alpha_text}")


print("Reference solution at T=1:")
print(f"theta(1) = {reference[0]:.15f}")
print(f"omega(1) = {reference[1]:.15f}")

print_convergence_table(heun, "Heun method")
print_convergence_table(rk4, "RK4 method")
```

**실행 결과**

```text
Reference solution at T=1:
theta(1) = 0.600085366127506
omega(1) = -0.754963713953130

Heun method
N       h           E_N                 alpha
100     0.01        9.215314300750e-06   ---
200     0.005       2.303371625687e-06   2.000286
400     0.0025      5.757846819650e-07   2.000146
800     0.00125     1.439388235167e-07   2.000074

RK4 method
N       h           E_N                 alpha
100     0.01        4.633082877220e-11   ---
200     0.005       2.890393539983e-12   4.002635
400     0.0025      1.804692712125e-13   4.001441
800     0.00125     1.100857650037e-14   4.035053
```

### Lab 3 — Exercise 3: Coupled Spring-Mass Systems

*Source: Lab-3.pdf, Exercise 3; Lec-2_3_scheme.pdf, High-Order Equations & Systems*

두 질량의 위치와 속도를 묶어 **4차원 상태 벡터**를 구성한다. 각각의 2계 방정식을 1계 미분방정식 네 개로 바꾸고, 일반 벡터 RK4 함수로 적분한 뒤 두 변위의 시간 변화를 같은 그래프에 나타낸다. 위상차와 서로 영향을 주고받는 진동 패턴을 확인할 수 있다.

**Python Implementation**

```python
import numpy as np
import matplotlib.pyplot as plt


# Convert two general second-order ODEs into a first-order system.
# State vector: z = [u, u_prime, v, v_prime]
def make_first_order_system(f1, f2):
    def rhs(t, z):
        u, u_prime, v, v_prime = z

        return np.array([
            u_prime,
            f1(t, u, u_prime, v, v_prime),
            v_prime,
            f2(t, u, u_prime, v, v_prime)
        ], dtype=float)

    return rhs


# General vector-form RK4 solver
def rk4(f, t, z0):
    z = np.zeros((len(t), len(z0)), dtype=float)
    z[0] = z0

    for n in range(len(t) - 1):
        h = t[n + 1] - t[n]

        k1 = f(t[n], z[n])
        k2 = f(t[n] + h / 2, z[n] + h * k1 / 2)
        k3 = f(t[n] + h / 2, z[n] + h * k2 / 2)
        k4 = f(t[n] + h, z[n] + h * k3)

        z[n + 1] = z[n] + (h / 6) * (
            k1 + 2*k2 + 2*k3 + k4
        )

    return z


# Model parameters
m1 = 1.0
m2 = 1.0
k1 = 1.0
k2 = 1.0
k3 = 1.0


# Accelerations of the two masses
def f1(t, u, u_prime, v, v_prime):
    return (-(k1 + k2) * u + k2 * v) / m1


def f2(t, u, u_prime, v, v_prime):
    return (k2 * u - (k2 + k3) * v) / m2


# Construct the first-order system
rhs = make_first_order_system(f1, f2)


# Initial conditions: [u(0), u'(0), v(0), v'(0)]
z0 = np.array([1.0, 0.0, 0.0, 0.0])

# Time interval and grid
T = 20.0
N = 4000
t = np.linspace(0.0, T, N + 1)

# Solve with RK4
solution = rk4(rhs, t, z0)

u = solution[:, 0]
v = solution[:, 2]


# Plot both displacement curves
plt.figure(figsize=(9, 5))
plt.plot(t, u, label="Mass 1: u(t)")
plt.plot(t, v, label="Mass 2: v(t)")
plt.xlabel("Time")
plt.ylabel("Displacement")
plt.title("Coupled Oscillations of Two Masses")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
```

![Lab 3 — Exercise 3: Coupled Spring-Mass Systems — Figure 1](/assets/images/ode_lab3_ex3_1.png)

*Figure. Lab 3 — Exercise 3: Coupled Spring-Mass Systems (1/1).*

## Summary and Takeaways

Taylor 전개는 방법의 정확도 조건을 제공하고, RK 방법은 다중 기울기 평가로 그 정확도를 효율적으로 구현한다. **Euler 1차, Heun 2차, RK4 4차**라는 수렴차수는 그래프의 기울기와 연결된다.

## References

- Raj Dhara, *Lec-2_3_scheme*, 2026-10-12, pp. 5–9.
- Raj Dhara, *Numerical Method Lab: ODE Numerical Schemes & Convergence*, Exercise 4.
