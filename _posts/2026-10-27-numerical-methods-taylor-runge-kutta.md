---
title: "Taylor and Runge-Kutta Methods"
date: 2026-10-27 18:20:00 +0100
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

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7Dy%28t_n%2Bh%29%3Dy%28t_n%29%2Bhy%5Cprime%28t_n%29%2B%5Cfrac%7Bh%5E2%7D%7B2%7Dy%5Cprime%5Cprime%28t_n%29%2B%5Ccdots%2B%5Cfrac%7Bh%5Em%7D%7Bm%21%7Dy%5E%7B%28m%29%7D%28t_n%29%2BO%28h%5E%7Bm%2B1%7D%29" alt="Taylor expansion" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

미분방정식 `y′=f(t,y)`를 대입하면 첫 도함수는 이미 알고 있다. 2차 미분은 합성함수 미분을 통해 구한다.

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7Dy%5Cprime%5Cprime%3Df_t%2Bf_yf" alt="Second derivative from ODE" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

**m차 Taylor 방법**은 이 전개에서 m차까지 남기고 나머지항을 버리는 방식이다. Forward Euler는 m=1인 특수한 사례이다. 다만 고차 도함수를 반복 계산해야 하므로 실제로는 구현 부담이 크다.

### 1.1 Lecture Example — Third-Order Taylor Method

강의의 `y′ = y − t³ + t + 1`, `y(0)=0.5`에서 `f=y−t³+t+1`이다.

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7Df%5Cprime%3Df_t%2Bf_y%20f%3Dy-t%5E3-3t%5E2%2Bt%2B2" alt="First total derivative" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7Df%5Cprime%5Cprime%3Dy-t%5E3-3t%5E2-5t%2B2" alt="Second total derivative" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

따라서 `T₃ = f + (h/2)f′ + (h²/6)f″`이며 `yₙ₊₁ = yₙ + h T₃(tₙ,yₙ)`으로 업데이트한다.

## 2. Runge-Kutta Methods

Runge–Kutta 방법은 Taylor와 유사한 정확도를 얻되 고차 도함수 대신 **여러 지점에서 계산한 기울기의 가중평균**을 사용한다.

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7Dy_%7Bn%2B1%7D%3Dy_n%2Bh%5Csum_%7Bi%3D1%7D%5E%7Bm%7Dw_iK_i" alt="Runge Kutta weighted slopes" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

기울기를 어느 지점에서 평가하고 어떤 비율로 섞을지 선택하여 Taylor 전개 계수와 일치시키는 것이 핵심이다.

### 2.1 Second-Order Runge-Kutta (RK2)

일반적인 2단계 RK 방법에서 `K₁=f(tₙ,yₙ)`, `K₂=f(tₙ+αh,yₙ+βhK₁)`라 두면 2차 정확도를 위해 `w₁+w₂=1`, `w₂α=1/2`, `w₂β=1/2`가 필요하다.

### 2.2 Heun's Method — Lab Handout Exercise 4

`w₁=w₂=1/2`, `α=β=1`로 두면 Heun의 Predictor–Corrector 공식이 된다.

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cwidetilde%20y_%7Bn%2B1%7D%3Dy_n%2Bhf%28t_n%2Cy_n%29" alt="Heun predictor" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7Dy_%7Bn%2B1%7D%3Dy_n%2B%5Cfrac%20h2%5Cleft%5Bf%28t_n%2Cy_n%29%2Bf%28t_%7Bn%2B1%7D%2C%5Cwidetilde%20y_%7Bn%2B1%7D%29%5Cright%5D" alt="Heun corrector" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

Lab 문제는 `y′=2y`, `y(0)=1`, `T=1`, `h=1/n`이다. 먼저 예측값을 구하면 `ỹₙ₊₁=(1+2h)yₙ`이고, 교정 단계에 대입하면 다음과 같다.

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7Dy_%7Bn%2B1%7D%3D%281%2B2h%2B2h%5E2%29y_n" alt="Heun update factor" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7Dy_n%3D%5Cleft%281%2B%5Cfrac2n%2B%5Cfrac2%7Bn%5E2%7D%5Cright%29%5En%5Clongrightarrow%20e%5E2" alt="Heun limit" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

로그를 취하면 `n ln(1+2/n+2/n²) → 2`이므로 정확해 `y(1)=e²`로 수렴한다. Forward Euler보다 한 단계 안에서 두 기울기를 사용하므로 이 문제의 전역 오차는 2차이다.

### 2.3 Modified Euler (Midpoint RK2)

다른 RK2 선택 `w₁=0`, `w₂=1`, `α=β=1/2`는 구간 중간점의 기울기를 사용하는 방법이다.

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7Dy_%7Bn%2B1%7D%3Dy_n%2Bh%20f%5Cleft%28t_n%2B%5Cfrac%20h2%2C%5C%2Cy_n%2B%5Cfrac%20h2%20f%28t_n%2Cy_n%29%5Cright%29" alt="Midpoint RK2" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

Heun과 Midpoint RK2는 세부 업데이트는 다르지만 둘 다 일반적인 조건 아래 **2차 방법**이다.

## 3. Fourth-Order Runge-Kutta (RK4)

RK4는 구간 시작, 중간 두 번, 끝에서 계산한 기울기 네 개를 조합한다.

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7DK_1%3Df%28t_n%2Cy_n%29%2C%5Cquad%20K_2%3Df%28t_n%2B%5Cfrac%20h2%2Cy_n%2B%5Cfrac%20h2K_1%29" alt="RK4 first stages" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7DK_3%3Df%28t_n%2B%5Cfrac%20h2%2Cy_n%2B%5Cfrac%20h2K_2%29%2C%5Cquad%20K_4%3Df%28t_n%2Bh%2Cy_n%2BhK_3%29" alt="RK4 last stages" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7Dy_%7Bn%2B1%7D%3Dy_n%2B%5Cfrac%20h6%28K_1%2B2K_2%2B2K_3%2BK_4%29" alt="RK4 update" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

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



## Detailed Lecture Notes — Taylor Series, RK Derivations, and ODE Systems

강의 `Lec-2_3_scheme.pdf`의 **Higher-Order Taylor Methods, Runge–Kutta Methods, Accuracy Comparison, Multi-Step Preview, High-Order Equations, Exercises, Big O and Small o Appendix**를 보다 정확한 계산 순서로 정리한다. Lab 3의 Python 코드와 시각화는 위에 이미 배치되어 있으며 원래 경로를 유지한다.

### 1.2 Full Taylor Expansion and the m-th Order Method


해가 충분히 매끄러울 때

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20t_n" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>


에서

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20t_%7Bn%2B1%7D%3Dt_n%2Bh" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

로 전개하면

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y%28t_n%2Bh%29%3Dy%28t_n%29%2Bhy%27%28t_n%29%2B%5Cfrac%7Bh%5E2%7D%7B2%21%7Dy%27%27%28t_n%29%2B%5Ccdots%2B%5Cfrac%7Bh%5Em%7D%7Bm%21%7Dy%5E%7B%28m%29%7D%28t_n%29%2B%5Cfrac%7Bh%5E%7Bm%2B1%7D%7D%7B%28m%2B1%29%21%7Dy%5E%7B%28m%2B1%29%7D%28%5Cxi_n%29" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>


이다.

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y%27%3Df%28t%2Cy%29" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

를 이용하면 시간에 대한 **전체 미분(total derivative)**은

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cfrac%7Bd%7D%7Bdt%7Df%28t%2Cy%28t%29%29%3Df_t%28t%2Cy%29%2Bf_y%28t%2Cy%29f%28t%2Cy%29" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

이므로

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y%27%27%28t%29%3Df_t%2Bf_y%20f" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

이다. 2차 Taylor 방법은

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y_%7Bn%2B1%7D%3Dy_n%2Bh%20f%28t_n%2Cy_n%29%2B%5Cfrac%7Bh%5E2%7D%7B2%7D%28f_t%2Bf_y%20f%29%28t_n%2Cy_n%29" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>


이다. 일반적으로

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20T_m%3Df%2B%5Cfrac%20h%7B2%21%7Df%5E%7B%281%29%7D%2B%5Ccdots%2B%5Cfrac%7Bh%5E%7Bm-1%7D%7D%7Bm%21%7Df%5E%7B%28m-1%29%7D" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

로 두면

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y_%7Bn%2B1%7D%3Dy_n%2BhT_m%28t_n%2Cy_n%29" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>


이며

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20m%3D1" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>


이면 Euler로 돌아간다. Taylor법은 정확한 차수를 달성할 수 있지만 **

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

의 반복 미분 계산이 복잡**하다는 약점이 있다.

### 1.3 Detailed Lecture Example — Compute T3

강의 예제는

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y%27%3Df%28t%2Cy%29%3Dy-t%5E3%2Bt%2B1%2C%5Cqquad%20y%280%29%3D0.5" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>


이다. 먼저

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y%27%3Df" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

를 그대로 적는다.

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f%28t%2Cy%28t%29%29%3Dy-t%5E3%2Bt%2B1" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

한 번 시간미분하면

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f%27%3Dy%27-3t%5E2%2B1" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%3Df-3t%5E2%2B1%3Dy-t%5E3-3t%5E2%2Bt%2B2" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

이다. 다시 한 번 미분하면

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f%27%27%3Dy%27-3t%5E2-6t%2B1" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%3Dy-t%5E3-3t%5E2-5t%2B2" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

이다. 따라서

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20T_3%28t%2Cy%29%3Df%2B%5Cfrac%20h2%20f%27%2B%5Cfrac%7Bh%5E2%7D%7B6%7Df%27%27" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%3Dy-t%5E3%2Bt%2B1%2B%5Cfrac%20h2%28y-t%5E3-3t%5E2%2Bt%2B2%29%2B%5Cfrac%7Bh%5E2%7D%7B6%7D%28y-t%5E3-3t%5E2-5t%2B2%29" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>


를 얻는다.

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20h%5E2%2F6" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>


을 곱하는 마지막 항은

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y%5E%7B%283%29%7D%3Df%27%27" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>


에 대응한다. 강의 Exercise에서는 같은 예제의

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20T_4" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>


및

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20h%3D1%2F2" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

에서 2회의 Taylor 업데이트를 요구한다.

#### Additional Calculation — T4 and Two Taylor Steps

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f%5E%7B%283%29%7D%3Dy%27-3t%5E2-6t-5%3Dy-t%5E3-3t%5E2-5t-4" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

이므로

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20T_4%3DT_3%2B%5Cfrac%7Bh%5E3%7D%7B24%7D%28y-t%5E3-3t%5E2-5t-4%29" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

이다. 아래 코드는 **강의 Exercise 1의 계산을 재현하는 추가 코드**다. 새로운 그림을 요구하지 않으므로 기존 PNG는 건드리지 않는다.

```python
import numpy as np

def taylor_step(t, y, h, order):
    f0 = y - t**3 + t + 1
    f1 = y - t**3 - 3*t**2 + t + 2
    f2 = y - t**3 - 3*t**2 - 5*t + 2
    f3 = y - t**3 - 3*t**2 - 5*t - 4
    update = y + h*f0 + h**2*f1/2 + h**3*f2/6
    if order == 4:
        update += h**4*f3/24
    return update

for order in (3, 4):
    t, y = 0., 0.5
    for _ in range(2):
        y = taylor_step(t, y, 0.5, order)
        t += 0.5
    # Lecture exact solution: 4 + 5*t + 3*t**2 + t**3 - 3.5*np.exp(t)
    exact = 4 + 5*t + 3*t**2 + t**3 - 3.5*np.exp(t)
    print(order, y, abs(y-exact))
```


**Formula check:** 초기값

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y%280%29%3D0.5" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>


는 강의에 적힌 정확해

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y%3D4%2B5t%2B3t%5E2%2Bt%5E3-%287%2F2%29e%5Et" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>


와 일치한다. 이 표현을 미분해 원래

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

와 비교하여 사용할 수 있다.

### 2.4 Why RK Methods Avoid High Derivatives


일반적인 Runge–Kutta 방법은 구간 내 여러 지점에서 기울기

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20K_i" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

를 평가한 뒤 가중평균으로 한 단계 전진한다.

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y_%7Bn%2B1%7D%3Dy_n%2Bh%5Csum_%7Bi%3D1%7D%5Em%20w_iK_i%2C%5Cqquad%20%5Csum_iw_i%3D1" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>


Taylor 방법처럼

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f_t%2Cf_y%2Cf_%7Btt%7D" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>


를 직접 구하는 대신 **어디에서

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

를 평가할지와 어떤 가중치를 쓸지**를 정하여 차수를 확보한다.

### 2.5 Complete RK2 Coefficient Derivation

다음 2단계 방법을 생각한다.

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20K_1%3Df%28t_n%2Cy_n%29" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20K_2%3Df%28t_n%2B%5Calpha%20h%2Cy_n%2B%5Cbeta%20hK_1%29" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y_%7Bn%2B1%7D%3Dy_n%2Bh%28w_1K_1%2Bw_2K_2%29" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

한편 정확해의 Taylor 전개는

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y%28t%2Bh%29%3Dy%2Bhy%27%2B%5Cfrac%7Bh%5E2%7D%7B2%7Dy%27%27%2BO%28h%5E3%29" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%3Dy%2Bhf%2B%5Cfrac%7Bh%5E2%7D%7B2%7D%28f_t%2Bf_yf%29%2BO%28h%5E3%29" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>


이다.

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20K_2" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

를 두 변수에 대해 Taylor 전개하면

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20K_2%3Df%2B%5Calpha%20h%20f_t%2B%5Cbeta%20h%20f_yK_1%2BO%28h%5E2%29" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%3Df%2B%5Calpha%20h%20f_t%2B%5Cbeta%20h%20f_yf%2BO%28h%5E2%29" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

가 된다. 이를 RK식에 대입하면

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y_%7Bn%2B1%7D%3Dy%2Bh%28w_1%2Bw_2%29f%2Bh%5E2%28w_2%5Calpha%20f_t%2Bw_2%5Cbeta%20f_yf%29%2BO%28h%5E3%29" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>


이다. 정확해 전개에서

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>


,

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f_t" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>


,

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f_yf" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

의 계수를 비교해

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20w_1%2Bw_2%3D1%2C%5Cqquad%20w_2%5Calpha%3D%5Cfrac12%2C%5Cqquad%20w_2%5Cbeta%3D%5Cfrac12" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>


라는 **강의 식 (32)의 조건**을 얻는다. 이 세 관계를 만족하면 RK2의 국소 한 단계 결손량은

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20O%28h%5E3%29" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

이며 적절한 안정성 조건 아래 전역 2차 수렴을 기대한다.

### 2.6 Heun and Modified Euler as Two RK2 Choices


**Heun:**

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20w_1%3Dw_2%3D1%2F2" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>


,

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Calpha%3D%5Cbeta%3D1" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

로 설정하면

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20K_1%3Df%28t_n%2Cy_n%29%2C%5Cqquad%20K_2%3Df%28t_n%2Bh%2Cy_n%2BhK_1%29" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y_%7Bn%2B1%7D%3Dy_n%2B%5Cfrac%20h2%28K_1%2BK_2%29" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>


이다.

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y_n%2BhK_1" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

로 먼저 예측하고, 시작점과 끝점 기울기를 평균하는 predictor-corrector 구조다.


**Modified Euler / Midpoint:**

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20w_1%3D0%2Cw_2%3D1" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>


,

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Calpha%3D%5Cbeta%3D1%2F2" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

로 설정하면

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y_%7Bn%2B1%7D%3Dy_n%2Bh%20f%5Cleft%28t_n%2B%5Cfrac%20h2%2Cy_n%2B%5Cfrac%20h2f%28t_n%2Cy_n%29%5Cright%29" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>


가 된다. 둘 모두 전역 2차이지만 중간에

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

를 평가하는 방식이 다르다.

#### Lab Handout Exercise 4 — Full Limit Calculation

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y%27%3D2y%2Cy%280%29%3D1%2Ch%3D1%2Fn" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>


에 Heun을 적용한다. 첫 기울기는

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20K_1%3D2y_n" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>


이고 예측값은

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cwidetilde%20y_%7Bn%2B1%7D%3Dy_n%2B2hy_n%3D%281%2B2h%29y_n" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>


이다. 두 번째 기울기는

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20K_2%3D2%5Cwidetilde%20y_%7Bn%2B1%7D%3D2%281%2B2h%29y_n" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

이므로

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y_%7Bn%2B1%7D%3Dy_n%2B%5Cfrac%20h2%282y_n%2B2%281%2B2h%29y_n%29" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%3D%281%2B2h%2B2h%5E2%29y_n" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

이다. 따라서

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y_n%3D%281%2B2%2Fn%2B2%2Fn%5E2%29%5En" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

이다. 로그를 취하면

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Clog%20y_n%3Dn%5Clog%281%2B2%2Fn%2B2%2Fn%5E2%29" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%3Dn%5Cleft%28%5Cfrac2n%2B%5Cfrac2%7Bn%5E2%7D-%5Cfrac12%28%5Cfrac2n%2B%5Cfrac2%7Bn%5E2%7D%29%5E2%2BO%28n%5E%7B-3%7D%29%5Cright%29" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%3D2%2BO%28n%5E%7B-2%7D%29" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>


이므로

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Clim_%7Bn%5Cto%5Cinfty%7Dy_n%3De%5E2" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>


이다. 정확해

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y%28t%29%3De%5E%7B2t%7D" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

와 일치한다.

### 3.1 Classical RK4: Four Stages and Their Roles

강의에 나오는 대표적인 4차 방법은

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20K_1%3Df%28t_n%2Cy_n%29" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20K_2%3Df%28t_n%2Bh%2F2%2Cy_n%2BhK_1%2F2%29" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20K_3%3Df%28t_n%2Bh%2F2%2Cy_n%2BhK_2%2F2%29" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20K_4%3Df%28t_n%2Bh%2Cy_n%2BhK_3%29" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y_%7Bn%2B1%7D%3Dy_n%2B%5Cfrac%20h6%28K_1%2B2K_2%2B2K_3%2BK_4%29" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>


이다. 시작 기울기와 끝 기울기에 가중치 1, 가운데 두 기울기에 가중치 2를 준다. 충분히 매끄러운 문제에서 **한 단계 결손량

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20O%28h%5E5%29" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>


**, **고정 구간 전역오차

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20O%28h%5E4%29" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>


**라는 차수를 갖는다. 강의에는 잔여항을

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20h%5E5y%5E%7B%285%29%7D%28%5Cxi%29%2F5%21" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

형태로 표기하는 설명이 있으며, 실제 RK4의 정확한 잔여항 계수는 문제에 따라 달라질 수 있으므로 그 표현은 **크기 차수 설명**으로 이해해야 한다.

### 3.2 From Local Error to Global Error


각 단계의 오차가

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20O%28h%5E5%29" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>


이고

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20n%3D%28T-t_0%29%2Fh" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

번의 단계가 있다면 단순한 차수 계산에서는

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20n%5Ccdot%20O%28h%5E5%29%3D%5Cfrac%7BT-t_0%7D%7Bh%7DO%28h%5E5%29%3DO%28h%5E4%29" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

이다. 다만 엄밀한 전역 상계에는 오차 전달의 안정성 또는 Lipschitz 상계가 필요하다. 예를 들어

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%7Ce_%7Bn%2B1%7D%7C%5Cle%281%2BCh%29%7Ce_n%7C%2BC_1h%5E5" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>


를 반복하면 등비합이

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20O%281%2Fh%29" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>


이므로

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%7Ce_n%7C%3DO%28h%5E4%29" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

를 얻는다. 강의 Exercise 2가 묻는 local-to-global 논리의 핵심이다.

### 4.1 Accuracy Comparison from the Lecture Table


강의는 강제진동계에 대해

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20t%3D1" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

에서의 오차를 다음과 같이 비교한다. 괄호는 두 격자 크기를 비교한 경험적 수렴차수다.

| 1/h | Euler | Heun | RK4 |
|---:|---:|---:|---:|
| 100 | 1.19 | 3.31E-2 | 2.61E-5 |
| 200 | 4.83E-1 (1.3) | 8.27E-3 (2.0) | 1.63E-6 (4.0) |
| 400 | 2.18E-1 (1.1) | 2.07E-3 (2.0) | 1.02E-7 (4.0) |
| 800 | 1.04E-1 (1.1) | 5.17E-4 (2.0) | 6.38E-9 (4.0) |


벡터 오차는

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20E%28h%29%3D%5Csqrt%7B%28y_%7B1%2Ch%7D%281%29-y_1%281%29%29%5E2%2B%28y_%7B2%2Ch%7D%281%29-y_2%281%29%29%5E2%7D" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>


와 같이

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cell%5E2" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>


노름으로 측정한다. 예를 들어 Heun의

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20h%3D1%2F200%2C1%2F400" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

을 비교하면

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20p%5Capprox%5Cfrac%7B%5Clog%28%288.27%5Ctimes10%5E%7B-3%7D%29%2F%282.07%5Ctimes10%5E%7B-3%7D%29%29%7D%7B%5Clog2%7D%5Capprox2" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>


이다. RK4에서도 같은 방식으로

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20p%5Capprox4" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

가 나온다. 이 표는 단지 오차가 작다는 사실뿐 아니라, **격자를 절반으로 줄였을 때 감소 비율이 이론적 차수를 따른다**는 사실을 보여준다.

### 4.2 Forced Spring-Mass System: Exact Solution and First-Order Form

강의의 강제진동 방정식은

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y%27%27%28t%29%2B%5Cfrac%5Ckappa%20m%20y%28t%29%3D%5Cfrac%7BF_0%7Dm%5Ccos%28%5Cmu%20t%29%2C%5Cqquad%20y%280%29%3Dc_0%2C%5C%20y%27%280%29%3D0" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>


이다.

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Comega%3D%5Csqrt%7B%5Ckappa%2Fm%7D" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>


를 정의하고

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cmu%5Cne%5Comega" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

라고 하면

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y%28t%29%3DA%5Ccos%28%5Comega%20t%29%2B%5Cfrac%7BF_0%7D%7Bm%28%5Comega%5E2-%5Cmu%5E2%29%7D%5Ccos%28%5Cmu%20t%29" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>


가 해다. 초기조건으로

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20c_0%3DA%2BF_0%2F%5Bm%28%5Comega%5E2-%5Cmu%5E2%29" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>


를 정한다. 강의는

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20m%3D1%2CF_0%3D40%2CA%3D1%2C%5Comega%3D4%5Cpi%2C%5Cmu%3D2%5Cpi" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

를 선택한다. 따라서 두 진동수의 비가 유리수이고 공통 기본주기가 1이라는 설명과 연결된다.

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y_1%3Dy%2C%5C%20y_2%3D-y_1%27%2F%5Comega" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

로 놓으면

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y_1%27%3D-%5Comega%20y_2" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y_2%27%3D%5Comega%20y_1-%5Cfrac%7BF_0%7D%7Bm%5Comega%7D%5Ccos%28%5Cmu%20t%29" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>


이고 초기조건은

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y_1%280%29%3Dc_0%2C%5C%20y_2%280%29%3D0" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

이다. **Lab 3 Exercise 1의 RK4 벡터 구현**은 바로 이 연립계에 대한 수치해다.

### 5.1 Adaptive Methods and RKF45


시간 간격

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20h" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>


를 작게 만들면 일반적으로 정확도는 증가하지만 함수 평가 횟수도 증가한다. **Adaptive method**는 해가 급격히 변화하는 구간에서는

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20h" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>


를 줄이고, 완만한 구간에서는

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20h" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

를 키워 전체 계산 비용을 줄인다. 강의는 RK4와 RK5를 짝지어 지역 오차를 추정하는 **Runge–Kutta–Fehlberg (RKF45)**를 예로 든다. 여기서는 강의에 제시되지 않은 RKF45의 모든 계수표를 임의로 추가하지 않는다.

### 6. High-Order ODEs as First-Order Systems


2계 미분방정식

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y%27%27%3Df%28t%2Cy%2Cy%27%29%2C%5C%20y%28t_0%29%3Dy_0%2C%5C%20y%27%28t_0%29%3Du_0" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>


는

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20u%3Dy%27" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

를 두면

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y%27%3Du%2C%5Cqquad%20u%27%3Df%28t%2Cy%2Cu%29" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>


인 1계 연립계로 변환된다. 상태벡터

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20Y%3D%28y%2Cu%29%5E%5Ctop" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>


를 정의하면

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20Y%27%3DF%28t%2CY%29" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

형식이 되어 vector RK4를 바로 적용할 수 있다. Lab 3의 **비선형 진자와 결합 스프링계**도 같은 원리를 사용한다.

#### Lecture Exercise 3 — Coupled Second-Order System


강의의 연립식은

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x%27%27%2Bx%27y-2y%27%27%3Dt" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>


및

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20-2y%2By%27%27%2Bx%3De%5E%7B-t%7D" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

이다. 두 번째 식에서

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y%27%27%3De%5E%7B-t%7D%2B2y-x" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

를 얻고 첫 번째에 대입하면

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x%27%27%3Dt-x%27y%2B2%28e%5E%7B-t%7D%2B2y-x%29" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>


이다.

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20u%3Dx%27%2Cv%3Dy%27" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

로 설정하면 네 개의 1계 방정식

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x%27%3Du%2C%5Cquad%20u%27%3Dt-uy%2B2e%5E%7B-t%7D%2B4y-2x" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y%27%3Dv%2C%5Cquad%20v%27%3De%5E%7B-t%7D%2B2y-x" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>


을 얻는다. 상태

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%28x%2Cu%2Cy%2Cv%29" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

를 vector RK4의 입력으로 사용하면 된다. 이는 **Lecture 2–3의 Exercise 3**을 계산 과정까지 전개한 것이다.

### Appendix: Big O and Small o


강의 마지막 Appendix는 두 표기의 차이를 정리한다.

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f%28h%29%3DO%28g%28h%29%29" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>


는 충분히 작은

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20h" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>


에서

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%7Cf%28h%29%7C%5Cle%20C%7Cg%28h%29%7C" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>


인 상수

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20C" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>


가 존재한다는 뜻이다. 반면

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f%28h%29%3Do%28g%28h%29%29" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

는

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Clim_%7Bh%5Cto0%7D%5Cfrac%7Bf%28h%29%7D%7Bg%28h%29%7D%3D0" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>


인 더 강한 조건이다(

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20g%28h%29%5Cne0" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

인 범위에서). 예를 들어

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%205h%2B2h%5E2%3DO%28h%29" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>


이고

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20h%5E2%3Do%28h%29" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

이다. 지수함수 Taylor 전개에서는

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20e%5Eh%3D1%2Bh%2B%5Cfrac%7Bh%5E2%7D2%2B%5Cfrac%7Bh%5E3%7D6%2B%5Ccdots" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20e%5Eh%3D1%2Bh%2BO%28h%5E2%29%3D1%2Bh%2Bo%28h%29" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

라고 둘 다 쓸 수 있다. **O는 상수배로 제어되는 오차 크기, o는 비교 대상보다 엄격히 빠른 소멸**을 나타낸다.

## Source-by-Source Verification — Numerical Calculations and Exercises

### T3 and T4: Differentiating Along the ODE

강의 p.5의 문제 `y′=y−x³+x+1`에서는 `f′`가 단순한 편미분 `∂f/∂x`가 아니라 **해 곡선을 따라 취하는 전미분**이라는 점이 중요하다.

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f%3Dy-x%5E3%2Bx%2B1" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f%5Cprime%3Df_x%2Bf_y%20y%5Cprime%3D%28-3x%5E2%2B1%29%2Bf%3Dy-x%5E3-3x%5E2%2Bx%2B2" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f%5Cprime%5Cprime%3Dy%5Cprime-3x%5E2-6x%2B1%3Dy-x%5E3-3x%5E2-5x%2B2" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f%5E%7B%283%29%7D%3Dy%5Cprime-3x%5E2-6x-5%3Dy-x%5E3-3x%5E2-5x-4" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

따라서

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20T_3%3Df%2B%5Cfrac%20h2%20f%5Cprime%2B%5Cfrac%7Bh%5E2%7D%7B6%7Df%5Cprime%5Cprime%2C%5Cqquad%20T_4%3DT_3%2B%5Cfrac%7Bh%5E3%7D%7B24%7Df%5E%7B%283%29%7D" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

이다. 여기서 `f^3`는 해를 따라 세 번 전미분한 값이다. **추가 검산:** `h=1/2`, `y₀=1/2`에서 3차 Taylor의 두 스텝을 정확한 유리수 연산으로 계산한 결과는 아래와 같다.

```text
Taylor order 3: step 1 = 1.6145833333; step 2 = 3.5193142361
Taylor order 4: step 1 = 1.6054687500; step 2 = 3.4892883301
Exact y(1) = 3.4860136004
```

**독립 검산:** 각 단계에서 다항식을 정확한 유리수로 계산한 뒤 소수 열 자리로 변환했다.

### RK4 Global Error: Conditions Behind the Bound

강의 p.8의 RK4 업데이트는 고전적 4단계 기울기 평균이고, 매 스텝 **국소 업데이트 오차가 O(h⁵)**, 구간 길이 T 고정 시 **전역 오차가 O(h⁴)**라는 표준 결과로 이어진다. 단순히 `N=T/h`개를 곱하는 설명만으로는 완전한 증명은 아니며, 오차가 증폭되는 정도를 제어하는 Lipschitz 안정성 상계까지 함께 사용해야 한다. 또한 강의의 `y⁽⁵⁾(ξ)` 꼴 식은 해당 도함수의 매끈함을 전제한 오차 표현으로 이해한다.

### Lecture Accuracy Table: Explicit Order Checks

강의 p.9에 제시된 수치는 소수 세 자리로 반올림되어 있으므로 경험적 차수 역시 근삿값이다. 예를 들어

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20p_%7B%5Cmathrm%7BHeun%7D%7D%3D%5Cfrac%7B%5Clog%288.27%5Ctimes10%5E%7B-3%7D%2F2.07%5Ctimes10%5E%7B-3%7D%29%7D%7B%5Clog%202%7D%5Capprox%201.998" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20p_%7B%5Cmathrm%7BRK4%7D%7D%3D%5Cfrac%7B%5Clog%281.63%5Ctimes10%5E%7B-6%7D%2F1.02%5Ctimes10%5E%7B-7%7D%29%7D%7B%5Clog%202%7D%5Capprox%203.998" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

를 얻는다. Euler는 비교 격자 수준에서 아직 점근구간으로 완전히 들어가지 않아 추정 차수가 1과 차이날 수 있다.

### Lab Handout Heun: Closed Form

`y′=2y`, `y(0)=1`, `h=1/n`에서 예측값은 `(1+2h)yₙ`, 보정값은

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y_%7Bn%2B1%7D%3Dy_n%2B%5Cfrac%20h2%282y_n%2B2%281%2B2h%29y_n%29%3D%281%2B2h%2B2h%5E2%29y_n" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

이므로

<div style="display:block;clear:both;width:100%;max-width:100%;text-align:center;overflow-x:auto;margin:0.9rem 0;">
<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y_n%3D%281%2B2%2Fn%2B2%2Fn%5E2%29%5En" alt="Equation" style="display:inline-block;max-width:100%;height:auto;vertical-align:middle;" />
</div>

이다. `u_n=2/n+2/n²`를 두고 `log(1+u_n)=u_n−u_n²/2+O(u_n³)`를 적용하면 `n log(1+u_n)=2+O(n^(−2))`, 따라서 `y_n→e²`이다. **1/n 항이 상쇄되는 계산**이 Heun의 2차 정확도를 보여준다.


Taylor 전개는 방법의 정확도 조건을 제공하고, RK 방법은 다중 기울기 평가로 그 정확도를 효율적으로 구현한다. **Euler 1차, Heun 2차, RK4 4차**라는 수렴차수는 그래프의 기울기와 연결된다.

## References

- Raj Dhara, *Lec-2_3_scheme*, 2026-10-12, pp. 5–9.
- Raj Dhara, *Numerical Method Lab: ODE Numerical Schemes & Convergence*, Exercise 4.


---
