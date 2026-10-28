---
title: "Multistep Methods and Stability"
date: 2026-10-27 18:30:00 +0100
categories:
  - Numerical Analysis
tags:
  - Multistep-Methods
  - Zero-Stability
  - Absolute-Stability
  - Crank-Nicolson
  - Stiffness
math: false
---

## Introduction

이번 글에서는 이전 값을 여러 개 사용하는 Linear Multistep Methods를 정리하고, Consistency와 Zero-Stability가 수렴에 함께 필요하다는 것을 확인한다. 또한 Lab 4의 Adams, Stiffness, Shooting 실습을 코드와 결과 그래프로 검증한다. `Lec-2_3_scheme.pdf`의 10–11쪽과 `Lec-4_scheme.pdf`의 10–14쪽을 연결하고, Lab Handout Exercise 5의 Crank–Nicolson 스킴을 안정성 함수로 해석한다.

## 1. Linear Multistep Methods

**One-step method**는 다음 값을 계산할 때 바로 이전 근삿값만 사용하지만, **Multistep method**는 `yₙ`, `yₙ₋₁` 등 여러 이전 상태를 활용한다. 일반적인 k-step 선형 다단계법은 다음과 같다.

![Linear multistep method](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Csum_%7Bj%3D0%7D%5E%7Bk%7D%5Calpha_jy_%7Bn%2Bj%7D%3Dh%5Csum_%7Bj%3D0%7D%5E%7Bk%7D%5Cbeta_j%20f%28t_%7Bn%2Bj%7D%2Cy_%7Bn%2Bj%7D%29)

현재 단계의 우변 함수 평가 `f(tₙ₊ₖ,yₙ₊ₖ)`가 없다면 Explicit, 포함하면 Implicit이다.

### 1.1 Adams-Bashforth (Explicit)

대표적인 2-step Adams–Bashforth는 이전 두 구간의 도함숫값을 이용한다.

![AB2](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7Dy_%7Bn%2B2%7D%3Dy_%7Bn%2B1%7D%2B%5Cfrac%20h2%283f_%7Bn%2B1%7D-f_n%29)

### 1.2 Adams-Moulton (Implicit)

1-step Adams–Moulton은 Trapezoidal Rule과 같다.

![Adams Moulton trapezoidal](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7Dy_%7Bn%2B1%7D%3Dy_n%2B%5Cfrac%20h2%28f_n%2Bf_%7Bn%2B1%7D%29)

### 1.3 Backward Differentiation Formula (BDF)

BDF-1은 바로 Backward Euler이다.

![BDF1 backward Euler](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7Dy_%7Bn%2B1%7D%3Dy_n%2Bh%20f%28t_%7Bn%2B1%7D%2Cy_%7Bn%2B1%7D%29)

BDF-2는 다음과 같다.

![BDF2](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7Dy_%7Bn%2B2%7D-%5Cfrac43%20y_%7Bn%2B1%7D%2B%5Cfrac13%20y_n%3D%5Cfrac%7B2h%7D%7B3%7D%20f%28t_%7Bn%2B2%7D%2Cy_%7Bn%2B2%7D%29)

### 1.4 Predictor-Corrector Methods

계산하기 쉬운 Explicit 방법으로 예측값을 얻고, 그 예측값을 사용해 Implicit 보정 공식을 평가한다. 강의에 제시된 Adams–Bashforth–Moulton 방식은 이 원리를 따른다.

## 2. Consistency and Zero-Stability

다단계법의 특성다항식을 다음과 같이 정의한다.

![Characteristic polynomials](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Crho%28%5Cxi%29%3D%5Csum_%7Bj%3D0%7D%5E%7Bk%7D%5Calpha_j%5Cxi%5Ej%2C%5Cqquad%5Csigma%28%5Cxi%29%3D%5Csum_%7Bj%3D0%7D%5E%7Bk%7D%5Cbeta_j%5Cxi%5Ej)

1차 이상 Consistency를 가지려면 다음 두 조건이 필요하다.

![Consistency conditions](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Crho%281%29%3D0%2C%5Cqquad%5Crho%5Cprime%281%29%3D%5Csigma%281%29)

**Zero-Stability의 Root Condition**은 `ρ(ξ)=0`의 모든 근이 `|ξ|≤1`이며, 절댓값이 1인 근은 단순근이어야 한다는 것이다.

**Dahlquist Equivalence Theorem:** 적절한 선형 다단계 방법에서 `Consistency + Zero-Stability`는 수렴에 필요충분한 조건이다.

### 2.1 Connection to Lab 1 — Leap-Frog

앞서 실험한 `yₙ₊₁=yₙ₋₁−2hyₙ`은 `y′=−y`의 Leap-Frog 차분식이다. `h`를 고정하고 이 ODE에 대입하면 특성근이 `−h±√(1+h²)`이다. 음수 근의 절댓값이 1보다 커서 장기 계산에서 기생 모드가 커질 수 있다.

다만 **엄밀하게 구별해야 할 점**이 있다. Leap-Frog 자체의 첫 번째 특성다항식 `ρ(ξ)=ξ²−1`은 근 `±1`이 모두 단순근이므로 **Zero-Stability의 Root Condition을 만족한다.** 따라서 이 현상을 “Leap-Frog가 zero-stable하지 않다”로 설명하면 잘못이다. 여기서 관찰한 것은 주어진 유한한 `h`에서 발생하는 절대 안정성 문제와 기생 모드의 증폭이다. `T=1`에서 `h→0`일 때의 수렴과 `T`가 길어졌을 때의 불안정성은 구별해야 한다.

![Long time leap frog instability](/assets/images/ode_lab1_leapfrog_long.png)

*Figure 1. Lab Notebook 실험. 동일한 h를 오래 사용하면 기생 진동 모드가 증폭되어 정확해와 크게 달라질 수 있다.*

## 3. Absolute Stability and the Test Equation

수치 방법의 안정성을 살펴볼 때 흔히 다음 테스트 방정식을 사용한다.

![Dahlquist test equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7Dy%5Cprime%3D%5Clambda%20y%2C%5Cqquad%5CRe%28%5Clambda%29%3C0%2C%5Cqquad%20z%3Dh%5Clambda)

정확해는 감소하지만 수치해는 한 단계마다 `R(z)`배 변하므로 일반적으로 `|R(z)|≤1` 조건이 안정성 판단의 기준이 된다(경계점 조건은 방법에 따라 따로 살펴야 한다).

### 3.1 Amplification Factors

| Method | Amplification Factor | Stability for Real Negative z |
|---|---|---|
| Forward Euler | `(1 + z)` | `−2 ≤ z ≤ 0` |
| Backward Euler | `1 / (1 − z)` | 모든 `z < 0` |
| Trapezoidal / Crank–Nicolson | `(1 + z/2) / (1 − z/2)` | 모든 `z < 0` |

![Euler and trapezoidal stability comparison](/assets/images/ode_schemes_stability_comparison.png)

*Figure 2. 추가 Python 시각화. 음의 실수축에서 Forward Euler는 안정 구간이 제한되지만 Backward Euler와 Trapezoidal은 전체 구간에서 증폭률이 1을 넘지 않는다.*

**Python Implementation**

```python
import numpy as np
import matplotlib.pyplot as plt

z = np.linspace(-5, 0, 500)
R_forward = 1 + z
R_backward = 1 / (1 - z)
R_trapezoidal = (1 + z/2) / (1 - z/2)

plt.plot(z, np.abs(R_forward), label="Forward Euler")
plt.plot(z, np.abs(R_backward), label="Backward Euler")
plt.plot(z, np.abs(R_trapezoidal), label="Trapezoidal")
plt.axhline(1, color="gray", linestyle="--")
plt.xlabel("z = h lambda")
plt.ylabel("|R(z)|")
plt.ylim(0, 4)
plt.legend()
plt.grid(True)
plt.show()
```

## 4. Crank-Nicolson / Trapezoidal Scheme — Lab Handout Exercise 5

### Lab Handout Exercise 5 — Crank–Nicolson / Trapezoidal Scheme

마지막 문제는 `y'=λy`, `y(0)=y₀`, 최종시간 `T`, 간격 `h=T/n`이다. 양 끝 시점의 변화율을 평균해 업데이트한다.

![Trapezoidal scheme](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7Dy_%7Bj%2B1%7D%3Dy_j%2B%5Cfrac%20h2%5Cleft%5B%5Clambda%20y_j%2B%5Clambda%20y_%7Bj%2B1%7D%5Cright%5D)

![Trapezoidal amplification factor](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7DR%28h%5Clambda%29%3D%5Cfrac%7B1%2Bh%5Clambda%2F2%7D%7B1-h%5Clambda%2F2%7D)

![Trapezoidal closed form](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7Dy_n%3Dy_0%5Cleft%28%5Cfrac%7B1%2B%5Clambda%20T%2F%282n%29%7D%7B1-%5Clambda%20T%2F%282n%29%7D%5Cright%29%5En)

`n`이 충분히 클 때 로그를 전개하면 `log(yₙ/y₀) → λT`에 해당하는 증폭인자 극한을 얻는다. 따라서 `yₙ → y₀e^(λT)`이며 정확한 해와 일치한다. 분모 `1−hλ/2=0`인 경우에는 업데이트를 정의할 수 없다는 점도 확인해야 한다.

이는 `R(z)=(1+z/2)/(1-z/2)` 형태이며 `n=T/h`번 적용하면 `yₙ = R(λT/n)ⁿ y₀`이다. `n→∞`에서 로그 전개를 사용하면 `n ln R(λT/n)→λT`이고, 결국 정확해 `y(T)=e^(λT)y₀`에 수렴한다.

## 5. A-Stability, L-Stability and Stiff Equations

**A-Stability**는 복소평면에서 실수부가 음수인 모든 `z`를 절대안정영역에 포함한다는 뜻이다. Backward Euler와 Trapezoidal 방법은 A-stable이다. 하지만 **L-Stability**는 여기에 더해 `z→−∞`에서 `R(z)→0`이라는 조건까지 요구한다.

- Backward Euler: `R(z)=1/(1−z)→0`이므로 L-stable.
- Trapezoidal: `R(z)=(1+z/2)/(1−z/2)→−1`이므로 A-stable이지만 L-stable은 아님.

Stiff 문제에서는 매우 빠르게 감소하는 모드 때문에 Explicit 방법의 시간 간격이 정확도보다 **안정성 제한**에 의해 매우 작아질 수 있다. 따라서 안정성 성질은 계산 효율과 직접 연결된다.

## 6. Further Topics — Boundary Value Problems

`Lec-4_scheme.pdf`의 15–16쪽에서는 Boundary Value Problem, Shooting Method, Finite Difference Method 및 삼대각행렬을 푸는 Thomas Algorithm을 소개한다. 이번 글의 Section 6에서 BVP와 Shooting, FDM 및 Thomas Algorithm을 상세히 정리하고, 아래 Lab 4 실습에서 Shooting을 구현한다.


### Lab 2 Connection — Leap-Frog Stability

`LAB_2.ipynb`의 Exercise 3에서는 Leap-Frog의 국소절단오차와 특성근을 직접 확인한다. `T = 1`에서는 `h`를 작게 할수록 전역오차가 줄어드는 결과가 관찰된다. 그러나 `h`를 고정하고 장시간 적분하면 `|−h − sqrt(1+h²)| > 1`인 모드가 증폭된다. 이는 **고정된 구간에서 h → 0인 수렴 문제**와 **고정된 h에서 장시간 적분할 때의 절대 안정성 문제**가 다르다는 점을 강조한다. 실제 실행 코드와 세 그래프는 앞선 Euler Schemes 글의 Lab 2 Exercise 3에 배치했다.

### Lab 3 Connection — Runge-Kutta for ODE Systems

`LAB_3.ipynb`에서는 벡터 RK4로 강제진동계 및 결합된 스프링-질량계를 계산한다. 이들은 **one-step Runge–Kutta를 연립 ODE에 적용하는 실습**이므로 자세한 구현과 그래프는 앞선 Taylor and Runge-Kutta Methods 글에 배치했다. 이 글에서는 반복적으로 사용되는 안정성 및 수렴차수 개념과의 연결만 정리한다.



## Lab 4 — Computational Experiments: Multistep Methods, Stiffness, and Shooting

이번 실습에서는 앞에서 정리한 선형 다단계법, 강성 문제의 절대안정성, 경계값 문제의 Shooting Method를 실제 수치 계산으로 연결한다. 기본 공식은 Lecture 2–3과 Lecture 4를 따르며, **Exercise 1은 강의에 주어진 IVP**, Exercise 2·3은 **Lab 4 원문에 구체적인 방정식과 경계조건이 없어 검증 가능한 예제를 선정**했다. 그래프와 출력값은 `LAB_4.ipynb`의 각 코드 셀에서 생성했다.

### Lab 4 — Exercise 1(a): AB4와 AM4 구현

> **Lab 4 · Exercise 1(a) — Problem**
>
> Implement k-step explicit Adams-Bashforth (AB) and implicit Adams-Moulton (AM) methods to solve a given first-order Initial Value Problem (IVP).
{: .prompt-tip }

Lecture 2–3의 예제 `y'(t)=t-y(t)`, `y(0)=1`, `0 ≤ t ≤ 1`, `h=0.2`를 사용한다. 정확해는 `y(t)=t-1+2exp(-t)`이다. AB4는 이전 네 도함숫값만 사용하는 explicit 방법이며, AM4는 새 시점 도함수를 포함하는 implicit 방법이다. 시작값 `y₁,y₂,y₃`은 RK4로 계산했다. 4차 AM 공식은 **3-step, order-4** 방법이라는 점에서 AB4의 4-step과 구분된다.

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y_%7Bn%2B1%7D%5E%7BAB4%7D%3Dy_n%2B%5Cfrac%7Bh%7D%7B24%7D%2855f_n-59f_%7Bn-1%7D%2B37f_%7Bn-2%7D-9f_%7Bn-3%7D%29)

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y_%7Bn%2B1%7D%5E%7BAM4%7D%3Dy_n%2B%5Cfrac%7Bh%7D%7B24%7D%289f_%7Bn%2B1%7D%2B19f_n-5f_%7Bn-1%7D%2Bf_%7Bn-2%7D%29)

![Figure 3. AB4와 AM4의 수치해 및 정확해 비교.](/assets/images/ode_lab4_1a_1.png)

*Figure 3. AB4와 AM4의 수치해 및 정확해 비교.*

**Python 구현**

```python
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar

# IVP from Lecture 2_3, page 2
# y' = t - y, y(0) = 1, 0 <= t <= 1, h = 0.2

def f(t, y):
    return t - y

def exact(t):
    return t - 1 + 2 * np.exp(-t)

t0, T = 0.0, 1.0
y0 = 1.0
h = 0.2

N = int(round((T - t0) / h))
t = np.linspace(t0, T, N + 1)

# RK4 starting method
def rk4_step(t, y, h):
    k1 = f(t, y)
    k2 = f(t + h/2, y + h*k1/2)
    k3 = f(t + h/2, y + h*k2/2)
    k4 = f(t + h, y + h*k3)
    return y + h*(k1 + 2*k2 + 2*k3 + k4)/6

y_ab = np.zeros(N + 1)
y_am = np.zeros(N + 1)

y_ab[0] = y_am[0] = y0

# Compute y1, y2, y3 using RK4
for n in range(3):
    y_ab[n+1] = rk4_step(t[n], y_ab[n], h)
    y_am[n+1] = y_ab[n+1]

# Adams-Bashforth 4-step (explicit)
for n in range(3, N):
    y_ab[n+1] = y_ab[n] + h/24 * (
        55*f(t[n], y_ab[n])
        - 59*f(t[n-1], y_ab[n-1])
        + 37*f(t[n-2], y_ab[n-2])
        - 9*f(t[n-3], y_ab[n-3])
    )

# Adams-Moulton, fourth-order (implicit)
for n in range(3, N):
    known = y_am[n] + h/24 * (
        19*f(t[n], y_am[n])
        - 5*f(t[n-1], y_am[n-1])
        + f(t[n-2], y_am[n-2])
    )

    def residual(z):
        return z - known - 9*h/24*f(t[n+1], z)

    sol = root_scalar(
        residual,
        bracket=[-10, 10],
        method="brentq"
    )

    y_am[n+1] = sol.root

# Exact solution and errors
y_exact = exact(t)

print(" t       Exact        AB4          AM4       AB Error    AM Error")
for i in range(N + 1):
    print(
        f"{t[i]:.1f}   {y_exact[i]:.8f}   "
        f"{y_ab[i]:.8f}   {y_am[i]:.8f}   "
        f"{abs(y_ab[i]-y_exact[i]):.2e}   "
        f"{abs(y_am[i]-y_exact[i]):.2e}"
    )

print("\nMaximum errors:")
print(f"AB4: {np.max(np.abs(y_ab - y_exact)):.6e}")
print(f"AM4: {np.max(np.abs(y_am - y_exact)):.6e}")

# Plot numerical and exact solutions
plt.figure(figsize=(9, 5))
plt.plot(t, y_exact, "k-", label="Exact solution")
plt.plot(t, y_ab, "o--", label="Adams-Bashforth 4")
plt.plot(t, y_am, "s--", label="Adams-Moulton 4")

plt.xlabel("t")
plt.ylabel("y(t)")
plt.title("Exercise 1(a): AB4 and AM4 Methods")
plt.legend()
plt.grid(True)
plt.show()
```

### Lab 4 — Exercise 1(b): 국소 절단오차와 수렴차수

> **Lab 4 · Exercise 1(b) — Problem**
>
> Compute the Local Truncation Error (LTE) empirically for various time step sizes (h) and verify the theoretical order of accuracy by analyzing the error ratio between grid spacings h and 2h.
{: .prompt-tip }

정확해를 차분식에 직접 대입해 `h`로 정규화한 LTE를 계산했다. 각 격자에서 최대 절댓값을 구하고 `h`를 절반으로 줄여 오차비를 비교한다. **국소 잔차의 차수**와 전체 구간에서 누적되는 **전역오차의 차수**는 개념적으로 구별해야 한다.

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20p_%7B%5Cmathrm%7Bobs%7D%7D%3D%5Cfrac%7B%5Clog%28E%282h%29%2FE%28h%29%29%7D%7B%5Clog%202%7D%2C%5Cqquad%20E%28h%29%3D%5Cmax_n%7C%5Ctau_%7Bn%2B1%7D%28h%29%7C)

AB4와 AM4의 정규화된 LTE는 이론적으로 `O(h⁴)`이며, 작은 `h`에서 오차비는 16, 관측차수는 4에 가까워져야 한다.

![Figure 4. 시간 간격별 AB4·AM4 국소절단오차와 4차 기준선.](/assets/images/ode_lab4_1b_1.png)

*Figure 4. 시간 간격별 AB4·AM4 국소절단오차와 4차 기준선.*

**Python 구현**

```python
import numpy as np
import matplotlib.pyplot as plt

# IVP from Lecture 2_3
def f(t, y):
    return t - y

def exact(t):
    return t - 1 + 2*np.exp(-t)

# Lecture example h = 0.2, followed by grid refinement
h_values = np.array([0.2, 0.1, 0.05, 0.025, 0.0125])

ab_errors = []
am_errors = []

for h in h_values:
    t = np.arange(0, 1 + h/2, h)
    y = exact(t)
    fy = f(t, y)

    lte_ab = []
    lte_am = []

    for n in range(3, len(t)-1):

        # AB4 local truncation error
        ab_residual = (
            (y[n+1] - y[n])/h
            - (
                55*fy[n] - 59*fy[n-1]
                + 37*fy[n-2] - 9*fy[n-3]
              )/24
        )

        # AM4 local truncation error
        am_residual = (
            (y[n+1] - y[n])/h
            - (
                9*fy[n+1] + 19*fy[n]
                - 5*fy[n-1] + fy[n-2]
              )/24
        )

        lte_ab.append(abs(ab_residual))
        lte_am.append(abs(am_residual))

    ab_errors.append(max(lte_ab))
    am_errors.append(max(lte_am))

ab_errors = np.array(ab_errors)
am_errors = np.array(am_errors)

print("h        AB4 LTE       AB Ratio   AB Order    "
      "AM4 LTE       AM Ratio   AM Order")

for i, h in enumerate(h_values):
    if i == 0:
        print(f"{h:<8.4f} {ab_errors[i]:.4e}   ---        ---       "
              f"{am_errors[i]:.4e}   ---        ---")
    else:
        r_ab = ab_errors[i-1] / ab_errors[i]
        r_am = am_errors[i-1] / am_errors[i]
        p_ab = np.log2(r_ab)
        p_am = np.log2(r_am)

        print(f"{h:<8.4f} {ab_errors[i]:.4e}   "
              f"{r_ab:7.3f}    {p_ab:7.3f}    "
              f"{am_errors[i]:.4e}   "
              f"{r_am:7.3f}    {p_am:7.3f}")

plt.figure(figsize=(8,5))
plt.loglog(h_values, ab_errors, 'o-', label='AB4 LTE')
plt.loglog(h_values, am_errors, 's-', label='AM4 LTE')
plt.loglog(
    h_values,
    ab_errors[0]*(h_values/h_values[0])**4,
    'k--', label='Reference: O(h^4)'
)
plt.xlabel('Step size h')
plt.ylabel('Maximum absolute LTE')
plt.title('Exercise 1(b): Local Truncation Error')
plt.legend()
plt.grid(True, which='both')
plt.show()
```

### Lab 4 — Exercise 1(c): PECE Predictor–Corrector

> **Lab 4 · Exercise 1(c) — Problem**
>
> Program a Predictor-Corrector (PC) sequence in PECE mode by pairing an explicit AB predictor to calculate an initial guess with an implicit AM corrector.
{: .prompt-tip }

**Predict → Evaluate → Correct → Evaluate** 순서로 AB4 예측값을 만들고 AM4 공식을 한 번 적용한다. 완전한 implicit AM 해를 매 단계 구하는 것과는 다른 계산이며, 마지막 Evaluate는 보정된 새 상태의 도함수를 평가하는 단계다.

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y_%7Bn%2B1%7D%5E%7B%2A%7D%3Dy_n%2B%5Cfrac%7Bh%7D%7B24%7D%2855f_n-59f_%7Bn-1%7D%2B37f_%7Bn-2%7D-9f_%7Bn-3%7D%29)

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y_%7Bn%2B1%7D%3Dy_n%2B%5Cfrac%7Bh%7D%7B24%7D%289f%28t_%7Bn%2B1%7D%2Cy_%7Bn%2B1%7D%5E%7B%2A%7D%29%2B19f_n-5f_%7Bn-1%7D%2Bf_%7Bn-2%7D%29)

![Figure 5. AB4, AM4 및 PECE 수치해 비교.](/assets/images/ode_lab4_1c_1.png)

*Figure 5. AB4, AM4 및 PECE 수치해 비교.*

![Figure 6. 각 방법의 절대오차 비교.](/assets/images/ode_lab4_1c_2.png)

*Figure 6. 각 방법의 절대오차 비교.*

**Python 구현**

```python
import numpy as np
import matplotlib.pyplot as plt

# IVP from Lecture 2_3
def f(t, y):
    return t - y

def exact(t):
    return t - 1 + 2*np.exp(-t)

# RK4 starting method
def rk4_step(t, y, h):
    k1 = f(t, y)
    k2 = f(t + h/2, y + h*k1/2)
    k3 = f(t + h/2, y + h*k2/2)
    k4 = f(t + h, y + h*k3)
    return y + h*(k1 + 2*k2 + 2*k3 + k4)/6

# Initial conditions from the lecture
t0, T = 0.0, 1.0
y0 = 1.0
h = 0.2

N = int(round((T-t0)/h))
t = np.linspace(t0, T, N+1)

y_ab = np.zeros(N+1)
y_am = np.zeros(N+1)
y_pc = np.zeros(N+1)

y_ab[0] = y_am[0] = y_pc[0] = y0

# Starting values y1, y2, y3 using RK4
for n in range(3):
    y_pc[n+1] = rk4_step(t[n], y_pc[n], h)

y_ab[:4] = y_pc[:4]
y_am[:4] = y_pc[:4]

for n in range(3, N):

    # Adams-Bashforth 4-step
    y_ab[n+1] = y_ab[n] + h/24 * (
        55*f(t[n], y_ab[n])
        -59*f(t[n-1], y_ab[n-1])
        +37*f(t[n-2], y_ab[n-2])
        -9*f(t[n-3], y_ab[n-3])
    )

    # Implicit Adams-Moulton (solved analytically
    # for this particular linear IVP)
    known = y_am[n] + h/24 * (
        19*f(t[n], y_am[n])
        -5*f(t[n-1], y_am[n-1])
        +f(t[n-2], y_am[n-2])
    )

    y_am[n+1] = (
        known + (9*h/24)*t[n+1]
    ) / (1 + 9*h/24)

    # P: Predict with AB4
    y_predict = y_pc[n] + h/24 * (
        55*f(t[n], y_pc[n])
        -59*f(t[n-1], y_pc[n-1])
        +37*f(t[n-2], y_pc[n-2])
        -9*f(t[n-3], y_pc[n-3])
    )

    # E: Evaluate predicted derivative
    f_predict = f(t[n+1], y_predict)

    # C: Correct with AM4
    y_correct = y_pc[n] + h/24 * (
        9*f_predict
        +19*f(t[n], y_pc[n])
        -5*f(t[n-1], y_pc[n-1])
        +f(t[n-2], y_pc[n-2])
    )

    y_pc[n+1] = y_correct

    # E: Evaluate corrected derivative
    f_correct = f(t[n+1], y_pc[n+1])

# Compare results
y_exact = exact(t)

print(" t      Exact        AB4          AM4          PECE")
for i in range(N+1):
    print(
        f"{t[i]:.1f}   {y_exact[i]:.8f}   "
        f"{y_ab[i]:.8f}   {y_am[i]:.8f}   "
        f"{y_pc[i]:.8f}"
    )

print("\nMaximum absolute errors:")
for name, values in [
    ("AB4", y_ab),
    ("AM4", y_am),
    ("PECE", y_pc)
]:
    error = np.max(np.abs(values-y_exact))
    print(f"{name}: {error:.6e}")

# Plot solutions
plt.figure(figsize=(9, 5))
plt.plot(t, y_exact, 'k-', label='Exact')
plt.plot(t, y_ab, 'o--', label='AB4')
plt.plot(t, y_am, 's--', label='AM4')
plt.plot(t, y_pc, 'd--', label='PECE')

plt.xlabel('t')
plt.ylabel('y(t)')
plt.title('Exercise 1(c): AB4-AM4 Predictor-Corrector')
plt.legend()
plt.grid(True)
plt.show()

# Plot absolute errors
plt.figure(figsize=(9, 5))
plt.semilogy(t, np.abs(y_ab-y_exact), 'o-', label='AB4')
plt.semilogy(t, np.abs(y_am-y_exact), 's-', label='AM4')
plt.semilogy(t, np.abs(y_pc-y_exact), 'd-', label='PECE')

plt.xlabel('t')
plt.ylabel('Absolute error')
plt.title('Exercise 1(c): Error Comparison')
plt.legend()
plt.grid(True)
plt.show()
```

### Lab 4 — Exercise 2(a): Jacobian과 Stiffness Ratio

> **Lab 4 · Exercise 2(a) — Problem**
>
> Write a small function or program to calculate the stiffness ratio of a given system of ODEs by evaluating the eigenvalues of its Jacobian matrix.
{: .prompt-tip }

Lab 4에서 특정 강성 시스템을 지정하지 않아, 느린 감쇠와 빠른 감쇠를 동시에 가진 선형 연립계를 선택했다.

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cmathbf%7By%7D%5Cprime%3DA%5Cmathbf%7By%7D%2C%5Cqquad%20A%3D%5Cbegin%7Bpmatrix%7D-1%260%5C%5C0%26-1000%5Cend%7Bpmatrix%7D%2C%5Cqquad%20%5Cmathbf%7By%7D%280%29%3D%5Cbegin%7Bpmatrix%7D1%5C%5C1%5Cend%7Bpmatrix%7D)

Jacobian은 상수행렬 `A`이고 고유값은 `−1`, `−1000`이다. 따라서 강성비는 `S=1000`이다. 서로 다른 시간척도를 가진다는 의미이며, 이것만으로 특정 수치해법의 안정성을 대신 판단할 수는 없다.

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20S%3D%5Cfrac%7B%5Cmax_i%7C%5Coperatorname%7BRe%7D%28%5Clambda_i%29%7C%7D%7B%5Cmin_i%7C%5Coperatorname%7BRe%7D%28%5Clambda_i%29%7C%7D%3D%5Cfrac%7B1000%7D%7B1%7D%3D1000)

**Python 구현**

```python
import numpy as np

# Stiff ODE system: y' = A y
# Example selected for Exercise 2(a)-(c)
A = np.array([
    [-1.0, 0.0],
    [0.0, -1000.0]
])

y0 = np.array([1.0, 1.0])

def f(t, y):
    return A @ y

# Jacobian of the system
def jacobian(t, y):
    return A

# Calculate stiffness ratio using Jacobian eigenvalues
def stiffness_ratio(J):
    eigenvalues = np.linalg.eigvals(J)

    real_parts = np.abs(np.real(eigenvalues))

    if np.any(real_parts == 0):
        return np.inf, eigenvalues

    ratio = np.max(real_parts) / np.min(real_parts)

    return ratio, eigenvalues

J = jacobian(0.0, y0)

ratio, eigenvalues = stiffness_ratio(J)

print("Jacobian matrix:")
print(J)

print("\nEigenvalues:")
print(eigenvalues)

print(f"\nStiffness ratio: {ratio:.2f}")
```

### Lab 4 — Exercise 2(b): Explicit AB2의 시간 간격 제한

> **Lab 4 · Exercise 2(b) — Problem**
>
> Attempt to solve a highly stiff system using an explicit LMM and observe the severe step-size restrictions (h < maxC|λi | ) dictated by absolute stability requirements.
{: .prompt-tip }

같은 시스템에 AB2를 적용하고 `h=0.0005`와 `h=0.01`을 비교했다. 빠른 성분의 고유값이 `−1000`이고, AB2의 음의 실수축 안정구간은 `−1 ≤ hλ ≤ 0`이므로 안정성 경계는 `h=0.001`이다. 따라서 `h=0.01`은 안정영역 밖이다. 초기 한 단계는 강의에서 소개한 RK4로 시작한다.

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cmathbf%7By%7D_%7Bn%2B1%7D%3D%5Cmathbf%7By%7D_n%2B%5Cfrac%20h2%283A%5Cmathbf%7By%7D_n-A%5Cmathbf%7By%7D_%7Bn-1%7D%29)

![Figure 7. AB2의 시간 간격에 따른 빠른 성분 안정성 비교.](/assets/images/ode_lab4_2b_1.png)

*Figure 7. AB2의 시간 간격에 따른 빠른 성분 안정성 비교.*

정확해 `y₂(t)=exp(−1000t)`는 감소하지만, 큰 시간 간격에서는 수치적인 기생 모드가 성장한다. 이는 절대안정성 제약** 때문이다.

**Python 구현**

```python
import numpy as np
import matplotlib.pyplot as plt

# Same stiff system as Exercise 2(a)
A = np.array([[-1.0, 0.0],
              [0.0, -1000.0]])
y0 = np.array([1.0, 1.0])

def f(t, y):
    return A @ y

def exact(t):
    return np.column_stack((np.exp(-t), np.exp(-1000*t)))

# RK4 starting step, following Lecture 2_3
def rk4_step(t, y, h):
    k1 = f(t, y)
    k2 = f(t+h/2, y+h*k1/2)
    k3 = f(t+h/2, y+h*k2/2)
    k4 = f(t+h, y+h*k3)
    return y + h*(k1+2*k2+2*k3+k4)/6

# Explicit Adams-Bashforth 2-step method
def ab2(h, T):
    N = int(round(T/h))
    t = np.linspace(0, T, N+1)
    y = np.zeros((N+1, 2))
    y[0] = y0

    # One starting value from RK4
    y[1] = rk4_step(t[0], y[0], h)

    for n in range(1, N):
        y[n+1] = y[n] + h/2 * (
            3*f(t[n], y[n]) - f(t[n-1], y[n-1])
        )
    return t, y

# Stability restriction for AB2
eigenvalues = np.linalg.eigvals(A)
h_limit = 1 / max(abs(eigenvalues))

print(f"Jacobian eigenvalues: {eigenvalues}")
print(f"AB2 stability boundary: h = {h_limit:.6f}")

# Stable and unstable step sizes
h_values = [0.0005, 0.01]
T = 0.1

fig, axes = plt.subplots(1, 2, figsize=(12, 4))

for ax, h in zip(axes, h_values):
    t, y = ab2(h, T)
    y_true = exact(t)

    ax.plot(t, y[:, 1], 'o-', markersize=3,
            label=f'AB2 (h={h})')
    ax.plot(t, y_true[:, 1], 'k--',
            label='Exact solution')

    ax.set_xlabel('t')
    ax.set_ylabel('y2(t)')
    ax.set_title(f'AB2: h = {h}')
    ax.grid(True)
    ax.legend()

    error = np.max(np.abs(y - y_true))
    print(f"h = {h:.4f}: maximum absolute error = {error:.6e}")
    print(f"         final y2 = {y[-1, 1]:.6e}")

plt.suptitle('Exercise 2(b): Explicit AB2 Stability')
plt.tight_layout()
plt.show()
```

### Lab 4 — Exercise 2(c): BDF-1 및 BDF-2 안정성 비교

> **Lab 4 · Exercise 2(c) — Problem**
>
> Implement Backward Differentiation Formulas (BDF-1 and BDF-2) and evaluate their performance on the stiff system, validating their theoretical advantages as A-stable and L-stable implicit methods.
{: .prompt-tip }

같은 연립계와 시간 간격에서 implicit BDF 방법을 적용했다. BDF-1은 Backward Euler이며, BDF-2는 2차 다단계법이다. BDF-2의 시작값은 한 번의 BDF-1 단계로 계산한다.

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%28I-hA%29%5Cmathbf%7By%7D_%7Bn%2B1%7D%3D%5Cmathbf%7By%7D_n)

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%283I-2hA%29%5Cmathbf%7By%7D_%7Bn%2B1%7D%3D4%5Cmathbf%7By%7D_n-%5Cmathbf%7By%7D_%7Bn-1%7D)

![Figure 8. 작은 시간 간격에서 AB2, BDF-1, BDF-2 비교.](/assets/images/ode_lab4_2c_1.png)

*Figure 8. 작은 시간 간격에서 AB2, BDF-1, BDF-2 비교.*

![Figure 9. 큰 시간 간격에서 명시적 AB2와 암시적 BDF 방법의 비교.](/assets/images/ode_lab4_2c_2.png)

*Figure 9. 큰 시간 간격에서 명시적 AB2와 암시적 BDF 방법의 비교.*

![Figure 10. 느린 모드에서 BDF-1과 BDF-2 비교.](/assets/images/ode_lab4_2c_3.png)

*Figure 10. 느린 모드에서 BDF-1과 BDF-2 비교.*

BDF-1과 BDF-2는 모두 A-stable이므로 해당 음의 실수 고유값에 대해 큰 스텝에서도 안정적이다. BDF-1의 단일 증폭인자와 달리 BDF-2의 안정성은 **두 특성근**으로 분석해야 한다. 두 근의 크기는 강성 극한에서 0으로 향하지만, 안정하다는 사실이 **모든 `h`에서 작은 오차를 보장하지는 않는다.** 특히 BDF-2 시작값을 BDF-1로 구했기 때문에 관측 오차에 시작 절차의 영향이 남을 수 있다.

**Python 구현**

```python
import numpy as np
import matplotlib.pyplot as plt

# Same stiff system as Exercise 2(a) and 2(b)
A = np.array([
    [-1.0, 0.0],
    [0.0, -1000.0]
])

y0 = np.array([1.0, 1.0])
I = np.eye(2)

def f(t, y):
    return A @ y

def exact(t):
    return np.column_stack((
        np.exp(-t),
        np.exp(-1000*t)
    ))

# BDF-1 (Backward Euler)
def bdf1(h, T):
    N = int(round(T/h))
    t = np.linspace(0, T, N+1)
    y = np.zeros((N+1, 2))
    y[0] = y0

    M = I - h*A

    for n in range(N):
        y[n+1] = np.linalg.solve(M, y[n])

    return t, y

# BDF-2
def bdf2(h, T):
    N = int(round(T/h))
    t = np.linspace(0, T, N+1)
    y = np.zeros((N+1, 2))
    y[0] = y0

    # Starting value computed using BDF-1
    if N >= 1:
        y[1] = np.linalg.solve(I-h*A, y[0])

    M = 3*I - 2*h*A

    for n in range(1, N):
        rhs = 4*y[n] - y[n-1]
        y[n+1] = np.linalg.solve(M, rhs)

    return t, y

# Explicit AB2 from Exercise 2(b), for comparison
def rk4_step(t, y, h):
    k1 = f(t, y)
    k2 = f(t+h/2, y+h*k1/2)
    k3 = f(t+h/2, y+h*k2/2)
    k4 = f(t+h, y+h*k3)
    return y + h*(k1+2*k2+2*k3+k4)/6

def ab2(h, T):
    N = int(round(T/h))
    t = np.linspace(0, T, N+1)
    y = np.zeros((N+1, 2))
    y[0] = y0

    if N >= 1:
        y[1] = rk4_step(t[0], y[0], h)

    for n in range(1, N):
        y[n+1] = y[n] + h/2 * (
            3*f(t[n], y[n])
            - f(t[n-1], y[n-1])
        )

    return t, y

# Numerical experiments
T = 0.1
h_values = [0.0005, 0.01]

for h in h_values:
    t, y_ab = ab2(h, T)
    _, y_bdf1 = bdf1(h, T)
    _, y_bdf2 = bdf2(h, T)

    y_true = exact(t)

    print(f"\nStep size h = {h}")
    print("-" * 68)
    print(f"{'Method':<10} {'Max Error':<15} "
          f"{'Final y1':<15} {'Final y2':<15}")

    for name, y in [
        ("AB2", y_ab),
        ("BDF-1", y_bdf1),
        ("BDF-2", y_bdf2)
    ]:
        error = np.max(
            np.abs(y - y_true)
        )

        print(
            f"{name:<10} {error:<15.6e} "
            f"{y[-1,0]:<15.6e} "
            f"{y[-1,1]:<15.6e}"
        )

    print(
        f"{'Exact':<10} {'---':<15} "
        f"{y_true[-1,0]:<15.6e} "
        f"{y_true[-1,1]:<15.6e}"
    )

    # Plot fast component
    plt.figure(figsize=(9, 5))

    plt.plot(t, y_true[:,1], 'k-',
             label='Exact')
    plt.plot(t, y_ab[:,1], 'o--',
             markersize=3, label='AB2')
    plt.plot(t, y_bdf1[:,1], 's-',
             markersize=3, label='BDF-1')
    plt.plot(t, y_bdf2[:,1], 'd-',
             markersize=3, label='BDF-2')

    plt.xlabel('t')
    plt.ylabel('y2(t)')
    plt.title(f'Exercise 2(c): Stability Comparison (h={h})')
    plt.legend()
    plt.grid(True)
    plt.show()

# Compare slow component for a large step size
h = 0.01

t, y_bdf1 = bdf1(h, T)
_, y_bdf2 = bdf2(h, T)
y_true = exact(t)

plt.figure(figsize=(9, 5))
plt.plot(t, y_true[:,0], 'k-', label='Exact')
plt.plot(t, y_bdf1[:,0], 's--', label='BDF-1')
plt.plot(t, y_bdf2[:,0], 'd--', label='BDF-2')

plt.xlabel('t')
plt.ylabel('y1(t)')
plt.title('Slow Component: BDF-1 vs BDF-2')
plt.legend()
plt.grid(True)
plt.show()

# Stability amplification analysis
print("\nStiff-mode amplification for h=0.01")

z = -1000 * 0.01

# BDF-1 amplification factor
R_bdf1 = 1 / (1-z)

# BDF-2 characteristic roots:
# (3 - 2z) r^2 - 4r + 1 = 0
roots_bdf2 = np.roots([3-2*z, -4, 1])

print(f"z = {z}")
print(f"BDF-1 amplification magnitude: {abs(R_bdf1):.6f}")
print(f"BDF-2 root magnitudes: {np.abs(roots_bdf2)}")
print("Both methods damp the stiff mode.")
```

### Lab 4 — Exercise 3(a): Linear Shooting

> **Lab 4 · Exercise 3(a) — Problem**
>
> Solve a linear second-order Dirichlet BVP by expressing the solution as a linear combination of two independent IVPs (one inhomogeneous and one homogeneous).
{: .prompt-tip }

선형 2계 Dirichlet BVP는 **비제차 IVP와 제차 IVP**를 각각 풀어 선형결합한다. 이번 예제는 `y′′=y+x`, `y(0)=0`, `y(1)=1`이고 정확해는 `y(x)=2sinh(x)/sinh(1)−x`다.

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y_1%5Cprime%5Cprime%3Dy_1%2Bx%2C%5Cquad%20y_1%280%29%3D0%2C%5Cquad%20y_1%5Cprime%280%29%3D0)

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y_2%5Cprime%5Cprime%3Dy_2%2C%5Cquad%20y_2%280%29%3D0%2C%5Cquad%20y_2%5Cprime%280%29%3D1)

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y%28x%29%3Dy_1%28x%29%2Bcy_2%28x%29%2C%5Cqquad%20c%3D%5Cfrac%7B1-y_1%281%29%7D%7By_2%281%29%7D)

![Figure 11. 두 IVP를 결합한 Linear Shooting 해와 정확해.](/assets/images/ode_lab4_3a_1.png)

*Figure 11. 두 IVP를 결합한 Linear Shooting 해와 정확해.*

![Figure 12. Linear Shooting의 절대오차.](/assets/images/ode_lab4_3a_2.png)

*Figure 12. Linear Shooting의 절대오차.*

**Python 구현**

```python
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Linear BVP selected for Exercise 3(a)
# y'' = y + x, y(0) = 0, y(1) = 1

a, b = 0.0, 1.0
alpha, beta = 0.0, 1.0

x_grid = np.linspace(a, b, 201)

# IVP 1: Inhomogeneous equation
# y1'' = y1 + x
# y1(0) = 0, y1'(0) = 0
def ivp1(x, Y):
    y, yp = Y
    return [yp, y + x]

sol1 = solve_ivp(
    ivp1, (a, b), [alpha, 0.0],
    t_eval=x_grid,
    rtol=1e-10, atol=1e-12
)

# IVP 2: Homogeneous equation
# y2'' = y2
# y2(0) = 0, y2'(0) = 1
def ivp2(x, Y):
    y, yp = Y
    return [yp, y]

sol2 = solve_ivp(
    ivp2, (a, b), [0.0, 1.0],
    t_eval=x_grid,
    rtol=1e-10, atol=1e-12
)

assert sol1.success and sol2.success

y1 = sol1.y[0]
y2 = sol2.y[0]

# Linear combination coefficient
c = (beta - y1[-1]) / y2[-1]

# Shooting solution
y_shoot = y1 + c*y2

# Exact solution
def exact(x):
    return 2*np.sinh(x)/np.sinh(1) - x

y_exact = exact(x_grid)

# Errors
errors = np.abs(y_shoot - y_exact)

print(f"Shooting coefficient c = {c:.10f}")
print(f"Initial boundary y(0) = {y_shoot[0]:.10f}")
print(f"Final boundary y(1) = {y_shoot[-1]:.10f}")
print(f"Maximum absolute error = {np.max(errors):.6e}")

# Plot the two IVPs and the combined solution
plt.figure(figsize=(9, 5))

plt.plot(x_grid, y1, '--', label='IVP 1: Inhomogeneous')
plt.plot(x_grid, y2, '--', label='IVP 2: Homogeneous')
plt.plot(x_grid, y_shoot, 'b-', label='Linear Shooting')
plt.plot(x_grid, y_exact, 'k:', label='Exact Solution')

plt.xlabel('x')
plt.ylabel('y(x)')
plt.title('Exercise 3(a): Linear Shooting Method')
plt.legend()
plt.grid(True)
plt.show()

# Error plot
plt.figure(figsize=(9, 4))
plt.plot(x_grid, errors)

plt.xlabel('x')
plt.ylabel('Absolute Error')
plt.title('Exercise 3(a): Error of Linear Shooting')
plt.grid(True)
plt.show()
```

### Lab 4 — Exercise 3(b): Nonlinear Shooting과 Newton–Raphson

> **Lab 4 · Exercise 3(b) — Problem**
>
> Solve a non-linear BVP by implementing the Nonlinear Shooting Method paired with a Newton-Raphson root-finding iteration.
{: .prompt-tip }

비선형 예제로 `y′′=2y³`, `y(0)=0.5`, `y(1)=1`을 사용했다. 미지의 초기 기울기 `s=y′(0)`를 Newton 반복으로 찾아 오른쪽 경계잔차 `F(s)=y(1;s)−1`을 0으로 만든다. 비교할 정확해는 `y(x)=1/(2−x)`이므로 정확한 초기 기울기는 `0.25`다. 이 단계에서는 `F′(s)`를 중앙차분으로 근사했다.

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20s_%7Bk%2B1%7D%3Ds_k-%5Cfrac%7BF%28s_k%29%7D%7BF%5Cprime%28s_k%29%7D%2C%5Cqquad%20F%28s%29%3Dy%281%3Bs%29-1)

![Figure 13. 비선형 Shooting 해와 정확해.](/assets/images/ode_lab4_3b_1.png)

*Figure 13. 비선형 Shooting 해와 정확해.*

![Figure 14. Newton 반복에 따른 경계잔차 감소.](/assets/images/ode_lab4_3b_2.png)

*Figure 14. Newton 반복에 따른 경계잔차 감소.*

**Python 구현**

```python
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Nonlinear BVP selected for Exercises 3(b) and 3(c)
# y'' = 2y^3, y(0) = 0.5, y(1) = 1

a, b = 0.0, 1.0
alpha, beta = 0.5, 1.0
x_grid = np.linspace(a, b, 201)

# Parameterized IVP: initial slope = s
def nonlinear_ivp(x, Y):
    y, yp = Y
    return [yp, 2*y**3]

def shoot(s, t_eval=None):
    sol = solve_ivp(
        nonlinear_ivp, (a, b), [alpha, s],
        t_eval=t_eval,
        rtol=1e-10, atol=1e-12
    )
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol

# Boundary residual F(s) = y(1;s) - beta
def F(s):
    return shoot(s).y[0, -1] - beta

# Central finite difference approximation of F'(s)
def dF(s, epsilon=1e-5):
    return (F(s + epsilon) - F(s - epsilon))/(2*epsilon)

# Newton-Raphson iteration
s = 0.2
tol = 1e-10
max_iter = 20

s_history = [s]
residual_history = []

print("Iteration     Initial slope       Boundary residual")

for k in range(max_iter):
    residual = F(s)
    residual_history.append(abs(residual))

    print(f"{k:<13} {s:<19.12f} {residual:.6e}")

    if abs(residual) < tol:
        break

    derivative = dF(s)

    if abs(derivative) < 1e-14:
        raise RuntimeError("Newton derivative is too small.")

    s = s - residual/derivative
    s_history.append(s)
else:
    raise RuntimeError("Newton iteration did not converge.")

# Final numerical solution
sol = shoot(s, t_eval=x_grid)
y_shoot = sol.y[0]

# Exact solution for verification
def exact(x):
    return 1/(2-x)

y_exact = exact(x_grid)

error = np.max(np.abs(y_shoot - y_exact))

print("\nFinal results:")
print(f"Computed initial slope = {s:.12f}")
print(f"Exact initial slope    = {0.25:.12f}")
print(f"Computed y(1)          = {y_shoot[-1]:.12f}")
print(f"Boundary residual      = {abs(y_shoot[-1]-beta):.6e}")
print(f"Maximum absolute error = {error:.6e}")

# Plot numerical and exact solutions
plt.figure(figsize=(9, 5))
plt.plot(x_grid, y_shoot, 'b-', label='Nonlinear Shooting')
plt.plot(x_grid, y_exact, 'k--', label='Exact Solution')
plt.xlabel('x')
plt.ylabel('y(x)')
plt.title('Exercise 3(b): Nonlinear Shooting Method')
plt.legend()
plt.grid(True)
plt.show()

# Plot Newton-Raphson convergence
plt.figure(figsize=(8, 4))
plt.semilogy(
    range(len(residual_history)),
    residual_history, 'o-'
)
plt.xlabel('Newton Iteration')
plt.ylabel('Absolute Boundary Residual')
plt.title('Exercise 3(b): Newton-Raphson Convergence')
plt.grid(True)
plt.show()
```

### Lab 4 — Exercise 3(c): Variational Equation과 정확한 Newton 도함수

> **Lab 4 · Exercise 3(c) — Problem**
>
> Formulate and simultaneously integrate the sensitivity function (variational equation) alongside the parameterized IVP to compute the exact derivative needed for the Newton-Raphson update steps.
{: .prompt-tip }

3(b)와 동일한 BVP에서 민감도 `z(x;s)=∂y(x;s)/∂s`를 정의한다. 원래 방정식을 `s`에 대해 미분하면 **변분방정식**이 얻어진다. 네 개의 1차 ODE를 함께 적분해 `F′(s)=z(1;s)`를 구하므로 Newton 업데이트에서는 유한차분을 사용하지 않는다.

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20z%5Cprime%5Cprime%3D6y%5E2z%2C%5Cqquad%20z%280%29%3D0%2C%5Cquad%20z%5Cprime%280%29%3D1)

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cbegin%7Baligned%7Dy%5Cprime%26%3Dv%2C%26v%5Cprime%26%3D2y%5E3%2C%5C%5Cz%5Cprime%26%3Dw%2C%26w%5Cprime%26%3D6y%5E2z%5Cend%7Baligned%7D)

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20s_%7Bk%2B1%7D%3Ds_k-%5Cfrac%7By%281%3Bs_k%29-1%7D%7Bz%281%3Bs_k%29%7D)

![Figure 15. Sensitivity를 이용한 Nonlinear Shooting과 정확해.](/assets/images/ode_lab4_3c_1.png)

*Figure 15. Sensitivity를 이용한 Nonlinear Shooting과 정확해.*

![Figure 16. 계산된 Sensitivity Function.](/assets/images/ode_lab4_3c_2.png)

*Figure 16. 계산된 Sensitivity Function.*

![Figure 17. Sensitivity 기반 Newton 반복의 수렴.](/assets/images/ode_lab4_3c_3.png)

*Figure 17. Sensitivity 기반 Newton 반복의 수렴.*

유한차분으로 계산한 경계잔차 도함수와 민감도 방정식으로 계산한 도함수를 비교해 검산한다. 여기서 “exact derivative”는 **수학적으로 올바르게 유도된 도함수**를 뜻하며, IVP 자체의 수치 적분 오차가 사라진다는 의미는 아니다.

**Python 구현**

```python
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Nonlinear BVP from Exercise 3(b)
# y'' = 2y^3, y(0) = 0.5, y(1) = 1
# Example selected for the lab demonstration

a, b = 0.0, 1.0
alpha, beta = 0.5, 1.0
x_grid = np.linspace(a, b, 201)

# Augmented system:
# Y = [y, y', z, z']
# z = partial derivative of y with respect to s
#
# Original equation: y'' = 2y^3
# Sensitivity equation: z'' = 6y^2*z

def augmented_system(x, Y):
    y, yp, z, zp = Y

    return [
        yp,
        2*y**3,
        zp,
        6*y**2*z
    ]

# Solve the parameterized IVP and sensitivity equation
# simultaneously for a given initial slope s
def solve_augmented(s, t_eval=None):
    Y0 = [alpha, s, 0.0, 1.0]

    sol = solve_ivp(
        augmented_system,
        (a, b),
        Y0,
        t_eval=t_eval,
        rtol=1e-10,
        atol=1e-12
    )

    if not sol.success:
        raise RuntimeError(sol.message)

    return sol

# Newton-Raphson iteration using sensitivity derivative
s = 0.2
tol = 1e-10
max_iter = 20

s_history = [s]
residual_history = []

print("Iteration     Slope s          F(s)           F'(s)")

for k in range(max_iter):
    sol = solve_augmented(s)

    y_b = sol.y[0, -1]
    z_b = sol.y[2, -1]

    # Boundary residual F(s) = y(1;s) - beta
    F = y_b - beta

    # Derivative F'(s) = z(1;s), from the variational equation
    dF = z_b

    residual_history.append(abs(F))

    print(
        f"{k:<13} {s:<16.12f} "
        f"{F:<14.6e} {dF:.6e}"
    )

    if abs(F) < tol:
        break

    if abs(dF) < 1e-14:
        raise RuntimeError("Sensitivity derivative is too small.")

    # Newton-Raphson update
    s = s - F/dF
    s_history.append(s)

else:
    raise RuntimeError("Newton iteration did not converge.")

# Final numerical solution and sensitivity function
sol = solve_augmented(s, t_eval=x_grid)

y_num = sol.y[0]
z_num = sol.y[2]

# Exact solution for numerical verification
def exact(x):
    return 1/(2-x)

y_exact = exact(x_grid)

# Maximum absolute solution error
max_error = np.max(np.abs(y_num - y_exact))

# Verify the sensitivity derivative at the final shooting parameter.
# The variational equation provides F'(s) = z(1;s).
# A central finite difference is used only for verification,
# not for the Newton-Raphson updates.

epsilon = 1e-5

y_plus = solve_augmented(s + epsilon).y[0, -1]
y_minus = solve_augmented(s - epsilon).y[0, -1]

dF_fd = (y_plus - y_minus)/(2*epsilon)
dF_sensitivity = sol.y[2, -1]

print("\nFinal results:")
print(f"Computed initial slope       = {s:.12f}")
print(f"Exact initial slope          = {0.25:.12f}")
print(f"Computed y(1)                = {y_num[-1]:.12f}")
print(f"Boundary residual            = {abs(y_num[-1]-beta):.6e}")
print(f"Maximum absolute error       = {max_error:.6e}")
print(f"Derivative by sensitivity    = {dF_sensitivity:.12f}")
print(f"Derivative by finite diff.   = {dF_fd:.12f}")
print(f"Derivative difference        = {abs(dF_sensitivity-dF_fd):.6e}")

# Plot 1: Numerical and exact solutions
plt.figure(figsize=(9, 5))

plt.plot(
    x_grid, y_num, 'b-',
    label='Shooting with Sensitivity'
)
plt.plot(
    x_grid, y_exact, 'k--',
    label='Exact Solution'
)

plt.xlabel('x')
plt.ylabel('y(x)')
plt.title('Exercise 3(c): Nonlinear Shooting with Sensitivity')
plt.legend()
plt.grid(True)
plt.show()

# Plot 2: Sensitivity function
plt.figure(figsize=(9, 5))

plt.plot(
    x_grid, z_num, 'r-',
    label='Sensitivity z(x)'
)

plt.xlabel('x')
plt.ylabel('z(x)')
plt.title('Exercise 3(c): Sensitivity Function')
plt.legend()
plt.grid(True)
plt.show()

# Plot 3: Newton-Raphson convergence
plt.figure(figsize=(8, 4))

plt.semilogy(
    range(len(residual_history)),
    residual_history,
    'o-'
)

plt.xlabel('Newton Iteration')
plt.ylabel('Absolute Boundary Residual')
plt.title('Exercise 3(c): Newton-Raphson Convergence')
plt.grid(True)
plt.show()
```

### Lab 4 — 실험 결과에서 확인한 점

- **Adams Methods:** 같은 IVP에서 AB4, AM4 및 PECE의 동작을 비교했고, 정규화한 LTE의 수렴차수가 4에 가까워지는지 확인했다.
- **Stiffness:** 강성비 `1000`인 시스템에서 AB2는 큰 스텝으로 불안정해졌지만 BDF-1과 BDF-2는 강한 감쇠 모드를 안정적으로 처리했다. **안정성과 정확성은 별도로 검증**해야 한다.
- **Shooting Methods:** 선형 문제는 두 IVP의 결합으로, 비선형 문제는 초기 기울기에 대한 Newton 업데이트로 해결했다. 변분방정식을 함께 적분하면 도함수에 대한 유한차분 근사 없이 Newton 반복이 가능하다.


## Detailed Lecture Notes — Linear Multistep Methods, Stability, and BVPs

이 절은 `Lec-2_3_scheme.pdf`의 multistep 도입과 `Lec-4_scheme.pdf`의 **Section 4–6 전체**를 대응시켜 정리한다. 기존 설명에서 간략하게 지나갔던 **계수의 정확한 정의, 일관성 조건의 Taylor 유도, Root Condition, Dahlquist 정리, Shooting, FDM, Thomas 알고리즘**을 포함한다.

### 1.5 General Linear Multistep Formulation

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20k)

-step linear multistep method(LMM)는 앞 단계 값과 도함수 평가값 여러 개를 결합한다. 강의는

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Csum_%7Bj%3D0%7D%5Ek%5Calpha_j%20y_%7Bn%2Bj%7D%3Dh%5Csum_%7Bj%3D0%7D%5Ek%5Cbeta_j%20f%28t_%7Bn%2Bj%7D%2Cy_%7Bn%2Bj%7D%29)


로 정의하며

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Calpha_k%3D1)


로 정규화한다.

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cbeta_k%3D0)


이면 새로운

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f_%7Bn%2Bk%7D)


가 등장하지 않아 **explicit**,

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cbeta_k%5Cne0)


이면 **implicit**이다. 구별해야 할 점은

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20k)


는 한 단계 업데이트에서 참조하는 과거 정보의 수이며,

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20p)

는 오차의 수렴차수라는 것이다.

### 1.6 Characteristic Polynomials and Consistency

LMM의 특성다항식은

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Crho%28%5Czeta%29%3D%5Csum_%7Bj%3D0%7D%5Ek%5Calpha_j%5Czeta%5Ej%2C%5Cqquad%5Csigma%28%5Czeta%29%3D%5Csum_%7Bj%3D0%7D%5Ek%5Cbeta_j%5Czeta%5Ej)

이다. 정확해를 대입한 강의의 정규화된 local truncation error는

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cmathcal%20T_%7Bn%2Bk%7D%3D%5Cfrac1h%5Cleft%28%5Csum_%7Bj%3D0%7D%5Ek%5Calpha_j%20y%28t_%7Bn%2Bj%7D%29-h%5Csum_%7Bj%3D0%7D%5Ek%5Cbeta_j%20y%27%28t_%7Bn%2Bj%7D%29%5Cright%29)


이다. 공통 기준점

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20t_n)

에서 Taylor 전개하면

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y%28t_n%2Bjh%29%3Dy%2Bjhy%27%2B%5Cfrac%7Bj%5E2h%5E2%7D%7B2%21%7Dy%27%27%2B%5Ccdots)

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y%27%28t_n%2Bjh%29%3Dy%27%2Bjhy%27%27%2B%5Cfrac%7Bj%5E2h%5E2%7D%7B2%21%7Dy%5E%7B%283%29%7D%2B%5Ccdots)

이므로

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cmathcal%20T_%7Bn%2Bk%7D%3D%5Cfrac%7BC_0%7D%7Bh%7Dy%2BC_1y%27%2BC_2h%20y%27%27%2B%5Ccdots%2BC_p%20h%5E%7Bp-1%7Dy%5E%7B%28p%29%7D%2BO%28h%5Ep%29)

이다. 계수는

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20C_0%3D%5Csum_%7Bj%3D0%7D%5Ek%5Calpha_j%3D%5Crho%281%29)

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20C_1%3D%5Csum_%7Bj%3D0%7D%5Ek%28j%5Calpha_j-%5Cbeta_j%29%3D%5Crho%27%281%29-%5Csigma%281%29)

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20C_q%3D%5Csum_%7Bj%3D0%7D%5Ek%5Cleft%28%5Cfrac%7Bj%5Eq%7D%7Bq%21%7D%5Calpha_j-%5Cfrac%7Bj%5E%7Bq-1%7D%7D%7B%28q-1%29%21%7D%5Cbeta_j%5Cright%29%2C%5Cquad%20q%5Cge2)

이다. 모든 적절한 매끄러운 해에 대해 consistency를 얻기 위한 조건은

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20C_0%3DC_1%3D0)

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Crho%281%29%3D0%2C%5Cqquad%5Crho%27%281%29%3D%5Csigma%281%29)


이다.

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20C_0%3D%5Ccdots%3DC_p%3D0)


이고

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20C_%7Bp%2B1%7D%5Cne0)


이면 강의에서 정의하는 정확도 차수

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20p)

가 된다.

### 1.7 Adams-Bashforth Schemes


Adams–Bashforth는 explicit multistep 방법이다. 도함수

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f)

를 과거 격자에서 보간한 뒤 적분해 만든다. 강의에 제시된 대표적인 공식은 다음과 같다.

**AB1 / Forward Euler**

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y_%7Bn%2B1%7D%3Dy_n%2Bh%20f_n)

**AB2**

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y_%7Bn%2B2%7D%3Dy_%7Bn%2B1%7D%2B%5Cfrac%20h2%283f_%7Bn%2B1%7D-f_n%29)

**AB4**

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y_%7Bn%2B4%7D%3Dy_%7Bn%2B3%7D%2B%5Cfrac%20h%7B24%7D%2855f_%7Bn%2B3%7D-59f_%7Bn%2B2%7D%2B37f_%7Bn%2B1%7D-9f_n%29)


여기서

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f_j%3Df%28t_j%2Cy_j%29)


이다. 높은 차수의 방법은 이전 값이 여러 개 필요하므로

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y_1%2Cy_2%2Cy_3)

등을 RK4 같은 단일 단계 방법으로 먼저 계산할 수 있다.

### 1.8 Adams-Moulton Schemes


Adams–Moulton은 새 시점의

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f_%7Bn%2Bk%7D)

까지 포함하는 implicit 방법이다.

**AM1 / Trapezoidal Rule**

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y_%7Bn%2B1%7D%3Dy_n%2B%5Cfrac%20h2%28f_%7Bn%2B1%7D%2Bf_n%29)

**AM2 (order 3)**

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y_%7Bn%2B2%7D%3Dy_%7Bn%2B1%7D%2B%5Cfrac%20h%7B12%7D%285f_%7Bn%2B2%7D%2B8f_%7Bn%2B1%7D-f_n%29)

**AM3 (order 4)**

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y_%7Bn%2B3%7D%3Dy_%7Bn%2B2%7D%2B%5Cfrac%20h%7B24%7D%289f_%7Bn%2B3%7D%2B19f_%7Bn%2B2%7D-5f_%7Bn%2B1%7D%2Bf_n%29)

강의는 AB와 AM을 predictor–corrector로 묶는 **PECE** 절차도 제시한다: Predict, Evaluate, Correct, Evaluate. 먼저 explicit 예측값을 만들고, 그 값에서 도함수를 계산한 뒤 implicit 교정식을 평가하는 순서다.

### 1.9 Backward Differentiation Formulas


BDF는 새 시점의 도함수

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f_%7Bn%2Bk%7D)

하나만 우변에 사용하는 implicit LMM이다.

**BDF1 / Backward Euler**

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y_%7Bn%2B1%7D-y_n%3Dhf_%7Bn%2B1%7D)

**BDF2**

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y_%7Bn%2B2%7D-%5Cfrac43y_%7Bn%2B1%7D%2B%5Cfrac13y_n%3D%5Cfrac%7B2h%7D3f_%7Bn%2B2%7D)

특히 BDF1이 Lab Handout의 Backward Euler와 동일한 공식이라는 점을 확인한다. BDF 방법은 강의의 stiffness 논의와 연결된다.

### 2.2 Zero-Stability: The Root Condition

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20h%5Cto0)


일 때 한 단계에서 발생한 작은 오차가 반복 과정에서 무제한 증폭되지 않아야 한다. LMM의 동차 오차방정식은

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Crho%28%5Czeta%29%3D0)

의 근으로 설명할 수 있다. 강의의 Root Condition은 다음과 같다.


1. 모든 근

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Czeta_i)


가

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%7C%5Czeta_i%7C%5Cle1)

을 만족한다.

2.

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%7C%5Czeta_i%7C%3D1)

인 모든 근은 **중근이 아닌 단순근**이어야 한다.


왜 중근을 제한할까? 특성다항식의 중근

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Czeta%3D1)


에 대응하는 일반해에는

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20n%5Czeta%5En%3Dn)

과 같은 항이 포함되어 오차가 선형으로 증가할 수 있기 때문이다.

### 2.3 Dahlquist Equivalence and First Barrier

강의 Theorem 5.2(**Dahlquist Equivalence Theorem**)에 따르면, 적절한 가정의 선형 다단계법에서 **consistency와 zero-stability를 함께 만족하는 것**이 수렴의 필요충분조건이 된다. 일관성이 있다고 해서 안정성까지 자동으로 따라오지 않는다는 뜻이다.


강의 Theorem 5.3(**Dahlquist First Barrier**)은 zero-stable

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20k)


-step LMM의 달성 가능한 차수에 제약을 둔다. 일반적 상한은

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20k)


가 홀수이면

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20p%5Cle%20k%2B1)


, 짝수이면

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20p%5Cle%20k%2B2)


이며, explicit인 경우

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20p%5Cle%20k)

이다. 따라서 과거 값을 많이 참조한다고 정확도를 무한히 높일 수 있는 것은 아니다.

#### Connecting Leap-Frog to the Root Condition

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y_%7Bn%2B1%7D%3Dy_%7Bn-1%7D%2B2hf_n)


의

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20h%3D0)


동차식은

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y_%7Bn%2B1%7D-y_%7Bn-1%7D%3D0)

이며

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Crho%28%5Czeta%29%3D%5Czeta%5E2-1%3D%28%5Czeta-1%29%28%5Czeta%2B1%29)


이다. 근

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%2B1%2C-1)


이 모두 단순하고 단위원 위에 있으므로 **zero-stable**이다. 하지만

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y%27%3D-y)


에 대해 유한한

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20h)


에서의

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20r_-%3D-h-%5Csqrt%7B1%2Bh%5E2%7D)


는

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%7Cr_-%7C%3E1)


이므로 긴 시간 동안 parasitic mode가 증폭될 수 있다. 이 두 사실은 모순이 아니다. **Zero-stability는

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20h%5Cto0)

의 고정 구간 분석이고, 장기 고정-step 증폭은 절대안정성 문제**이기 때문이다.

### 3.2 Absolute Stability and Stability Polynomials

강의는 Dahlquist test equation

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y%27%3D%5Clambda%20y%2C%5Cqquad%5CRe%28%5Clambda%29%3C0)


을 사용한다.

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20z%3Dh%5Clambda)

라 두고 LMM에 대입하면

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Csum_%7Bj%3D0%7D%5Ek%28%5Calpha_j-z%5Cbeta_j%29y_%7Bn%2Bj%7D%3D0)

이다. 따라서 안정성 다항식은

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cpi%28%5Czeta%3Bz%29%3D%5Crho%28%5Czeta%29-z%5Csigma%28%5Czeta%29)


이고, 주어진

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20z)

에서 모든 근의 절댓값이 1보다 작으면 해당 증폭 모드가 감쇠한다. 강의의 안정성 영역은 이 근 조건으로 기술한다. **단위원 경계에서의 처리는 zero-stability 조건과 혼동하지 않아야 한다.**

### 3.3 Deriving Three Amplification Factors


**Forward Euler:**

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y_%7Bn%2B1%7D%3D%281%2Bz%29y_n)

이므로

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20R_%7BFE%7D%28z%29%3D1%2Bz)


이며 실수

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20z%3C0)


에서는

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%7C1%2Bz%7C%3C1)

,

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20-2%3Cz%3C0)

가 절대안정 조건이다.


**Backward Euler:**

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y_%7Bn%2B1%7D%3Dy_n%2Bzy_%7Bn%2B1%7D)

이므로

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20R_%7BBE%7D%28z%29%3D%5Cfrac1%7B1-z%7D)


이다.

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5CRe%28z%29%3C0)


이면

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%7C1-z%7C%3E1)

로 감쇠한다.


**Trapezoidal / Crank–Nicolson:**

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y_%7Bn%2B1%7D%3Dy_n%2B%28z%2F2%29%28y_n%2By_%7Bn%2B1%7D%29)

이므로

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20R_%7BTR%7D%28z%29%3D%5Cfrac%7B1%2Bz%2F2%7D%7B1-z%2F2%7D)


이다.

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5CRe%28z%29%3C0)


에서

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%7CR_%7BTR%7D%28z%29%7C%3C1)


이지만

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20z%5Cto-%5Cinfty)


일 때

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20R_%7BTR%7D%28z%29%5Cto-1)

이므로 강한 감쇠가 발생하지는 않는다. 기존 PNG `ode_schemes_stability_comparison.png`는 이 세 함수의 비교 그림이다.

### 4.1 Crank–Nicolson: Full Closed Form and Limit


Lab Handout Exercise 5는

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y%27%3D%5Clambda%20y%2C%5C%20y%280%29%3Dy_0%2C%5C%20h%3DT%2Fn)

이다. 업데이트를 정리하면

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y_%7Bj%2B1%7D%3Dy_j%2B%5Cfrac%7Bh%5Clambda%7D%7B2%7D%28y_j%2By_%7Bj%2B1%7D%29)

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cleft%281-%5Cfrac%7Bh%5Clambda%7D%7B2%7D%5Cright%29y_%7Bj%2B1%7D%3D%5Cleft%281%2B%5Cfrac%7Bh%5Clambda%7D%7B2%7D%5Cright%29y_j)

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20R%28h%5Clambda%29%3D%5Cfrac%7B1%2Bh%5Clambda%2F2%7D%7B1-h%5Clambda%2F2%7D)

이고, 분모가 0이 아닌 경우

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y_n%3Dy_0%5Cleft%28%5Cfrac%7B1%2B%5Clambda%20T%2F%282n%29%7D%7B1-%5Clambda%20T%2F%282n%29%7D%5Cright%29%5En)


이다.

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20a%3D%5Clambda%20T%2F2)

라 놓으면

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Clog%5Cfrac%7B1%2Ba%2Fn%7D%7B1-a%2Fn%7D%3D2%5Cleft%28%5Cfrac%20an%2B%5Cfrac%7Ba%5E3%7D%7B3n%5E3%7D%2BO%28n%5E%7B-5%7D%29%5Cright%29)

이므로

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20n%5Clog%5Cfrac%7B1%2Ba%2Fn%7D%7B1-a%2Fn%7D%3D2a%2BO%28n%5E%7B-2%7D%29%3D%5Clambda%20T%2BO%28n%5E%7B-2%7D%29)


이다. 따라서

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Clim_%7Bn%5Cto%5Cinfty%7Dy_n%3Dy_0e%5E%7B%5Clambda%20T%7D)

를 얻는다. 이것이 **Lab Handout Exercise 5의 수렴성 증명**이다.

### 5.1 Stiff Differential Equations


강의는 서로 다른 속도로 감쇠하는 성분이 공존하여 explicit 방법의 스텝 크기가 **정확도보다 안정성에 의해 심하게 제한되는 문제**를 stiff라고 설명한다. Jacobian

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20J%3D%5Cpartial%20f%2F%5Cpartial%20y)


의 고유값

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Clambda_i)

가 모두 왼쪽 반평면에 있을 때, 감쇠 속도의 비율

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20S%3D%5Cfrac%7B%5Cmax_i%7C%5CRe%28%5Clambda_i%29%7C%7D%7B%5Cmin_i%7C%5CRe%28%5Clambda_i%29%7C%7D%5Cgg1)


은 stiffness의 한 지표다. explicit 적분은 가장 빠른 모드의 제약 때문에

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20h%5Clesssim%20C%2F%5Cmax_i%7C%5Clambda_i%7C)

가 필요할 수 있다. 따라서 실제 관심 대상인 느린 모드가 완만해도 매우 작은 step을 사용해야 한다.

### 5.2 A-Stability and L-Stability


**A-stability**는 안정성 영역이 왼쪽 열린 반평면 전체를 포함한다는 뜻이다. **L-stability**는 A-stability와 더불어

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Clim_%7Bz%5Cto-%5Cinfty%7D%7CR%28z%29%7C%3D0)

을 요구한다. 결론은:

- Forward Euler: A-stable이 아니므로 스텝 크기 제한이 있다.
- Backward Euler(BDF1): A-stable이면서 L-stable이다.

- Trapezoidal: A-stable이지만

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20R%28z%29%5Cto-1)

이므로 L-stable은 아니다.

- BDF2: L-stable 방법의 예로 제시한다. 단, 다단계 방법에서는 특정 모드 증폭을 **단일한 스칼라

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20R%28z%29)

하나로 표현하는 방식에 주의**해야 한다.

**Dahlquist Second Barrier:** 강의는 explicit LMM이 A-stable할 수 없고, A-stable한 implicit LMM은 차수가 최대 2이며, 그 한계를 달성하는 대표적 방법이 Trapezoidal Rule이라고 설명한다.

## 6. Boundary Value Problems — Full Lecture Coverage

여기부터는 기존 글에서 매우 짧게 소개했던 `Lec-4_scheme.pdf` **Section 6 전체**를 상세히 정리한다. BVP는 한 시각의 초기조건 대신 서로 다른 경계 위치의 조건을 지정한다.

### 6.1 Linear Shooting Method

2계 선형 BVP가

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y%27%27%3Dp%28x%29y%27%2Bq%28x%29y%2Br%28x%29%2C%5Cqquad%20y%28a%29%3D%5Calpha%2C%5C%20y%28b%29%3D%5Cbeta)

일 때, 두 개의 IVP를 풀어 선형결합한다.

**IVP 1 — Inhomogeneous:**

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y_1%27%27%3Dp%28x%29y_1%27%2Bq%28x%29y_1%2Br%28x%29%2C%5Cquad%20y_1%28a%29%3D%5Calpha%2C%5Cquad%20y_1%27%28a%29%3D0)

**IVP 2 — Homogeneous:**

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y_2%27%27%3Dp%28x%29y_2%27%2Bq%28x%29y_2%2C%5Cquad%20y_2%28a%29%3D0%2C%5Cquad%20y_2%27%28a%29%3D1)


선형성을 이용하면

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y%28x%29%3Dy_1%28x%29%2Bc%20y_2%28x%29)


이고,

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y%28b%29%3D%5Cbeta)

를 대입하면

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y_1%28b%29%2Bcy_2%28b%29%3D%5Cbeta)

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20c%3D%5Cfrac%7B%5Cbeta-y_1%28b%29%7D%7By_2%28b%29%7D)


을 얻는다(

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y_2%28b%29%5Cne0)

). 미지의 초기 기울기를 경계조건에서 결정해 BVP를 IVP로 푸는 방법이다.

### 6.2 Nonlinear Shooting and the Variational Equation

비선형 BVP는

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y%27%27%3Df%28x%2Cy%2Cy%27%29%2C%5Cqquad%20y%28a%29%3D%5Calpha%2C%5C%20y%28b%29%3D%5Cbeta)


이다. 미지의 초기기울기

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20t%3Dy%27%28a%29)

를 매개변수로 두고

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y%27%27%3Df%28x%2Cy%2Cy%27%29%2C%5Cquad%20y%28a%29%3D%5Calpha%2C%5Cquad%20y%27%28a%29%3Dt)


를 푼다. 경계 잔차

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20F%28t%29%3Dy%28b%3Bt%29-%5Cbeta)

가 0이 되도록 Newton 반복

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20t_%7Bk%2B1%7D%3Dt_k-%5Cfrac%7BF%28t_k%29%7D%7BF%27%28t_k%29%7D)


을 적용한다. 여기서

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20F%27%28t%29%3D%5Cpartial%20y%28b%3Bt%29%2F%5Cpartial%20t)


를 구하려고 민감도 함수

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20z%28x%3Bt%29%3D%5Cpartial%20y%28x%3Bt%29%2F%5Cpartial%20t)


를 도입하면, 양변을

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20t)

에 대해 미분한 결과

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20z%27%27%3Df_y%28x%2Cy%2Cy%27%29z%2Bf_%7By%27%7D%28x%2Cy%2Cy%27%29z%27) 


이고 초기조건은

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20z%28a%29%3D0%2C%5C%20z%27%28a%29%3D1)


이다. 따라서 원래 IVP와 이 **variational equation**을 함께 풀어

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20F%27%28t%29%3Dz%28b%3Bt%29)

를 구할 수 있다.

### 6.3 Finite Difference Method: Mesh and Central Differences


구간

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20a%5Cle%20x%5Cle%20b)


를

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20N)


개의 간격으로 분할하여

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20h%3D%28b-a%29%2FN)


,

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x_i%3Da%2Bih)

로 둔다. 내부 격자점에서는 Taylor 전개로

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y%28x_i%2Bh%29%3Dy_i%2Bhy_i%27%2B%5Cfrac%7Bh%5E2%7D2y_i%27%27%2B%5Cfrac%7Bh%5E3%7D6y_i%5E%7B%283%29%7D%2B%5Cfrac%7Bh%5E4%7D%7B24%7Dy_i%5E%7B%284%29%7D%2BO%28h%5E5%29)

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y%28x_i-h%29%3Dy_i-hy_i%27%2B%5Cfrac%7Bh%5E2%7D2y_i%27%27-%5Cfrac%7Bh%5E3%7D6y_i%5E%7B%283%29%7D%2B%5Cfrac%7Bh%5E4%7D%7B24%7Dy_i%5E%7B%284%29%7D%2BO%28h%5E5%29)

를 얻는다. 두 식을 **빼면**

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y%27%28x_i%29%3D%5Cfrac%7By_%7Bi%2B1%7D-y_%7Bi-1%7D%7D%7B2h%7D%2BO%28h%5E2%29)

이고, 두 식을 **더하면**

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y%27%27%28x_i%29%3D%5Cfrac%7By_%7Bi%2B1%7D-2y_i%2By_%7Bi-1%7D%7D%7Bh%5E2%7D%2BO%28h%5E2%29)

이다. 이 중심차분식이 Lecture 4의 FDM 기본 도구다.

### 6.4 Deriving the Tridiagonal Matrix System


선형 BVP

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y%27%27-p%28x%29y%27-q%28x%29y%3Dr%28x%29)

에 중심차분을 대입하면

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cfrac%7By_%7Bi%2B1%7D-2y_i%2By_%7Bi-1%7D%7D%7Bh%5E2%7D-p_i%5Cfrac%7By_%7Bi%2B1%7D-y_%7Bi-1%7D%7D%7B2h%7D-q_i%20y_i%3Dr_i)


이다. 양변에

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20h%5E2)


를 곱하고

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y_%7Bi-1%7D%2Cy_i%2Cy_%7Bi%2B1%7D)

의 계수를 모으면

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cleft%281%2B%5Cfrac%7Bhp_i%7D2%5Cright%29y_%7Bi-1%7D-%282%2Bh%5E2q_i%29y_i%2B%5Cleft%281-%5Cfrac%7Bhp_i%7D2%5Cright%29y_%7Bi%2B1%7D%3Dh%5E2r_i)


를 얻는다. 여기서 하대각

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cell_i%3D1%2Bhp_i%2F2)


, 대각

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20d_i%3D-%282%2Bh%5E2q_i%29)


, 상대각

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20u_i%3D1-hp_i%2F2)


이며, 경계값

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y_0%3D%5Calpha%2Cy_N%3D%5Cbeta)

는 첫 행과 마지막 행의 우변으로 이동시킨다.


따라서 미지수 벡터

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%28y_1%2C%5Cldots%2Cy_%7BN-1%7D%29%5E%5Ctop)

에 대해

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20Ay%3Db)

인 **삼중대각(tridiagonal) 선형시스템**이 만들어진다. 일반적인 밀집행렬 연산보다 훨씬 효율적으로 풀 수 있다.

### 6.5 Neumann / Robin Conditions and Ghost Points


Dirichlet 조건처럼 함수값이 주어지는 경우와 달리

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y%27%28a%29%3D%5Cgamma)


처럼 도함수 경계조건이 주어지면 가상의 점

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x_%7B-1%7D%3Da-h)

를 도입한다.

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cfrac%7By_1-y_%7B-1%7D%7D%7B2h%7D%3D%5Cgamma)

이므로

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y_%7B-1%7D%3Dy_1-2h%5Cgamma)


이다. 이 값을 경계점의 중심차분식에 대입하면 가상 미지수

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y_%7B-1%7D)


를 제거하면서 2차 정확도를 유지하는 경계방정식을 얻는다. Robin 조건

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20a_0%20y%28a%29%2Bb_0%20y%27%28a%29%3Dc_0)

도 같은 원리로 처리한다.

### 6.6 Thomas Algorithm: Forward Elimination and Back Substitution

삼중대각 시스템을 다음과 같이 적는다.

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cell_i%20y_%7Bi-1%7D%2Bd_i%20y_i%2Bu_i%20y_%7Bi%2B1%7D%3Db_i)

**Forward elimination:** 첫 행을 정규화한다.

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20u%27_1%3D%5Cfrac%7Bu_1%7D%7Bd_1%7D%2C%5Cqquad%20b%27_1%3D%5Cfrac%7Bb_1%7D%7Bd_1%7D)


이후

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20i%3D2%2C%5Cldots%2CN-1)

에 대해

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20u%27_i%3D%5Cfrac%7Bu_i%7D%7Bd_i-%5Cell_i%20u%27_%7Bi-1%7D%7D)

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20b%27_i%3D%5Cfrac%7Bb_i-%5Cell_i%20b%27_%7Bi-1%7D%7D%7Bd_i-%5Cell_i%20u%27_%7Bi-1%7D%7D)


를 계산한다(마지막 행에서는

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20u%27_%7BN-1%7D)

계산을 생략할 수 있음). 각 분모가 0이 아니어야 한다.

**Back substitution:** 마지막 미지수부터

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y_%7BN-1%7D%3Db%27_%7BN-1%7D)

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y_i%3Db%27_i-u%27_i%20y_%7Bi%2B1%7D%2C%5Cqquad%20i%3DN-2%2C%5Cldots%2C1)


로 구한다. 각 행에서 상수 횟수의 계산만 필요하므로 연산 복잡도는

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20O%28N%29)

이다.

```python
import numpy as np

def thomas(lower, diag, upper, rhs):
    """Solve a tridiagonal system (nonzero elimination pivots assumed)."""
    d = np.asarray(diag, dtype=float).copy()
    b = np.asarray(rhs, dtype=float).copy()
    l = np.asarray(lower, dtype=float)
    u = np.asarray(upper, dtype=float)
    n = len(d)
    for i in range(1, n):
        multiplier = l[i-1] / d[i-1]
        d[i] -= multiplier * u[i-1]
        b[i] -= multiplier * b[i-1]
    y = np.empty(n)
    y[-1] = b[-1] / d[-1]
    for i in range(n-2, -1, -1):
        y[i] = (b[i] - u[i] * y[i+1]) / d[i]
    return y

# Independent numerical check: [2,-1,0; -1,2,-1; 0,-1,2] y = [1,0,1]
A = np.array([[2., -1., 0.], [-1., 2., -1.], [0., -1., 2.]])
y = thomas([-1., -1.], [2., 2., 2.], [-1., -1.], [1., 0., 1.])
assert np.allclose(A @ y, [1., 0., 1.])
print('Thomas check passed:', y)
```

이 코드는 Lecture 4 Section 6.3의 알고리즘을 수치적으로 확인하기 위한 **추가 구현**이다. 

### Complete Lecture Coverage Check

`Lec-4_scheme.pdf`의 **Section 4: Linear Multistep, Section 5: Stability/Stiffness, Section 6: BVP/Shooting/FDM/Thomas**를 모두 이 포스트에서 다뤘다. Lecture 2–3의 **Adams 도입, 오차비교 표, 연립계 변환**은 Post 3과 교차 연결된다. 기존 PNG 파일 경로는 삭제하거나 이름을 바꾸지 않았다.

## Source-by-Source Verification — Stability and Boundary Conditions

### Characteristic Polynomials: Explicit Consistency Check

`k`-step 방법의 특성다항식을

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Crho%28%5Czeta%29%3D%5Csum_%7Bj%3D0%7D%5E%7Bk%7D%5Calpha_j%5Czeta%5Ej%2C%5Cquad%5Csigma%28%5Czeta%29%3D%5Csum_%7Bj%3D0%7D%5E%7Bk%7D%5Cbeta_j%5Czeta%5Ej)

로 쓰면 강의 p.11의 일관성 조건은

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Crho%281%29%3D0%2C%5Cqquad%5Crho%5Cprime%281%29%3D%5Csigma%281%29)

이다. Forward Euler의 경우 `ρ(ζ)=ζ−1`, `σ(ζ)=1`이므로 각각 `0`과 `1=1`로 검증된다. Trapezoidal Rule에서는 `ρ(ζ)=ζ−1`, `σ(ζ)=(ζ+1)/2`이므로 역시 `0`, `1=1`이 된다.

### BDF2 and the Test Equation

강의 p.12의 BDF2는

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y_%7Bn%2B2%7D-%5Cfrac43y_%7Bn%2B1%7D%2B%5Cfrac13y_n%3D%5Cfrac%7B2h%7D%7B3%7Df_%7Bn%2B2%7D)

이므로 테스트 방정식 `y′=λy`, `z=hλ`를 대입하면

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%281-%5Cfrac23z%29%5Czeta%5E2-%5Cfrac43%5Czeta%2B%5Cfrac13%3D0)

을 얻는다. 다단계 방법의 안정성은 **이 근들을 동시에** 살펴야 한다. `R(z)` 하나로 판정하는 일단계 방법과 구분해야 한다.

### Trapezoidal Rule: Amplification, A-Stability and Its Limitation

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20R%28z%29%3D%5Cfrac%7B1%2Bz%2F2%7D%7B1-z%2F2%7D)

`Re(z)<0`이면 `|1+z/2| < |1−z/2|`이므로 `|R(z)|<1`이고, 따라서 A-stable이다. 그러나 음의 실수축 `z→−∞`에서 `R(z)→−1`이라 감쇠가 강한 모드가 진동하며 사라지지 않는다.

### Crank–Nicolson Lab Exercise 5: Closed Form and Limit

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y_n%3Dy_0%5Cleft%28%5Cfrac%7B1%2B%5Clambda%20T%2F%282n%29%7D%7B1-%5Clambda%20T%2F%282n%29%7D%5Cright%29%5En)

분모가 0이 되는 `n`에서는 해당 암시적 한 스텝을 풀 수 없으므로 이 식은 유한한 값으로 정의되지 않는다. 충분히 큰 `n`에서는 분모가 0이 아니고,

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Clog%5Cfrac%7B1%2Bw%2F2%7D%7B1-w%2F2%7D%3Dw%2B%5Cfrac%7Bw%5E3%7D%7B12%7D%2BO%28w%5E5%29%2C%5Cquad%20w%3D%5Clambda%20T%2Fn)

따라서 `n log R(λT/n) → λT`이고 `yₙ→y₀e^(λT)`이다. 여기서 로그의 전개는 `n`이 충분히 커서 인수가 양수인 실수 구간에서 적용한다.

### BVP Matrix: Boundary-Value Contribution Check

`y''−p(x)y'−q(x)y=r(x)`의 중앙차분을 쓰면 첫 내부 격자점 `i=1`의 계수는

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cell_1y_0%2Bd_1y_1%2Bu_1y_2%3Dh%5E2r_1)

이다. `y₀=α`를 알려진 값으로 옮기면 **첫 우변은 `h²r₁−ℓ₁α`**이다. 마지막 내부 점에서는 `y_N=β`를 옮기므로 **마지막 우변은 `h²r_(N−1)−u_(N−1)β`**다. 이 경계항의 부호를 확인하지 않으면 중앙차분식 자체는 맞더라도 전체 선형계가 잘못된다.

### Lecture 4: Scope of Proofs

강의 p.13–14의 Dahlquist Equivalence, First Barrier, Second Barrier는 **정리의 진술이 제시된 부분**이며 본문에 일반적인 완전 증명이 제공되지는 않는다. 이 블로그의 특정 방법에 대한 특성다항식 대입과 root check는 **추가 검산**이고, 그 자체로 Dahlquist 장벽정리 전체를 증명하지는 않는다.


다단계법의 정확도를 평가할 때는 국소 절단오차와 **Zero-Stability와 Absolute Stability를 구분**해야 한다. Crank–Nicolson은 2차 정확도와 A-Stability를 얻지만, 매우 빠른 감쇠 모드를 강하게 없애는 L-Stability는 갖지 않는다.

## References

- Raj Dhara, *Lec-2_3_scheme*, 2026-10-12, pp. 10–11.
- Raj Dhara, *Lec-4_scheme*, pp. 10–14.
- Raj Dhara, *Numerical Method Lab: ODE Numerical Schemes & Convergence*, Exercise 5.
- Raj Dhara, *Numerical Method Lab 4: Group-2 — ODE Numerical Schemes & Convergence*, 2026-10-28, Exercises 1–3.
- `LAB_4.ipynb`, Exercises 1(a)–3(c): numerical experiments and generated figures.
---

