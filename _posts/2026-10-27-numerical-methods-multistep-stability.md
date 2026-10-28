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

이번 글에서는 이전 값을 여러 개 사용하는 Linear Multistep Methods를 정리하고, Consistency와 Zero-Stability가 수렴에 함께 필요하다는 것을 확인한다. `Lec-2_3_scheme.pdf`의 10–11쪽과 `Lec-4_scheme.pdf`의 10–14쪽을 연결하고, Lab Handout Exercise 5의 Crank–Nicolson 스킴을 안정성 함수로 해석한다.

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

### 9.5 Lab Handout Exercise 5 — Crank–Nicolson / Trapezoidal Scheme

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

`Lec-4_scheme.pdf`의 15–16쪽에서는 Boundary Value Problem, Shooting Method, Finite Difference Method 및 삼대각행렬을 푸는 Thomas Algorithm을 소개한다. 이 내용은 현재 Lab 1의 IVP 문제와는 별도 주제이므로 다음 글에서 구체적으로 다룬다.


### Lab 2 Connection — Leap-Frog Stability

`LAB_2.ipynb`의 Exercise 3에서는 Leap-Frog의 국소절단오차와 특성근을 직접 확인한다. `T = 1`에서는 `h`를 작게 할수록 전역오차가 줄어드는 결과가 관찰된다. 그러나 `h`를 고정하고 장시간 적분하면 `|−h − sqrt(1+h²)| > 1`인 모드가 증폭된다. 이는 **고정된 구간에서 h → 0인 수렴 문제**와 **고정된 h에서 장시간 적분할 때의 절대 안정성 문제**가 다르다는 점을 강조한다. 실제 실행 코드와 세 그래프는 앞선 Euler Schemes 글의 Lab 2 Exercise 3에 배치했다.

### Lab 3 Connection — Runge-Kutta for ODE Systems

`LAB_3.ipynb`에서는 벡터 RK4로 강제진동계 및 결합된 스프링-질량계를 계산한다. 이들은 Linear Multistep Method를 구현한 실습이 아니라 **one-step Runge–Kutta를 연립 ODE에 적용하는 실습**이므로 자세한 구현과 그래프는 앞선 Taylor and Runge-Kutta Methods 글에 배치했다. 이 글에서는 반복적으로 사용되는 안정성 및 수렴차수 개념과의 연결만 정리한다.



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

즉,

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


, 즉

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

을 요구한다. 강의의 결론은:

- Forward Euler: A-stable이 아니므로 스텝 크기 제한이 있다.
- Backward Euler(BDF1): A-stable이면서 L-stable이다.

- Trapezoidal: A-stable이지만

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20R%28z%29%5Cto-1)

이므로 L-stable은 아니다.

- BDF2: 강의에서 L-stable 방법의 예로 제시한다. 단, 다단계 방법에서는 특정 모드 증폭을 **단일한 스칼라

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

). 즉, 미지의 초기 기울기를 경계조건에서 결정해 BVP를 IVP로 푸는 방법이다.

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
A = np.array(![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%202.%2C%20-1.%2C%200.%5D%2C%20%5B-1.%2C%202.%2C%20-1.%5D%2C%20%5B0.%2C%20-1.%2C%202.))
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

`Re(z)<0`이면 `|1+z/2| < |1−z/2|`이므로 `|R(z)|<1`이고, 따라서 A-stable이다. 그러나 음의 실수축 `z→−∞`에서 `R(z)→−1`이라 감쇠가 강한 모드가 진동하며 사라지지 않는다. 즉 **A-stable이지만 L-stable은 아니다.**

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


다단계법의 정확도를 평가할 때는 국소 절단오차뿐 아니라 **Zero-Stability와 Absolute Stability를 구분**해야 한다. Crank–Nicolson은 2차 정확도와 A-Stability를 얻지만, 매우 빠른 감쇠 모드를 강하게 없애는 L-Stability는 갖지 않는다.

## References

- Raj Dhara, *Lec-2_3_scheme*, 2026-10-12, pp. 10–11.
- Raj Dhara, *Lec-4_scheme*, pp. 10–14.
- Raj Dhara, *Numerical Method Lab: ODE Numerical Schemes & Convergence*, Exercise 5.


---

