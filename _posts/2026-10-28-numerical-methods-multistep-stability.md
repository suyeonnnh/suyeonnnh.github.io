---
title: "Numerical Methods : Multistep Methods and Stability"
date: 2026-10-28 18:30:00 +0100
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

## Summary and Takeaways

다단계법의 정확도를 평가할 때는 국소 절단오차뿐 아니라 **Zero-Stability와 Absolute Stability를 구분**해야 한다. Crank–Nicolson은 2차 정확도와 A-Stability를 얻지만, 매우 빠른 감쇠 모드를 강하게 없애는 L-Stability는 갖지 않는다.

## References

- Raj Dhara, *Lec-2_3_scheme*, 2026-10-12, pp. 10–11.
- Raj Dhara, *Lec-4_scheme*, pp. 10–14.
- Raj Dhara, *Numerical Method Lab: ODE Numerical Schemes & Convergence*, Exercise 5.
