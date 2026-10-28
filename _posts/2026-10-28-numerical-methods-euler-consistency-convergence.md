---
title: "Numerical Methods : Euler Schemes, Consistency and Convergence"
date: 2026-10-28 18:10:00 +0100
categories:
  - Numerical Analysis
tags:
  - Euler-Method
  - Consistency
  - Convergence
  - Leap-Frog
  - Python-Lab
math: false
---

## Introduction

`Lec-2_3_scheme.pdf`의 1–4쪽을 중심으로 Forward/Backward Euler, Local Truncation Error, Consistency, Global Convergence 및 Leap-Frog의 불안정성을 공부한다. `LAB_1.ipynb`의 Python Exercise 1–3은 설명 직후에 배치했고, Lab Handout의 Forward/Backward Euler 문제도 함께 풀이했다.

## 8. Numerical Solutions of ODEs

실제 ODE는 정확한 해를 초등함수 형태로 나타내기 어려운 경우가 많다. 따라서 `t₀, t₁, …, tₙ`의 격자점에서 해를 근사한다.

### 8.1 Forward Euler Method

Forward Euler는 현재 위치에서의 기울기를 이용해 다음 지점을 예측하는 가장 기본적인 명시적 수치 스킴이다.

![Forward Euler update](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7Dy_%7Bn%2B1%7D%3Dy_n%2Bh%5C%2Cf%28t_n%2Cy_n%29)

계산 과정은 다음과 같다.

1. 초기값 `y₀`과 시간 간격 `h`를 지정한다.
2. 현재 위치에서 기울기 `f(tₙ,yₙ)`를 계산한다.
3. `yₙ₊₁ = yₙ + h·f(tₙ,yₙ)`로 한 단계 진행한다.
4. 목표 시간까지 반복한다.

기울기를 현재 시점에서 한 번만 계산하므로 구현이 간단하다. 반면 오차가 비교적 크고 시간 간격에 따라 안정성 문제가 발생할 수 있다.

### 8.2 Lecture Example — Two Euler Steps

강의의 예제는 `y'=t−y`, `y(0)=1`, `h=0.2`이다.

- 첫 단계: `f(0,1)=−1`이므로 `y₁=1+0.2(−1)=0.8`.
- 두 번째 단계: `f(0.2,0.8)=−0.6`이므로 `y₂=0.8+0.2(−0.6)=0.68`.

이처럼 **미분방정식의 우변을 매번 현재 근삿값에 대입한다**는 점이 Euler의 핵심이다.

### 8.3 Local Truncation Error, Consistency, and Global Error

여기서 세 개념을 구분해야 한다.

- **Local Truncation Error (LTE):** 이전 단계에 정확한 값을 넣었다고 가정할 때 한 단계의 수치식이 만족하지 못하는 차이.
- **Consistency:** 시간 간격 `h → 0`일 때 적절히 정규화한 LTE가 0으로 수렴하는 성질.
- **Global Error:** 시작점부터 반복한 근삿값 `yₙ`과 정확한 값 `y(tₙ)` 사이의 차이.

![Forward Euler normalized LTE](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Ctau_%7Bn%2B1%7D%3D%5Cfrac%7By%28t_%7Bn%2B1%7D%29-y%28t_n%29%7D%7Bh%7D-f%28t_n%2Cy%28t_n%29%29)

![Global error](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7De_n%3D%5Cleft%7Cy%28t_n%29-y_n%5Cright%7C)

여기서는 **`h`로 나눈 정규화 LTE**를 사용한다. 이 정의에 따르면 Forward Euler의 LTE는 보통 `O(h)`이고, 한 단계 업데이트 자체의 결함은 `O(h²)`이다. 두 표기가 혼용되므로 정의를 확인해야 한다.

### 8.4 Notebook Experiment 1 — Forward Euler LTE and Consistency

업로드한 `LAB_1.ipynb`의 첫 번째 실험에서는 `y'=−y`, 정확한 해 `y(t)=e^(−t)`를 사용한다. 고정된 지점 `t=1`에서 `h`를 변화시키며 정규화 LTE를 계산한다. 이 실험은 *Lab-1.pdf*의 Exercise 1과는 별개의 실험이다.

![Euler LTE for exponential decay](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Ctau%3D%5Cfrac%7Be%5E%7B-1%7D-e%5E%7B-%281-h%29%7D%7D%7Bh%7D%2Be%5E%7B-%281-h%29%7D)

Taylor 전개를 적용하면 LTE의 선도항 크기는 `e^(−1)·h/2`이므로 `h`가 충분히 작아지면 로그-로그 그래프에서 기울기 약 1을 기대한다.

**Python Implementation**

```python
import numpy as np
import matplotlib.pyplot as plt

# Exercise 1: Forward Euler LTE and Consistency

# Fixed evaluation point
x_n = 1.0

# Number of grid intervals
n_values = [10, 100, 1000, 10000]

# Lists to store results
h_values = []
lte_values = []

# Calculate LTE for different step sizes
for n in n_values:
    h = 1 / n
    x_prev = x_n - h

    # Exact solution y(x) = exp(-x)
    y_n = np.exp(-x_n)
    y_prev = np.exp(-x_prev)

    # Forward Euler LTE
    lte = (y_n - y_prev) / h + y_prev

    h_values.append(h)
    lte_values.append(abs(lte))

    print(f"n = {n:5d}, h = {h:.5f}, LTE = {lte:.10f}")

# Plot the results
plt.figure(figsize=(7, 5))
plt.loglog(h_values, lte_values, 'o-', label='LTE')
plt.loglog(h_values, np.exp(-x_n) * np.array(h_values) / 2,
           '--', label='Reference: exp(-1)*h/2')

plt.xlabel('Step size (h)')
plt.ylabel('Absolute LTE')
plt.title('Forward Euler: LTE and Consistency')
plt.legend()
plt.grid(True, which='both')
plt.tight_layout()
plt.show()
```

**실행 결과**

```text
n =    10, h = 0.10000, LTE = 0.0196674740
n =   100, h = 0.01000, LTE = 0.0018517060
n =  1000, h = 0.00100, LTE = 0.0001840624
n = 10000, h = 0.00010, LTE = 0.0000183952
```

![Forward Euler LTE](/assets/images/ode_lab1_euler_lte.png)

*Figure 1. Forward Euler의 정규화 LTE와 이론적 기준선. 시간 간격이 작아질수록 LTE가 0에 가까워진다.*

**Results and Analysis.** `h=0.1`일 때 LTE가 약 `0.019667`이었으나 `h=0.0001`에서 약 `0.0000184`로 감소했다. 수치 LTE가 기준선 `e^(−1)·h/2`에 점차 가까워지므로 이 예제에서 Forward Euler의 **Consistency**를 확인할 수 있다.

### 8.5 Notebook Experiment 2 — Forward Euler Global Convergence

두 번째 실험은 `y'=−2y`, `y(0)=1`이며 `T=1`에서 정확한 해 `e^(−2)`와 수치해를 비교한다.

![Euler convergence update](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7Dy_%7Bn%2B1%7D%3D%281-2h%29y_n%2C%5Cqquad%20y_n%3D%281-2%2Fn%29%5En)

Global Error를 여러 `h`에 대해 계산한 다음, `ln(Error)`와 `ln(h)`의 선형 회귀 기울기 `p`로 수렴 차수를 추정한다.

![Convergence order regression](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Clog%28E%28h%29%29%5Capprox%20p%5Clog%28h%29%2B%5Clog%28C%29)

기울기가 약 `1`이면 오차가 `h`에 비례한다는 뜻이다.

**Python Implementation**

```python
import numpy as np
import matplotlib.pyplot as plt

# Exercise 2: Forward Euler Convergence

# Initial condition and final time
y0 = 1.0
T = 1.0

# Exact solution: y(x) = exp(-2x)
exact_y = np.exp(-2 * T)

# Number of grid intervals
n_values = [10, 20, 40, 80, 160]

# Lists to store results
h_values = []
error_values = []

# Solve the IVP using Forward Euler
for n in n_values:
    h = T / n
    y = y0

    for i in range(n):
        y = y - 2 * h * y

    # Global error at T = 1
    error = abs(exact_y - y)

    h_values.append(h)
    error_values.append(error)

    print(f"n = {n:3d}, h = {h:.5f}, "
          f"y_n = {y:.10f}, Error = {error:.10f}")

# Estimate the empirical order of convergence
log_h = np.log(h_values)
log_error = np.log(error_values)

p, intercept = np.polyfit(log_h, log_error, 1)

print(f"\nEstimated convergence order p = {p:.6f}")

# Plot ln(Error) versus ln(h)
plt.figure(figsize=(7, 5))
plt.plot(log_h, log_error, 'o-', label='Global Error')
plt.plot(log_h, p * log_h + intercept,
         '--', label=f'Fitted line (p = {p:.4f})')

plt.xlabel('ln(h)')
plt.ylabel('ln(Error)')
plt.title('Forward Euler: Global Convergence')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
```

**실행 결과**

```text
n =  10, h = 0.10000, y_n = 0.1073741824, Error = 0.0279611008
n =  20, h = 0.05000, y_n = 0.1215766546, Error = 0.0137586286
n =  40, h = 0.02500, y_n = 0.1285121566, Error = 0.0068231267
n =  80, h = 0.01250, y_n = 0.1319378054, Error = 0.0033974778
n = 160, h = 0.00625, y_n = 0.1336400679, Error = 0.0016952153

Estimated convergence order p = 1.010556
```

![Forward Euler convergence](/assets/images/ode_lab1_euler_convergence.png)

*Figure 2. Forward Euler의 Global Error를 이용한 경험적 수렴차수 추정.*

**Results and Analysis.** 시간 간격을 절반으로 줄이면 Global Error도 대략 절반으로 감소했다. 선형 회귀로 구한 수렴차수는 `p ≈ 1.010556`이므로 이 문제에서 Forward Euler의 **1차 수렴(First-order Convergence)**을 확인할 수 있다.

### 8.6 Stability and the Difference Between Consistency and Convergence

Consistency는 한 단계의 수치식이 미분방정식을 잘 근사하는지를 확인하는 조건이다. 반면 Stability는 수치 연산 중 생긴 작은 오차가 반복 계산을 거치며 어떻게 전파되는지를 다룬다.

특히 다단계 방법에서는 정확한 미분방정식에서 유래하지 않은 **parasitic mode**가 추가로 나타날 수 있다. 따라서 Consistency만 확인해서는 전체 알고리즘의 수렴성을 일반적으로 보장할 수 없다.

### 8.7 Notebook Experiment 3 — Leap-Frog and Long-Time Instability

세 번째 실험은 다음 ODE와 다단계 스킴을 사용한다.

![Leap frog recurrence](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7Dy%27%3D-y%2C%5Cquad%20y%280%29%3D1%2C%5Cqquad%20y_%7Bn%2B1%7D%3Dy_%7Bn-1%7D-2hy_n)

Notebook에서는 두 초기 데이터를 `y₀=1`, `y₁=1`로 지정했다. 하지만 정확한 해의 첫 번째 다음 값은 `y(h)=e^(−h)`이므로 **시작값 `y₁=1`에는 `O(h)`의 시작 오차가 있다.** 이 때문에 고정 시간 `T=1`에서 관찰하는 수렴차수가 스킴 자체의 형식적 정확도와 다를 수 있다.

특성방정식을 세우면 다음과 같다.

![Leap frog characteristic roots](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Crho%5E2%2B2h%5Crho-1%3D0%2C%5Cqquad%5Crho_%7B%5Cpm%7D%3D-h%5Cpm%5Csqrt%7B1%2Bh%5E2%7D)

한 근은 약 `1−h`로 감소하는 해를 따라가지만, 다른 근은 `−1−h`에 가까워 절댓값이 1보다 크다. **후자의 근이 불안정한 기생 모드를 생성한다.** 고정된 `T=1`의 수치 오차가 감소하는 현상과 긴 시간에 진동이 증폭되는 현상은 모순이 아니다.

**Python Implementation**

```python
import numpy as np
import matplotlib.pyplot as plt

# Exercise 3: Consistency Does Not Imply Convergence
# Method: Leap-Frog
# ODE: y' = -y, y(0) = 1
# Scheme: y[n+1] = y[n-1] - 2*h*y[n]

# Initial values and final time
y0 = 1.0
y1 = 1.0
T = 1.0

# Different grid sizes
n_values = [10, 20, 40, 80, 160]

# Exact solution at T = 1
exact_y = np.exp(-T)

h_values = []
numerical_values = []
error_values = []

# Apply the Leap-Frog scheme
for n in n_values:
    h = T / n

    y = np.zeros(n + 1)
    y[0] = y0
    y[1] = y1

    for i in range(1, n):
        y[i + 1] = y[i - 1] - 2 * h * y[i] 

    # Numerical solution and global error at T = 1
    numerical_y = y[-1]
    error = abs(exact_y - numerical_y)

    h_values.append(h)
    numerical_values.append(numerical_y)
    error_values.append(error)

    print(f"n = {n:3d}, h = {h:.5f}, "
          f"y(T) = {numerical_y:.10f}, "
          f"Error = {error:.10f}")


# Plot 1: Numerical vs Exact Solution at T = 1
plt.figure(figsize=(9, 5))

plt.plot(n_values, numerical_values, 'o-',
         label='Leap-Frog numerical solution')

plt.plot(n_values, [exact_y] * len(n_values), 'o-',
         label='Exact solution')

plt.xlabel('Grid Size (n)')
plt.ylabel('Solution at T = 1')
plt.title('Leap-Frog Method: Numerical vs Exact Solution')
plt.xticks(n_values)
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()


# Plot 2: Global Error vs Step Size
plt.figure(figsize=(8, 5))

plt.loglog(h_values, error_values, 'o-',
           label='Global Error')

# Decreasing step size from left to right
plt.gca().invert_xaxis()

plt.xlabel('Step Size (h)')
plt.ylabel('Global Error at T = 1')
plt.title('Leap-Frog Method: Global Error')
plt.legend()
plt.grid(True, which='both', alpha=0.3)
plt.tight_layout()
plt.show()


# Plot 3: Long-Time Instability (Additional Analysis)
# Fixed step size h = 0.1 and extended time T = 20

h_long = 0.1
T_long = 20.0
n_long = int(T_long / h_long)

x_long = np.linspace(0, T_long, n_long + 1)
y_long = np.zeros(n_long + 1)

y_long[0] = 1.0
y_long[1] = 1.0

for i in range(1, n_long):
    y_long[i + 1] = y_long[i - 1] - 2 * h_long * y_long[i]

plt.figure(figsize=(9, 5))

plt.plot(x_long, y_long,
         label='Leap-Frog (h = 0.1)')

plt.plot(x_long, np.exp(-x_long), 'k--',
         linewidth=2, label='Exact Solution')

plt.xlabel('Time (t)')
plt.ylabel('y(t)')
plt.title('Leap-Frog Method: Long-Time Instability')
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()
```

**실행 결과**

```text
n =  10, h = 0.10000, y(T) = 0.2576276480, Error = 0.1102517932
n =  20, h = 0.05000, y(T) = 0.3108433549, Error = 0.0570360863
n =  40, h = 0.02500, y(T) = 0.3389179443, Error = 0.0289614969
n =  80, h = 0.01250, y(T) = 0.3532924543, Error = 0.0145869868
n = 160, h = 0.00625, y(T) = 0.3605599872, Error = 0.0073194539
```

![Leap-Frog short horizon](/assets/images/ode_lab1_leapfrog_short.png)

*Figure 3. T=1에서 격자 수에 따른 Leap-Frog 근삿값과 정확한 해의 비교.*

![Leap-Frog global error](/assets/images/ode_lab1_leapfrog_error.png)

*Figure 4. T=1에서 시간 간격 감소에 따른 Global Error.*

![Leap-Frog long-time instability](/assets/images/ode_lab1_leapfrog_long.png)

*Figure 5. 고정된 h=0.1, T=20에서 Leap-Frog의 불안정한 진동 모드가 성장하는 모습.*

**Results and Analysis.** `T=1`에서는 `n`이 커질수록 오차가 감소했다. 그러나 `h=0.1`로 고정하고 `T=20`까지 계산하면 정확한 해 `e^(−t)`는 0으로 감소하는 반면, 수치해에는 커지는 교대 진동이 나타난다. 따라서 이 실험은 **짧은 구간에서의 오차 감소가 장기 안정성을 의미하지는 않는다**는 점을 보여준다.

또한 Leap-Frog의 매끄러운 정확해에 대한 한 단계 잔차는 2차로 작아지더라도, 두 단계 시작값의 정확도와 특성근에 대한 안정성 조건을 별도로 검토해야 한다.

## Lab Handout Problems — Euler Schemes

### 9.1 Lab Handout Exercise 1 — Forward Euler with Negative Growth Rate

문제는 `y'=−4y`, `y(0)=3`이고 최종시간은 `T=1`, 간격은 `h=1/n`이다.

Forward Euler를 적용하면 다음과 같다.

![Forward Euler negative growth](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7Dy_%7Bj%2B1%7D%3Dy_j-4hy_j%3D%281-4h%29y_j)

![Closed form Euler negative growth](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7Dy_n%3D3%5Cleft%281-%5Cfrac4n%5Cright%29%5En)

정확한 해는 `y(t)=3e^(−4t)`이다. 지수함수의 극한 공식에서 `(1−4/n)^n → e^(−4)`이므로 `n → ∞`에서 수치해가 `y(1)=3e^(−4)`로 수렴한다. 다만 작은 `n`에서는 증폭인자 `1−4h`가 음수이거나 절댓값이 클 수 있으므로 수렴 극한과 개별 시간 간격의 안정성은 구분해야 한다.

### 9.3 Lab Handout Exercise 3 — Backward / Implicit Euler

문제는 `y'=5y`, `y(0)=2`, `T=1`이다. Backward Euler는 다음 시간의 미지값 `yₙ₊₁`에서 변화율을 평가한다.

![Backward Euler equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7Dy_%7Bj%2B1%7D%3Dy_j%2Bh%5C%2C5y_%7Bj%2B1%7D)

![Backward Euler closed form](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7Dy_n%3D2%5Cleft%281-%5Cfrac5n%5Cright%29%5E%7B-n%7D)

업데이트식이 `yₙ₊₁ = yₙ/(1−5h)`로 정리되므로 단계를 반복하면 위의 식을 얻는다. 정확한 해는 `2e^(5t)`이며 `n → ∞`에서 `yₙ → 2e^5`이다. **단, `n=5`이면 분모가 0이 되어 스킴이 정의되지 않는다.** 극한을 다룰 때는 충분히 큰 `n`을 고려한다.

### Important Note — Zero-Stability vs Long-Time Instability

Leap-Frog의 근 `r₋=−h−√(1+h²)`는 고정된 `h>0`에서 절댓값이 1보다 커 장기적으로 기생 모드가 증폭될 수 있다. **그러나 Zero-Stability를 판정하는 다항식은 h→0의 `ρ(ξ)=ξ²−1`이며 근 ±1이 단순근이므로 Zero-Stability 자체는 만족한다.** 따라서 유한한 시간 구간의 h→0 수렴과 고정 h에서의 장기 불안정을 혼동해서는 안 된다. 상세한 Root Condition은 Post 4에서 다룬다.


## Further Validation — Lab 2

Lab 1에서 확인한 세 가지 개념을 Lab 2에서 다시 검증한다. 같은 초기값 문제를 사용하지만, 두 번째 실습은 국소절단오차의 절댓값, 로그 척도에서의 수렴차수, Leap-Frog의 두 특성근과 불안정 모드의 증폭을 더 분명히 드러낸다. 아래 출력 수치는 제출한 `LAB_2.ipynb`에 저장된 결과이다.

### Lab 2 — Exercise 1: LTE and Consistency

*Source: Lab-2.pdf, Exercise 1; Lec-2_3_scheme.pdf, Euler Method*

`y′ = −y`, `y(0) = 1`에서 `y(x) = exp(−x)`를 이용한다. 고정된 `x = 1`에서 `n = 10, 100, 1000, 10000`으로 그리드를 세분화하면 LTE의 절댓값이 감소한다. 1차 일치성이 의미하는 것은 `h → 0`에서 미분방정식의 잔차가 사라진다는 사실이지, 이것만으로 모든 수치해가 수렴한다는 뜻은 아니다.

**Python Implementation**

```python
import numpy as np
import matplotlib.pyplot as plt

# Exercise 1: Forward Euler LTE and Consistency

# Fixed evaluation point
x_n = 1.0
n_values = [10, 100, 1000, 10000]

h_values = []
lte_values = []

# Calculate LTE using the exact solution y(x) = exp(-x)
for n in n_values:
    h = 1 / n
    x_prev = x_n - h
    y_n = np.exp(-x_n)
    y_prev = np.exp(-x_prev)

    # Residual LTE for Forward Euler applied to y' = -y
    lte = (y_n - y_prev) / h + y_prev
    h_values.append(h)
    lte_values.append(abs(lte))
    print(f"n = {n:5d}, h = {h:.5f}, |LTE| = {abs(lte):.10f}")

plt.figure(figsize=(7, 5))
plt.loglog(h_values, lte_values, 'o-', label='Absolute LTE')
plt.loglog(h_values, np.exp(-x_n) * np.array(h_values) / 2,
           '--', label='Reference: exp(-1) h / 2')
plt.xlabel('Step size (h)')
plt.ylabel('Absolute LTE')
plt.title('Exercise 1: Forward Euler Consistency')
plt.legend()
plt.grid(True, which='both')
plt.tight_layout()
plt.show()
```

**실행 결과**

```text
n =    10, h = 0.10000, |LTE| = 0.0196674740
n =   100, h = 0.01000, |LTE| = 0.0018517060
n =  1000, h = 0.00100, |LTE| = 0.0001840624
n = 10000, h = 0.00010, |LTE| = 0.0000183952
```

![Lab 2 — Exercise 1: LTE and Consistency — Figure 1](/assets/images/ode_lab2_ex1_1.png)

*Figure. Lab 2 — Exercise 1: LTE and Consistency (1/1).*

### Lab 2 — Exercise 2: Empirical Convergence Order

*Source: Lab-2.pdf, Exercise 2; Lec-2_3_scheme.pdf, Convergence*

`y′ = −2y`, `y(0) = 1`의 정해는 `exp(−2t)`이다. `T = 1`에서 전역오차를 측정하고 `ln(Error)`와 `ln(h)`의 선형 회귀 기울기를 추정한다. 저장된 결과는 `p = 1.0106`으로, Forward Euler의 이론적 1차 전역 수렴과 일치한다.

**Python Implementation**

```python
import numpy as np
import matplotlib.pyplot as plt

# Exercise 2: Forward Euler Global Convergence and Order

T = 1.0
n_values = [10, 20, 40, 80, 160]

# Exact solution to y' = -2y, y(0) = 1
def exact_solution(x):
    return np.exp(-2 * x)

# Forward Euler scheme
def forward_euler(y0, T, n):
    h = T / n
    x = np.linspace(0, T, n + 1)
    y = np.zeros(n + 1)
    y[0] = y0
    for i in range(n):
        y[i + 1] = y[i] + h * (-2 * y[i])
    return x, y, h

h_values = []
error_values = []

for n in n_values:
    x, y, h = forward_euler(1.0, T, n)
    error = abs(exact_solution(T) - y[-1])
    h_values.append(h)
    error_values.append(error)
    print(f"n = {n:5d}, h = {h:.5f}, y(T) = {y[-1]:.8f}, Error = {error:.10f}")

# Estimate order from the slope of log(error) vs log(h)
p = np.polyfit(np.log(h_values), np.log(error_values), 1)[0]
print(f"Estimated convergence order p = {p:.4f}")

plt.figure(figsize=(7, 5))
plt.loglog(h_values, error_values, 'o-', label='Global error at T = 1')
plt.loglog(h_values, error_values[0] * np.array(h_values) / h_values[0],
           '--', label='Reference: O(h)')
plt.xlabel('Step size (h)')
plt.ylabel('Absolute global error')
plt.title('Exercise 2: Forward Euler Global Convergence')
plt.legend()
plt.grid(True, which='both')
plt.tight_layout()
plt.show()
```

**실행 결과**

```text
n =    10, h = 0.10000, y(T) = 0.10737418, Error = 0.0279611008
n =    20, h = 0.05000, y(T) = 0.12157665, Error = 0.0137586286
n =    40, h = 0.02500, y(T) = 0.12851216, Error = 0.0068231267
n =    80, h = 0.01250, y(T) = 0.13193781, Error = 0.0033974778
n =   160, h = 0.00625, y(T) = 0.13364007, Error = 0.0016952153
Estimated convergence order p = 1.0106
```

![Lab 2 — Exercise 2: Empirical Convergence Order — Figure 1](/assets/images/ode_lab2_ex2_1.png)

*Figure. Lab 2 — Exercise 2: Empirical Convergence Order (1/1).*

### Lab 2 — Exercise 3: Leap-Frog and the Parasitic Mode

*Source: Lab-2.pdf, Exercise 3; Lec-2_3_scheme.pdf, Leap-Frog*

`y[n+1] = y[n−1] − 2h y[n]`, `y[0] = y[1] = 1`에서 특성근은 `−h ± sqrt(1+h²)`이다. 고정된 짧은 종점 `T = 1`에서 `h`를 줄이면 출력된 전역오차가 작아지지만, 고정된 `h`로 오래 적분하면 절댓값 1보다 큰 음의 근에 대응하는 *parasitic mode*가 증폭된다. 따라서 이 장기 불안정성과 유한 시간 구간에서의 수렴을 구분해야 한다.

**Python Implementation**

```python

import numpy as np
import matplotlib.pyplot as plt

# Exercise 3: Consistency Does Not Imply Convergence (The Leap-Frog Trap)

# Final time
T = 1.0

# Number of grid intervals
n_values = [10, 20, 40, 80, 160]

# Exact solution y(x) = exp(-x)
def exact_solution(x):
    return np.exp(-x)

# Leap-Frog method
def leap_frog(y0, y1, T, n):
    h = T / n
    x = np.linspace(0, T, n + 1)
    y = np.zeros(n + 1)

    # Initial values
    y[0] = y0
    y[1] = y1

    # Leap-Frog scheme
    for i in range(1, n):
        y[i + 1] = y[i - 1] - 2 * h * y[i]

    return x, y, h

# Lists to store results
h_values = []
lte_values = []
error_values = []

# Plot 1: Numerical solutions
plt.figure(figsize=(8, 5))

for n in n_values:
    x, y, h = leap_frog(1.0, 1.0, T, n)

    # Exact solution
    y_exact = exact_solution(x)

    # Local truncation error at x = T - h
    lte = ((exact_solution(T) - exact_solution(T - 2*h))
           / (2*h) + exact_solution(T - h))

    # Global error at T = 1
    error = abs(y[-1] - y_exact[-1])

    h_values.append(h)
    lte_values.append(abs(lte))
    error_values.append(error)

    print(f"n = {n:5d}, h = {h:.5f}, "
          f"LTE = {abs(lte):.10f}, "
          f"Global Error = {error:.10f}")

    plt.plot(x, y, label=f'n = {n}')

# Exact solution
x_exact = np.linspace(0, T, 1000)
plt.plot(x_exact, exact_solution(x_exact),
         'k--', linewidth=2, label='Exact solution')

plt.xlabel('Time (x)')
plt.ylabel('Solution y(x)')
plt.title('Leap-Frog: Numerical vs Exact Solution')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

# Plot 2: LTE and global error
plt.figure(figsize=(7, 5))
plt.loglog(h_values, lte_values, 'o-', label='LTE')
plt.loglog(h_values, error_values, 's-', label='Global Error')

plt.xlabel('Step size (h)')
plt.ylabel('Absolute Error')
plt.title('Leap-Frog: Consistency and Global Error')
plt.legend()
plt.grid(True, which='both')
plt.tight_layout()
plt.show()

# Plot 3: Growth of the oscillatory computational mode
n = 10
h = T / n
x = np.linspace(0, T, n + 1)

# Characteristic roots
r_plus = -h + np.sqrt(1 + h**2)
r_minus = -h - np.sqrt(1 + h**2)

# Coefficient determined by y0 = y1 = 1
B = (1.0 - r_plus) / (r_minus - r_plus)

# Actual oscillatory component
k = np.arange(n + 1)
oscillatory_mode = B * r_minus**k

plt.figure(figsize=(7, 5))
plt.plot(x, oscillatory_mode, 'o-',
         label='Oscillatory Mode')
plt.plot(x, np.abs(oscillatory_mode), '--',
         label='Mode Amplitude')

plt.xlabel('Time (x)')
plt.ylabel('Mode Value')
plt.title('Leap-Frog: Oscillatory Mode Amplification')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

print(f"\nUnstable root: {r_minus:.6f}")
print(f"Initial mode amplitude: {abs(B):.6f}")
print(f"Final mode amplitude: {abs(oscillatory_mode[-1]):.6f}")
print(f"Amplification factor: {abs(r_minus)**n:.4f}")
```

**실행 결과**

```text
n =    10, h = 0.10000, LTE = 0.0006779550, Global Error = 0.1102517932
n =    20, h = 0.05000, LTE = 0.0001611622, Global Error = 0.0570360863
n =    40, h = 0.02500, LTE = 0.0000392921, Global Error = 0.0289614969
n =    80, h = 0.01250, LTE = 0.0000097008, Global Error = 0.0145869868
n =   160, h = 0.00625, LTE = 0.0000024101, Global Error = 0.0073194539

Unstable root: -1.104988
Initial mode amplitude: 0.047270
Final mode amplitude: 0.128281
Amplification factor: 2.7138
```

![Lab 2 — Exercise 3: Leap-Frog and the Parasitic Mode — Figure 1](/assets/images/ode_lab2_ex3_1.png)

*Figure. Lab 2 — Exercise 3: Leap-Frog and the Parasitic Mode (1/3).*

![Lab 2 — Exercise 3: Leap-Frog and the Parasitic Mode — Figure 2](/assets/images/ode_lab2_ex3_2.png)

*Figure. Lab 2 — Exercise 3: Leap-Frog and the Parasitic Mode (2/3).*

![Lab 2 — Exercise 3: Leap-Frog and the Parasitic Mode — Figure 3](/assets/images/ode_lab2_ex3_3.png)

*Figure. Lab 2 — Exercise 3: Leap-Frog and the Parasitic Mode (3/3).*

## Summary and Takeaways

**Consistency는 국소 오차가 사라진다는 조건이고, Convergence는 실제 해에 접근한다는 조건**이다. Leap-Frog는 작은 LTE만으로 충분하지 않다는 사실을 보여 준다. Long-time Instability의 원인은 Post 4의 Root Condition과 연결해 다시 정리한다.

## References

- Raj Dhara, *Lec-2_3_scheme*, 2026-10-12, pp. 1–4.
- Raj Dhara, *Numerical Method Lab: ODE Numerical Schemes & Convergence*, Exercises 1 and 3.
- `LAB_1.ipynb`, Exercises 1–3 (Python implementation and stored results).
