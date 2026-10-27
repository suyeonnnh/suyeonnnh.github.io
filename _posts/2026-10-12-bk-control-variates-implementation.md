---
title: "Control Variates for Coupon-Bond Pricing under the Black–Karasinski Model 2"
date: 2026-10-12
categories:
  - Computational Finance
tags:
  - Black-Karasinski
  - Vasicek
  - Monte-Carlo
  - Control-Variates
math: false
---

## Introduction

지난 글에서는 Black–Karasinski(BK)와 Vasicek 모델의 수학적 구조, 할인계수, 쿠폰 채권의 Monte Carlo 가격 추정과 Control Variates의 기본 원리를 정리했다.

이번 글에서는 이론을 Python으로 구현하여 **Monte Carlo 가격 추정의 정확도와 분산 감소 효과**를 비교한다. 먼저 해석적 가격을 계산할 수 있는 Vasicek 모델에서 구현을 검증하고, Black–Karasinski 모델로 확장해 두 가지 변동성에서 쿠폰 채권 가격을 추정한다. 이후 Vasicek 경로를 이용한 통제변량과 BK 금리 적분값으로 구성한 통제변량을 비교하여 가격, 표준오차(Standard Error, SE), 분산 감소 배율(Variance Reduction Factor, VRF)을 살펴본다.

프로젝트를 하며 얻은 중요한 점은 **추정량의 표준오차가 작다고 해서 가격이 반드시 정확해지는 것은 아니라는 점**이었다. 통제변량의 기대값이 시뮬레이션의 이산화 방식과 일치하지 않거나 잘못 지정되면, 작은 SE를 유지하면서도 가격 추정에 편향이 발생할 수 있다. 이를 확인하기 위해 독립 표본을 이용한 기대값 검증, 시간 간격과 적분 방식 비교, 의도적인 기대값 오류 실험, 자동화된 테스트를 차례로 수행한다.

## Parameters

| 구분 | 파라미터 | 값 | 설명 |
|---|---|---|---|
| Black–Karasinski | κ (kappa) | 0.8 | 로그 금리의 평균회귀 속도 |
| Black–Karasinski | θ (theta) | ln(0.03) | 로그 금리의 장기 평균 수준 |
| Black–Karasinski | r₀ | 0.03 | 초기 단기금리 |
| Black–Karasinski | σ (sigma) | **0.2, 0.5** | 비교할 두 변동성 시나리오 |
| Vasicek (Baseline) | α (alpha) | 0.8 | 평균회귀 속도 |
| Vasicek (Baseline) | β (beta) | 0.03 | 장기 평균금리 |
| Vasicek (Baseline) | η (eta) | 0.01 | 금리 변동성 |
| Vasicek (Baseline) | r₀ᵛ | 0.03 | 초기 단기금리 |
| Coupon Bond | T | 5 years | 만기 |
| Coupon Bond | C | 0.04 | 매년 지급하는 쿠폰, 액면가 1 기준 |
| Coupon Bond | Principal | 1.0 | 만기 원금 |
| Simulation | Δt | 1/52 | 기본 주간 시간 간격 |
| Time-step Comparison | Δt | 1/252 | 더 촘촘한 시간 간격 |
| Pricing Sample | N | 20,000 | 채권 가격 평가 경로 수 |
| Pilot Sample | N_pilot | 2,000 | 통제변량 계수 추정용 독립 경로 수 |
| Random Seed | Seed | 2026 | 기본 가격 평가 난수 시드 |
| Independent Validation | N_test | 30,000 | 통제변량 기대값 검증용 경로 수 |

BK 모델에서는 **σ = 0.2와 σ = 0.5를 각각 고정**하고 동일한 추정법의 가격과 SE를 비교한다. 통제변량을 적용할 때는 가격 평가에 쓰는 경로와 계수를 추정하는 pilot 경로를 분리하며, 비교 대상에는 가능한 한 같은 난수 증분(Common Random Numbers)을 사용한다. 이때 목적은 변동성을 인위적으로 줄이는 것이 아니라, **같은 BK 모형과 같은 변동성 조건에서 Monte Carlo 추정량의 분산을 줄이는 것**이다.

경로별 채권 가격 Y의 Monte Carlo 추정량과 표준오차는 다음과 같다.

- **MC price:** P̂_MC = (1/N) Σᵢ Yᵢ
- **Standard error:** SE = s_Y / √N
- **Sample variance:** s_Y² = [1/(N−1)] Σᵢ (Yᵢ − Ȳ)²

기대값 μ_X = E[X]를 알고 있는 통제변량 X에 대해서는 다음과 같이 가격을 보정한다.

- **Control variate:** Y_CV = Y − ν̂ (X − μ_X)
- **Variance Reduction Factor:** VRF = Var(Y) / Var(Y_CV) (표본분산 기준)

**VRF가 클수록 같은 표본 수에서 더 작은 분산을 얻었다**는 뜻이며, 기대값의 정확성은 별도로 검증해야 한다. 여기서 ν̂는 독립 pilot 표본으로 추정한 통제변량 계수이다.

---

# Part A: Plain Monte Carlo in the Vasicek Model

## Task 1. Vasicek Plain Monte Carlo

Vasicek 모델은 해석적 채권 가격을 구할 수 있으므로 Monte Carlo 구현을 검증하는 기준으로 사용한다. 먼저 Euler 방법으로 금리 경로를 생성하고, 사다리꼴 적분으로 할인계수를 계산한다. 연간 쿠폰과 만기 원금을 각 지급 시점에 할인하여 경로별 채권 가치를 구한다.

**Python Implementation**

```python
import numpy as np

# 1. Parameters
alpha = 0.8
beta = 0.03
eta = 0.01
r0 = 0.03

T = 5
C = 0.04
principal = 1.0

dt = 1 / 52
N = 20000
seed = 2026

# 2. Vasicek Monte Carlo simulation
def vasicek_mc(N, dt, seed=2026):
    rng = np.random.default_rng(seed)

    M = round(T / dt)
    dt = T / M
    sqrt_dt = np.sqrt(dt)

    rates = np.full(N, r0, dtype=float)
    integrated_rates = np.zeros(N)
    payoffs = np.zeros(N)

    steps_per_year = M // T

    for k in range(1, M + 1):

        # Euler scheme
        Z = rng.standard_normal(N)

        next_rates = (
            rates
            + alpha * (beta - rates) * dt
            + eta * sqrt_dt * Z
        )

        # Trapezoidal rule
        integrated_rates += (
            0.5 * (rates + next_rates) * dt
        )

        rates = next_rates

        # Annual coupon payments
        if k % steps_per_year == 0:
            year = k // steps_per_year

            cashflow = C

            if year == T:
                cashflow += principal

            discount_factor = np.exp(-integrated_rates)

            payoffs += cashflow * discount_factor

    # Plain Monte Carlo price
    mc_price = np.mean(payoffs)

    # Estimated standard error
    sample_std = np.std(payoffs, ddof=1)
    se = sample_std / np.sqrt(N)

    return mc_price, se, payoffs


mc_price, se, payoffs = vasicek_mc(N, dt, seed)

print("=== Task 1: Vasicek Plain Monte Carlo ===")
print(f"Number of paths: {N:,}")
print(f"Time step: {dt:.8f}")
print(f"MC bond price: {mc_price:.8f}")
print(f"Standard error: {se:.8f}")
```

**저장된 실행 결과**

```text
=== Task 1: Vasicek Plain Monte Carlo ===
Number of paths: 20,000
Time step: 0.01923077
MC bond price: 1.04383592
Standard error: 0.00015148
```

## Task 2. Vasicek Analytical Comparison

Vasicek의 해석적 영이표채권 가격을 이용해 쿠폰 채권의 기준 가격을 구한다. MC 가격과 해석적 가격의 차이를 SE 단위로 비교하고, 시간 간격을 1/52에서 1/252로 변경하여 이산화 오차의 영향을 확인한다.

**Python Implementation**

```python

# 3. Analytical Vasicek zero-coupon bond price
def vasicek_zero_coupon(t):
    B = -np.expm1(-alpha * t) / alpha

    A = (
        (beta - eta**2 / (2 * alpha**2)) * (B - t)
        - (eta**2 / (4 * alpha)) * B**2
    )

    return np.exp(A - B * r0)


# 4. Analytical coupon-bond price
analytical_price = 0.0

for year in range(1, T + 1):
    cashflow = C

    if year == T:
        cashflow += principal

    analytical_price += (
        cashflow * vasicek_zero_coupon(year)
    )


# 5. Comparison with Monte Carlo
difference = mc_price - analytical_price
absolute_error = abs(difference)
error_in_se = absolute_error / se

print("=== Task 2: Analytical Comparison ===")
print(f"Analytical price: {analytical_price:.8f}")
print(f"MC price:         {mc_price:.8f}")
print(f"Standard error:   {se:.8f}")
print(f"Absolute error:   {absolute_error:.8f}")
print(f"Error / SE:       {error_in_se:.4f}")


# 6. Finer time step
dt_fine = 1 / 252

mc_fine, se_fine, _ = vasicek_mc(
    N, dt_fine, seed
)

print("\n=== Time Step Comparison ===")
print(f"dt = 1/52  : {mc_price:.8f}, SE = {se:.8f}")
print(f"dt = 1/252 : {mc_fine:.8f}, SE = {se_fine:.8f}")
print(f"Analytical: {analytical_price:.8f}")
```

**저장된 실행 결과**

```text
=== Task 2: Analytical Comparison ===
Analytical price: 1.04389268
MC price:         1.04383592
Standard error:   0.00015148
Absolute error:   0.00005676
Error / SE:       0.3747

=== Time Step Comparison ===
dt = 1/52  : 1.04383592, SE = 0.00015148
dt = 1/252 : 1.04382200, SE = 0.00015169
Analytical: 1.04389268
```

**Python Implementation**

```python
# Basic numerical validation

# Test 1: The simulation generates N payoffs
assert len(payoffs) == N

# Test 2: All simulated prices are finite
assert np.all(np.isfinite(payoffs))

# Test 3: The discounting produces positive path values
assert np.all(payoffs > 0)

# Test 4: The estimated standard error is positive
assert se > 0

# Test 5: The analytic zero-coupon price at t=0 equals 1
assert np.isclose(
    vasicek_zero_coupon(0),
    1.0,
    atol=1e-12
)

print("All basic validation tests passed.")
```

**저장된 실행 결과**

```text
All basic validation tests passed.
```

# Part B: Plain Monte Carlo in the BK Model

## Task 3. BK Plain Monte Carlo

BK에서는 로그 금리를 Euler 방법으로 생성하고 실제 금리를 지수변환하여 계산한다. 두 변동성 σ=0.2와 σ=0.5에서 Plain MC 채권 가격과 표준오차를 비교한다. 또한 σ=0인 경우 결정론적 할인 가격과 일치하는지 검증한다.

**Python Implementation**

```python
import numpy as np

# 1. BK parameters
kappa = 0.8
theta = np.log(0.03)
r0_bk = 0.03
sigmas = [0.2, 0.5]

# 2. Bond and simulation parameters
T = 5
C = 0.04
principal = 1.0

N = 20000
dt = 1 / 52
seed = 2026


def bk_mc(sigma, N=20000, dt=1/52, seed=2026):
    rng = np.random.default_rng(seed)

    M = round(T / dt)
    dt = T / M
    sqrt_dt = np.sqrt(dt)
    steps_per_year = M // T

    # Initial log short rates
    x = np.full(N, np.log(r0_bk), dtype=float)

    integrated_rates = np.zeros(N)
    payoffs = np.zeros(N)

    for k in range(1, M + 1):

        # Generate Brownian increments
        Z = rng.standard_normal(N)

        # Euler scheme for log short rates
        x_next = (
            x
            + kappa * (theta - x) * dt
            + sigma * sqrt_dt * Z
        )

        # Convert log rates to short rates
        r_current = np.exp(x)
        r_next = np.exp(x_next)

        # Trapezoidal integration
        integrated_rates += (
            0.5 * (r_current + r_next) * dt
        )

        x = x_next

        # Annual coupon payments
        if k % steps_per_year == 0:
            year = k // steps_per_year

            cashflow = C

            if year == T:
                cashflow += principal

            discount_factor = np.exp(-integrated_rates)

            payoffs += cashflow * discount_factor

    # Monte Carlo estimate
    mc_price = np.mean(payoffs)

    # Estimated standard error
    sample_std = np.std(payoffs, ddof=1)
    se = sample_std / np.sqrt(N)

    return mc_price, se, payoffs


# 3. Results for both volatility values
results_bk = {}

for sigma in sigmas:
    price, se, payoffs = bk_mc(
        sigma=sigma,
        N=N,
        dt=dt,
        seed=seed
    )

    results_bk[sigma] = {
        "price": price,
        "se": se,
        "payoffs": payoffs
    }

    print(f"=== BK Model: sigma = {sigma} ===")
    print(f"MC bond price: {price:.8f}")
    print(f"Standard error: {se:.8f}")
    print()
```

**저장된 실행 결과**

```text
=== BK Model: sigma = 0.2 ===
MC bond price: 1.04212403
Standard error: 0.00009197

=== BK Model: sigma = 0.5 ===
MC bond price: 1.03400908
Standard error: 0.00024446
```

![BK 변동성별 Plain MC 가격과 ±1.96 SE](/assets/images/bk_plain_mc_price.png)
*Figure. BK 변동성별 Plain MC 가격과 ±1.96 SE.*

**Python Implementation**

```python
# Validation for Part B

# Test 1: Correct number of simulated payoffs
for sigma in sigmas:
    payoffs = results_bk[sigma]["payoffs"]
    assert len(payoffs) == N

# Test 2: All payoffs are finite and positive
for sigma in sigmas:
    payoffs = results_bk[sigma]["payoffs"]
    assert np.all(np.isfinite(payoffs))
    assert np.all(payoffs > 0)

# Test 3: Standard errors are finite and positive
for sigma in sigmas:
    se = results_bk[sigma]["se"]
    assert np.isfinite(se)
    assert se > 0

# Test 4: Deterministic case (sigma = 0)
# Since theta = log(r0), the rate stays at r0.
price_zero, se_zero, _ = bk_mc(
    sigma=0.0,
    N=100,
    dt=dt,
    seed=seed
)

expected_price = (
    sum(C * np.exp(-r0_bk * year)
        for year in range(1, T + 1))
    + principal * np.exp(-r0_bk * T)
)

assert np.isclose(
    price_zero,
    expected_price,
    atol=1e-12
)
assert se_zero < 1e-12

print("All Part B validation tests passed.")
```

**저장된 실행 결과**

```text
All Part B validation tests passed.
```

# Part C: Vasicek Control Variates

## Part C Overview

Part C에서는 Vasicek 모델을 BK 가격 추정의 통제변량으로 사용한다. **동일한 Gaussian increments**로 두 모델을 함께 시뮬레이션하여 경로 간 상관관계를 높이고, 회귀계수는 독립적인 pilot 표본에서 추정한 후 가격 표본에 적용한다.

**Python Implementation**

```python
import numpy as np
import matplotlib.pyplot as plt

# BK parameters
KAPPA = 0.8
THETA = np.log(0.03)
R0 = 0.03
SIGMAS = (0.2, 0.5)

# Bond and simulation parameters
T = 5
COUPON = 0.04
N = 20000
N_PILOT = 2000
DT = 1 / 52
SEED = 2026

# Initial Vasicek parameters
BASE_PARAMS = (0.8, 0.03, 0.01)


def generate_normals(n, dt, seed):
    m = round(T / dt)

    assert np.isclose(m * dt, T)
    assert m % T == 0

    rng = np.random.default_rng(seed)
    return rng.standard_normal((m, n))


def simulate_joint(z, sigma, alpha, beta, eta, dt=DT):
    """
    Simulate BK and Vasicek models using identical Gaussian increments.
    Return BK coupon payoff, Vasicek zero-coupon control,
    and Vasicek coupon-bond control.
    """
    m, n = z.shape
    steps_per_year = m // T

    x = np.full(n, np.log(R0), dtype=float)
    v = np.full(n, R0, dtype=float)

    I_bk = np.zeros(n)
    I_vas = np.zeros(n)

    Y = np.zeros(n)
    X_coupon = np.zeros(n)

    sqrt_dt = np.sqrt(dt)

    for k in range(1, m + 1):
        zk = z[k - 1]

        # BK Euler scheme: simulate the log short rate
        x_next = (
            x
            + KAPPA * (THETA - x) * dt
            + sigma * sqrt_dt * zk
        )

        # Vasicek Euler scheme: same Gaussian increments
        v_next = (
            v
            + alpha * (beta - v) * dt
            + eta * sqrt_dt * zk
        )

        # Integrated rates using the trapezoidal rule
        I_bk += 0.5 * (np.exp(x) + np.exp(x_next)) * dt
        I_vas += 0.5 * (v + v_next) * dt

        x = x_next
        v = v_next

        if k % steps_per_year == 0:
            cashflow = COUPON + (1.0 if k == m else 0.0)

            Y += cashflow * np.exp(-I_bk)
            X_coupon += cashflow * np.exp(-I_vas)

    X_zcb = np.exp(-I_vas)

    return Y, X_zcb, X_coupon


# Separate Gaussian samples for coefficient estimation and pricing
z_pilot = generate_normals(N_PILOT, DT, SEED + 1)
z_pricing = generate_normals(N, DT, SEED)

print("Part C simulation functions ready.")
print(f"Pilot paths: {N_PILOT}")
print(f"Pricing paths: {N}")
```

**저장된 실행 결과**

```text
Part C simulation functions ready.
Pilot paths: 2000
Pricing paths: 20000
```

## Task 4. Two Vasicek Control Variates

통제변량은 Vasicek의 만기 할인계수(ZCB)와 쿠폰 채권 가격 두 종류다. 기대값은 **실제 시뮬레이션에 사용한 Euler·Trapezoidal 이산화와 일치하도록** 계산한다. 가격, SE, 추정 계수 ν̂, VRF를 비교한다.

**Python Implementation**

```python
def vasicek_control_expectations(alpha, beta, eta, dt=DT):
    """
    Exact expectations of the discretized Vasicek controls.

    Euler rates and trapezoidal integrated rates are Gaussian.
    Their means and covariance matrices are propagated exactly.
    """
    m = round(T / dt)
    steps_per_year = m // T

    # State = [short rate, integrated short rate]
    mean = np.array([R0, 0.0])
    cov = np.zeros((2, 2))

    a = 1 - alpha * dt
    b = alpha * beta * dt
    q = eta * np.sqrt(dt)

    A = np.array([
        [a, 0.0],
        [0.5 * dt * (1 + a), 1.0]
    ])

    shift = np.array([b, 0.5 * dt * b])
    noise = np.array([q, 0.5 * dt * q])

    coupon_expectation = 0.0

    for k in range(1, m + 1):
        mean = A @ mean + shift
        cov = A @ cov @ A.T + np.outer(noise, noise)

        if k % steps_per_year == 0:
            cashflow = COUPON + (1.0 if k == m else 0.0)

            # For Gaussian I: E[exp(-I)] = exp(-E[I] + Var(I)/2)
            expected_discount = np.exp(
                -mean[1] + 0.5 * cov[1, 1]
            )

            coupon_expectation += cashflow * expected_discount

    zcb_expectation = np.exp(
        -mean[1] + 0.5 * cov[1, 1]
    )

    return zcb_expectation, coupon_expectation


def control_variate(y_pilot, x_pilot, y_price, x_price, mu):
    """
    Fit the control coefficient on the pilot sample
    and evaluate the CV estimator on independent pricing paths.
    """
    x_pilot = np.asarray(x_pilot)
    x_price = np.asarray(x_price)

    # Sample covariance / sample variance
    nu = (
        np.cov(x_pilot, y_pilot, ddof=1)[0, 1]
        / np.var(x_pilot, ddof=1)
    )

    y_cv = y_price - nu * (x_price - mu)

    price = np.mean(y_cv)
    se = np.std(y_cv, ddof=1) / np.sqrt(len(y_cv))

    vrf = (
        np.var(y_price, ddof=1)
        / np.var(y_cv, ddof=1)
    )

    return {
        "price": float(price),
        "se": float(se),
        "nu": float(nu),
        "vrf": float(vrf)
    }


def evaluate_controls(sigma, params, z_pilot, z_pricing):
    alpha, beta, eta = params

    Y_pilot, Xz_pilot, Xc_pilot = simulate_joint(
        z_pilot, sigma, alpha, beta, eta
    )

    Y, Xz, Xc = simulate_joint(
        z_pricing, sigma, alpha, beta, eta
    )

    mu_z, mu_c = vasicek_control_expectations(
        alpha, beta, eta
    )

    zcb_result = control_variate(
        Y_pilot, Xz_pilot, Y, Xz, mu_z
    )

    coupon_result = control_variate(
        Y_pilot, Xc_pilot, Y, Xc, mu_c
    )

    plain_price = np.mean(Y)
    plain_se = np.std(Y, ddof=1) / np.sqrt(len(Y))

    return {
        "plain_price": float(plain_price),
        "plain_se": float(plain_se),
        "zcb": zcb_result,
        "coupon": coupon_result
    }


print("Control-variate functions ready.")
```

**저장된 실행 결과**

```text
Control-variate functions ready.
```

**Python Implementation**

```python
task4_results = {}

for sigma in SIGMAS:
    result = evaluate_controls(
        sigma,
        BASE_PARAMS,
        z_pilot,
        z_pricing
    )

    task4_results[sigma] = result

    print(f"\n=== BK sigma = {sigma} ===")
    print(f"Plain MC price: {result['plain_price']:.8f}")
    print(f"Plain MC SE:    {result['plain_se']:.8f}")

    for name in ["zcb", "coupon"]:
        r = result[name]

        print(f"\nControl: {name}")
        print(f"CV price: {r['price']:.8f}")
        print(f"CV SE:    {r['se']:.8f}")
        print(f"nu:       {r['nu']:.6f}")
        print(f"VRF:      {r['vrf']:.4f}")
```

**저장된 실행 결과**

```text

=== BK sigma = 0.2 ===
Plain MC price: 1.04212403
Plain MC SE:    0.00009197

Control: zcb
CV price: 1.04216099
CV SE:    0.00000884
nu:       0.672563
VRF:      108.3568

Control: coupon
CV price: 1.04215894
CV SE:    0.00000816
nu:       0.604759
VRF:      127.1633

=== BK sigma = 0.5 ===
Plain MC price: 1.03400908
Plain MC SE:    0.00024446

Control: zcb
CV price: 1.03410593
CV SE:    0.00004926
nu:       1.762282
VRF:      24.6262

Control: coupon
CV price: 1.03410052
CV SE:    0.00004851
nu:       1.584209
VRF:      25.3936
```

### Task 4. Results and Interpretation

Task 4에서는 두 통제변량 모두 Plain MC보다 SE를 낮췄다. 이번 실험에서 Vasicek 쿠폰 채권 통제변량은 동일한 현금흐름 구조를 반영하기 때문에 ZCB 통제변량보다 조금 높은 VRF를 보였다. 다만 결과는 선택된 모형 파라미터에 의존한다.

## Task 5. Vasicek Parameter Matching

BK의 단기금리 동역학을 장기 평균 금리 부근에서 1차로 근사하여 Vasicek 파라미터 대응 규칙을 유도한다. 이는 국소적인 근사 규칙이며 두 모형의 정확한 동등성을 의미하지 않는다.

**Python Implementation**

```python
def match_vasicek_parameters(sigma):
    r_bar = np.exp(THETA)

    alpha = KAPPA
    beta = r_bar
    eta = sigma * r_bar

    return alpha, beta, eta


matched_params = {}

for sigma in SIGMAS:
    params = match_vasicek_parameters(sigma)
    matched_params[sigma] = params

    alpha, beta, eta = params

    print(f"BK sigma = {sigma}")
    print(f"alpha = {alpha:.6f}")
    print(f"beta  = {beta:.6f}")
    print(f"eta   = {eta:.6f}")
    print()
```

**저장된 실행 결과**

```text
BK sigma = 0.2
alpha = 0.800000
beta  = 0.030000
eta   = 0.006000

BK sigma = 0.5
alpha = 0.800000
beta  = 0.030000
eta   = 0.015000
```

## Task 6. Numerical Search for Maximum VRF

α와 η를 격자에서 변경하면서 통제변량별 VRF를 계산한다. 모든 격자점에서 같은 Gaussian 난수를 사용해 비교의 표본잡음을 줄인다. 결과는 **탐색한 격자 안에서의 최대값**이며 전역 최적성은 보장하지 않는다.

**Python Implementation**

```python
# Search ranges
alpha_grid = np.linspace(0.35, 1.65, 9)
eta_grid = np.linspace(0.002, 0.032, 11)

beta_search = 0.03

# Common random numbers for the entire grid
z_search = generate_normals(
    n=2000,
    dt=DT,
    seed=123456
)


def calculate_grid_vrf(sigma):
    vrf_zcb = np.zeros(
        (len(eta_grid), len(alpha_grid))
    )

    vrf_coupon = np.zeros_like(vrf_zcb)

    for i, eta in enumerate(eta_grid):
        for j, alpha in enumerate(alpha_grid):

            Y, Xz, Xc = simulate_joint(
                z_search,
                sigma,
                alpha,
                beta_search,
                eta
            )

            # Optimal scalar CV variance reduction
            corr_z = np.corrcoef(Y, Xz)[0, 1]
            corr_c = np.corrcoef(Y, Xc)[0, 1]

            vrf_zcb[i, j] = 1 / (1 - corr_z**2)
            vrf_coupon[i, j] = 1 / (1 - corr_c**2)

    return vrf_zcb, vrf_coupon


grid_results = {}

for sigma in SIGMAS:
    grid_results[sigma] = calculate_grid_vrf(sigma)

    print(f"Grid search completed: sigma = {sigma}")


# Report the best parameters found within the grid
for sigma in SIGMAS:
    vrf_zcb, vrf_coupon = grid_results[sigma]

    for name, vrf_matrix in [
        ("ZCB control", vrf_zcb),
        ("Coupon-bond control", vrf_coupon)
    ]:
        i, j = np.unravel_index(
            np.argmax(vrf_matrix),
            vrf_matrix.shape
        )

        best_alpha = alpha_grid[j]
        best_eta = eta_grid[i]
        best_vrf = vrf_matrix[i, j]

        print(f"\nBK sigma = {sigma}")
        print(f"Control: {name}")
        print(f"Best alpha: {best_alpha:.6f}")
        print(f"Best beta: {beta_search:.6f}")
        print(f"Best eta: {best_eta:.6f}")
        print(f"Maximum grid VRF: {best_vrf:.4f}")
```

**저장된 실행 결과**

```text
Grid search completed: sigma = 0.2
Grid search completed: sigma = 0.5

BK sigma = 0.2
Control: ZCB control
Best alpha: 0.675000
Best beta: 0.030000
Best eta: 0.002000
Maximum grid VRF: 183.9261

BK sigma = 0.2
Control: Coupon-bond control
Best alpha: 0.837500
Best beta: 0.030000
Best eta: 0.002000
Maximum grid VRF: 188.3937

BK sigma = 0.5
Control: ZCB control
Best alpha: 0.675000
Best beta: 0.030000
Best eta: 0.002000
Maximum grid VRF: 32.0218

BK sigma = 0.5
Control: Coupon-bond control
Best alpha: 0.837500
Best beta: 0.030000
Best eta: 0.002000
Maximum grid VRF: 32.1381
```

**Python Implementation**

```python
# Additional check: is the best eta at the grid boundary?
# Compare smaller eta values using the same random numbers.

eta_check = [0.0005, 0.001, 0.002, 0.004]

for sigma in SIGMAS:
    print(f"\nBK sigma = {sigma}")

    for control_name, alpha in [
        ("ZCB", 0.675),
        ("Coupon", 0.8375)
    ]:
        print(f"\n{control_name} control")

        for eta in eta_check:
            Y, Xz, Xc = simulate_joint(
                z_search, sigma, alpha, beta_search, eta
            )

            X = Xz if control_name == "ZCB" else Xc
            rho = np.corrcoef(Y, X)[0, 1]
            vrf = 1 / (1 - rho**2)

            print(f"eta={eta:.4f}, VRF={vrf:.4f}")
```

**저장된 실행 결과**

```text

BK sigma = 0.2

ZCB control
eta=0.0005, VRF=193.9270
eta=0.0010, VRF=190.5250
eta=0.0020, VRF=183.9261
eta=0.0040, VRF=171.5209

Coupon control
eta=0.0005, VRF=197.0514
eta=0.0010, VRF=194.1145
eta=0.0020, VRF=188.3937
eta=0.0040, VRF=177.5445

BK sigma = 0.5

ZCB control
eta=0.0005, VRF=32.7107
eta=0.0010, VRF=32.4792
eta=0.0020, VRF=32.0218
eta=0.0040, VRF=31.1295

Coupon control
eta=0.0005, VRF=32.7160
eta=0.0010, VRF=32.5220
eta=0.0020, VRF=32.1381
eta=0.0040, VRF=31.3863
```

**Python Implementation**

```python
%matplotlib inline

import matplotlib.pyplot as plt

for sigma in SIGMAS:
    vrf_zcb, vrf_coupon = grid_results[sigma]

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

    controls = [
        ("Vasicek ZCB control", vrf_zcb),
        ("Vasicek coupon-bond control", vrf_coupon)
    ]

    for ax, (title, values) in zip(axes, controls):
        im = ax.imshow(
            values,
            origin="lower",
            aspect="auto",
            extent=[
                alpha_grid[0],
                alpha_grid[-1],
                eta_grid[0],
                eta_grid[-1]
            ]
        )

        ax.set_title(title)
        ax.set_xlabel("alpha")
        ax.set_ylabel("eta")

        fig.colorbar(im, ax=ax, label="VRF")

    fig.suptitle(f"BK sigma = {sigma}")

    plt.tight_layout()

    filename = f"vrf_heatmap_sigma_{sigma}.png"
    plt.savefig(filename, dpi=200, bbox_inches="tight")

    plt.show()

    print(f"Saved: {filename}")
```

**저장된 실행 결과**

```text
Saved: vrf_heatmap_sigma_0.2.png

Saved: vrf_heatmap_sigma_0.5.png
```

![BK σ=0.2의 Vasicek 통제변량 VRF 히트맵](/assets/images/vrf_heatmap_sigma_0.2.png)
*Figure. BK σ=0.2의 Vasicek 통제변량 VRF 히트맵.*

**Python Implementation**

```python
# Investigate the effect of beta
beta_values = [0.01, 0.03, 0.05]

alpha = 0.8
eta = 0.015
sigma = 0.5

z_beta = generate_normals(
    n=20000,
    dt=DT,
    seed=SEED + 5
)

for beta in beta_values:
    Y, Xz, Xc = simulate_joint(
        z_beta,
        sigma,
        alpha,
        beta,
        eta
    )

    corr_z = np.corrcoef(Y, Xz)[0, 1]
    corr_c = np.corrcoef(Y, Xc)[0, 1]

    vrf_z = 1 / (1 - corr_z**2)
    vrf_c = 1 / (1 - corr_c**2)

    print(f"beta = {beta:.2f}")
    print(f"ZCB VRF:    {vrf_z:.4f}")
    print(f"Coupon VRF: {vrf_c:.4f}")
    print()
```

**저장된 실행 결과**

```text
beta = 0.01
ZCB VRF:    24.8902
Coupon VRF: 25.6588

beta = 0.03
ZCB VRF:    24.8902
Coupon VRF: 25.6574

beta = 0.05
ZCB VRF:    24.8902
Coupon VRF: 25.6528
```

![BK σ=0.5의 Vasicek 통제변량 VRF 히트맵](/assets/images/vrf_heatmap_sigma_0.5.png)
*Figure. BK σ=0.5의 Vasicek 통제변량 VRF 히트맵.*

### Task 6. Results and Interpretation

히트맵을 비교하면 η가 기존 격자의 하한에 있을 때 VRF가 크게 나오는 경우가 있다. 따라서 추가로 더 작은 η 값(0.0005, 0.001, 0.002, 0.004)을 같은 난수로 평가했다. 탐색 범위가 결과에 영향을 주므로 경계점의 최대 VRF를 전역 최적값으로 해석하지 않는다. β를 변화시킬 때 ZCB의 최적 스칼라 VRF는 거의 변하지 않지만 쿠폰 채권 VRF는 조금 달라질 수 있다.

## Task 7. Comparison of Vasicek Parameter Choices

Part C의 기본값(Baseline), Task 5의 모형 근사값(Matched), Task 6의 탐색값(Grid best)을 **같은 독립 가격 표본**에서 비교한다. 가격·SE·ν̂·VRF를 함께 확인하고 단순 근사 규칙이 언제나 가장 높은 VRF를 주지는 않는다는 점을 살펴본다.

**Python Implementation**

```python
task7_results = []

for sigma in SIGMAS:

    # Baseline and Task 5 matching
    baseline = BASE_PARAMS
    matched = matched_params[sigma]

    for control_index, control_name in [
        (0, "zcb"),
        (1, "coupon")
    ]:
        # Best parameter combination from Task 6 grid
        vrf_grid = grid_results[sigma][control_index]

        best_i, best_j = np.unravel_index(
            np.argmax(vrf_grid),
            vrf_grid.shape
        )

        best_params = (
            float(alpha_grid[best_j]),
            beta_search,
            float(eta_grid[best_i])
        )

        choices = {
            "Baseline": baseline,
            "Matched": matched,
            "Grid best": best_params
        }

        print(f"\nBK sigma = {sigma}, control = {control_name}")

        for choice_name, params in choices.items():
            results = evaluate_controls(
                sigma,
                params,
                z_pilot,
                z_pricing
            )

            r = results[control_name]

            task7_results.append({
                "sigma": sigma,
                "control": control_name,
                "choice": choice_name,
                "alpha": params[0],
                "beta": params[1],
                "eta": params[2],
                "price": r["price"],
                "se": r["se"],
                "nu": r["nu"],
                "vrf": r["vrf"]
            })

            print(
                f"{choice_name:10s} "
                f"params={params} "
                f"price={r['price']:.8f} "
                f"SE={r['se']:.8f} "
                f"nu={r['nu']:.5f} "
                f"VRF={r['vrf']:.2f}"
            )
```

**저장된 실행 결과**

```text

BK sigma = 0.2, control = zcb
Baseline   params=(0.8, 0.03, 0.01) price=1.04216099 SE=0.00000884 nu=0.67256 VRF=108.36
Matched    params=(0.8, 0.029999999999999995, 0.005999999999999999) price=1.04216098 SE=0.00000837 nu=1.12178 VRF=120.78
Grid best  params=(0.6749999999999999, 0.03, 0.002) price=1.04215872 SE=0.00000730 nu=2.98761 VRF=158.51

BK sigma = 0.2, control = coupon
Baseline   params=(0.8, 0.03, 0.01) price=1.04215894 SE=0.00000816 nu=0.60476 VRF=127.16
Matched    params=(0.8, 0.029999999999999995, 0.005999999999999999) price=1.04215895 SE=0.00000766 nu=1.00866 VRF=144.01
Grid best  params=(0.8374999999999999, 0.03, 0.002) price=1.04215954 SE=0.00000721 nu=3.13084 VRF=162.75

BK sigma = 0.5, control = zcb
Baseline   params=(0.8, 0.03, 0.01) price=1.03410593 SE=0.00004926 nu=1.76228 VRF=24.63
Matched    params=(0.8, 0.029999999999999995, 0.014999999999999998) price=1.03410588 SE=0.00005083 nu=1.17259 VRF=23.13
Grid best  params=(0.6749999999999999, 0.03, 0.002) price=1.03410006 SE=0.00004620 nu=7.83555 VRF=28.00

BK sigma = 0.5, control = coupon
Baseline   params=(0.8, 0.03, 0.01) price=1.03410052 SE=0.00004851 nu=1.58421 VRF=25.39
Matched    params=(0.8, 0.029999999999999995, 0.014999999999999998) price=1.03410040 SE=0.00005007 nu=1.05416 VRF=23.84
Grid best  params=(0.8374999999999999, 0.03, 0.002) price=1.03410222 SE=0.00004607 nu=8.21159 VRF=28.15
```

### Task 7. Results and Interpretation

σ=0.2에서는 탐색값의 VRF가 기본값보다 높았다. σ=0.5에서는 Matched가 Baseline보다 낮은 VRF를 보이기도 했다. 따라서 모형을 국소적으로 잘 근사하는 파라미터가 통제변량 성능까지 최적으로 만든다고 볼 수 없다. ν̂의 크기도 통제변량의 스케일에 영향을 받으므로 단독 성능지표로 쓰지 않는다.

### Part C Validation

독립 pilot 표본과 가격 표본을 분리했는지, 통제변량의 기대값과 시뮬레이션 평균이 통계적 허용 범위 내에 있는지, 각 가격·SE·VRF가 유한한지 확인한다. 기본 테스트만으로 모든 계산 오류를 배제할 수는 없다.

**Python Implementation**

```python
for sigma in SIGMAS:
    alpha, beta, eta = BASE_PARAMS

    Y, Xz, Xc = simulate_joint(
        z_pricing,
        sigma,
        alpha,
        beta,
        eta
    )

    mu_z, mu_c = vasicek_control_expectations(
        alpha,
        beta,
        eta
    )

    # Finite and positive simulated payoffs
    assert np.all(np.isfinite(Y))
    assert np.all(np.isfinite(Xz))
    assert np.all(np.isfinite(Xc))

    assert np.all(Y > 0)
    assert np.all(Xz > 0)
    assert np.all(Xc > 0)

    # Control expectation checks
    se_z = np.std(Xz, ddof=1) / np.sqrt(N)
    se_c = np.std(Xc, ddof=1) / np.sqrt(N)

    assert abs(np.mean(Xz) - mu_z) < 4 * se_z
    assert abs(np.mean(Xc) - mu_c) < 4 * se_c

    # CV output checks
    results = task4_results[sigma]

    for control_name in ["zcb", "coupon"]:
        r = results[control_name]

        assert np.isfinite(r["price"])
        assert np.isfinite(r["se"])
        assert np.isfinite(r["nu"])
        assert np.isfinite(r["vrf"])

        assert r["se"] > 0
        assert r["vrf"] > 0

    print(f"BK sigma={sigma}: validation passed")

print("All Part C basic validation tests passed.")
```

**저장된 실행 결과**

```text
BK sigma=0.2: validation passed
BK sigma=0.5: validation passed
All Part C basic validation tests passed.
```

### Part C Summary

Part C에서는 두 Vasicek 통제변량의 분산 감소를 확인하고, BK 동역학 근사 및 수치 탐색에 따른 민감도를 비교했다. 마지막 선택은 VRF뿐 아니라 기대값의 정확성과 추가 계산 비용까지 고려해야 한다.

# Part D: Controls Built from the BK Path

## Task 8. BK Control Expectations

Part D에서는 별도의 금리 모형을 도입하지 않고 **BK 경로 자체**에서 통제변량을 만든다. 특히 할인계수의 Taylor 근사에 연결되는 적분 금리의 1차·2차 항을 이용한다.

### Theoretical Derivation

BK 로그 금리는 Gaussian OU 과정이며 실제 금리는 로그정규 구조를 가진다. 먼저 연속시간의 금리 모멘트를 설명하고, 실제 계산에는 **Euler 로그 금리와 사다리꼴 누적 금리의 이산시간 모멘트**를 사용한다. 통제변량의 기대값은 반드시 시뮬레이션된 통제변량의 기대값과 같아야 한다.

BK 로그 금리 x(t)는 Gaussian 과정이고 실제 금리는 exp(x(t))이므로 금리의 1차·2차 모멘트를 계산할 수 있다. 먼저 연속시간 기대값을 유도하고 실제 CV 구현에서는 Euler·Trapezoidal 방식과 일치하는 **이산시간 기대값**을 사용한다.

![BK stochastic differential equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7Ddx_t%3D%5Ckappa%28%5Ctheta-x_t%29dt%2B%5Csigma%5C%2CdW_t%2C%5Cqquad%20r_t%3De%5E%7Bx_t%7D)

The log short rate follows an Ornstein–Uhlenbeck process, which is Gaussian.

Assuming that the initial log rate x_0 is deterministic, its mean is:

![Mean of log short rate](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7Dm%28t%29%3DE%5Bx_t%5D%3D%5Ctheta%2B%28x_0-%5Ctheta%29e%5E%7B-%5Ckappa%20t%7D)

Its covariance function is:

![Covariance of log short rate](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7DK%28s%2Ct%29%3D%5Cfrac%7B%5Csigma%5E2%7D%7B2%5Ckappa%7D%5Cleft%28e%5E%7B-%5Ckappa%5Clvert%20t-s%5Crvert%7D-e%5E%7B-%5Ckappa%28t%2Bs%29%7D%5Cright%29)

Since r_t = exp(x_t), the short rate follows a lognormal distribution at each fixed time.

Therefore,

![Expected short rate](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7DE%5Br_t%5D%3D%5Cexp%5Cleft%28m%28t%29%2B%5Cfrac12%20K%28t%2Ct%29%5Cright%29)

The integrated short rate is:

![Integrated short rate](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7DI%28t%29%3D%5Cint_0%5Et%20r_s%5C%2Cds)

Its first moment is:

![First moment of integrated rate](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7DE%5BI%28t%29%5D%3D%5Cint_0%5Et%20E%5Br_s%5D%5C%2Cds)

To derive its second moment, we write:

![Second moment of integrated rate](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7DE%5BI%28t%29%5E2%5D%3D%5Cint_0%5Et%5Cint_0%5Et%20E%5Br_s%20r_u%5D%5C%2Cds%5C%2Cdu)

Because the log rates are jointly Gaussian,

![Joint short-rate moment](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7DE%5Br_s%20r_u%5D%3D%5Cexp%5Cleft%28m%28s%29%2Bm%28u%29%2B%5Cfrac12%20K%28s%2Cs%29%2B%5Cfrac12%20K%28u%2Cu%29%2BK%28s%2Cu%29%5Cright%29)

Thus, the expectations of the two coupon-weighted controls are:

![Expected X1](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7DE%5BX_1%5D%3D%5Csum_%7Bi%3D1%7D%5E%7Bn%7Dc_iE%5BI%28t_i%29%5D)

![Expected X2](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7DE%5BX_2%5D%3D%5Csum_%7Bi%3D1%7D%5E%7Bn%7Dc_iE%5BI%28t_i%29%5E2%5D)

These expressions provide the theoretical expectations under the continuous-time BK model.

However, our Monte Carlo implementation uses the Euler scheme and trapezoidal integration. Therefore, the simulated integrated rate is a discretized approximation.

For the control-variate calculation, we will derive its exact discrete-time moments using the Gaussian structure of the Euler log-rate process.


This ensures that the control expectations correspond to the quantities actually simulated.

**Python Implementation**

```python
import numpy as np
import matplotlib.pyplot as plt

# Part D parameters
D_KAPPA = 0.8
D_THETA = np.log(0.03)
D_R0 = 0.03

D_SIGMAS = (0.2, 0.5)

D_T = 5
D_COUPON = 0.04
D_PRINCIPAL = 1.0

D_DT = 1 / 52
D_N = 20000
D_N_PILOT = 2000
D_SEED = 2026

# Annual coupon payments, including principal at maturity
D_CASHFLOWS = np.full(D_T, D_COUPON)
D_CASHFLOWS[-1] += D_PRINCIPAL


def d_generate_normals(n, dt, seed):
    m = round(D_T / dt)

    assert np.isclose(m * dt, D_T)
    assert m % D_T == 0

    rng = np.random.default_rng(seed)
    return rng.standard_normal((m, n))


def d_simulate_bk(z, sigma, dt=D_DT):
    """
    Simulate BK using Euler and trapezoidal integration.

    Returns:
        Y: coupon-bond value on each path
        I: integrated rate at annual payment dates
    """
    m, n = z.shape
    steps_per_year = m // D_T

    x = np.full(n, np.log(D_R0))
    integrated_rate = np.zeros(n)

    annual_integrals = np.zeros((n, D_T))

    for k in range(1, m + 1):
        x_next = (
            x
            + D_KAPPA * (D_THETA - x) * dt
            + sigma * np.sqrt(dt) * z[k - 1]
        )

        r_current = np.exp(x)
        r_next = np.exp(x_next)

        integrated_rate += (
            0.5 * (r_current + r_next) * dt
        )

        x = x_next

        if k % steps_per_year == 0:
            year_index = k // steps_per_year - 1
            annual_integrals[:, year_index] = integrated_rate

    # Discount every coupon and principal payment
    discount_factors = np.exp(-annual_integrals)

    Y = discount_factors @ D_CASHFLOWS

    return Y, annual_integrals


# Independent samples for coefficient estimation and pricing
d_z_pilot = d_generate_normals(
    D_N_PILOT, D_DT, D_SEED + 1
)

d_z_pricing = d_generate_normals(
    D_N, D_DT, D_SEED
)

print("Part D BK simulation ready.")
print(f"Pilot paths: {D_N_PILOT}")
print(f"Pricing paths: {D_N}")
```

**저장된 실행 결과**

```text
Part D BK simulation ready.
Pilot paths: 2000
Pricing paths: 20000
```

**Python Implementation**

```python
def d_integrated_rate_moments(sigma, dt=D_DT):
    """
    Exact first and second moments of trapezoidal
    integrated BK rates on the Euler time grid.

    The Euler log rates are jointly Gaussian.
    Their exponentials have known first and
    second cross moments.
    """
    M = round(D_T / dt)

    a = 1 - D_KAPPA * dt
    q2 = sigma**2 * dt

    indices = np.arange(M + 1)

    # Mean and variance of the Euler log rates
    mean_x = np.zeros(M + 1)
    var_x = np.zeros(M + 1)

    mean_x[0] = np.log(D_R0)

    for k in range(1, M + 1):
        mean_x[k] = (
            a * mean_x[k - 1]
            + D_KAPPA * D_THETA * dt
        )

        var_x[k] = a**2 * var_x[k - 1] + q2

    # Covariance of two Euler log-rate observations
    k_min = np.minimum.outer(indices, indices)
    k_diff = np.abs(np.subtract.outer(indices, indices))

    cov_x = var_x[k_min] * a**k_diff

    # First moment of the BK short rate
    mean_r = np.exp(mean_x + 0.5 * var_x)

    # Second cross moment E[r_k * r_l]
    second_r = np.exp(
        mean_x[:, None]
        + mean_x[None, :]
        + 0.5 * (var_x[:, None] + var_x[None, :])
        + cov_x
    )

    # Trapezoidal weights for annual integrated rates
    weights = np.zeros((D_T, M + 1))

    for year in range(1, D_T + 1):
        end = round(year / dt)

        weights[year - 1, 0] = dt / 2
        weights[year - 1, 1:end] = dt
        weights[year - 1, end] = dt / 2

    # E[I(t_i)]
    mean_I = weights @ mean_r

    # E[I(t_i)^2]
    mean_I2 = np.sum(
        (weights @ second_r) * weights,
        axis=1
    )

    return mean_I, mean_I2


def d_control_expectations(sigma, dt=D_DT):
    mean_I, mean_I2 = d_integrated_rate_moments(
        sigma, dt
    )

    # Controls based on maturity T
    mu_IT = mean_I[-1]
    mu_IT2 = mean_I2[-1]

    # Coupon-weighted controls
    mu_X1 = D_CASHFLOWS @ mean_I
    mu_X2 = D_CASHFLOWS @ mean_I2

    return {
        "IT": float(mu_IT),
        "IT2": float(mu_IT2),
        "X1": float(mu_X1),
        "X2": float(mu_X2)
    }


for sigma in D_SIGMAS:
    mu = d_control_expectations(sigma)

    print(f"\nBK sigma = {sigma}")
    for name, value in mu.items():
        print(f"E[{name}] = {value:.8f}")
```

**저장된 실행 결과**

```text

BK sigma = 0.2
E[IT] = 0.15166495
E[IT2] = 0.02318686
E[X1] = 0.16984816
E[X2] = 0.02522634

BK sigma = 0.5
E[IT] = 0.16073316
E[IT2] = 0.02719000
E[X1] = 0.17991247
E[X2] = 0.02957599
```

## Task 9. BK Path-Based Control Variates

BK 만기 적분금리 I(T), I(T)와 제곱항, 쿠폰 가중 적분금리 X₁, 두 쿠폰 가중 항 X₁·X₂의 네 구성을 비교한다. 여러 통제변량을 동시에 사용하면 계수 ν̂는 벡터가 되며, 독립 pilot 표본에 대한 최소제곱으로 추정한다.

**Python Implementation**

```python
def d_make_controls(annual_integrals):
    """
    Build four quantities from the same BK path.
    """
    I_T = annual_integrals[:, -1]
    I_T2 = I_T**2

    X1 = annual_integrals @ D_CASHFLOWS
    X2 = (annual_integrals**2) @ D_CASHFLOWS

    return {
        "IT": I_T,
        "IT2": I_T2,
        "X1": X1,
        "X2": X2
    }


def d_fit_cv(
    y_pilot,
    x_pilot,
    y_pricing,
    x_pricing,
    mu
):
    """
    Works with one or multiple control variates.

    Pilot sample: estimate nu.
    Pricing sample: estimate corrected BK price.
    """
    x_pilot = np.asarray(x_pilot)
    x_pricing = np.asarray(x_pricing)
    mu = np.atleast_1d(np.asarray(mu, dtype=float))

    if x_pilot.ndim == 1:
        x_pilot = x_pilot[:, None]

    if x_pricing.ndim == 1:
        x_pricing = x_pricing[:, None]

    # Center pilot observations for regression
    centered_x = x_pilot - np.mean(
        x_pilot, axis=0
    )

    centered_y = y_pilot - np.mean(y_pilot)

    # Fit optimal linear-control coefficients on pilot paths
    nu = np.linalg.lstsq(
        centered_x,
        centered_y,
        rcond=None
    )[0]

    # Apply the correction to independent pricing paths
    y_cv = (
        y_pricing
        - (x_pricing - mu) @ nu
    )

    plain_price = np.mean(y_pricing)
    plain_se = np.std(
        y_pricing, ddof=1
    ) / np.sqrt(len(y_pricing))

    cv_price = np.mean(y_cv)
    cv_se = np.std(
        y_cv, ddof=1
    ) / np.sqrt(len(y_cv))

    vrf = (
        np.var(y_pricing, ddof=1)
        / np.var(y_cv, ddof=1)
    )

    return {
        "plain_price": float(plain_price),
        "plain_se": float(plain_se),
        "price": float(cv_price),
        "se": float(cv_se),
        "nu": nu,
        "vrf": float(vrf)
    }


print("BK path-based CV functions ready.")
```

**저장된 실행 결과**

```text
BK path-based CV functions ready.
```

**Python Implementation**

```python
d_task9_results = {}

for sigma in D_SIGMAS:

    # Independent pilot and pricing simulations
    Y_pilot, I_pilot = d_simulate_bk(
        d_z_pilot, sigma
    )

    Y_pricing, I_pricing = d_simulate_bk(
        d_z_pricing, sigma
    )

    controls_pilot = d_make_controls(I_pilot)
    controls_pricing = d_make_controls(I_pricing)

    mu = d_control_expectations(sigma)

    configurations = {
        "I(T)": ["IT"],
        "I(T), I(T)^2": ["IT", "IT2"],
        "X1": ["X1"],
        "X1, X2": ["X1", "X2"]
    }

    d_task9_results[sigma] = {}

    print(f"\n=== BK sigma = {sigma} ===")

    for name, variables in configurations.items():

        X_pilot = np.column_stack([
            controls_pilot[v] for v in variables
        ])

        X_pricing = np.column_stack([
            controls_pricing[v] for v in variables
        ])

        mu_vector = np.array([
            mu[v] for v in variables
        ])

        result = d_fit_cv(
            Y_pilot,
            X_pilot,
            Y_pricing,
            X_pricing,
            mu_vector
        )

        d_task9_results[sigma][name] = result

        print(f"\nControl: {name}")
        print(f"CV price: {result['price']:.8f}")
        print(f"SE:       {result['se']:.8f}")
        print(f"nu:       {np.round(result['nu'], 6)}")
        print(f"VRF:      {result['vrf']:.4f}")

    plain = result

    print("\nPlain Monte Carlo")
    print(f"Price: {plain['plain_price']:.8f}")
    print(f"SE:    {plain['plain_se']:.8f}")
```

**저장된 실행 결과**

```text

=== BK sigma = 0.2 ===

Control: I(T)
CV price: 1.04215152
SE:       0.00000345
nu:       [-0.953406]
VRF:      710.8564

Control: I(T), I(T)^2
CV price: 1.04215154
SE:       0.00000333
nu:       [-1.098813  0.472347]
VRF:      763.3863

Control: X1
CV price: 1.04215044
SE:       0.00000093
nu:       [-0.860525]
VRF:      9873.6842

Control: X1, X2
CV price: 1.04215038
SE:       0.00000001
nu:       [-0.989913  0.430068]
VRF:      57026630.2968

Plain Monte Carlo
Price: 1.04212403
SE:    0.00009197

=== BK sigma = 0.5 ===

Control: I(T)
CV price: 1.03404968
SE:       0.00001145
nu:       [-0.932645]
VRF:      456.1278

Control: I(T), I(T)^2
CV price: 1.03404802
SE:       0.00000892
nu:       [-1.094778  0.459406]
VRF:      751.8014

Control: X1
CV price: 1.03405042
SE:       0.00000710
nu:       [-0.842361]
VRF:      1186.3044

Control: X1, X2
CV price: 1.03404842
SE:       0.00000019
nu:       [-0.985351  0.414391]
VRF:      1722574.8421

Plain Monte Carlo
Price: 1.03400908
SE:    0.00024446
```

### Task 9. Results and Interpretation

이번 실험에서는 쿠폰 지급 시점을 반영한 X₁, X₂ 구성이 특히 큰 분산 감소 효과를 보였다. 그러나 VRF가 수백만 이상으로 커지는 경우 **기대값의 작은 불일치가 추정 SE보다 훨씬 커질 수 있으므로**, 높은 VRF만으로 추정량의 정확성을 판단하면 안 된다.

| BK σ | Control | SE | VRF |
|---|---|---:|---:|
| 0.2 | I(T) | 0.00000345 | 710.86 |
| 0.2 | X₁, X₂ | 0.00000001 | 57,026,630.30 |
| 0.5 | I(T) | 0.00001145 | 456.13 |
| 0.5 | X₁, X₂ | 0.00000019 | 1,722,574.84 |

## Task 10. Comparison with Part C

Part C의 기본 Vasicek 통제변량과 Part D의 BK 경로 기반 통제변량을 동일한 BK 경로·난수·표본수 조건에서 비교한다. Part C의 Baseline 설정을 기준으로 한 비교임을 분명히 한다.

**Python Implementation**

```python
# Baseline Vasicek parameters from Part C
D_ALPHA = 0.8
D_BETA = 0.03
D_ETA = 0.01


def d_simulate_vasicek_controls(
    z,
    alpha=D_ALPHA,
    beta=D_BETA,
    eta=D_ETA,
    dt=D_DT
):
    m, n = z.shape
    steps_per_year = m // D_T

    rates = np.full(n, D_R0)
    integrated = np.zeros(n)

    X_coupon = np.zeros(n)

    for k in range(1, m + 1):
        next_rates = (
            rates
            + alpha * (beta - rates) * dt
            + eta * np.sqrt(dt) * z[k - 1]
        )

        integrated += (
            0.5 * (rates + next_rates) * dt
        )

        rates = next_rates

        if k % steps_per_year == 0:
            year_index = k // steps_per_year - 1
            cashflow = D_CASHFLOWS[year_index]

            X_coupon += (
                cashflow * np.exp(-integrated)
            )

    X_zcb = np.exp(-integrated)

    return X_zcb, X_coupon


def d_vasicek_expectations(
    alpha=D_ALPHA,
    beta=D_BETA,
    eta=D_ETA,
    dt=D_DT
):
    """
    Exact expectations of Euler-trapezoidal
    Vasicek controls.
    """
    m = round(D_T / dt)
    steps_per_year = m // D_T

    mean = np.array([D_R0, 0.0])
    cov = np.zeros((2, 2))

    a = 1 - alpha * dt
    b = alpha * beta * dt
    q = eta * np.sqrt(dt)

    A = np.array([
        [a, 0.0],
        [0.5 * dt * (1 + a), 1.0]
    ])

    shift = np.array([b, 0.5 * dt * b])
    noise = np.array([q, 0.5 * dt * q])

    coupon_expectation = 0.0

    for k in range(1, m + 1):
        mean = A @ mean + shift
        cov = A @ cov @ A.T + np.outer(noise, noise)

        if k % steps_per_year == 0:
            year_index = k // steps_per_year - 1

            expected_discount = np.exp(
                -mean[1] + 0.5 * cov[1, 1]
            )

            coupon_expectation += (
                D_CASHFLOWS[year_index]
                * expected_discount
            )

    zcb_expectation = np.exp(
        -mean[1] + 0.5 * cov[1, 1]
    )

    return zcb_expectation, coupon_expectation


print("Part C comparison functions ready.")
```

**저장된 실행 결과**

```text
Part C comparison functions ready.
```

**Python Implementation**

```python
d_task10_results = {}

mu_vas_zcb, mu_vas_coupon = d_vasicek_expectations()

for sigma in D_SIGMAS:

    # The same BK paths are reused in Parts C and D
    Y_pilot, _ = d_simulate_bk(
        d_z_pilot, sigma
    )

    Y_pricing, _ = d_simulate_bk(
        d_z_pricing, sigma
    )

    # Vasicek controls with the same Gaussian increments
    Xz_pilot, Xc_pilot = d_simulate_vasicek_controls(
        d_z_pilot
    )

    Xz_pricing, Xc_pricing = d_simulate_vasicek_controls(
        d_z_pricing
    )

    c_zcb = d_fit_cv(
        Y_pilot,
        Xz_pilot,
        Y_pricing,
        Xz_pricing,
        [mu_vas_zcb]
    )

    c_coupon = d_fit_cv(
        Y_pilot,
        Xc_pilot,
        Y_pricing,
        Xc_pricing,
        [mu_vas_coupon]
    )

    d_task10_results[sigma] = {
        "Part C: Vasicek ZCB": c_zcb,
        "Part C: Vasicek Coupon": c_coupon,
        "Part D: I(T)": d_task9_results[sigma]["I(T)"],
        "Part D: I(T), I(T)^2": d_task9_results[sigma]["I(T), I(T)^2"],
        "Part D: X1": d_task9_results[sigma]["X1"],
        "Part D: X1, X2": d_task9_results[sigma]["X1, X2"]
    }

    print(f"\n=== BK sigma = {sigma} ===")

    for name, result in d_task10_results[sigma].items():
        print(
            f"{name:26s} "
            f"SE={result['se']:.8f} "
            f"VRF={result['vrf']:.2f} "
            f"nu={np.round(result['nu'], 5)}"
        )
```

**저장된 실행 결과**

```text

=== BK sigma = 0.2 ===
Part C: Vasicek ZCB        SE=0.00000884 VRF=108.36 nu=[0.67256]
Part C: Vasicek Coupon     SE=0.00000816 VRF=127.16 nu=[0.60476]
Part D: I(T)               SE=0.00000345 VRF=710.86 nu=[-0.95341]
Part D: I(T), I(T)^2       SE=0.00000333 VRF=763.39 nu=[-1.09881  0.47235]
Part D: X1                 SE=0.00000093 VRF=9873.68 nu=[-0.86052]
Part D: X1, X2             SE=0.00000001 VRF=57026630.30 nu=[-0.98991  0.43007]

=== BK sigma = 0.5 ===
Part C: Vasicek ZCB        SE=0.00004926 VRF=24.63 nu=[1.76228]
Part C: Vasicek Coupon     SE=0.00004851 VRF=25.39 nu=[1.58421]
Part D: I(T)               SE=0.00001145 VRF=456.13 nu=[-0.93265]
Part D: I(T), I(T)^2       SE=0.00000892 VRF=751.80 nu=[-1.09478  0.45941]
Part D: X1                 SE=0.00000710 VRF=1186.30 nu=[-0.84236]
Part D: X1, X2             SE=0.00000019 VRF=1722574.84 nu=[-0.98535  0.41439]
```

![Vasicek 및 BK 통제변량의 VRF 비교](/assets/images/bk_vrf_comparison.png)
*Figure. Vasicek 및 BK 통제변량의 VRF 비교.*

### Task 10. Results and Interpretation

이번 계산에서는 BK의 쿠폰 가중 적분금리 기반 통제변량이 더 높은 VRF를 기록했다. 이는 할인계수의 1차·2차 Taylor 근사와 관련되지만, 모든 변동성이나 파라미터에서 같은 우위가 보장되는 것은 아니다.

### Part D Validation

Part D의 이산 모멘트가 기초 부등식을 만족하는지, 독립 표본평균이 계산된 기대값과 일치하는지, CV 추정 결과가 유한한지 검사한다. 또한 σ=0일 때 결정론적인 채권 가격으로 돌아가는지 확인한다.

**Python Implementation**

```python
# Part D validation

for sigma in D_SIGMAS:

    Y, I = d_simulate_bk(
        d_z_pricing, sigma
    )

    controls = d_make_controls(I)
    mu = d_control_expectations(sigma)

    # Test 1: Finite, positive bond payoffs
    assert np.all(np.isfinite(Y))
    assert np.all(Y > 0)

    # Test 2: Moment consistency
    mean_I, mean_I2 = d_integrated_rate_moments(
        sigma
    )

    assert np.all(
        mean_I2 >= mean_I**2 - 1e-12
    )

    # Test 3: Control expectations versus MC means
    for name in ["IT", "IT2", "X1", "X2"]:
        X = controls[name]

        se_mean = np.std(
            X, ddof=1
        ) / np.sqrt(D_N)

        difference = abs(
            np.mean(X) - mu[name]
        )

        assert difference < 4 * se_mean, (
            f"Control expectation check failed: "
            f"sigma={sigma}, control={name}"
        )

    # Test 4: CV outputs are finite
    for name, result in d_task9_results[sigma].items():

        assert np.isfinite(result["price"])
        assert np.isfinite(result["se"])
        assert np.all(np.isfinite(result["nu"]))
        assert np.isfinite(result["vrf"])

        assert result["se"] > 0
        assert result["vrf"] > 0

    print(f"BK sigma={sigma}: validation passed")


# Test 5: Deterministic zero-volatility case
z_test = d_generate_normals(
    n=100,
    dt=D_DT,
    seed=9876
)

Y_zero, I_zero = d_simulate_bk(
    z_test, sigma=0.0
)

mu_zero = d_control_expectations(sigma=0.0)
controls_zero = d_make_controls(I_zero)

years = np.arange(1, D_T + 1)

expected_bond_price = np.sum(
    D_CASHFLOWS * np.exp(-D_R0 * years)
)

assert np.allclose(
    Y_zero,
    expected_bond_price,
    atol=1e-12
)

for name in ["IT", "IT2", "X1", "X2"]:
    assert np.allclose(
        controls_zero[name],
        mu_zero[name],
        atol=1e-12
    )

print("Zero-volatility validation passed.")
print("All Part D validation tests passed.")
```

**저장된 실행 결과**

```text
BK sigma=0.2: validation passed
BK sigma=0.5: validation passed
Zero-volatility validation passed.
All Part D validation tests passed.
```

### Part D Summary

Part D에서는 BK 경로를 직접 사용하는 통제변량을 구성하고 네 조합의 성능을 비교했다. 특히 낮은 표준오차와 **편향이 없다는 사실은 별개의 검증 대상**이라는 점을 확인했다.

# Part E: Testing the Known Expectations & Summary

## Task 11. Independent Validation of Control Expectations

통제변량의 기대값이 정확한지를 독립 Monte Carlo 표본 30,000개로 검증한다. 기본 Vasicek 통제변량뿐 아니라 **Matched와 Grid best까지 포함**하여 Part C에서 실제 사용한 모든 설정을 확인한다. Part D의 I(T), I(T)², X₁, X₂도 두 BK 변동성별로 검증한다. 허용 기준은 **3 SE**로 설정한다.

**Python Implementation**

```python
import numpy as np
import matplotlib.pyplot as plt
import time

# Part E simulation settings
E_SIGMAS = (0.2, 0.5)
E_DT = 1 / 52
E_N_TEST = 30000
E_SEED = 314159

E_ALPHA = 0.8
E_BETA = 0.03
E_ETA = 0.01

# Independent validation sample
z_test = d_generate_normals(
    E_N_TEST, E_DT, E_SEED
)

# Vasicek control expectations
mu_vas_z, mu_vas_c = d_vasicek_expectations(
    E_ALPHA, E_BETA, E_ETA, E_DT
)

# Simulate Vasicek controls
Xz, Xc = d_simulate_vasicek_controls(
    z_test, E_ALPHA, E_BETA, E_ETA, E_DT
)


def e_expectation_check(name, values, expected_mu):
    mean_mc = np.mean(values)
    se_mc = np.std(values, ddof=1) / np.sqrt(len(values))

    difference = mean_mc - expected_mu
    tolerance = 3.0 * se_mc

    passed = abs(difference) <= tolerance

    z_score = difference / se_mc

    # Executable validation
    assert passed, (
        f"{name}: expectation check failed. "
        f"Difference={difference:.8f}, "
        f"Tolerance={tolerance:.8f}"
    )

    return {
        "control": name,
        "expected": float(expected_mu),
        "mc_mean": float(mean_mc),
        "se": float(se_mc),
        "difference": float(difference),
        "tolerance": float(tolerance),
        "z_score": float(z_score),
        "passed": bool(passed)
    }


e_task11_results = []

# Part C: Vasicek controls
e_task11_results.append(
    e_expectation_check(
        "Vasicek ZCB", Xz, mu_vas_z
    )
)

e_task11_results.append(
    e_expectation_check(
        "Vasicek coupon", Xc, mu_vas_c
    )
)


# Part D: BK path-based controls
for sigma in E_SIGMAS:
    Y_test, I_test = d_simulate_bk(
        z_test, sigma, E_DT
    )

    controls = d_make_controls(I_test)
    expectations = d_control_expectations(sigma, E_DT)

    for name in ["IT", "IT2", "X1", "X2"]:
        result = e_expectation_check(
            f"BK sigma={sigma}, {name}",
            controls[name],
            expectations[name]
        )

        e_task11_results.append(result)


# Print validation table
print("=== Task 11: Control Expectation Validation ===")

print(
    f"{'Control':27s} "
    f"{'Expected':>11s} "
    f"{'Observed':>11s} "
    f"{'SE':>11s} "
    f"{'Tolerance':>11s} "
    f"{'Z-score':>9s} "
    f"{'Pass':>6s}"
)

for r in e_task11_results:
    print(
        f"{r['control']:27s} "
        f"{r['expected']:11.7f} "
        f"{r['mc_mean']:11.7f} "
        f"{r['se']:11.7f} "
        f"{r['tolerance']:11.7f} "
        f"{r['z_score']:9.3f} "
        f"{str(r['passed']):>6s}"
    )

print("\nAll Task 11 expectation checks passed.")
```

**저장된 실행 결과**

```text
=== Task 11: Control Expectation Validation ===
Control                        Expected    Observed          SE   Tolerance   Z-score   Pass
Vasicek ZCB                   0.8609221   0.8610819   0.0001099   0.0003298     1.454   True
Vasicek coupon                1.0438936   1.0440688   0.0001226   0.0003679     1.428   True
BK sigma=0.2, IT              0.1516649   0.1515368   0.0000778   0.0002334    -1.648   True
BK sigma=0.2, IT2             0.0231869   0.0231450   0.0000240   0.0000719    -1.749   True
BK sigma=0.2, X1              0.1698482   0.1697081   0.0000864   0.0002593    -1.621   True
BK sigma=0.2, X2              0.0252263   0.0251813   0.0000260   0.0000780    -1.733   True
BK sigma=0.5, IT              0.1607332   0.1603438   0.0002108   0.0006323    -1.847   True
BK sigma=0.5, IT2             0.0271900   0.0270429   0.0000748   0.0002244    -1.967   True
BK sigma=0.5, X1              0.1799125   0.1794867   0.0002340   0.0007021    -1.819   True
BK sigma=0.5, X2              0.0295760   0.0294177   0.0000811   0.0002433    -1.951   True

All Task 11 expectation checks passed.
```

**Python Implementation**

```python
# =========================================================
# TASK 11 - ADDITIONAL VASICEK EXPECTATION VALIDATION
# Matched and Grid Best parameter configurations
# =========================================================

e_task11_extra_results = []

for sigma in E_SIGMAS:

    for control_name in ["zcb", "coupon"]:

        # Retrieve exactly the parameters used in Task 7
        records = [
            r for r in task7_results
            if r["sigma"] == sigma
            and r["control"] == control_name
            and r["choice"] in ["Matched", "Grid best"]
        ]

        assert len(records) == 2, (
            f"Expected Matched and Grid best records "
            f"for sigma={sigma}, control={control_name}"
        )

        for record in records:

            alpha = record["alpha"]
            beta = record["beta"]
            eta = record["eta"]

            # Exact expectation under Euler-trapezoidal
            # Vasicek discretization
            mu_z, mu_c = d_vasicek_expectations(
                alpha, beta, eta, E_DT
            )

            # Independent MC sample from Task 11
            Xz_test, Xc_test = d_simulate_vasicek_controls(
                z_test, alpha, beta, eta, E_DT
            )

            if control_name == "zcb":
                values = Xz_test
                expected_mu = mu_z
            else:
                values = Xc_test
                expected_mu = mu_c

            label = (
                f"Vasicek {control_name}, "
                f"sigma={sigma}, {record['choice']}"
            )

            result = e_expectation_check(
                label,
                values,
                expected_mu
            )

            e_task11_extra_results.append(result)


# Verify that all 8 additional configurations
# have been checked.
assert len(e_task11_extra_results) == 8

assert all(
    result["passed"]
    for result in e_task11_extra_results
)

print("=== Task 11: Additional Vasicek Checks ===")

print(
    f"{'Control':42s} "
    f"{'Expected':>12s} "
    f"{'Observed':>12s} "
    f"{'Z-score':>10s} "
    f"{'Pass':>6s}"
)

for r in e_task11_extra_results:
    print(
        f"{r['control']:42s} "
        f"{r['expected']:12.8f} "
        f"{r['mc_mean']:12.8f} "
        f"{r['z_score']:10.3f} "
        f"{str(r['passed']):>6s}"
    )

print(
    "\nAdditional checks passed: "
    f"{sum(r['passed'] for r in e_task11_extra_results)}/8"
)

print(
    "Total Task 11 checks: "
    f"{len(e_task11_results) + len(e_task11_extra_results)}"
)
```

**저장된 실행 결과**

```text
=== Task 11: Additional Vasicek Checks ===
Control                                        Expected     Observed    Z-score   Pass
Vasicek zcb, sigma=0.2, Matched              0.86078504   0.86088177      1.467   True
Vasicek zcb, sigma=0.2, Grid best            0.86071889   0.86075523      1.463   True
Vasicek coupon, sigma=0.2, Matched           1.04374331   1.04384930      1.441   True
Vasicek coupon, sigma=0.2, Grid best         1.04366754   1.04370208      1.457   True
Vasicek zcb, sigma=0.5, Matched              0.86118973   0.86142681      1.437   True
Vasicek zcb, sigma=0.5, Grid best            0.86071889   0.86075523      1.463   True
Vasicek coupon, sigma=0.5, Matched           1.04418733   1.04444725      1.412   True
Vasicek coupon, sigma=0.5, Grid best         1.04366754   1.04370208      1.457   True

Additional checks passed: 8/8
Total Task 11 checks: 18
```

### Task 11. Results and Interpretation

추가 Vasicek 설정 8개도 독립 표본과 비교한 결과 모두 3 SE 기준을 통과했다. 기본 Vasicek 2개와 BK 통제변량 8개를 합하면 **총 18개 검사에서 18개 모두 통과**했다. 단, 통계적 허용 범위 통과는 수학적인 무오류 증명을 의미하지 않는다.

| Validation category | Checks | Passed |
|---|---:|---:|
| Vasicek Baseline | 2 | 2 |
| Vasicek Matched | 4 | 4 |
| Vasicek Grid best | 4 | 4 |
| BK path-based | 8 | 8 |
| **Total** | **18** | **18** |

## Task 12. Discretization Effects on Control Expectations

시간 간격 1/52와 1/252에서 계산한 기대값을 비교한다. Vasicek은 연속시간 가격과 이산시간 기대값을 대조하고, BK는 두 시간격자에서 통제변량 기대값이 어떻게 달라지는지 확인한다. **사다리꼴 적분과 왼쪽 직사각형 적분의 차이도 Vasicek과 BK 모두 비교**한다.

**Python Implementation**

```python
# Task 12: time-step effects

E_STEPS = (1 / 52, 1 / 252)

def e_continuous_vasicek_zcb(t, alpha, beta, eta, r0):
    B = -np.expm1(-alpha * t) / alpha
    A = (
        (beta - eta**2 / (2 * alpha**2)) * (B - t)
        - eta**2 * B**2 / (4 * alpha)
    )
    return np.exp(A - B * r0)


vas_cont_zcb = e_continuous_vasicek_zcb(
    D_T, E_ALPHA, E_BETA, E_ETA, D_R0
)

vas_cont_coupon = sum(
    D_CASHFLOWS[year - 1]
    * e_continuous_vasicek_zcb(
        year, E_ALPHA, E_BETA, E_ETA, D_R0
    )
    for year in range(1, D_T + 1)
)

e_task12_expectations = {}

print("=== Vasicek: discrete vs continuous ===")

for dt in E_STEPS:
    mu_z, mu_c = d_vasicek_expectations(
        E_ALPHA, E_BETA, E_ETA, dt
    )

    e_task12_expectations[("Vasicek", dt)] = (mu_z, mu_c)

    print(f"\ndt = {dt:.8f}")
    print(f"ZCB discrete:    {mu_z:.10f}")
    print(f"ZCB continuous:  {vas_cont_zcb:.10f}")
    print(f"ZCB difference:  {mu_z - vas_cont_zcb:+.10f}")
    print(f"Coupon discrete:   {mu_c:.10f}")
    print(f"Coupon continuous: {vas_cont_coupon:.10f}")
    print(f"Coupon difference: {mu_c - vas_cont_coupon:+.10f}")


print("\n=== BK: discrete control expectations ===")

for sigma in E_SIGMAS:
    print(f"\nBK sigma = {sigma}")

    for dt in E_STEPS:
        mu = d_control_expectations(sigma, dt)
        e_task12_expectations[(sigma, dt)] = mu

        print(f"dt={dt:.8f}: {mu}")

    coarse = e_task12_expectations[(sigma, E_STEPS[0])]
    fine = e_task12_expectations[(sigma, E_STEPS[1])]

    print("Differences (coarse - fine):")
    for name in ["IT", "IT2", "X1", "X2"]:
        print(f"{name}: {coarse[name] - fine[name]:+.10f}")
```

**저장된 실행 결과**

```text
=== Vasicek: discrete vs continuous ===

dt = 0.01923077
ZCB discrete:    0.8609220553
ZCB continuous:  0.8609212015
ZCB difference:  +0.0000008538
Coupon discrete:   1.0438936455
Coupon continuous: 1.0438926820
Coupon difference: +0.0000009634

dt = 0.00396825
ZCB discrete:    0.8609213775
ZCB continuous:  0.8609212015
ZCB difference:  +0.0000001760
Coupon discrete:   1.0438928806
Coupon continuous: 1.0438926820
Coupon difference: +0.0000001986

=== BK: discrete control expectations ===

BK sigma = 0.2
dt=0.01923077: {'IT': 0.1516649462365563, 'IT2': 0.023186855000636046, 'X1': 0.16984816437304046, 'X2': 0.025226337976935968}
dt=0.00396825: {'IT': 0.15165327717174187, 'IT2': 0.02318269486844777, 'X1': 0.16983512578044205, 'X2': 0.025221793731939134}
Differences (coarse - fine):
IT: +0.0000116691
IT2: +0.0000041601
X1: +0.0000130386
X2: +0.0000045442

BK sigma = 0.5
dt=0.01923077: {'IT': 0.16073316165592397, 'IT2': 0.027189996095222484, 'X1': 0.1799124735617631, 'X2': 0.029575994136731204}
dt=0.00396825: {'IT': 0.16065572462433808, 'IT2': 0.02715914601134587, 'X1': 0.17982597968038702, 'X2': 0.029542299413942}
Differences (coarse - fine):
IT: +0.0000774370
IT2: +0.0000308501
X1: +0.0000864939
X2: +0.0000336947
```

**Python Implementation**

```python
def e_vasicek_expectations_with_rule(
    alpha, beta, eta, dt, rule="trapezoid"
):
    m = round(D_T / dt)
    steps_per_year = m // D_T

    mean = np.array([D_R0, 0.0])
    cov = np.zeros((2, 2))

    a = 1 - alpha * dt
    b = alpha * beta * dt
    q = eta * np.sqrt(dt)

    if rule == "trapezoid":
        w_old = dt * (1 + a) / 2
        w_noise = dt * q / 2
        w_shift = dt * b / 2
    elif rule == "left":
        w_old = dt
        w_noise = 0.0
        w_shift = 0.0
    else:
        raise ValueError("Unknown integration rule")

    A = np.array([
        [a, 0.0],
        [w_old, 1.0]
    ])

    shift = np.array([b, w_shift])
    noise = np.array([q, w_noise])

    coupon_mu = 0.0

    for k in range(1, m + 1):
        mean = A @ mean + shift
        cov = A @ cov @ A.T + np.outer(noise, noise)

        if k % steps_per_year == 0:
            year_index = k // steps_per_year - 1

            expected_discount = np.exp(
                -mean[1] + 0.5 * cov[1, 1]
            )

            coupon_mu += (
                D_CASHFLOWS[year_index] * expected_discount
            )

    zcb_mu = np.exp(
        -mean[1] + 0.5 * cov[1, 1]
    )

    return zcb_mu, coupon_mu


for dt in E_STEPS:
    trap = e_vasicek_expectations_with_rule(
        E_ALPHA, E_BETA, E_ETA, dt, "trapezoid"
    )

    left = e_vasicek_expectations_with_rule(
        E_ALPHA, E_BETA, E_ETA, dt, "left"
    )

    print(f"\nTime step = {dt:.8f}")

    for name, trap_mu, left_mu in zip(
        ["ZCB", "Coupon"], trap, left
    ):
        print(
            f"{name}: trapezoid={trap_mu:.10f}, "
            f"left={left_mu:.10f}, "
            f"difference={left_mu - trap_mu:+.10f}"
        )
```

**저장된 실행 결과**

```text

Time step = 0.01923077
ZCB: trapezoid=0.8609220553, left=0.8609214338, difference=-0.0000006215
Coupon: trapezoid=1.0438936455, left=1.0438929251, difference=-0.0000007204

Time step = 0.00396825
ZCB: trapezoid=0.8609213775, left=0.8609212490, difference=-0.0000001285
Coupon: trapezoid=1.0438928806, left=1.0438927316, difference=-0.0000001490
```

**Python Implementation**

```python
# BK trapezoidal vs left-point integration
# Compare exact expectations of discretized controls.

def e_bk_expectations_with_rule(sigma, dt, rule):
    M = round(D_T / dt)
    a = 1 - D_KAPPA * dt

    mean_x = np.zeros(M + 1)
    var_x = np.zeros(M + 1)

    mean_x[0] = np.log(D_R0)

    for k in range(1, M + 1):
        mean_x[k] = (
            a * mean_x[k - 1]
            + D_KAPPA * D_THETA * dt
        )
        var_x[k] = (
            a**2 * var_x[k - 1]
            + sigma**2 * dt
        )

    indices = np.arange(M + 1)
    k_min = np.minimum.outer(indices, indices)
    k_diff = np.abs(np.subtract.outer(indices, indices))

    cov_x = var_x[k_min] * a**k_diff

    # Lognormal first and second moments
    mean_r = np.exp(mean_x + 0.5 * var_x)

    second_r = np.exp(
        mean_x[:, None]
        + mean_x[None, :]
        + 0.5 * (
            var_x[:, None] + var_x[None, :]
        )
        + cov_x
    )

    # Integration weights at annual payment dates
    weights = np.zeros((D_T, M + 1))

    for year in range(1, D_T + 1):
        end = round(year / dt)

        if rule == "trapezoid":
            weights[year - 1, 0] = dt / 2
            weights[year - 1, 1:end] = dt
            weights[year - 1, end] = dt / 2

        elif rule == "left":
            weights[year - 1, :end] = dt

        else:
            raise ValueError("Unknown integration rule")

    mean_I = weights @ mean_r

    mean_I2 = np.sum(
        (weights @ second_r) * weights,
        axis=1
    )

    return {
        "IT": float(mean_I[-1]),
        "IT2": float(mean_I2[-1]),
        "X1": float(D_CASHFLOWS @ mean_I),
        "X2": float(D_CASHFLOWS @ mean_I2)
    }


print("=== Task 12: BK Integration Rule Comparison ===")

for sigma in E_SIGMAS:
    for dt in E_STEPS:

        trap = e_bk_expectations_with_rule(
            sigma, dt, "trapezoid"
        )

        left = e_bk_expectations_with_rule(
            sigma, dt, "left"
        )

        # Cross-check against the original implementation
        original = d_control_expectations(sigma, dt)

        for name in ["IT", "IT2", "X1", "X2"]:
            assert np.isclose(
                trap[name],
                original[name],
                rtol=1e-10,
                atol=1e-12
            )

        print(f"\nBK sigma = {sigma}, dt = {dt:.8f}")
        print(
            f"{'Control':8s} "
            f"{'Trapezoid':>14s} "
            f"{'Left-point':>14s} "
            f"{'Difference':>14s}"
        )

        for name in ["IT", "IT2", "X1", "X2"]:
            difference = left[name] - trap[name]

            print(
                f"{name:8s} "
                f"{trap[name]:14.10f} "
                f"{left[name]:14.10f} "
                f"{difference:+14.10f}"
            )

print("\nAll BK integration-rule checks passed.")
```

**저장된 실행 결과**

```text
=== Task 12: BK Integration Rule Comparison ===

BK sigma = 0.2, dt = 0.01923077
Control       Trapezoid     Left-point     Difference
IT         0.1516649462   0.1516612907  -0.0000036555
IT2        0.0231868550   0.0231852098  -0.0000016452
X1         0.1698481644   0.1698438141  -0.0000043503
X2         0.0252263380   0.0252244812  -0.0000018567

BK sigma = 0.2, dt = 0.00396825
Control       Trapezoid     Left-point     Difference
IT         0.1516532772   0.1516525275  -0.0000007497
IT2        0.0231826949   0.0231823566  -0.0000003383
X1         0.1698351258   0.1698342337  -0.0000008920
X2         0.0252217937   0.0252214119  -0.0000003818

BK sigma = 0.5, dt = 0.01923077
Control       Trapezoid     Left-point     Difference
IT         0.1607331617   0.1607095407  -0.0000236210
IT2        0.0271899961   0.0271784454  -0.0000115507
X1         0.1799124736   0.1798843700  -0.0000281036
X2         0.0295759941   0.0295629588  -0.0000130353

BK sigma = 0.5, dt = 0.00396825
Control       Trapezoid     Left-point     Difference
IT         0.1606557246   0.1606508815  -0.0000048431
IT2        0.0271591460   0.0271567722  -0.0000023738
X1         0.1798259797   0.1798202181  -0.0000057616
X2         0.0295422994   0.0295396205  -0.0000026789

All BK integration-rule checks passed.
```

**Python Implementation**

```python
# Task 12: expectation mismatch relative to standard error

e_task12_bias_results = []

for sigma in E_SIGMAS:
    Y_pilot, I_pilot = d_simulate_bk(
        d_z_pilot, sigma, D_DT
    )

    Y_price, I_price = d_simulate_bk(
        d_z_pricing, sigma, D_DT
    )

    cp = d_make_controls(I_pilot)
    cf = d_make_controls(I_price)

    mu_correct = d_control_expectations(sigma, D_DT)
    mu_fine = d_control_expectations(sigma, 1 / 252)

    for control_name, keys in {
        "I(T)": ["IT"],
        "I(T), I(T)^2": ["IT", "IT2"],
        "X1": ["X1"],
        "X1, X2": ["X1", "X2"]
    }.items():

        xp = np.column_stack([cp[k] for k in keys])
        xf = np.column_stack([cf[k] for k in keys])

        correct_vector = np.array([mu_correct[k] for k in keys])
        wrong_vector = np.array([mu_fine[k] for k in keys])

        result = d_fit_cv(
            Y_pilot, xp, Y_price, xf, correct_vector
        )

        shift = float(
            result["nu"] @ (wrong_vector - correct_vector)
        )

        ratio = abs(shift) / result["se"]

        e_task12_bias_results.append(
            (sigma, control_name, shift, ratio)
        )


print("=== Task 12: Expectation mismatch ===")
print("Mismatch: fine-grid mu used for coarse-grid controls")

for sigma, name, shift, ratio in e_task12_bias_results:
    print(
        f"sigma={sigma}, {name:14s}: "
        f"price shift={shift:+.9f}, "
        f"|shift|/SE={ratio:.3f}"
    )


# Repeat for both Vasicek controls
Xzp, Xcp = d_simulate_vasicek_controls(d_z_pilot)
Xzf, Xcf = d_simulate_vasicek_controls(d_z_pricing)

mu_vas_coarse = d_vasicek_expectations(
    E_ALPHA, E_BETA, E_ETA, D_DT
)
mu_vas_fine = d_vasicek_expectations(
    E_ALPHA, E_BETA, E_ETA, 1 / 252
)

for sigma in E_SIGMAS:
    Yp, _ = d_simulate_bk(d_z_pilot, sigma)
    Yf, _ = d_simulate_bk(d_z_pricing, sigma)

    for name, xp, xf, correct_mu, wrong_mu in zip(
        ["Vasicek ZCB", "Vasicek coupon"],
        [Xzp, Xcp],
        [Xzf, Xcf],
        mu_vas_coarse,
        mu_vas_fine
    ):
        result = d_fit_cv(
            Yp, xp, Yf, xf, [correct_mu]
        )

        shift = float(result["nu"][0] * (wrong_mu - correct_mu))

        print(
            f"sigma={sigma}, {name}: "
            f"price shift={shift:+.9f}, "
            f"|shift|/SE={abs(shift)/result['se']:.3f}"
        )
```

**저장된 실행 결과**

```text
=== Task 12: Expectation mismatch ===
Mismatch: fine-grid mu used for coarse-grid controls
sigma=0.2, I(T)          : price shift=+0.000011125, |shift|/SE=3.225
sigma=0.2, I(T), I(T)^2  : price shift=+0.000010857, |shift|/SE=3.262
sigma=0.2, X1            : price shift=+0.000011220, |shift|/SE=12.123
sigma=0.2, X1, X2        : price shift=+0.000010953, |shift|/SE=899.334
sigma=0.5, I(T)          : price shift=+0.000072221, |shift|/SE=6.310
sigma=0.5, I(T), I(T)^2  : price shift=+0.000070604, |shift|/SE=7.919
sigma=0.5, X1            : price shift=+0.000072859, |shift|/SE=10.265
sigma=0.5, X1, X2        : price shift=+0.000071264, |shift|/SE=382.602
sigma=0.2, Vasicek ZCB: price shift=-0.000000456, |shift|/SE=0.052
sigma=0.2, Vasicek coupon: price shift=-0.000000463, |shift|/SE=0.057
sigma=0.5, Vasicek ZCB: price shift=-0.000001195, |shift|/SE=0.024
sigma=0.5, Vasicek coupon: price shift=-0.000001212, |shift|/SE=0.025
```

### Task 12. Results and Interpretation

시간 간격을 줄이면 이 실험에서 이산화로 인한 기대값 차이가 작아졌다. BK에서는 적분 규칙을 왼쪽 직사각형으로 바꾸면 사다리꼴 기대값과 차이가 나고, 높은 변동성에서 차이가 더 컸다. 특히 [X₁, X₂]처럼 SE가 매우 작은 추정량에서는 작은 기대값 불일치도 수백 SE 규모의 가격 이동을 만들었다.

| BK σ | dt | Left-point − Trapezoid: E[I(T)] | Left-point − Trapezoid: E[X₁] |
|---|---|---:|---:|
| 0.2 | 1/52 | −0.0000036555 | −0.0000043503 |
| 0.2 | 1/252 | −0.0000007497 | −0.0000008920 |
| 0.5 | 1/52 | −0.0000236210 | −0.0000281036 |
| 0.5 | 1/252 | −0.0000048431 | −0.0000057616 |

## Task 13. Incorrect Control Expectations

BK 통제변량 X₁의 기대값에 의도적으로 δ=0.005를 더한다. 잘못된 기대값은 모든 보정 가격에 동일한 상수를 더하므로 **가격 추정값을 이동시키지만 표본분산과 SE를 바꾸지 않는다.** 95% 신뢰구간의 포함률도 독립 반복 실험으로 비교한다.

**Python Implementation**

```python
E_DELTA = 0.005

e_task13_results = {}

for sigma in E_SIGMAS:
    # BK pilot and pricing paths
    Yp, Ip = d_simulate_bk(d_z_pilot, sigma)
    Y, I = d_simulate_bk(d_z_pricing, sigma)

    # Use the BK control X1 for this illustration
    Xp = d_make_controls(Ip)["X1"]
    X = d_make_controls(I)["X1"]

    correct_mu = d_control_expectations(sigma)["X1"]
    incorrect_mu = correct_mu + E_DELTA

    correct = d_fit_cv(
        Yp, Xp, Y, X, [correct_mu]
    )

    incorrect = d_fit_cv(
        Yp, Xp, Y, X, [incorrect_mu]
    )

    e_task13_results[sigma] = {
        "correct": correct,
        "incorrect": incorrect
    }

    actual_shift = incorrect["price"] - correct["price"]
    predicted_shift = correct["nu"][0] * E_DELTA

    print(f"\n=== BK sigma = {sigma} ===")
    print(f"Correct CV price:   {correct['price']:.8f}")
    print(f"Incorrect CV price: {incorrect['price']:.8f}")
    print(f"Correct SE:         {correct['se']:.8f}")
    print(f"Incorrect SE:       {incorrect['se']:.8f}")
    print(f"Price shift:        {actual_shift:+.8f}")
    print(f"Predicted shift:    {predicted_shift:+.8f}")

    # Changing only mu must not change the SE
    assert np.isclose(
        correct["se"],
        incorrect["se"],
        atol=1e-12,
        rtol=1e-12
    )

    assert np.isclose(
        actual_shift,
        predicted_shift,
        atol=1e-12
    )

print("\nIncorrect-expectation shift tests passed.")
```

**저장된 실행 결과**

```text

=== BK sigma = 0.2 ===
Correct CV price:   1.04215044
Incorrect CV price: 1.03784782
Correct SE:         0.00000093
Incorrect SE:       0.00000093
Price shift:        -0.00430262
Predicted shift:    -0.00430262

=== BK sigma = 0.5 ===
Correct CV price:   1.03405042
Incorrect CV price: 1.02983861
Correct SE:         0.00000710
Incorrect SE:       0.00000710
Price shift:        -0.00421180
Predicted shift:    -0.00421180

Incorrect-expectation shift tests passed.
```

![잘못된 기대값이 CV 가격을 이동시키는 현상](/assets/images/bk_wrong_mu_shift.png)
*Figure. 잘못된 기대값이 CV 가격을 이동시키는 현상.*

**Python Implementation**

```python
# Task 13: Confidence interval coverage experiment
# Uses the Part D X1 control and independent simulations.

E_REPETITIONS = 80
E_N_COVERAGE = 2000
E_N_REFERENCE = 120000

e_coverage_results = {}

for sigma in E_SIGMAS:
    mu = d_control_expectations(sigma)["X1"]

    # Estimate nu using independent pilot paths
    Yp, Ip = d_simulate_bk(d_z_pilot, sigma)
    Xp = d_make_controls(Ip)["X1"]

    nu = np.cov(Xp, Yp, ddof=1)[0, 1] / np.var(
        Xp, ddof=1
    )

    # Independent, larger reference calculation.
    # Process in batches to avoid excessive memory use.
    reference_cv_sum = 0.0
    reference_cv_sumsq = 0.0
    processed = 0

    batch_size = 10000

    for batch in range(E_N_REFERENCE // batch_size):
        zb = d_generate_normals(
            batch_size,
            D_DT,
            700000 + 1000 * int(sigma * 10) + batch
        )

        Yb, Ib = d_simulate_bk(zb, sigma)
        Xb = d_make_controls(Ib)["X1"]

        cvb = Yb - nu * (Xb - mu)

        reference_cv_sum += np.sum(cvb)
        reference_cv_sumsq += np.sum(cvb**2)
        processed += len(cvb)

    reference_price = reference_cv_sum / processed

    reference_var = (
        reference_cv_sumsq
        - processed * reference_price**2
    ) / (processed - 1)

    reference_se = np.sqrt(
        max(reference_var, 0) / processed
    )

    covered_correct = 0
    covered_incorrect = 0

    widths_correct = []
    widths_incorrect = []

    for rep in range(E_REPETITIONS):
        z = d_generate_normals(
            E_N_COVERAGE,
            D_DT,
            900000 + 10000 * int(sigma * 10) + rep
        )

        Y, I = d_simulate_bk(z, sigma)
        X = d_make_controls(I)["X1"]

        cv_correct = Y - nu * (X - mu)
        cv_incorrect = Y - nu * (X - (mu + E_DELTA))

        price_correct = np.mean(cv_correct)
        price_incorrect = np.mean(cv_incorrect)

        se_correct = np.std(
            cv_correct, ddof=1
        ) / np.sqrt(E_N_COVERAGE)

        se_incorrect = np.std(
            cv_incorrect, ddof=1
        ) / np.sqrt(E_N_COVERAGE)

        # Approximate 95% normal confidence intervals
        lower_correct = price_correct - 1.96 * se_correct
        upper_correct = price_correct + 1.96 * se_correct

        lower_incorrect = price_incorrect - 1.96 * se_incorrect
        upper_incorrect = price_incorrect + 1.96 * se_incorrect

        covered_correct += (
            lower_correct <= reference_price <= upper_correct
        )

        covered_incorrect += (
            lower_incorrect <= reference_price <= upper_incorrect
        )

        widths_correct.append(upper_correct - lower_correct)
        widths_incorrect.append(upper_incorrect - lower_incorrect)

    e_coverage_results[sigma] = {
        "reference": reference_price,
        "reference_se": reference_se,
        "correct_coverage": covered_correct / E_REPETITIONS,
        "incorrect_coverage": covered_incorrect / E_REPETITIONS,
        "correct_width": np.mean(widths_correct),
        "incorrect_width": np.mean(widths_incorrect),
        "nu": nu
    }

    print(f"\n=== BK sigma = {sigma} ===")
    print(f"Reference price: {reference_price:.8f}")
    print(f"Reference SE:    {reference_se:.8f}")
    print(
        f"Correct mu coverage: "
        f"{covered_correct / E_REPETITIONS:.1%}"
    )
    print(
        f"Incorrect mu coverage: "
        f"{covered_incorrect / E_REPETITIONS:.1%}"
    )
    print(
        f"Correct mean CI width: "
        f"{np.mean(widths_correct):.8f}"
    )
    print(
        f"Incorrect mean CI width: "
        f"{np.mean(widths_incorrect):.8f}"
    )
```

**저장된 실행 결과**

```text

=== BK sigma = 0.2 ===
Reference price: 1.04215049
Reference SE:    0.00000037
Correct mu coverage: 92.5%
Incorrect mu coverage: 0.0%
Correct mean CI width: 0.00001127
Incorrect mean CI width: 0.00001127

=== BK sigma = 0.5 ===
Reference price: 1.03404736
Reference SE:    0.00000280
Correct mu coverage: 97.5%
Incorrect mu coverage: 0.0%
Correct mean CI width: 0.00008549
Incorrect mean CI width: 0.00008549
```

![올바른 기대값과 잘못된 기대값의 신뢰구간 포함률](/assets/images/bk_ci_coverage.png)
*Figure. 올바른 기대값과 잘못된 기대값의 신뢰구간 포함률.*

### Task 13. Results and Interpretation

σ=0.2와 σ=0.5 모두 잘못된 기대값 때문에 채권 가격이 크게 이동했지만 SE는 동일했다. 80회 반복한 실험에서 올바른 기대값을 썼을 때 포함률은 각각 **92.5%, 97.5%**, 잘못된 기대값에서는 둘 다 **0%**였다. 참고 가격도 분석해가 아닌 큰 독립 MC 기반 수치 추정값이다.

| BK σ | Correct CV price | Incorrect CV price | SE |
|---|---:|---:|---:|
| 0.2 | 1.04215044 | 1.03784782 | 0.00000093 |
| 0.5 | 1.03405042 | 1.02983861 | 0.00000710 |

## Task 14. Final Comparison of Estimators

마지막으로 Plain MC, Vasicek 통제변량의 기본·Matched·Grid 설정, BK 경로 기반 통제변량을 가격·SE·VRF·실행 시간으로 비교한다. 실행 시간은 하드웨어와 구현 순서에 의존하므로 절대적인 성능이라고 해석하지 않는다.

**Python Implementation**

```python
import time
import numpy as np

E_SUMMARY = []

# Previously generated independent samples
# d_z_pilot: coefficient estimation
# d_z_pricing: final pricing

def e_add_summary(sigma, method, result, elapsed):
    E_SUMMARY.append({
        "sigma": sigma,
        "method": method,
        "price": result["price"],
        "se": result["se"],
        "vrf": result["vrf"],
        "seconds": elapsed
    })


for sigma in E_SIGMAS:

    # Plain Monte Carlo
    start = time.perf_counter()

    Yp, Ip = d_simulate_bk(d_z_pilot, sigma)
    Y, I = d_simulate_bk(d_z_pricing, sigma)

    elapsed = time.perf_counter() - start

    plain_price = float(np.mean(Y))
    plain_se = float(
        np.std(Y, ddof=1) / np.sqrt(len(Y))
    )

    e_add_summary(
        sigma,
        "Plain MC",
        {
            "price": plain_price,
            "se": plain_se,
            "vrf": 1.0
        },
        elapsed
    )

    # Vasicek parameter configurations
    vasicek_choices = {
        "Baseline": (0.8, 0.03, 0.01),
        "Matched": (
            D_KAPPA,
            np.exp(D_THETA),
            sigma * np.exp(D_THETA)
        )
    }

    # Optionally include Part C grid-search choices
    if all(
        name in globals()
        for name in ["grid_results", "alpha_grid", "eta_grid"]
    ):
        for index, label in [
            (0, "Grid ZCB"),
            (1, "Grid Coupon")
        ]:
            vrf_grid = grid_results[sigma][index]

            iy, ix = np.unravel_index(
                np.argmax(vrf_grid),
                vrf_grid.shape
            )

            vasicek_choices[label] = (
                float(alpha_grid[ix]),
                0.03,
                float(eta_grid[iy])
            )

    # Part C: Vasicek controls
    for choice_name, params in vasicek_choices.items():
        alpha, beta, eta = params

        for control_name in ["ZCB", "Coupon"]:

            # A grid optimum is specific to its own control
            if choice_name == "Grid ZCB" and control_name != "ZCB":
                continue

            if choice_name == "Grid Coupon" and control_name != "Coupon":
                continue

            start = time.perf_counter()

            Xzp, Xcp = d_simulate_vasicek_controls(
                d_z_pilot, alpha, beta, eta
            )

            Xz, Xc = d_simulate_vasicek_controls(
                d_z_pricing, alpha, beta, eta
            )

            mu_z, mu_c = d_vasicek_expectations(
                alpha, beta, eta
            )

            if control_name == "ZCB":
                xp, xf, mu = Xzp, Xz, mu_z
            else:
                xp, xf, mu = Xcp, Xc, mu_c

            # Include target BK simulation cost
            Yp_local, _ = d_simulate_bk(
                d_z_pilot, sigma
            )
            Y_local, _ = d_simulate_bk(
                d_z_pricing, sigma
            )

            result = d_fit_cv(
                Yp_local, xp, Y_local, xf, [mu]
            )

            elapsed = time.perf_counter() - start

            e_add_summary(
                sigma,
                f"Vasicek {control_name} ({choice_name})",
                result,
                elapsed
            )

    # Part D: BK path-based controls
    configurations = {
        "I(T)": ["IT"],
        "I(T), I(T)^2": ["IT", "IT2"],
        "X1": ["X1"],
        "X1, X2": ["X1", "X2"]
    }

    for control_name, keys in configurations.items():
        start = time.perf_counter()

        Yp_local, Ip_local = d_simulate_bk(
            d_z_pilot, sigma
        )
        Y_local, I_local = d_simulate_bk(
            d_z_pricing, sigma
        )

        cp = d_make_controls(Ip_local)
        cf = d_make_controls(I_local)

        mu_all = d_control_expectations(sigma)

        xp = np.column_stack([cp[k] for k in keys])
        xf = np.column_stack([cf[k] for k in keys])
        mu = np.array([mu_all[k] for k in keys])

        result = d_fit_cv(
            Yp_local, xp, Y_local, xf, mu
        )

        elapsed = time.perf_counter() - start

        e_add_summary(
            sigma,
            f"BK {control_name}",
            result,
            elapsed
        )


# Print final comparison table
print("=== Task 14: Final Estimator Comparison ===")

header = (
    f"{'Sigma':>5} {'Method':36s} "
    f"{'Price':>12} {'SE':>12} "
    f"{'VRF':>9} {'Time(s)':>9}"
)
print(header)
print("-" * len(header))

for r in E_SUMMARY:
    print(
        f"{r['sigma']:5.1f} "
        f"{r['method']:36s} "
        f"{r['price']:12.8f} "
        f"{r['se']:12.8f} "
        f"{r['vrf']:9.2f} "
        f"{r['seconds']:9.4f}"
    )

# Identify highest-VRF method for each BK volatility
for sigma in E_SIGMAS:
    candidates = [
        r for r in E_SUMMARY
        if r["sigma"] == sigma
    ]

    best = max(candidates, key=lambda r: r["vrf"])

    print(f"\nBest VRF for BK sigma={sigma}:")
    print(f"Method: {best['method']}")
    print(f"VRF:    {best['vrf']:.2f}")
    print(f"SE:     {best['se']:.8f}")
```

**저장된 실행 결과**

```text
=== Task 14: Final Estimator Comparison ===
Sigma Method                                      Price           SE       VRF   Time(s)
----------------------------------------------------------------------------------------
  0.2 Plain MC                               1.04212403   0.00009197      1.00    0.1738
  0.2 Vasicek ZCB (Baseline)                 1.04216099   0.00000884    108.36    0.2027
  0.2 Vasicek Coupon (Baseline)              1.04215894   0.00000816    127.16    0.2297
  0.2 Vasicek ZCB (Matched)                  1.04216098   0.00000837    120.78    0.2395
  0.2 Vasicek Coupon (Matched)               1.04215895   0.00000766    144.01    0.2201
  0.2 Vasicek ZCB (Grid ZCB)                 1.04215872   0.00000730    158.51    0.1946
  0.2 Vasicek Coupon (Grid Coupon)           1.04215954   0.00000721    162.75    0.1975
  0.2 BK I(T)                                1.04215152   0.00000345    710.86    0.1435
  0.2 BK I(T), I(T)^2                        1.04215154   0.00000333    763.39    0.1511
  0.2 BK X1                                  1.04215044   0.00000093   9873.68    0.1466
  0.2 BK X1, X2                              1.04215038   0.00000001 57026630.30    0.1508
  0.5 Plain MC                               1.03400908   0.00024446      1.00    0.1385
  0.5 Vasicek ZCB (Baseline)                 1.03410593   0.00004926     24.63    0.1926
  0.5 Vasicek Coupon (Baseline)              1.03410052   0.00004851     25.39    0.1953
  0.5 Vasicek ZCB (Matched)                  1.03410588   0.00005083     23.13    0.2036
  0.5 Vasicek Coupon (Matched)               1.03410040   0.00005007     23.84    0.2662
  0.5 Vasicek ZCB (Grid ZCB)                 1.03410006   0.00004620     28.00    0.2173
  0.5 Vasicek Coupon (Grid Coupon)           1.03410222   0.00004607     28.15    0.1844
  0.5 BK I(T)                                1.03404968   0.00001145    456.13    0.1465
  0.5 BK I(T), I(T)^2                        1.03404802   0.00000892    751.80    0.1563
  0.5 BK X1                                  1.03405042   0.00000710   1186.30    0.1854
  0.5 BK X1, X2                              1.03404842   0.00000019 1722574.84    0.1450

Best VRF for BK sigma=0.2:
Method: BK X1, X2
VRF:    57026630.30
SE:     0.00000001

Best VRF for BK sigma=0.5:
Method: BK X1, X2
VRF:    1722574.84
SE:     0.00000019
```

### Task 14. Final Interpretation

실험 결과에서는 **BK X₁, X₂** 조합이 두 변동성 모두에서 가장 높은 VRF를 나타냈다. 다만 값이 극단적으로 큰 경우에도 기대값이 잘못되면 편향이 SE에 나타나지 않는다. 따라서 최종 방법을 선택할 때 VRF와 함께 기대값 검증, 이산화 오차, 연산 비용을 함께 고려한다.

# Final Validation

## 1. Executable Test Functions

검증 함수를 별도로 구성하고, 실행 결과를 **Test Register**로 정리한다. 각 테스트는 기대 결과, 관측 결과, 허용오차와 근거를 함께 남긴다.

**Python Implementation**

```python
# =========================================================
# FINAL VALIDATION - STEP 1: EXECUTABLE TEST FUNCTIONS
# =========================================================

final_test_results = []


def test_F01_bk_zero_volatility():
    """
    Validation ladder 2: Reduction to a known case.

    With sigma=0 and x0=theta, the BK rate is
    constant at 3%, so the bond price is known.
    """
    test_n = 1000

    price, se, simulated_payoffs = bk_mc(
        sigma=0.0,
        N=test_n,
        dt=1/52,
        seed=2026
    )

    expected = sum(
        (C + (principal if year == T else 0.0))
        * np.exp(-0.03 * year)
        for year in range(1, T + 1)
    )

    tolerance = 1e-10
    observed_error = abs(price - expected)

    assert observed_error <= tolerance
    assert np.allclose(
        simulated_payoffs, expected,
        rtol=0,
        atol=tolerance
    )
    assert se < tolerance

    return {
        "check": "F01",
        "property": "BK zero-volatility limit",
        "expected": expected,
        "observed": price,
        "tolerance": tolerance,
        "pass": True
    }


def test_F02_mc_standard_error():
    """
    Verify that the Monte Carlo standard error
    equals sample standard deviation / sqrt(N).
    """
    data = results_bk[0.2]

    simulated_payoffs = data["payoffs"]

    expected_se = (
        np.std(simulated_payoffs, ddof=1)
        / np.sqrt(len(simulated_payoffs))
    )

    observed_se = data["se"]

    tolerance = 1e-12

    assert abs(observed_se - expected_se) <= tolerance

    return {
        "check": "F02",
        "property": "MC standard error formula",
        "expected": expected_se,
        "observed": observed_se,
        "tolerance": tolerance,
        "pass": True
    }


def test_F03_control_expectations():
    """
    Validate all 18 independent control-expectation
    checks from Task 11 using the three-SE criterion.
    """

    assert len(e_task11_results) == 10
    assert len(e_task11_extra_results) == 8

    all_results = (
        e_task11_results + e_task11_extra_results
    )

    assert len(all_results) == 18

    for result in all_results:
        difference = abs(
            result["mc_mean"] - result["expected"]
        )

        tolerance = 3.0 * result["se"]

        assert np.isfinite(tolerance)
        assert tolerance > 0

        assert difference <= tolerance, (
            f"Control expectation failed: "
            f"{result['control']}, "
            f"error={difference:.8g}, "
            f"tolerance={tolerance:.8g}"
        )

        assert result["passed"]

    return {
        "check": "F03",
        "property": "Independent control expectations",
        "expected": "All 18 checks within 3 SE",
        "observed": "All 18 checks passed",
        "tolerance": "3 SE per control",
        "pass": True
    }

# Execute validation tests
for test_function in [
    test_F01_bk_zero_volatility,
    test_F02_mc_standard_error,
    test_F03_control_expectations
]:
    result = test_function()
    final_test_results.append(result)
    print(f"{result['check']}: PASS - {result['property']}")

print("\nAll Step 1 validation tests passed.")
```

**저장된 실행 결과**

```text
F01: PASS - BK zero-volatility limit
F02: PASS - MC standard error formula
F03: PASS - Independent control expectations

All Step 1 validation tests passed.
```

## 2. Control-Variate Formula Validation

수작업으로 계산할 수 있는 작은 데이터에 대해 CV 가격을 검증하고, 잘못된 기대값을 넣었을 때 가격 이동이 이론적 항등식과 맞는지 확인한다.

**Python Implementation**

```python
# =========================================================
# FINAL VALIDATION - STEP 2: CV FORMULA TESTS
# =========================================================

def test_F04_hand_computable_cv():
    """
    Validation ladder 1:
    Check the CV formula with a hand-computable dataset.
    """

    x_pilot = np.array([0.0, 1.0, 2.0])
    y_pilot = np.array([1.0, 3.0, 5.0])

    x_price = np.array([0.0, 1.0, 2.0])
    y_price = np.array([4.0, 8.0, 12.0])

    mu = 1.0

    result = control_variate(
        y_pilot, x_pilot,
        y_price, x_price, mu
    )

    # Hand-computed coefficient: nu = 2
    expected_nu = 2.0

    # Corrected payoffs: [6, 8, 10]
    expected_payoffs = np.array([6.0, 8.0, 10.0])

    expected_price = 8.0
    expected_se = (
        np.std(expected_payoffs, ddof=1)
        / np.sqrt(3)
    )

    tolerance = 1e-12

    assert abs(result["nu"] - expected_nu) <= tolerance
    assert abs(result["price"] - expected_price) <= tolerance
    assert abs(result["se"] - expected_se) <= tolerance

    return {
        "check": "F04",
        "property": "Hand-computed CV price",
        "expected": expected_price,
        "observed": result["price"],
        "tolerance": tolerance,
        "pass": True
    }


def test_F05_expectation_shift_identity():
    """
    Validation ladder 3:
    Changing mu by delta shifts the CV price by nu*delta,
    without changing SE or VRF.
    """

    x_pilot = np.array([0.0, 1.0, 2.0])
    y_pilot = np.array([1.0, 3.0, 5.0])

    x_price = np.array([0.0, 1.0, 2.0])
    y_price = np.array([4.0, 8.0, 12.0])

    mu = 1.0
    delta = 0.1

    correct = control_variate(
        y_pilot, x_pilot,
        y_price, x_price, mu
    )

    shifted = control_variate(
        y_pilot, x_pilot,
        y_price, x_price, mu + delta
    )

    expected_shift = correct["nu"] * delta
    observed_shift = shifted["price"] - correct["price"]

    tolerance = 1e-12

    assert abs(observed_shift - expected_shift) <= tolerance

    assert abs(
        shifted["se"] - correct["se"]
    ) <= tolerance

    assert abs(
        shifted["vrf"] - correct["vrf"]
    ) <= tolerance

    return {
        "check": "F05",
        "property": "CV expectation-shift identity",
        "expected": expected_shift,
        "observed": observed_shift,
        "tolerance": tolerance,
        "pass": True
    }


# Execute Step 2 tests
for test_function in [
    test_F04_hand_computable_cv,
    test_F05_expectation_shift_identity
]:
    result = test_function()
    final_test_results.append(result)

    print(
        f"{result['check']}: PASS - "
        f"{result['property']}"
    )

print("\nAll Step 2 validation tests passed.")
```

**저장된 실행 결과**

```text
F04: PASS - Hand-computed CV price
F05: PASS - CV expectation-shift identity

All Step 2 validation tests passed.
```

## 3. Deliberate-Bug Test

BK 통제변량 X₁의 기대값에 의도적으로 0.005를 더한 후 독립 기대값 검증을 실행한다. 올바른 기대값은 통과하고 잘못된 기대값은 **AssertionError로 탐지되어야 한다.** 예외는 테스트 안에서 기록하므로 Notebook 전체 재실행을 중단시키지 않는다.

**Python Implementation**

```python
# =========================================================
# FINAL VALIDATION - STEP 3: DELIBERATE-BUG TEST
# =========================================================

def test_F06_deliberate_wrong_expectation():
    """
    Verify that the independent expectation test
    detects an intentionally incorrect BK X1 mean.

    The expected AssertionError is caught so that
    Restart & Run All can finish successfully.
    """

    sigma = 0.2
    delta = 0.005

    # Independent simulation already used in Task 11
    _, I_test = d_simulate_bk(z_test, sigma, E_DT)
    controls = d_make_controls(I_test)

    x1_samples = np.asarray(controls["X1"])

    correct_mu = d_control_expectations(
        sigma, E_DT
    )["X1"]

    wrong_mu = correct_mu + delta

    # First verify that the CORRECT expectation passes.
    correct_result = e_expectation_check(
        "BK X1: correct expectation",
        x1_samples,
        correct_mu
    )

    assert correct_result["passed"]

    # Now deliberately introduce the incorrect mean.
    bug_detected = False
    failure_message = ""

    try:
        e_expectation_check(
            "BK X1: deliberately wrong expectation",
            x1_samples,
            wrong_mu
        )

    except AssertionError as error:
        bug_detected = True
        failure_message = str(error)

    # This final assertion verifies detection itself.
    assert bug_detected, (
        "The independent expectation test failed "
        "to detect the deliberately incorrect mean."
    )

    observed_mean = float(np.mean(x1_samples))
    observed_se = float(
        np.std(x1_samples, ddof=1)
        / np.sqrt(len(x1_samples))
    )

    tolerance = 3.0 * observed_se
    observed_error = abs(observed_mean - wrong_mu)

    print("=== DELIBERATE-BUG TEST EVIDENCE ===")
    print(f"BK sigma:         {sigma}")
    print(f"Correct mu:       {correct_mu:.10f}")
    print(f"Incorrect mu:     {wrong_mu:.10f}")
    print(f"Independent mean: {observed_mean:.10f}")
    print(f"Tolerance (3 SE): {tolerance:.10f}")
    print(f"Observed error:   {observed_error:.10f}")
    print()

    print("Expected failure: AssertionError")
    print("Actual result:    AssertionError caught")
    print(f"Failure message:  {failure_message}")
    print()
    print("Correct expectation: PASS")
    print("Incorrect expectation: DETECTED")
    print("Original pricing implementation: UNCHANGED")

    return {
        "check": "F06",
        "property": "Deliberately incorrect BK X1 expectation",
        "expected": "AssertionError on incorrect mu",
        "observed": "AssertionError caught",
        "tolerance": f"3 SE = {tolerance:.10g}",
        "pass": True
    }


f06_result = test_F06_deliberate_wrong_expectation()
final_test_results.append(f06_result)

print(
    f"\n{f06_result['check']}: PASS - "
    "Deliberate bug successfully detected"
)
```

**저장된 실행 결과**

```text
=== DELIBERATE-BUG TEST EVIDENCE ===
BK sigma:         0.2
Correct mu:       0.1698481644
Incorrect mu:     0.1748481644
Independent mean: 0.1697080592
Tolerance (3 SE): 0.0002593218
Observed error:   0.0051401052

Expected failure: AssertionError
Actual result:    AssertionError caught
Failure message:  BK X1: deliberately wrong expectation: expectation check failed. Difference=-0.00514011, Tolerance=0.00025932

Correct expectation: PASS
Incorrect expectation: DETECTED
Original pricing implementation: UNCHANGED

F06: PASS - Deliberate bug successfully detected
```

## 4. Test Register

F01부터 F06까지 여섯 가지 테스트의 결과를 실제 실행된 검증 함수에서 수집해 기록한다. 반복 실행 중 기록이 중복되지 않도록 정리하며 허용오차와 검출 가능한 오류 범위를 명시한다.

**Python Implementation**

```python
# =========================================================
# FINAL VALIDATION - STEP 4: TEST REGISTER
# =========================================================

import pandas as pd

test_justifications = {
    "F01": (
        "Deterministic BK limit: numerical tolerance "
        "1e-10 accounts for floating-point rounding."
    ),
    "F02": (
        "Standard-error identity: tolerance 1e-12 "
        "accounts for floating-point rounding."
    ),
    "F03": (
        "Independent MC expectation checks: 3 SE "
        "allows for sampling fluctuations. "
        "This is a statistical, not absolute, criterion."
    ),
    "F04": (
        "Hand-computable CV example: tolerance 1e-12 "
        "checks numerical agreement."
    ),
    "F05": (
        "CV expectation-shift identity: tolerance "
        "1e-12 checks the mathematical relationship."
    ),
    "F06": (
        "The deliberately wrong expectation must "
        "trigger an AssertionError in the independent "
        "3-SE expectation test."
    )
}

test_function_names = {
    "F01": "test_F01_bk_zero_volatility",
    "F02": "test_F02_mc_standard_error",
    "F03": "test_F03_control_expectations",
    "F04": "test_F04_hand_computable_cv",
    "F05": "test_F05_expectation_shift_identity",
    "F06": "test_F06_deliberate_wrong_expectation"
}

expected_checks = [f"F{i:02d}" for i in range(1, 7)]

# Keep the latest result for each test ID.
# This avoids duplicate entries when validation cells
# have been executed more than once.
unique_results = {}

for result in final_test_results:
    unique_results[result["check"]] = result

final_test_results = [
    unique_results[check_id]
    for check_id in expected_checks
    if check_id in unique_results
]

# Check for missing or unexpected tests.
recorded_checks = [
    result["check"] for result in final_test_results
]

assert len(recorded_checks) == len(set(recorded_checks)), (
    "Duplicate test IDs found in final_test_results."
)

assert set(recorded_checks) == set(expected_checks), (
    f"Missing or unexpected tests: {recorded_checks}"
)

# Verify that all registered test functions exist.
for check_id in expected_checks:
    function_name = test_function_names[check_id]

    assert callable(globals().get(function_name)), (
        f"Missing test function: {function_name}"
    )

# Generate the register from actual test results.
register_rows = []

for result in final_test_results:
    check_id = result["check"]

    register_rows.append({
        "Test ID": check_id,
        "Test function": test_function_names[check_id],
        "Property": result["property"],
        "Expected": result["expected"],
        "Observed": result["observed"],
        "Tolerance": result["tolerance"],
        "Tolerance justification": (
            test_justifications[check_id]
        ),
        "Pass": bool(result["pass"])
    })

test_register = pd.DataFrame(register_rows)

test_register = test_register.sort_values(
    "Test ID"
).reset_index(drop=True)

# All registered tests must have passed.
assert test_register["Pass"].all(), (
    "At least one validation test did not pass."
)

print("=== FINAL TEST REGISTER ===")

display(test_register)

print(
    f"\nTests passed: "
    f"{int(test_register['Pass'].sum())}"
    f"/{len(test_register)}"
)

print("Test Register generated successfully.")
```

**저장된 실행 결과**

```text
=== FINAL TEST REGISTER ===

  Test ID                          Test function  \
0     F01            test_F01_bk_zero_volatility   
1     F02             test_F02_mc_standard_error   
2     F03          test_F03_control_expectations   
3     F04            test_F04_hand_computable_cv   
4     F05    test_F05_expectation_shift_identity   
5     F06  test_F06_deliberate_wrong_expectation   

                                   Property                        Expected  \
0                  BK zero-volatility limit                        1.043659   
1                 MC standard error formula                        0.000092   
2          Independent control expectations       All 18 checks within 3 SE   
3                    Hand-computed CV price                             8.0   
4             CV expectation-shift identity                             0.2   
5  Deliberately incorrect BK X1 expectation  AssertionError on incorrect mu   

                Observed               Tolerance  \
0               1.043659                     0.0   
1               0.000092                     0.0   
2   All 18 checks passed        3 SE per control   
3                    8.0                     0.0   
4                    0.2                     0.0   
5  AssertionError caught  3 SE = 0.0002593218362   

                             Tolerance justification  Pass  
0  Deterministic BK limit: numerical tolerance 1e...  True  
1  Standard-error identity: tolerance 1e-12 accou...  True  
2  Independent MC expectation checks: 3 SE allows...  True  
3  Hand-computable CV example: tolerance 1e-12 ch...  True  
4  CV expectation-shift identity: tolerance 1e-12...  True  
5  The deliberately wrong expectation must trigge...  True


Tests passed: 6/6
Test Register generated successfully.
```

## 5. One Bug I Found

Test Register를 만드는 과정에서 검증 셀을 여러 번 실행하자 결과 목록에 같은 Test ID가 중복 저장되는 문제가 있었다. 원인은 Jupyter가 이전 실행의 리스트 상태를 유지한다는 점이었다. Test ID별 가장 최근 결과만 남기도록 수정한 뒤 **6/6 테스트 통과**를 확인했다. 실제 채권 가격이나 통제변량 계산식은 변경하지 않았다. 생성형 AI는 코드 작성·구조화·설명 정리의 보조 수단으로 사용했으며 최종 결과는 실행 결과를 기준으로 검토했다.

### AI Declaration
코드 구성과 설명을 정리하는 데 생성형 AI를 보조적으로 활용했다. 수식의 적용 방식, 구현 결과 및 검증 책임은 작성자에게 있으며, 가격 추정 결과는 실제 코드 실행과 테스트를 통해 확인했다.

---

## Final Summary

Part A–E에서 가격 추정과 통제변량의 분산 감소 효과를 비교했으며, **독립적인 기대값 검증 18개, BK 적분 규칙 비교, 의도적 오류 검출, 테스트 등록표 F01–F06**를 통해 구현 결과의 신뢰성을 함께 확인했다. 특히 작은 표준오차만으로는 편향을 발견할 수 없으므로, 통제변량의 기대값을 시뮬레이션 방식과 일치시키고 독립적으로 검증하는 절차가 중요하다.
