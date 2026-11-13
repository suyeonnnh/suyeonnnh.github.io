---
title: "Computational Finance : Milstein Scheme, Stochastic Volatility and Bond Pricing"
date: 2026-11-11 12:30:00 +0100
categories:
  - Computational Finance
tags:
  - Euler-Maruyama
  - Milstein Scheme
  - Stochastic Volatility
  - Correlated Brownian Motion
  - Zero-Coupon Bond
  - Coupon Bond
  - Vasicek
  - Black-Karasinski
  - Control Variates
math: true
---

이번 글은 Euler scheme과 Milstein scheme의 수치적 근사, correlated Brownian motion을 사용하는 stochastic volatility model, 채권의 할인 현금흐름 및 Monte Carlo pricing을 정리한다.

---

## 1. Euler Scheme for Stochastic Differential Equations

### 1.1 Time Discretization

1차원 stochastic differential equation(SDE)을 고려한다.

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BdX_t%3Da%28t%2CX_t%29%5C%2Cdt%2Bb%28t%2CX_t%29%5C%2CdW_t)

시간구간을 동일한 길이로 나눈다.

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3Bt_n%3Dn%5Cdelta%2C%5Cquad%20%5Cdelta%3D%5Cfrac%7BT%7D%7BN%7D%2C%5Cquad%20n%3D0%2C1%2C%5Cldots%2CN)

한 단계의 정확한 적분식은 다음과 같다.

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BX_%7Bt_%7Bn%2B1%7D%7D%3DX_%7Bt_n%7D%2B%5Cint_%7Bt_n%7D%5E%7Bt_%7Bn%2B1%7D%7Da%28s%2CX_s%29%5C%2Cds%2B%5Cint_%7Bt_n%7D%5E%7Bt_%7Bn%2B1%7D%7Db%28s%2CX_s%29%5C%2CdW_s)

각 적분항의 계수를 왼쪽 끝점에서 고정하면 Euler–Maruyama 갱신식을 얻는다.

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BX_%7Bt_%7Bn%2B1%7D%7D%5E%7B%5Cdelta%7D%3DX_%7Bt_n%7D%5E%7B%5Cdelta%7D%2Ba%28t_n%2CX_%7Bt_n%7D%5E%7B%5Cdelta%7D%29%5Cdelta%2Bb%28t_n%2CX_%7Bt_n%7D%5E%7B%5Cdelta%7D%29%28W_%7Bt_%7Bn%2B1%7D%7D-W_%7Bt_n%7D%29)

Brownian increment의 분포는 다음과 같다.

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BW_%7Bt_%7Bn%2B1%7D%7D-W_%7Bt_n%7D%5Csim%5Cmathcal%7BN%7D%280%2C%5Cdelta%29)

### 1.2 Strong and Weak Errors

Euler 근사의 오차는 strong error와 weak error로 구분한다.

**Strong error**는 같은 Brownian motion으로 정의된 정확한 해와 근사해 사이의 경로 오차를 평가한다.

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Csup_%7Bt%5Cin%5B0%2CT%5D%7D%5Cmathbb%7BE%7D%5B%7CX_t-X_t%5E%7B%5Cdelta%7D%7C%5D)

**Weak error**는 시험함수 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cvarphi)를 적용한 기댓값의 차이를 본다.

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%7C%5Cmathbb%7BE%7D%5B%5Cvarphi%28X_T%29%5D-%5Cmathbb%7BE%7D%5B%5Cvarphi%28X_T%5E%7B%5Cdelta%7D%29%5D%7C)

Euler scheme의 대표 수렴차수는 strong order ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B1%2F2), weak order ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B1)이다. 계수와 시험함수에 대한 정칙성 조건하에서 성립한다. Brownian increment의 표준편차는 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Csqrt%7B%5Cdelta%7D)이다.

## 2. Derivation of the Milstein Scheme

### 2.1 First-Order Expansion of the Diffusion Coefficient

시간에 명시적으로 의존하지 않는 1차원 SDE를 고려한다.

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BdX_t%3Da%28X_t%29%5C%2Cdt%2Bb%28X_t%29%5C%2CdW_t)

한 단계에서 drift 적분은 Euler 방식과 같이 근사한다.

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cint_%7Bt_n%7D%5E%7Bt_%7Bn%2B1%7D%7Da%28X_s%29%5C%2Cds%5Capprox%20a%28X_%7Bt_n%7D%29%28t_%7Bn%2B1%7D-t_n%29)

Diffusion coefficient는 Taylor expansion으로 근사한다. ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3Bb%28X_s%29)를 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BX_%7Bt_n%7D) 주변에서 일차 전개한다.

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3Bb%28X_s%29%5Capprox%20b%28X_%7Bt_n%7D%29%2Bb%5Cprime%28X_%7Bt_n%7D%29%28X_s-X_%7Bt_n%7D%29)

짧은 구간에서 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BX_s-X_%7Bt_n%7D)의 주된 확률항만 남기면

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BX_s-X_%7Bt_n%7D%5Capprox%20b%28X_%7Bt_n%7D%29%28W_s-W_%7Bt_n%7D%29)

이므로

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3Bb%28X_s%29%5Capprox%20b%28X_%7Bt_n%7D%29%2Bb%5Cprime%28X_%7Bt_n%7D%29b%28X_%7Bt_n%7D%29%28W_s-W_%7Bt_n%7D%29)

을 얻는다. 이 근사를 stochastic integral에 대입한다.

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cint_%7Bt_n%7D%5E%7Bt_%7Bn%2B1%7D%7Db%28X_s%29%5C%2CdW_s%5Capprox%20b%28X_%7Bt_n%7D%29%5CDelta%20W_n%2Bb%5Cprime%28X_%7Bt_n%7D%29b%28X_%7Bt_n%7D%29%5Cint_%7Bt_n%7D%5E%7Bt_%7Bn%2B1%7D%7D%28W_s-W_%7Bt_n%7D%29%5C%2CdW_s)

여기서 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5CDelta%20W_n%3DW_%7Bt_%7Bn%2B1%7D%7D-W_%7Bt_n%7D)이다.

### 2.2 Evaluation of the Stochastic Integral

다음 Itô 적분을 계산한다.

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cint_%7Bt_n%7D%5E%7Bt_%7Bn%2B1%7D%7D%28W_s-W_%7Bt_n%7D%29%5C%2CdW_s)

Itô formula에서 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3Bd%28W_t%5E2%29%3D2W_t%5C%2CdW_t%2Bdt)를 사용하면

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cint_0%5Et%20W_s%5C%2CdW_s%3D%5Cfrac12%28W_t%5E2-t%29)

가 된다. 이를 구간 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Bt_n%2Ct_%7Bn%2B1%7D%5D)에 적용하면

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cbegin%7Baligned%7D%5Cint_%7Bt_n%7D%5E%7Bt_%7Bn%2B1%7D%7D%28W_s-W_%7Bt_n%7D%29%5C%2CdW_s%26%3D%5Cfrac12%28W_%7Bt_%7Bn%2B1%7D%7D%5E2-W_%7Bt_n%7D%5E2-%5Cdelta%29-W_%7Bt_n%7D%28W_%7Bt_%7Bn%2B1%7D%7D-W_%7Bt_n%7D%29%5C%5C%26%3D%5Cfrac12%28%28W_%7Bt_%7Bn%2B1%7D%7D-W_%7Bt_n%7D%29%5E2-%5Cdelta%29.%5Cend%7Baligned%7D)

따라서

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cint_%7Bt_n%7D%5E%7Bt_%7Bn%2B1%7D%7D%28W_s-W_%7Bt_n%7D%29%5C%2CdW_s%3D%5Cfrac12%28%28%5CDelta%20W_n%29%5E2-%5Cdelta%29)

이다. **![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B-%5Cdelta) 항은 Itô calculus에서 나오는 보정항**이다.

### 2.3 Milstein Update

이 결과를 diffusion 근사식에 대입하면 최종 Milstein scheme은 다음과 같다.

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cbegin%7Baligned%7DX_%7Bt_%7Bn%2B1%7D%7D%5E%7B%5Cdelta%7D%3D%7B%7D%26X_%7Bt_n%7D%5E%7B%5Cdelta%7D%2Ba%28X_%7Bt_n%7D%5E%7B%5Cdelta%7D%29%5Cdelta%2Bb%28X_%7Bt_n%7D%5E%7B%5Cdelta%7D%29%5CDelta%20W_n%5C%5C%26%2B%5Cfrac12b%5Cprime%28X_%7Bt_n%7D%5E%7B%5Cdelta%7D%29b%28X_%7Bt_n%7D%5E%7B%5Cdelta%7D%29%28%28%5CDelta%20W_n%29%5E2-%5Cdelta%29.%5Cend%7Baligned%7D)

Euler scheme와 달리 diffusion coefficient의 도함수 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3Bb%27)와 제곱증분 항이 추가된다.

### 2.4 Constant Diffusion and Convergence Orders

Diffusion coefficient가 상수인 경우를 고려한다. ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3Bb%28x%29%3D%5Cbeta)로 놓으면 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3Bb%27%28x%29%3D0)이므로 Milstein의 추가 보정항이 사라지고 **Euler와 Milstein의 갱신식이 같아진다.**

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3Bb%5Cprime%28x%29%3D0%5Cquad%5CLongrightarrow%5Cquad%20%5Cmathrm%7BMilstein%7D%3D%5Cmathrm%7BEuler%7D)

Milstein scheme의 대표적인 수렴차수는 **strong order = 1**, **weak order = 1**이다. 충분한 정칙성과 표준 수렴 가정을 전제로 한다.

## 3. Multidimensional SDEs

### 3.1 Vector Form of the Euler Scheme

1차원 SDE를 벡터 형태로 일반화한다. 상태변수 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BX_t%5Cin%5Cmathbb%7BR%7D%5E%7Bd%7D)와 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3Br)차원 Brownian motion ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BW_t%5Cin%5Cmathbb%7BR%7D%5E%7Br%7D)에 대해

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BdX_t%3Da%28t%2CX_t%29%5C%2Cdt%2Bb%28t%2CX_t%29%5C%2CdW_t)

로 쓴다. 이때 계수들의 차원은 다음과 같다.

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3Ba%3A%5B0%2CT%5D%5Ctimes%5Cmathbb%7BR%7D%5E%7Bd%7D%5Cto%5Cmathbb%7BR%7D%5E%7Bd%7D%2C%5Cquad%20b%3A%5B0%2CT%5D%5Ctimes%5Cmathbb%7BR%7D%5E%7Bd%7D%5Cto%5Cmathbb%7BR%7D%5E%7Bd%5Ctimes%20r%7D)

따라서 Euler 갱신식도 벡터 형태를 그대로 유지한다.

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BX_%7Bt_%7Bn%2B1%7D%7D%5E%7B%5Cdelta%7D%3DX_%7Bt_n%7D%5E%7B%5Cdelta%7D%2Ba%28t_n%2CX_%7Bt_n%7D%5E%7B%5Cdelta%7D%29%5Cdelta%2Bb%28t_n%2CX_%7Bt_n%7D%5E%7B%5Cdelta%7D%29%5CDelta%20W_n)

여기서 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5CDelta%20W_n)은 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3Br)차원 벡터이고, 일반적인 독립 좌표 Brownian motion에서는 각 성분이 독립적인 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cmathcal%20N%280%2C%5Cdelta%29)를 따른다.

## 4. Stochastic Volatility with Correlated Brownian Motions

### 4.1 Stock and Variance Processes

주가 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BS_t)와 분산 상태변수 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BV_t)로 구성된 stochastic volatility model을 고려한다.

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BdS_t%3D%5Cmu%20S_t%5C%2Cdt%2BS_t%5Csqrt%7BV_t%7D%5C%2CdW_t%5E%7B%281%29%7D)

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BdV_t%3D%28%5Csigma_0%5E2-V_t%29%5C%2Cdt%2B%5Ckappa%5Csqrt%7BV_t%7D%5C%2CdW_t%5E%7B%282%29%7D)

두 Brownian motion은 상관계수 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Crho)를 가지며, quadratic covariation은 다음과 같다.

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Clangle%20W%5E%7B%281%29%7D%2CW%5E%7B%282%29%7D%5Crangle_t%3D%5Crho%20t)

### 4.2 Construction from Independent Brownian Motions

상관된 두 Brownian motion을 만들기 위해, 서로 독립인 좌표를 갖는 2차원 Brownian motion ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%28%5Cwidetilde%20W%5E%7B%281%29%7D%2C%5Cwidetilde%20W%5E%7B%282%29%7D%29)를 먼저 생성한다.

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BW_t%5E%7B%281%29%7D%3D%5Cwidetilde%7BW%7D_t%5E%7B%281%29%7D)

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BW_t%5E%7B%282%29%7D%3D%5Crho%5Cwidetilde%7BW%7D_t%5E%7B%281%29%7D%2B%5Csqrt%7B1-%5Crho%5E2%7D%5C%2C%5Cwidetilde%7BW%7D_t%5E%7B%282%29%7D)

독립성이 있으므로

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cbegin%7Baligned%7D%5Coperatorname%7BVar%7D%28W_t%5E%7B%282%29%7D%29%26%3D%5Crho%5E2t%2B%281-%5Crho%5E2%29t%3Dt%2C%5C%5C%5Coperatorname%7BCov%7D%28W_t%5E%7B%281%29%7D%2CW_t%5E%7B%282%29%7D%29%26%3D%5Crho%20t.%5Cend%7Baligned%7D)

이다. 따라서 두 번째 과정은 표준 Brownian motion이면서 첫 번째 과정과 원하는 상관관계를 가진다. 여기서 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B-1%5Cle%5Crho%5Cle1)이다.

### 4.3 Diffusion Matrix

독립인 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cwidetilde%20W)를 기준으로 위 2차원 SDE를 다시 쓰면 diffusion matrix는 다음과 같다.

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3Bb%28S_t%2CV_t%29%3D%5Cbegin%7Bpmatrix%7DS_t%5Csqrt%7BV_t%7D%260%5C%5C%5Ckappa%5Csqrt%7BV_t%7D%5Crho%26%5Ckappa%5Csqrt%7BV_t%7D%5Csqrt%7B1-%5Crho%5E2%7D%5Cend%7Bpmatrix%7D)

첫 번째 열은 공통 Brownian shock가 주가와 분산에 함께 영향을 주는 부분이고, 두 번째 열은 분산 과정에만 들어가는 독립 shock의 부분이다.

## 5. Zero-Coupon and Coupon Bond Pricing

### 5.1 Zero-Coupon Bond

Zero-coupon bond는 만기 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BT)에 **1을 지급**한다. 시간은 연 단위로 계산한다. 확률적 단기이자율 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3Br_t)를 사용할 때, 위험중립측도에서 시점 0의 가격을 다음처럼 나타낸다.

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BP%280%2CT%29%3D%5Cmathbb%7BE%7D%5E%7B%5Cmathbb%7BQ%7D%7D%5Cleft%5B%5Cexp%5Cleft%28-%5Cint_0%5ET%20r_s%5C%2Cds%5Cright%29%5Cright%5D)

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cmathbb%20E)

### 5.2 Coupon Bond

Coupon bond는 지급 시점 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3Bt_1%2C%5Cldots%2Ct_n)에 금액 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BC_1%2C%5Cldots%2CC_n)을 지급한다. 일반적인 연 1회 지급 구조는 다음과 같다.

- ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BC_1%3D%5Ccdots%3DC_%7Bn-1%7D%3DC): 중간 지급액은 동일한 coupon
- ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BC_n%3D1%2BC): 마지막에는 원금 1과 coupon을 함께 지급
- ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3Bt_i%3Di): 연 1회 지급하는 예시

시점 0의 가격은 모든 현금흐름의 할인 기댓값의 합이다.

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cbegin%7Baligned%7DV_0%26%3D%5Cmathbb%7BE%7D%5E%7B%5Cmathbb%7BQ%7D%7D%5Cleft%5B%5Csum_%7Bi%3D1%7D%5En%20C_i%5Cexp%5Cleft%28-%5Cint_0%5E%7Bt_i%7Dr_s%5C%2Cds%5Cright%29%5Cright%5D%5C%5C%26%3D%5Csum_%7Bi%3D1%7D%5En%20C_i%5C%2C%5Cmathbb%7BE%7D%5E%7B%5Cmathbb%7BQ%7D%7D%5Cleft%5B%5Cexp%5Cleft%28-%5Cint_0%5E%7Bt_i%7Dr_s%5C%2Cds%5Cright%29%5Cright%5D%5C%5C%26%3D%5Csum_%7Bi%3D1%7D%5En%20C_iP%280%2Ct_i%29.%5Cend%7Baligned%7D)

이 등식은 기댓값의 선형성에 의해 성립한다.

## 6. Monte Carlo Pricing and Control Variates

### 6.1 Discount-Integral Notation

할인인자의 Taylor expansion을 위해 적분값을 정의한다.

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BI_t%3D%5Cint_0%5Et%20r_u%5C%2Cdu)

그러면

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3Be%5E%7B-I_t%7D%5Capprox1-I_t%2B%5Cfrac12I_t%5E2)

이다. 다음 변수들을 control variate 후보로 고려한다.

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BI_t%2C%5Cqquad%20I_t%5E2)

Coupon bond에서는 여러 지급 시점의 할인적분을 선형결합한다.

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Csum_i%20c_iI_%7Bt_i%7D%2C%5Cqquad%5Csum_i%20c_iI_%7Bt_i%7D%5E2)

이 항들의 기댓값을 알고 있거나 별도로 계산할 수 있을 때 통제변량으로 활용할 수 있다. 

### 6.2 Why the Control Variate Is Useful

할인인자 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BY%3De%5E%7B-I_t%7D)를 Monte Carlo로 계산할 때, 기댓값이 알려진 변수 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BC)를 사용하면

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BY_%7B%5Cmathrm%7BCV%7D%7D%3DY-%5Ctheta%28C-%5Cmathbb%7BE%7D%5BC%5D%29)

처럼 평균을 유지하면서 분산을 줄이도록 설계할 수 있다. ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Ctheta)의 선택에 따라 분산이 달라진다.

## 7. Computational Finance Project Stages

Vasicek와 Black–Karasinski 모형을 이용한 채권 가격 계산의 실험 단계는 다음과 같다.

### 7.1 Stage 1 — Analytical Vasicek Bond Price

먼저 Vasicek 모형에서 zero-coupon bond의 가격 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BP_%7B%5Cmathrm%7BVas%7D%7D%280%2CT%29)에 대한 해석적 공식을 다룬다. 가격 검증 조건은 다음과 같다.

- 만기가 바로 현재이면 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BP_%7B%5Cmathrm%7BVas%7D%7D%280%2C0%29%3D1)
- 이자율이 상수 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3Br_0)가 되는 특수한 경우 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BP%280%2CT%29%3De%5E%7B-r_0T%7D)
- 교재의 공식과 비교

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BP_%7B%5Cmathrm%7BVas%7D%7D%280%2C0%29%3D1%2C%5Cqquad%20r_t%5Cequiv%20r_0%5C%20%5CLongrightarrow%5C%20P%280%2CT%29%3De%5E%7B-r_0T%7D)

### 7.2 Stage 2 — Euler Scheme and Discounted Payoffs

Vasicek 및 Black–Karasinski(BK) 모형의 단기이자율 경로를 Euler scheme으로 생성하고, 채권 payoff와 그 moment를 확인한다. 한 단계의 계산을 직접 수행하고 Vasicek/BK의 결과를 비교한다.

확률적인 이자율 경로가 주어졌을 때 할인적분을 시간격자에서 근사하는 일반적인 형태는 다음과 같다.

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BI_T%5Capprox%5Csum_%7Bn%3D0%7D%5E%7BN-1%7Dr_%7Bt_n%7D%5Cdelta%2C%5Cquad%20e%5E%7B-I_T%7D%5Capprox%5Cexp%5Cleft%28-%5Csum_%7Bn%3D0%7D%5E%7BN-1%7Dr_%7Bt_n%7D%5Cdelta%5Cright%29)

### 7.3 Stage 3 — Plain Monte Carlo for Vasicek

Vasicek의 zero-coupon bond 가격을 일반 Monte Carlo로 계산한다. Stage 1의 해석가격이 구현 점검 기준이 된다.

### 7.4 Stage 4 — Plain Monte Carlo for Black–Karasinski

BK 모형에서도 같은 만기·지급 구조에 대해 plain Monte Carlo를 실행한다.

### 7.5 Stage 5 — Vasicek as a Control Variate

Vasicek에서 계산한 값을 BK 가격 추정의 control variate로 활용한다. 실제 사용 시에는 두 추정량을 같은 시뮬레이션에서 어떻게 결합하고, 어떤 기대값을 기준값으로 넣을지 정해야 한다.

---
