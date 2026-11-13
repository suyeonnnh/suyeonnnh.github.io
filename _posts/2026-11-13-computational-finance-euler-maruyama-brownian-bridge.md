---
title: "Computational Finance : Euler–Maruyama and Brownian Bridge"
date: 2026-11-13 12:00:00 +0100
categories:
  - Computational Finance
tags:
  - Monte Carlo
  - Brownian Motion
  - Euler-Maruyama
  - Stochastic Differential Equations
  - Strong Convergence
  - Weak Convergence
  - Brownian Bridge
math: true
---

 Monte Carlo에서 확률변수를 생성하는 문제부터 Brownian motion의 이산화, strong/weak convergence, Euler–Maruyama scheme, SDE의 해의 조건, Brownian bridge의 조건부분포 유도까지 필기 순서와 수식을 중심으로 정리했다.

---

## 1. Monte Carlo and the Simulation Problem

Monte Carlo에서는 확률변수 X의 기댓값을 구하려고 한다.

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cmathbb%7BE%7D%5BX)]

앞선 수업에서는 같은 기댓값을 유지하면서 분산을 줄이는 변량감소 방법을 다뤘다.

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cmathbb%7BE%7D%5BY%5D%3D%5Cmathbb%7BE%7D%5BX%5D%2C%5Cqquad%20%5Coperatorname%7BVar%7D%28Y%29%3C%5Coperatorname%7BVar%7D%28X%29)

이번 필기의 문제는 **X의 분포로부터 표본을 어떻게 생성하는가**이다. 특히 X가 확률미분방정식(SDE)의 해에 의해 정의되는 경우를 생각한다. SDE의 해를 정확히 시뮬레이션하기 어려우면 이산 시간격자 위에서 경로를 근사한다.

## 2. Simulation of Brownian Motion

### 2.1 Time Grid

시간구간 [0,T]와 분할 수 N을 정하고 시간간격을 정의한다.

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5CDelta%3D%5Cfrac%7BT%7D%7BN%7D%2C%5Cqquad%20t_i%3Di%5CDelta%2C%5Cquad%20i%3D0%2C%5Cldots%2CN)

필기의 도식처럼 0에서 T까지 N개의 구간을 만든다.

### 2.2 Independent Gaussian Increments

표준 Brownian motion은 0에서 시작한다.

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BW_0%3D0)

임의의 t_i와 t_{i+1}에서 증분은 다음 분포를 따른다.

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BW_%7Bt_%7Bi%2B1%7D%7D-W_%7Bt_i%7D%5Csim%5Cmathcal%7BN%7D%280%2Ct_%7Bi%2B1%7D-t_i%29%3D%5Cmathcal%7BN%7D%280%2C%5CDelta%29)

겹치지 않는 구간의 증분은 서로 독립이다. 따라서 매 단계 표준정규 확률변수를 생성한다.

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BZ_i%5Coverset%7B%5Cmathrm%7Bi.i.d.%7D%7D%7B%5Csim%7D%5Cmathcal%7BN%7D%280%2C1%29)

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cwidehat%20W_%7Bt_0%7D%3D0%2C%5Cqquad%20%5Cwidehat%20W_%7Bt_%7Bi%2B1%7D%7D%3D%5Cwidehat%20W_%7Bt_i%7D%2B%5Csqrt%7B%5CDelta%7D%5C%2CZ_i)

이렇게 구한 값은 격자점에서 Brownian motion의 분포를 정확히 재현한다. **중간 시점의 값**은 뒤의 Brownian bridge에서 따로 생성한다.

## 3. Strong and Weak Convergence

정확한 과정 X와 근사 과정 X^Δ가 있을 때, 시간간격 Δ를 0으로 보내면서 오차를 비교한다.

### 3.1 Strong Convergence

Strong convergence는 같은 확률적 입력으로 정의한 두 경로가 얼마나 가까운지를 살핀다.

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Csup_%7Bt%5Cin%5B0%2CT%5D%7D%5Cmathbb%7BE%7D%5Cbig%5B%7CX_t-X_t%5E%7B%5CDelta%7D%7C%5Cbig%5D%5Cxrightarrow%5B%5CDelta%5Cto0%5D%7B%7D0)

**Strong order γ**는 아래 오차 상계로 표현한다.

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Csup_%7Bt%5Cin%5B0%2CT%5D%7D%5Cmathbb%7BE%7D%5Cbig%5B%7CX_t-X_t%5E%7B%5CDelta%7D%7C%5Cbig%5D%5Cle%20K%5CDelta%5E%7B%5Cgamma%7D)

상수 K는 Δ에 의존하지 않는다. 필기에서 strong convergence는 **경로 사이 오차**를 기준으로 정의했다.

### 3.2 Weak Convergence

Weak convergence는 각 경로의 오차 대신 시험함수 φ를 적용한 기댓값을 비교한다.

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cleft%7C%5Cmathbb%7BE%7D%5B%5Cvarphi%28X_T%29%5D-%5Cmathbb%7BE%7D%5B%5Cvarphi%28X_T%5E%7B%5CDelta%7D%29%5D%5Cright%7C%5Cxrightarrow%5B%5CDelta%5Cto0%5D%7B%7D0)

필기에는 최종 시점 T에서의 기댓값이 적혀 있고, 아래처럼 **weak order γ**를 정의한다.

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cleft%7C%5Cmathbb%7BE%7D%5B%5Cvarphi%28X_T%29%5D-%5Cmathbb%7BE%7D%5B%5Cvarphi%28X_T%5E%7B%5CDelta%7D%29%5D%5Cright%7C%5Cle%20K_%7B%5Cvarphi%7D%5CDelta%5E%7B%5Cgamma%7D)

함수 φ는 임의의 함수 전체가 아니라, 정리에서 지정한 함수 종류에 속해야 한다. 필기에는 continuous, bounded 함수가 예시로 적혀 있다. 수렴차수 정리에서는 추가적인 매끄러움 조건이 필요할 수 있다.

### 3.3 Euler–Maruyama Convergence Orders

칠판에 쓰인 대표적인 결과는 다음과 같다.

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Ctext%7BStrong%20convergence%20order%7D%3D%5Cfrac12)

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Ctext%7BWeak%20convergence%20order%7D%3D1)

이는 필요한 계수·초기값·시험함수 조건이 성립하는 표준적인 경우의 결과다. Strong order 1/2는 경로 오차가 통상 Δ의 제곱근 차수로 감소한다는 뜻이고, weak order 1은 기댓값 오차가 통상 Δ 차수로 감소한다는 뜻이다.

### 3.4 Problems with Euler–Maruyama

칠판에는 다음 문제가 적혀 있다.

- **Exploding coefficients:** 계수가 급격하게 커지는 경우에는 단순한 Euler–Maruyama scheme이 불안정하거나 발산할 수 있다.
- **State-space constraints:** 해가 (0,∞) 또는 [0,1] 같은 구간에 있어야 하더라도 이산 근사값이 그 구간 밖으로 나갈 수 있다. 예를 들어 양수여야 할 근사값이 음수가 될 수 있다.

## 4. Derivation of the Euler–Maruyama Scheme

### 4.1 Stochastic Differential Equation

필기의 SDE는 다음과 같다.

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BdX_t%3Da%28t%2CX_t%29%5C%2Cdt%2Bb%28t%2CX_t%29%5C%2CdW_t%2C%5Cqquad%20X_0%3Dx_0)

여기서 a는 drift, b는 diffusion coefficient다. 한 시간구간 [t_i,t_{i+1}]에서 적분형으로 쓰면

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BX_%7Bt_%7Bi%2B1%7D%7D%3DX_%7Bt_i%7D%2B%5Cint_%7Bt_i%7D%5E%7Bt_%7Bi%2B1%7D%7Da%28s%2CX_s%29%5C%2Cds%2B%5Cint_%7Bt_i%7D%5E%7Bt_%7Bi%2B1%7D%7Db%28s%2CX_s%29%5C%2CdW_s)

이다.

### 4.2 Freezing the Coefficients

구간 길이가 작을 때, 두 적분의 계수를 구간 왼쪽 끝점에서 평가하여 근사한다.

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cint_%7Bt_i%7D%5E%7Bt_%7Bi%2B1%7D%7Da%28s%2CX_s%29%5C%2Cds%5Capprox%20a%28t_i%2CX_%7Bt_i%7D%29%28t_%7Bi%2B1%7D-t_i%29)

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cint_%7Bt_i%7D%5E%7Bt_%7Bi%2B1%7D%7Db%28s%2CX_s%29%5C%2CdW_s%5Capprox%20b%28t_i%2CX_%7Bt_i%7D%29%28W_%7Bt_%7Bi%2B1%7D%7D-W_%7Bt_i%7D%29)

### 4.3 Update Formula

따라서 격자점에서의 Euler–Maruyama 근사식은

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BX_%7Bt_%7Bi%2B1%7D%7D%5E%7B%5CDelta%7D%3DX_%7Bt_i%7D%5E%7B%5CDelta%7D%2Ba%28t_i%2CX_%7Bt_i%7D%5E%7B%5CDelta%7D%29%28t_%7Bi%2B1%7D-t_i%29%2Bb%28t_i%2CX_%7Bt_i%7D%5E%7B%5CDelta%7D%29%28W_%7Bt_%7Bi%2B1%7D%7D-W_%7Bt_i%7D%29)

이다. 등간격 시간격자에서

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5CDelta%20W_i%3DW_%7Bt_%7Bi%2B1%7D%7D-W_%7Bt_i%7D%3D%5Csqrt%7B%5CDelta%7D%5C%2CZ_i)

를 사용하면 실제 갱신식은 다음과 같다.

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cboxed%7BX_%7Bt_%7Bi%2B1%7D%7D%5E%7B%5CDelta%7D%3DX_%7Bt_i%7D%5E%7B%5CDelta%7D%2Ba%28t_i%2CX_%7Bt_i%7D%5E%7B%5CDelta%7D%29%5CDelta%2Bb%28t_i%2CX_%7Bt_i%7D%5E%7B%5CDelta%7D%29%5Csqrt%7B%5CDelta%7D%5C%2CZ_i%7D)

초기 근사값은 정확한 초기값과 동일하게 둔다.

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BX_%7Bt_0%7D%5E%7B%5CDelta%7D%3Dx_0)

## 5. Vectorized Simulation and Memory Usage

### 5.1 Multiple Paths at Once

필기에는 **특정 시점 T에서의 기댓값**을 구하는 계산 예시가 있고, 그 뒤 여러 경로를 동시에 생성하는 방법이 적혀 있다.

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cmathbb%7BE%7D%5Bf%28X_T%29)]

Monte Carlo 표본 경로를 M개 만들 때, 각 경로의 값을 하나의 M차원 벡터에 저장한다.

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cmathbf%20X_i%3D%5Cbegin%7Bpmatrix%7DX_%7Bt_i%7D%5E%7B%281%29%7D%5C%5C%5Cvdots%5C%5CX_%7Bt_i%7D%5E%7B%28M%29%7D%5Cend%7Bpmatrix%7D%5Cin%5Cmathbb%20R%5EM)

각 시간단계에서 M개의 독립 표준정규변수를 생성한다.

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cboldsymbol%5Cxi_i%5Csim%5Cmathcal%20N%28%5Cmathbf%200%2CI_M%29%2C%5Cqquad%20%5Csqrt%7B%5CDelta%7D%5Cboldsymbol%5Cxi_i%5Csim%5Cmathcal%20N%28%5Cmathbf%200%2C%5CDelta%20I_M%29)

그러면 Euler–Maruyama update를 성분별로 한꺼번에 수행할 수 있다.

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cmathbf%20X_%7Bi%2B1%7D%3D%5Cmathbf%20X_i%2B%5Cmathbf%20a%28t_i%2C%5Cmathbf%20X_i%29%5CDelta%2B%5Cmathbf%20b%28t_i%2C%5Cmathbf%20X_i%29%5Csqrt%7B%5CDelta%7D%5Cboldsymbol%5Cxi_i)

여기서 굵은 a와 b는 각 성분에 같은 계수함수를 적용한다는 의미로 썼다.

### 5.2 Simulation Procedure

1. M개 경로를 초기값으로 설정한다.
2. i=0부터 N−1까지 반복한다.
3. 각 단계에서 M차원 Gaussian vector를 하나 생성한다.
4. Euler–Maruyama 갱신식을 전체 벡터에 적용한다.
5. 마지막 시점의 벡터 성분에 f를 적용해 표본평균을 계산한다.

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cwidehat%5Cmu_M%3D%5Cfrac1M%5Csum_%7Bj%3D1%7D%5E%7BM%7Df%5Cbigl%28X_T%5E%7B%5CDelta%2C%28j%29%7D%5Cbigr%29)

### 5.3 Computational Advice

필기 마지막에는 **가능한 한 적은 메모리를 사용하는 것**이 좋다고 강조했다. 최종 시점 X_T만 필요한 계산이라면 각 경로의 전체 시간 이력을 모두 저장할 필요가 없다. 현재 시점의 벡터를 다음 시점의 벡터로 덮어쓰면 된다.

필기에는 CPU의 L1/L2/L3 cache와 RAM 구조도 그려져 있다. 요점은 전체 경로 행렬을 불필요하게 유지하지 않고, 메모리 사용량과 데이터 이동을 줄이는 것이다. 반대로 경로 전체에 의존하는 payoff를 계산한다면 필요한 경로 정보는 보존해야 한다.

## 6. Strong Solutions of SDEs

### 6.1 Integral Form and Adaptedness

필기에서는 위 SDE를 적분형으로 다시 적고, 해 X가 주어진 filtration에 적합(adapted)해야 한다고 표시했다.

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BX_t%3DX_0%2B%5Cint_0%5Et%20a%28s%2CX_s%29%5C%2Cds%2B%5Cint_0%5Et%20b%28s%2CX_s%29%5C%2CdW_s%2C%5Cqquad%20t%5Cin%5B0%2CT)]

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BX_t%5Ctext%7B%20is%20adapted%20to%20%7D%28%5Cmathcal%20F_t%29_%7Bt%5Cin%5B0%2CT%5D%7D)

또한 확률적분이 정의되도록 다음 제곱적분가능성 조건이 적혀 있다.

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cmathbb%20E%5Cleft%5B%5Cint_0%5ET%7Cb%28s%2CX_s%29%7C%5E2%5C%2Cds%5Cright%5D%3C%5Cinfty)

필기 옆에는 적응성 조건과 함께 **almost surely**라는 메모가 있다. 위 적분 등식은 모든 t에 대해 거의 확실하게 성립하는 해의 표현이다.

### 6.2 Lipschitz Condition

Strong solution의 존재성과 유일성을 보장하는 표준적인 충분조건 중 하나는 다음 전역 Lipschitz 조건이다.

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%7Ca%28t%2Cx%29-a%28t%2Cy%29%7C%2B%7Cb%28t%2Cx%29-b%28t%2Cy%29%7C%5Cle%20L%7Cx-y%7C)

어떤 L>0에 대해 모든 x,y 및 t∈[0,T]에서 성립하도록 가정한다.

### 6.3 Linear Growth Condition

다른 가정은 다음 linear growth condition이다.

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%7Ca%28t%2Cx%29%7C%2B%7Cb%28t%2Cx%29%7C%5Cle%20K%281%2B%7Cx%7C%29)

여기서 K>0이며, 필기에는 a,b의 measurability도 함께 적혀 있다. 적절한 초기값과 통상적인 확률공간 조건 아래 이러한 조건은 강한 해의 존재성과 경로별 유일성에 사용된다.

## 7. Brownian Bridge: Intermediate Time Points

### 7.1 The Conditional Sampling Problem

Brownian motion을 격자점에서 생성한 뒤, 이미 알고 있는 두 값 사이에서 시점 s의 값을 생성하려고 한다.

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3Bt_i%3Cs%3Ct_%7Bi%2B1%7D)

그때 알고 있는 값은 W_{t_i}, W_{t_{i+1}}이다. Brownian bridge의 목표는 이 두 값이 주어졌을 때 **W_s의 조건부 정규분포**를 구하는 것이다.

### 7.2 Independent Increments

다음처럼 서로 인접한 두 증분을 정의한다.

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BX%3DW_s-W_%7Bt_i%7D%2C%5Cqquad%20Y%3DW_%7Bt_%7Bi%2B1%7D%7D-W_s)

독립증분 성질에 의해 X와 Y는 독립이며 함께 정규분포를 따른다.

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BX%5Csim%5Cmathcal%20N%280%2Cs-t_i%29%2C%5Cqquad%20Y%5Csim%5Cmathcal%20N%280%2Ct_%7Bi%2B1%7D-s%29)

또한

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BX%2BY%3DW_%7Bt_%7Bi%2B1%7D%7D-W_%7Bt_i%7D)

이다. 따라서 **X+Y가 주어졌을 때 X의 조건부분포**를 구하면 된다.

### 7.3 Finding the Projection Coefficient

필기에서는 다음과 같은 선형결합 R을 만든다.

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BR%3DX-%5Calpha%28X%2BY%29%3D%281-%5Calpha%29X-%5Calpha%20Y)

R이 X+Y와 독립이 되도록 α를 고른다. (X,Y)는 jointly Gaussian이므로 선형결합 (R,X+Y)도 jointly Gaussian이다. 따라서 두 변수의 공분산이 0이면 독립이다.

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Coperatorname%7BCov%7D%28R%2CX%2BY%29%3D0)

독립성으로 Cov(X,Y)=0이므로

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cbegin%7Baligned%7D0%26%3D%5Coperatorname%7BCov%7D%28X-%5Calpha%28X%2BY%29%2CX%2BY%29%5C%5C%26%3D%5Coperatorname%7BVar%7D%28X%29-%5Calpha%5Coperatorname%7BVar%7D%28X%2BY%29%5C%5C%26%3D%281-%5Calpha%29%28s-t_i%29-%5Calpha%28t_%7Bi%2B1%7D-s%29%5Cend%7Baligned%7D)

을 얻는다. 이를 정리하면

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%28s-t_i%29%3D%5Calpha%28t_%7Bi%2B1%7D-t_i%29)

따라서

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cboxed%7B%5Calpha%3D%5Cfrac%7Bs-t_i%7D%7Bt_%7Bi%2B1%7D-t_i%7D%7D)

이다.

### 7.4 Independence and the Conditional Mean

위에서 구한 α를 사용하면 R과 X+Y는 독립이다. 이 부분이 필기에서 계산한 **독립성 증명**의 목적이다.

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BX%3D%5Calpha%28X%2BY%29%2BR)

R의 평균은 0이고 X+Y와 독립이므로

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cmathbb%20E%5BX%5Cmid%20X%2BY%5D%3D%5Calpha%28X%2BY%29)

가 된다. 따라서

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cmathbb%20E%5BW_s%5Cmid%20W_%7Bt_i%7D%2CW_%7Bt_%7Bi%2B1%7D%7D%5D%3DW_%7Bt_i%7D%2B%5Cfrac%7Bs-t_i%7D%7Bt_%7Bi%2B1%7D-t_i%7D%28W_%7Bt_%7Bi%2B1%7D%7D-W_%7Bt_i%7D%29)

이다. 조건부 평균은 두 끝점을 연결하는 **선형보간**이다.

### 7.5 Deriving the Conditional Variance

R은 X+Y와 독립이므로 조건부 분산은 Var(R)이다.

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cbegin%7Baligned%7D%5Coperatorname%7BVar%7D%28R%29%26%3D%5Coperatorname%7BVar%7D%28%281-%5Calpha%29X-%5Calpha%20Y%29%5C%5C%26%3D%281-%5Calpha%29%5E2%28s-t_i%29%2B%5Calpha%5E2%28t_%7Bi%2B1%7D-s%29%5Cend%7Baligned%7D)

이때

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B1-%5Calpha%3D%5Cfrac%7Bt_%7Bi%2B1%7D-s%7D%7Bt_%7Bi%2B1%7D-t_i%7D)

이므로

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cbegin%7Baligned%7D%5Coperatorname%7BVar%7D%28R%29%26%3D%5Cfrac%7B%28t_%7Bi%2B1%7D-s%29%5E2%28s-t_i%29%2B%28s-t_i%29%5E2%28t_%7Bi%2B1%7D-s%29%7D%7B%28t_%7Bi%2B1%7D-t_i%29%5E2%7D%5C%5C%26%3D%5Cfrac%7B%28s-t_i%29%28t_%7Bi%2B1%7D-s%29%7D%7Bt_%7Bi%2B1%7D-t_i%7D%5Cend%7Baligned%7D)

을 얻는다.

### 7.6 Brownian Bridge Distribution

따라서 두 끝점이 주어졌을 때 중간값의 조건부분포는 다음과 같다.

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cboxed%7BW_s%5Cmid%28W_%7Bt_i%7D%2CW_%7Bt_%7Bi%2B1%7D%7D%29%5Csim%5Cmathcal%20N%5Cleft%28W_%7Bt_i%7D%2B%5Cfrac%7Bs-t_i%7D%7Bt_%7Bi%2B1%7D-t_i%7D%28W_%7Bt_%7Bi%2B1%7D%7D-W_%7Bt_i%7D%29%2C%5C%20%5Cfrac%7B%28s-t_i%29%28t_%7Bi%2B1%7D-s%29%7D%7Bt_%7Bi%2B1%7D-t_i%7D%5Cright%29%7D)

### 7.7 Brownian Bridge Sampling Formula

위 조건부분포에서 Z∼N(0,1)을 독립적으로 생성하면

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cboxed%7BW_s%3DW_%7Bt_i%7D%2B%5Cfrac%7Bs-t_i%7D%7Bt_%7Bi%2B1%7D-t_i%7D%28W_%7Bt_%7Bi%2B1%7D%7D-W_%7Bt_i%7D%29%2B%5Csqrt%7B%5Cfrac%7B%28s-t_i%29%28t_%7Bi%2B1%7D-s%29%7D%7Bt_%7Bi%2B1%7D-t_i%7D%7D%5C%2CZ%7D)

를 사용해 중간값을 생성할 수 있다. 선형보간만 사용하면 조건부 분산이 0이 되어 올바른 Brownian motion의 중간값 분포를 재현하지 못한다.

---
