---
layout: post
title: "Computational Finance : Monte Carlo Greeks and Binomial Trees"
date: 2026-11-13
categories: [Computational Finance]
tags: [Monte-Carlo, Greeks, Malliavin, Binomial-Tree, Black-Scholes, CRR]
math: false
---

# Monte Carlo Greeks and Binomial Trees

옵션 가격의 기초자산 민감도(Delta)를 Monte Carlo로 계산하는 방법을 정리하고, 이어서 Black–Scholes 모델을 이항트리(binomial tree)로 근사하여 옵션 가격을 구하는 과정을 살펴본다. 핵심은 **기대값을 미분하는 문제**, **Monte Carlo 추정량의 분산**, **이항트리의 모멘트 일치**, **만기 격자와 행사가격의 위치**다.

## 1. Black–Scholes Model and Option Pricing

위험중립측도 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cmathbb%20Q) 아래에서 기초자산 가격은


![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BdS_t%3DrS_t%5C%2Cdt%2B%5Csigma%20S_t%5C%2CdW_t%2C%0A%5Cqquad%20B_t%3De%5E%7Brt%7D)


을 만족한다. 따라서 만기 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BT)의 가격은


![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BS_T%3DS%5Cexp%5Cleft%28%5Cleft%28r-%5Cfrac%7B%5Csigma%5E2%7D%7B2%7D%5Cright%29T%2B%5Csigma%20W_T%5Cright%29)


이다. 여기서 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BS%3DS_0)이고 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BW_T%5Csim%20N%280%2CT%29)다. 만기 지급액이 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3Bh%28S_T%29)인 유럽형 경로독립 옵션의 현재 가격은


![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cboxed%7Bv%28S%29%3D%5Cmathbb%20E_%7B%5Cmathbb%20Q%7D%5Cleft%5Be%5E%7B-rT%7Dh%28S_T%29%5Cright%5D%0A%3De%5E%7B-rT%7D%5Cmathbb%20E_%7B%5Cmathbb%20Q%7D%5Cleft%5Bh%5Cleft%28S%20e%5E%7B%28r-%5Csigma%5E2%2F2%29T%2B%5Csigma%20W_T%7D%5Cright%29%5Cright%5D.%7D)


Monte Carlo에서는 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BW_T%3D%5Csqrt%20T%20Z), ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BZ%5Csim%20N%280%2C1%29)을 생성하고 지급액을 평균한다. 그러나 **Delta**를 계산하려면 가격 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3Bv%28S%29)를 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BS)에 대해 미분해야 한다.


![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5CDelta%28S%29%3D%5Cfrac%7B%5Cpartial%20v%28S%29%7D%7B%5Cpartial%20S%7D.)


## 2. Likelihood-Ratio Derivative: Full Derivation

지급액 함수 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3Bh)가 매끄럽지 않아도 적용할 수 있는 기대값 미분 공식을 유도한다. 변수


![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3Ba%3D%5Cleft%28r-%5Cfrac%7B%5Csigma%5E2%7D%7B2%7D%5Cright%29T%2C%0A%5Cqquad%20U%3D%5Clog%20S_T%3D%5Clog%20S%2Ba%2B%5Csigma%20W_T)


를 두면 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BU%5Csim%20N%28%5Clog%20S%2Ba%2C%5Csigma%5E2T%29)이므로


![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3Bv%28S%29%3De%5E%7B-rT%7D%5Cint_%7B-%5Cinfty%7D%5E%7B%5Cinfty%7Dh%28e%5Eu%29%0A%5Cfrac%7B1%7D%7B%5Csqrt%7B2%5Cpi%7D%5Csigma%5Csqrt%20T%7D%0A%5Cexp%5Cleft%5B-%5Cfrac%7B%28u-%5Clog%20S-a%29%5E2%7D%7B2%5Csigma%5E2T%7D%5Cright%5Ddu.)


![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3Bh%28e%5Eu%29) 자체를 미분하는 대신 **정규분포의 밀도에서 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BS)에 의존하는 부분**을 미분한다. ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3Bm%3D%5Clog%20S%2Ba)로 놓으면 지수 부분은


![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3Bg%28S%29%3D-%5Cfrac%7B%28u-m%29%5E2%7D%7B2%5Csigma%5E2T%7D.)


연쇄법칙을 사용하면


![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cfrac%7B%5Cpartial%20g%7D%7B%5Cpartial%20S%7D%0A%3D-%5Cfrac%7B2%28u-m%29%7D%7B2%5Csigma%5E2T%7D%5Cleft%28-%5Cfrac1S%5Cright%29%0A%3D%5Cfrac%7Bu-%5Clog%20S-a%7D%7BS%5Csigma%5E2T%7D.)


따라서 적절한 적분 가능성 조건 아래


![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cbegin%7Baligned%7D%0Av%27%28S%29%0A%26%3De%5E%7B-rT%7D%5Cint_%7B-%5Cinfty%7D%5E%7B%5Cinfty%7Dh%28e%5Eu%29%0A%5Cfrac%7Be%5E%7B-%28u-%5Clog%20S-a%29%5E2%2F%282%5Csigma%5E2T%29%7D%7D%7B%5Csqrt%7B2%5Cpi%7D%5Csigma%5Csqrt%20T%7D%0A%5Cfrac%7Bu-%5Clog%20S-a%7D%7BS%5Csigma%5E2T%7D%5C%2Cdu.%0A%5Cend%7Baligned%7D)


이제 표준화 변수


![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cboxed%7Bz%3D%5Cfrac%7Bu-%5Clog%20S-a%7D%7B%5Csigma%5Csqrt%20T%7D%7D%2C%0A%5Cqquad%20du%3D%5Csigma%5Csqrt%20T%5C%2Cdz)


를 도입한다. 이때 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3Bu%3D%5Clog%20S%2Ba%2B%5Csigma%5Csqrt%20T%20z)이므로


![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cbegin%7Baligned%7D%0Av%27%28S%29%0A%26%3D%5Cfrac%7Be%5E%7B-rT%7D%7D%7BS%5Csigma%5Csqrt%20T%7D%0A%5Cint_%7B-%5Cinfty%7D%5E%7B%5Cinfty%7D%0A%20h%5Cleft%28e%5E%7B%5Clog%20S%2Ba%2B%5Csigma%5Csqrt%20Tz%7D%5Cright%29%0A%20z%5Cfrac%7Be%5E%7B-z%5E2%2F2%7D%7D%7B%5Csqrt%7B2%5Cpi%7D%7D%5C%2Cdz%5C%5C%0A%26%3D%5Cfrac%7Be%5E%7B-rT%7D%7D%7BS%5Csigma%5Csqrt%20T%7D%0A%5Cmathbb%20E%5Cleft%5Bh%5Cleft%28S%20e%5E%7Ba%2B%5Csigma%5Csqrt%20T%20Z%7D%5Cright%29Z%5Cright%5D.%0A%5Cend%7Baligned%7D)


![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BW_T%3D%5Csqrt%20T%20Z)이므로 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BZ%3DW_T%2F%5Csqrt%20T)이고, 최종적으로


![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cboxed%7B%5Cdisplaystyle%0Av%27%28S%29%3D%5Cfrac%7Be%5E%7B-rT%7D%7D%7BS%5Csigma%20T%7D%0A%5Cmathbb%20E_%7B%5Cmathbb%20Q%7D%5Cleft%5B%0A%20h%5Cleft%28S%20e%5E%7B%28r-%5Csigma%5E2%2F2%29T%2B%5Csigma%20W_T%7D%5Cright%29W_T%0A%5Cright%5D.%7D)


이 식은 **likelihood-ratio / score-function 방식**의 Delta 공식이다. 지급액 함수의 도함수 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3Bh%27)가 등장하지 않으므로 디지털 옵션처럼 불연속인 지급액에도 활용할 수 있다. 단, 실제 Monte Carlo 분산이 작다는 보장은 없다.

## 3. Derivatives via Monte Carlo

### 3.1 Plain finite differences

중앙차분으로 Delta를 근사하면


![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5CDelta%28S%29%5Capprox%5Cfrac%7Bv%28S%2B%5Cdelta%29-v%28S-%5Cdelta%29%7D%7B2%5Cdelta%7D.)


각 가격을 별개의 Monte Carlo 표본으로 계산하면 두 추정오차가 차분에서 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B1%2F%282%5Cdelta%29)만큼 증폭된다. 특히 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cdelta)가 작을수록 문제가 심해질 수 있다.

### 3.2 Variance reduction with common random numbers

같은 정규난수 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BZ_i)를 두 가격에 공통으로 사용하면


![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cwidehat%5CDelta_%7B%5Cmathrm%7BCRN%7D%7D%0A%3D%5Cfrac%7Be%5E%7B-rT%7D%7D%7BN%7D%5Csum_%7Bi%3D1%7D%5EN%0A%5Cfrac%7Bh%28%28S%2B%5Cdelta%29Y_i%29-h%28%28S-%5Cdelta%29Y_i%29%7D%7B2%5Cdelta%7D%2C%0A%5Cqquad%20Y_i%3De%5E%7B%28r-%5Csigma%5E2%2F2%29T%2B%5Csigma%5Csqrt%20T%20Z_i%7D.)


독립 추정량 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BX_%2B), ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BX_-)의 차이와 달리 공통난수를 사용하면


![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Coperatorname%7BVar%7D%28X_%2B-X_-%29%0A%3D%5Coperatorname%7BVar%7D%28X_%2B%29%2B%5Coperatorname%7BVar%7D%28X_-%29%0A-2%5Coperatorname%7BCov%7D%28X_%2B%2CX_-%29)


에서 **양의 공분산**을 활용할 수 있다. 지급액의 형태에 따라 실제 효과는 달라진다.

### 3.3 Malliavin estimator

칠판의 **Malliavin estimator (ML)**는 앞에서 유도한 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BW_T) 가중치 공식을 Monte Carlo로 계산하는 방식이다.


![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cboxed%7B%5Cdisplaystyle%0A%5Cwidehat%5CDelta_%7B%5Cmathrm%7BML%7D%7D%0A%3D%5Cfrac%7Be%5E%7B-rT%7D%7D%7BS%5Csigma%20T%7D%5Cfrac1N%0A%5Csum_%7Bi%3D1%7D%5EN%20h%28SY_i%29W_T%5E%7B%28i%29%7D.%7D)


비교 대상으로 사용한 지급액은 다음 두 가지다.

- **Binary call:** ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3Bh%28S_T%29%3D%5Cmathbf1_%7B%5C%7BS_T%5Cge%20K%5C%7D%7D). 행사가격 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BK)를 넘으면 1, 아니면 0을 지급한다.
- **Vanilla call:** ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3Bh%28S_T%29%3D%28S_T-K%29%5E%2B%3D%5Cmax%28S_T-K%2C0%29).

### 3.4 Numerical results

칠판의 실험은 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BN%3D1%7B%2C%7D000%7B%2C%7D000)을 사용한다. VR은 variance reduction을 적용한 차분 추정량, ML은 Malliavin 추정량을 나타낸다. 괄호 안의 두 값은 판서 그대로 기록한 **두 수치**이며, 각 값의 정확한 통계적 의미(반복 결과인지, 추정치와 다른 지표인지)는 별도로 명시되어 있지 않다.

| Payoff | VR | ML | Variance comparison |
|---|---|---|---|
| Binary call | ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%280.1147%2C%5C%200.1423%29) | ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%280.1233%2C%5C%200.1301%29) | ML에서 분산 약 **1,400배 감소**로 표기 |
| Vanilla call | ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%280.4550%2C%5C%200.4631%29) | ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%280.4454%2C%5C%200.4612%29) | ML에서 분산 약 **4배 증가**로 표기 |

결과의 중요한 특징은 **ML이 모든 지급액에서 더 좋은 방법은 아니라는 점**이다. 디지털 옵션처럼 지급액이 불연속인 경우에는 미분이 없는 공식이 특히 유용할 수 있지만, 바닐라 콜에서는 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BW_T)를 곱하는 가중치가 오히려 분산을 키울 수 있다. 표에 적힌 배수는 판서의 실험 결과이지 보편적인 이론적 배수가 아니다.

## 4. Portfolio Immunisation

여러 위험요인 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3Bx_1%2C%5Cldots%2Cx_m)에 따라 포트폴리오 가치가 변한다고 하자.


![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BV%3DV%28x_1%2C%5Cldots%2Cx_m%29.)


위험요인이 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5CDelta%20x)만큼 변할 때 다변수 Taylor 전개는


![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cbegin%7Baligned%7D%0AV%28x%2B%5CDelta%20x%29-V%28x%29%0A%26%5Capprox%20%5Csum_%7Bi%3D1%7D%5E%7Bm%7D%5Cfrac%7B%5Cpartial%20V%7D%7B%5Cpartial%20x_i%7D%5CDelta%20x_i%5C%5C%0A%26%5Cquad%2B%5Cfrac12%5Csum_%7Bi%3D1%7D%5E%7Bm%7D%5Csum_%7Bj%3D1%7D%5E%7Bm%7D%0A%5Cfrac%7B%5Cpartial%5E2V%7D%7B%5Cpartial%20x_i%5Cpartial%20x_j%7D%5CDelta%20x_i%5CDelta%20x_j.%0A%5Cend%7Baligned%7D)


**Portfolio immunisation**은 포트폴리오가 작은 위험요인 변화에 덜 민감하도록 편미분값을 0 또는 작은 값으로 만드는 접근이다. 1차 편미분을 0으로 맞추면 작은 변화에 대한 1차 영향을 제거하고, 필요에 따라 2차 민감도까지 관리한다. Delta, Gamma 등 Greeks가 이런 위험 관리에 연결된다.

## 5. Binomial Model under the Risk-Neutral Measure

Black–Scholes 모형에서는


![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BdS_t%3DrS_t%5C%2Cdt%2B%5Csigma%20S_t%5C%2CdW_t%2C%0A%5Cqquad%20S_t%3DS_0e%5E%7B%28r-%5Csigma%5E2%2F2%29t%2B%5Csigma%20W_t%7D%2C%0A%5Cqquad%20B_t%3De%5E%7Brt%7D.)


만기 지급액이 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3Bh%28S_T%29)인 유럽형 옵션의 가격은


![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cboxed%7BV_0%3D%5Cmathbb%20E_%7B%5Cmathbb%20Q%7D%5Be%5E%7B-rT%7Dh%28S_T%29%5D%7D.)


이 기대값을 빠르게 근사하기 위해 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5B0%2CT%5D)를 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BN)개의 기간으로 나누고 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5CDelta%20t%3DT%2FN)으로 설정한다. 매 기간 주가가 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3Bu)배 상승하거나 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3Bd)배 하락하는 **이항트리**를 구성한다.

![recombining binomial tree]({{ "/assets/img/cf-post4-recombining-binomial-tree.svg" | relative_url }})

위 트리는 네 기간의 재결합 이항트리다. 예를 들어 두 번 움직인 뒤에는 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BS_0u%5E2), ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BS_0ud), ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BS_0d%5E2)가 가능하다. ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3Bud%3Ddu)이므로 상승 후 하락과 하락 후 상승은 같은 노드로 합쳐진다. 일반적으로 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3Bi)기간 후 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3Bj)번 하락한 노드는


![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cboxed%7BS_%7Bi%2Cj%7D%3DS_0u%5E%7Bi-j%7Dd%5Ej%2C%5Cqquad%20j%3D0%2C1%2C%5Cldots%2Ci.%7D)


단순히 모든 경로를 나열하면 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B2%5EN)개의 경로가 생기지만, 재결합 트리의 서로 다른 노드 수는


![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B1%2B2%2B%5Ccdots%2B%28N%2B1%29%3D%5Cboxed%7B%5Cfrac%7B%28N%2B1%29%28N%2B2%29%7D%7B2%7D%7D)


개다. 칠판에서 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BN%3D100) 또는 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BN%3D1000)처럼 기간 수가 커질 때 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B2%5EN) 방식이 비현실적인 이유를 보여준 부분이다. 재결합 트리를 이용하면 저장과 계산을 크게 줄일 수 있다.

## 6. Constructing a Binomial Tree by Moment Matching

한 기간에서 이항트리 수익률은


![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cfrac%7B%5Cwidehat%20S_%7Bt%2B%5CDelta%20t%7D%7D%7B%5Cwidehat%20S_t%7D%0A%3D%5Cbegin%7Bcases%7Du%2C%26%5Ctext%7Bwith%20probability%20%7Dp%2C%5C%5Cd%2C%26%5Ctext%7Bwith%20probability%20%7D1-p.%5Cend%7Bcases%7D)


위험중립 Black–Scholes 모형의 정확한 한 기간 수익률은


![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cfrac%7BS_%7Bt%2B%5CDelta%20t%7D%7D%7BS_t%7D%0A%3D%5Cexp%5Cleft%5B%5Cleft%28r-%5Cfrac%7B%5Csigma%5E2%7D%7B2%7D%5Cright%29%5CDelta%20t%0A%2B%5Csigma%28W_%7Bt%2B%5CDelta%20t%7D-W_t%29%5Cright%5D.)


**Moment matching**은 이 두 수익률의 첫 번째 모멘트와 두 번째 모멘트를 같게 만드는 방법이다.

### 6.1 First moment

이항트리의 기대 수익률은


![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cmathbb%20E%5Cleft%5B%5Cfrac%7B%5Cwidehat%20S_%7Bt%2B%5CDelta%20t%7D%7D%7B%5Cwidehat%20S_t%7D%5Cright%5D%0A%3Dpu%2B%281-p%29d.)


Black–Scholes에서는 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BW_%7Bt%2B%5CDelta%20t%7D-W_t%5Csim%20N%280%2C%5CDelta%20t%29)이고 정규분포의 지수적률 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cmathbb%20E%5Be%5E%7BaZ%7D%5D%3De%5E%7Ba%5E2%2F2%7D)를 사용하면


![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cbegin%7Baligned%7D%0A%5Cmathbb%20E%5Cleft%5B%5Cfrac%7BS_%7Bt%2B%5CDelta%20t%7D%7D%7BS_t%7D%5Cright%5D%0A%26%3De%5E%7B%28r-%5Csigma%5E2%2F2%29%5CDelta%20t%7D%0A%5Cmathbb%20E%5Be%5E%7B%5Csigma%28W_%7Bt%2B%5CDelta%20t%7D-W_t%29%7D%5D%5C%5C%0A%26%3De%5E%7B%28r-%5Csigma%5E2%2F2%29%5CDelta%20t%7De%5E%7B%5Csigma%5E2%5CDelta%20t%2F2%7D%0A%3De%5E%7Br%5CDelta%20t%7D.%0A%5Cend%7Baligned%7D)


따라서 첫 번째 조건은


![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cboxed%7Bpu%2B%281-p%29d%3De%5E%7Br%5CDelta%20t%7D.%7D%5Ctag%7B1%7D)


### 6.2 Second moment

이항트리의 두 번째 모멘트는


![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cmathbb%20E%5Cleft%5B%5Cleft%28%5Cfrac%7B%5Cwidehat%20S_%7Bt%2B%5CDelta%20t%7D%7D%7B%5Cwidehat%20S_t%7D%5Cright%29%5E2%5Cright%5D%0A%3Dpu%5E2%2B%281-p%29d%5E2.)


Black–Scholes 수익률을 제곱하면


![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cleft%28%5Cfrac%7BS_%7Bt%2B%5CDelta%20t%7D%7D%7BS_t%7D%5Cright%29%5E2%0A%3De%5E%7B2%28r-%5Csigma%5E2%2F2%29%5CDelta%20t%2B2%5Csigma%28W_%7Bt%2B%5CDelta%20t%7D-W_t%29%7D.)


따라서


![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cbegin%7Baligned%7D%0A%5Cmathbb%20E%5Cleft%5B%5Cleft%28%5Cfrac%7BS_%7Bt%2B%5CDelta%20t%7D%7D%7BS_t%7D%5Cright%29%5E2%5Cright%5D%0A%26%3De%5E%7B%282r-%5Csigma%5E2%29%5CDelta%20t%7D%0A%5Cmathbb%20E%5Be%5E%7B2%5Csigma%28W_%7Bt%2B%5CDelta%20t%7D-W_t%29%7D%5D%5C%5C%0A%26%3De%5E%7B%282r-%5Csigma%5E2%29%5CDelta%20t%7De%5E%7B2%5Csigma%5E2%5CDelta%20t%7D%5C%5C%0A%26%3De%5E%7B%282r%2B%5Csigma%5E2%29%5CDelta%20t%7D.%0A%5Cend%7Baligned%7D)


두 번째 조건은


![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cboxed%7Bpu%5E2%2B%281-p%29d%5E2%3De%5E%7B%282r%2B%5Csigma%5E2%29%5CDelta%20t%7D.%7D%5Ctag%7B2%7D)


### 6.3 Why the solution is not unique

식 (1), (2)는 **두 개의 방정식**이지만 미지수는 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3Bp%2Cu%2Cd)의 **세 개**다. 따라서 추가 조건을 지정하지 않으면 해가 유일하게 정해지지 않는다.

두 식에 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BS_t)가 나타나지 않는다는 점도 중요하다. 선택한 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3Bp%2Cu%2Cd)를 현재 주가에 따라 다시 바꾸지 않고 모든 노드에 동일하게 적용할 수 있어 재결합 트리를 만들기 쉽다.

첫 번째 조건을 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3Bp)에 대해 풀면


![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3Bpu%2B%281-p%29d%3De%5E%7Br%5CDelta%20t%7D%0A%5Cquad%5CLongrightarrow%5Cquad%0Ap%28u-d%29%3De%5E%7Br%5CDelta%20t%7D-d)


이므로


![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cboxed%7Bp%3D%5Cfrac%7Be%5E%7Br%5CDelta%20t%7D-d%7D%7Bu-d%7D.%7D)


이 확률이 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B0%5Cle%20p%5Cle1)이 되려면 보통 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3Bd%5Cle%20e%5E%7Br%5CDelta%20t%7D%5Cle%20u)가 필요하다. 엄격한 무차익 조건을 적용하는 비퇴화 트리에서는 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3Bd%3Ce%5E%7Br%5CDelta%20t%7D%3Cu)를 사용한다.

## 7. Jarrow–Rudd and Cox–Ross–Rubinstein Trees

### 7.1 Jarrow–Rudd

**Jarrow–Rudd**에서는


![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cboxed%7Bp%3D%5Cfrac12%7D)


로 놓는다. 두 모멘트를 정확히 일치시키려면 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3Bm%3De%5E%7Br%5CDelta%20t%7D), ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3Bq%3De%5E%7B%282r%2B%5Csigma%5E2%29%5CDelta%20t%7D)에 대해


![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cfrac%7Bu%2Bd%7D%7B2%7D%3Dm%2C%0A%5Cqquad%20%5Cfrac%7Bu%5E2%2Bd%5E2%7D%7B2%7D%3Dq.)


![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3Bu%3Dm%2Ba), ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3Bd%3Dm-a)로 놓으면


![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cfrac%7B%28m%2Ba%29%5E2%2B%28m-a%29%5E2%7D%7B2%7D%3Dm%5E2%2Ba%5E2%3Dq%2C)


따라서


![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cboxed%7Bu%3Dm%2B%5Csqrt%7Bq-m%5E2%7D%2C%5Cqquad%20d%3Dm-%5Csqrt%7Bq-m%5E2%7D.%7D)


또한 널리 쓰는 대칭 로그수익률 형태는


![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3Bu%3De%5E%7B%28r-%5Csigma%5E2%2F2%29%5CDelta%20t%2B%5Csigma%5Csqrt%7B%5CDelta%20t%7D%7D%2C%0A%5Cqquad%0Ad%3De%5E%7B%28r-%5Csigma%5E2%2F2%29%5CDelta%20t-%5Csigma%5Csqrt%7B%5CDelta%20t%7D%7D%2C%0A%5Cqquad%20p%3D%5Ctfrac12%2C)


이며 이는 작은 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5CDelta%20t)에서 모멘트를 근사한다. 두 매개화는 구별해야 한다.

### 7.2 Cox–Ross–Rubinstein (CRR)

**CRR**에서는


![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cboxed%7Bd%3D%5Cfrac1u%7D)


라는 조건을 둔다. 표준적으로 사용하는 근사는


![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cboxed%7Bu%3De%5E%7B%5Csigma%5Csqrt%7B%5CDelta%20t%7D%7D%2C%5Cqquad%20d%3De%5E%7B-%5Csigma%5Csqrt%7B%5CDelta%20t%7D%7D%2C%5Cqquad%0Ap%3D%5Cfrac%7Be%5E%7Br%5CDelta%20t%7D-d%7D%7Bu-d%7D.%7D)


이 선택은 첫 번째 위험중립 모멘트를 정확히 만족시키며, 두 번째 모멘트는 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5CDelta%20t%5Cto0)에서 점근적으로 일치한다. 즉 **표준 CRR의 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3Bu%2Cd)가 위 두 모멘트를 유한한 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5CDelta%20t)에서 모두 정확히 맞추는 것은 아니다.**

## 8. Backward Induction for European Calls

만기에서 각 노드의 옵션 가치는 지급액 자체다.


![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BV_%7BN%2Cj%7D%3Dh%28S_%7BN%2Cj%7D%29%2C%0A%5Cqquad%20j%3D0%2C%5Cldots%2CN.)


유럽형 콜이라면


![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cboxed%7BV_%7BN%2Cj%7D%3D%28S_%7BN%2Cj%7D-K%29%5E%2B.%7D)


한 기간 이전의 가치는 위험중립 기대값을 한 기간 할인하여 계산한다.


![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cboxed%7B%5Cdisplaystyle%0AV_%7Bi%2Cj%7D%3De%5E%7B-r%5CDelta%20t%7D%5Cleft%5BpV_%7Bi%2B1%2Cj%7D%2B%281-p%29V_%7Bi%2B1%2Cj%2B1%7D%5Cright%5D%2C%0A%5Cqquad%20i%3DN-1%2C%5Cldots%2C0.%7D)


여기서 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3Bj)는 하락 횟수이므로 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BV_%7Bi%2B1%2Cj%7D)가 상승 자식 노드이고 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BV_%7Bi%2B1%2Cj%2B1%7D)이 하락 자식 노드다. 마지막 시점부터 시작해 한 단계씩 거슬러 올라가면 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BV_%7B0%2C0%7D)을 얻는다. 모든 경로를 개별 계산하거나 거대한 행렬을 만들 필요가 없다.

## 9. Oscillations in Binomial Call Prices

기간 수 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BN)을 늘리면 일반적으로 이항트리 가격이 Black–Scholes 가격으로 가까워지지만, 수렴이 매끄럽게 한쪽 방향으로만 진행되지는 않을 수 있다. 칠판의 그래프는 **옵션 가격(세로축)**과 **기간 수(가로축)**의 관계에서 가격이 위아래로 진동하면서 안정되는 현상을 보여준다.

![binomial price oscillations]({{ "/assets/img/cf-post4-binomial-price-oscillations.svg" | relative_url }})

*그래프는 칠판에 그려진 진동·수렴 형태를 설명하기 위한 개념도다. 실제 실험 데이터나 정확한 가격 수치를 재현한 그래프는 아니다.*

유럽형 콜의 만기 지급액은


![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3Bh%28S_T%29%3D%28S_T-K%29%5E%2B)


이고, ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BS_T%3DK)에서 기울기가 바뀐다. ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BN)이 달라지면 만기 트리의 노드들이 행사가격 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BK)를 기준으로 배치되는 위치도 달라진다. 그에 따라 근사오차가 번갈아 변하며 진동이 나타날 수 있다.

## 10. Centering the Strike to Reduce Oscillations

진동을 줄이기 위한 방법으로 **짝수 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BN)에서 행사가격 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BK)가 만기 노드의 중앙에 정확히 놓이도록** 상승·하락 계수를 조정한다. 아래 그림은 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BN%3D4)의 예로, 강조된 만기 중앙 노드가 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BK%3DS_0u%5E2d%5E2)를 만족하도록 만든 구조다.

![strike centered binomial tree]({{ "/assets/img/cf-post4-strike-centered-binomial-tree.svg" | relative_url }})

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BN)이 짝수이면 만기 중앙 노드는 상승 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BN%2F2)번, 하락 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BN%2F2)번을 거친 가격이므로


![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cboxed%7BS_0u%5E%7BN%2F2%7Dd%5E%7BN%2F2%7D%3DK.%7D)


표준 CRR의 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3Bu%3De%5E%7B%5Csigma%5Csqrt%7B%5CDelta%20t%7D%7D), ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3Bd%3De%5E%7B-%5Csigma%5Csqrt%7B%5CDelta%20t%7D%7D)를 그대로 쓰면 중앙 노드는 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BS_0)가 된다. 따라서 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BK%5Cne%20S_0)인 경우 중앙 노드와 행사가격이 어긋난다.

이를 조정하기 위해


![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cboxed%7Bu%3D%5Calpha%20e%5E%7B%5Csigma%5Csqrt%7B%5CDelta%20t%7D%7D%2C%5Cqquad%0Ad%3D%5Calpha%20e%5E%7B-%5Csigma%5Csqrt%7B%5CDelta%20t%7D%7D%7D)


로 놓는다. 그러면


![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cbegin%7Baligned%7D%0AS_0u%5E%7BN%2F2%7Dd%5E%7BN%2F2%7D%0A%26%3DS_0%5Cleft%28%5Calpha%20e%5E%7B%5Csigma%5Csqrt%7B%5CDelta%20t%7D%7D%5Cright%29%5E%7BN%2F2%7D%0A%5Cleft%28%5Calpha%20e%5E%7B-%5Csigma%5Csqrt%7B%5CDelta%20t%7D%7D%5Cright%29%5E%7BN%2F2%7D%5C%5C%0A%26%3DS_0%5Calpha%5E%7BN%2F2%2BN%2F2%7D%0A%20e%5E%7B%28N%2F2%29%5Csigma%5Csqrt%7B%5CDelta%20t%7D-%28N%2F2%29%5Csigma%5Csqrt%7B%5CDelta%20t%7D%7D%5C%5C%0A%26%3DS_0%5Calpha%5EN.%0A%5Cend%7Baligned%7D)


이 값을 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BK)와 같게 만들면


![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BS_0%5Calpha%5EN%3DK%0A%5Cquad%5CLongrightarrow%5Cquad%0A%5Cboxed%7B%5Calpha%3D%5Cleft%28%5Cfrac%20K%7BS_0%7D%5Cright%29%5E%7B1%2FN%7D%7D.)


로그를 취하면


![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cboxed%7B%5Clog%5Calpha%3D%5Cfrac1N%5Clog%5Cfrac%20K%7BS_0%7D%7D.)


따라서 조정된 계수는


![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cbegin%7Baligned%7D%0Au%26%3D%5Cexp%5Cleft%28%5Cfrac1N%5Clog%5Cfrac%20K%7BS_0%7D%2B%5Csigma%5Csqrt%7B%5CDelta%20t%7D%5Cright%29%2C%5C%5C%0Ad%26%3D%5Cexp%5Cleft%28%5Cfrac1N%5Clog%5Cfrac%20K%7BS_0%7D-%5Csigma%5Csqrt%7B%5CDelta%20t%7D%5Cright%29.%0A%5Cend%7Baligned%7D)


![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5CDelta%20t%3DT%2FN)을 대입하면 같은 식을


![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cboxed%7B%5Cdisplaystyle%0Au%3D%5Cexp%5Cleft%28%5Cfrac%7B%5CDelta%20t%7D%7BT%7D%5Clog%5Cfrac%20K%7BS_0%7D%2B%5Csigma%5Csqrt%7B%5CDelta%20t%7D%5Cright%29%2C%5Cqquad%0Ad%3D%5Cexp%5Cleft%28%5Cfrac%7B%5CDelta%20t%7D%7BT%7D%5Clog%5Cfrac%20K%7BS_0%7D-%5Csigma%5Csqrt%7B%5CDelta%20t%7D%5Cright%29.%7D)


로 쓸 수 있다. 이 계수를 사용한 뒤에도 위험중립 확률은


![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3Bp%3D%5Cfrac%7Be%5E%7Br%5CDelta%20t%7D-d%7D%7Bu-d%7D)


로 계산하며, **![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B0%3Cp%3C1)인지 확인**해야 한다. 행사가격을 격자에 맞추는 조정은 진동을 완화하기 위한 방법이지만, 모든 설정에서 진동을 완전히 제거한다는 뜻은 아니다.

## 11. Key Connections

- **Likelihood-ratio Delta:** 지급액을 직접 미분하지 않고 확률밀도를 미분하여 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3Bv%27%28S%29)를 기대값으로 표현한다.
- **Monte Carlo variance:** 공통난수와 Malliavin 가중치는 모두 Delta를 계산하는 방법이지만, 어느 쪽의 분산이 작은지는 지급액에 따라 달라진다.
- **Portfolio immunisation:** 편미분으로 측정하는 민감도를 줄여 포트폴리오의 위험요인 노출을 관리한다.
- **Moment matching:** 이항트리 수익률의 첫 두 모멘트를 Black–Scholes 수익률에 맞추는 조건을 세운다.
- **Binomial recursion:** 만기 지급액에서 시작하여 할인된 위험중립 기대값으로 이전 노드의 가격을 계산한다.
- **Strike alignment:** 콜옵션의 꺾이는 지점 ![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3BK)를 만기 격자의 중앙 노드에 맞춰 가격 진동을 완화한다.
