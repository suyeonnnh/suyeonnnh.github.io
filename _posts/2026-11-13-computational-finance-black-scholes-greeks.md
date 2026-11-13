---
layout: post
title: "Computational Finance : Black–Scholes Replication and Monte Carlo Greeks"
date: 2026-11-13 13:00:00 +0100
categories:
  - Computational Finance
tags:
  - Black-Scholes
  - Replicating Portfolio
  - Delta Hedging
  - Monte Carlo
  - Finite Difference
  - Common Random Numbers
  - Greeks
math: false
---

이 글은 Black–Scholes 모형에서 파생상품의 가격과 복제 포트폴리오를 구하고, Monte Carlo로 가격의 민감도인 Delta를 계산하는 방법을 정리한다. 유한차분, common random numbers, 경로별 미분(pathwise differentiation), likelihood-ratio 방식으로 이어지는 계산을 다룬다.

---

## 1. Black–Scholes Pricing and Replicating Portfolio

만기 `T`에 지급되는 payoff를 `h(S_T)`라고 하자. 시점 `t`의 파생상품 가격은 `v(t,S_t)`이며, 위험중립측도 `Q`에서

<div style="clear:both; overflow-x:auto; text-align:center; margin:20px 0 24px; line-height:1.2;">
  <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%20v%28t%2CS%29%3D%5Cmathbb%7BE%7D_%7B%5Cmathbb%7BQ%7D%7D%5Cleft%5Be%5E%7B-r%28T-t%29%7Dh%28S_T%29%5Cmid%20S_t%3DS%5Cright%5D" alt="2026-11-13-computational-finance-black-scholes-greeks equation 1" style="display:block; width:auto; max-width:100%; height:auto; margin:0 auto;" />
</div>

무위험 채권 `B_t`와 주식 `S_t`는 다음과 같이 주어진다.

<div style="clear:both; overflow-x:auto; text-align:center; margin:20px 0 24px; line-height:1.2;">
  <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%20B_t%3De%5E%7Brt%7D%2C%5Cqquad%20S_t%3DS_0%5Cexp%5Cleft%28%5Cleft%28r-%5Cfrac%7B%5Csigma%5E2%7D%7B2%7D%5Cright%29t%2B%5Csigma%20W_t%5Cright%29" alt="2026-11-13-computational-finance-black-scholes-greeks equation 2" style="display:block; width:auto; max-width:100%; height:auto; margin:0 auto;" />
</div>

주식 `φ_t`주와 채권 `ψ_t`단위로 이루어진 포트폴리오의 가치는

<div style="clear:both; overflow-x:auto; text-align:center; margin:20px 0 24px; line-height:1.2;">
  <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%20V_t%3D%5Cvarphi_t%20S_t%2B%5Cpsi_t%20B_t" alt="2026-11-13-computational-finance-black-scholes-greeks equation 3" style="display:block; width:auto; max-width:100%; height:auto; margin:0 auto;" />
</div>

이다. 복제 포트폴리오는 만기에

<div style="clear:both; overflow-x:auto; text-align:center; margin:20px 0 24px; line-height:1.2;">
  <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%20V_T%3Dh%28S_T%29" alt="2026-11-13-computational-finance-black-scholes-greeks equation 4" style="display:block; width:auto; max-width:100%; height:auto; margin:0 auto;" />
</div>

를 만족한다. 자기금융(self-financing) 조건은

<div style="clear:both; overflow-x:auto; text-align:center; margin:20px 0 24px; line-height:1.2;">
  <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%20dV_t%3D%5Cvarphi_t%5C%2CdS_t%2B%5Cpsi_t%5C%2CdB_t" alt="2026-11-13-computational-finance-black-scholes-greeks equation 5" style="display:block; width:auto; max-width:100%; height:auto; margin:0 auto;" />
</div>

이며, 복제 전략의 보유량은

<div style="clear:both; overflow-x:auto; text-align:center; margin:20px 0 24px; line-height:1.2;">
  <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%20%5Cboxed%7B%5Cvarphi_t%3D%5Cfrac%7B%5Cpartial%20v%7D%7B%5Cpartial%20S%7D%28t%2CS_t%29%7D" alt="2026-11-13-computational-finance-black-scholes-greeks equation 6" style="display:block; width:auto; max-width:100%; height:auto; margin:0 auto;" />
</div>

<div style="clear:both; overflow-x:auto; text-align:center; margin:20px 0 24px; line-height:1.2;">
  <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%20%5Cboxed%7B%5Cpsi_t%3D%5Cfrac%7Bv%28t%2CS_t%29-%5Cvarphi_t%20S_t%7D%7BB_t%7D%7D" alt="2026-11-13-computational-finance-black-scholes-greeks equation 7" style="display:block; width:auto; max-width:100%; height:auto; margin:0 auto;" />
</div>

로 표현된다. 주식 보유량 `φ_t`는 Delta에 해당한다.

## 2. Monte Carlo Approximation of Delta

시점 0의 가격을 초기 주가 `S`의 함수로 쓰면

<div style="clear:both; overflow-x:auto; text-align:center; margin:20px 0 24px; line-height:1.2;">
  <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%20v%28S%29%3D%5Cmathbb%7BE%7D_%7B%5Cmathbb%7BQ%7D%7D%5Cleft%5Be%5E%7B-rT%7Dh%5Cleft%28S%5Cexp%5Cleft%28%5Cleft%28r-%5Cfrac%7B%5Csigma%5E2%7D%7B2%7D%5Cright%29T%2B%5Csigma%20W_T%5Cright%29%5Cright%29%5Cright%5D" alt="2026-11-13-computational-finance-black-scholes-greeks equation 8" style="display:block; width:auto; max-width:100%; height:auto; margin:0 auto;" />
</div>

계산할 민감도는

<div style="clear:both; overflow-x:auto; text-align:center; margin:20px 0 24px; line-height:1.2;">
  <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%20%5CDelta%3D%5Cfrac%7B%5Cpartial%20v%7D%7B%5Cpartial%20S%7D%28S%29%3Dv%27%28S%29" alt="2026-11-13-computational-finance-black-scholes-greeks equation 9" style="display:block; width:auto; max-width:100%; height:auto; margin:0 auto;" />
</div>

이다. `v(S)`를 Monte Carlo로 추정할 수 있으므로, 유한차분으로 도함수를 근사할 수 있다.

### 2.1 Forward Difference

<div style="clear:both; overflow-x:auto; text-align:center; margin:20px 0 24px; line-height:1.2;">
  <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%20%5Cboxed%7Bv%27%28S%29%5Capprox%5Cfrac%7Bv%28S%2B%5CDelta%20S%29-v%28S%29%7D%7B%5CDelta%20S%7D%7D" alt="2026-11-13-computational-finance-black-scholes-greeks equation 10" style="display:block; width:auto; max-width:100%; height:auto; margin:0 auto;" />
</div>

Taylor 전개를 이용하면

<div style="clear:both; overflow-x:auto; text-align:center; margin:20px 0 24px; line-height:1.2;">
  <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%20v%28S%2B%5CDelta%20S%29%3Dv%28S%29%2Bv%27%28S%29%5CDelta%20S%2B%5Cfrac12v%27%27%28S%29%28%5CDelta%20S%29%5E2%2BO%28%28%5CDelta%20S%29%5E3%29" alt="2026-11-13-computational-finance-black-scholes-greeks equation 11" style="display:block; width:auto; max-width:100%; height:auto; margin:0 auto;" />
</div>

따라서

<div style="clear:both; overflow-x:auto; text-align:center; margin:20px 0 24px; line-height:1.2;">
  <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%20%5Cfrac%7Bv%28S%2B%5CDelta%20S%29-v%28S%29%7D%7B%5CDelta%20S%7D%3Dv%27%28S%29%2B%5Cfrac12v%27%27%28S%29%5CDelta%20S%2BO%28%28%5CDelta%20S%29%5E2%29" alt="2026-11-13-computational-finance-black-scholes-greeks equation 12" style="display:block; width:auto; max-width:100%; height:auto; margin:0 auto;" />
</div>

Forward difference의 절단오차는 일반적으로 `O(ΔS)`이다.

### 2.2 Central Difference

<div style="clear:both; overflow-x:auto; text-align:center; margin:20px 0 24px; line-height:1.2;">
  <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%20%5Cboxed%7Bv%27%28S%29%5Capprox%5Cfrac%7Bv%28S%2B%5CDelta%20S%29-v%28S-%5CDelta%20S%29%7D%7B2%5CDelta%20S%7D%7D" alt="2026-11-13-computational-finance-black-scholes-greeks equation 13" style="display:block; width:auto; max-width:100%; height:auto; margin:0 auto;" />
</div>

양쪽 Taylor 전개는

<div style="clear:both; overflow-x:auto; text-align:center; margin:20px 0 24px; line-height:1.2;">
  <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%20v%28S%2B%5CDelta%20S%29%3Dv%28S%29%2Bv%27%28S%29%5CDelta%20S%2B%5Cfrac12v%27%27%28S%29%28%5CDelta%20S%29%5E2%2B%5Cfrac16v%27%27%27%28S%29%28%5CDelta%20S%29%5E3%2BO%28%28%5CDelta%20S%29%5E4%29%2C" alt="2026-11-13-computational-finance-black-scholes-greeks equation 14" style="display:block; width:auto; max-width:100%; height:auto; margin:0 auto;" />
</div>

<div style="clear:both; overflow-x:auto; text-align:center; margin:20px 0 24px; line-height:1.2;">
  <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%20v%28S-%5CDelta%20S%29%3Dv%28S%29-v%27%28S%29%5CDelta%20S%2B%5Cfrac12v%27%27%28S%29%28%5CDelta%20S%29%5E2-%5Cfrac16v%27%27%27%28S%29%28%5CDelta%20S%29%5E3%2BO%28%28%5CDelta%20S%29%5E4%29" alt="2026-11-13-computational-finance-black-scholes-greeks equation 15" style="display:block; width:auto; max-width:100%; height:auto; margin:0 auto;" />
</div>

두 식을 빼고 `2ΔS`로 나누면

<div style="clear:both; overflow-x:auto; text-align:center; margin:20px 0 24px; line-height:1.2;">
  <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%20%5Cfrac%7Bv%28S%2B%5CDelta%20S%29-v%28S-%5CDelta%20S%29%7D%7B2%5CDelta%20S%7D%3Dv%27%28S%29%2B%5Cfrac16v%27%27%27%28S%29%28%5CDelta%20S%29%5E2%2BO%28%28%5CDelta%20S%29%5E3%29" alt="2026-11-13-computational-finance-black-scholes-greeks equation 16" style="display:block; width:auto; max-width:100%; height:auto; margin:0 auto;" />
</div>

충분히 매끄러운 `v`에 대해 central difference의 절단오차는 `O((ΔS)²)`이다.

## 3. Choosing the Finite-Difference Step

`ΔS`를 줄이면 절단오차가 작아진다. 반면 Monte Carlo 추정값을 서로 빼고 작은 `ΔS`로 나누는 과정에서는 표본오차가 확대될 수 있다. 매우 작은 간격에서는 부동소수점 연산의 정밀도 문제도 발생한다.

필요한 정확도에 맞추어 간격을 정한다. 한 가지 선택은

<div style="clear:both; overflow-x:auto; text-align:center; margin:20px 0 24px; line-height:1.2;">
  <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%20%5Cboxed%7B%5CDelta%20S%3DS%5Cvarepsilon%5E%7B1%2F4%7D%7D" alt="2026-11-13-computational-finance-black-scholes-greeks equation 17" style="display:block; width:auto; max-width:100%; height:auto; margin:0 auto;" />
</div>

이다. 여기서 `ε`은 계산에서 고려하는 오차 수준을 나타낸다. 실제 간격은 표본오차와 절단오차의 크기에 따라 조정한다.

## 4. Common Random Numbers for Delta Estimation

Monte Carlo 가격 추정량을 `v̂_MC(S)`라고 하면 central-difference Delta 추정량은

<div style="clear:both; overflow-x:auto; text-align:center; margin:20px 0 24px; line-height:1.2;">
  <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%20%5Cwidehat%5CDelta_%7B%5Cmathrm%7BFD%7D%7D%3D%5Cfrac%7B%5Cwidehat%20v_%7B%5Cmathrm%7BMC%7D%7D%28S%2B%5CDelta%20S%29-%5Cwidehat%20v_%7B%5Cmathrm%7BMC%7D%7D%28S-%5CDelta%20S%29%7D%7B2%5CDelta%20S%7D" alt="2026-11-13-computational-finance-black-scholes-greeks equation 18" style="display:block; width:auto; max-width:100%; height:auto; margin:0 auto;" />
</div>

두 가격을 계산할 때 **동일한 난수**를 사용한다. `Z_j ~ N(0,1)`을 생성하고

<div style="clear:both; overflow-x:auto; text-align:center; margin:20px 0 24px; line-height:1.2;">
  <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%20Y_j%3D%5Cexp%5Cleft%28%5Cleft%28r-%5Cfrac%7B%5Csigma%5E2%7D%7B2%7D%5Cright%29T%2B%5Csigma%5Csqrt%20T%5C%2CZ_j%5Cright%29" alt="2026-11-13-computational-finance-black-scholes-greeks equation 19" style="display:block; width:auto; max-width:100%; height:auto; margin:0 auto;" />
</div>

로 두면

<div style="clear:both; overflow-x:auto; text-align:center; margin:20px 0 24px; line-height:1.2;">
  <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%20%5Cboxed%7B%5Cwidehat%5CDelta_%7B%5Cmathrm%7BCRN%7D%7D%3D%5Cfrac%7Be%5E%7B-rT%7D%7D%7BN%7D%5Csum_%7Bj%3D1%7D%5EN%5Cfrac%7Bh%28%28S%2B%5CDelta%20S%29Y_j%29-h%28%28S-%5CDelta%20S%29Y_j%29%7D%7B2%5CDelta%20S%7D%7D" alt="2026-11-13-computational-finance-black-scholes-greeks equation 20" style="display:block; width:auto; max-width:100%; height:auto; margin:0 auto;" />
</div>

이다. 각 `j`에 대해 두 payoff에 같은 `Y_j`를 사용한다.

일반적으로 두 가격 추정량 `A, B`에 대해

<div style="clear:both; overflow-x:auto; text-align:center; margin:20px 0 24px; line-height:1.2;">
  <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%20%5Coperatorname%7BVar%7D%28A-B%29%3D%5Coperatorname%7BVar%7D%28A%29%2B%5Coperatorname%7BVar%7D%28B%29-2%5Coperatorname%7BCov%7D%28A%2CB%29" alt="2026-11-13-computational-finance-black-scholes-greeks equation 21" style="display:block; width:auto; max-width:100%; height:auto; margin:0 auto;" />
</div>

공통 난수로 양의 공분산을 유도할 수 있으면 차이의 분산이 감소한다. 동일한 계산 횟수에서도 난수 구성에 따라 두 추정량의 분산이 달라질 수 있다.

## 5. Lemma A: Pathwise Differentiation

다음과 같이 정의한다.

<div style="clear:both; overflow-x:auto; text-align:center; margin:20px 0 24px; line-height:1.2;">
  <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%20Y%3D%5Cexp%5Cleft%28%5Cleft%28r-%5Cfrac%7B%5Csigma%5E2%7D%7B2%7D%5Cright%29T%2B%5Csigma%20W_T%5Cright%29%2C%5Cqquad%20v%28S%29%3D%5Cmathbb%20E_%7B%5Cmathbb%20Q%7D%5Be%5E%7B-rT%7Dh%28SY%29%5D" alt="2026-11-13-computational-finance-black-scholes-greeks equation 22" style="display:block; width:auto; max-width:100%; height:auto; margin:0 auto;" />
</div>

`h ∈ C¹`이고 `h′`가 유계라면 미분과 기댓값을 교환하여

<div style="clear:both; overflow-x:auto; text-align:center; margin:20px 0 24px; line-height:1.2;">
  <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%20%5Cboxed%7Bv%27%28S%29%3D%5Cmathbb%20E_%7B%5Cmathbb%20Q%7D%5Cleft%5Be%5E%7B-rT%7Dh%27%28SY%29Y%5Cright%5D%7D" alt="2026-11-13-computational-finance-black-scholes-greeks equation 23" style="display:block; width:auto; max-width:100%; height:auto; margin:0 auto;" />
</div>

를 얻는다.

### 5.1 Proof Using Dominated Convergence

차분몫으로 시작한다.

<div style="clear:both; overflow-x:auto; text-align:center; margin:20px 0 24px; line-height:1.2;">
  <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%20v%27%28S%29%3D%5Clim_%7B%5Cepsilon%5Cto0%7D%5Cfrac%7Bv%28S%2B%5Cepsilon%29-v%28S%29%7D%7B%5Cepsilon%7D" alt="2026-11-13-computational-finance-black-scholes-greeks equation 24" style="display:block; width:auto; max-width:100%; height:auto; margin:0 auto;" />
</div>

<div style="clear:both; overflow-x:auto; text-align:center; margin:20px 0 24px; line-height:1.2;">
  <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%20%3D%5Clim_%7B%5Cepsilon%5Cto0%7D%5Cmathbb%20E_%7B%5Cmathbb%20Q%7D%5Cleft%5Be%5E%7B-rT%7D%5Cfrac%7Bh%28%28S%2B%5Cepsilon%29Y%29-h%28SY%29%7D%7B%5Cepsilon%7D%5Cright%5D" alt="2026-11-13-computational-finance-black-scholes-greeks equation 25" style="display:block; width:auto; max-width:100%; height:auto; margin:0 auto;" />
</div>

각 경로에서 `Y`를 고정하면 연쇄법칙으로

<div style="clear:both; overflow-x:auto; text-align:center; margin:20px 0 24px; line-height:1.2;">
  <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%20%5Clim_%7B%5Cepsilon%5Cto0%7D%5Cfrac%7Bh%28%28S%2B%5Cepsilon%29Y%29-h%28SY%29%7D%7B%5Cepsilon%7D%3Dh%27%28SY%29Y" alt="2026-11-13-computational-finance-black-scholes-greeks equation 26" style="display:block; width:auto; max-width:100%; height:auto; margin:0 auto;" />
</div>

평균값정리에 의해

<div style="clear:both; overflow-x:auto; text-align:center; margin:20px 0 24px; line-height:1.2;">
  <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%20%5Cleft%7C%5Cfrac%7Bh%28%28S%2B%5Cepsilon%29Y%29-h%28SY%29%7D%7B%5Cepsilon%7D%5Cright%7C%5Cleq%5ClVert%20h%27%5CrVert_%5Cinfty%20Y" alt="2026-11-13-computational-finance-black-scholes-greeks equation 27" style="display:block; width:auto; max-width:100%; height:auto; margin:0 auto;" />
</div>

또한

<div style="clear:both; overflow-x:auto; text-align:center; margin:20px 0 24px; line-height:1.2;">
  <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%20%5Cmathbb%20E%5BY%5D%3De%5E%7BrT%7D%3C%5Cinfty" alt="2026-11-13-computational-finance-black-scholes-greeks equation 28" style="display:block; width:auto; max-width:100%; height:auto; margin:0 auto;" />
</div>

따라서 dominated convergence theorem을 적용하여

<div style="clear:both; overflow-x:auto; text-align:center; margin:20px 0 24px; line-height:1.2;">
  <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%20%5Cboxed%7Bv%27%28S%29%3D%5Cmathbb%20E_%7B%5Cmathbb%20Q%7D%5Be%5E%7B-rT%7Dh%27%28SY%29Y%5D%7D" alt="2026-11-13-computational-finance-black-scholes-greeks equation 29" style="display:block; width:auto; max-width:100%; height:auto; margin:0 auto;" />
</div>

를 얻는다. 이를 이용한 Monte Carlo 추정량은

<div style="clear:both; overflow-x:auto; text-align:center; margin:20px 0 24px; line-height:1.2;">
  <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%20%5Cwidehat%5CDelta_%7B%5Cmathrm%7BPW%7D%7D%3D%5Cfrac%7Be%5E%7B-rT%7D%7DN%5Csum_%7Bj%3D1%7D%5EN%20h%27%28SY_j%29Y_j" alt="2026-11-13-computational-finance-black-scholes-greeks equation 30" style="display:block; width:auto; max-width:100%; height:auto; margin:0 auto;" />
</div>

## 6. Lemma B: Likelihood-Ratio Differentiation

`h`가 유계이고 측정가능한 경우에는 정규분포 밀도를 미분하는 방식으로 Delta를 구할 수 있다.

<div style="clear:both; overflow-x:auto; text-align:center; margin:20px 0 24px; line-height:1.2;">
  <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%20%5Cboxed%7Bv%27%28S%29%3D%5Cfrac%7Be%5E%7B-rT%7D%7D%7BS%5Csigma%20T%7D%5Cmathbb%20E_%7B%5Cmathbb%20Q%7D%5Cleft%5Bh%28SY%29W_T%5Cright%5D%7D" alt="2026-11-13-computational-finance-black-scholes-greeks equation 31" style="display:block; width:auto; max-width:100%; height:auto; margin:0 auto;" />
</div>

이 방식은 likelihood-ratio method 또는 score-function method로 불린다.

### 6.1 Proof Using the Lognormal Density

`W_T ~ N(0,T)`이므로 `U = log(SY)`는

<div style="clear:both; overflow-x:auto; text-align:center; margin:20px 0 24px; line-height:1.2;">
  <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%20U%5Csim%5Cmathcal%20N%5Cleft%28%5Clog%20S%2B%5Cleft%28r-%5Cfrac%7B%5Csigma%5E2%7D%7B2%7D%5Cright%29T%2C%5Csigma%5E2T%5Cright%29" alt="2026-11-13-computational-finance-black-scholes-greeks equation 32" style="display:block; width:auto; max-width:100%; height:auto; margin:0 auto;" />
</div>

따라서

<div style="clear:both; overflow-x:auto; text-align:center; margin:20px 0 24px; line-height:1.2;">
  <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%20v%28S%29%3De%5E%7B-rT%7D%5Cint_%7B-%5Cinfty%7D%5E%7B%5Cinfty%7Dh%28e%5Eu%29%5Cfrac%7B1%7D%7B%5Csqrt%7B2%5Cpi%7D%5Csigma%5Csqrt%20T%7D%5Cexp%5Cleft%28-%5Cfrac%7B%28u-%5Clog%20S-%28r-%5Csigma%5E2%2F2%29T%29%5E2%7D%7B2%5Csigma%5E2T%7D%5Cright%29du" alt="2026-11-13-computational-finance-black-scholes-greeks equation 33" style="display:block; width:auto; max-width:100%; height:auto; margin:0 auto;" />
</div>

밀도의 지수 부분을 `S`에 대해 미분하면

<div style="clear:both; overflow-x:auto; text-align:center; margin:20px 0 24px; line-height:1.2;">
  <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%20%5Cfrac%7B%5Cpartial%7D%7B%5Cpartial%20S%7D%5Cleft%5B-%5Cfrac%7B%28u-%5Clog%20S-%28r-%5Csigma%5E2%2F2%29T%29%5E2%7D%7B2%5Csigma%5E2T%7D%5Cright%5D%3D%5Cfrac%7Bu-%5Clog%20S-%28r-%5Csigma%5E2%2F2%29T%7D%7BS%5Csigma%5E2T%7D" alt="2026-11-13-computational-finance-black-scholes-greeks equation 34" style="display:block; width:auto; max-width:100%; height:auto; margin:0 auto;" />
</div>

변수변환

<div style="clear:both; overflow-x:auto; text-align:center; margin:20px 0 24px; line-height:1.2;">
  <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%20z%3D%5Cfrac%7Bu-%5Clog%20S-%28r-%5Csigma%5E2%2F2%29T%7D%7B%5Csigma%5Csqrt%20T%7D%2C%5Cqquad%20du%3D%5Csigma%5Csqrt%20T%5C%2Cdz" alt="2026-11-13-computational-finance-black-scholes-greeks equation 35" style="display:block; width:auto; max-width:100%; height:auto; margin:0 auto;" />
</div>

를 사용하면

<div style="clear:both; overflow-x:auto; text-align:center; margin:20px 0 24px; line-height:1.2;">
  <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%20v%27%28S%29%3D%5Cfrac%7Be%5E%7B-rT%7D%7D%7BS%5Csigma%5Csqrt%20T%7D%5Cint_%7B-%5Cinfty%7D%5E%7B%5Cinfty%7Dh%5Cleft%28S%20e%5E%7B%28r-%5Csigma%5E2%2F2%29T%2B%5Csigma%5Csqrt%20Tz%7D%5Cright%29z%5Cfrac%7Be%5E%7B-z%5E2%2F2%7D%7D%7B%5Csqrt%7B2%5Cpi%7D%7Ddz" alt="2026-11-13-computational-finance-black-scholes-greeks equation 36" style="display:block; width:auto; max-width:100%; height:auto; margin:0 auto;" />
</div>

`Z = W_T/√T`이므로

<div style="clear:both; overflow-x:auto; text-align:center; margin:20px 0 24px; line-height:1.2;">
  <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%20%5Cboxed%7Bv%27%28S%29%3D%5Cfrac%7Be%5E%7B-rT%7D%7D%7BS%5Csigma%20T%7D%5Cmathbb%20E_%7B%5Cmathbb%20Q%7D%5Bh%28SY%29W_T%5D%7D" alt="2026-11-13-computational-finance-black-scholes-greeks equation 37" style="display:block; width:auto; max-width:100%; height:auto; margin:0 auto;" />
</div>

Monte Carlo 추정량은

<div style="clear:both; overflow-x:auto; text-align:center; margin:20px 0 24px; line-height:1.2;">
  <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%20%5Cwidehat%5CDelta_%7B%5Cmathrm%7BLR%7D%7D%3D%5Cfrac%7Be%5E%7B-rT%7D%7D%7BNS%5Csigma%5Csqrt%20T%7D%5Csum_%7Bj%3D1%7D%5EN%20h%28SY_j%29Z_j" alt="2026-11-13-computational-finance-black-scholes-greeks equation 38" style="display:block; width:auto; max-width:100%; height:auto; margin:0 auto;" />
</div>


## 7. Discontinuous Payoffs and Other Greeks

### 7.1 Digital or Binary Option

Digital option은 만기 주가가 행사가격 K 이상이면 1을 지급하고, K보다 작으면 0을 지급하는 파생상품이다.

Payoff function은 다음과 같다.

<div style="clear:both; overflow-x:auto; text-align:center; margin:20px 0 24px; line-height:1.2;">
  <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%20h%28x%29%3D%5Cbegin%7Bcases%7D1%2C%26x%5Cgeq%20K%2C%5C%5C0%2C%26x%3CK.%5Cend%7Bcases%7D" alt="2026-11-13-computational-finance-black-scholes-greeks equation 39" style="display:block; width:auto; max-width:100%; height:auto; margin:0 auto;" />
</div>

즉, 주가가 K를 넘는 순간 payoff가 0에서 1로 바뀐다. 따라서 K에서 jump discontinuity가 발생한다.

이 옵션의 가격은 위험중립측도에서 다음과 같이 표현된다.

<div style="clear:both; overflow-x:auto; text-align:center; margin:20px 0 24px; line-height:1.2;">
  <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%20v%28S%29%3De%5E%7B-rT%7D%5Cmathbb%7BE%7D_%7B%5Cmathbb%20Q%7D%5B%5Cmathbf%7B1%7D_%7B%5C%7BSY%5Cgeq%20K%5C%7D%7D%5D%3De%5E%7B-rT%7D%5Cmathbb%7BQ%7D%28SY%5Cgeq%20K%29" alt="2026-11-13-computational-finance-black-scholes-greeks equation 40" style="display:block; width:auto; max-width:100%; height:auto; margin:0 auto;" />
</div>

여기서 Y는 Black–Scholes 모형의 주가 변화율을 나타낸다.

Digital option의 가격은 만기에 주가가 행사가격 이상일 위험중립확률에 할인계수를 곱한 값이다.

### 7.2 Why Pathwise Differentiation Fails

앞서 Lemma A에서 사용한 pathwise differentiation은 payoff function이 연속미분가능하고 도함수가 유계라는 조건을 필요로 한다.

그러나 digital option의 payoff는 K에서 불연속이다.

K를 제외한 모든 지점에서 payoff가 상수이므로

<div style="clear:both; overflow-x:auto; text-align:center; margin:20px 0 24px; line-height:1.2;">
  <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%20h%27%28x%29%3D0%5Cquad%20%28x%5Cneq%20K%29" alt="2026-11-13-computational-finance-black-scholes-greeks equation 41" style="display:block; width:auto; max-width:100%; height:auto; margin:0 auto;" />
</div>

이다.

만약 pathwise 공식을 그대로 적용한다면

<div style="clear:both; overflow-x:auto; text-align:center; margin:20px 0 24px; line-height:1.2;">
  <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%20%5Cwidehat%7B%5CDelta%7D_%7B%5Cmathrm%7BPW%7D%7D%3D%5Cfrac%7Be%5E%7B-rT%7D%7D%7BN%7D%5Csum_%7Bj%3D1%7D%5E%7BN%7Dh%27%28SY_j%29Y_j" alt="2026-11-13-computational-finance-black-scholes-greeks equation 42" style="display:block; width:auto; max-width:100%; height:auto; margin:0 auto;" />
</div>

에서 거의 모든 표본의 도함수가 0이므로 Delta 추정값도 0이 된다.

하지만 초기 주가 S가 증가하면 만기 주가가 K 이상일 확률도 증가한다. 따라서 digital option의 가격은 S에 따라 변하고 Delta도 일반적으로 0보다 크다.

이러한 차이가 발생하는 이유는 불연속 payoff에 대해 미분과 기댓값을 교환할 수 있는 조건이 성립하지 않기 때문이다.

### 7.3 Likelihood-Ratio Method

Likelihood-ratio method는 payoff function을 직접 미분하는 대신 확률밀도를 미분한다.

Lemma B에서 얻은 공식은 다음과 같다.

<div style="clear:both; overflow-x:auto; text-align:center; margin:20px 0 24px; line-height:1.2;">
  <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%20%5CDelta%3D%5Cfrac%7Be%5E%7B-rT%7D%7D%7BS%5Csigma%20T%7D%5Cmathbb%7BE%7D_%7B%5Cmathbb%20Q%7D%5Bh%28SY%29W_T%5D" alt="2026-11-13-computational-finance-black-scholes-greeks equation 43" style="display:block; width:auto; max-width:100%; height:auto; margin:0 auto;" />
</div>

이 공식에는 h'(x)가 포함되지 않는다. 따라서 payoff가 불연속이어도 Lemma B의 유계성과 측정가능성 조건을 만족하면 적용할 수 있다.

Digital option의 payoff를 대입하면 Monte Carlo 추정량은

<div style="clear:both; overflow-x:auto; text-align:center; margin:20px 0 24px; line-height:1.2;">
  <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%20%5Cwidehat%7B%5CDelta%7D_%7B%5Cmathrm%7BLR%7D%7D%3D%5Cfrac%7Be%5E%7B-rT%7D%7D%7BNS%5Csigma%5Csqrt%7BT%7D%7D%5Csum_%7Bj%3D1%7D%5E%7BN%7D%5Cmathbf%7B1%7D_%7B%5C%7BSY_j%5Cgeq%20K%5C%7D%7DZ_j" alt="2026-11-13-computational-finance-black-scholes-greeks equation 44" style="display:block; width:auto; max-width:100%; height:auto; margin:0 auto;" />
</div>

이다.

각 표본에서 만기 주가가 K 이상인지 확인하고, 해당 payoff에 표준정규변수 Z_j를 곱한 뒤 평균을 계산한다.

이 방식은 유한차분 간격 ΔS를 선택할 필요가 없다.

### 7.4 Analytical Delta of a Digital Option

Black–Scholes 모형에서는 digital option의 Delta를 해석적으로도 계산할 수 있다.

다음과 같이 정의하자.

<div style="clear:both; overflow-x:auto; text-align:center; margin:20px 0 24px; line-height:1.2;">
  <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%20d_2%3D%5Cfrac%7B%5Clog%28S%2FK%29%2B%28r-%5Csigma%5E2%2F2%29T%7D%7B%5Csigma%5Csqrt%7BT%7D%7D" alt="2026-11-13-computational-finance-black-scholes-greeks equation 45" style="display:block; width:auto; max-width:100%; height:auto; margin:0 auto;" />
</div>

그러면 digital option의 가격은

<div style="clear:both; overflow-x:auto; text-align:center; margin:20px 0 24px; line-height:1.2;">
  <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%20v%28S%29%3De%5E%7B-rT%7D%5CPhi%28d_2%29" alt="2026-11-13-computational-finance-black-scholes-greeks equation 46" style="display:block; width:auto; max-width:100%; height:auto; margin:0 auto;" />
</div>

이다.

여기서 Φ는 표준정규분포의 누적분포함수이다.

이를 S에 대해 미분하면

<div style="clear:both; overflow-x:auto; text-align:center; margin:20px 0 24px; line-height:1.2;">
  <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%20%5CDelta%3D%5Cfrac%7Be%5E%7B-rT%7D%5Cphi%28d_2%29%7D%7BS%5Csigma%5Csqrt%7BT%7D%7D" alt="2026-11-13-computational-finance-black-scholes-greeks equation 47" style="display:block; width:auto; max-width:100%; height:auto; margin:0 auto;" />
</div>

를 얻는다.

여기서 φ는 표준정규분포의 확률밀도함수이다.

따라서 digital option은 payoff가 불연속이더라도 가격 함수는 미분가능하며, Delta를 계산할 수 있다.

이 해석적 결과는 Monte Carlo로 구한 likelihood-ratio 추정값을 검증하는 기준으로도 사용할 수 있다.

### 7.5 Other Greeks: Gamma and Vega

Delta 외에도 파생상품 가격의 다양한 민감도를 계산할 수 있다.

**Gamma**는 주가 변화에 따른 Delta의 변화율이다. 가격 함수를 주가에 대해 두 번 미분한다.

<div style="clear:both; overflow-x:auto; text-align:center; margin:20px 0 24px; line-height:1.2;">
  <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%20%5CGamma%3D%5Cfrac%7B%5Cpartial%5E2%20v%7D%7B%5Cpartial%20S%5E2%7D" alt="2026-11-13-computational-finance-black-scholes-greeks equation 48" style="display:block; width:auto; max-width:100%; height:auto; margin:0 auto;" />
</div>

Central finite difference를 적용하면

<div style="clear:both; overflow-x:auto; text-align:center; margin:20px 0 24px; line-height:1.2;">
  <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%20%5CGamma%5Capprox%5Cfrac%7Bv%28S%2B%5CDelta%20S%29-2v%28S%29%2Bv%28S-%5CDelta%20S%29%7D%7B%28%5CDelta%20S%29%5E2%7D" alt="2026-11-13-computational-finance-black-scholes-greeks equation 49" style="display:block; width:auto; max-width:100%; height:auto; margin:0 auto;" />
</div>

로 근사할 수 있다.

즉, 주가가 증가하거나 감소했을 때의 가격을 현재 가격과 비교하여 가격 함수의 곡률을 계산한다.

**Vega**는 변동성 σ의 변화에 대한 파생상품 가격의 민감도이다.

<div style="clear:both; overflow-x:auto; text-align:center; margin:20px 0 24px; line-height:1.2;">
  <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%20%5Cmathrm%7BVega%7D%3D%5Cfrac%7B%5Cpartial%20v%7D%7B%5Cpartial%5Csigma%7D" alt="2026-11-13-computational-finance-black-scholes-greeks equation 50" style="display:block; width:auto; max-width:100%; height:auto; margin:0 auto;" />
</div>

Central finite difference를 적용하면

<div style="clear:both; overflow-x:auto; text-align:center; margin:20px 0 24px; line-height:1.2;">
  <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%20%5Cmathrm%7BVega%7D%5Capprox%5Cfrac%7Bv%28%5Csigma%2B%5CDelta%5Csigma%29-v%28%5Csigma-%5CDelta%5Csigma%29%7D%7B2%5CDelta%5Csigma%7D" alt="2026-11-13-computational-finance-black-scholes-greeks equation 51" style="display:block; width:auto; max-width:100%; height:auto; margin:0 auto;" />
</div>

로 근사할 수 있다. 여기서 다른 모형 변수들은 고정하고 변동성만 변화시킨다.

Delta와 Gamma는 주가 S에 대한 미분이고, Vega는 변동성 σ에 대한 미분이라는 차이가 있다.

유한차분은 가격 함수만 계산할 수 있어도 적용할 수 있지만 간격 선택에 따른 오차가 발생한다. 반면 pathwise differentiation과 likelihood-ratio method는 적절한 조건에서 미분을 기댓값 형태로 표현하므로 별도의 유한차분 간격 없이 Monte Carlo 추정량을 구성할 수 있다.

---
