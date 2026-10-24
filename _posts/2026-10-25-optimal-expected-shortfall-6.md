---
title: "Optimal Nonparametric Estimation of the Expected Shortfall Risk 6"
date: 2026-10-24 20:00:00 +0100
categories:
  - Computational Finance
tags:
  - Expected Shortfall
  - Beta Mixing
  - Time Series
  - Robust Statistics
  - Financial Risk
math: true
---
마지막 6편에서는 시간적으로 종속된 데이터(non-i.i.d. data)로 분석 범위를 넓힌다. 금융시장에서는 어제의 변동성과 오늘의 변동성이 관련될 수 있고, 큰 손실이 연속된 시점에 관측되기도 한다. 이런 환경에서 블록별 ES 추정량이 어떤 조건 아래 안정적인지 확인해보려 한다.

원문 Section 6에서는 **geometric β-mixing**, **Theorem 6.1·6.2**, **Lemma 6.3~6.6**을 다룬다. 독립성 가정을 완화하는 방법, 블록을 한 칸씩 건너뛰어 사용하는 이유, coupling을 이용한 확률 집중 증명에 초점을 맞춘다.

---

# 6. Extension to non i.i.d. data (p. 32–38)

## 6.1 종속성을 고려해야 하는 이유

앞선 장의 분석에서는 관측값 X₁, …, X_N이 독립이고 동일한 분포를 따른다고 가정했다. 이러한 가정은 실험과 증명을 구성하는 데 유용하지만, 실제 금융 시계열에는 연속된 관측값 사이의 상관관계와 변동성 군집이 존재할 수 있다.

논문은 strict stationarity(엄격 정상성)와 **geometric β-mixing**을 가정해 종속성을 다룬다. 엄격 정상성은 시간을 함께 이동시켜도 임의의 유한 구간에서 관측값들의 결합분포가 변하지 않는 성질이다. β-mixing은 시간 간격이 커질수록 과거와 미래의 의존성이 얼마나 약해지는지 정량화한다.

두 σ-대수 A, B에 대한 β 계수는 다음과 같이 정의한다.

![Beta mixing definition](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cbeta%28%5Cmathcal%20A%2C%5Cmathcal%20B%29%3D%5Csup_%7BC%5Cin%5Cmathcal%20A%5Cotimes%5Cmathcal%20B%7D%5Cleft%7CP_%7B%5Cmathcal%20A%5Cotimes%5Cmathcal%20B%7D%28C%29-%28P_%7B%5Cmathcal%20A%7D%5Cotimes%20P_%7B%5Cmathcal%20B%7D%29%28C%29%5Cright%7C)

P_{A⊗B}는 두 정보 집합의 실제 결합법칙을 나타내고, P_A⊗P_B는 독립인 경우에 대응하는 곱의 법칙이다. 두 법칙이 가까우면 β 계수도 작아진다.

시차 k에 대한 mixing 계수는 다음과 같다.

![Mixing coefficient](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cbeta%28k%29%3D%5Csup_%7Bt%5Cgeq1%7D%5Cbeta%5C%21%5Cleft%28%5Csigma%28X_1%2C%5Cldots%2CX_t%29%2C%5Csigma%28X_%7Bt%2Bk%7D%2C%5Cldots%2CX_N%29%5Cright%29)

논문의 **geometric β-mixing** 가정은 다음 부등식이다.

![Geometric beta mixing](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cbeta%28k%29%5Cleq%20C_%7B%5Cmathrm%7Bmix%7D%7D%5Crho%5Ek%2C%5Cqquad%20C_%7B%5Cmathrm%7Bmix%7D%7D%3E0%2C%5Cquad0%3C%5Crho%3C1)

ρ가 1보다 작으므로 시차가 커질수록 우변이 기하급수적으로 작아진다. 논문은 적절한 moment와 contraction 조건 아래 일부 GARCH 및 stochastic volatility 모형이 이 구조에 포함된다고 설명한다. 모든 금융 시계열이 자동으로 이 가정을 만족하는 것은 아니다.

> **핵심 해설**
>
> 서로 가까운 시점의 손실은 관련될 수 있다. 충분히 멀리 떨어진 시점들의 정보는 β-mixing 가정 아래에서 점차 독립에 가까워진다. 6장의 추정량은 이 성질을 블록 분할에 활용한다.
{: .prompt-info }

## 6.2 종속 데이터에서의 점근 분산

i.i.d. 데이터에서는 ES의 점근 분산을 같은 시점의 분포함수로 표현했다. 종속 데이터에서는 서로 다른 시점의 정보가 함께 작용하므로 **시차 공분산**을 포함한다.

먼저 시차 k와 임계값 t, s에 대해 다음과 같이 정의한다.

![Lagged indicator covariance](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cgamma_k%28t%2Cs%29%3D%5Coperatorname%7BCov%7D%5C%21%5Cleft%28%5Cmathbf1_%7B%5C%7BX_0%5Cleq%20t%5C%7D%7D%2C%5Cmathbf1_%7B%5C%7BX_k%5Cleq%20s%5C%7D%7D%5Cright%29)

원문은 이 시차 공분산을 합쳐 장기 공분산 함수 Γ를 정의한다.

![Long run covariance kernel](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5CGamma%28t%2Cs%29%3D%5Cgamma_0%28t%2Cs%29%2B2%5Csum_%7Bk%5Cgeq1%7D%5Cgamma_k%28t%2Cs%29)

그리고 q = VaR_{1−α}(X)로 두면 종속 데이터의 ES 관련 분산은 다음과 같다.

![Long run ES variance](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Coverline%5Csigma_%7B%5Cmathrm%7BES%7D%7D%5E2%3D%5Cint_q%5E%5Cinfty%5C%21%5Cint_q%5E%5Cinfty%5CGamma%28t%2Cs%29%5C%2Cds%5C%2Cdt)

i.i.d.에서 시차가 양수인 공분산은 0이므로 이전 장의 분산으로 돌아간다. 종속성이 존재하면 시차 공분산 항이 추가되며 그 부호와 크기에 따라 장기 분산이 달라질 수 있다.

**원문 식 (6.1)**은 다음과 같다.

![Equation 6.1](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Coverline%5Csigma_%7B%5Cmathrm%7BES%7D%7D%5E2%3D%5Csigma_%7B%5Cmathrm%7BES%7D%7D%5E2%2B2%5Csum_%7Bk%5Cgeq1%7D%5Cint_q%5E%5Cinfty%5C%21%5Cint_q%5E%5Cinfty%5Cgamma_k%28t%2Cs%29%5C%2Cds%5C%2Cdt)

Mixing 조건 아래에는 시차 공분산을 제어하는 부등식도 주어진다.

![Mixing covariance bound](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%7C%5Cgamma_k%28t%2Cs%29%7C%5Cleq%5Cmin%5C%7BC_%7B%5Cmathrm%7Bmix%7D%7D%5Crho%5Ek%2C1-F%28t%29%2C1-F%28s%29%5C%7D)

시차 k가 커지면 첫 번째 상한이 감소한다. 나머지 두 항은 꼬리확률과 연결된다. 이 부등식은 장기 공분산의 크기를 다룰 때 활용된다.

> **표기 확인**
>
> 원문 Section 6에는 독립 표본의 분산과 종속 표본의 장기 분산이 별개의 기호로 등장한다. 여기서는 종속 데이터의 분산에 윗줄이 있는 σ̄를 사용했다. 원문의 분산 정의는 적분식 자체를 기준으로 따라가면 된다.
{: .prompt-info }

## 6.3 Theorem 6.1 — 종속 데이터에서도 가능한 유한 표본 보장

저자들은 α ≥ 0.01, ρ < 0.99, C_mix ≥ 4를 추가로 가정한다. 장기 ES 관련 분산이 유한하고, 분위수 함수가 [1−2α, 1−α/2]에서 L-Lipschitz 연속이라고 하자.

ε < c₀√α에 더해, 정리에 인쇄된 표본 수 조건은 다음과 같다.

![Theorem 6.1 sample sizes](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20N%5Cgeq%20c_1%5Clog%28C_%7B%5Cmathrm%7Bmix%7D%7D%29%5Cvarepsilon%5E%7B-2%7D%2C%5Cqquad%20N%5Cgeq%20c_2%5Cvarepsilon%5E%7B-4%7D)

그러면 양의 절대상수 c₀, …, c₄에 대해 특정 추정량 S̃_N이 존재하여 아래 부등식을 만족한다.

![Theorem 6.1 concentration](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20P%5C%21%5Cleft%28%7C%5Cwidetilde%20S_N-%5Cmathrm%7BES%7D_%7B%5Calpha%7D%28X%29%7C%5Cgeq%5Cvarepsilon%5Coverline%5Csigma_%7B%5Cmathrm%7BES%7D%7D%2Bc_3L%5Cvarepsilon%5E2%5Cright%29%5Cleq2%5Cexp%28-c_4N%5Cvarepsilon%5E2%29)

증명에서 얻는 오차는 **ε × 장기 표준편차**와 **Lε²**로 구성된다. 지수형 확률 상한을 얻는 구조도 Theorem 1.4와 유사하다. 새로운 가정과 블록 구성에 맞춰 필요한 표본 크기 조건이 달라진다.

> **원문 표기 주의 — 표본 수 부등호**
>
> Theorem 6.1 본문에는 `N ≥ c₂ε⁻⁴`가 인쇄되어 있다. 바로 다음 문단에서는 같은 조건을 `N ≤ c₂ε⁻⁴`로 설명하고, 이를 상한이라고 부른다. 두 서술이 서로 일치하지 않는다. 아래 정리의 소개에는 정리 문구를 따랐으며, 이 조건을 적용하는 연구나 코드에서는 원문 저자에게 확인할 필요가 있다.
{: .prompt-info }

## 6.4 새로운 추정량의 블록 구성 — 한 칸씩 건너뛰는 이유

2편에서 구성한 추정량은 연속적인 블록을 사용했다. 연속된 블록 사이에도 시간적 의존성이 남을 수 있으므로 Section 6에서는 블록 사이에 **간격(gap)**을 둔다.

블록 크기와 사용할 블록의 수는 다음과 같이 정한다.

![Separated block size](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20m%3D%5Cleft%5Clceil%5Cfrac%7BC_0%5Clog%28C_%7B%5Cmathrm%7Bmix%7D%7D%29%7D%7B%5Cvarepsilon%5E2%7D%5Cright%5Crceil%2C%5Cqquad%20n%3D%5Cleft%5Clfloor%5Cfrac%7BN%7D%7B2m%7D%5Cright%5Crfloor)

전체 데이터를 길이 m인 2n개 블록으로 나눈 뒤, 추정량 계산에는 짝수 번호의 블록 I′₂, I′₄, …, I′₂ₙ만 사용한다. I_j = I′_{2j}로 쓴다.

예를 들어 m = 3일 때 배열은 다음과 같다.

```text
원본 관측값: [ X1 X2 X3 ] [ X4 X5 X6 ] [ X7 X8 X9 ] [ X10 X11 X12 ]
블록 구분:      I'1 (건너뜀)   I'2 (사용)     I'3 (건너뜀)    I'4 (사용)
```

같은 그룹 내부에서는 의존성이 남아 있어도, 사용한 블록 사이에는 최소 m개의 관측값 간격이 생긴다. 블록 ES 추정량 사이의 β 계수는 다음과 같이 제한된다.

![Block statistic mixing](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cbeta%5E%7BT%7D%28k%29%5Cleq%5Cbeta%28km%29%5Cleq%20C_%7B%5Cmathrm%7Bmix%7D%7D%5Crho%5E%7Bmk%7D)

m을 키우면 서로 다른 사용 블록 사이의 의존성이 더 약해지는 효과를 기대할 수 있다. 블록 개수 n은 감소하므로 두 요구 사이에서 균형을 잡아야 한다. 전체 표본의 plug-in T̂_N은 여전히 전체 데이터를 사용한다. 건너뛴 블록은 분위수 구간을 구성할 때 제외된다.

이 블록들에서 얻은 ES 추정값을 사용해 블록 분위수 구간을 만든 후, 전체 plug-in 추정값을 그 범위 안으로 제한하여 S̃_N을 정의한다.

---

## 6.5 Lemma 6.3과 6.4 — Coupling으로 의존성을 다루는 증명

Section 2에서 높은 확률의 집중 결과를 얻을 수 있었던 핵심은 블록 추정값들 사이의 독립성이었다. 이제 블록 통계량에 남은 의존성을 **coupling**으로 처리한다.

### Lemma 6.4 — 독립 복제열 구성

β-mixing인 수열 Z₁, Z₂, …가 주어졌다고 하자. 원문은 다음 세 조건을 만족하는 새로운 수열 Y₁, Y₂, …를 구성할 수 있음을 보인다.

1. 각 시점에서 Y_t와 Z_t의 주변분포가 같다.
2. Y₁, Y₂, …는 서로 독립이다.
3. 두 변수가 달라질 확률을 β(1)로 제어할 수 있다.

![Lemma 6.4 coupling mismatch](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20P%28Y_t%5Cne%20Z_t%29%5Cleq%5Cbeta%281%29)

> **Proof — Lemma 6.4**
>
> **Step 1.** 첫 시점에는 Y₁ = Z₁로 정한다.
>
> **Step 2.** 이미 구성한 Y₁, …, Y_{T−1}과 새로운 Z_T의 의존성을 β(1)로 제어한다. 보조 난수 U_T를 도입하고 Berbee의 coupling lemma를 적용한다.
>
> **Step 3.** Y_T는 Z_T와 같은 주변분포를 가지면서 앞에서 구성한 Y들로부터 독립이 되도록 선택한다. 이때 두 값이 일치하지 않을 확률이 β(1) 이하라는 조건도 얻는다.
>
> **Step 4.** 수학적 귀납법으로 모든 T에 대해 조건을 구성한다.
{: .prompt-tip }

이 Lemma 덕분에 **약하게 종속된 블록 통계량의 확률 집중 문제**를 독립 수열에 대한 문제와 연결할 수 있다.

### Lemma 6.3 — 나쁜 블록이 많이 나올 확률

각 블록의 plug-in 추정량 T̂_{I_j}가 실제 ES에서 δ 이상 벗어날 확률이 1/20 이하라고 가정한다.

![Equation 6.2](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20P%5C%21%5Cleft%28%7C%5Cwidehat%20T_%7BI_j%7D-%5Cmathrm%7BES%7D_%7B%5Calpha%7D%28X%29%7C%5Cgeq%5Cdelta%5Cright%29%5Cleq%5Cfrac1%7B20%7D)

δ는 다음과 같이 정한다.

![Delta radius](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cdelta%3D%5Cvarepsilon%5Coverline%5Csigma_%7B%5Cmathrm%7BES%7D%7D%2B10L%5Cvarepsilon%5E2)

그러면 블록 중 최소 30%가 이 범위를 벗어날 확률은 지수적으로 작아진다.

![Lemma 6.3 concentration](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20P%5C%21%5Cleft%28%5C%23%5C%7Bj%3A%7C%5Cwidehat%20T_%7BI_j%7D-%5Cmathrm%7BES%7D_%7B%5Calpha%7D%28X%29%7C%5Cgeq%5Cdelta%5C%7D%5Cgeq0.3n%5Cright%29%5Cleq2%5Cexp%28-n%2F50%29)

> **Proof — Lemma 6.3**
>
> **Step 1.** Lemma 6.4를 적용해 블록 ES 값 T̂_{I_j}를 독립적인 Y_j와 결합한다. 각 블록이 복제값과 달라질 확률은 최대 C_mixρ^m이다.
>
> **Step 2.** 결합된 값 중 하나라도 다른 사건의 확률을 union bound로 묶는다.
{: .prompt-tip }

![Coupling union bound](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20P%28%5Cexists%20j%3AY_j%5Cne%5Cwidehat%20T_%7BI_j%7D%29%5Cleq%20nC_%7B%5Cmathrm%7Bmix%7D%7D%5Crho%5Em)

> **Proof — Lemma 6.3 (계속)**
>
> **Step 3.** 복제 수열의 나쁜 블록 확률을 1/10 이하로 제한한 뒤, Section 2의 독립 블록 집중 부등식(Lemma 2.1)을 적용한다.
>
> **Step 4.** 결합 실패 사건과 독립 복제 수열의 집중 실패 사건을 합쳐 원하는 확률 상한을 얻는다.
{: .prompt-tip }

여기서 상한의 앞에 2가 붙는 것은 두 종류의 실패 사건을 함께 제어하는 구조와 연결된다.

---

## 6.6 Lemma 6.5 — 선형 오차항의 제어

2편의 증명과 동일하게 블록별 ES 오차를 **선형항 L_j**와 **나머지항 E_j**로 분해한다.

![Error decomposition](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cwidehat%20T_%7BI_j%7D-%5Cmathrm%7BES%7D_%7B%5Calpha%7D%28X%29%3DL_j%2BE_j)

Lemma 6.5는 선형항의 오차확률을 다음처럼 제한한다.

![Lemma 6.5](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20P%28%7CL_j%7C%5Cgeq%5Cvarepsilon%5Coverline%5Csigma_%7B%5Cmathrm%7BES%7D%7D%29%5Cleq%5Cfrac1%7B11%7D)

> **Proof — Lemma 6.5**
>
> 블록 경험분포함수의 공분산을 관측값의 시차 공분산으로 전개한다. 독립 표본에서는 서로 다른 시점의 공분산 항이 0이지만, 이 절에서는 γ_k가 합에 들어간다.
{: .prompt-tip }

![Block empirical covariance](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Coperatorname%7BCov%7D%28%5Cwidehat%20F_%7BI_j%7D%28t%29%2C%5Cwidehat%20F_%7BI_j%7D%28s%29%29%3D%5Cfrac1%7Bm%5E2%7D%5Csum_%7Bi%3D1%7D%5E%7Bm%7D%5Csum_%7B%5Cell%3D1%7D%5E%7Bm%7D%5Cgamma_%7B%7Ci-%5Cell%7C%7D%28t%2Cs%29)

이 합은 긴 시간 간격에 대한 공분산 가중합으로 표현할 수 있다.

![Equation 6.5](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Coperatorname%7BCov%7D%28%5Cwidehat%20F_%7BI_j%7D%28t%29%2C%5Cwidehat%20F_%7BI_j%7D%28s%29%29%5Cleq%5Cfrac1m%5CGamma%28t%2Cs%29)

원문은 이 상한을 이용해 선형항의 이차모멘트를 제어하고, Chebyshev 부등식을 적용한다.

![Linear term variance](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20E%5BL_j%5E2%5D%5Cleq%5Cfrac%7B%5Coverline%5Csigma_%7B%5Cmathrm%7BES%7D%7D%5E2%7D%7Bm%7D)

> **Proof — Lemma 6.5 (마무리)**
>
> m이 충분히 크도록 선택했기 때문에 Chebyshev 
부등식으로 ![Probability bound](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20P%5Cleft%28%5Cleft%7CL_j%5Cright%7C%5Cgeq%5Cvarepsilon%5Cbar%7B%5Csigma%7D_%7B%5Cmathrm%7BES%7D%7D%5Cright%29%5Cleq%5Cfrac%7B1%7D%7B11%7D)를 얻는다. 독립성 가정을 완화하면서도 블록 수준의 확률 제어를 유지하는 단계다.
{: .prompt-tip }

## 6.7 Lemma 6.6 — 분위수 근방의 나머지항

나머지항 E_j의 크기는 경험적 CDF가 **1−α 주변에서 실제 CDF를 얼마나 정확하게 따라가는지**에 달려 있다.

논문은 γ를 다음과 같이 선택한다.

![Gamma threshold](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cgamma%3D%5Cfrac54%5Cvarepsilon%5Csqrt%5Calpha)

Lemma 6.6은 다음 오차확률을 제공한다.

![Lemma 6.6](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20P%5C%21%5Cleft%28%7CE_j%7C%3E%5Cfrac%7B3%5Cgamma%7D%7B%5Calpha%7D%5Cleft%5BF%5E%7B-1%7D%281-%5Calpha%2B%5Cgamma%29-F%5E%7B-1%7D%281-%5Calpha-%5Cgamma%29%5Cright%5D%5Cright%29%5Cleq%5Cfrac1%7B110%7D)

> **Proof — Lemma 6.6**
>
> **Step 1.** 두 분위수 경계 t₀, t₁을 정한다. 각각 F(t₀) = 1−α−γ, F(t₁) = 1−α+γ이다.
>
> **Step 2.** 경험적 CDF가 이 두 지점에서 γ 이상 벗어날 확률을 제어한다. 시차 공분산은 mixing 계수와 꼬리확률 중 작은 값으로 제한한다.
{: .prompt-tip }

![Quantile covariance bound](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cgamma_k%28t%2Ct%29%5Cleq%5Cmin%5C%7BC_%7B%5Cmathrm%7Bmix%7D%7D%5Crho%5Ek%2C1-F%28t%29%5C%7D)

> **Proof — Lemma 6.6 (계속)**
>
> **Step 3.** k에 대해 공분산 상한을 합산하고, 블록 크기 m을 이용해 경험분포의 분산을 작게 만든다.
{: .prompt-tip }

![CDF variance bound](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Coperatorname%7BVar%7D%28%5Cwidehat%20F_%7BI_j%7D%28t%29%29%5Cleq%5Cfrac%7B%5Calpha%20C_1%5Clog%28C_%7B%5Cmathrm%7Bmix%7D%7D%29%7D%7Bm%7D)

Chebyshev 부등식과 γ의 선택에 따라 원문 식 (6.6)을 얻는다.

![Equation 6.6](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20P%5C%21%5Cleft%28%7C%5Cwidehat%20F_%7BI_j%7D%28t%29-F%28t%29%7C%5Cgeq%5Cgamma%5Cright%29%5Cleq%5Cfrac1%7B220%7D%2C%5Cqquad%20t%5Cin%5C%7Bt_0%2Ct_1%5C%7D)

> **Proof — Lemma 6.6 (마무리)**
>
> 두 분위수 경계에서 경험적 CDF의 오차가 작다는 사실을 얻었다. 나머지 증명은 Section 2의 Lemma 2.3에서 사용했던 분위수 구간과 나머지항 제어 논리를 적용한다.
{: .prompt-tip }

원문 Lemma 6.6의 계산에는 C_mix와 ρ의 크기에 관한 조건이 중요하다. 특히 본문의 보조 계산 중에는 ρ ≤ 0.9라는 더 강한 부등식이 쓰이지만, Theorem 6.1의 서두에는 ρ < 0.99로 쓰여 있다. 이 차이는 실제 상수를 적용할 때 별도로 검토할 필요가 있다.

---

## 6.8 Theorem 6.1과 6.2의 증명 완성

Lemma 6.5와 Lemma 6.6을 결합하면 개별 블록의 오차가 δ를 넘을 확률을 제어할 수 있다. Lemma 6.3에서는 **많은 블록이 동시에 실패할 확률**을 지수적으로 작게 제한했다.

> **Proof — Theorem 6.1**
>
> Section 2의 Theorem 1.4 증명 구조를 따른다. 독립 블록 집중에 사용했던 Lemma 2.1 자리에 Lemma 6.3을, 선형항 제어 Lemma 2.2 자리에 Lemma 6.5를, 나머지항 제어 Lemma 2.3 자리에 Lemma 6.6을 적용한다.
>
> 블록 분위수로 만든 구간이 높은 확률로 실제 ES 근처에 놓이고, 최종 추정값을 그 구간으로 제한하면서 Theorem 6.1의 지수형 오차 보장을 얻는다.
{: .prompt-tip }

### Theorem 6.2 — 종속 데이터와 적대적 변조

Theorem 6.1의 가정 아래 전체 N개 중 최대 C₀Nε²개가 악의적으로 변조되어도 다음 관계가 성립한다.

![Theorem 6.2](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20P%5C%21%5Cleft%28%7C%5Cwidetilde%20S_N-%5Cmathrm%7BES%7D_%7B%5Calpha%7D%28X%29%7C%5Cgeq%5Cvarepsilon%5Coverline%5Csigma_%7B%5Cmathrm%7BES%7D%7D%2BC_1L%5Cvarepsilon%5E2%5Cright%29%5Cleq2%5Cexp%28-C_2N%5Cvarepsilon%5E2%29)

> **Proof — Theorem 6.2**
>
> 원문은 Theorem 1.8에서 사용한 적대적 변조 논리를 적용한다. 변조된 관측값 수가 제한되어 있으므로 지나치게 많은 블록의 추정값이 동시에 망가지는 사건을 억제할 수 있다. 블록 분위수로 구성한 절단 구간은 결합 오차와 데이터 변조를 함께 고려한 확률 경계 안에서 안정성을 유지한다.
{: .prompt-tip }

> **6편에서 기억할 점**
>
> 블록을 한 칸씩 건너뛰어 사용하는 방식은 시간 간격을 확보한다. β-mixing과 coupling은 약하게 종속된 블록 통계량을 독립 블록의 집중 논리와 연결한다. Theorem 6.1은 종속된 데이터에 대한 오차 보장을, Theorem 6.2는 제한된 데이터 변조까지 고려한 보장을 제시한다.
{: .prompt-info }

---

## 6.9 금융 시계열 연구에 적용한다면

금융시장에서 손실의 극단적 크기와 발생 빈도를 함께 연구할 때, 관측값들이 시간에 따라 연결된다는 점은 중요하다. 특히 시장 변동성이 커지는 시기에 연속된 큰 손실이 관측될 수 있다. 이런 환경에서 위험 추정량을 평가하려면 **분포의 두꺼운 꼬리**, **시간적 종속성**, **데이터 품질**을 함께 살펴볼 필요가 있다.

이 논문을 Python 프로젝트로 확장한다면 독립 Pareto 표본과 종속적인 시계열 모형을 나란히 생성하고, plug-in과 블록 추정량의 오차확률을 비교할 수 있다. 다만 이론을 실험으로 검증할 때에는 정상성·mixing 가정·장기 분산의 유한성을 별도로 확인해야 한다.

Section 6의 논의는 5장에서 다뤘던 수치실험의 다음 단계와 연결된다. 현재 시계열 데이터의 종속성을 고려하면서 ES를 추정한다면 어떤 통계적 보장을 얻을 수 있는지에 관한 연구 질문이 남는다.

---

이 논문에서 다룬 주제는 **Expected Shortfall을 제한된 표본으로 얼마나 안정적으로 추정할 수 있는가**이다. 처음에는 경험적 분포함수를 사용하는 plug-in 추정량에서 출발했고, 강건한 블록 기반 추정량의 유한 표본 보장과 최적성을 살펴봤다. Pareto 분포 계산과 Numerical Experiments를 통해 극단적 손실의 영향도 확인했다. 마지막 Section 6에서는 시간적으로 종속된 데이터까지 방법을 확장했다.

실제 금융 데이터를 다룰 때는 오차의 크기뿐 아니라 **그 오차가 발생할 확률과 분석의 가정**까지 살펴봐야 한다. 이번 6편으로 논문 본문의 주요 흐름을 마무리한다.

---

> Bartl, D., & Eckstein, S. (2026). *Optimal Nonparametric Estimation of the Expected Shortfall Risk*. arXiv:2405.00357v2, Section 6, pp. 32–38. https://arxiv.org/abs/2405.00357
