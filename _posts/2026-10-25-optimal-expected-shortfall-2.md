---
title: "Optimal Nonparametric Estimation of the Expected Shortfall Risk — 2"
date: 2026-10-25 12:00:00 +0100
categories:
  - Computational Finance
tags:
  - Expected Shortfall
  - Robust Statistics
  - Proof
  - Concentration Inequalities
  - Adversarial Robustness
math: true
---

# 2. Proofs of Theorem 1.4 and Theorem 1.8 (원문 p. 12–18)

지난 1편에서는 Plug-in Estimator가 유한 표본과 heavy-tailed 분포에서 겪는 문제를 확인하고, 저자들이 제안한 강건 ES 추정량을 살펴보았다. 전체 표본으로 계산한 ES 추정값을 블록별 추정값의 분위수 구간에 제한하는 방법이었다.

이번 2편에서는 **그 방법이 어떤 이유로 높은 확률의 오차 보장을 갖는지** 증명 과정을 따라가려 한다. 특히 Theorem 1.4에서 주장한 finite-sample bound와 Theorem 1.8의 adversarial robustness를 중심으로 읽는다.

본문에는 Hoeffding's inequality, Chebyshev's inequality, Bernstein's inequality와 같은 확률 부등식이 등장한다. 각 부등식을 어디에 적용하는지, 앞에서 얻은 결과가 다음 증명에 어떻게 쓰이는지 순서대로 정리할 것이다. 


1. **Lemma 2.1**: 각 블록 추정이 대체로 정확하다면, 대부분의 블록이 정확한 범위 안에 들어올 확률이 높다는 사실을 보인다.
2. **ES의 적분 표현 (2.4)**: 분포함수의 차이를 이용해 추정 오차를 분석할 수 있도록 ES 공식을 바꾼다.
3. **Lemma 2.2**: 추정 오차의 선형항을 점근 분산으로 제어한다.
4. **Lemma 2.3**: 경험적 CDF가 분위수 경계 주변에서 얼마나 흔들리는지 분석해 나머지항을 제어한다.
5. **Theorem 1.4**: 블록별 오차와 Lemma 2.1을 합쳐 전체 추정량의 오차확률을 구한다.
6. **Theorem 1.8과 Corollary 2.4**: 일부 데이터가 변조된 상황에서도 강건성이 유지되는 조건과 확률수렴을 확인한다.

---

## 2.1. 증명을 위한 표기와 기본 가정 (p. 12)

증명에서는 표기를 단순하게 하기 위해 ES 관련 점근 분산을 σ²_ES, 손실분포의 CDF를 F로 적는다. 데이터 N개를 크기 m의 서로 겹치지 않는 n개 블록으로 분할한다.

<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20m%3D%5Cleft%5Clceil%5Cfrac%7B11%7D%7B%5Cvarepsilon%5E2%7D%5Cright%5Crceil%2C%5Cquad%20n%3D%5Cleft%5Clfloor%5Cfrac%7BN%7D%7Bm%7D%5Cright%5Crfloor" alt="block" loading="lazy">

논문은 ε ≤ √α/6, α < 1/2, Nε² ≥ 140을 가정한다. 또한 ε ≤ 1이다. 이 조건에서 다음을 얻는다.

**(2.1)**

<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cfrac%7B12%7D%7B%5Cvarepsilon%5E2%7D%5Cgeq%20m%5Cgeq11%2C%5Cqquad%20n%5Cgeq%5Cfrac%7BN%5Cvarepsilon%5E2%7D%7B14%7D%5Cgeq10" alt="21" loading="lazy">

m ≥ 11이라는 사실은 뒤에서 Chebyshev 부등식을 적용할 때 사용한다. n ≥ 10이라는 조건은 경험적 분위수를 구하기 위한 블록 수를 확보하고, 블록 중 일부가 부정확하더라도 나머지 블록을 이용할 수 있게 한다.

### 왜 F가 연속이고 strictly increasing이라고 가정할까?

저자들은 이 절의 증명을 진행하는 동안 F가 연속이며 엄격하게 증가한다고 가정한다. 뒤에서 분위수 F⁻¹(1−α±γ)를 다루기 편하게 만드는 조건이다.

일반적인 손실분포로 돌아가기 위해 다음과 같이 작은 정규분포 잡음을 더한다.

<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20Y%5E%7B%5Ceta%7D%3DX%2B%5Ceta%20G%2C%5Cquad%20G%5Csim%5Cmathcal%7BN%7D%280%2C1%29%2C%5Cquad%5Ceta%3E0" alt="gauss" loading="lazy">

η를 0으로 보내면 다음과 같이 ES와 관련 분산이 원래 분포의 값에 수렴한다.

<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cmathrm%7BES%7D_%7B%5Calpha%7D%28Y%5E%7B%5Ceta%7D%29%5Clongrightarrow%5Cmathrm%7BES%7D_%7B%5Calpha%7D%28X%29%2C%5Cquad%5Csigma%5E2_%7B%5Cmathrm%7BES%7D_%7B%5Calpha%7D%28Y%5E%7B%5Ceta%7D%29%7D%5Clongrightarrow%5Csigma%5E2_%7B%5Cmathrm%7BES%7D_%7B%5Calpha%7D%28X%29%7D" alt="gausslim" loading="lazy">

η를 충분히 작게 정하면 표본과 새 표본의 차이도 임의로 작게 만들 수 있다. 이렇게 연속성을 확보한 상황에서 얻은 결과를 원래 분포에 적용할 수 있다.

---

## 2.2. Lemma 2.1 — Binomial Concentration (p. 12–13)

첫 번째 Lemma는 **좋지 않은 블록이 너무 많이 발생할 확률**을 제어한다. 블록 하나의 추정값이 정답 주변의 허용구간을 벗어날 확률이 충분히 작다는 가정에서 시작한다.

독립인 Y₁, …, Yₙ이 같은 참값 y₀에 대해 다음을 만족한다고 하자.

<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20P%28%5Clvert%20Y_j-y_0%5Crvert%3E%5Cdelta%29%5Cleq%5Cfrac%7B1%7D%7B10%7D" alt="yass" loading="lazy">

그러면 **Lemma 2.1**은 다음을 보장한다.

<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20P%5C%21%5Cleft%28%5C%23%5C%7Bj%3A%5Clvert%20Y_j-y_0%5Crvert%3E%5Cdelta%5C%7D%5Cgeq0.3n%5Cright%29%5Cleq%5Cexp%28-n%2F50%29" alt="lemma21" loading="lazy">

각 블록이 잘못될 확률이 최대 10%라면, **전체 블록의 30% 이상이 잘못될 확률은 n에 대해 지수적으로 작아진다.**

<div markdown="1" style="background: rgba(34, 197, 94, 0.10); border-left: 4px solid #22c55e; border-radius: 8px; padding: 18px 20px; margin: 18px 0;">

**Proof**

j번째 블록이 허용구간을 벗어났을 때 1을 갖는 지시변수를 정의한다.

<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20B_j%3D%5Cmathbf%7B1%7D_%7B%5C%7B%5Clvert%20Y_j-y_0%5Crvert%3E%5Cdelta%5C%7D%7D%2C%5Cquad%20E%5BB_j%5D%5Cleq0.1" alt="bad" loading="lazy">

B_j는 0 또는 1을 갖는다. 전체 블록 중 문제가 발생한 개수는 B₁부터 Bₙ까지의 합이다. 이 변수들은 독립이므로 Hoeffding's inequality를 적용할 수 있다.

<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20P%5C%21%5Cleft%28%5Csum_%7Bj%3D1%7D%5En%20B_j%5Cgeq%5Csum_%7Bj%3D1%7D%5En%20E%5BB_j%5D%2B%5Clambda%5Cright%29%5Cleq%5Cexp%5C%21%5Cleft%28-%5Cfrac%7B%5Clambda%5E2%7D%7B2n%7D%5Cright%29" alt="hoeff" loading="lazy">

기댓값의 합은 최대 0.1n이다. 우리가 제어하려는 기준은 0.3n이므로, 기댓값에서 0.2n만큼 더 벗어나는 사건을 계산하면 된다.

<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Clambda%3D0.2n%5Cquad%5CLongrightarrow%5Cquad%5Cexp%5C%21%5Cleft%28-%5Cfrac%7B%280.2n%29%5E2%7D%7B2n%7D%5Cright%29%3De%5E%7B-n%2F50%7D" alt="hoeff2" loading="lazy">

이것으로 Lemma 2.1의 결론을 얻는다.

</div>

**이 Lemma가 중요한 이유**는 이후에 블록별 추정 오차를 먼저 계산한 다음, 그 결과를 전체 블록의 경험적 분위수에 적용할 수 있기 때문이다. 증명의 출발점에서 블록 사이의 독립성이 직접 사용된다.

---

## 2.3. 전체 추정량의 문제를 블록 하나의 문제로 바꾸기 (p. 13)

각 블록 I_j에서 계산한 Plug-in Estimator를 T̂_Ij로 적자. 블록들이 겹치지 않고 원래 관측치들이 i.i.d.이므로 블록 추정값도 서로 독립이다.

저자들은 다음 한 가지를 증명하는 데 집중한다.

**(2.2)**

<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20P%5C%21%5Cleft%28%5Clvert%5Cwidehat%7BT%7D_%7BI_j%7D-%5Cmathrm%7BES%7D_%7B%5Calpha%7D%28X%29%5Crvert%3E%5Cvarepsilon%5Csigma_%7B%5Cmathrm%7BES%7D%7D%2B10L%5Cvarepsilon%5E2%5Cright%29%5Cleq%5Cfrac%7B1%7D%7B10%7D" alt="22" loading="lazy">

σ_ES는 ES 관련 점근 표준편차이고, L은 분위수 함수의 Lipschitz 상수다. εσ_ES는 분산과 관련된 오차 규모를 나타내며, 10Lε²은 분위수 경계 주변의 영향을 반영한다.

여기에서 허용 오차 d를 다음처럼 쓰자.

<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20d%3D%5Cvarepsilon%5Csigma_%7B%5Cmathrm%7BES%7D%7D%2B10L%5Cvarepsilon%5E2%2C%5Cquad0.3n%5Cleq%20k%5Cleq%5Cell%5Cleq0.7n" alt="def_d" loading="lazy">

Lemma 2.1을 적용하면 **적어도 70%의 블록 추정값이 참값 ES에서 거리 d 이내**에 들어올 확률이 매우 높다. 따라서 블록 추정값을 정렬했을 때 중앙부의 순서통계량도 그 구간 안에 위치한다.

**(2.3)**

<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cwidehat%7BT%7D_%7B%28k%29%7D%2C%5Cwidehat%7BT%7D_%7B%28%5Cell%29%7D%5Cin%5B%5Cmathrm%7BES%7D_%7B%5Calpha%7D%28X%29-d%2C%5C%20%5Cmathrm%7BES%7D_%7B%5Calpha%7D%28X%29%2Bd%5D" alt="23" loading="lazy">

1편에서 정의한 최종 추정량을 다시 쓰면

<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cwidehat%7BS%7D_N%3D%5Cmin%5C%7B%5Cmax%5C%7B%5Cwidehat%7BT%7D_N%2C%5Cwidehat%7BQ%7D%28%5Cbeta_1%29%5C%7D%2C%5Cwidehat%7BQ%7D%28%5Cbeta_2%29%5C%7D" alt="clip" loading="lazy">

이다. β₁, β₂가 0.35에서 0.65 사이에 있으므로, 두 분위수는 블록 추정값을 정렬했을 때 중앙부에 해당한다. 중앙부 추정값들이 허용구간 안에 있으면, 그 사이로 제한된 Ŝ_N 역시 같은 구간 안에 놓인다.

따라서 앞으로의 핵심 과제는 **각 블록에서 (2.2)가 왜 성립하는지**를 보이는 것이다.

---

## 2.4. Expected Shortfall의 새로운 적분 표현 (2.4) (p. 13–14)

ES의 원래 정의는 분위수 함수의 적분이다. 이번 증명에서는 **경험적 누적분포함수 F̂와 실제 누적분포함수 F 사이의 차이**를 직접 분석해야 하므로, 분포함수를 적분하는 표현을 사용한다.

먼저 다음 함수를 정의한다.

<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cpsi%28x%29%3D%5Cmin%5Cleft%5C%7B%5Cfrac%7B1-x%7D%7B%5Calpha%7D%2C1%5Cright%5C%7D%2C%5Cquad%20x%5Cin%5B0%2C1%5D" alt="psi" loading="lazy">

이 함수는 F(t)가 1−α 아래에 있을 때 1이고, 그보다 커지면 선형으로 감소한다. 손실분포의 상위 α 구간을 구분해 주는 함수다.

이때 임의의 손실확률변수 Y에 대해 **(2.4)**가 성립한다.

<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cmathrm%7BES%7D_%7B%5Calpha%7D%28Y%29%3D%5Cint_%7B-%5Cinfty%7D%5E%7B0%7D%5B%5Cpsi%28F_Y%28t%29%29-1%5Ddt%2B%5Cint_0%5E%5Cinfty%5Cpsi%28F_Y%28t%29%29dt" alt="24" loading="lazy">

<div markdown="1" style="background: rgba(34, 197, 94, 0.10); border-left: 4px solid #22c55e; border-radius: 8px; padding: 18px 20px; margin: 18px 0;">

**Proof**

먼저 P(Y > 0) ≥ α인 경우를 생각한다. 이 조건에서는 ES 계산에 포함되는 분위수가 음수가 되지 않는다. 원래 ES 적분을 지시함수의 적분으로 표현하면

<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cmathrm%7BES%7D_%7B%5Calpha%7D%28Y%29%3D%5Cfrac1%5Calpha%5Cint_%7B1-%5Calpha%7D%5E1%5Cint_0%5E%5Cinfty%5Cmathbf%7B1%7D_%7B%5C%7Bt%5Cleq%20F_Y%5E%7B-1%7D%28u%29%5C%7D%7D%5C%2Cdt%5C%2Cdu" alt="positive" loading="lazy">

가 된다. Fubini's theorem으로 적분 순서를 바꾸면, 특정 t보다 큰 분위수 값이 차지하는 비율을 계산할 수 있다. 그 결과는 다음과 같다.

<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%3D%5Cint_0%5E%5Cinfty%5Cfrac%7B%5Cmin%5C%7B1-F_Y%28t%29%2C%5Calpha%5C%7D%7D%7B%5Calpha%7D%5C%2Cdt%3D%5Cint_0%5E%5Cinfty%5Cpsi%28F_Y%28t%29%29%5C%2Cdt" alt="positive2" loading="lazy">

이 경우 t가 음수인 구간에서는 ψ(F_Y(t))−1이 0이어서 (2.4)가 성립한다.

일반적인 경우에는 충분히 큰 상수 c를 더해 P(Y+c > 0) ≥ α가 되도록 한다. ES의 평행이동 성질과 CDF의 이동 관계를 이용한다.

<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cmathrm%7BES%7D_%7B%5Calpha%7D%28Y%2Bc%29%3D%5Cmathrm%7BES%7D_%7B%5Calpha%7D%28Y%29%2Bc%2C%5Cquad%20F_%7BY%2Bc%7D%28t%29%3DF_Y%28t-c%29" alt="shift" loading="lazy">

따라서 (2.4)는 일반적인 Y에도 적용된다.

</div>

### 이 표현을 왜 사용하는가?

ψ는 거의 모든 지점에서 미분 가능하며, 도함수는 다음과 같다.

<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cpsi%5E%7B%5Cprime%7D%28x%29%3D-%5Cfrac1%5Calpha%5Cmathbf%7B1%7D_%7B%281-%5Calpha%2C1%29%7D%28x%29%5Cquad%5Ctext%7Ba.e.%7D" alt="deriv" loading="lazy">

ES 추정값의 차이를 ψ에 대한 변화량으로 쓸 수 있으므로, **1차 선형항과 그 밖의 오차를 분리**해 분석할 수 있다. 다음 Lemma 2.2와 2.3이 각각 한 항씩 담당한다.

---

## 2.5. 추정 오차를 두 항으로 분해하기 (p. 14)

블록 I_j의 경험적 CDF를 F̂_Ij라고 하자. 원래 ES와 블록 ES의 차이는 (2.4)를 적용해 적분 형태로 나타낼 수 있다.

저자들은 이 추정 오차를 다음 두 항으로 분해한다.

<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cwidehat%7BT%7D_%7BI_j%7D-%5Cmathrm%7BES%7D_%7B%5Calpha%7D%28X%29%3D%5Cmathcal%7BL%7D_j%2B%5Cmathcal%7BE%7D_j" alt="split" loading="lazy">

**Linearized term: L_j**

<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cmathcal%7BL%7D_j%3D%5Cint_%7B%5Cmathbb%7BR%7D%7D%5Cpsi%5E%7B%5Cprime%7D%28F%28t%29%29%5C%2C%5B%5Cwidehat%20F_%7BI_j%7D%28t%29-F%28t%29%5D%5C%2Cdt" alt="lin" loading="lazy">

**Remainder term: E_j**

<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cmathcal%7BE%7D_j%3D%5Cint_%7B%5Cmathbb%7BR%7D%7D%5B%5Cpsi%28%5Cwidehat%20F_%7BI_j%7D%28t%29%29-%5Cpsi%28F%28t%29%29-%5Cpsi%5E%7B%5Cprime%7D%28F%28t%29%29%28%5Cwidehat%20F_%7BI_j%7D%28t%29-F%28t%29%29%5Ddt" alt="err" loading="lazy">

L_j는 CDF 오차 F̂_Ij−F에 선형적으로 반응하는 항이다. E_j에는 ψ의 기울기가 바뀌는 경계 근처의 효과가 포함된다.

ψ의 기울기는 F(t)가 1−α를 지나는 지점에서 달라진다. 그래서 나머지항은 **VaR 경계 주변에서 경험적 CDF가 어떻게 움직이는지** 살펴보아야 한다.

이제 두 항의 크기를 각각 확률적으로 제어한다.

---

## 2.6. Lemma 2.2 — 선형항의 분산과 Chebyshev's inequality (p. 14–15)

**Lemma 2.2**는 다음 결과를 보인다.

<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20P%28%5Clvert%5Cmathcal%7BL%7D_j%5Crvert%3E%5Cvarepsilon%5Csigma_%7B%5Cmathrm%7BES%7D%7D%29%5Cleq%5Cfrac1%7B11%7D" alt="lemma22" loading="lazy">

<div markdown="1" style="background: rgba(34, 197, 94, 0.10); border-left: 4px solid #22c55e; border-radius: 8px; padding: 18px 20px; margin: 18px 0;">

**Proof**

### Step 1. 선형항의 2차 모멘트 계산

경험적 CDF는 지시변수의 평균이다. 서로 독립인 m개의 관측값을 이용하므로 두 지점 t와 s에서의 공분산은 다음과 같다.

<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Coperatorname%7BCov%7D%28%5Cwidehat%20F_%7BI_j%7D%28t%29%2C%5Cwidehat%20F_%7BI_j%7D%28s%29%29%3D%5Cfrac%7BF%28%5Cmin%5C%7Bt%2Cs%5C%7D%29-F%28t%29F%28s%29%7D%7Bm%7D" alt="cov" loading="lazy">

왜 분모에 m이 등장할까? 각 관측값에서 생기는 공분산을 m개 더한 뒤, 평균을 계산하면서 m²으로 나누기 때문이다. 독립성 덕분에 서로 다른 관측치 사이의 교차 공분산은 0이 된다.

선형항 L_j의 제곱에 기댓값을 취한 뒤 적분 순서를 바꾸면, 1편에서 정의한 ES 점근 분산과 같은 형태가 나온다.

<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20E%5B%5Cmathcal%7BL%7D_j%5E2%5D%3D%5Cfrac%7B%5Csigma_%7B%5Cmathrm%7BES%7D%7D%5E2%7D%7Bm%7D" alt="var" loading="lazy">

**블록 크기 m이 커질수록 선형항의 분산이 1/m 비율로 감소**하는 것이다.

### Step 2. Chebyshev's inequality 적용

평균이 0인 L_j에 Chebyshev 부등식을 적용한다.

<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20P%28%5Clvert%5Cmathcal%7BL%7D_j%5Crvert%3E%5Cvarepsilon%5Csigma_%7B%5Cmathrm%7BES%7D%7D%29%5Cleq%5Cfrac%7B1%7D%7Bm%5Cvarepsilon%5E2%7D%5Cleq%5Cfrac1%7B11%7D" alt="cheb" loading="lazy">

마지막 부등식은 (2.1)의 m ≥ 11/ε²에서 얻는다.

이 결과로 **선형항이 εσ_ES를 초과할 확률은 최대 1/11**임을 확인했다.

</div>


---

## 2.7. Lemma 2.3 — 나머지항과 Bernstein's inequality (p. 15–16)

Lemma 2.3에서는 선형 근사로 설명되지 않는 E_j의 크기를 제어한다. 이 부분은 분위수 경계 주변을 직접 분석하므로 증명이 조금 더 길다.

먼저 다음과 같이 γ를 정의한다.

<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cgamma%3D%5Cfrac54%5Cvarepsilon%5Csqrt%7B%5Calpha%7D" alt="gamma" loading="lazy">

**Lemma 2.3**의 주장은 다음과 같다.

<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20P%5Cleft%28%5Clvert%5Cmathcal%7BE%7D_j%5Crvert%3E%5Cfrac%7B3%5Cgamma%7D%7B%5Calpha%7D%5BF%5E%7B-1%7D%281-%5Calpha%2B%5Cgamma%29-F%5E%7B-1%7D%281-%5Calpha-%5Cgamma%29%5D%5Cright%29%5Cleq%5Cfrac1%7B110%7D" alt="lemma23" loading="lazy">

여기에는 실제 분위수 함수 F⁻¹가 두 번 등장한다. 분위수 수준 1−α를 중심으로 γ만큼 이동했을 때 두 분위수 값 사이의 거리가 나머지항의 상한을 결정한다.

<div markdown="1" style="background: rgba(34, 197, 94, 0.10); border-left: 4px solid #22c55e; border-radius: 8px; padding: 18px 20px; margin: 18px 0;">

**Proof**

### Step 1. 분위수 경계 양옆의 두 지점

두 점 t₀, t₁을 다음과 같이 정한다.

<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20F%28t_0%29%3D1-%5Calpha-%5Cgamma%2C%5Cquad%20F%28t_1%29%3D1-%5Calpha%2B%5Cgamma" alt="quant" loading="lazy">

t₀와 t₁은 우리가 추정하려는 VaR 부근의 왼쪽과 오른쪽 경계다. 저자들은 이 두 점에서 경험적 CDF의 오차가 γ 이내인지 확인하는 사건 A를 정의한다.

<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20A%3D%5C%7B%5Clvert%5Cwidehat%20F_%7BI_j%7D%28s%29-F%28s%29%5Crvert%5Cleq%5Cgamma%5Ctext%7B%20for%20%7Ds%5Cin%5C%7Bt_0%2Ct_1%5C%7D%5C%7D" alt="event" loading="lazy">

이를 보이기 위해 **Bernstein's inequality**를 사용한다. F̂_Ij(s)는 0 또는 1인 지시변수 m개의 평균이며, 그 분산은 F(s)(1−F(s))이다.

λ = log(440)을 대입하면 다음과 같은 상한을 얻는다.

<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cvarphi%28s%29%3D%5Csqrt%7B%5Cfrac%7B2F%28s%29%281-F%28s%29%29%5Clog%28440%29%7D%7Bm%7D%7D%2B%5Cfrac%7B2%5Clog%28440%29%7D%7B3m%7D" alt="bern" loading="lazy">

논문은 m ≥ 11/ε² 및 ε ≤ √α/6을 사용해 s = t₀, t₁에서 φ(s) ≤ γ를 계산한다. 두 경계점에 대한 실패확률을 합하면

<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20P%28A%5Ec%29%5Cleq%5Cfrac1%7B110%7D%2C%5Cquad%5Cvarphi%28t_0%29%2C%5Cvarphi%28t_1%29%5Cleq%5Cgamma" alt="bern2" loading="lazy">

가 된다.

여기에서 얻은 것은 **경계점 두 곳에서 경험적 CDF가 실제 CDF와 크게 다르지 않을 확률이 높다**는 사실이다.

### Step 2. 나머지항이 경계 사이에만 남는 이유

사건 A가 성립한다고 하자. t가 t₀보다 작은 경우 경험적 CDF와 실제 CDF는 둘 다 1−α 아래에 있다. 예를 들어

<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cwidehat%20F_%7BI_j%7D%28t%29%5Cleq%5Cwidehat%20F_%7BI_j%7D%28t_0%29%5Cleq%20F%28t_0%29%2B%5Cgamma%3D1-%5Calpha%5Cquad%28t%3Ct_0%29" alt="outside" loading="lazy">

가 성립한다. t가 t₁보다 큰 경우 두 CDF는 모두 1−α 위에 있다.

두 영역에서는 ψ의 기울기가 일정하므로 1차 근사에서 남는 항이 0이 된다. 나머지항을 제어할 적분 구간이 [t₀, t₁]로 줄어든다.

**(2.5)**

<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Clvert%5Cmathcal%7BE%7D_j%5Crvert%5Cleq%5Cfrac1%5Calpha%5Cint_%7Bt_0%7D%5E%7Bt_1%7D%5Clvert%20F%28t%29-%5Cwidehat%20F_%7BI_j%7D%28t%29%5Crvert%5C%2Cdt" alt="25" loading="lazy">

### Step 3. 경계 내부의 오차 제어

F와 F̂_Ij는 모두 단조증가한다. 사건 A에서 경계점의 오차가 γ 이내라는 사실과 F(t₁)−F(t₀) = 2γ를 이용하면

**(2.6)**

<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Clvert%20F%28t%29-%5Cwidehat%20F_%7BI_j%7D%28t%29%5Crvert%5Cleq3%5Cgamma%5Cquad%28t_0%5Cleq%20t%5Cleq%20t_1%29" alt="26" loading="lazy">

을 얻는다. 따라서 (2.5)의 적분값도 제한할 수 있다.

<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Clvert%5Cmathcal%7BE%7D_j%5Crvert%5Cleq%5Cfrac%7B3%5Cgamma%7D%7B%5Calpha%7D%28t_1-t_0%29" alt="errbound" loading="lazy">

여기서 t₁−t₀는 두 분위수 사이의 거리다. Lemma 2.3은 이 결과와 사건 A의 확률을 합쳐 증명된다.

</div>


---

## 2.8. Proof of Theorem 1.4 — 최종 오차 보장 (p. 17)

<div markdown="1" style="background: rgba(34, 197, 94, 0.10); border-left: 4px solid #22c55e; border-radius: 8px; padding: 18px 20px; margin: 18px 0;">

**Proof**

이제 앞에서 증명한 결과들을 연결할 수 있다.

- **Lemma 2.2:** 선형항 L_j가 εσ_ES보다 클 확률은 최대 1/11이다.
- **Lemma 2.3:** 나머지항 E_j는 확률 1−1/110 이상에서 3γ(t₁−t₀)/α로 제한된다.

분위수 함수 F⁻¹는 해당 구간에서 L-Lipschitz 연속이므로

<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20t_1-t_0%3DF%5E%7B-1%7D%281-%5Calpha%2B%5Cgamma%29-F%5E%7B-1%7D%281-%5Calpha-%5Cgamma%29%5Cleq2L%5Cgamma" alt="lips" loading="lazy">

가 성립한다. γ = 5ε√α/4를 대입하면

<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Clvert%5Cmathcal%7BE%7D_j%5Crvert%5Cleq%5Cfrac%7B6L%5Cgamma%5E2%7D%7B%5Calpha%7D%3D%5Cfrac%7B75%7D%7B8%7DL%5Cvarepsilon%5E2%5Cleq10L%5Cvarepsilon%5E2" alt="lips2" loading="lazy">

를 얻는다. 이 계산이 **10Lε²라는 2차 오차항이 등장하는 이유**다.

두 오차항의 확률을 합집합 부등식으로 묶으면

<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20P%5C%21%5Cleft%28%5Clvert%5Cwidehat%20T_%7BI_j%7D-%5Cmathrm%7BES%7D_%7B%5Calpha%7D%28X%29%5Crvert%3E%5Cvarepsilon%5Csigma_%7B%5Cmathrm%7BES%7D%7D%2B10L%5Cvarepsilon%5E2%5Cright%29%5Cleq%5Cfrac1%7B11%7D%2B%5Cfrac1%7B110%7D%3D%5Cfrac1%7B10%7D" alt="oneblock" loading="lazy">

이 된다. 우리가 앞에서 목표로 설정한 **(2.2)**가 성립한 것이다.

마지막으로 Lemma 2.1을 적용한다. 블록 대부분이 참값 주변 구간에 있고, β₁, β₂가 지정된 범위에 있으므로 Ŝ_N도 그 구간에 들어간다.

**Theorem 1.4의 오차 보장**은 다음과 같은 형태로 정리된다.

<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20P%5C%21%5Cleft%28%5Clvert%5Cwidehat%20S_N-%5Cmathrm%7BES%7D_%7B%5Calpha%7D%28X%29%5Crvert%5Cgeq%5Cvarepsilon%5Csigma_%7B%5Cmathrm%7BES%7D%7D%2B10L%5Cvarepsilon%5E2%5Cright%29%5Cleq%5Cexp%28-N%5Cvarepsilon%5E2%2F700%29" alt="blockprob" loading="lazy">

이때 논문은 사용 가능한 상수의 한 예로 c₂ = 10, c₃ = 1/700을 제시한다. 정리의 상수값보다 중요한 부분은 **오차확률이 exp(−cNε²) 규모로 감소**하고, 오차의 크기는 εσ_ES + O(Lε²)로 제어된다는 사실이다.

</div>


블록 하나의 추정값은 10%까지 실패할 수 있다. 그래도 블록 수가 충분하고 서로 독립이면 대부분의 블록이 정확할 가능성이 높다. 최종 추정량은 그 블록들이 만든 중앙 분위수 구간에 제한된다. 이 구조에서 유한 표본의 확률적 보장이 나온다.

---

## 2.9. Proof of Theorem 1.8 — Adversarial Robustness (p. 17)

<div markdown="1" style="background: rgba(34, 197, 94, 0.10); border-left: 4px solid #22c55e; border-radius: 8px; padding: 18px 20px; margin: 18px 0;">

**Proof**

이제 관측값 일부가 악의적으로 바뀌는 상황을 분석한다. 블록별 오차에 관한 결과는 그대로 활용하고, 변조된 관측치가 **몇 개의 블록을 손상시킬 수 있는지** 계산한다.

먼저 허용 오차를 d = εσ_ES + 10Lε²로 정하자.

원래 표본에 대해서는 확률이 높은 사건에서 부정확한 블록의 개수가 다음처럼 제한된다.

**(2.7)**

<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5C%23%5C%7Bj%3A%5Clvert%5Cwidehat%20T_%7BI_j%7D-%5Cmathrm%7BES%7D_%7B%5Calpha%7D%28X%29%5Crvert%3Ed%5C%7D%3C%5Cfrac%7B2n%7D%7B10%7D" alt="27" loading="lazy">

이제 누군가 K개의 관측값을 임의로 변조한다고 하자. 한 관측값은 소속된 블록 하나에만 영향을 줄 수 있으므로 **K개 변조로 추가로 손상될 수 있는 블록 수는 최대 K개**다.

논문에서는

<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20K%5Cleq%5Cfrac%7BN%5Cvarepsilon%5E2%7D%7B140%7D%5Cleq%5Cfrac%7Bn%7D%7B10%7D" alt="kbound" loading="lazy">

를 가정한다. 따라서 원래 부정확한 블록이 20%보다 적었고, 변조 때문에 최대 10%의 블록이 추가로 손상되어도 남는 양호한 블록은 70%보다 많다.

<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20n-0.2n-0.1n%3D0.7n" alt="70" loading="lazy">

β₁, β₂가 0.35와 0.65 사이에 있으므로 중앙 분위수 구간이 여전히 참값 주변에 위치한다. 이 때문에 최종 추정량의 절단 구조가 작동한다.

**Theorem 1.8의 결론**은 다음과 같다.

<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20P%5C%21%5Cleft%28%5Clvert%5Cwidehat%20S_N-%5Cmathrm%7BES%7D_%7B%5Calpha%7D%28X%29%5Crvert%5Cgeq%5Cvarepsilon%5Csigma_%7B%5Cmathrm%7BES%7D%7D%2B10L%5Cvarepsilon%5E2%5Cright%29%5Cleq%5Cexp%28-N%5Cvarepsilon%5E2%2F2800%29" alt="robust" loading="lazy">

원문에서 사용 가능한 상수 중 하나는 C₀ = 1/140, C₁ = 10, C₂ = 1/2800이다.

이 정리의 핵심은 **원본 관측값 일부의 변경이 전체 추정 결과를 무한정 흔들지 못하도록 블록 분위수가 경계 역할을 한다**는 점이다. 특정 블록 추정값이 크게 변하더라도 전체 블록의 중앙부가 충분히 유지되면 최종 추정값의 오차를 제어할 수 있다.

</div>


---

## 2.10. Corollary 2.4 — 변조된 데이터에서도 Consistency가 성립하는가? (p. 17–18)

Section 2의 마지막은 **확률수렴(consistency)**에 관한 결과다.

ES 관련 점근 분산이 유한하고, 분위수 함수가 1−α 주변에서 Lipschitz 연속이며, 변조된 데이터의 수 K_N이 전체 N보다 작다고 하자. 특히

- K_N / N → 0
- K_N / N은 변조된 표본의 비율이다.

논문은 다음과 같이 ε_N을 설정한다.

<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cvarepsilon_N%3D12%5Csqrt%7BK_N%2FN%7D" alt="cor" loading="lazy">

그러면 **Corollary 2.4**에 의해

<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Clvert%5Cwidehat%20S_N-%5Cmathrm%7BES%7D_%7B%5Calpha%7D%28X%29%5Crvert%5Cxrightarrow%7BP%7D0%5Cqquad%28N%5Cto%5Cinfty%29" alt="conv" loading="lazy">

가 성립한다. 표본 수가 증가하면서 **변조된 관측치의 비율이 0으로 수렴**하면, 강건 추정량도 실제 ES에 확률적으로 수렴한다는 뜻이다.

다만 이 결과는 앞에서 사용한 점근 분산의 유한성과 분위수 함수의 국소 Lipschitz 조건을 전제로 한다.

---


| 증명 단계 | 사용한 도구 | 얻은 결과 |
|---|---|---|
| Lemma 2.1 | Hoeffding's inequality | 부정확한 블록이 30% 이상일 확률을 지수적으로 제어 |
| 식 (2.4) | Fubini's theorem | ES를 CDF에 관한 적분으로 표현 |
| Lemma 2.2 | 공분산 계산, Chebyshev's inequality | 선형항의 오차확률 ≤ 1/11 |
| Lemma 2.3 | Bernstein's inequality, CDF 단조성 | 나머지항을 VaR 경계 주변으로 국한 |
| Theorem 1.4 | 두 Lemma + 블록 집중 | 유한 표본에서 지수형 오차확률 상한 확보 |
| Theorem 1.8 | 변조된 블록 개수 계산 | 적대적 변조 아래에서도 오차 상한 유지 |
| Corollary 2.4 | 변조 비율에 관한 극한 | ES 추정량의 확률수렴 |

이번 절을 읽으면서 가장 중요하게 본 것은 **각 단계가 다음 단계의 가정을 확보해 준다는 점**이다. 선형항의 분산을 구하고, 경계 주변의 나머지항을 통제하고, 블록별 결과를 모아 전체 추정량의 확률을 제어한다. 마지막으로 이 블록 구조를 이용해 데이터 변조까지 다룬다.

다음 **3편 Section 3 — On the Optimality**에서는 이렇게 얻은 오차 상한이 왜 최적이라고 할 수 있는지 살펴볼 예정이다. 그 절에서는 다른 추정량을 선택하더라도 피하기 어려운 **하한(lower bound)**을 구성한다.

---

> Bartl, D., & Eckstein, S. (2026). *Optimal Nonparametric Estimation of the Expected Shortfall Risk*. Section 2, pp. 12–18. [arXiv:2405.00357v2](https://arxiv.org/abs/2405.00357)
