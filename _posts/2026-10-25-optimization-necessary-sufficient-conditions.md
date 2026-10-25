---
title: "Optimization Theory : Necessary and sufficient optimality conditions"
date: 2026-10-25 12:00:00 +0100
categories:
  - Optimization Theory
tags:
  - Optimization
  - Optimality Conditions
  - Hessian
  - Convexity
  - Coercivity
math: false
---

이번 글에서는 최적화의 기본 문제 설정부터 **필요조건(necessary condition)**, **충분조건(sufficient condition)**, **헤세 행렬(Hessian matrix)의 정부호성**, **볼록성(convexity)**, **강제성(coercivity)**까지 정리한다.

미적분에서 배운 극값 판정법을 먼저 복습한 후, 미분법만으로 풀기 어려운 문제와 전역 최적해를 판정하는 방법을 살펴본다.

**Definition — Optimization**

최적화란 가능한 선택지 가운데 정해진 기준에 따라 가장 좋은 것을 선택하는 과정이다. 기준을 함수로 나타낼 수 있다면 그 함수를 최소화 또는 최대화하는 문제로 표현된다.

---

## 1. Foundations of Optimization from Mathematical Analysis

### 1.1 Finding Local Minima and Maxima

미적분에서 배운 기본 절차는 다음과 같다.

1. 함수의 **1계 도함수(first-order derivatives)**를 계산한다. 다변수 함수라면 각 변수에 대한 편도함수를 구한다.
2. 도함수를 0으로 놓고 **임계점(critical points)**을 구한다.
3. 1변수에서는 2계 도함수의 부호를 확인하고, 2변수에서는 헤세 행렬식과 순수 2계 편도함수의 부호를 이용해 극값을 판정한다.

#### Theorem — Second Derivative Test for Functions of Two Variables

2변수 함수에서 사용하는 판별식은 다음과 같다.

![Discriminant](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20D%3Df_%7Bxx%7Df_%7Byy%7D-%28f_%7Bxy%7D%29%5E2)

임계점에서 다음 기준을 적용한다.

- D > 0이고 fₓₓ > 0이면 국소 최솟값이다.
- D > 0이고 fₓₓ < 0이면 국소 최댓값이다.
- D < 0이면 안장점(saddle point)이다.
- D = 0이면 이 판정법만으로 결론을 내릴 수 없다.

#### Example — Finding Critical Points of a Bivariate Polynomial

다음 함수를 살펴보자.

![Function](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f%28x%2Cy%29%3Dx%5E3%2Bx%5E2y-y%5E2-4y)

먼저 편도함수를 계산한다.

![First derivatives](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f_x%3D3x%5E2%2B2xy%2C%5Cqquad%20f_y%3Dx%5E2-2y-4)

따라서 임계점은 다음 연립방정식의 해이다.

![Stationary equations](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%203x%5E2%2B2xy%3D0%2C%5Cqquad%20x%5E2-2y-4%3D0)

첫 번째 식을 인수분해하면 다음과 같다.

![Factorization](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x%283x%2B2y%29%3D0)

따라서 두 경우를 구분한다.

**Case 1 — x = 0**

두 번째 식에 대입하면 y = −2를 얻는다.

임계점은 (0, −2)이다.

**Case 2 — y = −3x/2**

두 번째 식에 대입하면 다음과 같다.

![Polynomial](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x%5E2%2B3x-4%3D%28x%2B4%29%28x-1%29%3D0)

나머지 임계점은 (−4, 6)과 (1, −3/2)이다.

이제 2계 편도함수를 계산한다.

![Second derivatives](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f_%7Bxx%7D%3D6x%2B2y%2C%5Cquad%20f_%7Byy%7D%3D-2%2C%5Cquad%20f_%7Bxy%7D%3D2x)

따라서 판별식은 다음과 같다.

![Discriminant calculation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20D%3D%286x%2B2y%29%28-2%29-%282x%29%5E2)

각 임계점에 대입한 결과는 다음과 같다.

| Critical point | fₓₓ | D | Classification |
|---|---:|---:|---|
| (0, −2) | −4 | 8 | Local maximum |
| (−4, 6) | −12 | −40 | Saddle point |
| (1, −3/2) | 3 | −10 | Saddle point |

이 예제는 1계 조건으로 후보를 찾고, 2계 조건으로 후보의 성질을 분류하는 과정을 보여준다.

하지만 도함수를 0으로 만드는 방법만으로 모든 극값을 찾을 수는 없다.

#### Example — A Nondifferentiable Local Minimum

다음 함수를 생각하자.

![Nondifferentiable function](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f%28x%2Cy%29%3D%5Csqrt%7Bx%5E2%2By%5E2%7D)

이 함수는 원점 (0, 0)에서 최솟값 0을 가진다.

그러나 원점에서는 편도함수가 존재하지 않는다. 따라서 미분 가능한 점에 대한 정지조건만 확인하면 이 극솟점을 놓칠 수 있다.

---

### 1.2 Finding Extrema on Compact Sets

콤팩트 집합 위의 연속함수는 최댓값과 최솟값을 가진다.

전역 최댓값과 최솟값을 찾으려면 다음 후보를 모두 확인해야 한다.

1. 내부의 임계점
2. 미분 불가능한 점
3. 경계에서 발생하는 극값

#### Example — Global Extrema on the Closed Unit Disk

다음 최적화 문제를 생각하자.

![Constrained function](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f%28x%2Cy%29%3Dx%5E3-y%5E2%2Bxy%5E2)

제약조건은 다음과 같다.

![Disk constraint](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x%5E2%2By%5E2%5Cleq1)

**Step 1 — Interior Critical Points**

편도함수를 계산한다.

![Derivatives](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f_x%3D3x%5E2%2By%5E2%2C%5Cqquad%20f_y%3D2y%28x-1%29)

첫 번째 편도함수가 0이 되려면 x = 0, y = 0이어야 한다.

따라서 내부 임계점은 (0, 0)이며 함수값은 0이다.

**Step 2 — Boundary Analysis**

경계에서는 다음 관계가 성립한다.

![Boundary](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y%5E2%3D1-x%5E2%2C%5Cqquad%20-1%5Cleq%20x%5Cleq1)

이를 목적함수에 대입하면 다음 1변수 함수가 나온다.

![Boundary function](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20g%28x%29%3Dx%5E2%2Bx-1)

도함수를 계산한다.

![Boundary derivative](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20g%27%28x%29%3D2x%2B1)

따라서 경계 내부의 후보는 x = −1/2이고, 이때 y = ±√3/2이다.

**Step 3 — Comparing Candidates**

| Candidate | Function value |
|---|---:|
| (0, 0) | 0 |
| (−1/2, ±√3/2) | −5/4 |
| (−1, 0) | −1 |
| (1, 0) | 1 |

**결론:** 전역 최댓값은 1이며 (1, 0)에서 달성된다. 전역 최솟값은 −5/4이며 (−1/2, ±√3/2)에서 달성된다.

---

### 1.3 A Question Leading to Convexity

미적분에서 배운 방법은 국소 극값을 찾는 데 유용하다.

그렇다면 **국소 최솟값이 언제 전역 최솟값이 되는가?**

이 질문에 답하기 위해 볼록성(convexity)을 공부한다.

---

## 2. Why Optimization Theory Is Needed

미분법으로 최적화 문제를 푸는 것이 항상 가능한 것은 아니다.

실제 문제에서는 목적함수가 복잡하거나, 제약조건이 존재하거나, 변수가 정수여야 할 수도 있다.

최적화 이론은 이러한 문제들을 체계적으로 다룬다.

#### Example — A Maximization Problem with Nonlinear Equations

다음 목적함수를 최대화한다고 하자. 단, a > 0이다.

![Objective](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f%28p%2Ca%29%3D%282000%2B%5Cln%20a-20p%29%28p-2%29-a-20000)

편도함수를 0으로 놓으면 다음 조건을 얻는다.

![First condition](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f_p%3D2040%2B%5Cln%20a-40p%3D0)

![Second condition](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f_a%3D%5Cfrac%7Bp-2%7D%7Ba%7D-1%3D0)

두 번째 식으로부터 a = p − 2를 얻는다.

이를 첫 번째 식에 대입하면 다음 비선형 방정식이 나온다.

![Nonlinear equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%202040%2B%5Cln%28p-2%29%3D40p)

이 방정식은 단순한 대수적 조작만으로 풀기 어렵다.

실제 최적화에서는 이러한 비선형 문제를 수치적으로 해결하는 방법이 필요하다.

#### Example — Challenges from Constraints

다음 문제를 생각하자.

![Function](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f%28x%2Cy%29%3Dx%5E3-y%5E2%2Bxy)

![Constraint](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x%5E2%2By%5E2%5Cleq1)

내부 임계점을 찾으려면 다음 조건을 푼다.

![Stationary equations](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f_x%3D3x%5E2%2By%3D0%2C%5Cqquad%20f_y%3D-2y%2Bx%3D0)

두 번째 식에서 x = 2y를 얻는다.

이를 첫 번째 식에 대입하면 다음과 같다.

![Reduced equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%2012y%5E2%2By%3D0)

따라서 내부 임계점은 (0, 0)과 (−1/6, −1/12)이다.

하지만 전체 최적해를 구하려면 경계도 확인해야 한다.

경계에서는 다음 관계가 성립한다.

![Boundary relation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y%3D%5Cpm%5Csqrt%7B1-x%5E2%7D)

대입한 뒤 얻는 두 개의 경계 함수를 따로 분석해야 한다.

![First boundary function](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20g_1%28x%29%3Dx%5E3%2Bx%5E2-1%2Bx%5Csqrt%7B1-x%5E2%7D)

![Second boundary function](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20g_2%28x%29%3Dx%5E3%2Bx%5E2-1-x%5Csqrt%7B1-x%5E2%7D)

제약조건이 존재하면 내부의 정지조건만으로 전체 문제를 해결하기 어려울 수 있다.

최적화 이론에서는 이러한 상황을 라그랑주 승수법이나 KKT 조건 등을 통해 체계적으로 다룬다.

#### Example — A Retail Packaging Decision

실제 생산 문제에서는 포장 크기, 재료비, 운송비, 보관비 등을 함께 고려해야 한다.

이 경우 목적함수는 비용의 최소화일 수 있지만, 동시에 부피와 규격에 대한 제약조건이 존재한다.

더 나아가 포장 단위가 정수여야 한다면 연속적인 미분법만으로 문제를 해결할 수 없다.

### 2.1 Goals and Scope of Optimization Theory

최적화 이론의 주요 목표는 다음과 같다.

1. 최적해가 존재하는지 판단한다.
2. 최적해가 만족해야 할 필요조건을 찾는다.
3. 충분조건을 이용해 후보가 실제 최적해인지 확인한다.
4. 해석적으로 풀기 어려운 문제를 수치적으로 해결한다.
5. 국소 최적해와 전역 최적해를 구분한다.

---

## 3. Necessary and Sufficient Optimality Conditions

### 3.1 First- and Second-Order Necessary Conditions

#### Theorem — Necessary Optimality Condition

n차원 실수공간의 점 x*가 무제약 국소 최솟점 또는 국소 최댓점이라고 하자.

함수 f가 x*를 포함하는 열린 근방에서 연속 미분 가능하다면 다음 조건이 반드시 성립한다.

![First order necessary condition](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cnabla%20f%28x%5E%7B%5Cast%7D%29%3D0)

이 조건을 **1계 필요조건(first-order necessary condition)**이라고 한다.

추가로 함수가 두 번 연속 미분 가능하다면, 국소 최솟점에서 헤세 행렬은 양의 준정부호여야 한다.

![Second order necessary condition minimum](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20v%5ET%5Cnabla%5E2f%28x%5E%7B%5Cast%7D%29v%5Cgeq0)

반대로 국소 최댓점에서는 헤세 행렬이 음의 준정부호여야 한다.

![Second order necessary condition maximum](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20v%5ET%5Cnabla%5E2f%28x%5E%7B%5Cast%7D%29v%5Cleq0)

두 식은 모든 방향 벡터 v에 대해 성립해야 한다.

중요한 점은 이들이 **필요조건**이라는 사실이다. 조건을 만족한다고 해서 반드시 극값인 것은 아니다.

### 3.2 Second-Order Sufficient Conditions

#### Theorem — Second-Order Sufficient Optimality Condition

함수 f가 두 번 연속 미분 가능하다고 하자.

먼저 다음 정지조건을 만족해야 한다.

![Stationary condition](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cnabla%20f%28x%5E%7B%5Cast%7D%29%3D0)

이때 헤세 행렬이 양의 정부호이면 x*는 **엄격한 국소 최솟점(strict local minimum)**이다.

반대로 헤세 행렬이 음의 정부호이면 x*는 **엄격한 국소 최댓점(strict local maximum)**이다.

필요조건과 충분조건의 차이를 구분해야 한다.

- 필요조건: 최적해라면 반드시 만족해야 하는 조건
- 충분조건: 해당 조건을 만족하면 최적성 결론을 내릴 수 있는 조건

---

## 4. Definiteness of Matrices

### 4.1 Definitions

헤세 행렬을 이용한 최적성 판정에서는 행렬의 정부호성이 중요하다.

#### Definition — Definiteness of a Symmetric Real Matrix

실수 대칭행렬 A에 대해 다음과 같이 정의한다.

| Type | Definition |
|---|---|
| Positive definite | xᵀAx > 0 for all x ≠ 0 |
| Positive semidefinite | xᵀAx ≥ 0 for all x |
| Negative definite | xᵀAx < 0 for all x ≠ 0 |
| Negative semidefinite | xᵀAx ≤ 0 for all x |

양의 정부호는 모든 0이 아닌 방향에서 이차형식이 양수라는 뜻이다.

양의 준정부호에서는 이차형식이 0이 되는 방향이 존재할 수 있다.

이 차이가 2계 충분조건과 필요조건을 구분하는 데 중요하다.

### 4.2 Practical Tests for Definiteness

#### Theorem — Sylvester's Criterion

대칭행렬의 정부호성을 판단할 때는 **선도 주인수식(leading principal minors)**을 사용할 수 있다.

양의 정부호이기 위한 조건은 모든 선도 주인수식이 양수인 것이다.

음의 정부호이기 위한 조건은 선도 주인수식의 부호가 다음처럼 번갈아 나타나는 것이다.

- 홀수 차수: 음수
- 짝수 차수: 양수

#### Theorem — Eigenvalue Criterion

실수 대칭행렬의 고유값을 이용할 수도 있다.

- 모든 고유값이 양수이면 양의 정부호
- 모든 고유값이 음수이면 음의 정부호
- 모든 고유값이 0 이상이면 양의 준정부호
- 모든 고유값이 0 이하이면 음의 준정부호

주의할 점은 준정부호를 검사할 때 선도 주인수식만으로 충분하지 않을 수 있다는 것이다.

준정부호를 주인수식으로 판정하려면 **모든 주인수식(principal minors)**을 확인해야 한다.

#### Example — Testing Negative Definiteness of a 4×4 Matrix

다음 행렬을 생각하자.

![Matrix](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20A%3D%5Cbegin%7Bpmatrix%7D-2%26-1%261%261%5C%5C-1%26-4%260%26-2%5C%5C1%260%26-2%260%5C%5C1%26-2%260%26-4%5Cend%7Bpmatrix%7D)

선도 주인수식을 계산하면 다음과 같은 부호가 나온다.

![Principal minors](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5CDelta_1%3D-2%3C0%2C%5Cqquad%5CDelta_2%3D7%3E0)

![Third minor](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5CDelta_3%3D-10%3C0%2C%5Cqquad%5CDelta_4%3D12%3E0)

부호가 음수, 양수, 음수, 양수로 번갈아 나타난다.

따라서 이 행렬은 음의 정부호이다.

---

## 5. Applying Optimality Conditions to Functions

#### Example — A Local Minimum of a Three-Variable Function

다음 함수를 생각하자.

![Three-variable function](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f%28x%2Cy%2Cz%29%3D%28x%5E2-z%29%5E2%2B%28y-z%29%5E2%2B2z)

**Step 1 — Finding Stationary Points**

편도함수를 계산하고 모두 0으로 놓는다.

![Gradient system](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cbegin%7Bcases%7Df_x%3D4x%28x%5E2-z%29%3D0%5C%5Cf_y%3D2%28y-z%29%3D0%5C%5Cf_z%3D-2%28x%5E2-z%29-2%28y-z%29%2B2%3D0%5Cend%7Bcases%7D)

두 번째 식에서 y = z를 얻는다.

이를 세 번째 식에 대입하면 z = x² − 1이다.

첫 번째 식에 다시 대입하면 x = 0이므로 정지점은 (0, −1, −1)이다.

**Step 2 — Computing the Hessian**

헤세 행렬은 다음과 같다.

![Hessian](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cnabla%5E2f%3D%5Cbegin%7Bpmatrix%7D12x%5E2-4z%260%26-4x%5C%5C0%262%26-2%5C%5C-4x%26-2%264%5Cend%7Bpmatrix%7D)

정지점에 대입한다.

![Hessian at stationary point](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cnabla%5E2f%280%2C-1%2C-1%29%3D%5Cbegin%7Bpmatrix%7D4%260%260%5C%5C0%262%26-2%5C%5C0%26-2%264%5Cend%7Bpmatrix%7D)

선도 주인수식은 다음과 같다.

![Principal minors](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5CDelta_1%3D4%2C%5Cqquad%5CDelta_2%3D8%2C%5Cqquad%5CDelta_3%3D16)

모두 양수이므로 헤세 행렬은 양의 정부호이다.

따라서 (0, −1, −1)은 **엄격한 국소 최솟점**이다.

그렇다면 이 점은 전역 최솟점이기도 할까?

### 5.1 Why Local Optimality Is Not Enough

#### Example — A Global Minimum Without the Sufficient Condition

다음 함수를 보자.

![Function](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f%28x%2Cy%29%3D%28x%5E2%2By%5E2-1%29%5E2)

제곱은 항상 0 이상이므로 함수값도 0 이상이다.

또한 다음 조건을 만족하는 모든 점에서 함수값이 0이다.

![Unit circle](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x%5E2%2By%5E2%3D1)

따라서 단위원 위의 모든 점이 전역 최솟점이다.

이제 (1, 0)에서 2계 충분조건을 확인해 보자.

먼저 편도함수를 계산한다.

![First derivatives](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f_x%3D4x%28x%5E2%2By%5E2-1%29%2C%5Cquad%20f_y%3D4y%28x%5E2%2By%5E2-1%29)

이 점에서 그래디언트는 0이다.

하지만 헤세 행렬은 다음과 같다.

![Hessian](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cnabla%5E2f%281%2C0%29%3D%5Cbegin%7Bpmatrix%7D8%260%5C%5C0%260%5Cend%7Bpmatrix%7D)

이 행렬은 양의 준정부호이지만 양의 정부호는 아니다.

따라서 2계 충분조건을 만족하지 않는다.

그런데도 (1, 0)은 분명히 전역 최솟점이다.

**충분조건이 성립하지 않는다고 해서 최적해가 존재하지 않는 것은 아니다.**

#### Example — A Strict Local Minimum Without a Global Minimum

다음 함수를 생각하자.

![Function](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f%28x%2Cy%29%3Dx%5E2%2By%5E2%281-x%29%5E3)

원점 (0, 0)은 엄격한 국소 최솟점이다.

그러나 x가 충분히 커지고 y도 적절하게 선택되면 함수값은 매우 작아질 수 있다.

실제로 x = 2를 고정하면 다음과 같다.

![Counterexample](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f%282%2Cy%29%3D4-y%5E2)

y의 크기가 증가하면 함수값은 음의 무한대로 감소한다.

따라서 이 함수는 전역 최솟값을 가지지 않는다.

이 예제는 **엄격한 국소 최솟점이라도 전역 최솟점일 필요는 없다**는 사실을 보여준다.

---

## 6. Convexity

### 6.1 Convex and Concave Functions

#### Definition — Convex and Concave Functions

집합 X가 볼록집합이고, 그 위에 정의된 함수 f를 생각하자.

서로 다른 두 점 x, y와 0과 1 사이의 실수 λ에 대해 다음 관계가 성립하면 f를 볼록함수(convex function)라고 한다.

![Convexity](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f%28%5Clambda%20x%2B%281-%5Clambda%29y%29%5Cleq%5Clambda%20f%28x%29%2B%281-%5Clambda%29f%28y%29)

반대로 다음 부등식을 만족하면 오목함수(concave function)라고 한다.

![Concavity](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f%28%5Clambda%20x%2B%281-%5Clambda%29y%29%5Cgeq%5Clambda%20f%28x%29%2B%281-%5Clambda%29f%28y%29)

직관적으로 볼록함수에서는 그래프 위의 두 점을 잇는 선분이 함수 그래프보다 위에 놓인다.

부등식이 서로 다른 두 점에 대해 엄격하게 성립하면 엄격한 볼록함수(strictly convex function)라고 한다.

### 6.2 Testing Convexity Using the Hessian

#### Theorem — Characterizations of Convexity

함수가 두 번 연속 미분 가능하다고 하자.

볼록함수 여부를 확인하는 데 다음 조건들을 사용할 수 있다.

1. 헤세 행렬이 모든 점에서 양의 준정부호이면 함수는 볼록하다.
2. 헤세 행렬이 모든 점에서 양의 정부호이면 함수는 엄격하게 볼록하다.
3. 그래디언트는 볼록함수의 1계 최적성 조건과 연결된다.

대표적인 2계 판정 조건은 다음과 같다.

![Hessian convexity](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cnabla%5E2f%28x%29%5Csucceq0)

여기서 위 기호는 헤세 행렬이 양의 준정부호임을 뜻한다.

주의할 점은 엄격한 볼록함수라고 해서 헤세 행렬이 모든 점에서 반드시 양의 정부호인 것은 아니라는 사실이다.

### 6.3 Convexity and Optimality

#### Theorem — Convexity and Optimality

f가 볼록함수이고 정의역도 볼록하다고 하자.

그러면 다음 성질을 얻는다.

**(a)** 모든 국소 최솟점은 전역 최솟점이다.

**(b)** f가 미분 가능하다면, 정의역 내부에서 그래디언트가 0인 점은 전역 최솟점이다.

**(c)** f가 엄격하게 볼록하다면 전역 최솟점은 최대 하나이다.

볼록성은 국소적인 정보만으로 전역 최적성을 판단할 수 있게 하는 중요한 도구다.

---

## 7. Coercivity and the Existence of Global Minima

### 7.1 Definition of Coercivity

#### Definition — Coercive Function

함수 f의 입력 벡터 크기가 무한히 커질 때 함수값도 양의 무한대로 증가하면, f를 **강제적 함수(coercive function)**라고 한다.

![Coercivity](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Clim_%7B%5C%7Cx%5C%7C%5Cto%5Cinfty%7Df%28x%29%3D%2B%5Cinfty)

이는 입력 벡터가 어느 방향으로든 무한히 멀어지면 함수값도 증가한다는 뜻이다.

#### Theorem — Existence of a Global Minimizer

함수 f가 연속이고 강제적이라면 전역 최솟점이 존재한다.

강제성은 입력값이 무한대로 멀어지는 방향에서 함수값이 계속 감소하는 현상을 방지한다.

다만 **전역 최솟점의 존재와 유일성은 서로 다른 문제**다.

강제성은 존재성을 보장할 수 있지만, 유일성을 자동으로 보장하지 않는다.

### 7.2 Revisiting the Three-Variable Example

#### Example — A Nonconvex but Coercive Function

앞서 살펴본 다음 함수를 다시 생각하자.

![Function](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f%28x%2Cy%2Cz%29%3D%28x%5E2-z%29%5E2%2B%28y-z%29%5E2%2B2z)

앞에서는 (0, −1, −1)이 엄격한 국소 최솟점이라는 사실을 확인했다.

이제 전역 최솟점인지 확인한다.

**Step 1 — Testing Convexity**

헤세 행렬은 다음과 같다.

![Hessian](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cnabla%5E2f%3D%5Cbegin%7Bpmatrix%7D12x%5E2-4z%260%26-4x%5C%5C0%262%26-2%5C%5C-4x%26-2%264%5Cend%7Bpmatrix%7D)

예를 들어 (0, 0, 1)에서는 첫 번째 대각성분이 −4이다.

따라서 이 점에서 헤세 행렬은 양의 준정부호가 아니다.

그러므로 f는 전역적으로 볼록하지 않다.

**Step 2 — Testing Coercivity**

함수를 다음과 같은 형태로 다시 쓸 수 있다.

![Equivalent form 1](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f%3D%28x%5E2-z%29%5E2%2B%28y-z%29%5E2%2B2z)

![Equivalent form 2](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f%3D%28x%5E2-z%29%5E2%2B%28y-z-1%29%5E2%2B2y-1)

![Equivalent form 3](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f%3D%28x%5E2-z-1%29%5E2%2B2x%5E2-1%2B%28y-z%29%5E2)

이 세 가지 표현은 서로 동등하다.

원자료에서는 입력 벡터가 무한히 커지는 여러 경우를 구분한다.

- x의 크기가 증가하는 경우
- y의 크기가 증가하는 경우
- z의 크기가 증가하는 경우
- 여러 변수가 동시에 증가하는 경우

각 경우에 제곱항과 나머지 항의 증가를 비교하여 강제성을 확인할 수 있다.

따라서 함수는 강제적이며 전역 최솟점이 존재한다.

**Step 3 — Proving the Global Minimum Directly**

세 번째 표현을 다시 보자.

![Sum of squares](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f%28x%2Cy%2Cz%29%3D%28x%5E2-z-1%29%5E2%2B%28y-z%29%5E2%2B2x%5E2-1)

제곱항들은 모두 0 이상이다.

따라서 다음 하한을 얻는다.

![Lower bound](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f%28x%2Cy%2Cz%29%5Cgeq-1)

등호가 성립하려면 세 개의 비음수 항이 모두 0이어야 한다.

이를 만족하는 점은 (0, −1, −1)뿐이다.

따라서 다음 결론을 얻는다.

![Global minimum](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cmin%20f%3D-1)

전역 최솟점은 (0, −1, −1)이며 유일하다.

이 예제는 중요한 사실을 보여준다.

함수가 전역적으로 볼록하지 않아도 전역 최솟점이 존재할 수 있다.

또한 최적성 조건을 직접 검사하는 것보다 제곱완성을 통해 하한을 구하는 방식이 더 간단할 수도 있다.

---

## 8. Key Vocabulary

| English | Korean | Meaning |
|---|---|---|
| Optimization | 최적화 | 목적함수를 최소화 또는 최대화하는 과정 |
| Objective function | 목적함수 | 최적화 대상 함수 |
| Feasible set | 가능집합 | 제약조건을 만족하는 해의 집합 |
| Constraint | 제약조건 | 가능한 해를 제한하는 조건 |
| Critical point | 임계점 | 도함수가 0이거나 존재하지 않는 점 |
| Stationary point | 정지점 | 그래디언트가 0인 점 |
| Gradient | 그래디언트 | 1계 편도함수 벡터 |
| Hessian matrix | 헤세 행렬 | 2계 편도함수 행렬 |
| Necessary condition | 필요조건 | 최적점이라면 반드시 만족하는 조건 |
| Sufficient condition | 충분조건 | 만족하면 최적성 결론이 성립하는 조건 |
| Positive definite | 양의 정부호 | xᵀAx > 0 (x ≠ 0) |
| Positive semidefinite | 양의 준정부호 | xᵀAx ≥ 0 |
| Negative definite | 음의 정부호 | xᵀAx < 0 (x ≠ 0) |
| Negative semidefinite | 음의 준정부호 | xᵀAx ≤ 0 |
| Local minimizer | 국소 최솟점 | 근방에서 가장 작은 함수값을 갖는 점 |
| Global minimizer | 전역 최솟점 | 정의역 전체에서 가장 작은 함수값을 갖는 점 |
| Saddle point | 안장점 | 극솟점이나 극댓점이 아닌 정지점 |
| Convex function | 볼록함수 | 선분 위 함수값이 양 끝 함수값의 선형보간 이하 |
| Concave function | 오목함수 | 볼록함수의 반대 부등식을 만족 |
| Strictly convex | 엄격한 볼록함수 | 서로 다른 두 점 사이에서 엄격한 볼록부등식 성립 |
| Coercive | 강제적 | 입력 벡터가 무한히 멀어지면 함수값도 양의 무한대로 증가 |
| Leading principal minor | 선도 주인수식 | 왼쪽 위 정사각 부분행렬의 행렬식 |
| Eigenvalue criterion | 고유값 판정법 | 고유값의 부호로 정부호성을 판정 |
| Sylvester's criterion | 실베스터 판정법 | 선도 주인수식의 부호를 이용한 정부호 판정 |
| Lagrange multipliers | 라그랑주 승수법 | 등식 제약 최적화에 사용하는 방법 |
| KKT conditions | KKT 조건 | 부등식 제약 최적화의 최적성 조건 |

## Reference

Piotr Więcek, *Optimization Theory — Necessary and Sufficient Optimality Conditions*, Wrocław University of Science and Technology, 2026.
