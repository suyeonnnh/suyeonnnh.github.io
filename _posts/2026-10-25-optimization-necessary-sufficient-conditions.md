---
title: 'Optimization Theory : Necessary and sufficient optimality conditions'
date: 2026-10-25 12:00:00 +0200
categories: [Optimization Theory]
tags: [Optimization, Optimality Conditions, Hessian, Convexity, Coercivity]
math: false
---

# Optimization Theory : Necessary and sufficient optimality conditions

이번 글에서는 최적화의 기본 문제 설정부터 **necessary condition**, **sufficient condition)**, **헤세 행렬(Hessian matrix)의 정부호성**, **convexity**, **coercivity**까지 정리한다. 미적분에서 배운 극값 판정법을 먼저 복습한 후, 미분법만으로 풀기 어려운 문제와 전역 최적해를 판정하는 방법을 살펴본다.


> **Definition | Optimization**  
> 가능한 선택지 가운데 정해진 기준에 따라 가장 좋은 것을 선택하는 과정. 기준을 함수로 나타낼 수 있다면 그 함수를 **최소화 또는 최대화**하는 문제로 표현된다.
{: .prompt-tip }

## 1. Foundations of Optimization from Mathematical Analysis

### 1.1 Finding Local Minima and Maxima

미적분에서 배운 기본 절차는 다음과 같다.

1. 함수의 **1계 도함수(first-order derivatives)**, 다변수 함수라면 편도함수를 계산한다.
2. 도함수를 0으로 놓고 **임계점(critical points)**을 구한다.
3. 1변수에서는 2계 도함수의 부호를 확인하고, 2변수에서는 헤세 행렬식과 순수 2계 편도함수의 부호를 확인해 극값의 종류를 판정한다.

> **Theorem | Second Derivative Test for Functions of Two Variables**  
>
> **Formula:** `D=f_xxf_yy-(f_xy)^2.`
>
> 임계점에서 `D>0, f_xx>0`이면 **국소 최솟값**, `D>0, f_xx<0`이면 **국소 최댓값**이다. `D<0`이면 **안장점(saddle point)**이고, `D=0`이면 이 판정법만으로 결론을 내릴 수 없다.
{: .prompt-tip }

> **Example | Finding Critical Points of a Bivariate Polynomial**  
>
> **Formula:** `f(x,y)=x^3+x^2y-y^2-4y.`
>
> 먼저 편도함수를 계산한다.  
>
> **Formula:** `f_x=3x^2+2xy, f_y=x^2-2y-4.`
>
> 따라서 임계점은 다음 연립방정식의 해이다.  
>
> **Formula:** `3x^2+2xy=0, x^2-2y-4=0.`
>
> 첫 번째 식을 인수분해하면  
>
> **Formula:** `x(3x+2y)=0 ⇒ x=0 or y=-(3/2)x.`
>
> **Case 1: `x=0`.** 두 번째 식에 대입하면 `-2y-4=0`이므로 `(0,-2)`를 얻는다.  
> **Case 2: `y=-(3/2)x`.** 두 번째 식에 대입하면  
>
> **Formula:** `x^2+3x-4=(x+4)(x-1)=0.`
>
> 따라서 나머지 임계점은 `(-4,6)`과 `(1,-3/2)`이다.  
> 이제 2계 편도함수를 구한다.  
>
> **Formula:** `f_xx=6x+2y, f_yy=-2, f_xy=2x,`
>
>
> **Formula:** `D=(6x+2y)(-2)-(2x)^2=-12x-4y-4x^2.`
>
> 각 점에서의 결과는 다음과 같다.  
> 
>
> - **임계점** — `f_xx`, `D`, Classification
> - **`(0,-2)`** — `-4`, `8>0`, 국소 최댓값
> - **`(-4,6)`** — `-12`, `-40<0`, 안장점
> - **`(1,-3/2)`** — `3`, `-10<0`, 안장점
{: .prompt-info }

하지만 **도함수를 0으로 만드는 방식으로 모든 극값을 찾을 수는 없다.** 극값이 존재해도 해당 점에서 미분이 불가능할 수 있다.

> **Example | A Nondifferentiable Local Minimum**  
>
> **Formula:** `f(x,y)=√(x^2+y^2).`
>
> 이 함수는 `(0,0)`에서 최솟값 0을 가지지만, 원점에서 편도함수가 존재하지 않는다. 따라서 미분 가능한 점의 정지조건만 검사하면 이 극솟점을 놓칠 수 있다.
{: .prompt-info }

### 1.2 Finding Extrema on Compact Sets

콤팩트 집합 위에서 연속함수의 최댓값과 최솟값을 찾을 때는 **내부의 임계점, 미분 불가능한 점, 경계의 극값 후보**를 모두 조사하고 함수값을 비교한다. 다른 방법으로는 라그랑주 승수법(Lagrange multipliers)이 있다.

> **Example | Global Extrema on the Closed Unit Disk**  
>
> **Formula:** `f(x,y)=x^3-y^2+xy^2, x^2+y^2≤1.`
>
> **① Interior Critical Points**  
>
> **Formula:** `f_x=3x^2+y^2=0, f_y=-2y+2xy=2y(x-1)=0.`
>
> 첫 번째 식이 0이려면 실수 범위에서 `x=y=0`이어야 한다. 따라서 내부 임계점은 `(0,0)`이고 `f(0,0)=0`이다. 두 번째 식을 기준으로 나누어 보더라도 `y=0`이면 `x=0`, `x=1`이면 `3+y^2=0`이라는 모순이 나온다.  
> **② Boundary Analysis**  
> 경계 `x^2+y^2=1`에서는 `y^2=1-x^2`, `-1≤ x≤1`이다.  
>
> **Formula:** `g(x) =x^3-(1-x^2)+x(1-x^2); =x^2+x-1.`
>
> 도함수를 0으로 놓으면  
>
> **Formula:** `g'(x)=2x+1=0 ⇒ x=-1/2,`
>
>
> **Formula:** `y=±√(1-1/4)=±(\sqrt3)/(2).`
>
> 따라서 경계 임계점에서 `f(-1/2,±(\sqrt3)/(2))=-5/4`이다. 또한 구간 끝점에 대응하는 `(-1,0)`과 `(1,0)`에서 함수값은 각각 `-1`, `1`이다.  
> **③ Comparing Candidate Values**  
> 
>
> - **Candidate** — Function Value
> - **`(0,0)`** — `0`
> - **`(-1/2,±(\sqrt3)/(2))`** — `-5/4`
> - **`(-1,0)`** — `-1`
> - **`(1,0)`** — `1`
> 
> **전역 최댓값은 `1` (위치 `(1,0)`), 전역 최솟값은 `-5/4` (위치 `(-1/2,±(\sqrt3)/(2))`)**이다.
{: .prompt-info }

### 1.3 An Introduction to Convexity

미적분에서는 볼록함수의 정의와, 1변수에서 이를 실용적으로 판단하는 방법도 배운다. 이 개념은 나중에 **국소 최솟값을 전역 최솟값으로 확장할 수 있는지**를 판단하는 중요한 도구가 된다.

## 2. Why Do We Need Optimization Theory?

미분법으로 극값의 후보를 찾는 것은 출발점이지만, 실제 문제에서는 방정식이 어렵거나 제약조건이 있거나 변수가 정수여야 한다. 최적화 이론은 이런 경우를 체계적으로 다룬다.

> **Example | Maximization Leading to a Nonlinear Equation**  
> 다음 목적함수를 최대화한다고 하자. (`a>0`)  
>
> **Formula:** `f(p,a)=(2000+ln a-20p)(p-2)-a-20000.`
>
> 1계 조건은  
>
> **Formula:** `f_p=2000+ln a+40-40p=0,`
>
>
> **Formula:** `f_a=(p-2)/(a)-1=0.`
>
> 두 번째 조건으로부터 `a=p-2`를 얻고, 이를 대입하면  
>
> **Formula:** `2039+ln(p-2)=40p`
>
> 와 같은 비선형 방정식이 나온다. 이 식은 단순한 대수적 조작으로 풀기 어렵다. 일반적인 최적화에서는 수치적으로도 풀기 어려운 방정식을 자주 만난다.  
> *주:* 원자료의 중간 상수 표기에 산술상 차이가 있다. `f_p`를 직접 미분하면 상수항은 `2040`이므로, 정확한 치환식은 **`2040+ln(p-2)=40p`**이다. 위의 `2039`는 원자료에 제시된 식이다.
{: .prompt-info }

> **Example | Computational Challenges in Constrained Optimization**  
>
> **Formula:** `f(x,y)=x^3-y^2+xy, x^2+y^2≤1.`
>
> 내부 임계점은  
>
> **Formula:** `f_x=3x^2+y=0, f_y=-2y+x=0`
>
> 을 풀어 얻는다. `x=2y`를 대입하면  
>
> **Formula:** `12y^2+y=y(12y+1)=0,`
>
> 따라서 `(0,0)`과 `(-1/6,-\frac112)`가 나온다.  
> 경계 `x^2+y^2=1`에서는  
>
> **Formula:** `y=±√(1-x^2), -1≤ x≤1`
>
> 이므로 다음 두 함수의 임계점을 조사해야 한다.  
>
> **Formula:** `g_1(x)=x^3+x^2-1+x√(1-x^2),`
>
>
> **Formula:** `g_2(x)=x^3+x^2-1-x√(1-x^2).`
>
> 첫 번째 함수의 도함수는  
>
> **Formula:** `g_1'(x)=3x^2+2x+√(1-x^2)-(x^2)/(√(1-x^2)).`
>
> `g_1'(x)=0`에서 제곱근을 제거해 정리하면 다음과 같은 높은 차수의 식으로 이어진다.  
>
> **Formula:** `9x^6+12x^5-x^4-12x^3-8x^2+1=0.`
>
> 식을 수치적으로 풀 수 있더라도 부호 검증, 경계 끝점, 두 번째 곡선까지 고려해야 한다. 즉 **경계를 직접 매개화하는 방식은 계산이 매우 번거로워질 수 있다.**
{: .prompt-info }

> **Example | A Retailer’s Package Purchasing Problem**  
> 소매업자는 상품을 100개, 50개, 20개, 10개 묶음으로 구매하며 각 묶음 가격은 각각 USD 160, USD 90, USD 45, USD 28이다. 상품은 낱개당 USD 5에 판매하고, 주문량은 3,791개다. 각 포장 단위를 몇 묶음씩 구매해야 이익이 최대가 될까?  
> 포장별 구매 묶음 수를 `x_100,x_50,x_20,x_10`, 실제 판매 개수를 `u`라고 하면 목적함수는  
>
> **Formula:** `max 5u-160x_100-90x_50-45x_20-28x_10`
>
> 이며 제약조건은  
>
> **Formula:** `u≤3791,`
>
>
> **Formula:** `100x_100+50x_50+20x_20+10x_10≥ u,`
>
>
> **Formula:** `x_100,x_50,x_20,x_10≥0.`
>
> 구매 묶음 수는 실제로 **음이 아닌 정수**여야 한다. 선형 목적함수는 일반적으로 도함수가 0이 되지 않으므로 이 경우 미분만으로 답을 구할 수 없다. 연속 선형계획의 최적해는 적절한 조건하에서 가능영역의 꼭짓점에서 찾을 수 있지만, 변수 수에 따라 꼭짓점의 개수가 기하급수적으로 증가할 수 있다. 더구나 정수 제약은 별도로 처리해야 한다.
{: .prompt-info }

### 2.1 Goals and Scope of Optimization Theory

이 강의에서 강조하는 세 가지 목표는 다음과 같다.

1. 도함수를 활용해 최솟점과 최댓점을 찾는 **수치적 도구**를 제공한다.
2. **제약 최적화 문제**를 미분 기반 방법으로 다룰 수 있도록 변형하는 분석 도구를 제공한다.
3. 미분 기반 접근이 적합하지 않은 **선형·이산 최적화** 문제를 효율적으로 해결하는 방법을 제공한다.

전체 학습 주제는 (1) 최적성 필요·충분조건, (2) 무제약 최적화를 위한 경사 기반 방법, (3) 선형계획법과 정수 선형계획법, (4) 라그랑주 승수법과 Karush–Kuhn–Tucker(KKT) 조건, (5) 이차계획법, (6) 일반적인 제약 최적화 방법, (7) 이산 최적화를 위한 휴리스틱 방법이다.

> **Vocabulary | Basic Optimization Terminology**  
>
> **Formula:** `min_x∈ Df(x) or max_x∈ Df(x).`
>
> - **Objective function(목적함수)**: 최소화·최대화할 대상인 `f`; 최소화 시 *cost*, 최대화 시 *reward*라고도 한다.  
> - **Variable(변수)**: 선택하거나 조정하는 `x`.  
> - **Feasible set(가능집합)**: 허용되는 선택의 집합 `D`.  
> - **Constraint(제약조건)**: `D`를 정의하는 등식 또는 부등식.  
> - **`argmin` / `argmax`**: 가장 작은/큰 **함수값을 만들어 내는 변수(해)**.  
> - **Minimizer / Maximizer(최솟점 / 최댓점)**: 최적해를 이루는 점.  
> - **Constrained / Unconstrained**: 제약 / 무제약.  
> - **Local / Global**: 국소 / 전역.  
> 주의: `min f(x)`는 최적 **값**이고, `argminf(x)`는 그 값을 달성하는 **점 또는 점들의 집합**이다.
{: .prompt-tip }

## 3. Necessary and Sufficient Optimality Conditions

### 3.1 First- and Second-Order Necessary Conditions

> **Theorem | Necessary Optimality Conditions**  
> `x^*∈ℝ^n`이 **무제약 국소 최솟점 또는 최댓점**이고 `f`가 `x^*`를 포함하는 열린 근방에서 연속 미분 가능하다고 하자. 그러면 반드시  
>
> **Formula:** `∇ f(x^*)=0`
>
> 이다. 여기에 `f`가 2번 연속 미분 가능하다는 조건까지 추가하면,  
> - 국소 최솟점에서 `∇^2f(x^*)`는 **양의 준정부호(positive semidefinite)**이고,
> - 국소 최댓점에서 `∇^2f(x^*)`는 **음의 준정부호(negative semidefinite)**이다.  
> 즉, 국소 극값이라면 반드시 만족해야 하는 조건이지만, 이것만으로 극값임이 보장되지는 않는다.
{: .prompt-tip }

여기서 `∇ f`는 모든 1계 편도함수를 모은 **그래디언트(gradient)**이고 `∇^2f`는 모든 2계 편도함수로 구성한 **헤세 행렬(Hessian)**이다.

### 3.2 Second-Order Sufficient Conditions

> **Theorem | Second-Order Sufficient Optimality Conditions**  
> `f`가 열린 근방에서 두 번 연속 미분 가능하고 `∇ f(x^*)=0`이라고 하자.  
> - `∇^2 f(x^*)`가 **양의 정부호(positive definite)**이면, `x^*`는 **엄격한 국소 최솟점(strict local minimizer)**이다.
> - `∇^2 f(x^*)`가 **음의 정부호(negative definite)**이면, `x^*`는 **엄격한 국소 최댓점(strict local maximizer)**이다.  
> 준정부호인 경우에는 이 충분조건만으로 엄격한 극값을 확정할 수 없다.
{: .prompt-tip }

**필요조건**은 최적점이 *반드시 만족해야 할 조건*이고, **충분조건**은 해당 조건이 성립하면 *최적점임을 결론 내릴 수 있는 조건*이다. 필요조건을 통과해도 안장점일 수 있지만, 충분조건을 만족하면 그 점의 국소 최적성은 보장된다.

## 4. Definiteness of Matrices

### 4.1 Definition

> **Definition | Definiteness of a Symmetric Real Matrix**  
> `A`를 `n× n` 대칭 실수행렬이라고 하자.  
> - **양의 정부호(positive definite)**: `x≠0`인 모든 `x∈ℝ^n`에 대해 `x^TAx>0`.
> - **양의 준정부호(positive semidefinite)**: 모든 `x∈ℝ^n`에 대해 `x^TAx≥0`.
> - **음의 정부호(negative definite)**: `x≠0`인 모든 `x∈ℝ^n`에 대해 `x^TAx<0`.
> - **음의 준정부호(negative semidefinite)**: 모든 `x∈ℝ^n`에 대해 `x^TAx≤0`.
{: .prompt-tip }

### 4.2 Practical Tests for Definiteness

> **Theorem | Sylvester’s Criterion and the Eigenvalue Criterion**  
> **Sylvester’s Criterion — Strict Definiteness**  
> 선도 주인수식(leading principal minors) `Δ_k=det(A_1:k,1:k)`을 계산한다.  
> - **양의 정부호** `⇔` 모든 `Δ_k>0`.
> - **음의 정부호** `⇔` `(-1)^kΔ_k>0` (홀수 차수는 음수, 짝수 차수는 양수).  
> **Note on Semidefiniteness:** 원자료는 선도 주인수식의 부호만으로 준정부호를 판정할 수 있는 것처럼 적고 있지만, 일반적으로는 **모든 주인수식(principal minors)**을 확인해야 한다. 양의 준정부호는 모든 주인수식 `≥0`, 음의 준정부호는 모든 `k`차 주인수식에 대해 `(-1)^kdet≥0`이다.  
> **Eigenvalue Criterion**  
> - 양의 정부호: 모든 고유값 `0`; 양의 준정부호: 모든 고유값 `≥0`.
> - 음의 정부호: 모든 고유값 `<0`; 음의 준정부호: 모든 고유값 `≤0`.
{: .prompt-tip }

> **Example | Testing Negative Definiteness of a 4×4 Matrix**  
>
> **Formula:** `A=[-2 -1 1 1; -1 -4 0 -2; 1 0 -2 0; 1 -2 0 -4].`
>
> 선도 주인수식 네 개를 계산하면  
>
> **Formula:** `Δ_1=-2<0, Δ_2=det[-2 -1; -1 -4]=7>0,`
>
>
> **Formula:** `Δ_3=det[-2 -1 1; -1 -4 0; 1 0 -2]=-10<0,`
>
>
> **Formula:** `Δ_4=det(A)=12>0.`
>
> 부호가 음·양·음·양으로 번갈아 나타나므로 `A`는 **음의 정부호(negative definite)**이다.
{: .prompt-info }

어떤 선도 주인수식이 **0**이면 엄격한 정부호 판정은 통과하지 못하며, 준정부호를 판정하려면 빠진 주인수식을 추가로 확인해야 한다. 4×4 대칭행렬에서는 선도 주인수식 4개 외에 확인할 수 있는 비선도 주인수식이 **11개** 더 있다(전체 비어 있지 않은 주인수식 15개). 3×3에서는 1×1 비선도 주인수식 2개와 2×2 비선도 주인수식 2개가 더 있다. 2×2에서는 선도식 외에 다른 1×1 주인수식 하나도 별도로 고려해야 하므로, **원자료의 '2×2에서는 추가 확인이 필요 없다'는 문장은 일반적인 준정부호 판정에 대해 성립하지 않는다.** 예를 들어 `diag(0,-1)`의 두 선도식은 모두 0이지만 양의 준정부호는 아니다.

## 5. Applying Optimality Conditions to Functions

> **Example | A Local Minimum of a Three-Variable Function**  
>
> **Formula:** `f(x,y,z)=(x^2-z)^2+(y-z)^2+2z.`
>
> **① Finding Stationary Points**  
>
> **Formula:** `f_x=4x(x^2-z)=0,; f_y=2(y-z)=0,; f_z=-2(x^2-z)-2(y-z)+2=0.`
>
> `y=z`를 세 번째 식에 대입하면 `x^2-z=1`이다. 따라서 첫 번째 식으로부터 `x=0`, 그리고 `z=y=-1`을 얻는다. 정지점은 **`(0,-1,-1)` 하나**다.  
> **② Computing the Hessian Matrix**  
>
> **Formula:** `∇^2f(x,y,z)=[12x^2-4z 0 -4x; 0 2 -2; -4x -2 4].`
>
> 정지점에서  
>
> **Formula:** `∇^2f(0,-1,-1)=[4 0 0; 0 2 -2; 0 -2 4].`
>
> 이 행렬의 선도 주인수식은 `Δ_1=4`, `Δ_2=8`, `Δ_3=16`으로 모두 양수다. 따라서 양의 정부호이고, **`(0,-1,-1)`은 엄격한 국소 최솟점**이다.  
> 하지만 이것이 **전역 최솟점**인지 알려면 추가적인 논증이 필요하다.
{: .prompt-info }

### 5.1 Why Local Optimality Is Not Enough

> **Example | A Global Minimizer Without the Second-Order Sufficient Condition**  
>
> **Formula:** `f(x,y)=(x^2+y^2-1)^2.`
>
> 제곱은 항상 0 이상이므로 `f(x,y)≥0`. 그리고 `x^2+y^2=1`인 모든 점에서 `f=0`이다. 따라서 **단위원 위의 모든 점이 전역 최솟점**이다.  
> 편도함수는  
>
> **Formula:** `f_x=4x(x^2+y^2-1), f_y=4y(x^2+y^2-1),`
>
>
> **Formula:** `f_xx=4(x^2+y^2-1)+8x^2, f_yy=4(x^2+y^2-1)+8y^2, f_xy=8xy.`
>
> `(1,0)`에서 그래디언트는 0이고 헤세 행렬은  
>
> **Formula:** `∇^2f(1,0)=[8 0; 0 0].`
>
> 이는 **양의 준정부호이지만 양의 정부호는 아니다.** 따라서 2계 충분조건을 만족하지 않지만, 이미 전역 최솟점임을 확인했다.  
> 또한 `(x,y)=(1-ε,0)`, `0<ε<2`이면  
>
> **Formula:** `∇^2 f(1-ε,0)=[4(ε^2-2ε)+8(1-ε)^2 0; 0 4(ε^2-2ε)].`
>
> 오른쪽 아래 성분이 음수이므로 이 행렬은 양의 준정부호가 아니다. 즉 이 함수는 단위원 내부의 이런 점에서는 볼록함수가 아니다.  
> *주:* 원자료의 `(1,0)`에서 헤세 행렬 `(4,0;0,0)` 표기는 도함수에 비춰 볼 때 오타이며, 위처럼 `(8,0;0,0)`이 정확하다. 원자료의 근방 헤세 행렬 평가점도 `(1-ε,0)`으로 바로잡았다.
{: .prompt-info }

> **Example | A Strict Local Minimum Without a Global Minimum**  
>
> **Formula:** `f(x,y)=x^2+y^2(1-x)^3.`
>
> 정지점은  
>
> **Formula:** `f_x=2x-3y^2(1-x)^2=0, f_y=2y(1-x)^3=0`
>
> 으로 구한다. `x=1`이면 첫 번째 식에서 `2=0`이라는 모순이 생기고, `y=0`이면 `x=0`이므로 정지점은 `(0,0)` 하나다.  
> 헤세 행렬은  
>
> **Formula:** `∇^2 f(x,y)=[2+6y^2(1-x) -6y(1-x)^2; -6y(1-x)^2 2(1-x)^3].`
>
> `(0,0)`에서는  
>
> **Formula:** `∇^2 f(0,0)=[2 0; 0 2],`
>
> 따라서 **엄격한 국소 최솟점**이다. 그러나  
>
> **Formula:** `f(x,1)=x^2+(1-x)^3→-∞ (x→+∞).`
>
> 함수값이 아래로 한없이 내려갈 수 있으므로 **전역 최솟값은 존재하지 않는다.**
{: .prompt-info }

이 예제들은 **'국소 최적'과 '전역 최적'이 다르다**는 점, 그리고 2계 충분조건이 전역 최적성을 보장하지 않는다는 점을 보여 준다.

## 6. Convexity

### 6.1 Definitions of Convex and Concave Functions

> **Definition | Convex and Concave Functions**  
> 두 점 `x,y`와 `0<λ<1`에 대해  
>
> **Formula:** `f(λ x+(1-λ)y)≤λ f(x)+(1-λ)f(y)`
>
> 를 만족하면 `f`는 **볼록함수(convex function)**다.  
> 반대로  
>
> **Formula:** `f(λ x+(1-λ)y)≥λ f(x)+(1-λ)f(y)`
>
> 를 만족하면 **오목함수(concave function)**다.  
> 서로 다른 `x≠ y`에서 부등호가 항상 엄격하면 각각 **엄격한 볼록함수(strictly convex)**, **엄격한 오목함수(strictly concave)**라고 한다. 원자료의 '모든 `x,y`에서 strict'라는 표현은 `x=y`일 때 양변이 같으므로 `x≠ y`로 이해해야 한다.
{: .prompt-tip }

### 6.2 Characterizing Convexity Using the Hessian

> **Theorem | Characterization of Convexity**  
> `f:ℝ^n→ℝ`가 두 번 연속 미분 가능할 때  
> - `f`가 **볼록함수** `⇔` 모든 `x`에서 `∇^2 f(x)`가 **양의 준정부호**.
> - 모든 `x`에서 `∇^2 f(x)`가 **양의 정부호**라면 `f`는 **엄격한 볼록함수**.
> - `f`가 **오목함수** `⇔` 모든 `x`에서 `∇^2 f(x)`가 **음의 준정부호**.
> - 모든 `x`에서 `∇^2 f(x)`가 **음의 정부호**라면 `f`는 **엄격한 오목함수**.  
> 위 정의와 판정은 볼록집합 및 열린 볼록집합 위의 함수에도 적절히 확장된다.  
> **Important:** '헤세 행렬이 모든 점에서 양의 정부호'는 엄격한 볼록성을 보장하는 **충분조건**이지만, **필요조건은 아니다**. 예컨대 `f(x)=x^4`는 엄격하게 볼록하지만 `f''(0)=0`이다. 반대로 음의 정부호 조건도 엄격한 오목성의 충분조건일 뿐이다.
{: .prompt-tip }

### 6.3 Convexity and Optimality

> **Theorem | Convexity and Optimality**  
> `C⊆ℝ^n`가 볼록집합이고 `f:C→ℝ`가 볼록함수이면 다음이 성립한다.  
> **(a)** 모든 **국소 최솟점**은 **전역 최솟점**이다.  
> **(b)** `f`가 연속적으로 1계 미분 가능하고 적절한 내부 정지점 `x^*`에서 `∇ f(x^*)=0`이면 그 점은 **전역 최솟점**이다.  
> **(c)** `f`가 **엄격한 볼록함수**이면 전역 최솟점은 **최대 하나**다(존재 자체를 보장하지는 않음).  
> 오목함수의 경우에도 '최솟점'을 '최댓점'으로 바꾼 유사한 정리가 성립한다.  
> 따라서 **볼록집합에서 볼록함수를 최소화하는 문제(convex optimization problem)**는 최적화 이론에서 특히 중요하다.
{: .prompt-tip }

## 7. Coercivity and Existence of Global Minimizers

> **Definition | Coercive Function**  
> 함수 `f:ℝ^n→ℝ`가  
>
> **Formula:** `lim_(||x||→∞) f(x)=+∞`
>
> 를 만족하면 **강제적(coercive)**이라고 한다. 즉 입력 벡터가 어느 방향으로든 무한히 멀어지면 함수값도 `+∞`로 증가한다.
{: .prompt-tip }

> **Theorem | Existence of a Global Minimizer**  
> `f:ℝ^n→ℝ`가 **연속**이고 **강제적**이면 `f`는 적어도 하나의 **전역 최솟점**을 갖는다. 실무적으로는 국소 최솟점 후보를 찾아 그 함수값을 비교해 전역 최솟점을 판정하는 데 활용할 수 있다.
{: .prompt-tip }

### 7.1 Revisiting the Three-Variable Example

> **Example | A Nonconvex but Coercive Function**  
> 다시  
>
> **Formula:** `f(x,y,z)=(x^2-z)^2+(y-z)^2+2z`
>
> 를 생각하자. 앞서 `(0,-1,-1)`이 엄격한 국소 최솟점임을 확인했다. 이를 전역 최솟점으로 확정하기 위해 **볼록성** 또는 **강제성**을 확인해 보자.  
> **① Testing Convexity**  
>
> **Formula:** `∇^2f=[12x^2-4z 0 -4x; 0 2 -2; -4x -2 4].`
>
> `(x,y,z)=(0,0,1)`에서는 첫 번째 대각성분이 `-4<0`이다. 따라서 헤세 행렬이 양의 준정부호일 수 없고, **이 함수는 전역적으로 볼록하지 않다.**  
> **② Testing Coercivity: Three Equivalent Forms**  
>
> **Formula:** `f =(x^2-z)^2+(y-z)^2+2z (1); =(x^2-z)^2+(y-z-1)^2+2y-1 (2); =(x^2-z-1)^2+2x^2-1+(y-z)^2 (3)`
>
> 원자료에서는 다음 상황을 나누어 모든 무한대로 가는 경로에서 `f→+∞`임을 설명한다.  
> 1. `|x|→∞`이고 `y,z`가 유계라면 (1)에서 `f→+∞`.
> 2. `|y|→∞`이고 `z`가 유계이거나, `|z|→∞`이고 `y`가 유계라면 (1) 또는 (2)에서 `f→+∞`.
> 3. `y→+∞` 또는 `z→+∞`이면 (2) 또는 (1)에서 `f→+∞`.
> 4. `y→-∞`와 `z→-∞`가 함께 발생하면 (3)을 이용해 `f→+∞`임을 확인한다.  
> **Conclusion:** 함수는 강제적이며 전역 최솟점을 갖는다. 유일한 정지점 `(0,-1,-1)`이 바로 **전역 최솟점**이다.
{: .prompt-info }

여기서는 더 간결하게 **전역 최적값을 직접 증명**할 수도 있다. (3)에서
> **Formula:** `f(x,y,z)=(x^2-z-1)^2+(y-z)^2+2x^2-1≥-1.`
등호가 성립하려면 `x=0`, `y=z`, `x^2-z-1=0`이므로 `z=y=-1`이다. 따라서
> **Formula:** `argmin f = {(0,-1,-1)}, min f = -1`
> **Key Takeaways | Connecting the Optimality Tests**  
> **Stationarity Condition** `∇ f=0` → 최적해 *후보* 찾기.  
> **Hessian Definiteness** → 후보의 *국소* 최적성 분류.  
> **Convexity** → 국소 최적성을 *전역* 최적성으로 연결.  
> **Coercivity + Continuity** → *전역 최솟점의 존재* 보장.  
> 때로는 제곱완성과 같은 직접적인 하한 증명이 이 모든 방법보다 간단하다.
{: .prompt-tip }
 
>
> - **Optimization** — 최적화 · 목적함수를 최솟값 또는 최댓값으로 만드는 과정
> - **Objective function** — 목적함수 · 최적화 대상 함수
> - **Feasible set** — 가능집합 · 제약을 만족하는 해의 집합
> - **Constraint** — 제약조건 · 가능한 해를 제한하는 조건
> - **Critical point** — 임계점 · 도함수가 0이거나 존재하지 않는 점(맥락에 따라 정의에 차이)
> - **Stationary point** — 정지점 · 그래디언트가 0인 점
> - **Gradient** — 그래디언트, 기울기 벡터 · 1계 편도함수 벡터
> - **Hessian matrix** — 헤세 행렬 · 2계 편도함수 행렬
> - **Necessary condition** — 필요조건 · 최적점이 반드시 만족하는 조건
> - **Sufficient condition** — 충분조건 · 만족하면 최적성 결론이 따라오는 조건
> - **Positive definite** — 양의 정부호 · `x^TAx>0` (`x≠0`)
> - **Positive semidefinite** — 양의 준정부호 · `x^TAx≥0`
> - **Negative definite** — 음의 정부호 · `x^TAx<0` (`x≠0`)
> - **Negative semidefinite** — 음의 준정부호 · `x^TAx≤0`
> - **Local minimizer** — 국소 최솟점 · 근방에서 가장 작은 값의 점
> - **Global minimizer** — 전역 최솟점 · 가능집합 전체에서 가장 작은 값의 점
> - **Saddle point** — 안장점 · 극솟점·극댓점이 아닌 정지점
> - **Convex function** — 볼록함수 · 선분 위 함수값이 양 끝값의 선형보간 이하
> - **Concave function** — 오목함수 · 볼록함수의 반대 부등식을 만족
> - **Strictly convex** — 엄격한 볼록 · 서로 다른 두 점 사이에서 엄격한 볼록부등식
> - **Coercive** — 강제적, `||x||→∞`일 때 `f(x)→+∞`
> - **Leading principal minor** — 선도 주인수식 · 왼쪽 위 `k× k` 부분행렬의 행렬식
> - **Eigenvalue criterion** — 고유값 판정법 · 고유값의 부호를 통해 정부호성 판정
> - **Sylvester's criterion** — 실베스터 판정법 · 선도 주인수식 부호로 엄격한 정부호성 판정
> - **Lagrange multipliers** — 라그랑주 승수법 · 등식 제약 최적화 방법
> - **KKT conditions** — KKT 조건 · 제약 최적화의 주요 최적성 조건
> - **Linear / Integer programming** — 선형 / 정수계획법 · 선형 목적·제약 또는 정수 변수 최적화
{: .prompt-tip }

---

Piotr Więcek, *Optimization Theory*, Wrocław University of Science and Technology
