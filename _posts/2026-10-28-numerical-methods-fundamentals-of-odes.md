---
title: "Numerical Methods : Fundamentals of Ordinary Differential Equations"
date: 2026-10-27 18:00:00 +0100
categories:
  - Numerical Analysis
tags:
  - Ordinary-Differential-Equations
  - Initial-Value-Problems
  - Lipschitz-Continuity
  - Existence-and-Uniqueness
math: false
---

## Introduction

ODE의 해를 수치적으로 근사하기 전에 **해가 존재하는지, 유일한지, 초기조건의 변화에 안정적으로 반응하는지**부터 살펴본다. 강의 `Lec-1_ODE2026.pdf`의 Section 을 순서대로 정리하고, Peano 현상에 관한 Lab Handout Exercise 2를 해당 이론 뒤에 배치했다.

## 1. Framework and Problem Formulation

ODE는 알려진 초기 상태에서 미지 함수의 변화를 추적하는 방정식이다. n차원에서는 `x(t) ∈ ℝⁿ`인 함수를 구하며, 강의에서의 기본 조건은 연속 벡터장 `f: Ω → ℝⁿ`이다.

![Initial value problem](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7Dx%27%28t%29%3Df%28t%2Cx%28t%29%29%2C%5Cqquad%20x%28t_0%29%3Dx_0)

여기서 `Ω`는 `(t, x)`가 속하는 열린 영역이며, 해 `x`는 초기시간을 포함하는 구간에서 연속 미분 가능해야 한다. **미분방정식과 초기조건을 동시에 만족해야 초기값 문제의 해**라고 할 수 있다.

### 1.1 Integral Operator Reformulation

미분방정식의 양변을 초기시간 `t₀`부터 `t`까지 적분하면 다음과 같다.

![Volterra integral equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7Dx%28t%29%3Dx_0%2B%5Cint_%7Bt_0%7D%5E%7Bt%7Df%28s%2Cx%28s%29%29%5C%2Cds)

이는 **Volterra integral equation**이라고 한다. 미분형과 적분형은 적절한 연속성 조건 아래 동치이다. 적분형으로 바꾸면 함수를 입력받아 다시 함수를 반환하는 연산자를 정의할 수 있어, 뒤에서 고정점 정리로 존재성과 유일성을 증명할 때 활용한다.

![Picard operator](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%28%5Cmathcal%7BT%7Dx%29%28t%29%3Dx_0%2B%5Cint_%7Bt_0%7D%5E%7Bt%7Df%28s%2Cx%28s%29%29%5C%2Cds)

해를 구한다는 것은 결국 `𝒯x = x`를 만족하는 **고정점(Fixed Point)**을 찾는 문제로 해석할 수 있다.

## 2. Local Existence and Uniqueness

미분방정식에 해가 존재한다는 것과 그 해가 하나로 결정된다는 것은 다른 주장이다. 강의에서는 이를 구분하기 위해 Lipschitz 조건을 도입한다.

### 2.1 Lipschitz Continuity

함수 `f(t, x)`가 상태변수 `x`에 대해 Lipschitz라는 것은, 두 상태의 차이가 함수값에서 지나치게 증폭되지 않는다는 의미이다.

![Lipschitz condition](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5CVert%20f%28t%2Cx_1%29-f%28t%2Cx_2%29%5CVert%5Cle%20L%5CVert%20x_1-x_2%5CVert)

- **Locally Lipschitz:** 각 점 근처에서 어떤 상수 `L`을 선택할 수 있다.
- **Uniformly Lipschitz:** 지정한 영역에서 시간 `t`와 무관한 동일한 `L`을 사용할 수 있다.

`f`가 미분 가능하고 상태변수에 관한 도함수가 영역에서 유계이면 평균값정리를 사용해 Lipschitz 조건을 확인할 수 있다. 하지만 단순히 함수가 연속이라는 사실만으로는 이 조건을 보장할 수 없다.

### 2.2 Picard–Lindelöf Theorem

**Picard–Lindelöf 정리**는 `f`가 연속이고 `x`에 대해 국소 Lipschitz일 때 초기점 주변에서 유일한 국소해가 존재한다고 말한다.

강의의 증명 아이디어는 적분 연산자 `𝒯`를 적절한 함수공간에 정의하고 **Banach Fixed Point Theorem**을 적용하는 것이다. 적절한 노름과 구간의 길이를 선택해 `𝒯`가 수축사상(Contraction)이 됨을 보이면, 반복 적용으로 유일한 고정점에 도달한다.

![Contraction estimate](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5CVert%5Cmathcal%7BT%7Dx-%5Cmathcal%7BT%7Dy%5CVert_%7B%5Clambda%7D%5Cle%5Cfrac%7BL%7D%7B%5Clambda%7D%5CVert%20x-y%5CVert_%7B%5Clambda%7D%2C%5Cqquad%5Clambda%3EL)

위 식에서 `L/λ < 1`이므로 수축 조건이 성립한다. **유일성은 단순한 연속성보다 강한 조건으로부터 나온다**는 점이 핵심이다.

### 2.3 Lecture Example — A Unique Solution

강의에서 제시한 문제는 다음과 같다.

![Unique solution example](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7Dy%27%3D1%2Bt%5Csin%28ty%29%2C%5Cqquad%20y%280%29%3D0%2C%5Cqquad%200%5Cle%20t%5Cle2)

우변을 `y`로 편미분하면 `t² cos(ty)`이며, `0 ≤ t ≤ 2`에서 절댓값이 `4` 이하이다. 따라서 상태변수 `y`에 관해 Lipschitz 조건이 성립하여 해의 국소 존재성과 유일성을 확인할 수 있다. 이는 유일성이 바로 수치 오차를 줄여 준다는 뜻은 아니며, **수치해가 근사하려는 대상이 하나로 정해진다**는 수학적 출발점이다.

## 3. Existence Without Uniqueness: Peano’s Theorem

**Peano의 존재 정리**는 `f`가 연속이면 적어도 하나의 국소해가 존재함을 보장한다. 그러나 Lipschitz 조건이 없으면 여러 해가 같은 초기조건을 만족할 수 있다.

강의에서는 균일 유계성과 동등연속성을 확보한 근사해의 부분수열을 **Arzelà–Ascoli 정리**로 수렴시키는 증명 아이디어를 소개한다. 이는 수축사상을 통해 유일한 고정점을 얻는 Picard–Lindelöf 접근과 다르다.

### 3.1 Peano Phenomenon — Lecture Example / Lab Handout Exercise 2

![Nonunique initial value problem](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7Dx%27%3D3x%5E%7B2%2F3%7D%2C%5Cqquad%20x%280%29%3D0)

우변은 `x=0`에서 연속이지만 미분계수가 발산하므로 그 점 근처에서 Lipschitz가 아니다. 실제로 임의의 `c ≥ 0`에 대해 다음 함수가 해가 된다.

![Family of Peano solutions](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7Dx_c%28t%29%3D%5Cbegin%7Bcases%7D0%2C%26t%3Cc%2C%5C%5C%28t-c%29%5E3%2C%26t%5Cge%20c.%5Cend%7Bcases%7D)

`c`를 바꾸면 **0에 머무르다가 서로 다른 시간에 증가하기 시작하는 해**를 무한히 만들 수 있다.

**What happens with Euler’s method?**

Forward Euler를 `x₀=0`에서 시작하면 `x₁=0+h·3·0^(2/3)=0`이다. 이후에도 같은 계산이 반복되어 모든 `xₙ=0`이 된다. 즉, Euler 방법은 가능한 수많은 해 중에서 **항상 0인 해를 선택**한다. 수치 알고리즘이 하나의 결과를 출력했다고 해서 연속 문제의 해가 유일하다는 뜻은 아니다.

## 4. Well-Posed Problems

ODE의 수치해를 신뢰하려면 먼저 원래 문제가 **Well-posed**인지 살펴봐야 한다. 강의에서의 조건은 다음과 같다.

1. 해가 존재하고 유일해야 한다.
2. 초기값이나 방정식의 우변에 작은 변화가 생겼을 때 해도 작게 변해야 한다.

원래 문제와 교란된 문제를 다음과 같이 비교한다.

![Original IVP](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7Dy%27%3Df%28t%2Cy%29%2C%5Cquad%20y%28a%29%3D%5Calpha)

![Perturbed IVP](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7Dz%27%3Df%28t%2Cz%29%2B%5Cdelta%28t%29%2C%5Cquad%20z%28a%29%3D%5Calpha%2B%5Cdelta_0)

`|δ(t)| < ε`이고 `|δ₀| < ε`일 때, 해의 차이가 구간 전체에서 `|z(t)-y(t)| < kε`로 제어된다면 작은 입력 오차가 작은 출력 오차로 이어진다.

이 개념이 수치해석에서 중요한 이유는 **반올림 오차와 근사 오차를 완전히 피할 수 없기 때문**이다. 문제 자체가 작은 교란에 민감하면 계산 결과가 불안정해지기 쉽다.

### 4.1 Lecture Example — Well-Posed IVP

![Well posed IVP](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7Dy%27%3Dy-t%5E2%2B1%2C%5Cqquad%20y%280%29%3D0.5%2C%5Cqquad%200%5Cle%20t%5Cle2)

`f(t,y)=y−t²+1`이므로 `∂f/∂y=1`이다. 따라서 `y`에 대해 전역 Lipschitz 조건이 성립한다. 이 선형 방정식은 주어진 유한 구간에서 유일한 해를 가지며, 뒤의 Grönwall 부등식으로 교란에 따른 해의 차이도 제어할 수 있다.

## 5. Maximal Solutions and Continuation

국소적으로 해가 존재한다고 해서 그 해를 모든 시간으로 연장할 수 있는 것은 아니다. 더 이상 확장할 수 없는 해가 정의되는 구간을 **Maximal Interval**이라고 한다.

강의의 **Blow-up Alternative**는 연속성과 국소 Lipschitz 조건 아래, 유한 시간에 최대해가 끝난다면 해가 정의역의 모든 콤팩트 부분집합에서 벗어나야 한다는 내용을 담고 있다. 정의역이 전체 `ℝ × ℝⁿ`이라면 유한시간 종료는 해의 크기 발산과 연결된다.

### 5.1 Lecture Example — Finite-Time Blow-Up

![Blow up IVP](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7Dx%27%3D1%2Bx%5E2%2C%5Cqquad%20x%280%29%3D0)

![Tangent solution blow up](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7Dx%28t%29%3D%5Ctan%28t%29%2C%5Cqquad%20t%5Cuparrow%5Cfrac%7B%5Cpi%7D%7B2%7D%5C%20%5CLongrightarrow%5C%20x%28t%29%5Cto%2B%5Cinfty)

정확한 해가 `π/2`에 가까워질수록 급격히 커진다. 고정된 시간 간격으로 계산하는 Euler 스킴에서는 이 특이 시점을 정확하게 표현하기 어려우므로 **계산이 끝까지 실행됐다는 사실만으로 해가 전 시간에 존재한다고 판단할 수 없다.**

## 6. Continuous and Smooth Dependence on Parameters

초기값이나 모델 파라미터를 조금 바꾸면 해도 조금 바뀌는지 확인하는 것은 모델의 신뢰성을 이해하는 데 중요하다.

### 6.1 Grönwall–Bellman Inequality

Grönwall 부등식은 적분형 부등식을 만족하는 오차가 얼마나 증가할 수 있는지 상한을 제공한다.

![Gronwall hypothesis](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7Du%28t%29%5Cle%20C%2B%5Cint_%7Bt_0%7D%5E%7Bt%7D%5Cbeta%28s%29u%28s%29%5C%2Cds)

![Gronwall conclusion](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7Du%28t%29%5Cle%20C%5Cexp%5Cleft%28%5Cint_%7Bt_0%7D%5E%7Bt%7D%5Cbeta%28s%29%5C%2Cds%5Cright%29)

즉, 초기 오차 `C`가 시간에 따라 확대될 수 있지만, 그 증가를 지수함수 형태의 상한으로 제어할 수 있다. 이는 **초기조건에 대한 연속 의존성과 수치 오차의 누적 분석**에 연결된다.

## 7. Linear Systems of Ordinary Differential Equations

강의에서는 비자율 선형 시스템을 다음과 같이 정의한다.

![Linear ODE system](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7Dx%27%28t%29%3DA%28t%29x%28t%29%2Bg%28t%29)

### 7.1 Fundamental Solutions and Matrix Exponentials

동차방정식 `x'=A(t)x`의 선형독립인 해들을 열벡터로 모으면 **Fundamental Matrix** `Ψ(t)`를 얻는다.

![Fundamental matrix](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5CPsi%27%28t%29%3DA%28t%29%5CPsi%28t%29%2C%5Cqquad%5Cdet%5CPsi%28t%29%5Cneq0)

해 공간의 기저를 이룬다는 점에서 역행렬이 존재해야 한다. 강의에 제시된 Liouville–Abel–Jacobi 공식은 다음과 같다.

![Liouville formula](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdet%5CPsi%28t%29%3D%5Cdet%5CPsi%28t_0%29%5Cexp%5Cleft%28%5Cint_%7Bt_0%7D%5E%7Bt%7D%5Coperatorname%7Btr%7D%28A%28s%29%29%5C%2Cds%5Cright%29)

### 7.2 Variation of Constants — Duhamel’s Principle

비동차항 `g(t)`가 있으면 기본행렬을 이용해 해를 표현한다.

![Duhamel formula](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7Dx%28t%29%3D%5CPsi%28t%29%5CPsi%28t_0%29%5E%7B-1%7Dx_0%2B%5Cint_%7Bt_0%7D%5E%7Bt%7D%5CPsi%28t%29%5CPsi%28s%29%5E%7B-1%7Dg%28s%29%5C%2Cds)

행렬 `A`가 상수이면 행렬지수함수 `e^(tA)`를 사용할 수 있다. 첫 번째 항은 초기조건의 영향, 적분항은 외부 입력 `g`가 누적된 효과를 나타낸다.



## Detailed Lecture Notes — Definitions, Proofs, and Calculations

아래 내용은 앞에서 정리한 Lecture 1의 각 절을 **정의의 정확한 조건, 증명 과정, 예제 계산**을 중심으로 다시 전개한 것이다. 앞부분은 개념의 흐름을 빠르게 확인하기 위한 것이고, 이 절에서는 강의자료에 있는 수학적 논거를 순서대로 확인한다.

### 1.2 Precise Definition of an Initial Value Problem

강의에서는 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5COmega%5Csubset%5Cmathbb%20R%5Ctimes%5Cmathbb%20R%5En)을 열린 집합, ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f%3A%5COmega%5Cto%5Cmathbb%20R%5En)을 연속 벡터장으로 둔다. 초기값 문제의 해는 단순히 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x%27%3Df)를 형식적으로 만족하는 식이 아니라, 다음 세 조건을 동시에 충족하는 **연속 미분 가능한 함수** ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x%3AJ%5Cto%5Cmathbb%20R%5En)이다.

1. 초기시각 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20t_0)을 포함하는 열린 구간 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20J) 위에서 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%28t%2Cx%28t%29%29%5Cin%5COmega)이다.
2. 모든 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20t%5Cin%20J)에 대해 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x%27%28t%29%3Df%28t%2Cx%28t%29%29)이다.
3. 초기조건 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x%28t_0%29%3Dx_0)을 만족한다.

이때 적분형 표현의 동치는 미적분학의 기본정리에서 얻는다. 해가 존재한다면 양변 적분으로

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cint_%7Bt_0%7D%5E%7Bt%7Dx%27%28s%29ds%3D%5Cint_%7Bt_0%7D%5E%7Bt%7Df%28s%2Cx%28s%29%29ds)

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x%28t%29-x_0%3D%5Cint_%7Bt_0%7D%5E%7Bt%7Df%28s%2Cx%28s%29%29ds)

를 얻는다. 반대로 오른쪽 적분함수가 연속이므로 적분식을 미분하면 원래 미분방정식을 복원한다. 이 **양방향 논증**이 Volterra integral equation과 IVP가 동치인 이유다.

### 2.4 Deriving the Lipschitz Estimate

강의의 Theorem 2.3에서는 스칼라 함수에 대해 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%7C%5Cpartial%20f%2F%5Cpartial%20y%7C%5Cle%20L)인 경우 평균값정리를 적용한다. 고정된 시각 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20t)에서 두 점 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y_1%2Cy_2) 사이의 어떤 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cxi)에 대해

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f%28t%2Cy_1%29-f%28t%2Cy_2%29%3D%5Cfrac%7B%5Cpartial%20f%7D%7B%5Cpartial%20y%7D%28t%2C%5Cxi%29%28y_1-y_2%29)

이므로 절댓값을 취하면

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%7Cf%28t%2Cy_1%29-f%28t%2Cy_2%29%7C%5Cle%20L%7Cy_1-y_2%7C)

이다. 벡터값 문제에는 편미분 대신 상태변수에 관한 Jacobian의 연산자 노름을 사용한다. **연속성만으로는 Lipschitz가 따라오지 않는다.** 예를 들어 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f%28x%29%3D3x%5E%7B2%2F3%7D)은 0에서 연속이지만 미분의 크기가 발산한다.

### 2.5 Full Statement of the Picard–Lindelöf Theorem

초기점 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%28t_0%2Cx_0%29%5Cin%5COmega) 주변에 실린더

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20Q%3D%5Bt_0-a%2Ct_0%2Ba%5D%5Ctimes%5Coverline%20B_b%28x_0%29%5Csubset%5COmega)

가 존재한다고 하자. ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20M%3D%5Cmax_%7B%28t%2Cx%29%5Cin%20Q%7D%5C%7Cf%28t%2Cx%29%5C%7C)을 정의하고, ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20M%3E0)이면

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cdelta%3D%5Cmin%5Cleft%28a%2C%5Cfrac%20bM%5Cright%29)

을 선택한다. ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20M%3D0)인 경우에는 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cdelta%3Da)로 택할 수 있다. 이 범위에서 국소 Lipschitz 조건이 성립하면 초기값 문제의 **유일한 국소해**가 존재한다. 강의의 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20b%2FM) 조건은 해가 실린더의 공간 영역을 벗어나지 않도록 하기 위한 길이 제한이다.

### 2.6 Proof via Banach's Fixed Point Theorem

**Step 1 — Construct the function space.**

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20X%3D%5C%7Bx%5Cin%20C%28%5Bt_0-%5Cdelta%2Ct_0%2B%5Cdelta%5D%2C%5Cmathbb%20R%5En%29%3Ax%28t_0%29%3Dx_0%2C%5C%20%5Csup_t%5C%7Cx%28t%29-x_0%5C%7C%5Cle%20b%5C%7D)

이 공간에 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5C%7Cx%5C%7C_%5Cinfty%3D%5Csup_t%5C%7Cx%28t%29%5C%7C)를 사용한다. ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20X)는 연속함수 공간의 닫힌 부분집합이므로 완비다.

**Step 2 — Define the Picard operator.**

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%28%7B%5Ccal%20T%7Dx%29%28t%29%3Dx_0%2B%5Cint_%7Bt_0%7D%5E%7Bt%7Df%28s%2Cx%28s%29%29ds)

**Step 3 — Show that the operator maps X into itself.** 적분의 방향에 주의하고 적분구간을 작은 시간에서 큰 시간으로 잡으면

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5C%7C%28%7B%5Ccal%20T%7Dx%29%28t%29-x_0%5C%7C%5Cle%5Cleft%7C%5Cint_%7Bt_0%7D%5E%7Bt%7D%5C%7Cf%28s%2Cx%28s%29%29%5C%7Cds%5Cright%7C%5Cle%20M%7Ct-t_0%7C%5Cle%20M%5Cdelta%5Cle%20b)

이다. 또한 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%28%7B%5Ccal%20T%7Dx%29%28t_0%29%3Dx_0)이므로 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%7B%5Ccal%20T%7D%28X%29%5Csubset%20X)이다.

**Step 4 — Prove contraction.** 강의는 가중 Bielecki 노름

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5C%7Cx%5C%7C_%5Clambda%3D%5Csup_%7Bt%7De%5E%7B-%5Clambda%7Ct-t_0%7C%7D%5C%7Cx%28t%29%5C%7C)

을 사용한다. Lipschitz 조건으로

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5C%7C%28%7B%5Ccal%20T%7Dx%29%28t%29-%28%7B%5Ccal%20T%7Dy%29%28t%29%5C%7C%5Cle%20L%5Cint_%7B%5Cmin%28t%2Ct_0%29%7D%5E%7B%5Cmax%28t%2Ct_0%29%7D%5C%7Cx%28s%29-y%28s%29%5C%7Cds)

이고, 정의에 의해 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5C%7Cx%28s%29-y%28s%29%5C%7C%5Cle%20e%5E%7B%5Clambda%7Cs-t_0%7C%7D%5C%7Cx-y%5C%7C_%5Clambda)이므로

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20e%5E%7B-%5Clambda%7Ct-t_0%7C%7D%5C%7C%28%7B%5Ccal%20T%7Dx%29%28t%29-%28%7B%5Ccal%20T%7Dy%29%28t%29%5C%7C%5Cle%5Cfrac%20L%5Clambda%281-e%5E%7B-%5Clambda%7Ct-t_0%7C%7D%29%5C%7Cx-y%5C%7C_%5Clambda)

결국

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5C%7C%7B%5Ccal%20T%7Dx-%7B%5Ccal%20T%7Dy%5C%7C_%5Clambda%5Cle%5Cfrac%20L%5Clambda%5C%7Cx-y%5C%7C_%5Clambda)

이다. ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Clambda%3EL)로 선택하면 수축상수가 1보다 작아 Banach 고정점 정리를 적용할 수 있다. 따라서 유일한 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x%5E%2A%3D%7B%5Ccal%20T%7Dx%5E%2A)가 존재하고, 이 함수가 IVP의 유일한 해다. 이것이 Lecture 1에서 말하는 **존재성과 유일성의 동시에 성립하는 이유**다.

### 2.7 Worked Lipschitz Example

강의 예제는 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y%27%3D1%2Bt%5Csin%28ty%29%2C%5C%20y%280%29%3D0%2C%5C%200%5Cle%20t%5Cle2)이다.

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f%28t%2Cy%29%3D1%2Bt%5Csin%28ty%29)

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cfrac%7B%5Cpartial%20f%7D%7B%5Cpartial%20y%7D%3Dt%5E2%5Ccos%28ty%29)

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cleft%7C%5Cfrac%7B%5Cpartial%20f%7D%7B%5Cpartial%20y%7D%5Cright%7C%5Cle%20t%5E2%5Cle%204)

이므로 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20L%3D4)를 선택할 수 있다. 따라서 국소 유일성이 확보되며, 주어진 구간에서 우변의 크기가 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%7Cf%7C%5Cle%201%2B%7Ct%7C%5Cle3)으로 제한되므로 해의 유한시간 폭주 문제도 발생하지 않는다. 강의의 계산은 편미분 상계 하나로 Lipschitz 조건을 확인하는 전형적 예제다.

### 3.2 Peano's Theorem and the Arzelà–Ascoli Proof Strategy

**Peano 정리**는 연속인 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f)에 대해 적어도 하나의 국소해가 존재한다고 주장하지만, 유일성은 주장하지 않는다. 강의의 증명 전략은 다음과 같다.

1. Euler polygonal approximation 또는 연속적 regularization을 이용해 근사해 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x_%5Cvarepsilon)의 족을 만든다.
2. ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5C%7Cf%5C%7C%5Cle%20M)인 작은 영역에서 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5C%7Cx_%5Cvarepsilon%27%5C%7C%5Cle%20M)을 확보한다.
3. 따라서 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5C%7Cx_%5Cvarepsilon%28t%29-x_%5Cvarepsilon%28s%29%5C%7C%5Cle%20M%7Ct-s%7C)로 동등연속성(equicontinuity)을 얻는다.
4. 공간적 유계성까지 확보하면 **Arzelà–Ascoli 정리**에 따라 균등수렴하는 부분수열이 존재한다.
5. 근사식의 적분형 표현에서 극한을 취해 Volterra 방정식을 얻는다.

여기까지가 강의에 제시된 증명 전략이다. 강의의 질문에 있는 **Schauder Fixed Point Theorem**은 연속·콤팩트 연산자를 적절한 닫힌 볼록집합에 적용하는 대안적 접근으로 설명할 수 있다. 그러나 강의가 Peano 증명을 Schauder로 완성해 제시한 것은 아니므로, 이는 **Additional Explanation**으로 구분한다.

#### Non-Uniqueness Calculation

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x_c%28t%29%3D0%5C%20%28t%3Cc%29%2C%5Cquad%20x_c%28t%29%3D%28t-c%29%5E3%5C%20%28t%5Cge%20c%29)라고 하자. ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20t%3Cc)에서 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x_c%27%3D0%3D3x_c%5E%7B2%2F3%7D)이다. ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20t%3Ec)에서는

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x_c%27%28t%29%3D3%28t-c%29%5E2)

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%203x_c%28t%29%5E%7B2%2F3%7D%3D3%28%28t-c%29%5E3%29%5E%7B2%2F3%7D%3D3%28t-c%29%5E2)

이다. 접합점 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20t%3Dc)에서도 좌우 미분계수가 0으로 일치한다. 따라서 모든 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20c%5Cge0)이 동일한 초기조건을 만족하는 해를 준다.

**Euler와 연결:** ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x_0%3D0)에서 Forward Euler는 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x_%7Bn%2B1%7D%3Dx_n%2B3h%20x_n%5E%7B2%2F3%7D)이므로 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x_1%3D0), 이후에도 귀납적으로 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x_n%3D0)이다. 즉, 이 알고리즘은 무한히 많은 가능한 해 중 **항상 0에 머무르는 해**를 선택한다. 이것이 Lab Handout Exercise 2의 핵심이다.

### 4.2 Exact Definition of Well-Posedness

강의는 다음 세 요소를 강조한다: **존재**, **유일성**, **초기조건과 우변의 작은 섭동에 대한 연속적 안정성**. 원래 문제가

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y%27%3Df%28t%2Cy%29%2C%5Cqquad%20y%28a%29%3D%5Calpha)

이면 섭동 문제는

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20z%27%3Df%28t%2Cz%29%2B%5Cdelta%28t%29%2C%5Cqquad%20z%28a%29%3D%5Calpha%2B%5Cdelta_0)

이다. ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%7C%5Cdelta%28t%29%7C%3C%5Cvarepsilon) 및 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%7C%5Cdelta_0%7C%3C%5Cvarepsilon)에서 해가 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%7Cz%28t%29-y%28t%29%7C%3Ck%5Cvarepsilon)를 만족해야 well-posed다. 이는 수치해석에서 반올림 오차와 이산화 오차가 작은 교란으로 작용하기 때문에 중요하다.

#### Lecture Example: Checking Well-Posedness

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y%27%3Dy-t%5E2%2B1%2C%5Cqquad%20y%280%29%3D0.5%2C%5Cquad%200%5Cle%20t%5Cle2)

에서 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cpartial%20f%2F%5Cpartial%20y%3D1)이므로 전역 Lipschitz 상수 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20L%3D1)을 사용할 수 있다. 섭동 해와 원래 해의 차이 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20w%3Dz-y)는

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20w%27%3Dw%2B%5Cdelta%28t%29%2C%5Cqquad%20w%280%29%3D%5Cdelta_0)

를 만족한다. 적분인자 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20e%5E%7B-t%7D)를 사용하면

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%28e%5E%7B-t%7Dw%28t%29%29%27%3De%5E%7B-t%7D%5Cdelta%28t%29)

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20w%28t%29%3De%5Et%5Cdelta_0%2B%5Cint_0%5Et%20e%5E%7Bt-s%7D%5Cdelta%28s%29ds)

이므로

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%7Cw%28t%29%7C%5Cle%20e%5Et%5Cvarepsilon%2B%28e%5Et-1%29%5Cvarepsilon%3D%282e%5Et-1%29%5Cvarepsilon%5Cle%282e%5E2-1%29%5Cvarepsilon)

이다. 따라서 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20k%3D2e%5E2)처럼 더 큰 상수를 택하면 강의의 엄격 부등식 조건을 확보할 수 있다. 섭동이 작을수록 해의 차이도 선형적으로 작아지는 것을 계산으로 확인했다.

### 5.2 Maximal Interval and Blow-Up Alternative

해 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%28J_2%2Cx_2%29)가 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%28J_1%2Cx_1%29)의 extension이라는 것은 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20J_1%5Csubsetneq%20J_2)이며 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x_2%7C_%7BJ_1%7D%3Dx_1)임을 뜻한다. **Maximal solution**은 더 넓은 구간으로 연장할 수 없는 해다.

강의의 Blow-up Alternative는 유한한 최대 존재시간에 접근할 때 해의 그래프가 정의역의 모든 콤팩트 부분집합을 결국 벗어난다는 내용이다. 특별히 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5COmega%3D%5Cmathbb%20R%5Ctimes%5Cmathbb%20R%5En)이고 오른쪽 끝점 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20t_%2B%3C%5Cinfty)이면

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Clim_%7Bt%5Cuparrow%20t_%2B%7D%5C%7Cx%28t%29%5C%7C%3D%5Cinfty)

가 된다. **Additional Calculation:** ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x%27%3D1%2Bx%5E2%2C%5C%20x%280%29%3D0)에서는 변수분리로 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Carctan%28x%29%3Dt)이므로 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x%28t%29%3D%5Ctan%28t%29)이며 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20t_%2B%3D%5Cpi%2F2)이다. 수치 격자에서 발산 직전까지는 계산이 가능해도 실제 해는 유한시간에 폭주한다는 점과 연결된다.

### 6.2 Grönwall–Bellman Inequality: General and Constant Cases

강의의 일반형은 비음수 연속함수 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20u%2C%5Calpha%2C%5Cbeta)에 대해

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20u%28t%29%5Cle%5Calpha%28t%29%2B%5Cint_%7Bt_0%7D%5E%7Bt%7D%5Cbeta%28s%29u%28s%29ds)

이면

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20u%28t%29%5Cle%5Calpha%28t%29%2B%5Cint_%7Bt_0%7D%5E%7Bt%7D%5Calpha%28s%29%5Cbeta%28s%29%5Cexp%5Cleft%28%5Cint_s%5Et%5Cbeta%28r%29dr%5Cright%29ds)

라는 것이다. 특히 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Calpha%28t%29%3DC)인 경우

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20u%28t%29%5Cle%20C%5Cexp%5Cleft%28%5Cint_%7Bt_0%7D%5E%7Bt%7D%5Cbeta%28s%29ds%5Cright%29)

로 단순화된다.

**Derivation for constant forcing:** ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20v%28t%29%3DC%2B%5Cint_%7Bt_0%7D%5Et%5Cbeta%28s%29u%28s%29ds)로 정의하면 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20u%5Cle%20v)이고 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20v%27%28t%29%3D%5Cbeta%28t%29u%28t%29%5Cle%5Cbeta%28t%29v%28t%29)이다. ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20v%3E0)에서

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cfrac%7Bd%7D%7Bdt%7D%5Clog%20v%28t%29%3D%5Cfrac%7Bv%27%28t%29%7D%7Bv%28t%29%7D%5Cle%5Cbeta%28t%29)

이므로 적분하여 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Clog%20v%28t%29-%5Clog%20C%5Cle%5Cint_%7Bt_0%7D%5E%7Bt%7D%5Cbeta%28s%29ds)를 얻는다. 지수함수를 취하면 원래 결론이 나온다. 일반형에서는 적분인자를 적용해 위에 제시한 적분형 상계를 얻는다.

**Dependence on initial data:** 같은 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f)의 두 해 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x%2Cy)에 Lipschitz 상수 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20L)이 있으면

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5C%7Cx%28t%29-y%28t%29%5C%7C%5Cle%5C%7Cx%28t_0%29-y%28t_0%29%5C%7C%2BL%5Cint_%7Bt_0%7D%5Et%5C%7Cx%28s%29-y%28s%29%5C%7Cds)

이므로 Grönwall로

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5C%7Cx%28t%29-y%28t%29%5C%7C%5Cle%20e%5E%7BL%28t-t_0%29%7D%5C%7Cx%28t_0%29-y%28t_0%29%5C%7C)

이다. 여기서 **지수 상계는 연속적 의존성을 증명하지만, 작은 오차가 시간에 따라 전혀 증가하지 않는다는 뜻은 아니다.**

### 7.3 Fundamental Matrices and the Liouville Formula

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x%27%3DA%28t%29x)에서 fundamental matrix ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5CPsi%28t%29)는 열벡터가 해공간의 기저를 이루고, 모든 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20t)에서 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cdet%5CPsi%28t%29%5Cne0)인 행렬해다.

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5CPsi%27%28t%29%3DA%28t%29%5CPsi%28t%29)

강의의 Liouville–Abel–Jacobi 공식은

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cdet%5CPsi%28t%29%3D%5Cdet%5CPsi%28t_0%29%5Cexp%5Cleft%28%5Cint_%7Bt_0%7D%5E%7Bt%7D%5Cmathrm%7Btr%7D%28A%28s%29%29ds%5Cright%29)

이다. 이를 이해하기 위해 행렬식의 미분 공식(Jacobi's formula)을 적용하면

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cfrac%20d%7Bdt%7D%5Cdet%5CPsi%28t%29%3D%5Cdet%5CPsi%28t%29%5C%2C%5Cmathrm%7Btr%7D%28%5CPsi%28t%29%5E%7B-1%7D%5CPsi%27%28t%29%29)

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%3D%5Cdet%5CPsi%28t%29%5C%2C%5Cmathrm%7Btr%7D%28%5CPsi%5E%7B-1%7DA%5CPsi%29%3D%5Cdet%5CPsi%28t%29%5C%2C%5Cmathrm%7Btr%7D%28A%28t%29%29)

가 된다. 따라서 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20D%27%3D%5Cmathrm%7Btr%7D%28A%29D)인 스칼라 방정식을 적분하면 공식이 증명된다. 마지막 등호에는 행렬 trace의 순환성질을 사용했다.

### 7.4 Variation of Constants: Derivation of Duhamel's Formula

비동차 선형 시스템 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x%27%3DA%28t%29x%2Bg%28t%29)에 대해 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x%28t%29%3D%5CPsi%28t%29c%28t%29)를 가정한다.

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x%27%28t%29%3D%5CPsi%27%28t%29c%28t%29%2B%5CPsi%28t%29c%27%28t%29)

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%3DA%28t%29%5CPsi%28t%29c%28t%29%2B%5CPsi%28t%29c%27%28t%29)

이를 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20A%28t%29x%28t%29%2Bg%28t%29)와 비교하면

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5CPsi%28t%29c%27%28t%29%3Dg%28t%29)

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20c%27%28t%29%3D%5CPsi%28t%29%5E%7B-1%7Dg%28t%29)

이고, 적분하여

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20c%28t%29%3D%5CPsi%28t_0%29%5E%7B-1%7Dx_0%2B%5Cint_%7Bt_0%7D%5E%7Bt%7D%5CPsi%28s%29%5E%7B-1%7Dg%28s%29ds)

따라서

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x%28t%29%3D%5CPsi%28t%29%5CPsi%28t_0%29%5E%7B-1%7Dx_0%2B%5Cint_%7Bt_0%7D%5E%7Bt%7D%5CPsi%28t%29%5CPsi%28s%29%5E%7B-1%7Dg%28s%29ds)

를 얻는다. ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20A%28t%29%3DA)인 자율계에서는 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5CPsi%28t%29%5CPsi%28s%29%5E%7B-1%7D%3De%5E%7B%28t-s%29A%7D)이므로

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x%28t%29%3De%5E%7B%28t-t_0%29A%7Dx_0%2B%5Cint_%7Bt_0%7D%5E%7Bt%7De%5E%7B%28t-s%29A%7Dg%28s%29ds)

가 된다. 이것이 Lecture 1의 Duhamel's Principle을 계산으로 유도한 과정이다.

### 8. Numerical Solutions: Connection to the Next Post

Lecture 1의 마지막 절은 해석해가 초등함수 형태로 표현되지 않거나 계산하기 어렵기 때문에 수치해 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y_n)이 필요하다는 내용이다. 일정 간격 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20h)의 격자에서 접선 기울기 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f%28t_n%2Cy_n%29)로 다음 점을 예측하면

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y_%7Bn%2B1%7D%3Dy_n%2Bh%20f%28t_n%2Cy_n%29)

인 Euler 방법을 얻는다. 강의의 ![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y%27%3Dt-y%2C%5C%20y%280%29%3D1%2C%5C%20h%3D0.2) 예제에서는

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y_1%3D1%2B0.2%280-1%29%3D0.8)

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y_2%3D0.8%2B0.2%280.2-0.8%29%3D0.68)

이다. Euler가 구현하기 쉽지만 전역 1차 정확도와 안정성 제약을 갖는다는 논의는 **Post 2에서 전체 증명과 실험으로 연결**한다.

## Source-by-Source Verification — Lecture 1 Remaining Details

### Existence Proof: Self-Mapping and Contraction

강의 p.2–3의 Picard 증명은 단순히 수축상수가 작다고 말하는 것이 아니라, **사상 자체가 허용된 함수공간 안에 남는지** 먼저 확인한다. 반경 `b`인 닫힌 공과 시간 폭 `a`를 택하여 `Q` 안에서 우변의 최대 노름을 `M`이라 하자. `M=0`일 때는 시간 폭에 대한 `b/M` 조건이 필요하지 않고, 그렇지 않다면 다음과 같이 고른다.

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cdelta%3D%5Cmin%5C%7Ba%2Cb%2FM%5C%7D%2C%5Cqquad%20%5CVert%20%28%5Cmathcal%20T%20x%29%28t%29-x_0%5CVert%20%5Cle%20M%7Ct-t_0%7C%5Cle%20M%5Cdelta%5Cle%20b)

그러면 `𝒯(X) ⊆ X`가 성립한다. 이어 강의의 Bielecki 노름으로 `L/λ < 1`을 만들면 수축사상이 되어 Banach 고정점정리를 적용할 수 있다. **사상 불변성 → 수축성 → 유일한 고정점**이 증명의 정확한 논리 순서다.

### Peano Theorem: Fixed-Point Comparison

강의 p.3 질문에 대한 답은 다음과 같다. Picard–Lindelöf는 Lipschitz 조건으로 *수축*을 확보하므로 Banach 정리를 쓴다. Peano에서는 일반적으로 수축성이 없지만, 연속성·균등유계성·동등연속성에 의해 근사해 집합이 상대콤팩트함을 보이고 Arzelà–Ascoli로 수렴 부분수열을 잡는다. 추가 관점으로 Schauder 고정점정리를 쓰려면 **비어 있지 않은 닫힌 볼록집합**, 그 집합을 자기 자신으로 보내는 **연속이며 콤팩트한 사상**이 필요하다. 이 후자의 방식은 강의 p.3에 완성된 증명으로 적혀 있지는 않으므로 여기서는 추가 설명으로 구분한다.

### Grönwall Inequality: Nonconstant Bound

강의 p.4–5는 상수형만이 아니라 다음 일반형도 제시한다.

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20u%28t%29%5Cle%5Calpha%28t%29%2B%5Cint_%7Bt_0%7D%5E%7Bt%7D%5Cbeta%28s%29u%28s%29%5C%2Cds)

![Equation](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20u%28t%29%5Cle%5Calpha%28t%29%2B%5Cint_%7Bt_0%7D%5E%7Bt%7D%5Calpha%28s%29%5Cbeta%28s%29%5Cexp%5C%21%5Cleft%28%5Cint_s%5Et%5Cbeta%28r%29%5C%2Cdr%5Cright%29ds)

여기서 `α, β, u`는 음이 아닌 연속함수이고, 특히 `α(t)=C`이면 적분항을 평가해 `u(t) ≤ C exp(∫β)`로 단순화된다. 강의의 'Continuous and Smooth Dependence on Parameters' 절은 이 부등식을 매개변수나 초기조건 변화가 해에 미치는 영향을 제어하는 수단으로 소개한다. 구체적인 매개변수 미분방정식은 강의 원문에는 전개되지 않는다.

### Finite-Time Blow-Up and the Limit of Euler Simulation

강의 p.6 마지막 예제는 `x′=1+x²`, `x(0)=0`으로 정확해 `tan(t)`가 `t=π/2` 직전에 발산한다. Euler는 유한한 격자값만 출력하므로 **실제 해가 유한시간에 발산한다는 사실을 수치 출력만으로 보장하거나 부정할 수 없다.** 이는 최대해의 연장 정리(p.4)와 연결된다.



초기값 문제에서 **존재성**, **유일성**, **연속 의존성**은 서로 다른 성질이다. Picard–Lindelöf 정리는 유일한 국소해를 보장하는 반면, Peano 정리는 연속성만으로 존재성을 보장한다. 수치해법의 수렴을 논하기 전에 정확히 무엇을 근사하는지 이해하는 것이 중요하다.

## References

- Raj Dhara, *A reminder of the basic facts of Ordinary Differential Equations (ODE)*, 2026-10-05.
- Raj Dhara, *Numerical Method Lab: ODE Numerical Schemes & Convergence*, 2026-10-07, Exercise 2.
