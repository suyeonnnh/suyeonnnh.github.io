---
title: "Numerical Methods : Fundamentals of Ordinary Differential Equations"
date: 2026-10-28 18:00:00 +0100
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

ODE의 해를 수치적으로 근사하기 전에 **해가 존재하는지, 유일한지, 초기조건의 변화에 안정적으로 반응하는지**부터 살펴본다. 강의 `Lec-1_ODE2026.pdf`의 Section 1–7을 순서대로 정리하고, Peano 현상에 관한 Lab Handout Exercise 2를 해당 이론 뒤에 배치했다.

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

초기 오차 `C`가 시간에 따라 확대될 수 있지만, 그 증가를 지수함수 형태의 상한으로 제어할 수 있다. 이는 **초기조건에 대한 연속 의존성과 수치 오차의 누적 분석**에 연결된다.

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

## Summary and Takeaways

초기값 문제에서 **존재성**, **유일성**, **연속 의존성**은 서로 다른 성질이다. Picard–Lindelöf 정리는 유일한 국소해를 보장하는 반면, Peano 정리는 연속성만으로 존재성을 보장한다. 수치해법의 수렴을 논하기 전에 정확히 무엇을 근사하는지 이해하는 것이 중요하다.

## References

- Raj Dhara, *A reminder of the basic facts of Ordinary Differential Equations (ODE)*, 2026-10-05.
- Raj Dhara, *Numerical Method Lab: ODE Numerical Schemes & Convergence*, 2026-10-07, Exercise 2.
