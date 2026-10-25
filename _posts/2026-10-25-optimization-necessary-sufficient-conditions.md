---
title: 'Optimization Theory : Necessary and sufficient optimality conditions'
date: 2026-10-25 12:00:00 +0200
categories: [Optimization Theory]
tags: [Optimization, Optimality Conditions, Hessian, Convexity, Coercivity]
math: true
---

# Optimization Theory : Necessary and sufficient optimality conditions

이번 글에서는 최적화의 기본 문제 설정부터 **필요조건(necessary condition)**, **충분조건(sufficient condition)**, **헤세 행렬(Hessian matrix)의 정부호성**, **볼록성(convexity)**, **강제성(coercivity)**까지 정리한다. 미적분에서 배운 극값 판정법을 먼저 복습한 후, 미분법만으로 풀기 어려운 문제와 전역 최적해를 판정하는 방법을 살펴본다.


> <strong>Definition | Optimization</strong>  
> 가능한 선택지 가운데 정해진 기준에 따라 가장 좋은 것을 선택하는 과정. 기준을 함수로 나타낼 수 있다면 그 함수를 **최소화 또는 최대화**하는 문제로 표현된다.
{: .prompt-tip }

## 1. Foundations of Optimization from Mathematical Analysis

### 1.1 Finding Local Minima and Maxima

미적분에서 배운 기본 절차는 다음과 같다.

1. 함수의 **1계 도함수(first-order derivatives)**, 다변수 함수라면 편도함수를 계산한다.
2. 도함수를 0으로 놓고 **임계점(critical points)**을 구한다.
3. 1변수에서는 2계 도함수의 부호를 확인하고, 2변수에서는 헤세 행렬식과 순수 2계 편도함수의 부호를 확인해 극값의 종류를 판정한다.

> <strong>Theorem | Second Derivative Test for Functions of Two Variables</strong>  
>
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20D%3Df_%7Bxx%7Df_%7Byy%7D-%28f_%7Bxy%7D%29%5E2." alt="D=f_{xx}f_{yy}-(f_{xy})^2." style="display:block;max-width:100%;height:auto;max-height:16rem;margin:0.75rem auto;object-fit:contain;">
>
> 임계점에서 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20D%3E0%2C%5C%20f_%7Bxx%7D%3E0" alt="D&gt;0,\ f_{xx}&gt;0" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">이면 **국소 최솟값**, <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20D%3E0%2C%5C%20f_%7Bxx%7D%3C0" alt="D&gt;0,\ f_{xx}&lt;0" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">이면 **국소 최댓값**이다. <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20D%3C0" alt="D&lt;0" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">이면 **안장점(saddle point)**이고, <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20D%3D0" alt="D=0" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">이면 이 판정법만으로 결론을 내릴 수 없다.
{: .prompt-tip }

> <strong>Example | Finding Critical Points of a Bivariate Polynomial</strong>  
>
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f%28x%2Cy%29%3Dx%5E3%2Bx%5E2y-y%5E2-4y." alt="f(x,y)=x^3+x^2y-y^2-4y." style="display:block;max-width:100%;height:auto;max-height:16rem;margin:0.75rem auto;object-fit:contain;">
>
> 먼저 편도함수를 계산한다.  
>
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f_x%3D3x%5E2%2B2xy%2C%5Cqquad%20f_y%3Dx%5E2-2y-4." alt="f_x=3x^2+2xy,\qquad f_y=x^2-2y-4." style="display:block;max-width:100%;height:auto;max-height:16rem;margin:0.75rem auto;object-fit:contain;">
>
> 따라서 임계점은 다음 연립방정식의 해이다.  
>
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%203x%5E2%2B2xy%3D0%2C%5Cqquad%20x%5E2-2y-4%3D0." alt="3x^2+2xy=0,\qquad x^2-2y-4=0." style="display:block;max-width:100%;height:auto;max-height:16rem;margin:0.75rem auto;object-fit:contain;">
>
> 첫 번째 식을 인수분해하면  
>
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x%283x%2B2y%29%3D0%5Cquad%5CRightarrow%5Cquad%20x%3D0%5C%20%5Ctext%7B%EB%98%90%EB%8A%94%7D%5C%20y%3D-%5Cfrac32x." alt="수식" style="display:block;max-width:100%;height:auto;max-height:16rem;margin:0.75rem auto;object-fit:contain;">
>
> **경우 1: <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x%3D0" alt="x=0" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">.** 두 번째 식에 대입하면 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20-2y-4%3D0" alt="-2y-4=0" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">이므로 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%280%2C-2%29" alt="(0,-2)" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">를 얻는다.  
> **경우 2: <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y%3D-%5Cfrac32x" alt="y=-\frac32x" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">.** 두 번째 식에 대입하면  
>
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x%5E2%2B3x-4%3D%28x%2B4%29%28x-1%29%3D0." alt="x^2+3x-4=(x+4)(x-1)=0." style="display:block;max-width:100%;height:auto;max-height:16rem;margin:0.75rem auto;object-fit:contain;">
>
> 따라서 나머지 임계점은 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%28-4%2C6%29" alt="(-4,6)" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">과 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%281%2C-%5Cfrac32%29" alt="(1,-\frac32)" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">이다.  
> 이제 2계 편도함수를 구한다.  
>
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f_%7Bxx%7D%3D6x%2B2y%2C%5Cquad%20f_%7Byy%7D%3D-2%2C%5Cquad%20f_%7Bxy%7D%3D2x%2C" alt="f_{xx}=6x+2y,\quad f_{yy}=-2,\quad f_{xy}=2x," style="display:block;max-width:100%;height:auto;max-height:16rem;margin:0.75rem auto;object-fit:contain;">
>
>
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20D%3D%286x%2B2y%29%28-2%29-%282x%29%5E2%3D-12x-4y-4x%5E2." alt="D=(6x+2y)(-2)-(2x)^2=-12x-4y-4x^2." style="display:block;max-width:100%;height:auto;max-height:16rem;margin:0.75rem auto;object-fit:contain;">
>
> 각 점에서의 결과는 다음과 같다.  
> 

| 임계점 | <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f_%7Bxx%7D" alt="f_{xx}" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;"> | <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20D" alt="D" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;"> | 판정 |
|:--|--:|--:|:--|
| <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%280%2C-2%29" alt="(0,-2)" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;"> | <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20-4" alt="-4" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;"> | <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%208%3E0" alt="8&gt;0" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;"> | 국소 최댓값 |
| <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%28-4%2C6%29" alt="(-4,6)" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;"> | <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20-12" alt="-12" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;"> | <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20-40%3C0" alt="-40&lt;0" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;"> | 안장점 |
| <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%281%2C-%5Cfrac32%29" alt="(1,-\frac32)" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;"> | <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%203" alt="3" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;"> | <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20-10%3C0" alt="-10&lt;0" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;"> | 안장점 |


하지만 **도함수를 0으로 만드는 방식으로 모든 극값을 찾을 수는 없다.** 극값이 존재해도 해당 점에서 미분이 불가능할 수 있다.

> <strong>Example | A Nondifferentiable Local Minimum</strong>  
>
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f%28x%2Cy%29%3D%5Csqrt%7Bx%5E2%2By%5E2%7D." alt="f(x,y)=\sqrt{x^2+y^2}." style="display:block;max-width:100%;height:auto;max-height:16rem;margin:0.75rem auto;object-fit:contain;">
>
> 이 함수는 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%280%2C0%29" alt="(0,0)" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">에서 최솟값 0을 가지지만, 원점에서 편도함수가 존재하지 않는다. 따라서 미분 가능한 점의 정지조건만 검사하면 이 극솟점을 놓칠 수 있다.
{: .prompt-info }

### 1.2 Finding Extrema on Compact Sets

콤팩트 집합 위에서 연속함수의 최댓값과 최솟값을 찾을 때는 **내부의 임계점, 미분 불가능한 점, 경계의 극값 후보**를 모두 조사하고 함수값을 비교한다. 다른 방법으로는 라그랑주 승수법(Lagrange multipliers)이 있다.

> <strong>Example | Global Extrema on the Closed Unit Disk</strong>  
>
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f%28x%2Cy%29%3Dx%5E3-y%5E2%2Bxy%5E2%2C%5Cqquad%20x%5E2%2By%5E2%5Cle1." alt="f(x,y)=x^3-y^2+xy^2,\qquad x^2+y^2\le1." style="display:block;max-width:100%;height:auto;max-height:16rem;margin:0.75rem auto;object-fit:contain;">
>
> **Step 1 | Interior Critical Points**  
>
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f_x%3D3x%5E2%2By%5E2%3D0%2C%5Cqquad%20f_y%3D-2y%2B2xy%3D2y%28x-1%29%3D0." alt="f_x=3x^2+y^2=0,\qquad f_y=-2y+2xy=2y(x-1)=0." style="display:block;max-width:100%;height:auto;max-height:16rem;margin:0.75rem auto;object-fit:contain;">
>
> 첫 번째 식이 0이려면 실수 범위에서 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x%3Dy%3D0" alt="x=y=0" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">이어야 한다. 따라서 내부 임계점은 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%280%2C0%29" alt="(0,0)" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">이고 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f%280%2C0%29%3D0" alt="f(0,0)=0" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">이다. 두 번째 식을 기준으로 나누어 보더라도 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y%3D0" alt="y=0" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">이면 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x%3D0" alt="x=0" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">, <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x%3D1" alt="x=1" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">이면 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%203%2By%5E2%3D0" alt="3+y^2=0" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">이라는 모순이 나온다.  
> **Step 2 | Boundary Analysis**  
> 경계 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x%5E2%2By%5E2%3D1" alt="x^2+y^2=1" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">에서는 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y%5E2%3D1-x%5E2" alt="y^2=1-x^2" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">, <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20-1%5Cle%20x%5Cle1" alt="-1\le x\le1" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">이다.  
>
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cbegin%7Baligned%7Dg%28x%29%26%3Dx%5E3-%281-x%5E2%29%2Bx%281-x%5E2%29%5C%5C%26%3Dx%5E2%2Bx-1.%5Cend%7Baligned%7D" alt="수식" style="display:block;max-width:100%;height:auto;max-height:16rem;margin:0.75rem auto;object-fit:contain;">
>
> 도함수를 0으로 놓으면  
>
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20g%27%28x%29%3D2x%2B1%3D0%5Cquad%5CRightarrow%5Cquad%20x%3D-%5Cfrac12%2C" alt="g&#x27;(x)=2x+1=0\quad\Rightarrow\quad x=-\frac12," style="display:block;max-width:100%;height:auto;max-height:16rem;margin:0.75rem auto;object-fit:contain;">
>
>
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y%3D%5Cpm%5Csqrt%7B1-%5Cfrac14%7D%3D%5Cpm%5Cfrac%7B%5Csqrt3%7D%7B2%7D." alt="y=\pm\sqrt{1-\frac14}=\pm\frac{\sqrt3}{2}." style="display:block;max-width:100%;height:auto;max-height:16rem;margin:0.75rem auto;object-fit:contain;">
>
> 따라서 경계 임계점에서 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f%28-%5Cfrac12%2C%5Cpm%5Cfrac%7B%5Csqrt3%7D%7B2%7D%29%3D-%5Cfrac54" alt="f(-\frac12,\pm\frac{\sqrt3}{2})=-\frac54" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">이다. 또한 구간 끝점에 대응하는 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%28-1%2C0%29" alt="(-1,0)" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">과 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%281%2C0%29" alt="(1,0)" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">에서 함수값은 각각 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20-1" alt="-1" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">, <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%201" alt="1" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">이다.  
> **Step 3 | Comparing Candidates**  
> 

| 후보 | 함수값 |
|:--|--:|
| <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%280%2C0%29" alt="(0,0)" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;"> | <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%200" alt="0" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;"> |
| <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%28-%5Cfrac12%2C%5Cpm%5Cfrac%7B%5Csqrt3%7D%7B2%7D%29" alt="(-\frac12,\pm\frac{\sqrt3}{2})" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;"> | <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20-%5Cfrac54" alt="-\frac54" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;"> |
| <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%28-1%2C0%29" alt="(-1,0)" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;"> | <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20-1" alt="-1" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;"> |
| <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%281%2C0%29" alt="(1,0)" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;"> | <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%201" alt="1" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;"> |

> 
> **전역 최댓값은 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%201" alt="1" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;"> (위치 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%281%2C0%29" alt="(1,0)" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">), 전역 최솟값은 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20-%5Cfrac54" alt="-\frac54" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;"> (위치 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%28-%5Cfrac12%2C%5Cpm%5Cfrac%7B%5Csqrt3%7D%7B2%7D%29" alt="(-\frac12,\pm\frac{\sqrt3}{2})" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">)**이다.
{: .prompt-info }

### 1.3 A Question Leading to Convexity

미적분에서는 볼록함수의 정의와, 1변수에서 이를 실용적으로 판단하는 방법도 배운다. 이 개념은 나중에 **국소 최솟값을 전역 최솟값으로 확장할 수 있는지**를 판단하는 중요한 도구가 된다.

## 2. Why Optimization Theory Is Needed

미분법으로 극값의 후보를 찾는 것은 출발점이지만, 실제 문제에서는 방정식이 어렵거나 제약조건이 있거나 변수가 정수여야 한다. 최적화 이론은 이런 경우를 체계적으로 다룬다.

> <strong>Example | A Maximization Problem with Nonlinear Equations</strong>  
> 다음 목적함수를 최대화한다고 하자. (<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20a%3E0" alt="a&gt;0" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">)  
>
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f%28p%2Ca%29%3D%282000%2B%5Cln%20a-20p%29%28p-2%29-a-20000." alt="f(p,a)=(2000+\ln a-20p)(p-2)-a-20000." style="display:block;max-width:100%;height:auto;max-height:16rem;margin:0.75rem auto;object-fit:contain;">
>
> 1계 조건은  
>
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f_p%3D2000%2B%5Cln%20a%2B40-40p%3D0%2C" alt="f_p=2000+\ln a+40-40p=0," style="display:block;max-width:100%;height:auto;max-height:16rem;margin:0.75rem auto;object-fit:contain;">
>
>
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f_a%3D%5Cfrac%7Bp-2%7D%7Ba%7D-1%3D0." alt="f_a=\frac{p-2}{a}-1=0." style="display:block;max-width:100%;height:auto;max-height:16rem;margin:0.75rem auto;object-fit:contain;">
>
> 두 번째 조건으로부터 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20a%3Dp-2" alt="a=p-2" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">를 얻고, 이를 대입하면  
>
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%202039%2B%5Cln%28p-2%29%3D40p" alt="2039+\ln(p-2)=40p" style="display:block;max-width:100%;height:auto;max-height:16rem;margin:0.75rem auto;object-fit:contain;">
>
> 와 같은 비선형 방정식이 나온다. 이 식은 단순한 대수적 조작으로 풀기 어렵다. 일반적인 최적화에서는 수치적으로도 풀기 어려운 방정식을 자주 만난다.  
> *주:* 원자료의 중간 상수 표기에 산술상 차이가 있다. <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f_p" alt="f_p" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">를 직접 미분하면 상수항은 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%202040" alt="2040" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">이므로, 정확한 치환식은 **<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%202040%2B%5Cln%28p-2%29%3D40p" alt="2040+\ln(p-2)=40p" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">**이다. 위의 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%202039" alt="2039" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">는 원자료에 제시된 식이다.
{: .prompt-info }

> <strong>Example | Complications from Constraints</strong>  
>
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f%28x%2Cy%29%3Dx%5E3-y%5E2%2Bxy%2C%5Cqquad%20x%5E2%2By%5E2%5Cle1." alt="f(x,y)=x^3-y^2+xy,\qquad x^2+y^2\le1." style="display:block;max-width:100%;height:auto;max-height:16rem;margin:0.75rem auto;object-fit:contain;">
>
> 내부 임계점은  
>
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f_x%3D3x%5E2%2By%3D0%2C%5Cquad%20f_y%3D-2y%2Bx%3D0" alt="f_x=3x^2+y=0,\quad f_y=-2y+x=0" style="display:block;max-width:100%;height:auto;max-height:16rem;margin:0.75rem auto;object-fit:contain;">
>
> 을 풀어 얻는다. <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x%3D2y" alt="x=2y" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">를 대입하면  
>
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%2012y%5E2%2By%3Dy%2812y%2B1%29%3D0%2C" alt="12y^2+y=y(12y+1)=0," style="display:block;max-width:100%;height:auto;max-height:16rem;margin:0.75rem auto;object-fit:contain;">
>
> 따라서 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%280%2C0%29" alt="(0,0)" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">과 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%28-%5Cfrac16%2C-%5Cfrac1%7B12%7D%29" alt="(-\frac16,-\frac1{12})" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">가 나온다.  
> 경계 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x%5E2%2By%5E2%3D1" alt="x^2+y^2=1" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">에서는  
>
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y%3D%5Cpm%5Csqrt%7B1-x%5E2%7D%2C%5Cqquad%20-1%5Cle%20x%5Cle1" alt="y=\pm\sqrt{1-x^2},\qquad -1\le x\le1" style="display:block;max-width:100%;height:auto;max-height:16rem;margin:0.75rem auto;object-fit:contain;">
>
> 이므로 다음 두 함수의 임계점을 조사해야 한다.  
>
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20g_1%28x%29%3Dx%5E3%2Bx%5E2-1%2Bx%5Csqrt%7B1-x%5E2%7D%2C" alt="g_1(x)=x^3+x^2-1+x\sqrt{1-x^2}," style="display:block;max-width:100%;height:auto;max-height:16rem;margin:0.75rem auto;object-fit:contain;">
>
>
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20g_2%28x%29%3Dx%5E3%2Bx%5E2-1-x%5Csqrt%7B1-x%5E2%7D." alt="g_2(x)=x^3+x^2-1-x\sqrt{1-x^2}." style="display:block;max-width:100%;height:auto;max-height:16rem;margin:0.75rem auto;object-fit:contain;">
>
> 첫 번째 함수의 도함수는  
>
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20g_1%27%28x%29%3D3x%5E2%2B2x%2B%5Csqrt%7B1-x%5E2%7D-%5Cfrac%7Bx%5E2%7D%7B%5Csqrt%7B1-x%5E2%7D%7D." alt="수식" style="display:block;max-width:100%;height:auto;max-height:16rem;margin:0.75rem auto;object-fit:contain;">
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20g_1%27%28x%29%3D0" alt="g_1&#x27;(x)=0" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">에서 제곱근을 제거해 정리하면 다음과 같은 높은 차수의 식으로 이어진다.  
>
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%209x%5E6%2B12x%5E5-x%5E4-12x%5E3-8x%5E2%2B1%3D0." alt="9x^6+12x^5-x^4-12x^3-8x^2+1=0." style="display:block;max-width:100%;height:auto;max-height:16rem;margin:0.75rem auto;object-fit:contain;">
>
> 식을 수치적으로 풀 수 있더라도 부호 검증, 경계 끝점, 두 번째 곡선까지 고려해야 한다. 즉 **경계를 직접 매개화하는 방식은 계산이 매우 번거로워질 수 있다.**
{: .prompt-info }

> <strong>Example | A Retailer’s Package-Size Decision</strong>  
> 소매업자는 상품을 100개, 50개, 20개, 10개 묶음으로 구매하며 각 묶음 가격은 각각 \\$160, \\$90, \\$45, \\$28이다. 상품은 낱개당 \\$5에 판매하고, 주문량은 3,791개다. 각 포장 단위를 몇 묶음씩 구매해야 이익이 최대가 될까?  
> 포장별 구매 묶음 수를 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x_%7B100%7D%2Cx_%7B50%7D%2Cx_%7B20%7D%2Cx_%7B10%7D" alt="x_{100},x_{50},x_{20},x_{10}" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">, 실제 판매 개수를 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20u" alt="u" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">라고 하면 목적함수는  
>
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cmax%5C%205u-160x_%7B100%7D-90x_%7B50%7D-45x_%7B20%7D-28x_%7B10%7D" alt="수식" style="display:block;max-width:100%;height:auto;max-height:16rem;margin:0.75rem auto;object-fit:contain;">
>
> 이며 제약조건은  
>
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20u%5Cle3791%2C" alt="u\le3791," style="display:block;max-width:100%;height:auto;max-height:16rem;margin:0.75rem auto;object-fit:contain;">
>
>
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20100x_%7B100%7D%2B50x_%7B50%7D%2B20x_%7B20%7D%2B10x_%7B10%7D%5Cge%20u%2C" alt="100x_{100}+50x_{50}+20x_{20}+10x_{10}\ge u," style="display:block;max-width:100%;height:auto;max-height:16rem;margin:0.75rem auto;object-fit:contain;">
>
>
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x_%7B100%7D%2Cx_%7B50%7D%2Cx_%7B20%7D%2Cx_%7B10%7D%5Cge0." alt="x_{100},x_{50},x_{20},x_{10}\ge0." style="display:block;max-width:100%;height:auto;max-height:16rem;margin:0.75rem auto;object-fit:contain;">
>
> 구매 묶음 수는 실제로 **음이 아닌 정수**여야 한다. 선형 목적함수는 일반적으로 도함수가 0이 되지 않으므로 이 경우 미분만으로 답을 구할 수 없다. 연속 선형계획의 최적해는 적절한 조건하에서 가능영역의 꼭짓점에서 찾을 수 있지만, 변수 수에 따라 꼭짓점의 개수가 기하급수적으로 증가할 수 있다. 더구나 정수 제약은 별도로 처리해야 한다.
{: .prompt-info }

### 2.1 Goals and Scope of Optimization Theory

이 강의에서 강조하는 세 가지 목표는 다음과 같다.

1. 도함수를 활용해 최솟점과 최댓점을 찾는 **수치적 도구**를 제공한다.
2. **제약 최적화 문제**를 미분 기반 방법으로 다룰 수 있도록 변형하는 분석 도구를 제공한다.
3. 미분 기반 접근이 적합하지 않은 **선형·이산 최적화** 문제를 효율적으로 해결하는 방법을 제공한다.

전체 학습 주제는 (1) 최적성 필요·충분조건, (2) 무제약 최적화를 위한 경사 기반 방법, (3) 선형계획법과 정수 선형계획법, (4) 라그랑주 승수법과 Karush–Kuhn–Tucker(KKT) 조건, (5) 이차계획법, (6) 일반적인 제약 최적화 방법, (7) 이산 최적화를 위한 휴리스틱 방법이다.

> **Terminology | Basic Optimization Notation**  
>
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cmin_%7Bx%5Cin%20D%7Df%28x%29%5Cqquad%5Ctext%7B%EB%98%90%EB%8A%94%7D%5Cqquad%5Cmax_%7Bx%5Cin%20D%7Df%28x%29." alt="수식" style="display:block;max-width:100%;height:auto;max-height:16rem;margin:0.75rem auto;object-fit:contain;">
>
> - **Objective function(목적함수)**: 최소화·최대화할 대상인 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f" alt="f" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">; 최소화 시 *cost*, 최대화 시 *reward*라고도 한다.  
> - **Variable(변수)**: 선택하거나 조정하는 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x" alt="x" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">.  
> - **Feasible set(가능집합)**: 허용되는 선택의 집합 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20D" alt="D" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">.  
> - **Constraint(제약조건)**: <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20D" alt="D" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">를 정의하는 등식 또는 부등식.  
> - **<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Coperatorname%2A%7Bargmin%7D" alt="\operatorname*{argmin}" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;"> / <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Coperatorname%2A%7Bargmax%7D" alt="\operatorname*{argmax}" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">**: 가장 작은/큰 **함수값을 만들어 내는 변수(해)**.  
> - **Minimizer / Maximizer(최솟점 / 최댓점)**: 최적해를 이루는 점.  
> - **Constrained / Unconstrained**: 제약 / 무제약.  
> - **Local / Global**: 국소 / 전역.  
> 주의: <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cmin%20f%28x%29" alt="\min f(x)" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">는 최적 **값**이고, <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Coperatorname%2A%7Bargmin%7Df%28x%29" alt="\operatorname*{argmin}f(x)" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">는 그 값을 달성하는 **점 또는 점들의 집합**이다.
{: .prompt-tip }

## 3. Necessary and Sufficient Optimality Conditions

### 3.1 First- and Second-Order Necessary Conditions

> <strong>Theorem | Necessary Optimality Condition</strong>  
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x%5E%2A%5Cin%5Cmathbb%20R%5En" alt="x^*\in\mathbb R^n" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">이 **무제약 국소 최솟점 또는 최댓점**이고 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f" alt="f" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">가 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x%5E%2A" alt="x^*" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">를 포함하는 열린 근방에서 연속 미분 가능하다고 하자. 그러면 반드시  
>
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cboxed%7B%5Cnabla%20f%28x%5E%2A%29%3D0%7D" alt="\boxed{\nabla f(x^*)=0}" style="display:block;max-width:100%;height:auto;max-height:16rem;margin:0.75rem auto;object-fit:contain;">
>
> 이다. 여기에 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f" alt="f" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">가 2번 연속 미분 가능하다는 조건까지 추가하면,  
> - 국소 최솟점에서 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cnabla%5E2f%28x%5E%2A%29" alt="\nabla^2f(x^*)" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">는 **양의 준정부호(positive semidefinite)**이고,
> - 국소 최댓점에서 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cnabla%5E2f%28x%5E%2A%29" alt="\nabla^2f(x^*)" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">는 **음의 준정부호(negative semidefinite)**이다.  
> 즉, 국소 극값이라면 반드시 만족해야 하는 조건이지만, 이것만으로 극값임이 보장되지는 않는다.
{: .prompt-tip }

여기서 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cnabla%20f" alt="\nabla f" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">는 모든 1계 편도함수를 모은 **그래디언트(gradient)**이고 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cnabla%5E2f" alt="\nabla^2f" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">는 모든 2계 편도함수로 구성한 **헤세 행렬(Hessian)**이다.

### 3.2 Second-Order Sufficient Conditions

> <strong>Theorem | Second-Order Sufficient Optimality Condition</strong>  
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f" alt="f" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">가 열린 근방에서 두 번 연속 미분 가능하고 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cnabla%20f%28x%5E%2A%29%3D0" alt="\nabla f(x^*)=0" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">이라고 하자.  
> - <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cnabla%5E2%20f%28x%5E%2A%29" alt="\nabla^2 f(x^*)" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">가 **양의 정부호(positive definite)**이면, <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x%5E%2A" alt="x^*" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">는 **엄격한 국소 최솟점(strict local minimizer)**이다.
> - <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cnabla%5E2%20f%28x%5E%2A%29" alt="\nabla^2 f(x^*)" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">가 **음의 정부호(negative definite)**이면, <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x%5E%2A" alt="x^*" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">는 **엄격한 국소 최댓점(strict local maximizer)**이다.  
> 준정부호인 경우에는 이 충분조건만으로 엄격한 극값을 확정할 수 없다.
{: .prompt-tip }

**필요조건**은 최적점이 *반드시 만족해야 할 조건*이고, **충분조건**은 해당 조건이 성립하면 *최적점임을 결론 내릴 수 있는 조건*이다. 필요조건을 통과해도 안장점일 수 있지만, 충분조건을 만족하면 그 점의 국소 최적성은 보장된다.

## 4. Definiteness of Matrices

### 4.1 Definitions

> <strong>Definition | Definiteness of a Symmetric Matrix</strong>  
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20A" alt="A" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">를 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20n%5Ctimes%20n" alt="n\times n" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;"> 대칭 실수행렬이라고 하자.  
> - **양의 정부호(positive definite)**: <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x%5Cne0" alt="x\ne0" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">인 모든 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x%5Cin%5Cmathbb%20R%5En" alt="x\in\mathbb R^n" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">에 대해 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x%5ETAx%3E0" alt="x^TAx&gt;0" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">.
> - **양의 준정부호(positive semidefinite)**: 모든 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x%5Cin%5Cmathbb%20R%5En" alt="x\in\mathbb R^n" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">에 대해 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x%5ETAx%5Cge0" alt="x^TAx\ge0" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">.
> - **음의 정부호(negative definite)**: <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x%5Cne0" alt="x\ne0" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">인 모든 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x%5Cin%5Cmathbb%20R%5En" alt="x\in\mathbb R^n" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">에 대해 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x%5ETAx%3C0" alt="x^TAx&lt;0" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">.
> - **음의 준정부호(negative semidefinite)**: 모든 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x%5Cin%5Cmathbb%20R%5En" alt="x\in\mathbb R^n" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">에 대해 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x%5ETAx%5Cle0" alt="x^TAx\le0" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">.
{: .prompt-tip }

### 4.2 Practical Tests for Definiteness

> <strong>Theorem | Sylvester’s Criterion and Eigenvalue Tests</strong>  
> **Sylvester’s Criterion — Strict Definiteness**  
> 선도 주인수식(leading principal minors) <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5CDelta_k%3D%5Cdet%28A_%7B1%3Ak%2C1%3Ak%7D%29" alt="\Delta_k=\det(A_{1:k,1:k})" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">을 계산한다.  
> - **양의 정부호** <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5CLongleftrightarrow" alt="\Longleftrightarrow" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;"> 모든 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5CDelta_k%3E0" alt="\Delta_k&gt;0" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">.
> - **음의 정부호** <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5CLongleftrightarrow" alt="\Longleftrightarrow" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;"> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%28-1%29%5Ek%5CDelta_k%3E0" alt="(-1)^k\Delta_k&gt;0" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;"> (홀수 차수는 음수, 짝수 차수는 양수).  
> **Note on Semidefiniteness:** 원자료는 선도 주인수식의 부호만으로 준정부호를 판정할 수 있는 것처럼 적고 있지만, 일반적으로는 **모든 주인수식(principal minors)**을 확인해야 한다. 양의 준정부호는 모든 주인수식 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cge0" alt="\ge0" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">, 음의 준정부호는 모든 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20k" alt="k" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">차 주인수식에 대해 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%28-1%29%5Ek%5Cdet%5Cge0" alt="(-1)^k\det\ge0" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">이다.  
> **Eigenvalue Test**  
> - 양의 정부호: 모든 고유값 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%200" alt="0" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">; 양의 준정부호: 모든 고유값 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cge0" alt="\ge0" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">.
> - 음의 정부호: 모든 고유값 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%3C0" alt="&lt;0" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">; 음의 준정부호: 모든 고유값 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cle0" alt="\le0" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">.
{: .prompt-tip }

> <strong>Example | Testing Negative Definiteness of a 4×4 Matrix</strong>  
>
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20A%3D%5Cbegin%7Bpmatrix%7D-2%26-1%261%261%5C%5C-1%26-4%260%26-2%5C%5C1%260%26-2%260%5C%5C1%26-2%260%26-4%5Cend%7Bpmatrix%7D." alt="수식" style="display:block;max-width:100%;height:auto;max-height:16rem;margin:0.75rem auto;object-fit:contain;">
>
> 선도 주인수식 네 개를 계산하면  
>
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5CDelta_1%3D-2%3C0%2C%5Cqquad%5CDelta_2%3D%5Cbegin%7Bvmatrix%7D-2%26-1%5C%5C-1%26-4%5Cend%7Bvmatrix%7D%3D7%3E0%2C" alt="수식" style="display:block;max-width:100%;height:auto;max-height:16rem;margin:0.75rem auto;object-fit:contain;">
>
>
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5CDelta_3%3D%5Cbegin%7Bvmatrix%7D-2%26-1%261%5C%5C-1%26-4%260%5C%5C1%260%26-2%5Cend%7Bvmatrix%7D%3D-10%3C0%2C" alt="수식" style="display:block;max-width:100%;height:auto;max-height:16rem;margin:0.75rem auto;object-fit:contain;">
>
>
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5CDelta_4%3D%5Cdet%28A%29%3D12%3E0." alt="\Delta_4=\det(A)=12&gt;0." style="display:block;max-width:100%;height:auto;max-height:16rem;margin:0.75rem auto;object-fit:contain;">
>
> 부호가 음·양·음·양으로 번갈아 나타나므로 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20A" alt="A" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">는 **음의 정부호(negative definite)**이다.
{: .prompt-info }

어떤 선도 주인수식이 **0**이면 엄격한 정부호 판정은 통과하지 못하며, 준정부호를 판정하려면 빠진 주인수식을 추가로 확인해야 한다. 4×4 대칭행렬에서는 선도 주인수식 4개 외에 확인할 수 있는 비선도 주인수식이 **11개** 더 있다(전체 비어 있지 않은 주인수식 15개). 3×3에서는 1×1 비선도 주인수식 2개와 2×2 비선도 주인수식 2개가 더 있다. 2×2에서는 선도식 외에 다른 1×1 주인수식 하나도 별도로 고려해야 하므로, **원자료의 '2×2에서는 추가 확인이 필요 없다'는 문장은 일반적인 준정부호 판정에 대해 성립하지 않는다.** 예를 들어 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Coperatorname%7Bdiag%7D%280%2C-1%29" alt="\operatorname{diag}(0,-1)" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">의 두 선도식은 모두 0이지만 양의 준정부호는 아니다.

## 5. Applying Optimality Conditions to Functions

> <strong>Example | A Local Minimum of a Three-Variable Function</strong>  
>
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f%28x%2Cy%2Cz%29%3D%28x%5E2-z%29%5E2%2B%28y-z%29%5E2%2B2z." alt="f(x,y,z)=(x^2-z)^2+(y-z)^2+2z." style="display:block;max-width:100%;height:auto;max-height:16rem;margin:0.75rem auto;object-fit:contain;">
>
> **Step 1 | Stationary Points**  
>
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cbegin%7Bcases%7Df_x%3D4x%28x%5E2-z%29%3D0%2C%5C%5Cf_y%3D2%28y-z%29%3D0%2C%5C%5Cf_z%3D-2%28x%5E2-z%29-2%28y-z%29%2B2%3D0.%5Cend%7Bcases%7D" alt="수식" style="display:block;max-width:100%;height:auto;max-height:16rem;margin:0.75rem auto;object-fit:contain;">
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y%3Dz" alt="y=z" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">를 세 번째 식에 대입하면 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x%5E2-z%3D1" alt="x^2-z=1" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">이다. 따라서 첫 번째 식으로부터 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x%3D0" alt="x=0" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">, 그리고 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20z%3Dy%3D-1" alt="z=y=-1" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">을 얻는다. 정지점은 **<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%280%2C-1%2C-1%29" alt="(0,-1,-1)" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;"> 하나**다.  
> **Step 2 | Computing the Hessian**  
>
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cnabla%5E2f%28x%2Cy%2Cz%29%3D%5Cbegin%7Bpmatrix%7D12x%5E2-4z%260%26-4x%5C%5C0%262%26-2%5C%5C-4x%26-2%264%5Cend%7Bpmatrix%7D." alt="수식" style="display:block;max-width:100%;height:auto;max-height:16rem;margin:0.75rem auto;object-fit:contain;">
>
> 정지점에서  
>
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cnabla%5E2f%280%2C-1%2C-1%29%3D%5Cbegin%7Bpmatrix%7D4%260%260%5C%5C0%262%26-2%5C%5C0%26-2%264%5Cend%7Bpmatrix%7D." alt="수식" style="display:block;max-width:100%;height:auto;max-height:16rem;margin:0.75rem auto;object-fit:contain;">
>
> 이 행렬의 선도 주인수식은 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5CDelta_1%3D4" alt="\Delta_1=4" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">, <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5CDelta_2%3D8" alt="\Delta_2=8" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">, <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5CDelta_3%3D16" alt="\Delta_3=16" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">으로 모두 양수다. 따라서 양의 정부호이고, **<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%280%2C-1%2C-1%29" alt="(0,-1,-1)" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">은 엄격한 국소 최솟점**이다.  
> 하지만 이것이 **전역 최솟점**인지 알려면 추가적인 논증이 필요하다.
{: .prompt-info }

### 5.1 Why Local Optimality Is Not Enough

> <strong>Example | A Global Minimum Without the Sufficient Condition</strong>  
>
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f%28x%2Cy%29%3D%28x%5E2%2By%5E2-1%29%5E2." alt="f(x,y)=(x^2+y^2-1)^2." style="display:block;max-width:100%;height:auto;max-height:16rem;margin:0.75rem auto;object-fit:contain;">
>
> 제곱은 항상 0 이상이므로 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f%28x%2Cy%29%5Cge0" alt="f(x,y)\ge0" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">. 그리고 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x%5E2%2By%5E2%3D1" alt="x^2+y^2=1" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">인 모든 점에서 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f%3D0" alt="f=0" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">이다. 따라서 **단위원 위의 모든 점이 전역 최솟점**이다.  
> 편도함수는  
>
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f_x%3D4x%28x%5E2%2By%5E2-1%29%2C%5Cquad%20f_y%3D4y%28x%5E2%2By%5E2-1%29%2C" alt="f_x=4x(x^2+y^2-1),\quad f_y=4y(x^2+y^2-1)," style="display:block;max-width:100%;height:auto;max-height:16rem;margin:0.75rem auto;object-fit:contain;">
>
>
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f_%7Bxx%7D%3D4%28x%5E2%2By%5E2-1%29%2B8x%5E2%2C%5Cquad%20f_%7Byy%7D%3D4%28x%5E2%2By%5E2-1%29%2B8y%5E2%2C%5Cquad%20f_%7Bxy%7D%3D8xy." alt="수식" style="display:block;max-width:100%;height:auto;max-height:16rem;margin:0.75rem auto;object-fit:contain;">
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%281%2C0%29" alt="(1,0)" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">에서 그래디언트는 0이고 헤세 행렬은  
>
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cnabla%5E2f%281%2C0%29%3D%5Cbegin%7Bpmatrix%7D8%260%5C%5C0%260%5Cend%7Bpmatrix%7D." alt="수식" style="display:block;max-width:100%;height:auto;max-height:16rem;margin:0.75rem auto;object-fit:contain;">
>
> 이는 **양의 준정부호이지만 양의 정부호는 아니다.** 따라서 2계 충분조건을 만족하지 않지만, 이미 전역 최솟점임을 확인했다.  
> 또한 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%28x%2Cy%29%3D%281-%5Cvarepsilon%2C0%29" alt="(x,y)=(1-\varepsilon,0)" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">, <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%200%3C%5Cvarepsilon%3C2" alt="0&lt;\varepsilon&lt;2" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">이면  
>
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cnabla%5E2%20f%281-%5Cvarepsilon%2C0%29%3D%5Cbegin%7Bpmatrix%7D4%28%5Cvarepsilon%5E2-2%5Cvarepsilon%29%2B8%281-%5Cvarepsilon%29%5E2%260%5C%5C0%264%28%5Cvarepsilon%5E2-2%5Cvarepsilon%29%5Cend%7Bpmatrix%7D." alt="수식" style="display:block;max-width:100%;height:auto;max-height:16rem;margin:0.75rem auto;object-fit:contain;">
>
> 오른쪽 아래 성분이 음수이므로 이 행렬은 양의 준정부호가 아니다. 즉 이 함수는 단위원 내부의 이런 점에서는 볼록함수가 아니다.  
> *주:* 원자료의 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%281%2C0%29" alt="(1,0)" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">에서 헤세 행렬 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%284%2C0%3B0%2C0%29" alt="(4,0;0,0)" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;"> 표기는 도함수에 비춰 볼 때 오타이며, 위처럼 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%288%2C0%3B0%2C0%29" alt="(8,0;0,0)" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">이 정확하다. 원자료의 근방 헤세 행렬 평가점도 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%281-%5Cvarepsilon%2C0%29" alt="(1-\varepsilon,0)" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">으로 바로잡았다.
{: .prompt-info }

> <strong>Example | A Strict Local Minimum Without a Global Minimum</strong>  
>
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f%28x%2Cy%29%3Dx%5E2%2By%5E2%281-x%29%5E3." alt="f(x,y)=x^2+y^2(1-x)^3." style="display:block;max-width:100%;height:auto;max-height:16rem;margin:0.75rem auto;object-fit:contain;">
>
> 정지점은  
>
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f_x%3D2x-3y%5E2%281-x%29%5E2%3D0%2C%5Cqquad%20f_y%3D2y%281-x%29%5E3%3D0" alt="f_x=2x-3y^2(1-x)^2=0,\qquad f_y=2y(1-x)^3=0" style="display:block;max-width:100%;height:auto;max-height:16rem;margin:0.75rem auto;object-fit:contain;">
>
> 으로 구한다. <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x%3D1" alt="x=1" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">이면 첫 번째 식에서 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%202%3D0" alt="2=0" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">이라는 모순이 생기고, <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y%3D0" alt="y=0" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">이면 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x%3D0" alt="x=0" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">이므로 정지점은 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%280%2C0%29" alt="(0,0)" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;"> 하나다.  
> 헤세 행렬은  
>
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cnabla%5E2%20f%28x%2Cy%29%3D%5Cbegin%7Bpmatrix%7D2%2B6y%5E2%281-x%29%26-6y%281-x%29%5E2%5C%5C-6y%281-x%29%5E2%262%281-x%29%5E3%5Cend%7Bpmatrix%7D." alt="수식" style="display:block;max-width:100%;height:auto;max-height:16rem;margin:0.75rem auto;object-fit:contain;">
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%280%2C0%29" alt="(0,0)" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">에서는  
>
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cnabla%5E2%20f%280%2C0%29%3D%5Cbegin%7Bpmatrix%7D2%260%5C%5C0%262%5Cend%7Bpmatrix%7D%2C" alt="수식" style="display:block;max-width:100%;height:auto;max-height:16rem;margin:0.75rem auto;object-fit:contain;">
>
> 따라서 **엄격한 국소 최솟점**이다. 그러나  
>
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f%28x%2C1%29%3Dx%5E2%2B%281-x%29%5E3%5Clongrightarrow-%5Cinfty%5Cquad%28x%5Cto%2B%5Cinfty%29." alt="수식" style="display:block;max-width:100%;height:auto;max-height:16rem;margin:0.75rem auto;object-fit:contain;">
>
> 함수값이 아래로 한없이 내려갈 수 있으므로 **전역 최솟값은 존재하지 않는다.**
{: .prompt-info }

이 예제들은 **'국소 최적'과 '전역 최적'이 다르다**는 점, 그리고 2계 충분조건이 전역 최적성을 보장하지 않는다는 점을 보여 준다.

## 6. Convexity

### 6.1 Convex and Concave Functions

> <strong>Definition | Convex and Concave Functions</strong>  
> 두 점 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x%2Cy" alt="x,y" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">와 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%200%3C%5Clambda%3C1" alt="0&lt;\lambda&lt;1" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">에 대해  
>
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cboxed%7Bf%28%5Clambda%20x%2B%281-%5Clambda%29y%29%5Cle%5Clambda%20f%28x%29%2B%281-%5Clambda%29f%28y%29%7D" alt="수식" style="display:block;max-width:100%;height:auto;max-height:16rem;margin:0.75rem auto;object-fit:contain;">
>
> 를 만족하면 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f" alt="f" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">는 **볼록함수(convex function)**다.  
> 반대로  
>
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f%28%5Clambda%20x%2B%281-%5Clambda%29y%29%5Cge%5Clambda%20f%28x%29%2B%281-%5Clambda%29f%28y%29" alt="수식" style="display:block;max-width:100%;height:auto;max-height:16rem;margin:0.75rem auto;object-fit:contain;">
>
> 를 만족하면 **오목함수(concave function)**다.  
> 서로 다른 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x%5Cne%20y" alt="x\ne y" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">에서 부등호가 항상 엄격하면 각각 **엄격한 볼록함수(strictly convex)**, **엄격한 오목함수(strictly concave)**라고 한다. 원자료의 '모든 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x%2Cy" alt="x,y" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">에서 strict'라는 표현은 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x%3Dy" alt="x=y" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">일 때 양변이 같으므로 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x%5Cne%20y" alt="x\ne y" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">로 이해해야 한다.
{: .prompt-tip }

### 6.2 Testing Convexity Using the Hessian

> <strong>Theorem | Characterizations of Convexity</strong>  
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f%3A%5Cmathbb%20R%5En%5Cto%5Cmathbb%20R" alt="f:\mathbb R^n\to\mathbb R" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">가 두 번 연속 미분 가능할 때  
> - <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f" alt="f" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">가 **볼록함수** <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5CLongleftrightarrow" alt="\Longleftrightarrow" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;"> 모든 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x" alt="x" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">에서 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cnabla%5E2%20f%28x%29" alt="\nabla^2 f(x)" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">가 **양의 준정부호**.
> - 모든 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x" alt="x" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">에서 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cnabla%5E2%20f%28x%29" alt="\nabla^2 f(x)" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">가 **양의 정부호**라면 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f" alt="f" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">는 **엄격한 볼록함수**.
> - <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f" alt="f" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">가 **오목함수** <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5CLongleftrightarrow" alt="\Longleftrightarrow" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;"> 모든 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x" alt="x" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">에서 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cnabla%5E2%20f%28x%29" alt="\nabla^2 f(x)" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">가 **음의 준정부호**.
> - 모든 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x" alt="x" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">에서 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cnabla%5E2%20f%28x%29" alt="\nabla^2 f(x)" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">가 **음의 정부호**라면 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f" alt="f" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">는 **엄격한 오목함수**.  
> 위 정의와 판정은 볼록집합 및 열린 볼록집합 위의 함수에도 적절히 확장된다.  
> **Important:** '헤세 행렬이 모든 점에서 양의 정부호'는 엄격한 볼록성을 보장하는 **충분조건**이지만, **필요조건은 아니다**. 예컨대 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f%28x%29%3Dx%5E4" alt="f(x)=x^4" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">는 엄격하게 볼록하지만 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f%27%27%280%29%3D0" alt="f&#x27;&#x27;(0)=0" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">이다. 반대로 음의 정부호 조건도 엄격한 오목성의 충분조건일 뿐이다.
{: .prompt-tip }

### 6.3 Convexity and Optimality

> <strong>Theorem | Convexity and Optimality</strong>  
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20C%5Csubseteq%5Cmathbb%20R%5En" alt="C\subseteq\mathbb R^n" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">가 볼록집합이고 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f%3AC%5Cto%5Cmathbb%20R" alt="f:C\to\mathbb R" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">가 볼록함수이면 다음이 성립한다.  
> **(a)** 모든 **국소 최솟점**은 **전역 최솟점**이다.  
> **(b)** <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f" alt="f" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">가 연속적으로 1계 미분 가능하고 적절한 내부 정지점 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x%5E%2A" alt="x^*" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">에서 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cnabla%20f%28x%5E%2A%29%3D0" alt="\nabla f(x^*)=0" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">이면 그 점은 **전역 최솟점**이다.  
> **(c)** <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f" alt="f" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">가 **엄격한 볼록함수**이면 전역 최솟점은 **최대 하나**다(존재 자체를 보장하지는 않음).  
> 오목함수의 경우에도 '최솟점'을 '최댓점'으로 바꾼 유사한 정리가 성립한다.  
> 따라서 **볼록집합에서 볼록함수를 최소화하는 문제(convex optimization problem)**는 최적화 이론에서 특히 중요하다.
{: .prompt-tip }

## 7. Coercivity and the Existence of Global Minima

> <strong>Definition | Coercive Function</strong>  
> 함수 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f%3A%5Cmathbb%20R%5En%5Cto%5Cmathbb%20R" alt="f:\mathbb R^n\to\mathbb R" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">가  
>
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cboxed%7B%5Clim_%7B%5ClVert%20x%5CrVert%5Cto%5Cinfty%7Df%28x%29%3D%2B%5Cinfty%7D" alt="수식" style="display:block;max-width:100%;height:auto;max-height:16rem;margin:0.75rem auto;object-fit:contain;">
>
> 를 만족하면 **강제적(coercive)**이라고 한다. 즉 입력 벡터가 어느 방향으로든 무한히 멀어지면 함수값도 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%2B%5Cinfty" alt="+\infty" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">로 증가한다.
{: .prompt-tip }

> <strong>Theorem | Existence of a Global Minimizer</strong>  
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f%3A%5Cmathbb%20R%5En%5Cto%5Cmathbb%20R" alt="f:\mathbb R^n\to\mathbb R" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">가 **연속**이고 **강제적**이면 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f" alt="f" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">는 적어도 하나의 **전역 최솟점**을 갖는다. 실무적으로는 국소 최솟점 후보를 찾아 그 함수값을 비교해 전역 최솟점을 판정하는 데 활용할 수 있다.
{: .prompt-tip }

### 7.1 Revisiting the Three-Variable Example

> <strong>Example | A Nonconvex but Coercive Function</strong>  
> 다시  
>
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f%28x%2Cy%2Cz%29%3D%28x%5E2-z%29%5E2%2B%28y-z%29%5E2%2B2z" alt="f(x,y,z)=(x^2-z)^2+(y-z)^2+2z" style="display:block;max-width:100%;height:auto;max-height:16rem;margin:0.75rem auto;object-fit:contain;">
>
> 를 생각하자. 앞서 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%280%2C-1%2C-1%29" alt="(0,-1,-1)" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">이 엄격한 국소 최솟점임을 확인했다. 이를 전역 최솟점으로 확정하기 위해 **볼록성** 또는 **강제성**을 확인해 보자.  
> **Step 1 | Testing Convexity**  
>
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cnabla%5E2f%3D%5Cbegin%7Bpmatrix%7D12x%5E2-4z%260%26-4x%5C%5C0%262%26-2%5C%5C-4x%26-2%264%5Cend%7Bpmatrix%7D." alt="수식" style="display:block;max-width:100%;height:auto;max-height:16rem;margin:0.75rem auto;object-fit:contain;">
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%28x%2Cy%2Cz%29%3D%280%2C0%2C1%29" alt="(x,y,z)=(0,0,1)" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">에서는 첫 번째 대각성분이 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20-4%3C0" alt="-4&lt;0" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">이다. 따라서 헤세 행렬이 양의 준정부호일 수 없고, **이 함수는 전역적으로 볼록하지 않다.**  
> **Step 2 | Testing Coercivity: Three Equivalent Forms**  
>
> <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cbegin%7Baligned%7Df%26%3D%28x%5E2-z%29%5E2%2B%28y-z%29%5E2%2B2z%20%26%26%5Ctext%7B%281%29%7D%5C%5C%26%3D%28x%5E2-z%29%5E2%2B%28y-z-1%29%5E2%2B2y-1%20%26%26%5Ctext%7B%282%29%7D%5C%5C%26%3D%28x%5E2-z-1%29%5E2%2B2x%5E2-1%2B%28y-z%29%5E2%20%26%26%5Ctext%7B%283%29%7D%5Cend%7Baligned%7D" alt="수식" style="display:block;max-width:100%;height:auto;max-height:16rem;margin:0.75rem auto;object-fit:contain;">
>
> 원자료에서는 다음 상황을 나누어 모든 무한대로 가는 경로에서 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f%5Cto%2B%5Cinfty" alt="f\to+\infty" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">임을 설명한다.  
>
> 1. <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Clvert%20x%5Crvert%5Cto%5Cinfty" alt="\lvert x\rvert\to\infty" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">이고 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y%2Cz" alt="y,z" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">가 유계라면 (1)에서 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f%5Cto%2B%5Cinfty" alt="f\to+\infty" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">.
> 2. <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Clvert%20y%5Crvert%5Cto%5Cinfty" alt="\lvert y\rvert\to\infty" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">이고 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20z" alt="z" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">가 유계이거나, <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Clvert%20z%5Crvert%5Cto%5Cinfty" alt="\lvert z\rvert\to\infty" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">이고 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y" alt="y" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">가 유계라면 (1) 또는 (2)에서 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f%5Cto%2B%5Cinfty" alt="f\to+\infty" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">.
> 3. <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y%5Cto%2B%5Cinfty" alt="y\to+\infty" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;"> 또는 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20z%5Cto%2B%5Cinfty" alt="z\to+\infty" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">이면 (2) 또는 (1)에서 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f%5Cto%2B%5Cinfty" alt="f\to+\infty" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">.
> 4. <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y%5Cto-%5Cinfty" alt="y\to-\infty" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">와 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20z%5Cto-%5Cinfty" alt="z\to-\infty" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">가 함께 발생하면 (3)을 이용해 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f%5Cto%2B%5Cinfty" alt="f\to+\infty" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">임을 확인한다.  
> **결론:** 함수는 강제적이며 전역 최솟점을 갖는다. 유일한 정지점 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%280%2C-1%2C-1%29" alt="(0,-1,-1)" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">이 바로 **전역 최솟점**이다.
{: .prompt-info }

여기서는 더 간결하게 **전역 최적값을 직접 증명**할 수도 있다. (3)에서



<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f%28x%2Cy%2Cz%29%3D%28x%5E2-z-1%29%5E2%2B%28y-z%29%5E2%2B2x%5E2-1%5Cge-1." alt="f(x,y,z)=(x^2-z-1)^2+(y-z)^2+2x^2-1\ge-1." style="display:block;max-width:100%;height:auto;max-height:16rem;margin:0.75rem auto;object-fit:contain;">



등호가 성립하려면 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x%3D0" alt="x=0" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">, <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20y%3Dz" alt="y=z" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">, <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x%5E2-z-1%3D0" alt="x^2-z-1=0" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">이므로 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20z%3Dy%3D-1" alt="z=y=-1" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">이다. 따라서



<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cboxed%7B%5Coperatorname%2A%7Bargmin%7D_%7B%28x%2Cy%2Cz%29%5Cin%5Cmathbb%20R%5E3%7Df%28x%2Cy%2Cz%29%3D%5C%7B%280%2C-1%2C-1%29%5C%7D%2C%5Cqquad%5Cmin%20f%3D-1.%7D" alt="수식" style="display:block;max-width:100%;height:auto;max-height:16rem;margin:0.75rem auto;object-fit:contain;">



> **Key Takeaway | Connecting the Optimality Tests**  
> **정지조건** <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cnabla%20f%3D0" alt="\nabla f=0" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;"> → 최적해 *후보* 찾기.  
> **헤세 행렬의 정부호성** → 후보의 *국소* 최적성 분류.  
> **볼록성** → 국소 최적성을 *전역* 최적성으로 연결.  
> **강제성 + 연속성** → *전역 최솟점의 존재* 보장.  
> 때로는 제곱완성과 같은 직접적인 하한 증명이 이 모든 방법보다 간단하다.
{: .prompt-tip }

## 8. Key Vocabulary

> **Vocabulary**  

| 영어 | 한국어 | 의미 |
|:--|:--|:--|
| Optimization | 최적화 | 목적함수를 최솟값 또는 최댓값으로 만드는 과정 |
| Objective function | 목적함수 | 최적화 대상 함수 |
| Feasible set | 가능집합 | 제약을 만족하는 해의 집합 |
| Constraint | 제약조건 | 가능한 해를 제한하는 조건 |
| Critical point | 임계점 | 도함수가 0이거나 존재하지 않는 점(맥락에 따라 정의에 차이) |
| Stationary point | 정지점 | 그래디언트가 0인 점 |
| Gradient | 그래디언트, 기울기 벡터 | 1계 편도함수 벡터 |
| Hessian matrix | 헤세 행렬 | 2계 편도함수 행렬 |
| Necessary condition | 필요조건 | 최적점이 반드시 만족하는 조건 |
| Sufficient condition | 충분조건 | 만족하면 최적성 결론이 따라오는 조건 |
| Positive definite | 양의 정부호 | <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x%5ETAx%3E0" alt="x^TAx&gt;0" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;"> (<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x%5Cne0" alt="x\ne0" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">) |
| Positive semidefinite | 양의 준정부호 | <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x%5ETAx%5Cge0" alt="x^TAx\ge0" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;"> |
| Negative definite | 음의 정부호 | <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x%5ETAx%3C0" alt="x^TAx&lt;0" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;"> (<img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x%5Cne0" alt="x\ne0" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">) |
| Negative semidefinite | 음의 준정부호 | <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x%5ETAx%5Cle0" alt="x^TAx\le0" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;"> |
| Local minimizer | 국소 최솟점 | 근방에서 가장 작은 값의 점 |
| Global minimizer | 전역 최솟점 | 가능집합 전체에서 가장 작은 값의 점 |
| Saddle point | 안장점 | 극솟점·극댓점이 아닌 정지점 |
| Convex function | 볼록함수 | 선분 위 함수값이 양 끝값의 선형보간 이하 |
| Concave function | 오목함수 | 볼록함수의 반대 부등식을 만족 |
| Strictly convex | 엄격한 볼록 | 서로 다른 두 점 사이에서 엄격한 볼록부등식 |
| Coercive | 강제적 | <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5ClVert%20x%5CrVert%5Cto%5Cinfty" alt="\lVert x\rVert\to\infty" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;">일 때 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f%28x%29%5Cto%2B%5Cinfty" alt="f(x)\to+\infty" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;"> |
| Leading principal minor | 선도 주인수식 | 왼쪽 위 <img src="https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20k%5Ctimes%20k" alt="k\times k" style="display:inline-block;max-width:none;width:auto;height:1.2em;vertical-align:-0.23em;margin:0 0.12em;"> 부분행렬의 행렬식 |
| Eigenvalue criterion | 고유값 판정법 | 고유값의 부호를 통해 정부호성 판정 |
| Sylvester's criterion | 실베스터 판정법 | 선도 주인수식 부호로 엄격한 정부호성 판정 |
| Lagrange multipliers | 라그랑주 승수법 | 등식 제약 최적화 방법 |
| KKT conditions | KKT 조건 | 제약 최적화의 주요 최적성 조건 |
| Linear / Integer programming | 선형 / 정수계획법 | 선형 목적·제약 또는 정수 변수 최적화 |


---

Piotr Więcek, *Optimization Theory*, Wrocław University of Science and Technology
