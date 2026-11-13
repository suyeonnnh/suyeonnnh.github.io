---
title: "Optimization Theory : Gradient Methods 2"
date: 2026-11-08 00:00:00 +0100
categories:
  - Optimization Theory
tags:
  - Gradient Methods
  - Quasi-Newton
  - DFP
  - BFGS
  - Conjugate Gradient
  - Zeroth-Order Methods
math: true
---
>
> 이번 파트에서는 steepest descent와 Newton method의 비교, quasi-Newton methods, DFP와 BFGS, conjugate gradient method, Fletcher–Reeves method, zeroth-order methods를 다룬다.
{: .prompt-info }

---

## 1. Comparison of the Steepest Descent and Newton Methods (p.3)

Newton method는 stationary point 근처에서 steepest descent보다 훨씬 빠르게 수렴한다. 반면 steepest descent는 수렴에 필요한 조건이 더 약하고 구현이 단순하다.

Newton method의 기본 갱신식은 다음과 같다.

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x_%7Bk&plus;1%7D%3Dx_k-%5Cleft%28%5Cnabla%5E2f%28x_k%29%5Cright%29%5E%7B-1%7D%5Cnabla%20f%28x_k%29)

하지만 이 방법은 Hessian이 필요하고, Hessian의 역행렬 계산도 필요하다. 이계도함수를 구할 수 없거나, 행렬 차원이 크면 계산량과 수치오차가 커질 수 있다.

이 문제를 보완하기 위해 강의에서는 다음 두 계열의 방법을 소개한다.

- **Quasi-Newton methods**
- **Conjugate directions methods**

이들의 목적은 Newton method보다 계산은 단순하게 유지하면서도, steepest descent보다 빠른 수렴을 얻는 것이다.

## 2. Types of Alternative Methods (p.4)

### 2.1 Quasi-Newton Methods

Quasi-Newton methods에서는 Hessian의 역행렬을 직접 계산하지 않고, 그것을 근사하는 행렬 \( D_k \)를 반복적으로 갱신한다.

탐색 방향은 다음과 같이 잡는다.

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20d_k=-D_k%5Cnabla%20f%28x_k%29)

그리고 보정행렬 \( C_k \)를 이용해 \( D_k \)를 갱신한다.

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20D_%7Bk&plus;1%7D=D_k&plus;C_k)

갱신식은 \( D_k \), \( x_{k+1}-x_k \), \( \nabla f(x_{k+1})-\nabla f(x_k) \)를 사용해 구성된다. 또한 매 단계에서 \( D_k \)가 positive definite이면 \( D_{k+1} \)도 positive definite가 되도록 설계한다.

### 2.2 Conjugate Directions Methods

이 방법들은 현재 gradient뿐 아니라 이전 탐색 방향도 함께 사용한다.  
이차함수에서는 탐색 방향들이 \( Q \)-conjugate가 되도록 구성한다.

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20d_k%5ETQd_l=0%5Cqquad%28k%5Cneq%20l%29)

이는 일반적인 직교 조건과 다르며, \( Q \)가 정의하는 내적에 대한 직교성이다.

## 3. General Quasi-Newton Algorithm (p.5)

Quasi-Newton method의 일반적인 구조는 다음과 같다.

1. 초기화: \( D:=I \), \( x:=x_0 \)
2. 현재 gradient를 계산한다.
3. 탐색 방향 \( d=-Dg \)를 정한다.
4. line search로 \( f(x+\alpha d) \)를 최소화하는 \( \alpha^\* \)를 찾는다.
5. \( p=\alpha^\* d \), \( x:=x+p \)로 갱신한다.
6. gradient 변화 \( q \)를 계산한다.
7. \( D \)와 \( p,q \)를 이용해 correction matrix \( C \)를 구한다.
8. \( D:=D+C \)로 갱신한다.

핵심 변수는 다음과 같다.

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20p_k%3D%5Calpha_kd_k%3Dx_%7Bk&plus;1%7D-x_k)

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20q_k%3D%5Cnabla%20f%28x_%7Bk&plus;1%7D%29-%5Cnabla%20f%28x_k%29)

Quasi-Newton methods의 차이는 전부 \( C_k \)를 어떻게 정의하느냐에서 나온다.

## 4. Requirements for the Correction Matrix (p.6)

Correction matrix \( C_k \)는 다음 조건을 만족하는 것이 바람직하다.

1. \( D_k \)가 positive definite이면 \( D_k+C_k \)도 positive definite이어야 한다.
2. \( f \)가 convex quadratic function이면, 적절한 step size 하에서 \( n \)번의 반복 후 \( D_k \)가 정확히 \( Q^{-1} \)이 되도록 해야 한다.
3. 이 경우 탐색 방향들은 \( Q \)-conjugate 집합을 이루어야 한다.

즉, quadratic case에서는 Newton method와 비슷한 유한 단계 성질을 얻는 것이 목표다.

## 5. Davidon–Fletcher–Powell (DFP) Method (pp.7–9)

DFP method는 secant equation에서 출발한다.  
Gradient의 Taylor expansion을 이용하면

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cnabla%20f%28x_%7Bk&plus;1%7D%29%5Capprox%20%5Cnabla%20f%28x_k%29&plus;%5Cnabla%5E2f%28x_k%29%28x_%7Bk&plus;1%7D-x_k%29)

이고, 이를 정리하면

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cleft%28%5Cnabla%5E2f%28x_k%29%5Cright%29%5E%7B-1%7Dq_k%5Capprox%20p_k)

를 얻는다. 여기서

- \( p_k=x_{k+1}-x_k \)
- \( q_k=\nabla f(x_{k+1})-\nabla f(x_k) \)

이다.

Inverse Hessian을 \( D_{k+1} \)로 근사하면 secant equation은

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20D_%7Bk&plus;1%7Dq_k=p_k)

가 된다.

DFP correction matrix는 다음과 같다.

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20C_k%3D%5Cfrac%7Bp_kp_k%5ET%7D%7Bp_k%5ETq_k%7D-%5Cfrac%7BD_kq_kq_k%5ETD_k%7D%7Bq_k%5ETD_kq_k%7D)

따라서 DFP update는

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20D_%7Bk&plus;1%7D%3DD_k&plus;%5Cfrac%7Bp_kp_k%5ET%7D%7Bp_k%5ETq_k%7D-%5Cfrac%7BD_kq_kq_k%5ETD_k%7D%7Bq_k%5ETD_kq_k%7D)

이다.

## 6. Broyden–Fletcher–Goldfarb–Shanno (BFGS) Method (pp.7, 10–11)

BFGS는 Hessian 근사행렬 \( B_k \)를 갱신한 뒤, 이를 inverse form으로 바꿔 사용하는 방법이다.

Hessian approximation update는 다음과 같다.

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20B_%7Bk&plus;1%7D%3DB_k&plus;%5Cfrac%7Bq_kq_k%5ET%7D%7Bq_k%5ETp_k%7D-%5Cfrac%7BB_kp_kp_k%5ETB_k%7D%7Bp_k%5ETB_kp_k%7D)

이로부터 inverse-Hessian form \( D_k=B_k^{-1} \)를 얻을 때 Sherman–Morrison formula를 사용한다.

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%28A&plus;uv%5ET%29%5E%7B-1%7D=A%5E%7B-1%7D-%5Cfrac%7BA%5E%7B-1%7Duv%5ETA%5E%7B-1%7D%7D%7B1&plus;v%5ETA%5E%7B-1%7Du%7D)

BFGS inverse update는 보통 다음처럼 쓴다.

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20D_%7Bk&plus;1%7D%3D%28I-%5Crho_kp_kq_k%5ET%29D_k%28I-%5Crho_kq_kp_k%5ET%29&plus;%5Crho_kp_kp_k%5ET%2C%5Cqquad%20%5Crho_k%3D%5Cfrac1%7Bq_k%5ETp_k%7D)

강의에서는 DFP와 BFGS 모두 convex function과 적절한 step size rule 하에서 positive definiteness를 유지한다고 설명한다.

## 7. Other Quasi-Newton Methods (p.12)

### 7.1 Hoshino Method

Hoshino method는 \( p_k, q_k, D_k \)를 이용해 다른 형태의 correction matrix를 만든다.

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20C_k%3D%5CTheta%20p_kp_k%5ET-%5CPsi%28p_kq_k%5ETD_k&plus;D_kq_kp_k%5ET&plus;D_kq_kq_k%5ETD_k%29)

여기서 \( \Theta \)와 \( \Psi \)는 \( p_k,q_k,D_k \)를 이용해 정해지는 스칼라다.

### 7.2 Broyden Method

Broyden family는 DFP와 BFGS를 매개변수 \( \Phi \)로 연결한 형태다.

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20C_k%5E%7BB%7D=%281-%5CPhi%29C_k%5E%7BDFP%7D&plus;%5CPhi%20C_k%5E%7BBFGS%7D%2C%5Cqquad%20%5CPhi%5Cin%5B0%2C1%5D)

- \( \Phi=0 \): DFP
- \( \Phi=1 \): BFGS

### 7.3 Fletcher (Switch) Method

Fletcher switch method는 DFP와 BFGS 중 어느 쪽이 inverse Hessian 방향을 더 잘 근사하는지 검사하여, \( \Phi \)를 0 또는 1로 바꾸는 방식이다.

## 8. Conjugate Gradient Method (pp.13–16)

### 8.1 Quadratic Problem

Convex quadratic function

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20f%28x%29%3D%5Cfrac12x%5ETQx&plus;b%5ETx)

를 최소화한다고 하자. 여기서 \( Q \)는 symmetric positive definite matrix다.

Gradient는

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cnabla%20f%28x%29%3DQx&plus;b)

이다.

강의에서는 negative gradient를

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20g_k=-%28Qx_k&plus;b%29)

로 두고, 초기 방향을 \( d_0=g_0 \)로 잡는다.

### 8.2 Step Size

정확한 minimization rule을 쓰면 step size는

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Calpha_k%3D%5Cfrac%7Bg_k%5ETg_k%7D%7Bd_k%5ETQd_k%7D)

가 된다.

그 후 새 점은

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x_%7Bk&plus;1%7D=x_k&plus;%5Calpha_kd_k)

로 갱신한다.

### 8.3 Direction Update

탐색 방향은 이전 방향과 현재 negative gradient를 이용해

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20d_k=g_k&plus;%5Cbeta_kd_%7Bk-1%7D)

형태로 만든다.  
여기서 계수는

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cbeta_k%3D%5Cfrac%7Bg_k%5ETg_k%7D%7Bg_%7Bk-1%7D%5ETg_%7Bk-1%7D%7D)

이다.

이 방법에서는 탐색 방향들이 \( Q \)-conjugate가 되며,

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20d_k%5ETQd_l=0%5Cqquad%28k%5Cneq%20l%29)

정확한 계산을 가정하면 최대 \( n \)번 안에 최소점에 도달한다.

### 8.4 Orthogonality of Gradients

Conjugate gradient의 유도에서 중요한 관계는

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20g_k-g_%7Bk-1%7D=-%5Calpha_%7Bk-1%7DQd_%7Bk-1%7D)

이다.

또한 정확한 line search를 사용하면 current gradient는 이전 탐색 방향들과 직교한다.

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20g_k%5ETd_i=0%5Cqquad%28i=0%2C1%2C%5Cdots%2Ck-1%29)

이 성질을 이용해 \( \beta_k \) 공식을 단순화할 수 있다.

## 9. Fletcher–Reeves Method (p.17)

Fletcher–Reeves method는 conjugate gradient를 일반 함수로 확장한 것이다.

기본 구조는 다음과 같다.

- 한 묶음에서 \( n \)단계까지 conjugate gradient처럼 진행한다.
- 그 다음에는 새 점에서 다시 steepest descent 방향으로 시작한다.
- 이를 반복하면서 수렴을 노린다.

공식은 다음과 같다.

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20g_k=-%5Cnabla%20f%28x_k%29)

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cbeta_k%3D%5Cfrac%7Bg_k%5ETg_k%7D%7Bg_%7Bk-1%7D%5ETg_%7Bk-1%7D%7D)

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20d_k=g_k&plus;%5Cbeta_kd_%7Bk-1%7D)

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x_%7Bk&plus;1%7D=x_k&plus;%5Calpha_kd_k)

각 단계의 \( \alpha_k \)는 line search로 정한다.

## 10. Other Conjugate Directions Methods (p.18)

### 10.1 Powell’s Method

Powell’s method는 처음에 각 좌표축 방향으로 차례대로 line search를 수행한다.  
그 결과 얻은 점 \( x_n \)을 사용해 새로운 방향

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20d_%7Bn&plus;1%7D=x_n-x_0)

를 만든다.

이후에는 오래된 방향 하나를 제거하고, 새 방향을 추가하면서 같은 절차를 반복한다.

### 10.2 Partan Method

Partan method는 처음 두 번은 steepest descent를 쓰고, 그 결과를 바탕으로 새로운 방향을 만든다.

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20d_3=x_2-x_0)

그 뒤에는 최근 두 step 정보를 사용해 방향을

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20d_k=x_%7Bk-1%7D-x_%7Bk-3%7D)

로 잡고, step size는 매번 최적으로 정한다.

## 11. Zeroth-Order Methods (pp.20–22)

### 11.1 General Idea

이 방법들은 gradient나 Hessian을 모를 때, 함수값만 이용해 최소화를 수행한다.  
Steepest descent의 gradient를 추정량으로 바꾼다.

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20x_%7Bk&plus;1%7D%3Ax_k-%5Calpha_kG_f%28x_k%2Ch%29)

즉, \( G_f(x_k,h) \)는 gradient estimator 역할을 한다.

### 11.2 Finite-Difference Estimator

가장 단순한 gradient estimator는 각 좌표축 방향의 차분을 이용하는 것이다.

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20G_f%28x_k%2Ch%29%3D%5Csum_%7Bi=1%7D%5En%5Cfrac%7Bf%28x_k&plus;he_i%29-f%28x_k-he_i%29%7D%7B2h%7De_i)

여기서 \( e_i \)는 표준기저벡터다.

### 11.3 Random-Direction Estimator

더 적은 수의 방향만 써서 추정할 수도 있다. \( z_i \)를 random direction이라고 하면,

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20G_f%28x_k%2Ch%29%3D%5Cfrac%7Bn%7Dm%5Csum_%7Bi=1%7D%5Em%5Cfrac%7Bf%28x_k&plus;hz_i%29-f%28x_k-hz_i%29%7D%7B2h%7Dz_i)

이다.

\( z_i \)는 i.i.d.인 \( n \)-차원 확률벡터이며, 구면대칭 분포를 따르고

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5Cmathbb%7BE%7D%5B%5C%7Cz_i%5C%7C%5E2%5D=1)

을 만족한다고 가정한다.

강의에서는 실제로 \( m<n \)인 경우를 고려한다.

### 11.4 Convergence

강의자료에서는 zeroth-order methods에 대해 다음을 언급한다.

- uniform distribution을 사용하면 trajectory는 steepest descent와 비슷하지만 fluctuation이 있다.
- \( f \)가 smooth하고 strongly convex이면 선형수렴(linear convergence)을 얻는다.
- 차원 \( n \)이 커질수록 수렴비는 1에 가까워진다.
- convexity 가정이 없어도 gradient norm에 대한 점근적 결과를 얻을 수 있다.


![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%20%5C%7C%5Cnabla%20f%28x_k%29%5C%7C%5E2%5Csim%20O%28k%5E%7B-2%7D%29)

## 12. Recapitulation (p.19)

강의에서는 large-scale optimization에서 널리 쓰이는 방법으로 다음을 강조한다.

- Fletcher–Reeves conjugate gradient method
- Davidon–Fletcher–Powell (DFP)
- Broyden–Fletcher–Goldfarb–Shanno (BFGS)

이 방법들은 Hessian이나 inverse Hessian을 직접 계산하지 않아도 되고, Newton method보다 반복 횟수는 더 많을 수 있지만 계산시간 측면에서는 더 유리할 수 있다.

---

**Reference**  
Piotr Więcek, *Optimization Theory — Gradient Methods: Part II*, Wrocław University of Science and Technology, 30.10.2026.
