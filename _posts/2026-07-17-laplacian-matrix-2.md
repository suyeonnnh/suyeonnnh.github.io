---
title: "Laplacian Matrix (2): Eigenvalues and Optimization : The Courant-Fischer Theorem"
date: 2026-07-17 22:10:00 +0900
categories: [Quant I]
tags: [Spectral Graph Theory, Eigenvalues, Rayleigh Quotient, Courant-Fischer Theorem, Spectral Theorem, Singular Value Decomposition]
math: true
toc: true
---

행렬의 고유값은 보통 다음 방정식을 만족하는 값으로 정의한다.

$$
Mv=\mu v
$$

이 정의만 보면 고유값은 특성방정식

$$
\det(\mu I-M)=0
$$

의 근에 불과한 것처럼 보인다.

그러나 대칭행렬의 고유값은 자연스러운 최적화 문제의 해로 나타난다. 가장 큰 고유값은 이차형식을 가장 크게 만드는 방향과 연결되고, 가장 작은 고유값은 이차형식을 가장 작게 만드는 방향과 연결된다.

중간에 있는 모든 고유값도 적절한 부분공간에서 수행되는 최대화와 최소화 문제로 표현할 수 있다. 이러한 관계를 정식화한 결과가 **Courant–Fischer 정리**이다.

## 2.1 Rayleigh quotient

실수 대칭행렬

$$
M\in\mathbb{R}^{n\times n}
$$

과 영벡터가 아닌 벡터 $x\in\mathbb{R}^n$을 생각하자.

### 정의 2.1.1. Rayleigh quotient

행렬 $M$에 대한 벡터 $x\neq 0$의 **Rayleigh quotient**를 다음과 같이 정의한다.

$$
\boxed{
\mathcal{R}_M(x)=
\frac{x^TMx}{x^Tx}
}
$$

분모는 벡터 $x$의 길이의 제곱이다.

$$
x^Tx=\|x\|_2^2
$$

따라서 Rayleigh quotient는 벡터 $x$의 방향에서 행렬 $M$의 이차형식 $x^TMx$가 어느 정도의 값을 가지는지를 나타낸다.

### 2.1.1 Rayleigh quotient의 크기 불변성

영이 아닌 스칼라 $\alpha$에 대하여 다음이 성립한다.

$$
\mathcal{R}_M(\alpha x) = \frac{(\alpha x)^TM(\alpha x)} {(\alpha x)^T(\alpha x)} = \frac{\alpha^2x^TMx} {\alpha^2x^Tx} = \mathcal{R}_M(x)
$$

벡터의 크기를 바꾸어도 Rayleigh quotient는 변하지 않는다. 즉 Rayleigh quotient는 벡터의 크기가 아니라 **방향**에만 의존한다.

따라서 Rayleigh quotient를 최적화할 때는 일반성을 잃지 않고

$$
\|x\|_2=1
$$

이라고 둘 수 있다. 이 경우

$$
\mathcal{R}_M(x)=x^TMx
$$

가 된다.

### 2.1.2 고유벡터의 Rayleigh quotient

$v$가 고유값 $\mu$에 대응하는 $M$의 고유벡터라고 하자.

$$
Mv=\mu v
$$

그러면

$$
\mathcal{R}_M(v) = \frac{v^TMv}{v^Tv} = \frac{v^T(\mu v)}{v^Tv} = \mu\frac{v^Tv}{v^Tv} = \mu
$$

이다.

즉 고유벡터 방향에서 계산한 Rayleigh quotient는 그 고유벡터에 대응하는 고유값과 같다.

$$
\boxed{
Mv=\mu v
\quad\Longrightarrow\quad
\mathcal{R}_M(v)=\mu
}
$$

## 2.2 Courant–Fischer Theorem

$M$이 실수 대칭행렬이고, 고유값을 큰 순서대로 배열하자.

$$
\mu_1\geq\mu_2\geq\cdots\geq\mu_n
$$

각 고유값에 대응하는 정규직교 고유벡터를

$$
v_1,v_2,\ldots,v_n
$$

이라고 하자. 즉

$$
Mv_i=\mu_i v_i
$$

이고

$$
v_i^Tv_j=
\begin{cases}
1,&i=j,\\\\
0,&i\neq j
\end{cases}
$$

이다.

### 정리 2.2.1. Courant–Fischer Theorem

실수 대칭행렬 $M$의 고유값을

$$
\mu_1\geq\mu_2\geq\cdots\geq\mu_n
$$

이라고 하자.

가장 큰 고유값과 가장 작은 고유값은 각각 Rayleigh quotient의 최댓값과 최솟값이다.

$$
\boxed{
\mu_1=
\max_{x\neq 0}
\frac{x^TMx}{x^Tx}
}
$$

$$
\boxed{
\mu_n=
\min_{x\neq 0}
\frac{x^TMx}{x^Tx}
}
$$

또한 모든 $1\leq k\leq n$에 대하여 다음이 성립한다.

$$
\boxed{
\mu_k=
\max_{\substack{
S\subseteq\mathbb{R}^n\\\\
\dim(S)=k
}}
\;
\min_{\substack{
x\in S\\\\
x\neq 0
}}
\frac{x^TMx}{x^Tx}
}
$$

그리고 같은 고유값을 다음과 같이 나타낼 수도 있다.

$$
\boxed{
\mu_k=
\min_{\substack{
T\subseteq\mathbb{R}^n\\\\
\dim(T)=n-k+1
}}
\;
\max_{\substack{
x\in T\\\\
x\neq 0
}}
\frac{x^TMx}{x^Tx}
}
$$

첫 번째 식은 **max-min 표현**, 두 번째 식은 **min-max 표현**이라고 부른다.

### 2.2.1 가장 큰 고유값과 가장 작은 고유값

Rayleigh quotient는 벡터의 크기에 영향을 받지 않으므로 다음과 같이 표현할 수 있다.

$$
\mu_1=
\max_{\|x\|_2=1}x^TMx
$$

$$
\mu_n=
\min_{\|x\|_2=1}x^TMx
$$

가장 큰 고유값에 대응하는 고유벡터 $v_1$은 단위벡터 중에서 이차형식 $x^TMx$를 가장 크게 만드는 방향이다.

$$
v_1
\in
\mathrm{arg\,max}_{\|x\|_2=1}
x^TMx
$$

가장 작은 고유값에 대응하는 고유벡터 $v_n$은 이차형식을 가장 작게 만드는 방향이다.

$$
v_n
\in
\mathrm{arg\,min}_{\|x\|_2=1}
x^TMx
$$

여기서 $\mathrm{arg\,max}$는 최댓값 자체가 아니라 최댓값을 만들어 내는 벡터를 의미한다.

### 2.2.2 $k$번째 고유값의 max-min 표현

다음 식을 살펴보자.

$$
\mu_k=
\max_{\dim(S)=k}
\;
\min_{x\in S\setminus\lbrace 0 \rbrace}
\mathcal{R}_M(x)
$$

먼저 $k$차원 부분공간 $S$를 하나 선택한다. 그 공간 안에서 Rayleigh quotient가 가장 작은 벡터를 찾는다.

그다음 이 최솟값이 가능한 한 커지도록 부분공간 $S$를 선택한다.

즉 다음과 같은 문제이다.

> $k$차원 부분공간을 하나 선택할 때, 그 공간 안의 가장 불리한 방향조차 가능한 한 큰 Rayleigh quotient를 갖도록 공간을 선택한다.

이 문제의 최적 부분공간은

$$
S^{*}=
\mathrm{span}(v_1,\ldots,v_k)
$$

이며, 이 공간에서 Rayleigh quotient의 최솟값은 $\mu_k$이다.

### 2.2.3 $k$번째 고유값의 min-max 표현

두 번째 표현은 다음과 같다.

$$
\mu_k=
\min_{\dim(T)=n-k+1}
\;
\max_{x\in T\setminus\lbrace 0 \rbrace}
\mathcal{R}_M(x)
$$

이번에는 $n-k+1$차원 부분공간 $T$를 선택한다. 그 공간 안에서 Rayleigh quotient가 가장 큰 벡터를 찾는다.

그다음 이 최댓값이 가능한 한 작아지도록 부분공간 $T$를 선택한다.

이 문제의 최적 부분공간은

$$
T^{*}=
\mathrm{span}(v_k,\ldots,v_n)
$$

이다.

### 2.2.4 고유벡터의 순차적 최적화

가장 큰 고유벡터 $v_1$을 구한 다음에는 $v_1$에 수직인 방향만 남겨 놓고 다시 이차형식을 최대화할 수 있다.

그 결과가 두 번째 고유벡터 $v_2$이다.

일반적으로 다음이 성립한다.

$$
\boxed{
v_k
\in
\mathrm{arg\,max}_{
\substack{
\|x\|_2=1\\\\
x^Tv_j=0,\;j<k
}}
x^TMx
}
$$

즉 $v_1,\ldots,v_{k-1}$ 방향을 제외하고 남은 공간에서 이차형식을 가장 크게 만드는 방향이 $v_k$이다.

반대 방향으로는 다음과 같이 표현된다.

$$
\boxed{
v_k
\in
\mathrm{arg\,min}_{
\substack{
\|x\|_2=1\\\\
x^Tv_j=0,\;j>k
}}
x^TMx
}
$$

고유값이 중복되는 경우에는 최적화 문제의 해가 하나의 벡터로 유일하게 결정되지 않을 수 있다. 이때 해당 고유값에 대응하는 고유공간 안의 여러 벡터가 모두 최적해가 된다.

## 2.3 고유벡터 기저에서의 벡터 전개

Courant–Fischer 정리의 증명은 임의의 벡터를 $M$의 고유벡터 기저로 전개하는 데서 시작한다.

$M$은 실수 대칭행렬이므로 스펙트럴 정리에 따라 정규직교 고유벡터 기저

$$
v_1,\ldots,v_n
$$

을 가진다.

따라서 모든 벡터 $x\in\mathbb{R}^n$는 다음과 같이 표현된다.

$$
\boxed{
x=\sum_{i=1}^{n}c_iv_i
}
$$

각 계수는 다음과 같다.

$$
c_i=v_i^Tx
$$

고유벡터들을 열벡터로 갖는 행렬을

$$
V=
\begin{bmatrix}
v_1&v_2&\cdots&v_n
\end{bmatrix}
$$

라고 하자.

고유벡터들이 정규직교하므로 $V$는 직교행렬이다.

$$
V^TV=VV^T=I
$$

계수벡터를

$$
c=
\begin{pmatrix}
c_1\\\\
\vdots\\\\
c_n
\end{pmatrix}
$$

라고 하면

$$
c=V^Tx
$$

이고

$$
x=Vc
$$

이다.

실제로

$$
Vc = VV^Tx = Ix = x
$$

가 성립한다.

### 보조정리 2.3.1. 고유벡터 기저에서의 이차형식

$M$이 실수 대칭행렬이고 고유값과 정규직교 고유벡터를 각각

$$
\mu_1,\ldots,\mu_n
$$

과

$$
v_1,\ldots,v_n
$$

이라고 하자.

벡터 $x$가

$$
x=\sum_{i=1}^{n}c_iv_i
$$

로 표현된다면 다음이 성립한다.

$$
\boxed{
x^TMx=
\sum_{i=1}^{n}\mu_ic_i^2
}
$$

### 증명

벡터 $x$를 고유벡터 기저로 전개한다.

$$
x=\sum_{i=1}^{n}c_iv_i
$$

따라서

$$
x^TMx = \left( \sum_{i=1}^{n}c_iv_i \right)^T M \left( \sum_{j=1}^{n}c_jv_j \right) = \left( \sum_{i=1}^{n}c_iv_i \right)^T \left( \sum_{j=1}^{n}c_jMv_j \right)
$$

각 $v_j$는 고유값 $\mu_j$에 대응하는 고유벡터이므로

$$
Mv_j=\mu_jv_j
$$

이다. 따라서

$$
x^TMx = \left( \sum_{i=1}^{n}c_iv_i \right)^T \left( \sum_{j=1}^{n}c_j\mu_jv_j \right) = \sum_{i=1}^{n} \sum_{j=1}^{n} c_ic_j\mu_jv_i^Tv_j
$$

고유벡터들이 정규직교하므로

$$
v_i^Tv_j=
\begin{cases}
1,&i=j,\\\\
0,&i\neq j
\end{cases}
$$

이다.

따라서 $i\neq j$인 모든 교차항이 사라지고 다음만 남는다.

$$
x^TMx=
\sum_{i=1}^{n}\mu_ic_i^2
$$

그러므로

$$
\boxed{
x^TMx=
\sum_{i=1}^{n}\mu_ic_i^2
}
$$

이다.

같은 방법으로 벡터의 길이는 다음과 같이 표현된다.

$$
x^Tx=
\sum_{i=1}^{n}c_i^2
$$

따라서 Rayleigh quotient는 다음과 같다.

$$
\boxed{
\mathcal{R}_M(x)=
\frac{\sum_{i=1}^{n}\mu_ic_i^2}
     {\sum_{i=1}^{n}c_i^2}
}
$$

### 2.3.1 Rayleigh quotient는 고유값의 가중평균이다

다음과 같이 가중치를 정의하자.

$$
w_i=
\frac{c_i^2}
{\sum_{j=1}^{n}c_j^2} 
(w_i\geq 0 , 
\sum_{i=1}^{n}w_i=1)
$$

따라서 Rayleigh quotient는 다음과 같이 표현된다.

$$
\mathcal{R}_M(x)=
\sum_{i=1}^{n}w_i\mu_i
$$

즉 Rayleigh quotient는 고유값들의 가중평균이다.

가중평균은 가장 작은 값보다 작아질 수 없고 가장 큰 값보다 커질 수 없으므로

$$
\boxed{
\mu_n
\leq
\mathcal{R}_M(x)
\leq
\mu_1
}
$$

이 성립한다.

이 사실만으로도 다음을 확인할 수 있다.

$$
\mu_1=
\max_{x\neq 0}
\mathcal{R}_M(x)
$$

$$
\mu_n=
\min_{x\neq 0}
\mathcal{R}_M(x)
$$

## 2.4 Courant–Fischer 정리의 증명

### 정리 2.2.1의 증명: max-min 표현

다음 식을 증명한다.

$$
\mu_k=
\max_{\dim(S)=k}
\;
\min_{x\in S\setminus\lbrace 0 \rbrace}
\mathcal{R}_M(x)
$$

증명은 두 단계로 나뉜다.

1. 어떤 $k$차원 부분공간에서는 실제로 $\mu_k$를 얻을 수 있음을 보인다.
2. 다른 어떤 $k$차원 부분공간을 선택해도 $\mu_k$보다 큰 값을 얻을 수 없음을 보인다.

#### 1단계: $\mu_k$를 달성하는 부분공간

다음 부분공간을 선택한다.

$$
S^{*}=
\mathrm{span}(v_1,\ldots,v_k)
$$

이 공간의 차원은 $k$이다.

모든 $x\in S^{*}$는 다음과 같이 표현된다.

$$
x=
\sum_{i=1}^{k}c_iv_i
$$

보조정리 2.3.1에 따라 Rayleigh quotient는

$$
\mathcal{R}_M(x)=
\frac{\sum_{i=1}^{k}\mu_ic_i^2}
     {\sum_{i=1}^{k}c_i^2}
$$

이다.

고유값은 큰 순서대로 배열되어 있으므로

$$
\mu_i\geq\mu_k,
\qquad 1\leq i\leq k
$$

이다. 따라서

$$
\mathcal{R}_M(x) = \frac{\sum_{i=1}^{k}\mu_ic_i^2} {\sum_{i=1}^{k}c_i^2} \geq \frac{\sum_{i=1}^{k}\mu_kc_i^2} {\sum_{i=1}^{k}c_i^2} = \mu_k
$$

이다.

따라서

$$
\min_{x\in S^{*}\setminus\lbrace 0 \rbrace}
\mathcal{R}_M(x)
\geq
\mu_k
$$

이다.

한편 $v_k\in S^{*}$이고

$$
\mathcal{R}_M(v_k)=\mu_k
$$

이므로 실제 최솟값은 정확히 $\mu_k$이다.

$$
\min_{x\in S^{*}\setminus\lbrace 0 \rbrace}
\mathcal{R}_M(x)=
\mu_k
$$

따라서

$$
\max_{\dim(S)=k}
\;
\min_{x\in S\setminus\lbrace 0 \rbrace}
\mathcal{R}_M(x)
\geq
\mu_k
$$

이다.

#### 2단계: 다른 부분공간에서는 $\mu_k$를 넘을 수 없음

이제 $S$를 임의의 $k$차원 부분공간이라고 하자.

다음 부분공간을 정의한다.

$$
T^{*}=
\mathrm{span}(v_k,v_{k+1},\ldots,v_n)
$$

이 공간의 차원은

$$
\dim(T^{*})=n-k+1
$$

이다.

두 공간의 차원을 더하면

$$
\dim(S)+\dim(T^{*})=
k+(n-k+1)=
n+1
$$

이다.

전체 공간 $\mathbb{R}^n$의 차원은 $n$이므로 두 부분공간은 영벡터 이외의 공통 벡터를 반드시 가진다.

실제로 부분공간의 차원 공식에 따라

$$
\dim(S\cap T^{*}) \geq \dim(S)+\dim(T^{*})-n = k+(n-k+1)-n = 1
$$

이다.

따라서 다음을 만족하는 벡터를 선택할 수 있다.

$$
0\neq x\in S\cap T^{*}
$$

$x\in T^{*}$이므로

$$
x=
\sum_{i=k}^{n}c_iv_i
$$

로 표현된다.

따라서

$$
\mathcal{R}_M(x)=
\frac{\sum_{i=k}^{n}\mu_ic_i^2}
     {\sum_{i=k}^{n}c_i^2}
$$

이다.

$i\geq k$이면

$$
\mu_i\leq\mu_k
$$

이므로

$$
\mathcal{R}_M(x)
\leq
\mu_k
$$

이다.

또한 $x\in S$이므로

$$
\min_{y\in S\setminus\lbrace 0 \rbrace}
\mathcal{R}_M(y)
\leq
\mathcal{R}_M(x)
\leq
\mu_k
$$

이다.

이는 모든 $k$차원 부분공간 $S$에 대해 성립하므로

$$
\max_{\dim(S)=k}
\;
\min_{x\in S\setminus\lbrace 0 \rbrace}
\mathcal{R}_M(x)
\leq
\mu_k
$$

이다.

앞에서 반대 방향의 부등식도 증명했으므로

$$
\boxed{
\mu_k=
\max_{\dim(S)=k}
\;
\min_{x\in S\setminus\lbrace 0 \rbrace}
\mathcal{R}_M(x)
}
$$

를 얻는다.

### 정리 2.2.1의 증명: min-max 표현

다음 식을 증명한다.

$$
\mu_k=
\min_{\dim(T)=n-k+1}
\;
\max_{x\in T\setminus\lbrace 0 \rbrace}
\mathcal{R}_M(x)
$$

먼저 다음 부분공간을 선택한다.

$$
T^{*}=
\mathrm{span}(v_k,\ldots,v_n)
$$

모든 $x\in T^{*}$는

$$
x=
\sum_{i=k}^{n}c_iv_i
$$

로 표현되므로

$$
\mathcal{R}_M(x)=
\frac{\sum_{i=k}^{n}\mu_ic_i^2}
     {\sum_{i=k}^{n}c_i^2}
\leq
\mu_k
$$

이다.

한편 $v_k\in T^{*}$이고

$$
\mathcal{R}_M(v_k)=\mu_k
$$

이므로

$$
\max_{x\in T^{*}\setminus\lbrace 0 \rbrace}
\mathcal{R}_M(x)=
\mu_k
$$

이다.

따라서

$$
\min_{\dim(T)=n-k+1}
\;
\max_{x\in T\setminus\lbrace 0 \rbrace}
\mathcal{R}_M(x)
\leq
\mu_k
$$

이다.

이제 $T$를 임의의 $n-k+1$차원 부분공간이라고 하자.

다음 공간을 정의한다.

$$
S^{*}=
\mathrm{span}(v_1,\ldots,v_k)
$$

두 공간의 차원을 더하면

$$
\dim(S^{*})+\dim(T)=
k+(n-k+1)=
n+1
$$

이므로

$$
S^{*}\cap T
$$

에는 영벡터가 아닌 벡터가 존재한다.

따라서

$$
0\neq x\in S^{*}\cap T
$$

인 $x$를 선택할 수 있다.

$x\in S^{*}$이므로

$$
x=
\sum_{i=1}^{k}c_iv_i
$$

이고

$$
\mathcal{R}_M(x)=
\frac{\sum_{i=1}^{k}\mu_ic_i^2}
     {\sum_{i=1}^{k}c_i^2}
\geq
\mu_k
$$

이다.

또한 $x\in T$이므로

$$
\max_{y\in T\setminus\lbrace 0 \rbrace}
\mathcal{R}_M(y)
\geq
\mathcal{R}_M(x)
\geq
\mu_k
$$

이다.

이는 모든 $n-k+1$차원 부분공간 $T$에 대해 성립하므로

$$
\min_{\dim(T)=n-k+1}
\;
\max_{x\in T\setminus\lbrace 0 \rbrace}
\mathcal{R}_M(x)
\geq
\mu_k
$$

이다.

따라서

$$
\boxed{
\mu_k=
\min_{\dim(T)=n-k+1}
\;
\max_{x\in T\setminus\lbrace 0 \rbrace}
\mathcal{R}_M(x)
}
$$

를 얻는다.

## 2.5 최적화를 이용한 스펙트럴 정리의 증명

Courant–Fischer 정리는 스펙트럴 정리를 이용해 증명할 수 있다.

반대로 Rayleigh quotient의 최적화 성질을 이용하면 실수 대칭행렬의 스펙트럴 정리도 증명할 수 있다.

### 정리 2.5.1. Rayleigh quotient를 최대화하는 벡터

$M$이 실수 대칭행렬이라고 하자.

영벡터가 아닌 벡터 $x$가 Rayleigh quotient

$$
\frac{x^TMx}{x^Tx}
$$

를 최대화한다면 $x$는 $M$의 가장 큰 고유값 $\mu_1$에 대응하는 고유벡터이다.

즉

$$
Mx=\mu_1x
$$

이다.

마찬가지로 Rayleigh quotient의 최솟값은 가장 작은 고유값에 대응하는 고유벡터에서 달성된다.

### 증명

Rayleigh quotient는 동차적이다.

$$
\mathcal{R}_M(\alpha x)=\mathcal{R}_M(x)
$$

따라서 Rayleigh quotient를 최대화할 때는 단위벡터만 고려해도 충분하다.

단위구면

$$
\lbrace x\in\mathbb{R}^n:\|x\|_2=1 \rbrace
$$

은 닫혀 있고 유계이므로 콤팩트하다. Rayleigh quotient는 연속함수이므로 이 집합 위에서 최댓값을 실제로 가진다.

이제 Rayleigh quotient를 다음과 같이 두자.

$$
f(x)=
\frac{x^TMx}{x^Tx}
$$

$M$이 대칭행렬이므로

$$
\nabla_x(x^TMx)=2Mx
$$

이고

$$
\nabla_x(x^Tx)=2x
$$

이다.

몫의 미분법을 적용하면

$$
\nabla_xf(x)=
\frac{
(x^Tx)(2Mx)-(x^TMx)(2x)
}{
(x^Tx)^2
}
$$

이다.

최댓값을 만드는 벡터에서는 기울기가 영벡터가 되어야 하므로

$$
(x^Tx)(2Mx)-(x^TMx)(2x)=0
$$

이다.

양변을 $2$로 나누면

$$
(x^Tx)Mx=(x^TMx)x
$$

이다.

$x\neq 0$이므로 $x^Tx>0$이다. 따라서

$$
Mx=
\frac{x^TMx}{x^Tx}x
$$

를 얻는다.

이는 $x$가 $M$의 고유벡터이고, 그 고유값이 $x$의 Rayleigh quotient라는 뜻이다.

$$
Mx=\mathcal{R}_M(x)x
$$

$x$는 Rayleigh quotient를 최대화하므로 이에 대응하는 고유값은 가장 큰 고유값 $\mu_1$이다.

따라서

$$
Mx=\mu_1x
$$

이다.

최소화하는 경우에도 같은 논리를 적용하면 가장 작은 고유값에 대응하는 고유벡터를 얻는다.

### 따름정리 2.5.2. 영행렬이 아닌 대칭행렬의 고유값

영행렬이 아닌 모든 실수 대칭행렬 $M$은 영이 아닌 고유값에 대응하는 고유벡터를 적어도 하나 가진다.

### proof

먼저 어떤 벡터 $x$에 대해서는

$$
x^TMx\neq 0
$$

임을 보인다.

어떤 $i$에 대하여 대각성분이

$$
M(i,i)\neq 0
$$

이라면 $i$번째 표준단위벡터 $e_i$를 선택한다.

그러면

$$
e_i^TMe_i=M(i,i)\neq 0
$$

이다.

이제 모든 대각성분이 영이라고 하자. $M$은 영행렬이 아니므로 어떤 $i\neq j$에 대하여

$$
M(i,j)\neq 0
$$

인 성분이 존재한다.

다음 벡터를 선택한다.

$$
x=e_i+e_j
$$

$M$은 대칭행렬이고 대각성분은 모두 영이므로

$$
x^TMx = (e_i+e_j)^TM(e_i+e_j) = M(i,j)+M(j,i) = 2M(i,j) \neq 0
$$

이다.

따라서 어떤 벡터 $x$에 대해서는 반드시

$$
x^TMx\neq 0
$$

이다.

만약

$$
x^TMx>0
$$

이라면 Rayleigh quotient의 최댓값도 양수이다. 정리 2.5.1에 따라 $M$은 양의 고유값을 가진다.

반대로

$$
x^TMx<0
$$

이라면 $-M$에 정리 2.5.1을 적용한다. 그러면 $-M$은 양의 고유값 $\nu>0$을 가진다.

이는 $M$이

$$
-\nu<0
$$

라는 영이 아닌 고유값을 가진다는 뜻이다.

따라서 영행렬이 아닌 모든 대칭행렬은 영이 아닌 고유값을 적어도 하나 가진다.

## 2.6 대칭행렬의 스펙트럴 분해

### 정리 2.6.1. 대칭행렬의 스펙트럴 분해

랭크가 $r$인 실수 대칭행렬 $M$에 대하여 영이 아닌 실수

$$
\mu_1,\ldots,\mu_r
$$

과 정규직교 벡터

$$
v_1,\ldots,v_r
$$

가 존재하여 다음이 성립한다.

$$
\boxed{
M=
\sum_{i=1}^{r}
\mu_iv_iv_i^T
}
$$

이 식의 양변에 $v_j$를 곱하면

$$
Mv_j = \sum_{i=1}^{r} \mu_iv_iv_i^Tv_j = \mu_jv_j
$$

를 얻는다.

따라서 $v_j$는 고유값 $\mu_j$에 대응하는 고유벡터이다. 이 정리는 실수 대칭행렬에 대한 스펙트럴 정리와 동등하다.

### 정리 2.6.2. 대칭행렬의 열공간과 영공간

실수 대칭행렬 $M$의 열공간은 영공간과 서로 직교한다.

즉

$$
\boxed{
\mathrm{im}(M)
\perp
\ker(M)
}
$$

이다.

### 증명

$y$가 $M$의 열공간에 속한다고 하자. 그러면 어떤 벡터 $x$가 존재하여

$$
y=Mx
$$

로 표현된다.

또한 $z$가 $M$의 영공간에 속한다고 하자.

$$
Mz=0
$$

$M$은 대칭행렬이므로

$$
M^T=M
$$

이다.

따라서

$$
z^Ty = z^TMx = (M^Tz)^Tx = (Mz)^Tx = 0^Tx = 0
$$

이다.

그러므로 열공간에 속하는 모든 벡터는 영공간에 속하는 모든 벡터와 직교한다.

$$
\mathrm{im}(M)
\perp
\ker(M)
$$

### 보조정리 2.6.3. 한 고유방향의 제거

$M$이 실수 대칭행렬이고, $v$가 영이 아닌 고유값 $\mu$에 대응하는 단위고유벡터라고 하자.

$$
Mv=\mu v,
\qquad
v^Tv=1
$$

다음 행렬을 정의한다.

$$
\widetilde{M}=
M-\mu vv^T
$$

그러면 다음이 성립한다.

1. $\ker(M)\subseteq\ker(\widetilde{M})$
2. $v\in\ker(\widetilde{M})$
3. $\mathrm{rank}(\widetilde{M})=\mathrm{rank}(M)-1$

### 증명

먼저

$$
\ker(M)\subseteq\ker(\widetilde{M})
$$

임을 보인다.

$x\in\ker(M)$이라고 하자.

$$
Mx=0
$$

고유값 $\mu$는 영이 아니므로

$$
v=\frac{1}{\mu}Mv
$$

이다. 따라서 $v$는 $M$의 열공간에 속한다.

정리 2.6.2에 따라 $M$의 열공간과 영공간은 직교하므로

$$
v^Tx=0
$$

이다.

따라서

$$
\widetilde{M}x = \left(M-\mu vv^T\right)x = Mx-\mu v(v^Tx) = 0-\mu v\cdot 0 = 0
$$

이다.

그러므로

$$
x\in\ker(\widetilde{M})
$$

이고

$$
\ker(M)\subseteq\ker(\widetilde{M})
$$

이다.

이제 $v\in\ker(\widetilde{M})$임을 보인다.

$$
\widetilde{M}v = Mv-\mu vv^Tv = \mu v-\mu v(v^Tv) = \mu v-\mu v = 0
$$

여기서 $v^Tv=1$을 사용하였다.

따라서

$$
v\in\ker(\widetilde{M})
$$

이다.

또한 $\mu\neq 0$이므로 $v\notin\ker(M)$이다. 그런데 $\widetilde{M}$에서는 $v$ 방향이 새롭게 영공간에 포함되었다.

따라서 영공간의 차원은 정확히 하나 증가한다.

$$
\dim\ker(\widetilde{M})=
\dim\ker(M)+1
$$

랭크-널리티 정리에 따라

$$
\mathrm{rank}(M)+\dim\ker(M)=n
$$

이고

$$
\mathrm{rank}(\widetilde{M})
+
\dim\ker(\widetilde{M})=
n
$$

이다.

따라서

$$
\boxed{
\mathrm{rank}(\widetilde{M})=
\mathrm{rank}(M)-1
}
$$

이다.

### 정리 2.6.1의 증명

행렬 $M$의 랭크에 대한 수학적 귀납법을 사용한다.

#### 기본 단계

$M$의 랭크가 $0$이면 $M$은 영행렬이다.

이 경우

$$
M=0
$$

이므로 정리는 자명하게 성립한다.

#### 귀납가정

랭크가 $r$인 모든 실수 대칭행렬에 대하여 다음 분해가 존재한다고 가정한다.

$$
M=
\sum_{i=1}^{r}
\mu_iv_iv_i^T
$$

#### 귀납단계

이제 랭크가 $r+1$인 실수 대칭행렬 $M$을 생각하자.

따름정리 2.5.2에 따라 $M$에는 영이 아닌 고유값 $\mu$에 대응하는 단위고유벡터 $v$가 존재한다.

$$
Mv=\mu v
$$

다음 행렬을 정의한다.

$$
\widetilde{M}=
M-\mu vv^T
$$

보조정리 2.6.3에 따라

$$
\mathrm{rank}(\widetilde{M})=r
$$

이다.

귀납가정을 $\widetilde{M}$에 적용하면 정규직교 벡터

$$
v_1,\ldots,v_r
$$

과 영이 아닌 고유값

$$
\mu_1,\ldots,\mu_r
$$

가 존재하여

$$
\widetilde{M}=
\sum_{i=1}^{r}
\mu_iv_iv_i^T
$$

가 된다.

따라서

$$
M = \widetilde{M}+\mu vv^T = \sum_{i=1}^{r} \mu_iv_iv_i^T + \mu vv^T
$$

이다.

이제

$$
v_{r+1}=v
$$

와

$$
\mu_{r+1}=\mu
$$

라고 두면

$$
\boxed{
M=
\sum_{i=1}^{r+1}
\mu_iv_iv_i^T
}
$$

를 얻는다.

또한 $v\in\ker(\widetilde{M})$이고 $v_1,\ldots,v_r$은 $\widetilde{M}$의 열공간에 속한다.

정리 2.6.2에 따라

$$
v^Tv_i=0,
\qquad
1\leq i\leq r
$$

이다.

따라서

$$
v_1,\ldots,v_r,v_{r+1}
$$

은 정규직교 벡터들이다.

이로써 랭크가 $r+1$인 경우에도 정리가 성립한다.

## 2.7 비대칭행렬과 특이값

Rayleigh quotient를 이용한 고유값의 최적화 표현은 대칭행렬에서 성립한다.

비대칭행렬에서는 고유값과 고유벡터 대신 **특이값**과 **특이벡터**를 사용하는 것이 자연스럽다.

### 정의 2.7.1. 특이값 분해

임의의 실수행렬

$$
A\in\mathbb{R}^{m\times n}
$$

의 특이값 분해는 다음과 같은 표현이다.

$$
\boxed{
A=U\Sigma V^T
}
$$

여기서 $U$의 열벡터들은 서로 정규직교하고, $V$의 열벡터들도 서로 정규직교한다. $\Sigma$는 대각성분이 음이 아닌 대각행렬이다.

$\Sigma$의 대각성분을 $A$의 **특이값**이라고 한다.

$U$의 열벡터를 왼쪽 특이벡터, $V$의 열벡터를 오른쪽 특이벡터라고 한다.

특이값을 큰 순서대로 다음과 같이 배열하자.

$$
\sigma_1\geq\sigma_2\geq\cdots\geq\sigma_r\geq 0
$$

여기서

$$
r=\min(m,n)
$$

이다.

$U$와 $V$의 열벡터를 각각

$$
u_1,\ldots,u_r
$$

과

$$
v_1,\ldots,v_r
$$

이라고 하면 특이값 분해는 다음과 같이 표현된다.

$$
A=
\sum_{i=1}^{r}
\sigma_iu_iv_i^T
$$

또한

$$
Av_i=\sigma_i u_i
$$

이고

$$
A^Tu_i=\sigma_i v_i
$$

이다.

따라서

$$
A^TAv_i=
\sigma_i^2v_i
$$

이고

$$
AA^Tu_i=
\sigma_i^2u_i
$$

이다.

즉 오른쪽 특이벡터는 $A^TA$의 고유벡터이고, 왼쪽 특이벡터는 $AA^T$의 고유벡터이다.

특이값의 제곱은 두 행렬의 고유값이다.

$$
\boxed{
\sigma_i^2=
\lambda_i(A^TA)=
\lambda_i(AA^T)
}
$$

### 정리 2.7.2. 특이값의 최적화 표현

임의의 실수행렬

$$
A\in\mathbb{R}^{m\times n}
$$

의 특이값을

$$
\sigma_1\geq\sigma_2\geq\cdots\geq\sigma_r\geq 0
$$

이라고 하자.

가장 큰 특이값은 다음과 같이 표현된다.

$$
\boxed{
\sigma_1=
\max_{\|v\|_2=1}
\|Av\|_2
}
$$

또한

$$
\boxed{
\sigma_1=
\max_{\substack{
\|u\|_2=1\\\\
\|v\|_2=1
}}
u^TAv
}
$$

이다.

일반적인 $k$번째 특이값은 다음과 같이 표현된다.

$$
\boxed{
\sigma_k=
\max_{\substack{
S\subseteq\mathbb{R}^n\\\\
\dim(S)=k
}}
\;
\min_{\substack{
v\in S\\\\
\|v\|_2=1
}}
\|Av\|_2
}
$$

또한

$$
\boxed{
\sigma_k=
\min_{\substack{
T\subseteq\mathbb{R}^n\\\\
\dim(T)=n-k+1
}}
\;
\max_{\substack{
v\in T\\\\
\|v\|_2=1
}}
\|Av\|_2
}
$$

이다.

### 증명

먼저 다음 관계를 사용한다.

$$
\|Av\|_2^2=
(Av)^T(Av)=
v^TA^TAv
$$

행렬 $A^TA$는 대칭행렬이며 양의 준정부호이다.

$A^TA$의 고유값은 특이값의 제곱이다.

$$
\sigma_1^2
\geq
\sigma_2^2
\geq
\cdots
$$

Courant–Fischer 정리를 $A^TA$에 적용하면

$$
\sigma_k^2=
\max_{\dim(S)=k}
\;
\min_{\substack{
v\in S\\\\
v\neq 0
}}
\frac{v^TA^TAv}{v^Tv}
$$

이다.

그런데

$$
\frac{v^TA^TAv}{v^Tv}=
\frac{\|Av\|_2^2}{\|v\|_2^2}
$$

이므로

$$
\sigma_k^2=
\max_{\dim(S)=k}
\;
\min_{\substack{
v\in S\\\\
\|v\|_2=1
}}
\|Av\|_2^2
$$

이다.

양변에 제곱근을 취하면

$$
\boxed{
\sigma_k=
\max_{\dim(S)=k}
\;
\min_{\substack{
v\in S\\\\
\|v\|_2=1
}}
\|Av\|_2
}
$$

를 얻는다.

min-max 표현도 $A^TA$에 Courant–Fischer 정리를 적용하여 같은 방식으로 얻는다.

이제 가장 큰 특이값의 쌍선형 표현을 확인하자.

고정된 단위벡터 $v$에 대하여 코시–슈바르츠 부등식에 의해

$$
u^TAv
\leq
\|u\|_2\|Av\|_2=
\|Av\|_2
$$

이다.

$Av\neq 0$일 때

$$
u=\frac{Av}{\|Av\|_2}
$$

로 두면 등호가 성립한다.

따라서

$$
\max_{\|u\|_2=1}u^TAv=
\|Av\|_2
$$

이다.

이제 $v$에 대해서도 최대화하면

$$
\max_{\substack{ \|u\|_2=1\\\\ \|v\|_2=1 }} u^TAv = \max_{\|v\|_2=1} \|Av\|_2 = \sigma_1
$$

를 얻는다.
