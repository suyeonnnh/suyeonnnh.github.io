---
layout: post
title: "Black-Scholes-Merton PDE 유도 2"
date: 2026-07-06 09:00:00 +0900
categories: [Finance Mathematics]
tags: [Black-Scholes-Merton, PDE, Ito Lemma, Replicating Portfolio]
math: true
permalink: /posts/bsm-pde-replicating-portfolio/
---
## Black-Scholes-Merton PDE 유도 2 

이번 절에서는 Black-Scholes-Merton 모형을 이용해 옵션 가격이 만족해야 하는 PDE를 다시 유도한다.

핵심 아이디어는 간단하다. 같은 옵션 가격 변화 $df(t,S_t)$를 두 가지 방식으로 계산한 뒤, 두 식을 서로 비교하는 것이다.

첫 번째 방식은 <span style="color:red;">복제 포트폴리오</span>와 <span style="color:red;">자기금융 조건</span>을 이용하는 것이다. 이 방법에서는 옵션을 주식과 무위험 자산으로 복제할 수 있다고 보고, 복제 포트폴리오의 가치 변화로부터 $df(t,S_t)$를 구한다.

두 번째 방식은 옵션 가격 함수 $f(t,S_t)$에 <span style="color:red;">Ito 공식</span>을 직접 적용하는 것이다. 이 방법에서는 주식 가격 과정 $dS_t$를 이용해 옵션 가격 변화 $df(t,S_t)$를 계산한다.

두 방법은 모두 같은 옵션 가격 변화를 설명해야 하므로 결과가 서로 같아야 한다. 특히 확률항인 $dB_t$ 부분은 일치하고, drift 부분을 비교하면 Black-Scholes-Merton PDE가 도출된다.

앞 절에서는 옵션과 주식을 조합한 포트폴리오에서 확률항을 제거하는 방식으로 PDE를 유도했다면, 이번 절에서는 복제 포트폴리오와 자기금융 조건을 이용해 더 일반적인 방식으로 접근한다.

---

### 1. 모형의 기본 설정

Black-Scholes-Merton 모형에서는 두 개의 거래 가능한 자산을 생각한다.

하나는 무위험 자산인 money market account이고, 다른 하나는 주식이다.

무위험 자산의 가격을 $M_t$, 주식 가격을 $S_t$라고 하면 각각 다음의 SDE를 따른다.

$$
dM_t = rM_tdt
$$

$$
dS_t = S_t(\mu dt + \sigma dB_t)
$$

여기서 각 기호의 의미는 다음과 같다.

| 기호 | 의미 |
|---|---|
| $M_t$ | 무위험 자산의 가격 |
| $S_t$ | 주식 가격 |
| $r$ | 무위험 이자율 |
| $\mu$ | 주식의 기대수익률 또는 drift |
| $\sigma$ | 주식의 변동성 |
| $B_t$ | 브라운 운동 |

무위험 자산은 확률항이 없으므로 확정적으로 증가한다.  
초기값을 편의상 $M_0 = 1$로 두면,

$$
M_t = e^{rt}
$$

가 된다.

반면 주식 가격 $S_t$는 기하 브라운 운동(Geometric Brownian Motion)을 따른다. 따라서 Ito 공식에 의해 다음과 같이 쓸 수 있다.

$$
S_t = S_0 \exp \left( \left(\mu - \frac{\sigma^2}{2}\right)t + \sigma B_t \right)
$$

즉, 주식 가격은 확률항 $B_t$의 영향을 받기 때문에 미래 가격이 확정되어 있지 않다.

---

### 2. 파생상품의 만기 보상

이제 만기 $T$에서 하나의 payoff를 갖는 파생상품을 생각하자.

이 파생상품의 만기 payoff를 다음과 같이 둔다.

$$
g(S_T)
$$

예를 들어 European call option의 경우 만기 payoff는 다음과 같다.

$$
g(S_T) = (S_T - K)^+
$$

즉,

$$
(S_T - K)^+ = \max(S_T - K, 0)
$$

이다.

콜옵션은 만기 시점의 주가가 행사가격보다 높으면 차익을 얻고, 그렇지 않으면 행사하지 않으므로 payoff가 0이 된다.

---

### 3. 복제 포트폴리오의 설정

Black-Scholes-Merton 모형에서는 옵션의 payoff를 복제할 수 있는 포트폴리오가 존재한다고 가정한다.

즉, 주식과 무위험 자산을 적절히 조합해서 옵션의 payoff와 똑같은 현금흐름을 만들 수 있다고 보는 것이다.

이때 복제 포트폴리오의 가치를 $V_t$라고 하면,

$$
V_t = \theta_t^{(S)} S_t + \theta_t^{(M)} M_t
$$

로 쓸 수 있다.

여기서 각 기호의 의미는 다음과 같다.

| 기호 | 의미 |
|---|---|
| $\theta_t^{(S)}$ | 시점 $t$에서 보유한 주식의 수 |
| $\theta_t^{(M)}$ | 시점 $t$에서 보유한 무위험 자산의 수 |
| $V_t$ | 복제 포트폴리오의 가치 |

만기에는 이 복제 포트폴리오가 옵션 payoff와 같아야 한다.

$$
V_T = g(S_T)
$$

즉,

$$
g(S_T) = \theta_T^{(S)}S_T + \theta_T^{(M)}M_T
$$

이다.

---

### 4. 옵션 가격 함수와 무차익 조건

옵션의 가격은 시간 $t$와 주식 가격 $S$의 함수라고 가정한다.

$$
f(t, S)
$$

따라서 실제 주가 과정 $S_t$ 위에서 옵션 가격은

$$
f(t, S_t)
$$

로 쓸 수 있다.

무차익 원리에 의해 옵션 가격은 복제 포트폴리오의 가치와 같아야 한다.

$$
f(t, S_t) = V_t
$$

왜냐하면 두 자산이 만기 payoff가 같은데 현재 가격이 다르다면, 싼 것을 사고 비싼 것을 파는 방식으로 무위험 차익거래가 가능해지기 때문이다.

따라서

$$
df(t, S_t) = dV_t
$$

가 되어야 한다.

---

### 5. 자기금융 조건 적용

복제 포트폴리오가 자기금융적이라는 것은 외부에서 돈을 추가로 넣거나 빼지 않는다는 뜻이다.

즉, 포트폴리오 가치 변화는 오직 보유 중인 자산 가격 변화에서만 발생한다.

따라서

$$
dV_t = \theta_t^{(S)} dS_t + \theta_t^{(M)} dM_t
$$

이다.

여기에

$$
dS_t = S_t(\mu dt + \sigma dB_t)
$$

$$
dM_t = rM_tdt
$$

를 대입하면,

$$
\color{#2563eb}{dV_t = \theta_t^{(S)} S_t(\mu dt + \sigma dB_t) + \theta_t^{(M)} rM_tdt}
$$

이다.

이를 정리하면,

$$
\color{#2563eb}{dV_t = \left(\theta_t^{(S)} S_t \mu + \theta_t^{(M)} M_t r\right)dt + \theta_t^{(S)} S_t \sigma dB_t}
$$

이다.

무차익 조건에 의해 옵션 가격의 변화와 복제 포트폴리오의 가치 변화가 같아야 하므로,

$$
df(t,S_t)=dV_t
$$

가 성립한다.

따라서,

$$
\color{#2563eb}{df(t,S_t) = \left(\theta_t^{(S)} S_t \mu + \theta_t^{(M)} M_t r\right)dt + \theta_t^{(S)} S_t \sigma dB_t}
$$

가 된다.

---

### 6. 델타와 복제 포트폴리오

옵션 가격 함수 $f(t, S)$를 주가 $S$에 대해 미분한 값을 델타(delta)라고 한다.

$$
\Delta = f_S(t, S_t)
$$

여기서

$$
f_S(t, S_t) = \frac{\partial f}{\partial S}(t, S_t)
$$

이다.

복제 포트폴리오에서 보유해야 하는 주식 수는 옵션의 델타와 같다.

$$
\theta_t^{(S)} = f_S(t, S_t)
$$

이제 무차익 조건

$$
f(t, S_t) = V_t
$$

와

$$
V_t = \theta_t^{(S)}S_t + \theta_t^{(M)}M_t
$$

를 이용하면,

$$
f(t, S_t)=\theta_t^{(S)}S_t + \theta_t^{(M)}M_t
$$

이다.

여기에 $\theta_t^{(S)} = f_S(t, S_t)$를 대입하면,

$$
f(t, S_t)=f_S(t, S_t)S_t + \theta_t^{(M)}M_t
$$

따라서 무위험 자산에 들어간 금액은

$$
\theta_t^{(M)}M_t=f(t, S_t) - f_S(t, S_t)S_t
$$

가 된다.

---

### 7. 복제 포트폴리오에서 얻은 가격 변화식

앞에서 구한 식은 다음과 같았다.

$$
df(t, S_t)=\left(\theta_t^{(S)} S_t \mu + \theta_t^{(M)} M_t r\right)dt + \theta_t^{(S)} S_t \sigma dB_t
$$

여기에

$$
\theta_t^{(S)} = f_S(t, S_t)
$$

와

$$
\theta_t^{(M)}M_t=f(t, S_t) - f_S(t, S_t)S_t
$$

를 대입한다.

그러면 drift 부분은 다음과 같이 바뀐다.

$$
\color{#2563eb}{\theta_t^{(S)} S_t \mu + \theta_t^{(M)} M_t r = f_S(t, S_t)S_t\mu + \left(f(t, S_t) - f_S(t, S_t)S_t\right)r}
$$

정리하면,

$$
\color{#2563eb}{f_S(t, S_t)S_t\mu + \left(f(t, S_t) - f_S(t, S_t)S_t\right)r = rf(t, S_t) + (\mu - r)S_t f_S(t, S_t)}
$$

이다.

따라서 복제 포트폴리오 관점에서 옵션 가격 변화는 다음과 같다.

$$
\color{#2563eb}{df(t, S_t)=\left[rf(t, S_t) + (\mu - r)S_t f_S(t, S_t)\right]dt + S_t\sigma f_S(t, S_t)dB_t}
$$

이 식은 자기금융 복제 포트폴리오를 이용해서 얻은 옵션 가격 변화식이다.

---

### 8. Ito 공식으로 얻은 가격 변화식

이번에는 같은 $df(t,S_t)$를 Ito 공식으로 구해보자.

일반적으로 어떤 확률과정 $X_t$가 다음과 같이 주어진다고 하자.

$$
dX_t = \mu_t dt + \sigma_t dB_t
$$

이때 시간 $t$와 확률과정 $X_t$에 의존하는 함수 $f(t,X_t)$에 Ito 공식을 적용하면 다음과 같다.

$$
df(t,X_t) = \left(f_t(t,X_t) + \mu_t f_x(t,X_t) + \frac{1}{2}\sigma_t^2 f_{xx}(t,X_t)\right)dt + \sigma_t f_x(t,X_t)dB_t
$$

이제 이 일반 공식을 주식 가격 과정에 적용한다.

BSM 모형에서 주식 가격은 다음 과정을 따른다.

$$
dS_t = S_t(\mu dt+\sigma dB_t)
$$

즉,

$$
dS_t = \mu S_tdt+\sigma S_tdB_t
$$

이다.

따라서 일반 Ito 공식의 입력 과정 $X_t$를 주식 가격 $S_t$로 보면 다음과 같다.

$$
\color{#2563eb}{X_t = S_t}
$$

또한 일반식의 drift 항과 diffusion 항은 각각 다음과 같이 대응된다.

$$
\color{#2563eb}{\mu_t = \mu S_t}
$$

$$
\color{#2563eb}{\sigma_t = \sigma S_t}
$$

즉, 대응 관계를 정리하면 다음과 같다.

| 일반 Ito 공식 | 주식 가격 과정에 적용 |
|---|---|
| $X_t$ | $S_t$ |
| $\mu_t$ | $\mu S_t$ |
| $\sigma_t$ | $\sigma S_t$ |
| $f_x(t,X_t)$ | $f_S(t,S_t)$ |
| $f_{xx}(t,X_t)$ | $f_{SS}(t,S_t)$ |

따라서 일반 Ito 공식

$$
df(t,X_t) = \left(f_t(t,X_t) + \mu_t f_x(t,X_t) + \frac{1}{2}\sigma_t^2 f_{xx}(t,X_t)\right)dt + \sigma_t f_x(t,X_t)dB_t
$$

에

$$
X_t=S_t,\qquad \mu_t=\mu S_t,\qquad \sigma_t=\sigma S_t
$$

를 대입하면 다음과 같다.

$$
df(t,S_t) = \left(f_t(t,S_t) + \mu S_t f_S(t,S_t) + \frac{1}{2}(\sigma S_t)^2 f_{SS}(t,S_t)\right)dt + \sigma S_t f_S(t,S_t)dB_t
$$

이를 정리하면 다음 식을 얻는다.

$$
\color{#2563eb}{df(t,S_t) = \left(f_t(t,S_t) + \mu S_t f_S(t,S_t) + \frac{1}{2}\sigma^2S_t^2 f_{SS}(t,S_t)\right)dt + \sigma S_t f_S(t,S_t)dB_t}
$$

여기서 각 미분항의 의미는 다음과 같다.

| 기호 | 의미 |
|---|---|
| $f_t$ | 시간에 대한 옵션 가격의 변화 |
| $f_S$ | 주가에 대한 옵션 가격의 1차 미분 ; Delta |
| $f_{SS}$ | 주가에 대한 옵션 가격의 2차 미분 ; Gamma |
| $\frac{1}{2}\sigma^2S_t^2f_{SS}$ | Ito 공식에서 추가되는 2차 변동성 항 |

여기서 중요한 점은 일반적인 미분 공식과 달리 Ito 공식에는 다음 항이 추가된다는 것이다.

$$
\frac{1}{2}\sigma_t^2 f_{xx}(t,X_t)
$$

이 항은 브라운 운동의 2차 변동성 때문에 생기는 항이며, 옵션 가격 함수의 볼록성, 즉 Gamma 효과를 반영한다.

> **Delta와 Gamma**
>
> **Delta**는 $\Delta=f_S(t,S)$로 나타내며, 주가가 1단위 변할 때 옵션 가격이 얼마나 변하는지를 나타내는 **주가 민감도**이다. 복제 포트폴리오에서는 옵션을 복제하기 위해 보유해야 할 주식 수량으로 쓰인다.
>
> **Gamma**는 $\Gamma=f_{SS}(t,S)$로 나타내며, 주가가 변할 때 Delta가 얼마나 변하는지를 나타내는 **Delta의 민감도**이다. 옵션 가격의 곡률과 변동성 효과를 설명할 때 쓰인다.

---

### 9. 두 가격 변화식 비교

이제 같은 대상인 $df(t,S_t)$를 두 가지 방식으로 구했다.

첫 번째는 **복제 포트폴리오의 자기금융 조건**으로부터 얻은 식이다.

$$
\color{#2563eb}{df(t,S_t) = \left[rf(t,S_t) + (\mu-r)S_t f_S(t,S_t)\right]dt + \sigma S_t f_S(t,S_t)dB_t}
$$

두 번째는 **Ito 공식**으로부터 얻은 식이다.

$$
\color{#2563eb}{df(t,S_t) = \left[f_t(t,S_t) + \mu S_t f_S(t,S_t) + \frac{1}{2}\sigma^2S_t^2 f_{SS}(t,S_t)\right]dt + \sigma S_t f_S(t,S_t)dB_t}
$$

두 식은 모두 같은 $df(t,S_t)$를 나타내므로 서로 같아야 한다.

먼저 확률항인 $dB_t$ 부분을 비교하면 다음과 같다.

$$
\sigma S_t f_S(t,S_t)dB_t
$$

확률항은 두 식에서 동일하다.

따라서 이제 $dt$가 붙은 drift 부분만 비교하면 된다.

복제 포트폴리오에서 얻은 drift는 다음과 같다.

$$
rf(t,S_t)+(\mu-r)S_t f_S(t,S_t)
$$

Ito 공식에서 얻은 drift는 다음과 같다.

$$
f_t(t,S_t)+\mu S_t f_S(t,S_t)+\frac{1}{2}\sigma^2S_t^2 f_{SS}(t,S_t)
$$

따라서 두 drift를 같게 두면 다음과 같다.

$$
\color{#2563eb}{rf(t,S_t)+(\mu-r)S_t f_S(t,S_t) = f_t(t,S_t)+\mu S_t f_S(t,S_t)+\frac{1}{2}\sigma^2S_t^2 f_{SS}(t,S_t)}
$$

이제 좌변을 전개하면 다음과 같다.

$$
rf(t,S_t)+\mu S_t f_S(t,S_t)-rS_t f_S(t,S_t) = f_t(t,S_t)+\mu S_t f_S(t,S_t)+\frac{1}{2}\sigma^2S_t^2 f_{SS}(t,S_t)
$$

양변에 있는 다음 항이 서로 소거된다.

$$
\mu S_t f_S(t,S_t)
$$

따라서 다음 식이 남는다.

$$
rf(t,S_t)-rS_t f_S(t,S_t) = f_t(t,S_t)+\frac{1}{2}\sigma^2S_t^2 f_{SS}(t,S_t)
$$

모든 항을 한쪽으로 모으면 다음과 같다.

$$
\color{#2563eb}{f_t(t,S_t)+\frac{1}{2}\sigma^2S_t^2 f_{SS}(t,S_t)+rS_t f_S(t,S_t)-rf(t,S_t)=0}
$$

일반적인 변수 $S$에 대해 쓰면 다음과 같다.

$$
\color{#2563eb}{f_t(t,S)+\frac{1}{2}\sigma^2S^2 f_{SS}(t,S)+rS f_S(t,S)-rf(t,S)=0}
$$

이것이 Black-Scholes-Merton PDE이다.

원문처럼 묶어서 쓰면 다음과 같다.

$$
\color{#2563eb}{r\left(Sf_S(t,S)-f(t,S)\right)+f_t(t,S)+\frac{1}{2}\sigma^2S^2f_{SS}(t,S)=0}
$$

여기서 중요한 점은 주식의 기대수익률을 나타내는 $\mu$가 최종 PDE에서 사라지며, 옵션 가격을 결정할 때 주식의 실제 기대수익률 $\mu$는 직접적으로 필요하지 않다는 점을 알 수 있다.

위험을 더 싫어하는 투자자라면 위험자산을 보유하기 위해 더 높은 기대수익률을 요구할 것이고, 그 결과 $\mu$는 달라질 수 있다.

하지만 Black-Scholes-Merton 모형에서는 이런 위험선호가 옵션 가격 공식에 직접적으로 $\mu$의 형태로 들어가지 않는다.

 옵션 가격을 구할 때 필요한 것은 시장에서 직접 관찰 가능한 현재 주가 $S_0$, 변동성 $\sigma$, 무위험 이자율 $r$, 만기 $T$, 행사가격 $K$이다.

---

### 10. 전체 유도 과정 정리

이번 절의 흐름은 다음과 같다.

먼저 주식과 무위험 자산의 가격 과정을 설정한다.

$$
dM_t = rM_tdt
$$

$$
dS_t = S_t(\mu dt + \sigma dB_t)
$$

다음으로 옵션 payoff를 복제하는 포트폴리오를 가정한다.

$$
V_t = \theta_t^{(S)}S_t + \theta_t^{(M)}M_t
$$

무차익 원리에 의해 옵션 가격과 복제 포트폴리오의 가치는 같아야 한다.

$$
f(t, S_t) = V_t
$$

자기금융 조건을 이용해 $df(t, S_t)$를 구하면 다음과 같다.

$$
df(t, S_t)=\left[rf(t, S_t) + (\mu - r)S_t f_S(t, S_t)\right]dt+S_t\sigma f_S(t, S_t)dB_t
$$

한편 Ito 공식으로도 $df(t, S_t)$를 구할 수 있다.

$$
df(t, S_t)=\left(f_t+\mu S_t f_S+\frac{1}{2}\sigma^2 S_t^2 f_{SS}\right)dt+\sigma S_t f_S dB_t
$$

두 식의 drift 부분을 같게 놓고 정리하면 Black-Scholes-Merton PDE를 얻는다.

$$
f_t+\frac{1}{2}\sigma^2 S^2 f_{SS}+rS f_S-rf=0
$$

최종 PDE에는 주식의 기대수익률 $\mu$가 등장하지 않는다.


옵션은 주식에서 파생된 상품이므로 주식의 기대수익률 $\mu$가 옵션 가격에 중요할 것처럼 보인다.

하지만 Black-Scholes-Merton 모형에서는 옵션을 주식과 무위험 자산으로 복제할 수 있다고 본다.

복제가 가능하다면 옵션 가격은 투자자의 주관적 기대수익률이 아니라, 복제에 필요한 비용으로 결정된다.

따라서 옵션 가격은 

> 같은 payoff를 만드는 두 포트폴리오는 현재 가격도 같아야 한다.

라는 무차익 원리에 의해 결정되게 한다.

 Black-Scholes-Merton PDE는 옵션 가격이 만족해야 하는 무차익 조건을 수식으로 표현한 것이다.

---

## 참고자료

<div style="border-left: 4px solid #2563eb; padding: 12px 16px; background-color: #f8fafc; border-radius: 6px;">

**Yoo, Shiyong.**  
*Black-Scholes-Merton Option Pricing Formula.*  
Lecture Note.

본 글의 Black-Scholes-Merton PDE 유도는 위 강의자료의  
**0.6 The Black-Scholes-Merton PDE 2** 절을 참고하여 정리하였다.

</div>
