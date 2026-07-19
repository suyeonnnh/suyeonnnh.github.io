---
title: "Finite Difference Methods for Option Pricing"
date: 2026-07-20 00:00:00 +0900
categories: [Numerical Analysis]
tags: [Finite Difference Method, Explicit FDM, Implicit FDM, Black-Scholes, Stability, Consistency, Convergence] 
math: true
mermaid: true
---

## 1. 유한차분법이 필요한 이유

옵션가격 \(V(t,S)\)는 Black–Scholes–Merton 편미분방정식과 계약별 만기조건·경계조건에 의해 결정된다. 유럽형 콜옵션과 풋옵션에는 해석해가 알려져 있다. 계약조건이 복잡해지면 편미분방정식의 해를 닫힌 형태로 구하기 어려워지며, 이때 수치해법을 사용한다.[1][2]

유한차분법은 시간과 주가를 격자로 나누고 편미분을 격자점 사이의 차분으로 근사한다. 편미분방정식은 격자점의 옵션가격들이 만족하는 대수방정식으로 바뀐다. 계산은 만기 payoff에서 출발하여 현재시점까지 진행된다.[1][3]

```mermaid
flowchart LR
    A["Black–Scholes–Merton PDE"] --> B["시간·주가 격자"]
    B --> C["편미분의 차분 근사"]
    C --> D["Explicit FDM"]
    C --> E["Implicit FDM"]
    D --> F["일치성·안정성·수렴성"]
    E --> F
```

---

## 2. Black–Scholes–Merton 편미분방정식

배당수익률이 $q$이고 변동성이 $\sigma$인 기초자산을 생각하자. 옵션가격 $V(t,S)$는 다음 Black–Scholes–Merton 편미분방정식을 만족한다.[1][2]

$$
\displaystyle \frac{\partial V}{\partial t}+\frac{1}{2}\sigma^2S^2\frac{\partial^2V}{\partial S^2}+(r-q)S\frac{\partial V}{\partial S}-rV=0
$$

여기서 $r$은 무위험이자율이다. 옵션의 payoff는 만기시점에서 주어지므로 잔존만기 변수 $\tau=T-t$를 정의한다. 이때 다음 관계가 성립한다.

$$
\displaystyle \frac{\partial V}{\partial t}=-\frac{\partial V}{\partial\tau}
$$

이를 Black–Scholes–Merton 편미분방정식에 대입하면 다음과 같다.

$$
\displaystyle \frac{\partial V}{\partial\tau}=\frac{1}{2}\sigma^2S^2\frac{\partial^2V}{\partial S^2}+(r-q)S\frac{\partial V}{\partial S}-rV
$$

$\tau=0$은 만기시점이고 $\tau=T$는 현재시점이다.

유럽형 콜옵션과 풋옵션의 만기조건은 각각 다음과 같다.

$$V(0,S)=\max(S-K,0)$$

$$V(0,S)=\max(K-S,0)$$

유럽형 콜옵션의 대표적인 경계조건은 다음과 같다.

$$V(\tau,0)=0,\qquad V(\tau,S_{\max})\approx S_{\max}e^{-q\tau}-Ke^{-r\tau}$$

유럽형 풋옵션의 대표적인 경계조건은 다음과 같다.

$$V(\tau,0)=Ke^{-r\tau},\qquad V(\tau,S_{\max})\approx0$$

---

## 3. 시간과 주가의 격자

주가구간 $0\le S\le S_{\max}$를 $M$등분하고, 잔존만기 구간 $0\le\tau\le T$를 $N$등분한다.

$$\Delta S=\frac{S_{\max}}{M},\qquad \Delta\tau=\frac{T}{N}$$

각 격자점은 다음과 같다.

$$S_j=j\Delta S,\qquad j=0,1,\ldots,M$$

$$\tau_n=n\Delta\tau,\qquad n=0,1,\ldots,N$$

격자점 $(\tau_n,S_j)$에서의 옵션가격을 $V_j^n=V(\tau_n,S_j)$로 나타낸다.

<div class="text-center">
<svg viewBox="0 0 760 390" width="100%" role="img" aria-label="시간과 주가의 유한차분 격자">
  <defs>
    <marker id="arrow-grid" markerWidth="8" markerHeight="8" refX="6" refY="3" orient="auto">
      <path d="M0,0 L0,6 L7,3 z" fill="currentColor"></path>
    </marker>
  </defs>
  <line x1="90" y1="325" x2="705" y2="325" stroke="currentColor" stroke-width="2" marker-end="url(#arrow-grid)"></line>
  <line x1="90" y1="325" x2="90" y2="45" stroke="currentColor" stroke-width="2" marker-end="url(#arrow-grid)"></line>
  <g stroke="currentColor" stroke-width="1" opacity="0.55">
    <line x1="90" y1="325" x2="650" y2="325"></line>
    <line x1="90" y1="265" x2="650" y2="265"></line>
    <line x1="90" y1="205" x2="650" y2="205"></line>
    <line x1="90" y1="145" x2="650" y2="145"></line>
    <line x1="90" y1="85" x2="650" y2="85"></line>
    <line x1="90" y1="85" x2="90" y2="325"></line>
    <line x1="230" y1="85" x2="230" y2="325"></line>
    <line x1="370" y1="85" x2="370" y2="325"></line>
    <line x1="510" y1="85" x2="510" y2="325"></line>
    <line x1="650" y1="85" x2="650" y2="325"></line>
  </g>
  <g fill="currentColor">
    <circle cx="90" cy="325" r="4"></circle><circle cx="230" cy="325" r="4"></circle><circle cx="370" cy="325" r="4"></circle><circle cx="510" cy="325" r="4"></circle><circle cx="650" cy="325" r="4"></circle>
    <circle cx="90" cy="265" r="4"></circle><circle cx="230" cy="265" r="4"></circle><circle cx="370" cy="265" r="4"></circle><circle cx="510" cy="265" r="4"></circle><circle cx="650" cy="265" r="4"></circle>
    <circle cx="90" cy="205" r="4"></circle><circle cx="230" cy="205" r="4"></circle><circle cx="370" cy="205" r="4"></circle><circle cx="510" cy="205" r="4"></circle><circle cx="650" cy="205" r="4"></circle>
    <circle cx="90" cy="145" r="4"></circle><circle cx="230" cy="145" r="4"></circle><circle cx="370" cy="145" r="4"></circle><circle cx="510" cy="145" r="4"></circle><circle cx="650" cy="145" r="4"></circle>
    <circle cx="90" cy="85" r="4"></circle><circle cx="230" cy="85" r="4"></circle><circle cx="370" cy="85" r="4"></circle><circle cx="510" cy="85" r="4"></circle><circle cx="650" cy="85" r="4"></circle>
  </g>
  <g fill="currentColor" font-size="16">
    <text x="695" y="355">잔존만기 τ</text>
    <text x="25" y="55">주가 S</text>
    <text x="72" y="350">τ₀</text>
    <text x="208" y="350">τ₁</text>
    <text x="348" y="350">τ₂</text>
    <text x="488" y="350">τ₃</text>
    <text x="625" y="350">τₙ</text>
    <text x="35" y="330">S₀</text>
    <text x="28" y="270">Sⱼ₋₁</text>
    <text x="45" y="210">Sⱼ</text>
    <text x="23" y="150">Sⱼ₊₁</text>
    <text x="23" y="90">Sₘ</text>
    <text x="73" y="378">만기</text>
    <text x="630" y="378">현재</text>
  </g>
</svg>
</div>

주가축의 내부 격자점은 $j=1,\ldots,M-1$이다. 양쪽 끝점인 $j=0$과 $j=M$에서는 각각 주가의 하한과 상한에 해당하는 경계조건을 사용한다.

---

## 4. 편미분의 차분 근사

차분식은 Taylor 전개에서 얻어진다.[3][4]

### 4.1 시간에 대한 전진차분

$$u(x,t+\Delta t)=u(x,t)+\Delta t\,u_t(x,t)+\frac{(\Delta t)^2}{2}u_{tt}(x,t)+O((\Delta t)^3)$$

양변에서 \(u(x,t)\)를 빼고 \(\Delta t\)로 나누면 다음 식을 얻는다.

$$\frac{u(x,t+\Delta t)-u(x,t)}{\Delta t}=u_t(x,t)+\frac{\Delta t}{2}u_{tt}(x,t)+O((\Delta t)^2)$$

따라서 전진차분은 시간에 대해 1차 정확도를 갖는다.

$$u_t(x,t)=\frac{u(x,t+\Delta t)-u(x,t)}{\Delta t}+O(\Delta t)$$

### 4.2 시간에 대한 후진차분

$$u(x,t-\Delta t)=u(x,t)-\Delta t\,u_t(x,t)+\frac{(\Delta t)^2}{2}u_{tt}(x,t)+O((\Delta t)^3)$$

정리하면 다음 식을 얻는다.

$$\frac{u(x,t)-u(x,t-\Delta t)}{\Delta t}=u_t(x,t)-\frac{\Delta t}{2}u_{tt}(x,t)+O((\Delta t)^2)$$

후진차분도 시간에 대해 1차 정확도를 갖는다.

### 4.3 공간에 대한 중심차분

$$u(x+\Delta x,t)=u(x,t)+\Delta x\,u_x(x,t)+\frac{(\Delta x)^2}{2}u_{xx}(x,t)+\frac{(\Delta x)^3}{6}u_{xxx}(x,t)+O((\Delta x)^4)$$

$$u(x-\Delta x,t)=u(x,t)-\Delta x\,u_x(x,t)+\frac{(\Delta x)^2}{2}u_{xx}(x,t)-\frac{(\Delta x)^3}{6}u_{xxx}(x,t)+O((\Delta x)^4)$$

두 식을 빼면 다음 1차 중심차분을 얻는다.

$$u_x(x,t)=\frac{u(x+\Delta x,t)-u(x-\Delta x,t)}{2\Delta x}+O((\Delta x)^2)$$

두 식을 더하면 다음 2차 중심차분을 얻는다.

$$u_{xx}(x,t)=\frac{u(x+\Delta x,t)-2u(x,t)+u(x-\Delta x,t)}{(\Delta x)^2}+O((\Delta x)^2)$$

---

## 5. Explicit Finite Difference Method

### 5.1 Black–Scholes–Merton PDE에서의 도출

Explicit FDM은 공간미분을 이미 계산된 시간단계 \(n\)에서 평가한다.

$$\frac{V_j^{n+1}-V_j^n}{\Delta\tau}=\frac{1}{2}\sigma^2S_j^2\frac{V_{j+1}^n-2V_j^n+V_{j-1}^n}{(\Delta S)^2}+(r-q)S_j\frac{V_{j+1}^n-V_{j-1}^n}{2\Delta S}-rV_j^n$$

$S_j=j\Delta S$를 대입하면 $S_j/\Delta S=j$이고 $S_j^2/(\Delta S)^2=j^2$이다.

이를 $V_j^{n+1}$에 대해 정리하면 다음 식을 얻는다.


$$V_j^{n+1}=a_jV_{j-1}^n+b_jV_j^n+c_jV_{j+1}^n$$

$$a_j=\frac{\Delta\tau}{2}\left(\sigma^2j^2-(r-q)j\right)$$

$$b_j=1-\Delta\tau\left(\sigma^2j^2+r\right)$$

$$c_j=\frac{\Delta\tau}{2}\left(\sigma^2j^2+(r-q)j\right)$$

<div class="text-center">
<svg viewBox="0 0 720 250" width="100%" role="img" aria-label="Explicit FDM stencil">
  <defs>
    <marker id="arrow-explicit" markerWidth="8" markerHeight="8" refX="6" refY="3" orient="auto">
      <path d="M0,0 L0,6 L7,3 z" fill="currentColor"></path>
    </marker>
  </defs>
  <g stroke="currentColor" stroke-width="1" opacity="0.45">
    <line x1="80" y1="60" x2="640" y2="60"></line><line x1="80" y1="190" x2="640" y2="190"></line>
    <line x1="170" y1="35" x2="170" y2="215"></line><line x1="360" y1="35" x2="360" y2="215"></line><line x1="550" y1="35" x2="550" y2="215"></line>
  </g>
  <g stroke="currentColor" stroke-width="2" marker-end="url(#arrow-explicit)">
    <line x1="170" y1="60" x2="355" y2="185"></line>
    <line x1="360" y1="60" x2="360" y2="183"></line>
    <line x1="550" y1="60" x2="365" y2="185"></line>
  </g>
  <g fill="currentColor">
    <circle cx="170" cy="60" r="6"></circle><circle cx="360" cy="60" r="6"></circle><circle cx="550" cy="60" r="6"></circle>
    <circle cx="360" cy="190" r="9" fill="none" stroke="currentColor" stroke-width="3"></circle>
  </g>
  <g fill="currentColor" font-size="17">
    <text x="130" y="35">Vⱼ₋₁ⁿ</text><text x="330" y="35">Vⱼⁿ</text><text x="515" y="35">Vⱼ₊₁ⁿ</text>
    <text x="325" y="230">Vⱼⁿ⁺¹</text><text x="65" y="65">n</text><text x="35" y="195">n+1</text>
  </g>
</svg>
</div>

> **Explicit FDM 계산 구조**
>
> 1. $\tau=0$에서 만기 payoff를 입력한다.
> 2. 시간단계 $n$에서 아는 세 격자값 $V_{j-1}^n$, $V_j^n$, $V_{j+1}^n$을 사용하여 $V_j^{n+1}$을 계산한다.
> 3. 새로운 시간단계의 양쪽 경계조건을 입력한다.
> 4. $\tau=T$에 도달할 때까지 같은 계산을 반복한다.
{: .prompt-info }
---

## 6. Implicit Finite Difference Method

### 6.1 Black–Scholes–Merton PDE에서의 도출

Implicit FDM은 공간미분을 새로 구할 시간단계 \(n+1\)에서 평가한다.

$$\frac{V_j^{n+1}-V_j^n}{\Delta\tau}=\frac{1}{2}\sigma^2S_j^2\frac{V_{j+1}^{n+1}-2V_j^{n+1}+V_{j-1}^{n+1}}{(\Delta S)^2}+(r-q)S_j\frac{V_{j+1}^{n+1}-V_{j-1}^{n+1}}{2\Delta S}-rV_j^{n+1}$$

\(S_j=j\Delta S\)를 대입하고 시간단계 \(n+1\)의 항들을 왼쪽으로 모으면 다음 식을 얻는다.

$$\ell_jV_{j-1}^{n+1}+d_jV_j^{n+1}+u_jV_{j+1}^{n+1}=V_j^n$$

$$\ell_j=-\frac{\Delta\tau}{2}\left(\sigma^2j^2-(r-q)j\right)$$

$$d_j=1+\Delta\tau\left(\sigma^2j^2+r\right)$$

$$u_j=-\frac{\Delta\tau}{2}\left(\sigma^2j^2+(r-q)j\right)$$

<div class="text-center">
<svg viewBox="0 0 720 250" width="100%" role="img" aria-label="Implicit FDM stencil">
  <defs>
    <marker id="arrow-implicit" markerWidth="8" markerHeight="8" refX="6" refY="3" orient="auto">
      <path d="M0,0 L0,6 L7,3 z" fill="currentColor"></path>
    </marker>
  </defs>
  <g stroke="currentColor" stroke-width="1" opacity="0.45">
    <line x1="80" y1="60" x2="640" y2="60"></line><line x1="80" y1="190" x2="640" y2="190"></line>
    <line x1="170" y1="35" x2="170" y2="215"></line><line x1="360" y1="35" x2="360" y2="215"></line><line x1="550" y1="35" x2="550" y2="215"></line>
  </g>
  <g stroke="currentColor" stroke-width="2" marker-end="url(#arrow-implicit)">
    <line x1="360" y1="60" x2="175" y2="185"></line>
    <line x1="360" y1="60" x2="360" y2="183"></line>
    <line x1="360" y1="60" x2="545" y2="185"></line>
  </g>
  <g fill="currentColor">
    <circle cx="360" cy="60" r="6"></circle>
    <circle cx="170" cy="190" r="9" fill="none" stroke="currentColor" stroke-width="3"></circle>
    <circle cx="360" cy="190" r="9" fill="none" stroke="currentColor" stroke-width="3"></circle>
    <circle cx="550" cy="190" r="9" fill="none" stroke="currentColor" stroke-width="3"></circle>
  </g>
  <g fill="currentColor" font-size="17">
    <text x="330" y="35">Vⱼⁿ</text>
    <text x="125" y="230">Vⱼ₋₁ⁿ⁺¹</text><text x="325" y="230">Vⱼⁿ⁺¹</text><text x="505" y="230">Vⱼ₊₁ⁿ⁺¹</text>
    <text x="65" y="65">n</text><text x="35" y="195">n+1</text>
  </g>
</svg>
</div>

내부 격자값을 벡터로 모으면 다음과 같다.

$$
A\mathbf{V}^{n+1}=\mathbf{V}^n+\mathbf{g}^{n+1}
$$

$$
A=\begin{bmatrix}d_1&u_1&0&\cdots&0\\\ell_2&d_2&u_2&\ddots&\vdots\\0&\ell_3&d_3&\ddots&0\\\vdots&\ddots&\ddots&\ddots&u_{M-2}\\0&\cdots&0&\ell_{M-1}&d_{M-1}\end{bmatrix}
$$

여기서 $\mathbf{g}^{n+1}$은 주가축 양쪽 끝에서 주어지는 경계조건을 반영한 벡터이다. 행렬 $A$는 주대각선과 그 위아래의 대각선에만 값이 존재하는 삼중대각행렬이다. 따라서 일반적인 역행렬을 직접 계산하지 않고 Thomas algorithm을 사용하여 효율적으로 해를 구할 수 있다.

> **Implicit FDM 계산 구조**
>
> 1. $\tau=0$에서 만기 payoff를 입력한다.
> 2. 계수 $\ell_j$, $d_j$, $u_j$를 사용하여 삼중대각행렬 $A$를 구성한다.
> 3. 해당 시간단계의 경계조건을 벡터 $\mathbf{g}^{n+1}$에 반영한다.
> 4. 선형방정식 $A\mathbf{V}^{n+1}=\mathbf{V}^n+\mathbf{g}^{n+1}$을 풀어 새로운 시간단계의 옵션가격을 구한다.
> 5. $\tau=T$에 도달할 때까지 같은 계산을 반복한다.
{: .prompt-info }

---

## 7. 안정성 분석을 위한 열전도방정식 변환

Black–Scholes–Merton 편미분방정식은 로그주가와 시간변환을 거쳐 열전도방정식으로 바뀐다.[3][4]

다음 변수를 정의한다.

$$x=\ln\left(\frac{S}{K}\right),\qquad \theta=\frac{\sigma^2}{2}(T-t),\qquad v(x,\theta)=\frac{V(t,S)}{K}$$

연쇄법칙을 적용하면 다음 식을 얻는다.

$$V_t=-\frac{\sigma^2K}{2}v_\theta,\qquad V_S=\frac{K}{S}v_x,\qquad V_{SS}=\frac{K}{S^2}(v_{xx}-v_x)$$

이를 Black–Scholes–Merton 편미분방정식에 대입하면 다음 식을 얻는다.

$$v_\theta=v_{xx}+av_x-bv$$

$$a=\frac{2(r-q)}{\sigma^2}-1,\qquad b=\frac{2r}{\sigma^2}$$

이제 $v=e^{\alpha x+\beta\theta}u$로 두고 다음 값을 선택한다.

$$\alpha=-\frac{a}{2},\qquad \beta=\alpha^2+a\alpha-b=-\frac{a^2}{4}-b$$

그러면 1차 미분항과 0차항이 제거되어 열전도방정식을 얻는다.

$$u_\theta=u_{xx}$$

이후의 일치성·안정성·수렴성 증명은 이 열전도방정식을 기준으로 진행한다. 공간격자를 $\Delta x$, 시간격자를 $\Delta\theta$로 두고 다음 값을 정의한다.

$$\gamma=\frac{\Delta\theta}{(\Delta x)^2}$$

---

## 8. 일치성·안정성·수렴성의 정의

정확해의 격자값을 $U_i^n=u(x_i,\theta_n)$, 차분법으로 계산한 수치해를 $u_i^n$으로 나타낸다.

### 8.1 일치성

정확해 $U_i^n$을 차분방정식에 대입했을 때 남는 오차를 국소절단오차 $R_i^n$이라고 한다. 공간격자와 시간격자를 세분화할 때 국소절단오차가 0으로 수렴하면 해당 차분법은 일치성을 갖는다.

$$
\lim_{\Delta\theta\to0,\ \Delta x\to0}R_i^n=0
$$

### 8.2 안정성

초기조건이나 계산 과정에서 발생한 작은 오차가 시간단계를 반복하면서 계속 증폭되지 않으면 차분법은 안정적이다. von Neumann 안정성 분석에서는 오차를 다음과 같은 Fourier 모드로 나타낸다.

$$
\varepsilon_i^n=G^ne^{\mathrm{i}\xi i\Delta x}
$$

여기서 $G$는 한 시간단계가 진행될 때 오차가 얼마나 증폭되는지를 나타내는 증폭인자이고, $\xi$는 파수이다. 모든 파수 $\xi$에 대해 $|G|\le1$이면 해당 차분법은 안정적이다.
### 8.3 수렴성

정확한 해와 수치해의 차이를 전역오차라고 한다.

$$e_i^n=U_i^n-u_i^n$$

격자를 세분화할 때 전역오차가 0으로 수렴하면 수치해는 수렴성을 갖는다.

$$\lim_{\Delta\theta\to0,\ \Delta x\to0}e_i^n=0$$

---

## 9. Explicit FDM의 일치성·안정성·수렴성

### 9.1 Explicit FDM의 일치성 도출

열전도방정식의 Explicit 차분식은 다음과 같다.

$$\frac{u_i^{n+1}-u_i^n}{\Delta\theta}=\frac{u_{i+1}^n-2u_i^n+u_{i-1}^n}{(\Delta x)^2}$$

정확한 해 \(U_i^n\)을 대입하여 국소절단오차를 정의한다.

$$R_{i,\mathrm{E}}^n=\frac{U_i^{n+1}-U_i^n}{\Delta\theta}-\frac{U_{i+1}^n-2U_i^n+U_{i-1}^n}{(\Delta x)^2}$$

시간항을 \((x_i,\theta_n)\)에서 Taylor 전개하면 다음 식을 얻는다.

$$\frac{U_i^{n+1}-U_i^n}{\Delta\theta}=u_\theta+\frac{\Delta\theta}{2}u_{\theta\theta}+O((\Delta\theta)^2)$$

공간항을 같은 점에서 Taylor 전개하면 다음 식을 얻는다.

$$\frac{U_{i+1}^n-2U_i^n+U_{i-1}^n}{(\Delta x)^2}=u_{xx}+\frac{(\Delta x)^2}{12}u_{xxxx}+O((\Delta x)^4)$$

두 전개식을 국소절단오차에 대입한다.

$$R_{i,\mathrm{E}}^n=u_\theta-u_{xx}+\frac{\Delta\theta}{2}u_{\theta\theta}-\frac{(\Delta x)^2}{12}u_{xxxx}+O((\Delta\theta)^2)+O((\Delta x)^4)$$

정확한 해는 \(u_\theta=u_{xx}\)를 만족하므로 첫 두 항이 소거된다.

$$R_{i,\mathrm{E}}^n=\frac{\Delta\theta}{2}u_{\theta\theta}-\frac{(\Delta x)^2}{12}u_{xxxx}+O((\Delta\theta)^2)+O((\Delta x)^4)$$

따라서 Explicit FDM의 국소절단오차는 다음 차수를 갖는다.[4]

$$R_{i,\mathrm{E}}^n=O(\Delta\theta)+O((\Delta x)^2)$$

Explicit FDM은 시간에 대해 1차, 공간에 대해 2차 일치성을 갖는다.

### 9.2 Explicit FDM의 안정성 도출

Explicit 차분식을 $\gamma=\Delta\theta/(\Delta x)^2$를 사용하여 정리하면 다음과 같다.

$$
u_i^{n+1}=\gamma u_{i-1}^n+(1-2\gamma)u_i^n+\gamma u_{i+1}^n
$$

Fourier 오차모드 $\varepsilon_i^n=G^ne^{\mathrm{i}\xi i\Delta x}$를 대입한다.

$$
G^{n+1}e^{\mathrm{i}\xi i\Delta x}=\gamma G^ne^{\mathrm{i}\xi(i-1)\Delta x}+(1-2\gamma)G^ne^{\mathrm{i}\xi i\Delta x}+\gamma G^ne^{\mathrm{i}\xi(i+1)\Delta x}
$$

양변을 $G^ne^{\mathrm{i}\xi i\Delta x}$로 나누면 다음 식을 얻는다.

$$
G=\gamma e^{-\mathrm{i}\xi\Delta x}+(1-2\gamma)+\gamma e^{\mathrm{i}\xi\Delta x}
$$

$e^{\mathrm{i}z}+e^{-\mathrm{i}z}=2\cos z$와 $1-\cos z=2\sin^2(z/2)$를 사용하면 증폭인자는 다음과 같이 정리된다.

$$
G=1-2\gamma+2\gamma\cos(\xi\Delta x)=1-4\gamma\sin^2\left(\frac{\xi\Delta x}{2}\right)
$$

안정성을 위해서는 모든 파수 $\xi$에 대해 $|G|\le1$이 성립해야 한다. $\sin^2(\xi\Delta x/2)$의 범위가 $[0,1]$이므로 다음 조건이 필요하다.

$$
-1\le1-4\gamma\sin^2\left(\frac{\xi\Delta x}{2}\right)\le1
$$

오른쪽 부등식은 $\gamma\ge0$일 때 항상 성립한다. 왼쪽 부등식은 $\sin^2(\xi\Delta x/2)=1$일 때 가장 강한 조건을 준다.

$$
-1\le1-4\gamma
$$

이를 정리하면 Explicit FDM의 안정성 조건을 얻는다.

$$
0\le\gamma\le\frac{1}{2}
$$

이를 격자간격으로 나타내면 다음과 같다.

$$
\Delta\theta\le\frac{(\Delta x)^2}{2}
$$

Black–Scholes 방정식의 잔존만기 변수에서는 $\Delta\theta=\sigma^2\Delta\tau/2$이므로 다음 조건을 얻는다.

$$
\frac{\sigma^2\Delta\tau}{(\Delta x)^2}\le1
$$

> **Explicit FDM의 안정성 조건**
>
> 공간격자 $\Delta x$를 절반으로 줄이면 시간격자 $\Delta\theta$는 기존 값의 4분의 1 이하로 줄여야 한다. 이 조건을 위반하면 고주파 오차모드가 시간단계를 거치면서 증폭될 수 있다.
{: .prompt-warning }

### 9.3 Explicit FDM의 수렴성 도출

내부 격자값을 하나의 벡터로 모으면 Explicit FDM은 다음과 같이 표현된다.

$$
\mathbf{u}^{n+1}=A_{\mathrm{E}}\mathbf{u}^n
$$

$$
A_{\mathrm{E}}=\begin{bmatrix}1-2\gamma&\gamma&0&\cdots&0\\\gamma&1-2\gamma&\gamma&\ddots&\vdots\\0&\gamma&1-2\gamma&\ddots&0\\\vdots&\ddots&\ddots&\ddots&\gamma\\0&\cdots&0&\gamma&1-2\gamma\end{bmatrix}
$$

정확해의 격자값은 국소절단오차를 포함하여 다음 식을 만족한다.

$$
\mathbf{U}^{n+1}=A_{\mathrm{E}}\mathbf{U}^n+\Delta\theta\,\mathbf{R}_{\mathrm{E}}^n
$$

전역오차를 $\mathbf{e}^n=\mathbf{U}^n-\mathbf{u}^n$으로 정의한다. 정확해가 만족하는 식에서 수치해가 만족하는 식을 빼면 다음 오차방정식을 얻는다.

$$
\mathbf{e}^{n+1}=A_{\mathrm{E}}\mathbf{e}^n+\Delta\theta\,\mathbf{R}_{\mathrm{E}}^n
$$

안정성 조건 $0\le\gamma\le1/2$에서는 $\|A_{\mathrm{E}}\|_2\le1$이므로 다음 부등식이 성립한다.

$$
\|\mathbf{e}^{n+1}\|_2\le\|\mathbf{e}^n\|_2+\Delta\theta\|\mathbf{R}_{\mathrm{E}}^n\|_2
$$

이를 반복하여 적용하면 다음과 같다.

$$
\|\mathbf{e}^n\|_2\le\|\mathbf{e}^0\|_2+\Delta\theta\sum_{m=0}^{n-1}\|\mathbf{R}_{\mathrm{E}}^m\|_2
$$

$n\Delta\theta\le T_\theta$를 사용하면 다음 상계를 얻는다.

$$
\|\mathbf{e}^n\|_2\le\|\mathbf{e}^0\|_2+T_\theta\max_{0\le m<n}\|\mathbf{R}_{\mathrm{E}}^m\|_2
$$

초기조건을 격자에 정확히 입력하면 $\mathbf{e}^0=\mathbf{0}$이다. 일치성 결과 $\mathbf{R}_{\mathrm{E}}^m=O(\Delta\theta)+O((\Delta x)^2)$를 대입하면 다음 결론을 얻는다.

$$
\|\mathbf{e}^n\|_2=O(\Delta\theta)+O((\Delta x)^2)
$$

따라서 Explicit FDM은 안정성 조건 $0\le\gamma\le1/2$ 아래에서 정확해로 수렴한다.

---

## 10. Implicit FDM의 일치성·안정성·수렴성

### 10.1 Implicit FDM의 일치성 도출

열전도방정식의 Implicit 차분식은 다음과 같다.

$$\frac{u_i^{n+1}-u_i^n}{\Delta\theta}=\frac{u_{i+1}^{n+1}-2u_i^{n+1}+u_{i-1}^{n+1}}{(\Delta x)^2}$$

정확한 해 \(U_i^n\)을 대입하여 국소절단오차를 정의한다.

$$R_{i,\mathrm{I}}^{n+1}=\frac{U_i^{n+1}-U_i^n}{\Delta\theta}-\frac{U_{i+1}^{n+1}-2U_i^{n+1}+U_{i-1}^{n+1}}{(\Delta x)^2}$$

시간항을 새 시간점 \((x_i,\theta_{n+1})\)에서 후진 Taylor 전개한다.

$$\frac{U_i^{n+1}-U_i^n}{\Delta\theta}=u_\theta-\frac{\Delta\theta}{2}u_{\theta\theta}+O((\Delta\theta)^2)$$

공간항도 \((x_i,\theta_{n+1})\)에서 전개한다.

$$\frac{U_{i+1}^{n+1}-2U_i^{n+1}+U_{i-1}^{n+1}}{(\Delta x)^2}=u_{xx}+\frac{(\Delta x)^2}{12}u_{xxxx}+O((\Delta x)^4)$$

두 식을 국소절단오차에 대입한다.

$$R_{i,\mathrm{I}}^{n+1}=u_\theta-u_{xx}-\frac{\Delta\theta}{2}u_{\theta\theta}-\frac{(\Delta x)^2}{12}u_{xxxx}+O((\Delta\theta)^2)+O((\Delta x)^4)$$

정확한 해는 \(u_\theta=u_{xx}\)를 만족하므로 다음 식을 얻는다.

$$R_{i,\mathrm{I}}^{n+1}=-\frac{\Delta\theta}{2}u_{\theta\theta}-\frac{(\Delta x)^2}{12}u_{xxxx}+O((\Delta\theta)^2)+O((\Delta x)^4)$$

따라서 Implicit FDM의 국소절단오차는 다음 차수를 갖는다.[4]

$$R_{i,\mathrm{I}}^{n+1}=O(\Delta\theta)+O((\Delta x)^2)$$

Implicit FDM도 시간에 대해 1차, 공간에 대해 2차 일치성을 갖는다.

### 10.2 Implicit FDM의 안정성 도출

Implicit 차분식을 $\gamma=\Delta\theta/(\Delta x)^2$를 사용하여 정리하면 다음과 같다.

$$
-\gamma u_{i-1}^{n+1}+(1+2\gamma)u_i^{n+1}-\gamma u_{i+1}^{n+1}=u_i^n
$$

Fourier 오차모드 $\varepsilon_i^n=G^ne^{\mathrm{i}\xi i\Delta x}$를 대입한다.

$$
-\gamma G^{n+1}e^{\mathrm{i}\xi(i-1)\Delta x}+(1+2\gamma)G^{n+1}e^{\mathrm{i}\xi i\Delta x}-\gamma G^{n+1}e^{\mathrm{i}\xi(i+1)\Delta x}=G^ne^{\mathrm{i}\xi i\Delta x}
$$

양변을 $G^ne^{\mathrm{i}\xi i\Delta x}$로 나누면 다음 식을 얻는다.

$$
G\left[-\gamma e^{-\mathrm{i}\xi\Delta x}+(1+2\gamma)-\gamma e^{\mathrm{i}\xi\Delta x}\right]=1
$$

$e^{\mathrm{i}z}+e^{-\mathrm{i}z}=2\cos z$를 사용하면 다음과 같이 정리된다.

$$
G\left[1+2\gamma-2\gamma\cos(\xi\Delta x)\right]=1
$$

$1-\cos z=2\sin^2(z/2)$를 사용하면 증폭인자는 다음과 같다.

$$
G=\frac{1}{1+4\gamma\sin^2\left(\frac{\xi\Delta x}{2}\right)}
$$

$\gamma>0$이고 $0\le\sin^2(\xi\Delta x/2)\le1$이므로 다음 부등식이 성립한다.

$$
0<\frac{1}{1+4\gamma}\le G\le1
$$

따라서 모든 $\gamma>0$에 대해 $|G|\le1$이다. 그러므로 Implicit FDM은 von Neumann 안정성 기준에서 무조건 안정적이다.

> **Implicit FDM의 안정성**
>
> 안정성을 위해 $\Delta\theta$와 $\Delta x$ 사이에 추가적인 비율조건을 설정할 필요가 없다. 그러나 격자가 거칠면 절단오차가 커질 수 있으므로 정확도를 높이기 위해서는 두 격자간격을 충분히 줄여야 한다.
{: .prompt-tip }

### 10.3 Implicit FDM의 수렴성 도출

내부 격자값을 하나의 벡터로 모으면 Implicit FDM은 다음과 같이 표현된다.

$$
B_{\mathrm{I}}\mathbf{u}^{n+1}=\mathbf{u}^n
$$

$$
B_{\mathrm{I}}=\begin{bmatrix}1+2\gamma&-\gamma&0&\cdots&0\\-\gamma&1+2\gamma&-\gamma&\ddots&\vdots\\0&-\gamma&1+2\gamma&\ddots&0\\\vdots&\ddots&\ddots&\ddots&-\gamma\\0&\cdots&0&-\gamma&1+2\gamma\end{bmatrix}
$$

정확해의 격자값은 국소절단오차를 포함하여 다음 식을 만족한다.

$$
B_{\mathrm{I}}\mathbf{U}^{n+1}=\mathbf{U}^n+\Delta\theta\,\mathbf{R}_{\mathrm{I}}^{n+1}
$$

전역오차를 $\mathbf{e}^n=\mathbf{U}^n-\mathbf{u}^n$으로 정의한다. 정확해가 만족하는 식에서 수치해가 만족하는 식을 빼면 다음 오차방정식을 얻는다.

$$
B_{\mathrm{I}}\mathbf{e}^{n+1}=\mathbf{e}^n+\Delta\theta\,\mathbf{R}_{\mathrm{I}}^{n+1}
$$

양변에 $B_{\mathrm{I}}^{-1}$을 곱하면 다음과 같다.

$$
\mathbf{e}^{n+1}=B_{\mathrm{I}}^{-1}\mathbf{e}^n+\Delta\theta B_{\mathrm{I}}^{-1}\mathbf{R}_{\mathrm{I}}^{n+1}
$$

Implicit FDM의 증폭인자는 모든 파수에서 $0<G\le1$을 만족하므로 $\|B_{\mathrm{I}}^{-1}\|_2\le1$이다. 따라서 다음 부등식이 성립한다.

$$
\|\mathbf{e}^{n+1}\|_2\le\|\mathbf{e}^n\|_2+\Delta\theta\|\mathbf{R}_{\mathrm{I}}^{n+1}\|_2
$$

이를 반복하여 적용하고 $n\Delta\theta\le T_\theta$를 사용하면 다음 상계를 얻는다.

$$
\|\mathbf{e}^n\|_2\le\|\mathbf{e}^0\|_2+T_\theta\max_{1\le m\le n}\|\mathbf{R}_{\mathrm{I}}^m\|_2
$$

초기조건을 격자에 정확히 입력하면 $\mathbf{e}^0=\mathbf{0}$이다. 여기에 일치성 결과 $\mathbf{R}_{\mathrm{I}}^m=O(\Delta\theta)+O((\Delta x)^2)$를 대입하면 다음 결론을 얻는다.

$$
\|\mathbf{e}^n\|_2=O(\Delta\theta)+O((\Delta x)^2)
$$

따라서 Implicit FDM은 모든 $\gamma>0$에서 안정적이며, 격자간격이 0으로 수렴할 때 정확해로 수렴한다.

---

## 11. Lax 동등정리와 두 방법의 결론

선형이며 적절성 조건을 만족하는 초기값문제에 대해, 일치성을 갖는 유한차분법은 안정성과 수렴성이 서로 동등하다. 이를 Lax 동등정리라고 한다.

$$
\text{일치성을 가정할 때}\qquad\text{안정성}\Longleftrightarrow\text{수렴성}
$$

앞 절에서는 Explicit FDM과 Implicit FDM의 국소절단오차를 계산하여 일치성을 확인했다. 이어서 Fourier 오차모드를 이용해 안정성을 분석하고, 전역오차의 점화식을 통해 수렴성을 확인했다.

| 방법 | 일치성 | 안정성 | 수렴성 |
|---|---|---|---|
| Explicit FDM | $O(\Delta\theta)+O((\Delta x)^2)$ | $0\le\gamma\le1/2$에서 안정 | 안정성 조건 아래에서 수렴 |
| Implicit FDM | $O(\Delta\theta)+O((\Delta x)^2)$ | 모든 $\gamma>0$에서 안정 | 격자간격이 0으로 수렴할 때 비율조건 없이 수렴 |

---

## 12. Explicit FDM과 Implicit FDM 비교

| 구분 | Explicit FDM | Implicit FDM |
|---|---|---|
| 공간미분 평가시점 | 시간단계 $n$ | 시간단계 $n+1$ |
| 계산구조 | 새 격자값을 직접 계산 | 같은 시간단계의 미지수들을 함께 계산 |
| 선형계 | 필요 없음 | 삼중대각 선형계 |
| 시간정확도 | 1차 | 1차 |
| 공간정확도 | 2차 | 2차 |
| 일치성 | 있음 | 있음 |
| 안정성 | 조건부 안정 | 무조건 안정 |
| 안정성 조건 | $0\le\gamma\le1/2$ | 모든 $\gamma>0$ |
| 한 시간단계의 계산량 | 비교적 작음 | 선형계 풀이가 필요 |
| 격자 세분화 | $\Delta x$를 줄이면 $\Delta\theta$도 안정성 조건에 맞게 줄여야 함 | 안정성에 관한 비율조건 없이 세분화 가능 |

강현주의 유럽형 풋옵션 수치실험에서는 $\Delta t$를 고정한 상태에서 $\Delta S$를 줄였을 때 Explicit FDM이 불안정해졌다. 반면 Implicit FDM은 같은 격자에서도 안정적인 값을 계산했다.

| 방법 | $\Delta S=10$ | $\Delta S=5$ | $\Delta S=1$ | $\Delta S=0.5$ |
|---|---:|---:|---:|---:|
| Explicit FDM | 18.8787 | 18.8531 | 불안정 | 불안정 |
| Implicit FDM | 18.8836 | 18.8592 | 18.8492 | 18.8489 |

해당 실험에서 Black–Scholes 해는 18.8455이다.

---

## 13. 정리

유한차분법은 Black–Scholes–Merton 편미분방정식을 격자 위의 차분방정식으로 변환하여 옵션가격을 계산한다.

Explicit FDM은 시간단계 $n$의 세 격자값을 이용하여 시간단계 $n+1$의 값을 직접 계산한다. 시간에 대해 1차, 공간에 대해 2차의 정확도를 가지며, $0\le\gamma\le1/2$에서 안정적이고 수렴한다.

Implicit FDM에서는 시간단계 $n+1$의 격자값들이 하나의 삼중대각 선형계로 연결된다. 시간에 대해 1차, 공간에 대해 2차의 정확도를 가지며, 모든 $\gamma>0$에서 안정적이다. 또한 $\Delta\theta$와 $\Delta x$가 0으로 수렴하면 두 격자간격 사이의 추가적인 비율조건 없이 정확해로 수렴한다.

---

## References

[1] 이재석. (2006). *옵션 가치평가에 관한 연구: FDM과 Monte Carlo Simulation 기법의 비교*. 중앙대학교 대학원 석사학위논문.

[2] Yoo, Shiyong. (n.d.). *Black-Scholes-Merton Option Pricing Formula*. Course Notes.

[3] 강현주. (2010). *Numerical Methods for Pricing of European Option*. 이화여자대학교 대학원 석사학위논문.

[4] 최병선. (2007). *계산재무론: Computational Finance*. 세경사.
