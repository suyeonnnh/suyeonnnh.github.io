---
title: "Optimal Nonparametric Estimation of the Expected Shortfall Risk 5"
date: 2026-10-24 19:00:00 +0100
categories:
  - Computational Finance
tags:
  - Expected Shortfall
  - Numerical Experiments
  - Robust Statistics
  - Monte Carlo
  - Heavy Tails
math: true
---
이번 5편에서는 **Numerical Experiments**를 살펴본다. 저자들이 제시한 추정량이 실제 표본에서는 어떤 결과를 보이는지, 여러 손실분포와 오차 수준을 바꿔가며 비교한다. 특히 **작은 추정 오차가 중요한 상황**, **큰 오차를 최대한 피해야 하는 상황**, **데이터 일부가 변조된 상황**을 구분해서 읽어본다.

논문 Section 5에는 **Figure 2부터 Figure 7까지** 총 여섯 개의 그림이 있다. 본문에서는 각 그림의 실험 조건과 관찰 결과를 먼저 설명하고, 마지막에 핵심 알고리즘을 Python으로 재현할 수 있는 예제를 추가한다.

---

# 5. Numerical Experiments (p. 26–32)

## 5.1 Experimental Setup — 실험에서 무엇을 비교하는가?

저자들은 실험 전반에서 **α = 0.1**을 사용한다. 가장 나쁜 상위 10% 손실에 대한 ES를 추정하는 설정이다. α가 작을수록 꼬리에 포함되는 데이터가 줄어들기 때문에, 충분한 반복실험을 수행하는 계산 비용이 커질 수 있다.

세 가지 기본 추정량을 비교한다.

**① Plug-in estimator (T̂_N)**

전체 데이터를 이용하여 경험적 ES를 계산한다. 1편에서 구현한 것처럼, 표본의 상위 α 비율을 가중 평균하는 방식이다. 표본이 충분하고 손실분포가 안정적인 상황에서는 효율적으로 작동할 수 있다.

**② Median of Blocks estimator**

데이터를 여러 블록으로 나누고 블록마다 ES를 계산한 뒤 그 추정값들의 중앙값을 사용한다.

![Median of blocks](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%20%5Cwidehat%7BM%7D_N%3D%5Coperatorname%7Bmedian%7D%5C%21%5Cleft%5C%7B%5Cwidehat%7BT%7D_%7BI_1%7D%2C%5Cldots%2C%5Cwidehat%7BT%7D_%7BI_n%7D%5Cright%5C%7D)

블록별 추정값 일부가 극단적으로 커져도 중앙값에 미치는 영향을 제한할 수 있다. 블록 크기가 작을 때는 추정량에 **편향(bias)**이 나타날 수 있다.

**③ Proposed estimator (Ŝ_N)**

전체 plug-in 추정값을 블록별 추정값의 두 분위수 사이로 제한한다.

![Proposed estimator](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%20%5Cwidehat%7BS%7D_N%3D%5Cmin%5C%21%5Cleft%5C%7B%5Cmax%5C%21%5Cleft%5C%7B%5Cwidehat%7BT%7D_N%2C%5Cwidehat%7BQ%7D%28%5Cbeta_1%29%5Cright%5C%7D%2C%5Cwidehat%7BQ%7D%28%5Cbeta_2%29%5Cright%5C%7D)

논문에서 일반적으로 선택한 설정은 **β₁ = 0.5, β₂ = 0.6**이며, **m ≈ 250** 부근의 블록 크기에 맞춘 선택이다. Figure 5의 정규·로그정규 실험에서는 블록 크기를 **m = 125**로 바꾼다.

![Error probability](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%20p_N%28%5Cdelta%29%3DP%5C%21%5Cleft%28%5Cleft%7C%5Cwidehat%7BR%7D_N-%5Coperatorname%7BES%7D_%7B%5Calpha%7D%28X%29%5Cright%7C%5Cgeq%5Cdelta%5Cright%29)

여기서 R̂_N은 비교할 추정량이고, δ는 **허용 오차의 기준값**이다. 위 확률이 작을수록 실제 ES에서 δ 이상 벗어날 위험이 작다. 논문 그림의 각 점은 **100만 회(10⁶)의 반복실험**으로 추정한 오차확률이다.

> **실험 설정에서 주의할 점**
>
> 저자들은 α = 0.01로 바꾸면 비슷한 수준의 결과를 얻기 위해 대략 10배의 데이터가 필요할 것으로 예상한다. 이는 실험을 해석하기 위한 저자들의 전망이며, 모든 분포와 조건에서 성립하는 정리로 제시된 수치는 아니다. 현재 실험은 금융 실무에서 충분한 데이터가 확보된 **이상화된 환경**을 가정한다.
{: .prompt-info }

---

## 5.2 Figure 2 — Pareto Distribution: 큰 오차와 작은 오차

첫 실험에서는 **Pareto λ = 2.1, λ = 3.5** 두 분포를 비교한다. λ가 2에 가까우면 두꺼운 꼬리의 영향이 강해지고, 매우 큰 손실이 표본에 나타날 가능성이 상대적으로 높아진다.

![Pareto survival](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%20P%28X%5Cgeq%20x%29%3Dx%5E%7B-%5Clambda%7D%5Cqquad%28x%5Cgeq%201%29)

Figure 2에서는 **m = 250**, **β₁ = 0.5**, **β₂ = 0.6**을 사용한다. 그래프의 세로축은 오차확률이며, 가로축은 데이터 개수 N이다.

![Figure 2 Pareto estimators](/assets/images/es_part5_figure2.png)

*Figure 2. Pareto λ = 2.1 및 λ = 3.5에서 오차확률 비교. 출처: Bartl & Eckstein (2026), Figure 2, p. 28.*

### 큰 오차를 피하려는 경우

왼쪽 열은 **δ = 1**인 실험이다. ES 추정값이 참값으로부터 1 이상 벗어나는 사건을 본다. 여기서는 **proposed estimator와 median of blocks가 plug-in보다 빠르게 오차확률을 낮추는 경향**이 관찰된다.

논문의 설명에 따르면 제안 추정량과 블록 중앙값 방식의 집중확률은 지수적 감소와 연결되고, Pareto 분포에서 plug-in 추정량은 4편 Lemma 4.1에서 확인한 다항식 하한의 영향을 받는다.

![Rates comparison](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%20p_N%28%5Cdelta%29%5Capprox%20C_%5Cdelta%20N%5E%7B-%28%5Clambda-1%29%7D%5Cquad%5Ctext%7Bversus%7D%5Cquad%5Cexp%28-c_%5Cdelta%20N%29)

위 식은 각 감소 형태를 비교하는 **개략적 표현**이다. 그래프의 모든 지점이 일정한 상수 하나로 정확히 표현된다는 뜻은 아니다.

### 작은 오차를 구분하려는 경우

오른쪽 열은 **δ = 0.25 또는 δ = 0.10**인 설정이다. 이 영역에서는 proposed estimator의 곡선이 plug-in과 가까워진다. 극단적인 오차가 발생하지 않는 일반적인 표본에서는 전체 데이터를 사용하는 추정 방식의 장점을 상당 부분 유지한다는 해석과 연결된다.


>
> 오차 기준 δ에 따라 추정량들의 상대적 성능이 달라진다. 제안 추정량은 큰 오차를 억제하는 영역에서 블록 중앙값의 장점을 활용하고, 작은 오차를 평가하는 영역에서는 plug-in의 정확도에 가까운 결과를 얻는다.
{: .prompt-info }

---

## 5.3 Figure 3 — Student's t Distribution: 다른 Heavy Tail에서도 나타날까

Pareto 분포에서 얻은 결과가 특정한 분포에만 의존하는지 확인하기 위해, 저자들은 **Student's t 분포의 자유도 ν = 2.2**를 사용한다. 자유도가 작은 t 분포는 정규분포에 비해 두꺼운 꼬리를 갖는다. 본 실험에서도 블록 크기는 **m = 250**이다.

![Figure 3 Student t distribution](/assets/images/es_part5_figure3.png)

*Figure 3. Student-t ν = 2.2에 대한 추정량별 오차확률. 출처: Bartl & Eckstein (2026), Figure 3, p. 29.*

Figure 3은 **δ = 0.25와 δ = 0.10**을 비교한다. Figure 2에서 관찰한 특징이 Student-t 분포에서도 확인된다. 특히 원문은 이 결과를 통해 제안 방법의 장점이 Pareto라는 하나의 분포에 국한되지 않음을 설명한다.

여기서 Student-t는 양수 손실만 갖는 분포가 아니다. 그러나 손실 변수는 음수 값을 가질 수도 있으므로 ES의 정의와 알고리즘을 그대로 적용할 수 있다. 우리가 관심을 갖는 부분은 **상단 꼬리의 큰 손실값**이다.


> 제안 추정량의 성능을 평가할 때 꼬리모수 하나만 바꾸는 것으로 끝내지 않고, **분포의 형태까지 바꾸어 비교**한다는 점이 중요하다. 분포별로 추정값의 편향과 오차확률의 크기가 달라질 수 있기 때문이다.
{: .prompt-info }

---

## 5.4 Figure 4 — Trimmed Estimator와 비교

이번에는 기존 연구의 **Trimmed Estimator**와도 비교한다. 원문은 참고문헌 [27]에서 제시한 계열의 추정량을 사용하며, 그림에서는 **Hill estimator**, **proposed estimator**를 비교한다.

![Figure 4 trimmed estimator comparison](/assets/images/es_part5_figure4.png)

*Figure 4. Pareto λ = 2.1에서 trimmed 계열과 제안 추정량의 오차확률 비교. 출처: Bartl & Eckstein (2026), Figure 4, p. 29.*

Trimmed 방식은 표본의 극단적인 상위 관측값 일부를 제거하고 추정한다. 극단값의 영향을 줄이는 데 도움이 될 수 있지만, 제거하는 데이터가 실제 분포의 꼬리 정보를 포함할 수도 있다.

논문에서는 **가장 큰 ⌊0.25N^(1/3)⌋개 관측값**을 제거한다. 선행연구의 실험에서 사용한 **⌊0.25N^(2/3)⌋개**는 이 논문의 실험에서는 제거 정도가 지나치게 강하다고 판단해 변경했다.

그림은 **δ = 1**에서 표준 trimmed 버전이 유리한 특성을 보이는 반면, **δ = 0.5**에서는 그 장점이 줄어드는 모습을 보여준다. 편향 보정을 적용한 버전은 상대적으로 반대의 경향을 보인다.

> **실험 설계의 시사점**
>
> trimming 개수를 어떻게 선택하는지에 따라 결과가 달라진다. 서로 다른 추정량을 비교할 때는 **제거 비율, 블록 크기, 오차 기준 δ**를 함께 기록해야 한다. 이 설정을 생략하면 성능 차이의 원인을 해석하기 어렵다.
{: .prompt-info }

---

## 5.5 Figure 5 — Normal & Log-normal: Light Tail

앞선 실험들은 두꺼운 꼬리를 가진 분포를 주로 다뤘다. 이제 상대적으로 안정적인 손실분포에서 제안 추정량이 어떤 성능을 보이는지 확인한다. 저자들은 **표준정규분포와 표준 로그정규분포**를 사용하고, 블록 크기를 **m = 125**로 둔다.

![Figure 5 normal lognormal comparison](/assets/images/es_part5_figure5.png)

*Figure 5. 표준정규분포와 표준 로그정규분포에 대한 오차확률 및 추정값 히스토그램. 출처: Bartl & Eckstein (2026), Figure 5, p. 30.*

**왼쪽 그래프**는 N에 따른 오차확률을 보여주고, **오른쪽 히스토그램**은 N = 4,750일 때의 추정값 분포를 보여준다. 정상적인 조건에서 plug-in이 우수한 정확도를 보이며, proposed estimator도 거의 비슷한 결과를 얻는다. 반면 **median of blocks는 음의 편향이 눈에 띈다.**

특히 Figure 5는 1편의 Remark 1.6에서 설명한 **블록 중앙값의 편향 문제**를 시각적으로 확인하는 부분이다.

> **Median of Blocks가 아래쪽으로 치우치는 이유**
>
> 블록 크기가 제한되면 각 블록의 상단 꼬리 평균이 전체 모집단의 ES보다 작게 계산되는 경우가 반복될 수 있다. 그런 블록 추정값들의 중앙값을 사용하면 음의 편향이 나타날 수 있다. 그림의 히스토그램은 이 효과를 보여주며, 분포 및 블록 크기에 따라 편향의 정도는 달라진다.
{: .prompt-info }

---

## 5.6 Figure 6 — Adversarial Contamination: 데이터 3개만 바꾸었을 때

2편에서 살펴본 **Theorem 1.8**은 제안 추정량의 적대적 강건성을 다룬다. Figure 6은 이 결과를 실험으로 확인한다. 원래 데이터는 **Pareto λ = 2.2**, **N = 3,250**에서 생성한다.

그리고 전체 데이터 중 **처음 세 관측값**을 다음 규칙으로 변조한다.

![Adversarial contamination](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%20%5Cwidetilde%7BX%7D_i%3D%5Cmax%5C%21%5Cleft%5C%7BX_i%2CU_i%5Cright%5C%7D%2C%5Cquad%20U_i%5Csim%5Cmathcal%7BN%7D%285%2C250%5E2%29)

세 개의 U_i는 서로 독립이며 원래 데이터와도 독립이다. 정규분포의 평균은 5, 표준편차는 250이다. 각 값을 기존 관측값과 U_i 중 더 큰 값으로 교체한다.

![Figure 6 adversarial contamination](/assets/images/es_part5_figure6.png)

*Figure 6. 변조 전후의 추정값 히스토그램. 출처: Bartl & Eckstein (2026), Figure 6, p. 31.*

왼쪽은 변조 전, 오른쪽은 변조 후의 추정값 분포다. 그림의 파란 세로선은 **실제 ES−1, 실제 ES, 실제 ES+1**을 나타낸다. 단 세 개의 데이터가 변조된 후 **plug-in의 히스토그램이 크게 왜곡되는 반면 proposed estimator는 원래 분포에 가까운 모습을 유지**한다.

세 데이터는 전체 3,250개의 약 0.092%에 해당한다. 비율은 매우 작지만, plug-in 추정량에는 상당한 영향을 줄 수 있다.

> **Theorem 1.8과의 연결**
>
> 이 실험은 데이터 변조에 대한 추정량의 민감성을 직접 보여준다. Theorem 1.8의 엄밀한 보장은 변조된 관측값 수와 분포에 관한 조건에 의존하므로, Figure 6의 한 실험 결과만으로 모든 변조 상황에 대한 성능을 보장할 수는 없다.
{: .prompt-info }

---

## 5.7 Figure 7 — Infinite Variance: 이론 가정 밖

마지막 실험은 **Pareto λ = 1.5**이다. 파레토 분포에서 λ가 2 이하이면 분산이 무한대다. λ = 1.5에서는 ES 자체는 유한하지만, 논문의 Theorem 1.4에서 요구하는 **유한한 ES 점근 분산** 조건이 충족되지 않는다.

![Infinite variance](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%20%5Clambda%3D1.5%5Cquad%5CLongrightarrow%5Cquad%5Coperatorname%7BVar%7D%28X%29%3D%5Cinfty)

![Figure 7 infinite variance example](/assets/images/es_part5_figure7.png)

*Figure 7. Pareto λ = 1.5에서 추정량별 오차확률. 출처: Bartl & Eckstein (2026), Figure 7, p. 32.*

이 실험에서 proposed estimator는 plug-in보다 훨씬 안정적인 결과를 보이며, 논문은 특히 **이론적 보장이 적용되지 않는 경우에도 제안 방법이 좋은 실험 성능을 나타낼 수 있음**을 강조한다.


> Figure 7은 실험적으로 관찰된 현상이다. 분산이 무한한 모든 분포에서 동일한 오차 상한이 성립한다는 수학적 증명으로 해석해서는 안 된다. 이러한 조건에서 어떤 보장이 가능한지는 별도의 연구 문제다.
{: .prompt-info }

---

## 5.8 Python — 핵심 추정량 직접 비교

이번 글의 그림은 논문에서 발췌했다. 원문은 각 오차확률을 **100만 회 반복실험**으로 계산했으므로, 모든 설정을 그대로 재현하려면 상당한 계산이 필요하다. 아래 코드는 **원문의 핵심 추정량과 비교 논리를 이해하기 위한 소규모 실습**이다. Figure 2~7의 곡선을 그대로 재현한 결과로 해석하면 안 된다.

### Step 1. 세 가지 ES 추정량 구현

```python
import numpy as np

def plugin_es(x, alpha=0.1):
    """Empirical ES, including the fractional boundary weight."""
    x = np.sort(np.asarray(x, dtype=float))[::-1]
    if not len(x) or not (0 < alpha < 1):
        raise ValueError("Nonempty x and 0 < alpha < 1 are required")
    tail_count = alpha * len(x)
    k = int(np.floor(tail_count))
    result = x[:k].sum()
    if k < len(x):
        result += (tail_count - k) * x[k]
    return float(result / tail_count)

def three_estimators(x, alpha=0.1, m=250,
                     beta1=0.5, beta2=0.6):
    x = np.asarray(x, dtype=float)
    n_blocks = len(x) // m
    if n_blocks < 2:
        raise ValueError("At least 2 complete blocks are required")

    full = plugin_es(x, alpha)
    blocks = x[:n_blocks * m].reshape(n_blocks, m)
    values = np.array([plugin_es(b, alpha) for b in blocks])
    lower, upper = np.quantile(
        values, [beta1, beta2], method="linear"
    )
    return np.array([
        full,
        np.median(values),
        np.clip(full, lower, upper)
    ])
```

첫 번째 값은 **plug-in**, 두 번째 값은 **median of blocks**, 세 번째 값은 **proposed estimator**다. 전체 plug-in에는 모든 N개 관측값이 들어가며, 블록 계산에서는 마지막에 남는 관측값이 제외될 수 있다. 이 방식은 1편에서 설명한 Definition 1.3의 구조를 따른다.

### Step 2. Pareto 분포에서 오차확률 비교

```python
import matplotlib.pyplot as plt

rng = np.random.default_rng(2026)
alpha, lam = 0.1, 2.2
m, reps = 250, 1000
sizes = [1000, 2000, 3250, 5000]
threshold = 0.5
true_es = (lam / (lam - 1)) * alpha**(-1 / lam)

errors = np.zeros((len(sizes), 3))
for i, N in enumerate(sizes):
    estimates = np.empty((reps, 3))
    for j in range(reps):
        x = rng.pareto(lam, size=N) + 1
        estimates[j] = three_estimators(x, alpha, m)
    errors[i] = np.mean(
        np.abs(estimates - true_es) >= threshold, axis=0
    )

labels = ["Plug-in", "Median of blocks", "Proposed"]
for col, label in enumerate(labels):
    plt.plot(sizes, errors[:, col], marker="o", label=label)
plt.xlabel("Sample size N")
plt.ylabel("Estimated P(absolute error >= 0.5)")
plt.title("Illustrative Pareto ES comparison")
plt.legend()
plt.tight_layout()
plt.show()
```

![Python Monte Carlo experiment](/assets/images/es_part5_python_experiment.png)

*Figure 8. Python 소규모 재현 실험: Pareto 분포에서 세 가지 ES 추정량의 오차확률 비교 (α = 0.1, λ = 2.2, m = 250, 1,000회 반복).*

이 코드는 α = 0.1, λ = 2.2, m = 250에서 **실험 반복 횟수 1,000회**로 비교한다. 원문 Figure 2는 다른 λ와 δ, 100만 회 반복을 사용했다. 따라서 **곡선의 수치와 순위가 원문과 똑같이 나오리라고 기대하면 안 된다.** 반복 횟수가 작으면 아주 낮은 오차확률은 정확하게 측정하기 어렵다.

![Pareto theoretical ES](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%20%5Coperatorname%7BES%7D_%7B%5Calpha%7D%28X%29%3D%5Cfrac%7B%5Clambda%7D%7B%5Clambda-1%7D%5Calpha%5E%7B-1%2F%5Clambda%7D)

위의 이론적 ES는 최소 손실 x₀ = 1인 파레토 분포에서 사용할 수 있다.

### Step 3. 데이터 변조 전후 비교

```python
rng = np.random.default_rng(7)
N, alpha, m = 3250, 0.1, 250
x = rng.pareto(2.2, size=N) + 1

baseline = three_estimators(x, alpha, m)
corrupted = x.copy()
U = rng.normal(5, 250, size=3)
corrupted[:3] = np.maximum(corrupted[:3], U)
after = three_estimators(corrupted, alpha, m)

for name, a, b in zip(
    ["Plug-in", "Median of blocks", "Proposed"],
    baseline, after
):
    print(f"{name:18s}: before={a:.4f}, after={b:.4f}")
```

이 코드는 Figure 6의 변조 규칙을 한 번 적용해 추정값 변화를 비교한다. Figure 6처럼 **히스토그램의 형태까지 재현**하려면 데이터를 새로 생성하는 과정을 수많은 횟수 반복해야 한다.

---

## 5.9 Numerical Experiments

Section 5는 서로 다른 추정량의 성능을 **하나의 숫자로 평가하기 어렵다**는 사실을 잘 보여준다. 오차 기준 δ, 분포의 꼬리 형태, 블록 크기, 데이터 변조 여부에 따라 비교 결과가 달라진다.

| 실험 | 분포 및 설정 | 논문에서 확인한 핵심 결과 |
|---|---|---|
| Figure 2 | Pareto λ = 2.1, 3.5 | 큰 오차에서는 강건 추정량의 장점, 작은 오차에서는 plug-in에 가까운 정확도 |
| Figure 3 | Student-t ν = 2.2 | 다른 heavy-tailed 분포에서도 유사한 경향 |
| Figure 4 | Pareto λ = 2.1 | Trimmed 추정량은 보정 여부와 δ에 따라 성능 변화 |
| Figure 5 | Normal, Log-normal | Plug-in과 proposed는 비슷한 성능, median of blocks는 음의 편향 |
| Figure 6 | Pareto λ = 2.2, 3개 변조 | Proposed estimator의 데이터 변조 강건성 |
| Figure 7 | Pareto λ = 1.5 | 무한분산에서도 좋은 실험 성능 관찰, 정리의 보장 범위 밖 |

특히 Figure 2와 Figure 5를 함께 보면 제안 추정량의 설계 의도가 선명해진다. **일반적인 표본에서는 전체 데이터를 사용하는 추정의 정확도를 상당 부분 살리고, 극단적인 표본에서는 블록 분위수가 추정값의 폭주를 제어한다.**

Figure 6과 7은 이론적 조건을 어디까지 신뢰할 수 있는지 생각하게 만든다. 표본 일부가 변조된 환경에서 안정성을 확인할 수 있지만, 무한분산 실험에서 관찰한 좋은 결과를 Theorem 1.4의 보장으로 확대해서 해석할 수는 없다.

다음 6편에서는 **i.i.d. 가정**을 완화한다. 실제 금융 시계열에는 시간에 따른 종속성이 나타날 수 있으며, 원문 Section 6에서는 **β-mixing** 조건을 통해 오차 보장을 확장하는 방법을 다룬다.

---

> Bartl, D., & Eckstein, S. (2026). *Optimal Nonparametric Estimation of the Expected Shortfall Risk*. Section 5, pp. 26–32, Figures 2–7. [arXiv:2405.00357v2](https://arxiv.org/abs/2405.00357)
