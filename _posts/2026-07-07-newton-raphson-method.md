---
title: "Newton-Raphson Method"
date: 2026-07-07 20:00:00 +0900
categories: [Numerical Analysis]
tags: [Newton-Raphson Method, Newton Method, Numerical Analysis, Python, MATLAB, C++]
math: true
---

# Newton-Raphson Method


$$
f(x)=0
$$
의 해를 직접 구할 수 있으면 가장 좋겠지만, 실제 계산에서는 그렇지 않은 경우가 많다. 다항방정식이어도 차수가 높아지면 해를 깔끔하게 쓰기 어렵고, 삼각함수나 지수함수가 섞인 초월방정식에서는 손으로 해를 구하기가 더 어려워진다.

이럴 때 필요한 것이 근사해를 구하는 수치해법이다. 앞 장에서 설명했었던 이분법은 근이 들어 있는 구간을 잡고 그 구간을 계속 반으로 줄이는 방법이었다. 

반면 뉴턴-랩슨법은 구간보다 한 점에서 출발한다. 현재 점에서 함수의 접선을 그리고, 그 접선이 $x$축과 만나는 지점을 다음 근사값으로 삼는다.

뉴턴-랩슨법은 뉴턴법이라고도 부른다. 미분 정보를 사용하기 때문에 수렴이 빠른 편이고, 초기값이 해에 충분히 가까우면 매우 적은 반복으로도 정확한 근사해를 얻을 수 있다. 

다만 함수값만 보는 방법이 아니므로, 도함수가 존재해야 하고 도함수 값이 0에 가까워지는 경우에는 계산이 불안정해질 수 있다.


| 항목 | 이분법 | Newton 방법 |
|---|---|---|
| 수렴속도 | 느림 | 빠름 |
| 수렴차수 | 1 | 2 |
| 수렴까지에 필요한 점의 수 | 2점 | 1점 |
| 수렴조건 | 양 끝점 함수값의 부호가 다름 | 실제 해의 초기값이 가까이 있어야 함 |
| 함수의 미분값 | 필요 없음 | 필요함 |
| 평가 비용 | 수렴속도는 같음, 수렴차수 1 | M이 커짐에 따라 수렴속도가 이분법보다 더 느려짐 |
| 2차원 이상으로의 확대 | 불가능 | 가능 |

---

## Basic Principle of the Newton-Raphson Method

현재 근사값을 $x_n$이라고 하자. 이 점에서의 함수값은 $f(x_n)$이고, 
 접선의 기울기는 $f'(x_n)$이다.

점 $(x_n, f(x_n))$에서의 접선의 방정식은 $$
y = f'(x_n)(x-x_n)+f(x_n)
$$ 이다. 

뉴턴-랩슨법에서는 이 접선이 $x$축과 만나는 점을 다음 근사값 $x_{n+1}$로 정한다. $x$축과 만난다는 것은 $y=0$이라는 뜻이므로,

$$
0 = f'(x_n)(x_{n+1}-x_n)+f(x_n)
$$ 이다. 이를 $x_{n+1}$에 대해 정리하면

$$
x_{n+1}=x_n-\frac{f(x_n)}{f'(x_n)}
$$ 을 얻는다.

이 식이 뉴턴-랩슨법의 핵심 점화식이다.


> #### Newton-Raphson Iteration
>
> $$
> x_{n+1}=x_n-\frac{f(x_n)}{f'(x_n)}
> $$
>
> 현재 점 $x_n$에서 접선을 긋고, 그 접선의 $x$절편을 다음 근사값 $x_{n+1}$로 잡는다.
{: .prompt-warning }

뉴턴-랩슨법은 이분법처럼 양 끝점의 부호가 달라야 한다는 조건에서 비교적 자유롭다. 시작점 하나만 정하면 계산을 시작할 수 있다. 대신 그 시작점이 어디인지에 따라 수렴 여부가 달라진다.

---

## Algorithm

뉴턴-랩슨법의 계산 과정은 단순하다. 초기값을 정하고, 점화식을 반복하고, 두 근사값의 차이가 충분히 작아지면 멈춘다.

> #### Newton-Raphson Algorithm
>
> 1. 초기 추정값 $x_0$을 정한다.
> 2. 다음 식으로 새로운 근사값을 계산한다.
>
>    $$
>    x_{n+1}=x_n-\frac{f(x_n)}{f'(x_n)}
>    $$
>
> 3. 다음 조건 중 하나를 만족하면 반복을 멈춘다.
>
>    $$
>    |x_{n+1}-x_n|<\mathrm{TOL}
>    $$
>
>    또는
>
>    $$
>    |f(x_{n+1})|<\mathrm{TOL}
>    $$
>
> 4. 조건을 만족하지 않으면 $x_{n+1}$을 새로운 기준점으로 삼아 다시 계산한다.
{: .prompt-warning }

여기서 $\mathrm{TOL}$은 tolerance, 곧 허용 오차를 뜻한다. 계산을 무한히 반복할 수는 없기 때문에 어느 정도 가까워졌을 때 멈출지를 정해둔다.

---

## Calculation Process Through Tangents

뉴턴-랩슨법은 공식만 보면 간단하지만, 실제 의미는 접선 그림을 보면 훨씬 잘 보인다.

현재 점이 $x_n$이면 곡선 위의 점은

$$
(x_n, f(x_n))
$$

이다. 이 점에서 접선을 그은 뒤, 그 접선이 $x$축과 만나는 점을 $x_{n+1}$로 둔다. 다시 $x_{n+1}$에서 함수값을 계산하고, 또 접선을 그린다. 이 과정을 반복하면 점들이 근 쪽으로 이동한다.

예를 들어

$$
f(x)=x^2-2
$$

의 근은

$$
x=\sqrt{2}
$$

이다. 뉴턴-랩슨법으로 이 값을 구하면 $1.41421\cdots$에 빠르게 가까워진다.

아래 Python 코드는 $f(x)=x^2-2$에 뉴턴-랩슨법을 적용하고, 각 반복 단계에서의 접선까지 함께 시각화한다.

```python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

# =========================
# 한글 폰트 설정
# =========================
font_path = "C:/Windows/Fonts/malgun.ttf"
font_prop = fm.FontProperties(fname=font_path)

# =========================
# 함수 및 도함수 정의
# =========================
def f(x):
    return x**2 - 2

def df(x):
    return 2 * x

# =========================
# 접선을 시각화하는 뉴턴법
# =========================
def newton_method_with_tangents(x0, tol=1e-4, max_iter=10):
    steps = [(x0, f(x0))]
    tangents = []

    for _ in range(max_iter):
        fx = f(x0)
        dfx = df(x0)

        if dfx == 0:
            break

        x1 = x0 - fx / dfx
        steps.append((x1, f(x1)))

        # 접선 계산 : y = slope*x + intercept
        slope = dfx
        intercept = fx - slope * x0

        x_range = np.linspace(x0 - 1, x0 + 1, 2)
        y_range = slope * x_range + intercept

        tangents.append((x_range, y_range))

        if abs(x1 - x0) < tol:
            break

        x0 = x1

    return x1, steps, tangents

# =========================
# 실행
# =========================
x0 = 1.0
root, steps, tangents = newton_method_with_tangents(x0)

print(f"뉴턴법으로 구한 해: {root}")

# =========================
# 시각화
# =========================
x_vals = np.linspace(0, 2, 400)
y_vals = f(x_vals)

plt.figure(figsize=(8, 5))

plt.plot(
    x_vals,
    y_vals,
    label="함수 $f(x)=x^2-2$",
    color="blue"
)

plt.axhline(
    0,
    color="black",
    linewidth=0.8
)

plt.axvline(
    root,
    color="purple",
    linestyle="--",
    label=f"뉴턴법 근사 해 = {root:.5f}"
)

# 반복 과정의 점
for x, y in steps:
    plt.plot(x, y, "mo", markersize=4)

# 접선
for x_tangent, y_tangent in tangents:
    plt.plot(
        x_tangent,
        y_tangent,
        "r--",
        linewidth=1
    )

plt.title(
    "뉴턴법 반복 과정과 접선 시각화",
    fontproperties=font_prop
)

plt.xlabel(
    "x",
    fontproperties=font_prop
)

plt.ylabel(
    "f(x)",
    fontproperties=font_prop
)

plt.legend(prop=font_prop)
plt.grid(True)
plt.tight_layout()

plt.savefig("newton_method_with_tangents.png")
plt.show()
```

![뉴턴법 접선 시각화](/assets/img/posts/newton_method_with_tangents.png)

수직 점선은 뉴턴법으로 얻은 근사해의 위치를 나타낸다. 곡선 위의 점들은 각 반복 단계에서의 추정값 $x_n$과 함수값 $f(x_n)$이다. 붉은 점선은 각 점에서 그린 접선이다. 

접선의 $x$절편이 다음 근사값 $x_{n+1}$이 된다.

$x_0=1$에서 시작하면 $x^2-2=0$의 양의 해인 $\sqrt{2}$ 쪽으로 빠르게 이동한다. 

이분법이 구간을 하나씩 줄여가는 방식이라면, 뉴턴-랩슨법은 접선의 방향을 이용해 다음 지점을 한 번에 크게 이동시키는 방식이다.

---

## Example

다음 방정식을 생각하자.

$$
x^3-5x+3=0
$$

여기서

$$
f(x)=x^3-5x+3
$$

이고 도함수는

$$
f'(x)=3x^2-5
$$

이다. 뉴턴-랩슨법을 적용하려면 초기값 $x_0$을 정하고

$$
x_{n+1}=x_n-\frac{x_n^3-5x_n+3}{3x_n^2-5}
$$

를 반복하면 된다.

아래 코드는 `scipy.optimize.root_scalar`의 Newton 방법을 사용한 예제이다.

```python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from scipy.optimize import root_scalar

# =========================
# 한글 폰트 설정
# =========================
font_path = "C:/Windows/Fonts/malgun.ttf"
font_prop = fm.FontProperties(fname=font_path)

# =========================
# 함수 및 도함수 정의
# =========================
def f(x):
    return x**3 - 5*x + 3

def df(x):
    return 3*x**2 - 5

# =========================
# 뉴턴-랩슨 방법
# =========================
result = root_scalar(
    f,
    x0=1,
    fprime=df,
    method="newton"
)

print("뉴턴랩슨으로 구한 해:", result.root)

# =========================
# 그래프 데이터 생성
# =========================
x_vals = np.linspace(-3, 3, 400)
y_vals = f(x_vals)

# =========================
# 그래프 그리기
# =========================
plt.figure(figsize=(8, 5))

plt.plot(
    x_vals,
    y_vals,
    label="함수 $f(x)=x^3-5x+3$",
    color="blue"
)

plt.axhline(
    0,
    color="black",
    linewidth=0.8
)

plt.axvline(
    result.root,
    color="purple",
    linestyle="--",
    label=f"뉴턴법 근 = {result.root:.5f}"
)

plt.plot(
    result.root,
    0,
    "ro",
    label="근 위치"
)

plt.title(
    "뉴턴-랩슨을 이용한 비선형 방정식 풀이",
    fontproperties=font_prop
)

plt.xlabel(
    "x",
    fontproperties=font_prop
)

plt.ylabel(
    "f(x)",
    fontproperties=font_prop
)

plt.legend(prop=font_prop)
plt.grid(True)
plt.tight_layout()

plt.savefig("newton_raphson_cubic.png")
plt.show()
```

![뉴턴-랩슨 삼차방정식 예제](/assets/img/posts/newton_raphson_cubic.png)

## Convergence Rate of the Newton Method

뉴턴법은 적은 반복 횟수로 정확한 결과를 산출하며, 비선형 방정식의 근을 구하는 데 수렴 속도가 가장 빠른 방법 중 하나이다.

뉴턴법의 수렴 속도는 테일러 다항식으로 확인할 수 있다. 방정식

$$
f(x)=0
$$

의 근을 $\alpha$라고 하자. 그러면

$$
f(\alpha)=0
$$ 이다. $x_n$ 근방에서 $f(\alpha)$를 테일러 전개하면, 어떤 $\xi$에 대하여

$$
0=f(\alpha)=f(x_n)+(\alpha-x_n)f'(x_n)+\frac{(\alpha-x_n)^2}{2}f''(\xi)
$$ 가 된다.

여기에서 $f'(x_n)$으로 나누고 정리하면

$$
x_n-\frac{f(x_n)}{f'(x_n)}-\alpha = \frac{(\alpha-x_n)^2 f''(\xi)}{2f'(x_n)}
$$ 이다.

뉴턴법의 점화식은

$$
x_{n+1}=x_n-\frac{f(x_n)}{f'(x_n)}
$$

이므로 위 식은 다음과 같이 쓸 수 있다.

$$
x_{n+1}-\alpha = \frac{(\alpha-x_n)^2 f''(\xi)}{2f'(x_n)}
$$

오차를

$$
E_n=x_n-\alpha
$$

라고 두면,

$$
E_{n+1} = \frac{f''(\xi)}{2f'(x_n)}E_n^2
$$

의 형태가 된다. 따라서

$$
\frac{|E_{n+1}|}{|E_n|^2} = \left|\frac{f''(\xi)}{2f'(x_n)}\right|
$$

이다.

반복이 진행되어 $x_n$이 근 $\alpha$에 가까워지면 $\xi$도 $\alpha$에 가까워진다. 따라서

$$
\lim_{n\to\infty}
\frac{|E_{n+1}|}{|E_n|^2} = \left|\frac{f''(\alpha)}{2f'(\alpha)}\right|
$$

가 된다.

이 값이 유한하고 0이 아니면 뉴턴법은 2차 수렴한다고 말한다. 2차 수렴은 오차가 단순히 일정한 비율로 줄어드는 것이 아니라, 이전 오차의 제곱에 비례하여 줄어드는 것을 의미한다. 그래서 초기값이 근에 충분히 가까우면 뉴턴법은 매우 빠르게 수렴한다.

---

## Convergence Conditions and Precautions

뉴턴-랩슨법은 빠르지만, 항상 안정적인 방법은 아니다. 함수가 미분 가능해야 하고, 반복 과정에서 도함수 값이 0이 되거나 0에 가까워지면 계산이 크게 흔들릴 수 있다.

> #### 뉴턴-랩슨법이 잘 작동하기 위한 조건
>
> - $f(x)$가 해 근처에서 미분 가능해야 한다.
> - $f'(x)$가 해 근처에서 0에 너무 가까워지지 않아야 한다.
> - 초기값 $x_0$이 해에 어느 정도 가까워야 한다.
> - 반복 과정이 다른 근이나 발산 방향으로 이동하지 않아야 한다.
{: .prompt-warning }

초기값이 해에서 멀리 떨어져 있으면 접선이 엉뚱한 방향으로 향할 수 있다. 이 경우 수열 $\{x_n\}$이 발산하거나, 다른 근으로 수렴하거나, 특정 값들 사이를 반복할 수도 있다.

이분법은 수렴 속도는 느리지만, $f(a)f(b)<0$인 구간을 잡으면 근이 존재한다는 보장을 가지고 시작한다. 뉴턴-랩슨법은 그런 보장보다 속도를 택한 방법에 가깝다.

그래서 실제 계산에서는 처음부터 뉴턴법만 쓰기보다, 이분법으로 근이 있는 구간을 어느 정도 좁힌 다음 뉴턴법을 적용하는 방식도 자주 사용한다.


| 방법 | 필요한 정보 | 시작 조건 | 수렴 속도 | 특징 |
|---|---|---|---|---|
| 이분법 | 함수값의 부호 | $f(a)f(b)<0$인 구간 | 느림 | 안정적이고 근의 존재를 보장하기 쉽다 |
| 할선법 | 두 점의 함수값 | 초기값 두 개 | 이분법보다 빠른 편 | 도함수 없이 직선을 이용한다 |
| 뉴턴-랩슨법 | 함수값과 도함수 | 초기값 한 개 | 매우 빠른 편 | 도함수가 필요하고 초기값에 민감하다 |

---

## MATLAB Code

다음은 뉴턴-랩슨법을 MATLAB 함수로 구현한 코드이다.

```matlab
function [x, iter] = newton(f, df, x0, TOL, MaxIter)

% Newton 방법에 의해 f(x)=0의 해를 구함
% x0      : 초기 추측값
% TOL     : 오차한계
% MaxIter : 최대 반복횟수
% f       : 함수 이름
% df      : f의 도함수

if nargin < 5
    MaxIter = 100;
end

if nargin < 4
    TOL = 1.e-5;
end

x = x0;

for k = 1:MaxIter
    fval = feval(f, x0);
    dfval = feval(df, x0);

    x = x0 - fval / dfval;
    err = abs(x - x0);

    if err < TOL
        break;
    else
        x0 = x;
    end
end

iter = k;

end
```

이 코드는 초기값 $x_0$에서 시작하여

$$
x_{n+1}=x_n-\frac{f(x_n)}{f'(x_n)}
$$

을 반복한다. `err = abs(x - x0)`는 현재 근사값과 이전 근사값의 차이를 의미한다. 이 값이 `TOL`보다 작아지면 반복을 멈춘다.

![MATLAB 뉴턴-랩슨 함수 코드](/assets/img/posts/Matlab1.png)
---

## Example
다음은

$$
f(x)=x^5+3x-1
$$

의 근을 뉴턴법으로 구하는 MATLAB 코드이다. 도함수는

$$
f'(x)=5x^4+3
$$

이다.

```matlab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% newton.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; clf;

f = 'x^5+3*x-1';
df = '5*x^4+3';

n = 0;
x0 = 0.5;
tol = 10^(-10);
m = 1;

fprintf(' n          x_n          f(x_n)        x_n-x_(n-1)\n')

while m > tol
    n = n + 1;

    x = x0;
    y = eval(f);
    dy = eval(df);

    x1 = x0 - y / dy;
    m = abs(x1 - x0);

    x = x1;
    y = eval(f);

    fprintf('%2.0f   %2.4e   %2.4e   %2.4e \n', ...
            n, x1, y, x1 - x0)

    x0 = x1;
end
```

출력되는 값들은 반복 횟수, 현재 근사값, 현재 함수값, 이전 근사값과의 차이를 나타낸다. 반복이 진행될수록 $f(x_n)$이 0에 가까워지고, $x_n-x_{n-1}$의 크기도 작아진다.

![MATLAB 뉴턴-랩슨 함수 코드](/assets/img/posts/Matlab2.png)

## 비선형 연립방정식으로의 확장

뉴턴법은 일변수 방정식뿐 아니라 다변수 비선형 시스템에도 확장된다.

일변수에서는

$$
x_{n+1}=x_n-\frac{f(x_n)}{f'(x_n)}
$$

이었다. 다변수에서는 함수 $F$와 야코비안 행렬 $J$를 사용한다.

$$
p_{k+1}=p_k-J^{-1}F(p_k)
$$

여기서 $p_k$는 현재 근사 벡터이고, $F(p_k)$는 각 방정식의 함수값을 모은 벡터이다. $J$는 편미분값으로 이루어진 야코비안 행렬이다.

엄밀하게는 역행렬을 직접 구하기보다

$$
J(p_k)s_k = F(p_k)
$$

를 풀고,

$$
p_{k+1}=p_k-s_k
$$

처럼 계산하는 편이 더 안정적이다. 아래 코드는 원리를 보이기 위한 MATLAB 예제이다.

> #### 예제. 뉴턴 방법으로 비선형 연립방정식 풀기
>
> 다음 주어진 연립방정식을 뉴턴 방법으로 풀어보자.
>
> $$
> \begin{cases}
> x^2-\sin(y)+0.5\cos(z)=0.5,\\
> 3x-\cos(y)+\sin(z)=0,\\
> x^2+y^2+z^2=0.95.
> \end{cases}
> $$
>
> 초기 조건은 다음과 같이 둔다.
>
> $$
> x_0=0.5,\qquad y_0=0.5,\qquad z_0=0.5
> $$
>
> 이 문제는 미지수가 세 개인 비선형 연립방정식이므로, 일변수 뉴턴법의 점화식을 그대로 쓰지 않고 함수 벡터 $F$와 야코비안 행렬 $J$를 이용한다.
>
> $$
> p_{k+1}=p_k-J(p_k)^{-1}F(p_k)
> $$
>
> 여기서
>
> $$
> p_k=
> \begin{bmatrix}
> x_k\\
> y_k\\
> z_k
> \end{bmatrix}
> $$
>
> 이다.

```matlab
function F = f_func(x)

f1 = inline('x(1)^2-sin(x(2))+0.5*cos(x(3))-0.5','x');
f2 = inline('3*x(1)-cos(x(2))+sin(x(3))','x');
f3 = inline('x(1)^2+x(2)^2+x(3)^2-0.95','x');

F = [f1(x); f2(x); f3(x)];

end
```

```matlab
function J = j_func(x)

J = [ 2*x(1)          -cos(x(2))      -0.5*sin(x(3));
      3               sin(x(2))        cos(x(3));
      2*x(1)          2*x(2)           2*x(3) ];

end
```

```matlab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% newtonsys.m %%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc;

x = [0.5; 0.5; 0.5];

count = 0;
err = max(abs(f_func(x)));
tol = 1e-10;

fprintf('---------------------------\n');
fprintf('count      absolute\n');
fprintf('---------------------------\n');

while err > tol
    count = count + 1;

    x1 = x - inv(j_func(x)) * f_func(x);
    x = x1;

    err = max(abs(f_func(x)));

    fprintf('%d\t\t%f\n', count, err)
end

fprintf('---------------------------\n');
```
이 코드는 세 개의 미지수를 갖는 비선형 시스템에 뉴턴법을 적용한다.

![MATLAB 뉴턴-랩슨 함수 코드](/assets/img/posts/Matlab3.png)

 $F(x)$가 0에 가까워질 때까지 야코비안 행렬을 이용해 근사 벡터를 갱신한다.
 
---

## C++ Code

다음 C++ 코드는

$$
f(x)=1+\cos(3x)-x
$$

의 근을 뉴턴법으로 구한다. 도함수는

$$
f'(x)=-3\sin(3x)-1
$$

이다.

```cpp
#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

double f(double x) {
    return 1 + cos(3 * x) - x;
}

double df(double x) {
    return -3 * sin(3 * x) - 1;
}

int main() {
    double x_prev;
    double x = 1.0;
    double tol = 1e-12;
    int maxIter = 100;

    cout << fixed << setprecision(6);

    cout << "n\t"
         << "x_n\t\t"
         << "f(x_n)\t\t"
         << "|x_n-x_(n-1)|" << endl;

    cout << "0\t"
         << x << "\t"
         << f(x) << "\t"
         << "-" << endl;

    for (int n = 1; n <= maxIter; n++) {
        x_prev = x;

        x = x_prev - f(x_prev) / df(x_prev);

        cout << n << "\t"
             << x << "\t"
             << f(x) << "\t"
             << fabs(x - x_prev) << endl;

        if (fabs(x - x_prev) < tol) {
            break;
        }
    }

    return 0;
}
```

이 코드에서도 핵심은 같다. 현재값 `x_prev`에서 함수값과 도함수값을 계산하고,

```cpp
x = x_prev - f(x_prev) / df(x_prev);
```

를 통해 다음 근사값을 만든다.
![C++ 뉴턴-랩슨 함수 코드](/assets/img/posts/cpp1.png)
---

---

## Pros and Cons

뉴턴-랩슨법의 가장 큰 장점은 속도이다. 초기값이 해에 가까우면 수렴 속도가 매우 빠르다. 특히 단순근 근처에서는 보통 이차 수렴을 보인다. 이차 수렴이란 어느 정도 해에 가까워진 뒤에는 오차가 대략 제곱의 규모로 줄어드는 것을 말한다.

하지만 빠르다는 장점은 조건이 맞을 때의 이야기이다. 도함수가 0에 가까운 지점에서는

$$
\frac{f(x_n)}{f'(x_n)}
$$

이 지나치게 커질 수 있다. 그러면 다음 근사값이 해에서 멀리 벗어날 수 있다. 초기값이 좋지 않은 경우에도 같은 문제가 생긴다.

뉴턴-랩슨법은 접선의 방향을 믿고 이동하는 방법이다. 접선이 해 쪽을 향하면 빠르게 수렴하지만, 접선이 다른 방향을 향하면 오히려 계산이 망가진다.

> ####  Newton-Raphson Method
>
> - 수렴이 빠르다.
> - 도함수를 사용한다.
> - 초기값 하나에서 시작할 수 있다.
> - 초기값이 좋지 않으면 발산할 수 있다.
> - $f'(x_n)=0$이거나 0에 가까우면 계산이 불안정하다.
> - 이분법으로 구간을 먼저 좁힌 뒤 사용하면 더 안정적으로 쓸 수 있다.
{: .prompt-warning }

---

## References

- 이관수. 『공학도를 위한 수치해석』. 세화, 2014.
- 김학윤. 『수치해석 입문』. 영미디어, 1999.
- 권순학. 『수치해석 기초: MATLAB 활용』. 영남대학교 출판부, 2016.
- 이창영. 『수치해석의 기본』. 교문사(청문각), 2008.
- 박태희. 『MATLAB을 이용한 알기 쉬운 수치해석』. 생능출판, 2018.
