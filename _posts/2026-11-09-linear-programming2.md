---
title: "Optimization Theory : Linear Programming2"
date: 2026-11-09 11:00:00 +0100
categories:
  - Optimization Theory
tags:
  - Linear Programming
  - Simplex Algorithm
  - Simplex Tableau
  - Pivot
  - Phase I
  - Phase II
math: true
---

## 1. Geometrical Method and the Simplex Algorithm (p.3)

Linear programming의 최적해를 찾을 때 feasible set의 **vertex(꼭짓점)**를 이용할 수 있다. 최적값이 유한하고 feasible set에 꼭짓점이 존재하면 최적해 중 하나를 꼭짓점에서 찾을 수 있다.

Simplex algorithm은 한 꼭짓점에서 시작해, 목적함수 값을 개선할 수 있는 인접 꼭짓점으로 이동한다. 강의에서는 원점에서 시작할 수 있는 예를 먼저 살펴본다. 일반적인 LP에서는 원점이 feasible하지 않을 수도 있다.

## 2. Algebraic Interpretation of a Vertex (p.4)

부등식 제약에 **slack variable**을 추가한다.

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3Bz_i%3Db_i-A_i%20x%5Cquad%28i%3D1%2C%5Cldots%2Cm%29)

Slack variable을 도입하면 제약을 등식으로 쓰고, 원래 변수와 slack variable에 비음수 조건을 부여할 수 있다. 꼭짓점은 독립적인 제약이 충분히 활성화되는 점에 대응한다. 강의의 tableau 표현에서는 non-basic variable을 0으로 두고 basic variable의 값을 계산한다.

## 3. Initial Simplex Tableau (pp.5–6)

Short simplex tableau의 열은 non-basic variable, 상수항, basic variable 이름으로 구성된다. 처음에는 원래 변수 `x1, ..., xn`이 non-basic이고 slack variable `z1, ..., zm`이 basic이다.

| Basic variable | −x₁ | −x₂ | … | −xₙ | 1 (RHS) |
|---|---:|---:|---:|---:|---:|
| z₁ | a₁₁ | a₁₂ | … | a₁ₙ | b₁ |
| ⋮ | ⋮ | ⋮ | ⋮ | ⋮ | ⋮ |
| zₘ | aₘ₁ | aₘ₂ | … | aₘₙ | bₘ |
| v | −c₁ | −c₂ | … | −cₙ | 0 |

- **Non-basic variables:** 상단에 표시되며 현재 값이 0이다.
- **Basic variables:** 우측에 표시되며 RHS 열의 값을 갖는다.
- **v:** 목적함수와 현재 목적함수 값이다.

한 행은 다음과 같은 등식을 나타낸다.

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3Bz_i%3Db_i-%5Csum_%7Bj%3D1%7D%5E%7Bn%7Da_%7Bij%7Dx_j)

## 4. Pivot Operation (p.7)

한 꼭짓점에서 다음 꼭짓점으로 이동하려면 basic variable과 non-basic variable을 교환하는 **pivot**을 수행한다. Pivot 위치를 `(i₀, j₀)`로 놓고, 변경 전 원소를 `dᵢⱼ`, 변경 후 원소를 `d̂ᵢⱼ`라고 하자.

Pivot element는 역수로 바뀐다.

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cwidehat%20d_%7Bi_0j_0%7D%3D%5Cfrac%7B1%7D%7Bd_%7Bi_0j_0%7D%7D)

Pivot 행의 다른 원소는 pivot element로 나눈다.

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cwidehat%20d_%7Bi_0j%7D%3D%5Cfrac%7Bd_%7Bi_0j%7D%7D%7Bd_%7Bi_0j_0%7D%7D%5Cqquad%28j%5Cne%20j_0%29)

Pivot 열의 다른 원소는 부호를 바꾸고 pivot element로 나눈다.

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cwidehat%20d_%7Bij_0%7D%3D-%5Cfrac%7Bd_%7Bij_0%7D%7D%7Bd_%7Bi_0j_0%7D%7D%5Cqquad%28i%5Cne%20i_0%29)

나머지 원소는 다음과 같이 갱신한다.

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cwidehat%20d_%7Bij%7D%3Dd_%7Bij%7D-%5Cfrac%7Bd_%7Bi_0j%7Dd_%7Bij_0%7D%7D%7Bd_%7Bi_0j_0%7D%7D%5Cqquad%28i%5Cne%20i_0%2C%5C%20j%5Cne%20j_0%29)

## 5. Simplex Algorithm — Phase II (p.8)

Phase II는 **feasible basic solution**에서 출발해 목적함수를 최대화하는 과정이다.

1. 목적함수 행에 음수 계수가 남아 있는지 확인한다.
2. 가장 작은 음수 계수에 해당하는 열 `j₀`를 선택한다(entering variable).
3. 해당 열의 계수가 양수인 제약행에 대해서만 RHS / 계수의 비율을 계산한다.
4. 그 비율이 가장 작은 행 `i₀`를 선택한다(leaving variable).
5. `(i₀, j₀)`에서 pivot을 수행하고 반복한다.

Ratio test는 다음과 같다.

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3Bi_0%5Cin%5Coperatorname%2A%7Barg%5C%2Cmin%7D_%7Bi%3A%5C%2Cd_%7Bij_0%7D%3E0%7D%5Cfrac%7Bd_%7Bi%2Cn%2B1%7D%7D%7Bd_%7Bij_0%7D%7D)

들어올 변수를 제한하는 행이 없다면 목적함수는 위로 **unbounded**이다. 목적함수 행에 음수가 없으면 해당 단계에서 종료한다.

## 6. Why Phase I Is Needed (p.9)

원점이 항상 feasible한 것은 아니다. 강의에 나온 예는 다음과 같다.

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cmax%5C%3B%28x_1%2Bx_2%29)

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cbegin%7Baligned%7Dx_1%2B2x_2%26%5Cle6%2C%5C%5C%20x_1%2Bx_2%26%5Cge2%2C%5C%5C%20x_1%2Cx_2%26%5Cge0.%5Cend%7Baligned%7D)

`x₁ = x₂ = 0`은 두 번째 제약을 만족하지 않는다. 따라서 Phase II 이전에 feasible vertex를 찾는 **Phase I**이 필요하다.

## 7. Simplex Algorithm — Phase I (pp.10–11)

강의에서 제시한 Phase I 절차는 tableau의 RHS가 음수인 행을 대상으로 pivot을 반복한다.

1. RHS가 음수인 행이 있는지 확인한다.
2. 음수 RHS를 가진 행 중 강의의 규칙에 따라 행 `î`를 선택한다.
3. 해당 행에서 음수인 계수를 가진 열 `j₀`를 선택한다.
4. 강의의 ratio test로 pivot 행 `i₀`를 정한다.
5. Pivot 후 다시 RHS를 검사한다.

선택한 행에 필요한 음수 계수가 없다면 feasible set이 **empty**라고 판정한다. RHS가 모두 비음수가 되면 Phase II로 넘어간다.

강의자료 p.11은 **Phase I → Phase II** 순서로 두 절차를 이어서 제시한다.

## 8. Special Cases (pp.12–13)

### 8.1 Empty Feasible Set

제약을 동시에 만족하는 점이 없으면 Phase I에서 이를 판정한다.

![p.12 — Empty feasible set](/assets/images/ot_lp2_p12_empty_feasible.png)

*원본 강의 슬라이드: p.12 — Empty feasible set.*

### 8.2 Unbounded Objective Function

Feasible 방향으로 계속 이동하면서 목적함수 값이 무한히 증가할 수 있다면 Phase II에서 unbounded라고 판정한다.

![p.12 — Unbounded objective function](/assets/images/ot_lp2_p12_unbounded.png)

*원본 강의 슬라이드: p.12 — Unbounded objective function.*

### 8.3 Degeneracy and Cycling

두 변수에 비해 지나치게 많은 제약 직선이 하나의 꼭짓점에서 만나는 등 **degenerate vertex**가 있으면 pivot을 반복하더라도 점이 변하지 않을 수 있다. 특정 pivot 선택에서는 이전 tableau로 되돌아가는 **cycling**이 발생할 수 있다.

강의자료에서는 실제 구현 시 tableau 값에 perturbation을 주거나 변수의 label을 이용하는 방식을 언급한다.

![p.13 — Degenerate vertex and cycling](/assets/images/ot_lp2_p13_degeneracy.png)

*원본 강의 슬라이드: p.13 — Degenerate vertex and cycling.*

## 9. Example 1 — Initial Tableau (p.14)

이전 강의의 LP 예제를 다시 사용한다.

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cmax%5C%3B%2810x_1%2B9x_2%29)

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cbegin%7Baligned%7D6x_1%2B5x_2%26%5Cle60%2C%5C%5Cx_1%2B2x_2%26%5Cle15%2C%5C%5Cx_1%26%5Cle8%2C%5C%5Cx_1%2Cx_2%26%5Cge0.%5Cend%7Baligned%7D)

초기 RHS가 모두 비음수이므로 **Phase I 없이 Phase II부터** 시작할 수 있다.

| Basic variable | −x₁ | −x₂ | 1 (RHS) |
|---|---:|---:|---:|
| z₁ | 6 | 5 | 60 |
| z₂ | 1 | 2 | 15 |
| z₃ | 1 | 0 | 8 |
| v | −10 | −9 | 0 |

## 10. Example 1 — Pivot Sequence (pp.15–16)

처음에는 목적함수 행에서 `−10`을 선택한다. Ratio test는 `60/6 = 10`, `15/1 = 15`, `8/1 = 8`이므로 세 번째 행을 pivot 행으로 선택한다. 이후 pivot을 반복하면 강의자료의 최종 tableau를 얻는다.

각 tableau의 상단에는 non-basic variable이, 왼쪽에는 확인을 위한 basic variable 이름이 적혀 있다. Pivot을 수행할 때마다 두 변수의 역할이 바뀐다.

**Step 1. 초기 tableau**

| Basic variable | −x₁ | −x₂ | RHS |
|---|---:|---:|---:|
| z₁ | 6 | 5 | 60 |
| z₂ | 1 | 2 | 15 |
| z₃ | **1** | 0 | 8 |
| v | −10 | −9 | 0 |

**Step 2. x₁이 들어오고 z₃이 나간 후**

| Basic variable | −z₃ | −x₂ | RHS |
|---|---:|---:|---:|
| z₁ | −6 | **5** | 12 |
| z₂ | −1 | 2 | 7 |
| x₁ | 1 | 0 | 8 |
| v | 10 | −9 | 80 |

**Step 3. x₂가 들어오고 z₁이 나간 후**

| Basic variable | −z₃ | −z₁ | RHS |
|---|---:|---:|---:|
| x₂ | −6/5 | 1/5 | 12/5 |
| z₂ | **7/5** | −2/5 | 11/5 |
| x₁ | 1 | 0 | 8 |
| v | −4/5 | 9/5 | 508/5 |

**Step 4. z₃이 들어오고 z₂가 나간 후 (final tableau)**

| Basic variable | −z₂ | −z₁ | RHS |
|---|---:|---:|---:|
| x₂ | 6/7 | −1/7 | 30/7 |
| z₃ | 5/7 | −2/7 | 11/7 |
| x₁ | −5/7 | 2/7 | 45/7 |
| v | 4/7 | 11/7 | 720/7 |

마지막 목적함수 행에는 음수 계수가 없으므로 Phase II가 종료된다.

최종해와 목적함수 값은 다음과 같다.

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3Bx_1%3D%5Cfrac%7B45%7D%7B7%7D%2C%5Cqquad%20x_2%3D%5Cfrac%7B30%7D%7B7%7D%2C%5Cqquad%20v_%7B%5Cmax%7D%3D%5Cfrac%7B720%7D%7B7%7D)

강의자료는 tableau의 변화를 다음 꼭짓점 이동 경로와 연결한다.

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%280%2C0%29%5Clongrightarrow%288%2C0%29%5Clongrightarrow%5Cleft%288%2C%5Cfrac%7B12%7D%7B5%7D%5Cright%29%5Clongrightarrow%5Cleft%28%5Cfrac%7B45%7D%7B7%7D%2C%5Cfrac%7B30%7D%7B7%7D%5Cright%29)

![p.16 — Vertex-to-vertex path](/assets/images/ot_lp2_p16_vertex_path.png)

*원본 강의 슬라이드: p.16 — Vertex-to-vertex path.*

## 11. Example 2 — A Problem Requiring Phase I (p.17)

두 번째 예제는 다음과 같다.

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cmax%5C%3B%28x_1%2B2x_2%29)

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3B%5Cbegin%7Baligned%7Dx_1-x_2%26%5Cge1%2C%5C%5Cx_2%26%5Cge2%2C%5C%5Cx_1-x_2%26%5Cge-1%2C%5C%5Cx_1%26%5Cle4%2C%5C%5Cx_1%2Cx_2%26%5Cge0.%5Cend%7Baligned%7D)

이 문제는 초기 tableau의 RHS에 음수 값이 있으므로 Phase I이 필요하다.

| Basic variable | −x₁ | −x₂ | 1 (RHS) |
|---|---:|---:|---:|
| z₁ | −1 | 1 | −1 |
| z₂ | 0 | −1 | −2 |
| z₃ | −1 | 1 | 1 |
| z₄ | 1 | 0 | 4 |
| v | −1 | −2 | 0 |

## 12. Example 2 — Phase I (p.18)

RHS가 음수인 상태에서 pivot을 반복해 feasible tableau를 만든다. 아래는 강의 슬라이드의 중간 tableau와 pivot 위치다.

강의에서는 아래 네 개의 tableau를 순서대로 사용한다. 각 표의 열 제목은 non-basic variable이며, RHS가 모두 비음수가 되면 Phase I을 종료한다.

**Phase I — Initial tableau**

| Basic variable | −x₁ | −x₂ | RHS |
|---|---:|---:|---:|
| z₁ | −1 | 1 | −1 |
| z₂ | 0 | −1 | −2 |
| z₃ | −1 | **1** | 1 |
| z₄ | 1 | 0 | 4 |
| v | −1 | −2 | 0 |

**Pivot 1 — x₂ enters, z₃ leaves**

| Basic variable | −x₁ | −z₃ | RHS |
|---|---:|---:|---:|
| z₁ | 0 | −1 | −2 |
| z₂ | **−1** | 1 | −1 |
| x₂ | −1 | 1 | 1 |
| z₄ | 1 | 0 | 4 |
| v | −3 | 2 | 2 |

**Pivot 2 — x₁ enters, z₂ leaves**

| Basic variable | −z₂ | −z₃ | RHS |
|---|---:|---:|---:|
| z₁ | 0 | **−1** | −2 |
| x₁ | −1 | −1 | 1 |
| x₂ | −1 | 0 | 2 |
| z₄ | 1 | 1 | 3 |
| v | −3 | −1 | 5 |

**Pivot 3 — z₃ enters, z₁ leaves**

| Basic variable | −z₂ | −z₁ | RHS |
|---|---:|---:|---:|
| z₃ | 0 | −1 | 2 |
| x₁ | −1 | −1 | 3 |
| x₂ | −1 | 0 | 2 |
| z₄ | 1 | 1 | 1 |
| v | −3 | −1 | 7 |

최종 RHS가 모두 0 이상이므로 feasible basic solution을 얻었다.

## 13. Example 2 — Phase II (p.19)

Phase I이 끝나 모든 RHS가 비음수가 되면 Phase II를 수행한다. 마지막 tableau의 basic variable 값으로부터 다음 해를 읽을 수 있다.

![수식](https://latex.codecogs.com/svg.image?%5Ccolor%7Bwhite%7D%5Cdisplaystyle%5C%3Bx_1%3D4%2C%5Cqquad%20x_2%3D3%2C%5Cqquad%20v_%7B%5Cmax%7D%3D10)

Phase I 종료 시 tableau에서 시작해, 목적함수 행에 남은 음수 계수를 기준으로 pivot을 한 번 더 수행한다.

**Phase II — Initial tableau**

| Basic variable | −z₂ | −z₁ | RHS |
|---|---:|---:|---:|
| z₃ | 0 | −1 | 2 |
| x₁ | −1 | −1 | 3 |
| x₂ | −1 | 0 | 2 |
| z₄ | **1** | 1 | 1 |
| v | −3 | −1 | 7 |

**Phase II — Final tableau**

| Basic variable | −z₄ | −z₁ | RHS |
|---|---:|---:|---:|
| z₃ | 0 | −1 | 2 |
| x₁ | 1 | 0 | 4 |
| x₂ | 1 | 1 | 3 |
| z₂ | 1 | 1 | 1 |
| v | 3 | 2 | 10 |

따라서 `x₁ = 4`, `x₂ = 3`이고 최댓값은 `10`이다.

---

**Reference**  
Piotr Więcek, *Optimization Theory — Linear Programming Part II: The Simplex Algorithm*, Wrocław University of Science and Technology, 6 November 2026, pp. 1–19.
