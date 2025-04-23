# Nanocourse - Genomic analysis using sketching techniques


## Setup

Requirements:

 - gcc >= 12 or clang >=17
 - cmake >= 3.20
 - git

Checkout
```
git clone --recurse-submodules https://github.com/jonasschultemattler/nano.git
```

See https://docs.seqan.de/seqan3/main_user/setup.html for compiler setup. Compile with

```
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release .. -D CMAKE_CXX_COMPILER=g++-14
make
```

Test
```
./source/count
```

## Sketching

A *data sketch* of data $X$ is the output of a (randomized) function $f$ s.t.:
 - $|f(X)| \subseteq o(|X|)$
 - $f$ perserves some properties of $X$ e.g. approximation of the number of distinct elements
 - $f$ preserves certain similarity measures e.g. number of shared elements
 - $f(X)$ allows to be updated efficiently


## Counting distinct elements of a set

In an **online** setting (stream)


### Naive solution

Using a Hashmap or a Bitvector

TODO implement naive_couting() in count.cpp

<!-- TODO plot space consumption/number of k-mers in plot.py -->


#### Observation:

- Space consumption linear w.r.t. number distinct elements
- impractical for big data

-> approximate/probabilistic counting


### Flajolet-Martin’s algorithm

#### Recall:

Let $\mathcal{M}$ be a multiset of uniformly distributed random numbers.
 - The cardinality of $\mathcal{M}$ can be estimated by the maximum number of leading zeros in the binary representation of each number in $\mathcal{M}$.
 - If max leading zeros is $l$, one exepcts $2^l$ distinct elements
(the probability of observing a binary encoded number beginning with $k$ zeroes followed by a one is $1/2^{(k+1)}$ ).

#### Algorithm:

- Map each element $x$ to hash $h(x)$,
- remember the maximum number $l$ of leading 0-bits seen in any $h(x)$,
- estimate cardinality by $2^l$.

TODO:
- Implement flajolet_martin() in count.cpp.
- Test flajolet_martin() for different hash functions.
- Compare run time, space consumption and accuracy to exact solution.


#### Observations:

 - Hash $h(x) \rightarrow [0,L]$ requires $\log(L) \approx \log(n)$ space for $n$ distinct elements.
 - Large Variance


### HyperLogLog

#### Refinement:
- split $\mathcal{M}$ into $m$ subsets,
- estimate cardinalities of subsets,
- return mean.

The normalized version of the harmonic mean is the estimate
```math
E:=\frac{\alpha_m m^2}{\sum_{j=1}^m 2^{-M(j)}}.
```
for $m$ subsets $M(i)$ and normalization constant $\alpha_m \approx 0.7$.

TODO
- Implement hyperloglog() in count.cpp.
- Test hyperloglog() for different hash functions.
- Compare run time, space consumption and accuracy to Flajolet-Martin.


#### Observations:

- Smaller variance
- Less space concumption: $O(\log \log n)$



## Set Similarity

(Dis-)similarity of two sets $A$ and $B$ can be measured with Jaccard similarity
```math
J(A,B) := \frac{|A \cap B|}{|A \cup B|}
```


### Naive Algorithm

Hashmap or bitvector

TODO implement...

Space linear w.r.t. total number of distinct elements

Comparing $n$ sets requires $O(n^2)$ pairwise comarisons.


### MinHashing

Let MinHash $h_{\min}(A) = \min \lbrace h(x) \mid x \in A \rbrace$ and
```math
J_h(A,B) := \begin{cases}1, & \text{if } h_{\min}(A) = h_{\min}(B)\\ 0 & \text{otw.}\end{cases},
```
Then $E[J_h(A,B)] = J(A,B)$.

Algorithm:

- sample $h_{\min}$ from $k$ different hash functions
- let $l$ be the number of hash functions with $h_{\min}(A) = h_{\min}(B)$
- estimate $J(A,B)$ by $l/k$

TODO implement...

For any $\epsilon > 0$ there is a $k \in O(1/\epsilon^2)$ s.t. the expected error is at most $\epsilon$.

How many hash functions do you need to have an expected error at most $.05$?


### FracMinHashing

Let $h: \Omega \Rightarrow [0,H]$ for $H \in \mathbb{N}$ and $s \in [0,1]$ be a scaling factor. FracMinHash
```math
Ffrac_{s}(A) = \min \{ h(x) \mid x \in A h(x) \leq H s \}.
```

TODO implement...



### Evaluation

Plot time, space, gap/solution quality



## Take Aways

- Certain tasks for massive data as in molecular biology require data sketching.
- With a bit of randomness, measures for distinct elements in a set become tractable for big data.
- Especially, relying proximity measures that require pairwise comparisons.





