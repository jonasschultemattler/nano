# Nano Course - Data Sketching in Biology


## Setup

Requirements:

 - gcc >= 12 or clang >=17
 - cmake >= 3.20
 - git

Checkout
```
git clone --recurse-submodules https://github.com/jonasschultemattler/nano.git
```

Compile

See https://docs.seqan.de/seqan3/main_user/setup.html for compiler setup. Then compile with

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
 - $f(X)$ perserves some properties of $X$ e.g. approximation of the number of distinct elements
 <!-- - it can be updated efficiently -->
 - preserves some similarity measure


## Counting distinct elements of a set

**online** setting (stream)


### Naive solution

Hashmap or bitvector

TODO implement...

Space linear w.r.t. number distinct elements

-> approximate/probabilistic counting


### Flajolet-Martin’s algorithm

#### Recall:

Let $\mathcal{M}$ be a multiset of uniformly distributed random numbers.
 - the cardinality of $\mathcal{M}$ can be estimated by the maximum number of leading zeros in the binary representation of each number in $\mathcal{M}$
 - if max leading zeros is $l$, one exepcts $2^l$ distinct elements
(the probability of observing a binary encoded number beginning with $k$ zeroes followed by a one is $1/2^{(k+1)}$ )

#### Algorithm:

- map each element $x$ to hash $h(x)$
- remember the maximum number $l$ of leading 0-bits seen in any $h(x)$
- estimate cardinality by $2^l$ 

hash $h(x) \rightarrow [0,L]$ requires $\log(L) \approx \log(n)$ space for $n$ distinct elements.

TODO implement...


### HyperLogLog

#### Observation:
Large Variance

#### Refinement:
- split $\mathcal{M}$ into $m$ subsets
- estimate cardinalities of subsets
- return mean

The normalized version of the harmonic mean is the estimate
```math
E:=\frac{\alpha_m m^2}{\sum_{j=1}^m 2^{-M(j)}}.
```

It has smaller variance and requires $\log \log n$ space.


TODO implement...


### Evaluation

TODO: Plot time, space, gap/solution quality

TODO: Test different hash functions


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

Let MinHash $h_{\min}(A) = \min \{ h(x) \mid x \in A \}$ and
```math
J_h(A,B) := \begin{cases}1, & \text{if } h_{\min}(A) = h_{\min}(B)\\ 0 & \text{otw.}\end{cases},
```
then $E[J_h(A,B)] = J(A,B)$.

Algorithm:

- sample $h_{\min}$ for $k$ different hash functions
- let $l$ be the number of hash functions s.t. $h_{\min}(A) = h_{\min}(B)$
- estimate $J(A,B)$ by $l/k$

TODO implement...

For any $\epsilon > 0$ there is a $k \in O(1/\epsilon^2)$ s.t. the expected error is at most $\epsilon$.

How many hash functions do you need to have an expected error at most $.05$?


### FracMinHashing

<!-- Let FracMinHash $h_{frac}(A) = \min \{ s h(x) \mid x \in A \}$ for $s \in [0,1]$. -->

TODO implement...



### Evaluation

Plot time, space, gap/solution quality



## Take Aways

- Certain tasks on massive data as in molecular biology require data sketching.
- With a bit of randomness measures of distinct elements in a set become tractable for big data.
- Especially, relying proximity measures that require pairwise comparison.





