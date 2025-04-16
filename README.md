# Nano Course


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
 - the cardinality of $M$ can be estimated by the maximum number of leading zeros in the binary representation of each number in $M$
 - if max leading zeros is $l$, one exepcts $2^l$ distinct elements
(the probability of observing a binary encoded number beginning with $k$ zeroes followed by a one is $1/2^{(k+1)}$ )

#### Algorithm:

- map each element $x$ to hash $h(x)$
- remember the maximum number $l$ of leading 0-bits seen in any $h(x)$
- estimate cardinality by $2^l$ 

hash $h(x) \rightarrow [0,L]$ requires $\log(L) \approx \log(n)$ space for $n$ distinct elements.

TODO implement...

TODO: Test different hash functions


### HyperLogLog

#### Observation:
Large Variance

#### Refinement:
- split $\mathcal{M}$ into subsets
- estimate cardinalities of subsets
- return mean

The normalized version of the harmonic mean is the estimate $E:=\frac{\alpha_m m^2}{\sum_{j=1}^m 2^{-M(j)}}$.

smaller variance, $\log \log n$ space

TODO implement...


### Evaluation

Plot time, space, gap/solution quality



## Set Similarity

Jaccard similarity $J(A,B)$ of two sets $A$ and $B$: 


### Exact

TODO implement...

### FracMinHashing

TODO implement...

### HyperLogLog

TODO implement...



### Evaluation

Plot time, space, gap/solution quality


