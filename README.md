# Nano Course


## Setup

Checkout
```
git clone --recurse-submodules https://github.com/jonasschultemattler/nano.git
```

Compile
```
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release .. -D CMAKE_CXX_COMPILER=g++-14
make
```

Test
```
./source/hello
```

## Counting distinct elements of a set

**online** (stream)


### Naive solution

Hashmap or bitvector

TODO implement...

Space linear w.r.t. number distinct elements


-> approximate/probabilistic counting


### Flajolet-Martin’s algorithm

Recall:\
In a set of (uniformly distributed) hash values of cardinality $2^l$ , we expect to see one hash $h(x)$ prefixed by $l$ zeroes.

(the probability of observing a binary encoded hash beginning with $k$ zeroes followed by a one is $1/2^(k+1)$)

Algorithm:\
map each element $x$ to a $q$-bits hash $h(x)$,\
remember the maximum number $l = lb(h(x))$ of leading 0-bits seen in any $h(x)$,\
finally, return the estimate $2^l$ 


TODO implement...

TODO: use different hash functions


### HyperLogLog Sketching

Refinement:



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


