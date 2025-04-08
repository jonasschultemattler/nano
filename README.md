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
In a set of hash values of cardinality $2^l$ , we expect to see one hash $h(x_i)$ prefixed by $l$ zeroes.

Algorithm:\
map each element $x_i$ to a $q$-bits hash $h(x_i)$,\
remember the maximum number $l = lb(h(x_i))$ of leading 0-bits seen in any $h(x_i)$,\
finally, return the estimate $2^l$ 


TODO implement...

TODO: use different hash functions


### HyperLogLog Sketching

Refinement:



### Evaluation

Plot time, plot gap/solution quality



## Set Similarity

Jaccard similarity $J(A,B)$ of two sets $A$ and $B$: 


### Exact

TODO implement...

### FracMinHashing

TODO implement...

### HyperLogLog

TODO implement...



### Evaluation

Plot time, plot gap/solution quality


