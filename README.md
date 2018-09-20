# Discrete-Differential-Geometry

Work in progress. Symbolic math package for discrete differential geometry computations. Some goals:

- [x] abstract simplicial complex type
- [x] combinatorial *n*-manifold type
- [x] algorithms for simplicial complexes
    - [x] list connected components
    - [x] test for orientability
    - [x] test for n-manifold-ness for n=1,2,3
    - [x] return join of two complexes
    - [x] compute Euler characteristic
- [x] Metropolis-Hastings sampling of manifolds via bistellar moves
    - available objective components
        - [x] number of $n$-simplices
    - [ ] hinge moves
    - [ ] automatic temperature control

* Metropolis-Hastings sampling via Pachner moves and user-given objective
    * optimized implementations for special objective functions like total absolute angle-defect over hinges 
* vertex visibility determination in combinatorial manifolds
* pairwise discrete distance determination

