# Discrete-Differential-Geometry

Work in progress. Symbolic math package for discrete differential geometry computations. Some goals:

- [x] abstract simplicial complex type
- [x] combinatorial *n*-manifold type
- [x] algorithms for simplicial complexes
    - [x] list connected components
    - [x] test for orientability
    - [x] test if complex is a combinatorial *n*-manifold (*n* = 1,2,3)
    - [x] return join of two complexes
    - [x] compute Euler characteristic
    - [ ] edge diameter
    - [ ] randomly generated complexes
- [x] Metropolis-Hastings sampling of manifolds via bistellar moves
    - objective can target any combination of:
        - [x] number of *n*-simplices
        - [x] average degree of the codimension-2 simplices
        - [x] standard deviation in degree for codimension-2 simplices
        - [ ] more general user-defined local objectives
    - [ ] hinge moves
    - [ ] automatic temperature control
- [ ] vertex visibility determination in combinatorial manifolds
- [ ] pairwise discrete distance determination