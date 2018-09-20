# Discrete-Differential-Geometry

Work in progress. Symbolic math package for discrete differential geometry computations. The current focus is on efficient sampling of the space of combinatorial *n*-manifold triagulations with large numbers of *n*-simplices. 

Some goals:

- [x] abstract simplicial complex type
- [x] combinatorial *n*-manifold type
    - [ ] generate 2-manifolds (surfaces)
    - [ ] generate 3-manifolds from surfaces 
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
        - [ ] more general user-defined local properties
    - [ ] hinge moves
    - [ ] automatic temperature control
- [ ] vertex visibility determination in manifolds
- [ ] pairwise discrete distance determination