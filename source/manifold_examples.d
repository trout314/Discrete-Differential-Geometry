module manifold_examples;

import manifold : Manifold;
import unit_threaded : Name, shouldBeSameSetAs;
import std.range : iota;
import utility : subsetsOfSize;

auto trigonalBipyramid()
{
    return Manifold!2([[0,1,2],[0,1,3],[0,2,3],[1,2,4],[1,3,4],[2,3,4]]);
}

auto octahedron()
{
    return Manifold!2([[0,1,2], [0,2,3], [0,3,4], [0,1,4], [1,2,5],
        [2,3,5], [3,4,5], [1,4,5]]);
}

/******************************************************************************
Returns the "standard" triangulation of a sphere of the given dimension. This
is just the boundary of the simplex of one higher dimension.
*/
Manifold!(dim, int) standardSphere(int dim)()
{
    static assert(dim >= 1);
    return Manifold!(dim, int)(iota(dim + 2).subsetsOfSize(dim + 1));
}
///
@Name("standardSphere") unittest
{
    standardSphere!1.facets.shouldBeSameSetAs(
        [[0,1], [1,2], [0,2]]); 
    standardSphere!2.facets.shouldBeSameSetAs(
        [[0,1,2], [0,1,3], [0,2,3], [1,2,3]]);
    standardSphere!3.facets.shouldBeSameSetAs(
        [[0,1,2,3], [0,1,2,4], [0,1,3,4], [0,2,3,4],[1,2,3,4]]);
}

