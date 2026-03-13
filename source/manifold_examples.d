/// TO DO: Module description
module manifold_examples;

import std.range;
import unit_threaded;
import manifold, utility;

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

