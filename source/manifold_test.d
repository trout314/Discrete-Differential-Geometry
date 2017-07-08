import std.conv : to;
import std.algorithm : equal, map;
import simplex : Simplex, simplex;
import simplicial_complex;
import std.stdio : writeln;
import utility : throwsWithMsg;

auto test(alias Manifold)()
{
    // The octahedral surface with eight triangular facets
    auto octahedron = Manifold!2([[0, 1, 2], [0, 2, 3], [0, 3, 4], [0, 1, 4],
        [1, 2, 5], [2, 3, 5], [3, 4, 5], [1, 4, 5]]);

    static assert(octahedron.dimension == 2);
    static assert(is(octahedron.Vertex == int));
    static assert(is(octahedron.Facet == Simplex!2));
    static assert(is(octahedron.Facet == Simplex!(2, int)));

    assert(octahedron.fVector == [6,12,8]);
    assert(octahedron.eulerCharacteristic == 2);

    alias s = simplex;

    assert(octahedron.star(s(1,2)).equal([s(0,1,2), s(1,2,5)]));    
    static assert(!__traits(compiles, Manifold!2([["a", "bubba", "gump"]])));

    auto m1 = Manifold!(1, string)([["a", "b"], ["b", "c"], ["a", "c"]]);
    static assert(m1.dimension == 1);
    static assert(is(m1.Vertex == string));
    static assert(is(m1.Facet == Simplex!(1, string)));

    throwsWithMsg(Manifold!2([[1,2,3]]), "manifold constructor expects ridges "
        ~ "of degree 2, but found a ridge [1,2] with degree 1");

    // auto f = [[1,2,3], [2,3,4]].map!(v => Simplex!2(v));
    // auto y = Manifold!2(f);

    return true;
}