import std.algorithm : all, each, map;
import std.conv : to;
import std.exception : enforce;
import std.range : isInputRange, ElementType, walkLength;
import std.traits : isArray;
import simplex : Simplex, simplex;
import simplicial_complex : SimplicialComplex;

struct SmallManifold(int dimension_, Vertex_ = int)
{
    static immutable dimension = dimension_;
    alias Vertex = Vertex_;
    alias Facet = Simplex!(dimension, Vertex);    

    static assert(dimension >= 1, "dimension must be at least one, but got "
        ~ "dimension " ~ dimension.to!string);

    this(F)(F initialFacets)
    {
        static assert((isInputRange!F || isArray!F)
            && is(ElementType!F == Facet) || is(ElementType!F == Vertex[]),
            "a " ~ SmallManifold!(dimension, Vertex).stringof ~ " must be "
            ~ "constructed from a range or array of elements with type " 
            ~ Facet.stringof ~ " or " ~ Vertex[].stringof ~ " but got element "
            ~ "type " ~ ElementType!F.stringof);

        initialFacets.each!(f => simpComp_.insertFacet(f));

        // TO DO: Check links of codimension-1 simplices

        static if(dimension >= 2)
        {
            // TO DO: Check links of codimension-2 simplices
        }

        static if(dimension >= 3)
        {
            // TO DO: Check links of codimension-3 simplices
        }
    }

    auto star(int dim)(const Simplex!(dim, Vertex) s) const
    {
        return simpComp_.star(s).map!(verts => Facet(verts));
    }

    alias asSimplicialComplex this;
    ref const(SimplicialComplex!Vertex) asSimplicialComplex()
    {
        return simpComp_; 
    }


private:
    SimplicialComplex!Vertex simpComp_;
}

unittest
{
    import manifold_test : test;

    assert(test!SmallManifold);

    import std.stdio : writeln;

    alias s = simplex;

    auto sm = SmallManifold!1([s(1,2), s(2,3), s(1,3)]);
    sm.writeln;

    auto sm2 = SmallManifold!1([[1,2], [2,3], [1,3]]);
    assert(sm == sm2);

}