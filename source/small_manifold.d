import std.algorithm : all, each, find, map;
import std.conv : to;
import std.exception : enforce;
import std.range : array, isInputRange, ElementType, walkLength;
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

        static if(dimension >= 1)
        {
            auto badRidge = simpComp_.simplices!(dimension - 1).find!(
                s => degree(s) != 2);
            enforce(badRidge.empty, "manifold constructor expects ridges of "
            ~ "degree 2, but found a ridge " ~ badRidge.front.to!string 
            ~ " with degree " ~ degree(badRidge.front).to!string);
        }

        static if(dimension >= 2)
        {
            // this.simplices!(dimension - 2).map(
            //     s => link(s).to!SimplicialComplex)
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

    auto degree(int dim)(Simplex!(dim, Vertex) s) const
    {
        return star(s).walkLength;
    }


private:
    SimplicialComplex!Vertex simpComp_;
}

unittest
{
    import manifold_test : test;
    assert(test!SmallManifold);
}