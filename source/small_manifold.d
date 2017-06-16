import std.algorithm : each;
import std.conv : to;
import std.range : isInputRange, ElementType;
import std.traits : isArray;
import simplex : Simplex, simplex;
import simplicial_complex : SimplicialComplex;

struct SmallManifold(int dim, Vertex = int)
{
    static assert(dim >= 1, "dimension must be at least one, but got dimension " 
        ~ dim.to!string);

    alias Facet = Simplex!(dim, Vertex);    
    this(F)(F facets_)
    {
        static assert(isInputRange!F || isArray!F,
            "manifold must be constructed from a range or array");
        
        static assert(is(ElementType!F == Facet), "manifold must be constructed"
            ~ "from a range or array of elements with type " ~ Facet.stringof);

        facets_.each!(f => facets.insertFacet(f));

        static if(dim == 1)
        {
            // TO DO: Check manifold-ness
        }

        static if(dim == 2)
        {
            // TO DO: Check manifold-ness
        }

        static if(dim == 2)
        {
            // TO DO: Check manifold-ness
        }


    }
private:
    SimplicialComplex!Vertex facets;
}

unittest
{
    import std.stdio : writeln;

    alias s = simplex;

    auto sm = SmallManifold!1([s(1,2), s(2,3), s(1,3)]);

}