import std.algorithm : all, each, find, map, sort;
import std.conv : to;
import std.exception : enforce;
import std.range : array, front, isInputRange, ElementType, walkLength;
import std.traits : isArray;
import simplex : Simplex, simplex;
import simplicial_complex : SimplicialComplex, isCircle;
import utility : staticIota, subsets;

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

        enforce(initialFacets.all!(f => f.walkLength == dimension + 1),
            "expected facets of dimension " ~ dimension.to!string
            ~ " but got the facet " 
            ~ initialFacets.find!(f => f.walkLength != dimension + 1)
                .front.to!string);

        initialFacets.each!(f => simpComp_.insertFacet(f));

        static if(dimension >= 1)
        {
            auto badRidge = this.simplices!(dimension - 1).find!(
                s => this.degree(s) != 2);

            enforce(badRidge.empty, "manifold constructor expected ridges of "
            ~ "degree 2, but found a ridge " ~ badRidge.front.to!string 
            ~ " with degree " ~ this.degree(badRidge.front).to!string);
        }

        static if(dimension >= 2)
        {
            auto badHinge = this.simplices!(dimension - 2).find!(
                s => !(this.link(s).to!(SimplicialComplex!Vertex).isCircle));

            enforce(badHinge.empty, "manifold constructor expected hinges whose "
                ~ "links are circles but found a hinge " 
                ~ badHinge.front.to!string ~ " with link " 
                ~ this.link(badHinge.front).to!string);
        }

        static if(dimension >= 3)
        {
            // auto badCodim3 = this.simplices!(dimension - 3).find!(
            //     s => !(this.link(s).to!(SimplicialComplex!Vertex).eulerCharacteristic == 2));
        }

        enforce(this.isConnected, "manifold constructor expected a connected "
            ~ "simplicial complex but got one with "
            ~ this.connectedComponents.walkLength.to!string 
            ~ " connected components");
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

    auto pachnerMoves()
    {
        int[Vertex[]] degreeMap;

        // TO DO: Why do I need idup here?
        this.facets.each!(f => f.subsets.each!(s => ++degreeMap[s.idup]));

        Vertex[][] result;
        foreach(simp, deg; degreeMap)
        {

            /* TO DO: IMPORTANT! Check that for the replacement ball of the 
            form star(s) the simplex s is not already present. */

            if(deg == dimension + 2 - simp.walkLength)
            {
                result ~= simp.dup;
            }           
        }
        result.sort!mySort;

        return result;
    }


private:
    SimplicialComplex!Vertex simpComp_;
}

bool mySort(A, B)(A a, B b)
{
    if (a.walkLength < b.walkLength)
    {
        return true;
    }
    else if (a.walkLength == b.walkLength)
    {
        return a < b;
    }

    return false;
}

unittest
{
    import manifold_test : test;
    assert(test!SmallManifold);
}