import std.algorithm : all, canFind, each, find, joiner, map, maxElement, sort, uniq;
import std.conv : to;
import std.exception : enforce;
import std.range : array, empty, front, isInputRange, ElementType, popFront, walkLength;
import std.traits : isArray;
import simplex : Simplex, simplex;
import simplicial_complex : SimplicialComplex, isCircle, isPureOfDim, is2Sphere, join;
import utility : staticIota, subsets, subsetsOfSize;

import std.stdio : writeln;

struct SmallManifold(int dimension_, Vertex_ = int)
{
private:
    SimplicialComplex!Vertex simpComp_;
    Vertex[] verticesToReuse;
    static if(is(Vertex_ == int))
    {
        alias getNewVertex =
            sc => sc.simplices(0).map!(v => v[0]).maxElement + 1;
    }

    static if(is(Vertex_ == string))
    {
        alias getNewVertex = sc => "bubba!";
    }
public:
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
            ~ Facet.stringof ~ " or " ~ Vertex[].stringof);

        initialFacets.each!(f => simpComp_.insertFacet(f));

        enforce(this.isPureOfDim(dimension), "expected all facets to have " 
            ~ "dimension " ~ dimension.to!string ~ " but got the facet " 
            ~ facets.find!(f => f.walkLength != dimension + 1).front.to!string);

        enforce(this.isConnected, "manifold constructor expected a connected "
            ~ "simplicial complex but got one with "
            ~ this.connectedComponents.walkLength.to!string 
            ~ " connected components");

        static if(dimension >= 1)
        {{
            auto badRidge = this.simplices!(dimension - 1).find!(
                s => this.degree(s) != 2);

            enforce(badRidge.empty, "manifold constructor expected ridges of "
            ~ "degree 2, but found a ridge " ~ badRidge.front.to!string 
            ~ " with degree " ~ this.degree(badRidge.front).to!string);
        }}

        static if(dimension >= 2)
        {{
            auto badHinge = this.simplices!(dimension - 2).find!(
                s => !(this.link(s).to!(SimplicialComplex!Vertex).isCircle));

            enforce(badHinge.empty, "manifold constructor expected hinges whose"
                ~ " links are circles but found a hinge " 
                ~ badHinge.front.to!string ~ " with link " 
                ~ this.link(badHinge.front).to!string);
        }}

        static if(dimension >= 3)
        {{
            auto badCodim3 = this.simplices!(dimension - 3).find!(
                s => !this.link(s).to!(SimplicialComplex!Vertex).is2Sphere);

            enforce(badCodim3.empty, "manifold constructor expected the links "
                ~ "of all codimension-3 simplices to be 2-spheres but found "
                ~ "simplex " ~ badCodim3.front.to!string ~ " with link "
                ~ this.link(badCodim3.front).to!string);
        }}
    }

    auto star(int dim)(const Simplex!(dim, Vertex) s) const
    {
        return simpComp_.star(s).map!(verts => Facet(verts));
    }

    alias asSimplicialComplex this;
    ref const(SimplicialComplex!Vertex) asSimplicialComplex() const
    {
        return simpComp_; 
    }

    auto degree(int dim)(const Simplex!(dim, Vertex) s) const
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
            if(deg == dimension + 2 - simp.walkLength)
            {
                auto newSimp = this.link(simp).joiner.array.sort().uniq;
                if(newSimp.empty || !this.contains(newSimp))
                {
                    result ~= simp.dup;
                }
            }           
        }
        result.sort!mySort;
        return result;
    }

    void doPachner(Vertex[] s)
    {
        // TO DO: Better error
        enforce(this.pachnerMoves.canFind(s), "bad pachner move");

        auto coSimp = this.link(s).joiner.array.sort().uniq.array;

        if(coSimp.empty)
        {
            if(!verticesToReuse.empty)
            {
                coSimp = [verticesToReuse.front];
                verticesToReuse.popFront;
            }
            else
            {
                coSimp = [getNewVertex(this)];            
            }
        }

        simpComp_.star(s).each!(f => simpComp_.removeFacet(f));

        alias SComp = SimplicialComplex!Vertex;
        auto newPiece = join(
            SComp(s.subsetsOfSize(s.walkLength.to!int - 1)),
            SComp([coSimp])
        );

        assert(newPiece.isPureOfDim(this.dimension));
        newPiece.facets.each!(f => simpComp_.insertFacet(f));
    }
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