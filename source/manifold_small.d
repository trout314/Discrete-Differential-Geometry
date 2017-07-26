import std.algorithm : all, canFind, each, equal, find, joiner, map, maxElement,
    sort, uniq;
import std.conv : to;
import std.exception : enforce;
import std.range : array, empty, front, isInputRange, ElementType, popFront,
    walkLength;
import std.traits : isArray;
import simplex : Simplex, simplex;
import simplicial_complex : SimplicialComplex, fVector, eulerCharacteristic,
    isCircle, connectedComponents, isConnected, isPureOfDim, is2Sphere, join;
import utility : staticIota, subsets, subsetsOfSize, throwsWithMsg;
import std.stdio : writeln;

import fluent.asserts;
import unit_threaded : Name;

///
unittest
{
    alias Manifold = SmallManifold;

    auto octahedron = Manifold!2([[0,1,2], [0,2,3], [0,3,4], [0,1,4], [1,2,5],
        [2,3,5], [3,4,5], [1,4,5]]);

    static assert(octahedron.dimension == 2);
    static assert(is(octahedron.Vertex == int));
    static assert(is(octahedron.Facet == Simplex!2));
    static assert(is(octahedron.Facet == Simplex!(2, int)));

    assert(octahedron.fVector == [6,12,8]);
    assert(octahedron.eulerCharacteristic == 2);

    auto tetrahedron = Manifold!2([[1,2,3], [1,2,4], [1,3,4], [2,3,4]]);

    alias s = simplex;

    assert(octahedron.star(s(1,2)).equal([[0,1,2], [1,2,5]]));    

    octahedron.pachnerMoves.map!(s => s.idup).should.containOnly(
        [[0,1], [0,2], [0,3], [0,4], [1,2], [1,4],
         [1,5], [2,3], [2,5], [3,4], [3,5], [4,5],      // 1-simplices
         [0, 1, 2], [0, 1, 4], [0, 2, 3], [0, 3, 4],
         [1, 2, 5], [1, 4, 5], [2, 3, 5], [3, 4, 5]]);  // 2-simplices

    tetrahedron.pachnerMoves.map!(s => s.idup).should.containOnly(
        [[1, 2, 3], [1, 2, 4], [1, 3, 4], [2, 3, 4]]);

    // TO DO: Why do I need the map!(s => s.idup) part above?!

    tetrahedron.doPachner([1,2,3], 5);

    // TO DO: FINISH CHECKS!

    octahedron.writeln;
    octahedron.doPachner([1,2]);
    octahedron.writeln;

    // octahedron.doPachner([0,5]);
    octahedron.degree(s(0)).writeln;
}

// More tests
unittest
{
    alias Manifold = SmallManifold;

    static assert(!__traits(compiles, Manifold!2([["a", "bubba", "gump"]])));

    auto m1 = Manifold!(1, string)([["a", "b"], ["b", "c"], ["a", "c"]]);
    static assert(m1.dimension == 1);
    static assert(is(m1.Vertex == string));
    static assert(is(m1.Facet == Simplex!(1, string)));

    throwsWithMsg(Manifold!2([[1,2,3,4]]), "expected all facets to have "
        ~ "dimension 2 but got the facet [1, 2, 3, 4]");

    throwsWithMsg(Manifold!2([[1,2,3]]), "manifold constructor expected ridges "
        ~ "of degree 2, but found a ridge [1,2] with degree 1");

    throwsWithMsg(Manifold!2([[1,2,3], [1,2,4], [1,3,4], [2,3,4], [1,5,6], 
        [1,5,7], [1,6,7], [5,6,7]]), "manifold constructor expected hinges"
        ~ " whose links are circles but found a hinge [1] with link"
        ~ " [[2, 3], [2, 4], [3, 4], [5, 6], [5, 7], [6, 7]]");

    auto octahedron = Manifold!2([[0,1,2], [0,2,3], [0,3,4], [0,1,4], [1,2,5],
        [2,3,5], [3,4,5], [1,4,5]]);
    auto sphere = [[6,7,8], [6,7,9], [6,8,9], [7,8,9]];
    throwsWithMsg(Manifold!2(octahedron.facets ~ sphere), "manifold constructor"
        ~ " expected a connected simplicial complex but got one with"
        ~ " 2 connected components");

    auto sphere3 = Manifold!3(
        [[1,2,3,4], [1,2,3,5], [1,2,4,5], [1,3,4,5], [2,3,4,5]]);

    auto sphere3A = Manifold!3(
        [[1,6,7,8], [1,6,7,9], [1,6,8,9], [1,7,8,9], [6,7,8,9]]);

    throwsWithMsg(Manifold!3(sphere3.facets ~ sphere3A.facets),
        "manifold constructor expected the links of all codimension-3 "
        ~ "simplices to be 2-spheres but found simplex [1] with link "
        ~ "[[2, 3, 4], [2, 3, 5], [2, 4, 5], [3, 4, 5], [6, 7, 8], [6, 7, 9],"
        ~ " [6, 8, 9], [7, 8, 9]]");
}

struct SmallManifold(int dimension_, Vertex_ = int)
{
private:
    SimplicialComplex!Vertex simpComp_;
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

    auto degree(int dim)(const Simplex!(dim, Vertex) s) const
    {
        return star(s).walkLength;
    }

    /// We provide access to the manifold as a simplicial complex
    ref const(SimplicialComplex!Vertex) asSimplicialComplex() const
    {
        return simpComp_; 
    }

    /// We allow the manifold to be treated as a simplicial complex
    alias asSimplicialComplex this;
}

Vertex[][] pachnerMoves(Vertex, int dim)(const SmallManifold!(dim, Vertex) m)
{
    int[Vertex[]] degreeMap;

    // TO DO: Why do I need idup here?
    m.facets.each!(f => f.subsets.each!(s => ++degreeMap[s.idup]));

    Vertex[][] result;
    foreach(simp, deg; degreeMap)
    {
        if(deg == m.dimension + 2 - simp.walkLength)
        {
            auto newSimp = m.asSimplicialComplex.link(simp)
                .joiner.uniq.array.sort;            
            if(newSimp.empty || !m.contains(newSimp))
            {
                result ~= simp.dup;
            }
        }           
    }

    // TO DO: Should we even bother sorting?
    bool mySort(Vertex[] a, Vertex[] b)
    {
        if (a.length < b.length)
        {
            return true;
        }
        else if (a.length == b.length)
        {
            return a < b;
        }
        return false;
    }
    result.sort!mySort;

    return result;
}

void doPachner(Vertex, int dim)(
    SmallManifold!(dim, Vertex) manifold,
    Vertex[] centerFacet,
    Vertex newVertex)
{
    assert(centerFacet.walkLength == manifold.dimension + 1);
    manifold.doPachnerImpl(centerFacet, [newVertex]);
}

void doPachner(Vertex, int dim)(
    SmallManifold!(dim, Vertex) manifold,
    Vertex[] center)
{
    auto coCenter = manifold.link(center).joiner.array.sort().uniq.array;
    assert(!coCenter.empty);
    manifold.doPachnerImpl(center, coCenter);
}

private void doPachnerImpl(Vertex, int dim)(
    SmallManifold!(dim, Vertex) manifold,
    Vertex[] center,
    Vertex[] coCenter
)
{
    // TO DO: Better error
    enforce(manifold.pachnerMoves.canFind(center), "bad pachner move");

    manifold.star(center).each!(f => manifold.simpComp_.removeFacet(f));

    alias SComp = SimplicialComplex!Vertex;
    immutable cDim = center.walkLength.to!int - 1;
    auto newPiece = join(SComp(center.subsetsOfSize(cDim)), SComp([coCenter]));

    assert(newPiece.isPureOfDim(manifold.dimension));
    newPiece.facets.each!(f => manifold.simpComp_.insertFacet(f));

    // TO DO: Do sanity checking for manifold
}
