import std.algorithm : all, any, canFind, each, equal, filter, find, joiner,
    map, maxElement, sort, uniq;
import std.conv : to;
import std.exception : enforce;
import std.range : array, empty, enumerate, front, iota, isInputRange, ElementType,
    popFront, walkLength;
import std.traits : isArray;
import simplex : Simplex, simplex;
import simplicial_complex : SimplicialComplex, fVector;
import utility : staticIota, subsets, subsetsOfSize, throwsWithMsg;
import std.stdio : writeln;

import simplicial_complex_algorithms : connectedComponents, is2Sphere, isCircle,
    isConnected, join, eulerCharacteristic, isPureOfDim;

import fluent.asserts;
import unit_threaded : Name, writelnUt;

//dfmt off
///
@Name("SmallManifold doc tests")
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
         [1,5], [2,3], [2,5], [3,4], [3,5], [4,5],  // 1-simplices
         [0,1,2], [0,1,4], [0,2,3], [0,3,4],
         [1,2,5], [1,4,5], [2,3,5], [3,4,5]]);      // 2-simplices

    tetrahedron.pachnerMoves.map!(s => s.idup).should.containOnly(
        [[1, 2, 3], [1, 2, 4], [1, 3, 4], [2, 3, 4]]);

    // TO DO: Why do I need the map!(s => s.idup) part above?!

    tetrahedron.doPachner([1,2,3], 5);

    // TO DO: FINISH CHECKS!
    // octahedron.doPachner([1,2]);
    // octahedron.doPachner([0,5]);
}

// More tests
@Name("SmallManifold (additional)")
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

/*******************************************************************************
Returns a list of all the possible pachner moves for the given manifold
*/
Vertex[][] pachnerMoves(Vertex, int dim)(const ref SmallManifold!(dim, Vertex) m)
{
    int[Vertex[]] degreeMap;

    // TO DO: Why do I need idup here?
    m.facets.each!(f => f.subsets.each!(s => ++degreeMap[s.idup]));

    Vertex[][] result;
    foreach(simp, deg; degreeMap)
    {
        if(deg == m.dimension + 2 - simp.walkLength)
        {
            auto coCenter = m.getCoCenter(simp);
            if(coCenter.empty || (!m.contains(coCenter)))
            {
                result ~= simp.dup;
            }
        }           
    }
    return result;
}

/*******************************************************************************
Does a pachner move of type 1 -> (dim + 1) with the new vertex given by
the user.
*/
void doPachner(Vertex, int dim)(
    ref SmallManifold!(dim, Vertex) manifold,
    Vertex[] centerFacet,
    Vertex newVertex)
{
    assert(centerFacet.walkLength == manifold.dimension + 1);
    manifold.doPachnerImpl(centerFacet, [newVertex]);
}
///
@Name("doPachner 1 -> (dim + 1)")
unittest
{
    auto m = SmallManifold!2(4.iota.subsetsOfSize(3));
    m.doPachner([1,2,3], 4);
    assert(m.facets == [[0,1,2],[0,1,3], [0,2,3], [1,2,4], [1,3,4], [2,3,4]]);

    m.doPachner([0,2,3], 7);
    assert(m.facets == [[0,1,2], [0,1,3], [0,2,7], [0,3,7], [1,2,4], [1,3,4],
        [2,3,4], [2,3,7]]);

    // TO DO: More pachner move tests
}

/*******************************************************************************
Does a pachner move which replaces the star of the given simplex. This is a
move of type (dim + 2 - center.length) -> center.length
*/
void doPachner(Vertex, int dim)(
    ref SmallManifold!(dim, Vertex) manifold,
    Vertex[] center)
{
    auto coCenter = manifold.getCoCenter(center);
    assert(!coCenter.empty);
    manifold.doPachnerImpl(center, coCenter);
}
///
@Name("doPachner n -> (dim + 2 - n), 1 < n < dim + 2")
unittest
{
    auto octahedron = SmallManifold!2([[0,1,2], [0,2,3], [0,3,4], [0,1,4], [1,2,5],
        [2,3,5], [3,4,5], [1,4,5]]);
    octahedron.doPachner([1,2]);
    assert(octahedron.degree(simplex(0)) == 5);

    // We can undo the 2->2 move
    octahedron.doPachner([0,5]);
    assert(octahedron.simplices!0.all!(s => octahedron.degree(s) == 4));

    octahedron.doPachner([0,1,2], 99);
    octahedron.doPachner([99]);

    // Can't do 2->2 move on the boundary of a 3-simplex
    auto m = SmallManifold!2(4.iota.subsetsOfSize(3));
    m.doPachner([1,2]).throwsWithMsg("bad pachner move");

    // TO DO: More pachner move tests
}

// Factor out the common code for the two types of doPachner
private void doPachnerImpl(Vertex, int dim)(
    ref SmallManifold!(dim, Vertex) manifold,
    Vertex[] center,
    Vertex[] coCenter
)
{
    // TO DO: Better error
    enforce(manifold.pachnerMoves.canFind(center), "bad pachner move");
    manifold.star(center).each!(f => manifold.simpComp_.removeFacet(f));

    alias SComp = SimplicialComplex!Vertex;
    immutable cDim = center.walkLength.to!int - 1;
   
    auto newPiece = (cDim == 0)
        ? SComp([coCenter])
        : join(SComp(center.subsetsOfSize(cDim)), SComp([coCenter]));

    assert(newPiece.isPureOfDim(manifold.dimension));
    newPiece.facets.each!(f => manifold.simpComp_.insertFacet(f));

    // TO DO: Do sanity checking for manifold
}

private auto getCoCenter(Vertex, int dim)(
    const ref SmallManifold!(dim, Vertex) manifold,
    const(Vertex)[] center
)
{
    assert(manifold.contains(center));
    return manifold.link(center).joiner.array.sort().uniq.array;
}

