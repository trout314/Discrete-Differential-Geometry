import std.algorithm : all, any, canFind, each, equal, filter, find, joiner,
    map, maxElement, sort, uniq;
import std.conv : to;
import std.exception : assertThrown;
import std.range : array, chain, empty, enumerate, front, iota, isInputRange,
    ElementType, popFront, walkLength;
import std.traits : isArray;
import simplicial_complex : simplicialComplex, SimplicialComplex, fVector;
import utility : SmallMap, staticIota, subsets, subsetsOfSize, throwsWithMsg;
import std.stdio : writeln;

import simplicial_complex_algorithms : connectedComponents, is2Sphere, isCircle,
    isConnected, join, eulerCharacteristic, isPureOfDim;

import fluent.asserts;
import unit_threaded : Name, writelnUt;

//dfmt off
///
@Name("SmallManifold doc tests") unittest
{
    alias Manifold = SmallManifold;

    auto octahedron = Manifold!2([[0,1,2], [0,2,3], [0,3,4], [0,1,4], [1,2,5],
        [2,3,5], [3,4,5], [1,4,5]]);

    static assert(octahedron.dimension == 2);
    static assert(is(octahedron.Vertex == int));

    assert(octahedron.fVector == [6,12,8]);
    assert(octahedron.eulerCharacteristic == 2);

    auto tetrahedron = Manifold!2([[1,2,3], [1,2,4], [1,3,4], [2,3,4]]);

    assert(octahedron.star([1,2]).equal([[0,1,2], [1,2,5]]));    

    octahedron.pachnerMoves.should.containOnly(
        [[0,1], [0,2], [0,3], [0,4], [1,2], [1,4],
         [1,5], [2,3], [2,5], [3,4], [3,5], [4,5],  // 1-simplices
         [0,1,2], [0,1,4], [0,2,3], [0,3,4],
         [1,2,5], [1,4,5], [2,3,5], [3,4,5]]);      // 2-simplices

    tetrahedron.pachnerMoves.should.containOnly(
        [[1, 2, 3], [1, 2, 4], [1, 3, 4], [2, 3, 4]]);

    tetrahedron.doPachner([1,2,3], 5);

    // TO DO: FINISH CHECKS!
    // octahedron.doPachner([1,2]);
    // octahedron.doPachner([0,5]);
}

///
@Name("SmallManifold (errors)") pure @system unittest
{
    alias Manifold = SmallManifold;

    Manifold!2([[1,2,3,4]]).throwsWithMsg!Error(
        "not all facets have the correct dimension");

    Manifold!2([[1,2,3]]).throwsWithMsg!Error(
        "found a ridge with degree not equal to 2");

    Manifold!2([[1,2,3], [1,2,4], [1,3,4], [2,3,4], [1,5,6], [1,5,7], [1,6,7],
        [5,6,7]]).throwsWithMsg!Error("found a hinge whose link is not a circle");

    auto octahedron = [[0,1,2], [0,2,3], [0,3,4], [0,1,4], [1,2,5], [2,3,5],
        [3,4,5], [1,4,5]];

    auto sphere = [[6,7,8], [6,7,9], [6,8,9], [7,8,9]];

    Manifold!2(chain(octahedron, sphere)).throwsWithMsg!Error(
        "facets do not define a connected simplicial complex");

    // These facets separately define 3-spheres, but share the vertex 1
    auto sphere3 = [[1,2,3,4], [1,2,3,5], [1,2,4,5], [1,3,4,5], [2,3,4,5]];
    auto sphere3A =[[1,6,7,8], [1,6,7,9], [1,6,8,9], [1,7,8,9], [6,7,8,9]];

    Manifold!3(chain(sphere3, sphere3A)).throwsWithMsg!Error(
        "found a codimension-3 simplex whose link is not a 2-sphere");
}

// More tests
@Name("SmallManifold (additional)") unittest
{
    alias Manifold = SmallManifold;

    static assert(!__traits(compiles, Manifold!2([["a", "bubba", "gump"]])));

    auto m1 = Manifold!(1, string)([["a", "b"], ["b", "c"], ["a", "c"]]);
    static assert(m1.dimension == 1);
    static assert(is(m1.Vertex == string));
}

struct SmallManifold(int dimension_, Vertex_ = int)
{
private:
    SimplicialComplex!Vertex simpComp_;
public:
    static immutable dimension = dimension_;
    alias Vertex = Vertex_;

    static assert(dimension >= 1, "dimension must be at least one, but got "
        ~ "dimension " ~ dimension.to!string);

    this(F)(F initialFacets)
    {
        // TO DO: Put some nice constraints on F

        // static assert((isInputRange!F || isArray!F)
        //     && is(ElementType!F : Facet) || is(ElementType!F : Vertex[]),
        //     "a " ~ SmallManifold!(dimension, Vertex).stringof ~ " must be "
        //     ~ "constructed from a range or array of elements with type " 
        //     ~ Facet.stringof ~ " or " ~ Vertex[].stringof ~ " not " ~ F.stringof);

        initialFacets.each!(f => simpComp_.insertFacet(f));

        assert(this.isPureOfDim(dimension),
            "not all facets have the correct dimension");

        assert(this.isConnected,
            "facets do not define a connected simplicial complex");

        static if(dimension >= 1)
        {{
            auto badRidge = this.simplices(dimension - 1).find!(
                s => this.degree(s) != 2);
            assert(badRidge.empty, "found a ridge with degree not equal to 2");
        }}

        static if(dimension >= 2)
        {{
            auto badHinge = this.simplices(dimension - 2).find!(
                s => !this.link(s).to!(SimplicialComplex!Vertex).isCircle);
            assert(badHinge.empty, "found a hinge whose link is not a circle");
        }}

        static if(dimension >= 3)
        {{
            auto badCodim3 = this.simplices(dimension - 3).find!(
                s => !this.link(s).to!(SimplicialComplex!Vertex).is2Sphere);
            assert(badCodim3.empty,
                "found a codimension-3 simplex whose link is not a 2-sphere");
        }}
    }

    auto degree(V)(V vertices) const if (isInputRange!V)
    {
        return star(vertices).walkLength;
    }

    /// We provide access to the manifold as a simplicial complex
    ref const(SimplicialComplex!Vertex) asSimplicialComplex() const
    {
        return simpComp_; 
    }

    alias asSimplicialComplex this;
}

/*******************************************************************************
Returns a list of all the possible pachner moves for the given manifold
*/
const(Vertex)[][] pachnerMoves(Vertex, int dim)(const ref SmallManifold!(dim, Vertex) m)
{
    int[const(Vertex)[]] degreeMap;

    m.facets.each!(f => f.subsets.each!(s => ++degreeMap[s.array]));

    const(Vertex)[][] result;
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
///
@Name("pachnerMoves") unittest
{
    auto m = SmallManifold!2(
        [[1,2,3], [1,2,4], [1,3,4], [2,3,5], [2,4,5],[3,4,5]]);

    m.pachnerMoves.should.containOnly([[1], [2, 3], [2, 3, 5],
        [1, 3, 4], [1, 2, 3], [5], [2, 4, 5], [3, 4, 5], [1, 2, 4], [2, 4],
        [3, 4]]);

    auto octahedron = SmallManifold!2([[0,1,2], [0,2,3], [0,3,4], [0,1,4],
        [1,2,5], [2,3,5], [3,4,5], [1,4,5]]);
    
    assert(octahedron.simplices(0).all!(s => octahedron.degree(s) == 4));

    octahedron.pachnerMoves.should.containOnly(
        [[0,1,2], [0,2,3], [0,3,4], [0,1,4], [1,2,5], [2,3,5], [3,4,5], [1,4,5],
        [2, 3], [0, 1], [1, 5], [4, 5], [0, 3], [1, 4],
        [1, 2], [0, 4], [0, 2], [2, 5], [3, 5], [3, 4]]);

}

/*******************************************************************************
Does a pachner move of type 1 -> (dim + 1) with the new vertex given by
the user.
*/
void doPachner(Vertex, int dim)(
    ref SmallManifold!(dim, Vertex) manifold,
    const(Vertex)[] centerFacet,
    const(Vertex) newVertex)
{
    assert(centerFacet.walkLength == manifold.dimension + 1);
    manifold.doPachnerImpl(centerFacet, [newVertex]);
}
///
@Name("doPachner 1 -> (dim + 1)") unittest
{
    auto m = SmallManifold!2(4.iota.subsetsOfSize(3).map!(s => s.array.dup));
    m.doPachner([1,2,3], 4);
    m.facets.should.containOnly(
        [[0,1,2],[0,1,3], [0,2,3], [1,2,4], [1,3,4], [2,3,4]]);

    m.doPachner([0,2,3], 7);
    m.facets.should.containOnly([[0,1,2], [0,1,3], [0,2,7],
        [0,3,7], [1,2,4], [1,3,4], [2,3,4], [2,3,7]]);

    // TO DO: More pachner move tests
}

/*******************************************************************************
Does a pachner move which replaces the star of the given simplex. This is a
move of type (dim + 2 - center.length) -> center.length
*/
void doPachner(Vertex, int dim)(
    ref SmallManifold!(dim, Vertex) manifold,
    const(Vertex)[] center)
{
    // TO DO: Better error
    assert(manifold.pachnerMoves.canFind(center), "bad pachner move");

    auto coCenter = manifold.getCoCenter(center);
    assert(!coCenter.empty);

    manifold.doPachnerImpl(center, coCenter);
}
///
@Name("doPachner n -> (dim + 2 - n), 1 < n < dim + 2") unittest
{
    auto octahedron = SmallManifold!2([[0,1,2], [0,2,3], [0,3,4], [0,1,4],
        [1,2,5], [2,3,5], [3,4,5], [1,4,5]]);

    octahedron.doPachner([1,2]);
    assert(octahedron.degree([0]) == 5);
    assert(octahedron.degree([5]) == 5);
    assert(octahedron.degree([1]) == 3);
    assert(octahedron.degree([2]) == 3);

    // We can undo the 2->2 move
    octahedron.doPachner([0,5]);
    assert(octahedron.simplices(0).all!(s => octahedron.degree(s) == 4));

    octahedron.doPachner([0,1,2], 99);
    octahedron.doPachner([99]);

    // Can't do 2->2 move on the boundary of a 3-simplex
    auto m = SmallManifold!2([[1,2,3],[1,2,4], [1,3,4], [2,3,4]]);
    
    // m.doPachner([1,2]).assertThrown!Error;

    m.doPachner([1,2]).throwsWithMsg!Error("bad pachner move");

    // TO DO: More pachner move tests
}

int[int] degreeHistogram(Vertex, int dim)(
    const ref SmallManifold!(dim, Vertex) manifold,
    int histogramDim)
{
    // TO DO: Improve this

    assert(histogramDim >= 0);
    assert(histogramDim <= dim);

    int[int] result;
    int[Vertex[]] degreeMap;

    manifold.facets.each!(
        f => f.subsetsOfSize(histogramDim + 1).each!(s => ++degreeMap[s.array]));

    degreeMap.byValue.each!(deg => ++result[deg]);  
    return result;
}
///
@Name("degreeHistogram") unittest
{
    auto m = SmallManifold!2(
        [[1,2,3], [1,2,4], [1,3,4], [2,3,5], [2,4,5],[3,4,5]]);
    assert(m.degreeHistogram(0) == [4:3, 3:2]);
    assert(m.degreeHistogram(1) == [2:9]);
    assert(m.degreeHistogram(2) == [1:6]);

    // TO DO: Some more tests
}

auto pachnerMovesAndDegreeHistogram(Vertex, int dim)(
    const ref SmallManifold!(dim, Vertex) manifold,
    int histogramDim)
{
    assert(histogramDim >= 0);
    assert(histogramDim <= dim);

    int[Vertex[]] degreeMap;    // Stores degree of each simplex

    manifold.facets.each!(
        f => f.subsets.each!(s => ++degreeMap[s.array]));

    Vertex[][] moves_;
    int[int] histogram_;
    foreach(simp, deg; degreeMap)
    {
        // Add Pachner move with center `simp` if appropriate 
        if(deg == manifold.dimension + 2 - simp.walkLength)
        {
            auto coCenter = manifold.getCoCenter(simp);
            if(coCenter.empty || (!manifold.contains(coCenter)))
            {
                moves_ ~= simp.dup;
            }
        }

        // Increment histogram if needed
        if(simp.walkLength == histogramDim + 1)
        {
            ++histogram_[deg];
        }
    }
    
    static struct Result
    {
        Vertex[][] moves;
        int[int] histogram;
    }

    return Result(moves_, histogram_);
}
///
@Name("pachnerMovesAndDegreeHistogram") unittest
{
    auto m = SmallManifold!2(
        [[1,2,3], [1,2,4], [1,3,4], [2,3,5], [2,4,5],[3,4,5]]);

    auto r0 = m.pachnerMovesAndDegreeHistogram(0);
    auto r1 = m.pachnerMovesAndDegreeHistogram(1);
    auto r2 = m.pachnerMovesAndDegreeHistogram(2);

    foreach(r; [r0, r1, r2])
    {
        r.moves.should.containOnly([
            [1], [5],
            [2, 3], [2, 4], [3, 4],
            [2, 3, 5], [1, 3, 4], [1, 2, 3], [2, 4, 5], [3, 4, 5], [1, 2, 4]
        ]);
    }

    assert(r0.histogram == [4:3, 3:2]);
    assert(r1.histogram == [2:9]);
    assert(r2.histogram == [1:6]);

    // TO DO: More tests

    /* TO DO: Is this actually faster than two separate calls to pachnerMoves
    and degreeHistogram? */
}



// Factor out the common code for the two types of doPachner
private void doPachnerImpl(Vertex, int dim)(
    ref SmallManifold!(dim, Vertex) manifold,
    const(Vertex)[] center,
    const(Vertex)[] coCenter
)
{
    // TO DO: Better error
    assert(manifold.pachnerMoves.canFind(center), "bad pachner move");

    // need to make sure we have independent copies of the facets in the star   
    auto toRemove = manifold.star(center).map!(f => f.dup).array;
    toRemove.each!(f => manifold.simpComp_.removeFacet(f));

    alias SComp = SimplicialComplex!Vertex;
    immutable cDim = center.walkLength.to!int - 1;
   
    //  TO DO: Create SimplicialComplex constructor that can directly take subsetsOfSize
    auto newPiece = (cDim == 0)
        ? SComp([coCenter])
        : join(SComp(center.subsetsOfSize(cDim).map!(s => s.array.dup)), SComp([coCenter]));


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
    return manifold.link(center).joiner.array.dup.sort().uniq.array;
}

// TO DO: Separate unittesting for getCoCenter
