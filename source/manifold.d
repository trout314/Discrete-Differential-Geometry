import algorithms : eulerCharacteristic, is2Sphere, isCircle,
    isConnected, isPureOfDim, join;
import fluent.asserts;
import simplicial_complex : fVector, SimplicialComplex;
import std.algorithm : all, any, canFind, each, equal, filter, find, joiner,
    map, maxElement, sort, uniq;
import std.conv : to;
import std.range : array, chain, ElementType, empty, enumerate, front, iota,
    isInputRange, popFront, save, walkLength;
import unit_threaded : Name;
import utility : binomial, isSubsetOf, SmallMap, staticIota, subsets, subsetsOfSize, throwsWithMsg;

//dfmt off
///
@Name("Manifold doc tests") @system unittest
{
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
@Name("Manifold (errors)") /* pure */ @system unittest
{
    // TO DO: ldc doesn't like using "pure" above! Bugreport?

    Manifold!2([[1,2,3,4]]).throwsWithMsg(
        "not all facets have the correct dimension");

    Manifold!2([[1,2,3]]).throwsWithMsg(
        "found a ridge with degree not equal to 2");

    Manifold!2([[1,2,3], [1,2,4], [1,3,4], [2,3,4], [1,5,6], [1,5,7], [1,6,7],
        [5,6,7]]).throwsWithMsg("found a hinge whose link is not a circle");

    auto octahedron = [[0,1,2], [0,2,3], [0,3,4], [0,1,4], [1,2,5], [2,3,5],
        [3,4,5], [1,4,5]];

    auto sphere = [[6,7,8], [6,7,9], [6,8,9], [7,8,9]];

    Manifold!2(chain(octahedron, sphere)).throwsWithMsg(
        "facets do not define a connected simplicial complex");

    // These facets separately define 3-spheres, but share the vertex 1
    auto sphere3 = [[1,2,3,4], [1,2,3,5], [1,2,4,5], [1,3,4,5], [2,3,4,5]];
    auto sphere3A =[[1,6,7,8], [1,6,7,9], [1,6,8,9], [1,7,8,9], [6,7,8,9]];

    Manifold!3(chain(sphere3, sphere3A)).throwsWithMsg(
        "found a codimension-3 simplex whose link is not a 2-sphere");
}

// More tests
@Name("Manifold (additional)") /* pure */ @system unittest
{
    // TO DO: Why does this need to be @system? Make it @safe!
    // TO DO: ldc doesn't like using "pure" above! Bugreport?

    static assert(!__traits(compiles, Manifold!2([["a", "bubba", "gump"]])));

    auto m1 = Manifold!(1, string)([["a", "b"], ["b", "c"], ["a", "c"]]);
    static assert(m1.dimension == 1);
    static assert(is(m1.Vertex == string));
}

/*******************************************************************************
Manifold type... TO DO: More info here
*/
struct Manifold(int dimension_, Vertex_ = int)
{
private:
    SimplicialComplex!Vertex simpComp;
    size_t[dimension_ + 1] numSimplices;
    size_t[Vertex[]] degreeMap;
public:
    /// Dimension of the manifold
    static immutable dimension = dimension_;

    /// Vertex type used in the manifold
    alias Vertex = Vertex_;

    static assert(dimension >= 1, "dimension must be at least one, but got "
        ~ "dimension " ~ dimension.to!string);

    /// We can initialize the manifold from a range of ranges of vertices
    this(F)(F initialFacets) if (isInputRange!F)
    {
        // TO DO: Put some nice constraints on F
        foreach(f; initialFacets)
        {
            this.insertFacet(f);
        }

        assert(this.isPureOfDim(dimension),
            "not all facets have the correct dimension");

        assert(this.isConnected,
            "facets do not define a connected simplicial complex");

        static if(dimension >= 1)
        {
            assert(simplices(dimension - 1).find!(s => degree(s) != 2).empty,
                "found a ridge with degree not equal to 2");
        }

        static if(dimension >= 2)
        {
            // TO DO: Figure out why we need "this.link" below (low priority)
            assert(simplices(dimension - 2)
                .find!(s => !SimplicialComplex!Vertex(this.link(s)).isCircle)
                .empty, "found a hinge whose link is not a circle");
        }

        static if(dimension >= 3)
        {
            // TO DO: Figure out why we need "this.link" below (low priority)
            assert(simplices(dimension - 3).find!(
                s => !SimplicialComplex!Vertex(this.link(s)).is2Sphere).empty,
                "found a codimension-3 simplex whose link is not a 2-sphere");
        }

        numSimplices[] = simpComp.fVector[];
    }

    this(this)
    {
        degreeMap = degreeMap.dup;
    }

    /***************************************************************************
    Returns true if and only if the given simplex is in the manifold
    */
    bool contains(V)(V vertices) const if (isInputRange!V)
    {
        assert(this.simpComp.facets.any!(f => vertices.isSubsetOf(f))
            == !((vertices.array in degreeMap) is null));
        return !((vertices.array in degreeMap) is null);
    }


    /***************************************************************************
    Returns the degree of a simplex in the simplicial complex.
    */
    auto degree(V)(V vertices) const if (isInputRange!V)
    {
        assert(vertices.array in degreeMap);
        assert(degreeMap[vertices.array] == star(vertices).walkLength);
        return degreeMap[vertices.array];
    }

    /***************************************************************************
    Returns the degree of a simplex in the simplicial complex.
    */
    const(size_t)[] fVector() const
    {
        // import std.stdio : writeln;
        // writeln(numSimplices[], simpComp.fVector);
        assert(numSimplices[] == this.simpComp.fVector);
        return numSimplices[];
    }

    // Special version of insertFacet to update tracked info
    private void insertFacet(V)(V vertices) if (isInputRange!V)
    {
        assert(vertices.array !in degreeMap);
        foreach(s; vertices.subsets)
        {
            ++degreeMap[s.array.idup];
        }
        this.simpComp.insertFacet(vertices);
    }

    // Special version of removeFacet to update tracked info
    private void removeFacet(V)(V vertices) if (isInputRange!V)
    {
        assert(vertices.array in degreeMap);
        foreach(s; vertices.subsets)
        {
            --degreeMap[s.array.idup];
            if(degreeMap[s.array.idup] == 0)
            {
                degreeMap.remove(s.array.idup);
            }
        }
        assert(vertices.array !in degreeMap);
        this.simpComp.removeFacet(vertices);        
    }

    /// We provide access to the manifold as a simplicial complex
    ref const(SimplicialComplex!Vertex) asSimplicialComplex() const
    {
        return simpComp; 
    }

    alias asSimplicialComplex this;
}

/*******************************************************************************
Returns a list of all the possible pachner moves for the given manifold
*/
const(Vertex)[][] pachnerMoves(Vertex, int dim)(const ref Manifold!(dim, Vertex) m)
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
@Name("pachnerMoves") @system unittest
{
    // TO DO: Why does this need to be @system? Make it @safe!

    auto m = Manifold!2(
        [[1,2,3], [1,2,4], [1,3,4], [2,3,5], [2,4,5],[3,4,5]]);

    m.pachnerMoves.should.containOnly([[1], [2, 3], [2, 3, 5],
        [1, 3, 4], [1, 2, 3], [5], [2, 4, 5], [3, 4, 5], [1, 2, 4], [2, 4],
        [3, 4]]);

    auto octahedron = Manifold!2([[0,1,2], [0,2,3], [0,3,4], [0,1,4],
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
    ref Manifold!(dim, Vertex) manifold,
    const(Vertex)[] centerFacet,
    const(Vertex) newVertex)
{
    assert(centerFacet.walkLength == manifold.dimension + 1);
    manifold.doPachnerImpl(centerFacet, [newVertex]);
}
///
@Name("doPachner 1 -> (dim + 1)") @system unittest
{
    auto m = Manifold!2(4.iota.subsetsOfSize(3).map!(s => s.array.dup));
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
    ref Manifold!(dim, Vertex) manifold,
    const(Vertex)[] center)
{
    // TO DO: Better error
    assert(manifold.pachnerMoves.canFind(center), "bad pachner move");

    auto coCenter = manifold.getCoCenter(center);
    assert(!coCenter.empty);

    manifold.doPachnerImpl(center, coCenter);
}
///
@Name("doPachner n -> (dim + 2 - n), 1 < n < dim + 2") /* pure */ @system unittest
{
    auto octahedron = Manifold!2([[0,1,2], [0,2,3], [0,3,4], [0,1,4],
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
    auto m = Manifold!2([[1,2,3],[1,2,4], [1,3,4], [2,3,4]]);   
    m.doPachner([1,2]).throwsWithMsg("bad pachner move");

    // TO DO: More pachner move tests
}

/*******************************************************************************
Provides a historgram of simplex degrees for simplices of dimension 'dim' by
returning an associative array whose values are the number of simplices and 
whose keys are the degrees.
*/
int[int] degreeHistogram(Vertex, int dim)(
    const ref Manifold!(dim, Vertex) manifold,
    int histogramDim)
{
    // TO DO: Improve this

    assert(histogramDim >= 0);
    assert(histogramDim <= dim);

    int[int] result;
    int[Vertex[]] degreeMap;

    manifold.facets
        .each!(f => f.subsetsOfSize(histogramDim + 1)
        .each!(s => ++degreeMap[s.array]));

    degreeMap.byValue.each!(deg => ++result[deg]);  
    return result;
}
///
@Name("degreeHistogram") /* pure */ @system unittest
{
    // TO DO: ldc doesn't like using "pure" above! Bugreport?

    auto m = Manifold!2(
        [[1,2,3], [1,2,4], [1,3,4], [2,3,5], [2,4,5],[3,4,5]]);
    assert(m.degreeHistogram(0) == [4:3, 3:2]);
    assert(m.degreeHistogram(1) == [2:9]);
    assert(m.degreeHistogram(2) == [1:6]);

    // TO DO: Some more tests
}

///
auto pachnerMovesAndDegreeHistogram(Vertex, int dim)(
    const ref Manifold!(dim, Vertex) manifold,
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
@Name("pachnerMovesAndDegreeHistogram") @system unittest
{
    auto m = Manifold!2(
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
    ref Manifold!(dim, Vertex) manifold,
    const(Vertex)[] center,
    const(Vertex)[] coCenter
)
{
    assert(manifold.pachnerMoves.canFind(center), "bad pachner move");
    immutable centerDim = center.walkLength.to!int - 1;

    /* Need to ensure independent copies of the facets in the star since once a
    facet is removed, the range returned by star(center) becomes invalid! */   
    immutable toRemove = manifold.star(center).map!(f => f.dup).array;
    foreach(f; toRemove)
    {
        manifold.removeFacet(f);
    }
    // toRemove.each!(f => manifold.removeFacet(f));

    alias SC = SimplicialComplex!Vertex;
    auto newPiece = (centerDim == 0)
        ? SC([coCenter])
        : join(SC(center.subsetsOfSize(centerDim)), SC([coCenter]));

    assert(newPiece.isPureOfDim(dim));
    newPiece.facets.each!(f => manifold.insertFacet(f));

    manifold.numSimplices[].modifyFVector(center.length);

    // TO DO: Do sanity checking for manifold
}

// TO DO: Implement 4->4 moves (Important!)

private auto getCoCenter(Vertex, int dim)(
    const ref Manifold!(dim, Vertex) manifold,
    const(Vertex)[] center
)
{
    assert(manifold.contains(center));
    return manifold.link(center).joiner.array.dup.sort().uniq.array;
}

// TO DO: Separate unittesting for getCoCenter

// TO DO: Use standardSphereFacets below in appropriate places

/*******************************************************************************
Returns a lazy range that goes through the facets of the "standard" sphere
of dimension `dim`. This manifold is just the boundary if the "standard"
(dim+1)-simplex [0,1,2, ... (dim+1)].
*/
auto standardSphereFacets(int dim)
{
    assert(dim >= 1);
    return iota(dim + 2).subsetsOfSize(dim + 1).map!array;
}
///
@Name("standardSphereFacets") @system unittest
{
    import std.stdio : writeln;
    standardSphereFacets(2).should.containOnly(
        [[0,1,2], [0,1,3], [0,2,3], [1,2,3]]);

    auto s = Manifold!2(standardSphereFacets(2));
    s.facets.should.containOnly([[0,1,2], [0,1,3], [0,2,3], [1,2,3]]);

    foreach(dim; staticIota!(1,9))
    {
        auto m = Manifold!dim(standardSphereFacets(dim));
        assert(m.numFacets == dim + 2);
    }
}

///
@Name("facets(dim) (pure nothrow @nogc @safe)") @system unittest
{
    auto sc = Manifold!1([[0,1],[0,2],[1,2]]);
    int[2] edge01 = [0,1];
    int[2] edge02 = [0,2];
    int[2] edge12 = [1,2];

    () pure nothrow @nogc @safe {
        auto facetsRange = sc.facets(1);
        auto savedRange = facetsRange.save;

        int[2][3] facs;
        facs[0] = facetsRange.front;
        facetsRange.popFront;
        facs[1] = facetsRange.front;
        facetsRange.popFront;
        facs[2] = facetsRange.front;
        facetsRange.popFront;
        
        facs[].sort();

        assert(facs[0] == edge01);
        assert(facs[1] == edge02);
        assert(facs[2] == edge12);

        assert(facetsRange.empty);    
        assert(!savedRange.empty);    
    }();
}

///
@Name("facets() (pure nothrow @nogc @safe)") @system unittest
{
    auto sc = Manifold!1([[0,1],[0,2],[1,2]]);
    int[2] edge01 = [0,1];
    int[2] edge02 = [0,2];
    int[2] edge12 = [1,2];

    () pure nothrow @nogc @safe {
        auto facetsRange = sc.facets;
        auto savedRange = facetsRange.save;

        int[2][3] facs;
        facs[0] = facetsRange.front;
        facetsRange.popFront;
        facs[1] = facetsRange.front;
        facetsRange.popFront;
        facs[2] = facetsRange.front;
        facetsRange.popFront;
        
        facs[].sort();

        assert(facs[0] == edge01);
        assert(facs[1] == edge02);
        assert(facs[2] == edge12);

        assert(facetsRange.empty);    
        assert(!savedRange.empty);    
    }();
}

///
@Name("star(range) (pure nothrow @nogc @safe)") @system unittest
{
    auto sc = Manifold!1([[0,1],[0,2],[1,2]]);
    int[2] edge01 = [0,1];
    int[2] edge12 = [1,2];
    int[1] vertex1 = [1];

    () pure nothrow @nogc @safe {
        auto starRange = sc.star(vertex1[]);
        auto savedRange = starRange.save;

        int[2][2] edges;
        edges[0][] = starRange.front[];
        starRange.popFront;
        edges[1][] = starRange.front[];
        starRange.popFront;

        edges[].sort();        
        assert(edges[0] == edge01);
        assert(edges[1] == edge12);

        assert(starRange.empty);    
        assert(!savedRange.empty);    
    }();
}

///
@Name("link(range) (pure nothrow @nogc @safe)") @system unittest
{
    auto sc = Manifold!1([[0,1],[0,2],[1,2]]);
    int[1] vertex0 = [0];
    int[1] vertex1 = [1];
    int[1] vertex2 = [2];

    () pure nothrow @nogc @safe {
        auto linkRange = sc.link(vertex1[]);
        auto savedRange = linkRange.save;

        int[2] vertices;
        vertices[0] = linkRange.front.front;
        linkRange.popFront;
        vertices[1] = linkRange.front.front;
        linkRange.popFront;

        vertices[].sort();        
        assert(vertices[0] == 0);
        assert(vertices[1] == 2);

        assert(linkRange.empty);    
        assert(!savedRange.empty);    
    }();
}

///
@Name("contains (pure nothrow /* @nogc */ @safe)") @system unittest
{
    auto sc = Manifold!1([[0,1],[0,2],[1,2]]);
    int[2] edge01 = [0,1];
    int[2] edge07 = [0,7];

    () pure nothrow /* @nogc */ @safe {
        assert(sc.contains(edge01[]));
        assert(!sc.contains(edge07[]));
    }();  
}

///
// @Name("removeFacet unavailable") @system unittest
// {
//     auto sc = Manifold!1([[0,1],[0,2],[1,2]]);
//     static assert(!__traits(compiles, sc.removeFacet([0,1])));    
// }

///
// @Name("insertFacet unavailable") @system unittest
// {
//     auto sc = Manifold!1([[0,1],[0,2],[1,2]]);
//     static assert(!__traits(compiles, sc.insertFacet([0,3])));    
// }

/******************************************************************************
Modifies an fVector for a pachner move with center simplex that contains `n`
vertices.
*/
void modifyFVector(size_t[] fVector_, size_t centerLength)
{
    assert(centerLength >= 1, "center length must be at least one");
    assert(fVector_.length >= centerLength, "fVector must have length at least the center length");

    // We modify the fVector for the removal of the original star
    auto centerDim = centerLength - 1;
    auto dim = fVector_.length - 1;
    foreach(d; centerDim .. dim + 1)
    {
        fVector_[d] -= binomial(dim + 1 - centerDim, d - centerDim);
    }

    // Now modify it for the addition of the new star
    auto coDim = dim + 1 - centerLength;
    foreach(d; coDim .. dim + 1)
    {
        fVector_[d] += binomial(dim + 1 - coDim, d - coDim);
    }
}
///
@Name("modifyFVector") unittest
{
    size_t[4] fVec = [0,0,0,0];

    /* A 1 -> 4 move should give net: +1 vertex, +4 edges, +6 triangles,
    +3 tetrahdra */
    fVec[].modifyFVector(4);
    fVec.should.equal([1, 4, 6, 3]);
}

// TO DO: Adapt old code below for new manifold type!

// /******************************************************************************
// * Returns a manifold loaded from the file specified by fileName. If fileName
// * is the empty string, the returned manifold is the standard sphere.
// */
// Manifold!dim loadManifold(size_t dim = dimManifold)(string fileName)
// {
//     auto manifoldFile = File(fileName, "r"); // Open file in read-only mode

//     string facets;
//     foreach (string line; manifoldFile.lines)
//     {
//         // We allow comments starting with '#'
//         if (line.front != '#')
//         {
//             facets ~= line.strip;
//         }
//     }

//     try
//     {
//         auto loaded = Manifold!dim(facets.parse!(Simplex[]));
//         return loaded;
//     }
//     catch (Exception ex)
//     {
//         ex.msg ~= "\n\n\t ERROR: encountered malformed facet list in initial manifold triangulation file: "
//             ~ fileName ~ "\n";
//         throw ex;
//     }
// }
// ///
// unittest
// {

//     auto manifold = loadManifold!2("data/manifold_sampler_unittest_load.dat");
//     auto expectedManifold = Manifold!2([[0, 1, 2], [0, 2, 3], [0, 3, 4], [0, 4,
//             5], [0, 5, 6], [0, 6, 1], [1, 2, 6], [2, 3, 5], [2, 5, 6], [3, 4, 5]]);
//     assert(manifold == expectedManifold);

//     assertThrown(loadManifold!2("data/manifold_sampler_unittest_bad_load.dat"));
// }

// /******************************************************************************
// * Saves a manifold to file specified by fileName.
// */
// void saveTo(size_t dimension)(Manifold!dimension manifold, string fileName)
// {
//     auto saveFile = File(fileName, "w"); // Open in write-only mode

//     saveFile.writeln("# created ", Clock.currTime.to!DateTime);
//     saveFile.write("[");
//     foreach (record; manifold.facetRecords[0 .. $ - 1])
//     {
//         saveFile.write(record.facet, ",");
//     }
//     saveFile.writeln(manifold.facetRecords.back.facet, "]");
// }
// ///
// unittest
// {

//     auto fileName = "data/manifold_sampler_unittest_save.dat";
//     auto sphere = Manifold!4(standardSphere(4));
//     sphere.saveTo(fileName);
//     Manifold!4 loadedSphere = loadManifold!4(fileName);
//     assert(loadedSphere == sphere);
// }

// unittest
// {

//     auto octahedron = Manifold!2([[0, 1, 2], [0, 2, 3], [0, 3, 4], [0, 4, 1],
//             [5, 1, 2], [5, 2, 3], [5, 3, 4], [5, 4, 1]]);

//     assert(octahedron.fVector == [6, 12, 8]);

//     auto sphere7 = Manifold!7(standardSphere(7));

//     // fVector should be (9 choose 1), (9 choose 2), ... , (9 choose 8)
//     assert(sphere7.fVector == [9, 9 * 8 / 2, 9 * 8 * 7 / (3 * 2),
//             9 * 8 * 7 * 6 / (4 * 3 * 2), 9 * 8 * 7 * 6 / (4 * 3 * 2),
//             (9 * 8 * 7) / (3 * 2), (9 * 8) / 2, 9]);

//     auto hyperbolicDodecahedral = loadManifold!3("data/manifold_sampler_unittest_dodecahedral.dat");

//     // TO DO: Get reference for this...
//     assert(hyperbolicDodecahedral.fVector == [21, 190, 338, 169]);

// }
