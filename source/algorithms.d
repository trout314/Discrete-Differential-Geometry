import manifold : Manifold;
import simplicial_complex : fVector, simplicialComplex, SimplicialComplex;
import std.algorithm : all, any, canFind, chunkBy, equal, filter, find, joiner, map,
    setIntersection, sort, sum;
import std.conv : to;
import std.range : array, empty, enumerate, front, iota, popFront, save,
    walkLength;
import unit_threaded : Name;
import utility : subsetsOfSize, throwsWithMsg;


/*******************************************************************************
Returns true if the given manifold is orientable and false otherwise.
*/
bool isOrientable(Vertex, int dim)(Manifold!(dim, Vertex) manifold)
{
    /* We must choose a compatible orientation for each facet. Since the facets
    already come equipped with an ordering for the vertices, we need only
    indicate whether this orientation is the one we want (GivenOrder) or not
    (OppositeOrder). We let NotSet indicate that the facet has not yet been
    processed. */
    enum Orientation
    {
        NotSet,
        GivenOrder,
        OppositeOrder
    }

    static struct FacetRecord
    {
        const(Vertex)[] facet;
        Orientation label;
        bool done;
    }

    // An empty manifold is orientable
    if(manifold.numFacets == 0)
    {
        return true;
    }

    auto records = manifold.facets.map!(f => FacetRecord(f)).array;

    // We (arbitrarily) use the given vertex order to orient the first simplex
    assert(records.walkLength > 0);
    records.front.label = Orientation.GivenOrder;

    while(records.any!(r => !r.done))
    {
        /* If we're not done there must be some oriented facet that hasn't had
        its neighbors orientations set (or checked OK if already set) */
        auto toDo = records.find!(r => !r.done && r.label != Orientation.NotSet);
        assert(toDo.walkLength > 0);

        foreach(i, ridge; toDo.front.facet.subsetsOfSize(dim).map!array.enumerate)
        {
            auto oppFacet = manifold.star(ridge).filter!(
                f => f != toDo.front.facet);
            assert(!oppFacet.empty);

            auto j = oppFacet.front.subsetsOfSize(dim).map!array.enumerate
                .find!(p => p.value == ridge).front.index;

            Orientation oppFacetLabel;
            if(i % 2 != j % 2)
            {
                oppFacetLabel = toDo.front.label;
            }
            else
            {
                if(toDo.front.label == Orientation.GivenOrder)
                {
                    oppFacetLabel = Orientation.OppositeOrder;
                }
                else
                {
                    oppFacetLabel = Orientation.GivenOrder;
                }
            }

            auto oppRecord = records.find!(r => r.facet == oppFacet.front);
            if(oppRecord.front.label != Orientation.NotSet)
            {
                if (oppRecord.front.label != oppFacetLabel)
                {
                    return false;
                }
            }
            else
            {
                oppRecord.front.label = oppFacetLabel;
            }           
        }
        toDo.front.done = true;
    }
    return true;
}
///
@Name("isOrientable") /* pure */ @system unittest
{
    // TO DO: Why does this need to be @system? Make it @safe!
    // TO DO: ldc doesn't like using "pure" above! Bugreport?

    assert(Manifold!2().isOrientable);

    // http://page.math.tu-berlin.de/~lutz/stellar/manifolds_lex/manifolds_lex_d2_n10_o0_g5
    // Surface #4941 on non-orientable genus 5 list
    auto g5 = Manifold!2([[1, 2, 3], [1, 2, 4], [1, 3, 5], [1, 4, 6], [1,
            5, 7], [1, 6, 7], [2, 3, 6], [2, 4, 8], [2, 5, 7], [2, 5, 9], [2,
            6, 10], [2, 7, 8], [2, 9, 10], [3, 4, 6], [3, 4, 10], [3, 5, 8],
            [3, 7, 8], [3, 7, 10], [4, 5, 9], [4, 5, 10], [4, 8, 9], [5, 8,
            10], [6, 7, 9], [6, 8, 9], [6, 8, 10], [7, 9, 10]]);
    assert(!g5.isOrientable);

    // http://page.math.tu-berlin.de/~lutz/stellar/manifolds_lex/manifolds_lex_d2_n9_o1_g0
    // Surface #15 on orientable genus 0 list
    auto g0 = Manifold!2([[1, 2, 3], [1, 2, 4], [1, 3, 4], [2, 3, 5], [2,
            4, 5], [3, 4, 6], [3, 5, 7], [3, 6, 7], [4, 5, 8], [4, 6, 7], [4,
            7, 9], [4, 8, 9], [5, 7, 8], [7, 8, 9]]);
    assert(g0.isOrientable);

    // http://page.math.tu-berlin.de/~lutz/stellar/library_of_triangulations/poincare
    auto poincare = Manifold!3([[1, 2, 4, 9], [1, 2, 4, 15], [1, 2, 6,
            14], [1, 2, 6, 15], [1, 2, 9, 14], [1, 3, 4, 12], [1, 3, 4, 15],
            [1, 3, 7, 10], [1, 3, 7, 12], [1, 3, 10, 15], [1, 4, 9, 12], [1, 5,
            6, 13], [1, 5, 6, 14], [1, 5, 8, 11], [1, 5, 8, 13], [1, 5, 11,
            14], [1, 6, 13, 15], [1, 7, 8, 10], [1, 7, 8, 11], [1, 7, 11, 12],
            [1, 8, 10, 13], [1, 9, 11, 12], [1, 9, 11, 14], [1, 10, 13, 15],
            [2, 3, 5, 10], [2, 3, 5, 11], [2, 3, 7, 10], [2, 3, 7, 13], [2, 3,
            11, 13], [2, 4, 9, 13], [2, 4, 11, 13], [2, 4, 11, 15], [2, 5, 8,
            11], [2, 5, 8, 12], [2, 5, 10, 12], [2, 6, 10, 12], [2, 6, 10, 14],
            [2, 6, 12, 15], [2, 7, 9, 13], [2, 7, 9, 14], [2, 7, 10, 14], [2,
            8, 11, 15], [2, 8, 12, 15], [3, 4, 5, 14], [3, 4, 5, 15], [3, 4,
            12, 14], [3, 5, 10, 15], [3, 5, 11, 14], [3, 7, 12, 13], [3, 11,
            13, 14], [3, 12, 13, 14], [4, 5, 6, 7], [4, 5, 6, 14], [4, 5, 7,
            15], [4, 6, 7, 11], [4, 6, 10, 11], [4, 6, 10, 14], [4, 7, 11, 15],
            [4, 8, 9, 12], [4, 8, 9, 13], [4, 8, 10, 13], [4, 8, 10, 14], [4,
            8, 12, 14], [4, 10, 11, 13], [5, 6, 7, 13], [5, 7, 9, 13], [5, 7,
            9, 15], [5, 8, 9, 12], [5, 8, 9, 13], [5, 9, 10, 12], [5, 9, 10,
            15], [6, 7, 11, 12], [6, 7, 12, 13], [6, 10, 11, 12], [6, 12, 13,
            15], [7, 8, 10, 14], [7, 8, 11, 15], [7, 8, 14, 15], [7, 9, 14,
            15], [8, 12, 14, 15], [9, 10, 11, 12], [9, 10, 11, 16], [9, 10, 15,
            16], [9, 11, 14, 16], [9, 14, 15, 16], [10, 11, 13, 16], [10, 13,
            15, 16], [11, 13, 14, 16], [12, 13, 14, 15], [13, 14, 15, 16]]);
    assert(poincare.isOrientable);

    /* This should be an 5-sphere, hence orientable. See:
    http://page.math.tu-berlin.de/~lutz/stellar/5_manifolds
    http://page.math.tu-berlin.de/~lutz/stellar/5_manifolds.type */
    auto manifold_5_10_3_1 = Manifold!5([[1, 2, 3, 4, 5, 6], [1, 2, 3, 4,
            5, 10], [1, 2, 3, 4, 6, 7], [1, 2, 3, 4, 7, 8], [1, 2, 3, 4, 8, 9],
            [1, 2, 3, 4, 9, 10], [1, 2, 3, 5, 6, 10], [1, 2, 3, 6, 7, 10], [1,
            2, 3, 7, 8, 10], [1, 2, 3, 8, 9, 10], [1, 2, 4, 5, 6, 7], [1, 2, 4,
            5, 7, 8], [1, 2, 4, 5, 8, 9], [1, 2, 4, 5, 9, 10], [1, 2, 5, 6, 7,
            8], [1, 2, 5, 6, 8, 9], [1, 2, 5, 6, 9, 10], [1, 2, 6, 7, 8, 9],
            [1, 2, 6, 7, 9, 10], [1, 2, 7, 8, 9, 10], [1, 3, 4, 5, 6, 10], [1,
            3, 4, 6, 7, 10], [1, 3, 4, 7, 8, 10], [1, 3, 4, 8, 9, 10], [1, 4,
            5, 6, 7, 10], [1, 4, 5, 7, 8, 10], [1, 4, 5, 8, 9, 10], [1, 5, 6,
            7, 8, 10], [1, 5, 6, 8, 9, 10], [1, 6, 7, 8, 9, 10], [2, 3, 4, 5,
            6, 7], [2, 3, 4, 5, 7, 8], [2, 3, 4, 5, 8, 9], [2, 3, 4, 5, 9, 10],
            [2, 3, 5, 6, 7, 8], [2, 3, 5, 6, 8, 9], [2, 3, 5, 6, 9, 10], [2, 3,
            6, 7, 8, 9], [2, 3, 6, 7, 9, 10], [2, 3, 7, 8, 9, 10], [3, 4, 5, 6,
            7, 8], [3, 4, 5, 6, 8, 9], [3, 4, 5, 6, 9, 10], [3, 4, 6, 7, 8, 9],
            [3, 4, 6, 7, 9, 10], [3, 4, 7, 8, 9, 10], [4, 5, 6, 7, 8, 9], [4,
            5, 6, 7, 9, 10], [4, 5, 7, 8, 9, 10], [5, 6, 7, 8, 9, 10]]);
    assert(manifold_5_10_3_1.isOrientable);

    // http://page.math.tu-berlin.de/~lutz/stellar/manifolds_lex/manifolds_lex_d2_n9_o0_g1
    // surface #9 on the genus 1 non-orientable list
    auto g1 = Manifold!2([[1, 2, 3], [1, 2, 4], [1, 3, 4], [2, 3, 5], [2,
            4, 5], [3, 4, 6], [3, 5, 7], [3, 6, 8], [3, 7, 8], [4, 5, 8], [4,
            6, 9], [4, 7, 8], [4, 7, 9], [5, 6, 8], [5, 6, 9], [5, 7, 9]]);
    assert(!g1.isOrientable);

    // http://page.math.tu-berlin.de/~lutz/stellar/RP3
    auto rp3 = Manifold!3([[1, 2, 3, 7], [1, 2, 3, 11], [1, 2, 6, 9], [1,
            2, 6, 11], [1, 2, 7, 9], [1, 3, 5, 10], [1, 3, 5, 11], [1, 3, 7,
            10], [1, 4, 7, 9], [1, 4, 7, 10], [1, 4, 8, 9], [1, 4, 8, 10], [1,
            5, 6, 8], [1, 5, 6, 11], [1, 5, 8, 10], [1, 6, 8, 9], [2, 3, 4, 8],
            [2, 3, 4, 11], [2, 3, 7, 8], [2, 4, 6, 10], [2, 4, 6, 11], [2, 4,
            8, 10], [2, 5, 7, 8], [2, 5, 7, 9], [2, 5, 8, 10], [2, 5, 9, 10],
            [2, 6, 9, 10], [3, 4, 5, 9], [3, 4, 5, 11], [3, 4, 8, 9], [3, 5, 9,
            10], [3, 6, 7, 8], [3, 6, 7, 10], [3, 6, 8, 9], [3, 6, 9, 10], [4,
            5, 6, 7], [4, 5, 6, 11], [4, 5, 7, 9], [4, 6, 7, 10], [5, 6, 7, 8]]);
    assert(rp3.isOrientable);

    auto rp4 = Manifold!4([[1, 2, 4, 5, 11], [1, 2, 4, 5, 14], [1, 2, 4, 11,
        13], [1, 2, 4, 13, 14], [1, 2, 5, 11, 15], [1, 2, 5, 14, 15], [1, 2,
        7, 8, 13], [1, 2, 7, 8, 15], [1, 2, 7, 13, 14], [1, 2, 7, 14, 15], [1,
        2, 8, 11, 13], [1, 2, 8, 11, 15], [1, 3, 4, 10, 13], [1, 3, 4, 10, 16],
        [1, 3, 4, 11, 13], [1, 3, 4, 11, 16], [1, 3, 8, 9, 11], [1, 3, 8, 9,
        12], [1, 3, 8, 11, 13], [1, 3, 8, 12, 13], [1, 3, 9, 11, 16], [1, 3, 9,
        12, 16], [1, 3, 10, 12, 13], [1, 3, 10, 12, 16], [1, 4, 5, 11, 16], [1,
        4, 5, 14, 16], [1, 4, 10, 13, 14], [1, 4, 10, 14, 16], [1, 5, 6, 9,
        15], [1, 5, 6, 9, 16], [1, 5, 6, 14, 15], [1, 5, 6, 14, 16], [1, 5,
        9, 11, 15], [1, 5, 9, 11, 16], [1, 6, 7, 10, 12], [1, 6, 7, 10, 14],
        [1, 6, 7, 12, 15], [1, 6, 7, 14, 15], [1, 6, 9, 12, 15], [1, 6, 9, 12,
        16], [1, 6, 10, 12, 16], [1, 6, 10, 14, 16], [1, 7, 8, 12, 13], [1, 7,
        8, 12, 15], [1, 7, 10, 12, 13], [1, 7, 10, 13, 14], [1, 8, 9, 11, 15],
        [1, 8, 9, 12, 15], [2, 3, 5, 10, 12], [2, 3, 5, 10, 15], [2, 3, 5, 12,
        14], [2, 3, 5, 14, 15], [2, 3, 7, 9, 14], [2, 3, 7, 9, 16], [2, 3, 7,
        14, 15], [2, 3, 7, 15, 16], [2, 3, 9, 12, 14], [2, 3, 9, 12, 16], [2,
        3, 10, 12, 16], [2, 3, 10, 15, 16], [2, 4, 5, 11, 12], [2, 4, 5, 12,
        14], [2, 4, 6, 9, 12], [2, 4, 6, 9, 13], [2, 4, 6, 11, 12], [2, 4, 6,
        11, 13], [2, 4, 9, 12, 14], [2, 4, 9, 13, 14], [2, 5, 10, 11, 12], [2,
        5, 10, 11, 15], [2, 6, 8, 10, 11], [2, 6, 8, 10, 16], [2, 6, 8, 11,
        13], [2, 6, 8, 13, 16], [2, 6, 9, 12, 16], [2, 6, 9, 13, 16], [2, 6,
        10, 11, 12], [2, 6, 10, 12, 16], [2, 7, 8, 13, 16], [2, 7, 8, 15, 16],
        [2, 7, 9, 13, 14], [2, 7, 9, 13, 16], [2, 8, 10, 11, 15], [2, 8, 10,
        15, 16], [3, 4, 6, 7, 11], [3, 4, 6, 7, 15], [3, 4, 6, 11, 13], [3, 4,
        6, 13, 15], [3, 4, 7, 11, 16], [3, 4, 7, 15, 16], [3, 4, 10, 13, 15],
        [3, 4, 10, 15, 16], [3, 5, 6, 8, 13], [3, 5, 6, 8, 14], [3, 5, 6, 13,
        15], [3, 5, 6, 14, 15], [3, 5, 8, 12, 13], [3, 5, 8, 12, 14], [3, 5,
        10, 12, 13], [3, 5, 10, 13, 15], [3, 6, 7, 11, 14], [3, 6, 7, 14, 15],
        [3, 6, 8, 11, 13], [3, 6, 8, 11, 14], [3, 7, 9, 11, 14], [3, 7, 9, 11,
        16], [3, 8, 9, 11, 14], [3, 8, 9, 12, 14], [4, 5, 7, 8, 12], [4, 5, 7,
        8, 16], [4, 5, 7, 11, 12], [4, 5, 7, 11, 16], [4, 5, 8, 12, 14], [4, 5,
        8, 14, 16], [4, 6, 7, 11, 12], [4, 6, 7, 12, 15], [4, 6, 9, 12, 15],
        [4, 6, 9, 13, 15], [4, 7, 8, 12, 15], [4, 7, 8, 15, 16], [4, 8, 9, 10,
        14], [4, 8, 9, 10, 15], [4, 8, 9, 12, 14], [4, 8, 9, 12, 15], [4, 8,
        10, 14, 16], [4, 8, 10, 15, 16], [4, 9, 10, 13, 14], [4, 9, 10, 13,
        15], [5, 6, 8, 13, 16], [5, 6, 8, 14, 16], [5, 6, 9, 13, 15], [5, 6,
        9, 13, 16], [5, 7, 8, 12, 13], [5, 7, 8, 13, 16], [5, 7, 9, 10, 11],
        [5, 7, 9, 10, 13], [5, 7, 9, 11, 16], [5, 7, 9, 13, 16], [5, 7, 10, 11,
        12], [5, 7, 10, 12, 13], [5, 9, 10, 11, 15], [5, 9, 10, 13, 15], [6, 7,
        10, 11, 12], [6, 7, 10, 11, 14], [6, 8, 10, 11, 14], [6, 8, 10, 14,
        16], [7, 9, 10, 11, 14], [7, 9, 10, 13, 14], [8, 9, 10, 11, 14], [8, 9, 
        10, 11, 15]]);
    assert(!rp4.isOrientable);
}

/*******************************************************************************
Returns a range containing this simplicial complex's connected components
(returned as simplicial complexes of the same vertex type.)
*/
auto connectedComponents(Vertex)(const ref SimplicialComplex!Vertex sc)
{
    static struct FacetRecord
    {
        const(Vertex)[] facet;
        int label;      // labels which connected component facet is in
        bool seenNear;  // looked at facets in the star of the vertices?
    }

    /* We will label each facet with an integer 1, 2, 3 ... etc, marking its
    connected component */
    // auto labels = new int[numFacets];
    // auto allFacets = facets.array;

    /* We keep track of the facets we've yet to examine, the facets we've 
    labeled (but are not yet done with) and those that are completed (which
    means both labeled and all the facets in the star of its vertices are
    also labeled. */
    FacetRecord[] records = sc.facets.map!(f => FacetRecord(f)).array;

    int currentLabel = 1;

    // Deal with the facets that are (necessarily isolated) vertices
    auto nVertFac = sc.facets(0).walkLength;
    auto vertexFacets = records[0 .. nVertFac];

    foreach(indx, ref f; vertexFacets)
    {
        f.label = currentLabel;
        f.seenNear = true;
        ++currentLabel;
    }

    while(!records.find!(r => r.label == 0).empty)
    {
        records.find!(r => r.label == 0).front.label = currentLabel;

        // Now we propagate the current label as far as we can
        while(records.canFind!(r => r.label != 0 && !r.seenNear))
        {
            auto toDo = records.find!(r => r.label != 0 && !r.seenNear);

            foreach(f; toDo.front.facet.map!(v => sc.star([v])).map!array.joiner)
            {
                auto findIt = records.find!(r => r.facet == f);
                assert(!findIt.empty);
                findIt.front.label = currentLabel;
            }

            toDo.front.seenNear = true;           
        }
        currentLabel += 1;
    }

    return records.chunkBy!((r1, r2) => r1.label == r2.label)
        .map!(rList => SimplicialComplex!Vertex(rList.map!(r => r.facet).array));           
}

///
@Name("connectedComponents") /* pure */ @system unittest
{
    // TO DO: Why does this need to be @system? Make it @safe!
    // TO DO: ldc doesn't like using "pure" above! Bugreport?

    alias sComp = simplicialComplex;

    auto sc = sComp([[1], [9], [2,3], [3,4], [5,6], [6,7,8]]);

    auto c1 = sComp([[1]]);
    auto c2 = sComp([[9]]);
    auto c3 = sComp([[2,3], [3,4]]);
    auto c4 = sComp([[5,6], [6,7,8]]);

    assert(sc.connectedComponents.array == [c1, c2, c3, c4]);

    auto emptyComplex = SimplicialComplex!int();
    assert(emptyComplex.connectedComponents.empty); 
}

/*******************************************************************************
Returns the Euler characteristic of the simplicial complex
*/
int eulerCharacteristic(Vertex)(const ref SimplicialComplex!Vertex sc)
{
    return sc.fVector.enumerate.map!(f => (-1)^^f.index.to!int * f.value).sum;
}
///
@Name("eulerCharacteristic") pure @safe unittest
{
    SimplicialComplex!() sc;
    sc.insertFacet([1,2]);
    assert(sc.eulerCharacteristic == 1);
    sc.insertFacet([2,3]);
    sc.insertFacet([2,4]);
    assert(sc.eulerCharacteristic == 1);
    sc.insertFacet([3,4]);
    assert(sc.eulerCharacteristic == 0);
    sc.insertFacet([1,3]);
    assert(sc.eulerCharacteristic == -1);        
}

/*******************************************************************************
Decide if a simplicial complex is homeomorphic to a 2-sphere
*/
bool is2Sphere(Vertex)(SimplicialComplex!Vertex sc)
{

    return sc.isOrientableSurfaceOfGenus(0);
}
///
@Name("is2Sphere") /* pure */ @system unittest
{
    // TO DO: Why does this need to be @system? Make it @safe!
    // TO DO: ldc doesn't like using "pure" above! Bugreport?

    auto s1 = simplicialComplex([[1,2,3], [1,2,4], [1,3,4], [2,3,4]]);
    assert(s1.is2Sphere);

    auto disjoint2Spheres = SimplicialComplex!()(s1.facets.array
        ~ [[5,6,7], [5,6,8], [5,7,8], [6,7,8]].to!(const(int)[][]));
    assert(!disjoint2Spheres.is2Sphere);

    auto s2 = simplicialComplex(s1.facets.array
        ~ [[1,6,7], [1,6,8], [1,7,8], [6,7,8]].to!(const(int)[][]));
    assert(!s2.is2Sphere);

    // http://page.math.tu-berlin.de/~lutz/stellar/manifolds_lex/manifolds_lex_d2_n10_o1_g0
    // Surface #168 in the genus 0 (sphere) list
    auto g0 = simplicialComplex([[1,2,3],[1,2,4],[1,3,4],[2,3,5],[2,4,6],
        [2,5,6],[3,4,7],[3,5,6],[3,6,8],[3,7,9],[3,8,9],[4,6,10],[4,7,9],
        [4,9,10],[6,8,10],[8,9,10]]);
    assert(g0.is2Sphere);

    auto emptyComplex = SimplicialComplex!()();
    assert(!emptyComplex.is2Sphere);
}

/*******************************************************************************
Decide if a simplicial complex is homeomorphic to a 2-torus
*/
bool is2Torus(Vertex)(SimplicialComplex!Vertex sc)
{
    return sc.isOrientableSurfaceOfGenus(1);
}
///
@Name("is2Torus") /* pure */ @system unittest
{
    // TO DO: Why does this need to be @system? Make it @safe!
    // TO DO: ldc doesn't like using "pure" above! Bugreport?

    // http://page.math.tu-berlin.de/~lutz/stellar/manifolds_lex/manifolds_lex_d2_n10_o1_g1
    // Surface #2105 on genus 1 list
    auto g1 = simplicialComplex([[1,2,3],[1,2,4],[1,3,5],[1,4,6],[1,5,6],
        [2,3,7],[2,4,8],[2,5,8],[2,5,9],[2,7,9],[3,5,10],[3,7,10],[4,6,9],
        [4,8,10],[4,9,10],[5,6,8],[5,9,10],[6,7,8],[6,7,9],[7,8,10]]);
    assert(g1.is2Torus);

    auto emptyComplex = SimplicialComplex!()();
    assert(!emptyComplex.is2Torus);
}

/*******************************************************************************
Decide if a simplicial complex is homeomorphic to a 1-sphere (circle)
*/
bool isCircle(Vertex)(const SimplicialComplex!Vertex sc)
{
    return !sc.facets.empty
        && sc.facets.all!(f => f.walkLength == 2)
        && sc.simplices(0).all!(v => sc.star(v).walkLength == 2)
        && sc.isConnected;
}
///
@Name("isCircle") /* pure */ @system unittest
{
    // TO DO: Why does this need to be @system? Make it @safe!
    // TO DO: ldc doesn't like using "pure" above! Bugreport?

    // Start with a circle with 4 edges
    auto s = simplicialComplex([[1,2], [2,3], [3,4], [1,4]]);
    assert(s.isCircle);

    s.removeFacet([2,3]);
    assert(!s.isCircle);

    // Simplicial complex must be homeomorphic to a circle, not just homotopic
    s.insertFacet([2,3,5]);
    assert(!s.isCircle);

    // Must be a connected simplicial complex
    auto disjointCircles = simplicialComplex([
        [1,2], [2,3], [1,3],
        [4,5], [5,6], [4,6]]);

    assert(!disjointCircles.isCircle);

    auto emptyComplex = SimplicialComplex!()();
    assert(!emptyComplex.isCircle);
}

/***************************************************************************
Returns true if the simplicial complex is connected and false otherwise.
Note that an empty complex counts as connected.
*/
bool isConnected(Vertex)(const ref SimplicialComplex!Vertex sc)
{
    return sc.connectedComponents.walkLength <= 1;
}
///
@Name("isConnected") /* pure */ @system unittest
{
    // TO DO: Why does this need to be @system? Make it @safe!
    // TO DO: ldc doesn't like using "pure" above! Bugreport?

    auto disjointCircles = simplicialComplex([
        [1,2], [2,3], [1,3],
        [4,5], [5,6], [4,6]]);
        
    assert(!disjointCircles.isConnected);

    auto barbell = SimplicialComplex!()(disjointCircles.facets.array ~ [[3,4]].to!(const(int)[][]));
    assert(barbell.isConnected);

    auto emptyComplex = SimplicialComplex!()();
    assert(emptyComplex.isConnected);
}

/*******************************************************************************
Decide if a simplicial complex is pure of dimension `d`
*/
bool isPureOfDim(Vertex)(const ref SimplicialComplex!Vertex sc, int d)
{
    assert(d >= 0, "expected a non-negative dimension");
    return sc.facets(d).walkLength == sc.numFacets;
}
///
@Name("isPureOfDim") pure @system unittest
{
    // TO DO: Why does this need to be @system? Make it @safe!

    auto sc = SimplicialComplex!()();

    // An empty simplicial complex is pure of any dimension
    assert(iota(10).all!(d => sc.isPureOfDim(d)));

    sc.insertFacet([1,2]);
    sc.insertFacet([2,3]);
    assert(sc.isPureOfDim(1));
    assert(!sc.isPureOfDim(0));
    assert(iota(2, 10).all!(d => !sc.isPureOfDim(d)));
    
    sc.insertFacet([4,5,6]);
    assert(iota(10).all!(d => !sc.isPureOfDim(d)));
    
    throwsWithMsg!Error(sc.isPureOfDim(-2),
        "expected a non-negative dimension");
    
    auto emptyComplex = SimplicialComplex!()();
    assert(iota(16).all!(d => emptyComplex.isPureOfDim(d)));
}

/*******************************************************************************
Decide if a simplicial complex is homeomorphic to a surface of genus `g`
*/
bool isOrientableSurfaceOfGenus(Vertex)(const ref SimplicialComplex!Vertex sc, int g)
{
    alias SimpComp = SimplicialComplex!Vertex;

    // TO DO: Find a better way to do orientability check
    return !sc.facets.empty
        && sc.isConnected
        && sc.isPureOfDim(2)
        && sc.simplices(0).all!(v => sc.link(v).to!SimpComp.isCircle)
        && sc.eulerCharacteristic == 2 - 2 * g
        && Manifold!2(sc.facets).isOrientable;
}
///
@Name("isSurfaceOfGenus") /* pure */ @system unittest
{
    // TO DO: Why does this need to be @system? Make it @safe!
    // TO DO: ldc doesn't like using "pure" above! Bugreport?

    // http://page.math.tu-berlin.de/~lutz/stellar/manifolds_lex/manifolds_lex_d2_n10_o1_g2
    // Surface #514 in the genus 2 list 
    auto g2 = simplicialComplex([[1,2,3],[1,2,4],[1,3,5],[1,4,6],[1,5,6],
        [2,3,6],[2,4,7],[2,6,8],[2,7,9],[2,8,9],[3,5,8],[3,6,9],[3,7,8],
        [3,7,9],[4,5,7],[4,5,9],[4,6,10],[4,9,10],[5,6,9],[5,7,10],[5,8,10],
        [6,7,8],[6,7,10],[8,9,10]]);
    assert(g2.isOrientableSurfaceOfGenus(2));

    // http://page.math.tu-berlin.de/~lutz/stellar/manifolds_lex/manifolds_lex_d2_n12_o1_g6
    // Surface #30 in the genus 6 list
    auto g6 = simplicialComplex([[1,2,3],[1,2,4],[1,3,5],[1,4,6],[1,5,7],
        [1,6,8],[1,7,9],[1,8,10],[1,9,11],[1,10,12],[1,11,12],[2,3,6],[2,4,5],
        [2,5,9],[2,6,11],[2,7,8],[2,7,12],[2,8,11],[2,9,10],[2,10,12],[3,4,7],
        [3,4,9],[3,5,12],[3,6,10],[3,7,11],[3,8,11],[3,8,12],[3,9,10],[4,5,8],
        [4,6,12],[4,7,10],[4,8,9],[4,10,11],[4,11,12],[5,6,9],[5,6,12],[5,7,11],
        [5,8,10],[5,10,11],[6,7,8],[6,7,10],[6,9,11],[7,9,12],[8,9,12]]);
    assert(g6.isOrientableSurfaceOfGenus(6));

    auto emptyComplex = SimplicialComplex!()();
    assert(iota(0, 10).all!(g => !emptyComplex.isOrientableSurfaceOfGenus(g)));
}

/***************************************************************************
Returns the join of two simplicial complexes
*/
auto join(Vertex)(const SimplicialComplex!Vertex sc1, 
                  const SimplicialComplex!Vertex sc2)
{
    const(Vertex)[][] result;
    foreach(f1; sc1.facets)
    {
        foreach(f2; sc2.facets)
        {
            assert(setIntersection(f1, f2).empty);
            result ~= (f1 ~ f2).array.dup.sort().array;
        }
    }

    return SimplicialComplex!Vertex(result);
}
///
@Name("join") pure @safe unittest
{
    auto sc1 = simplicialComplex([[1,2], [2,3,4], [5]]);
    auto sc2 = simplicialComplex([[6,7], [8]]);

    assert(join(sc1, sc2).facets.equal([
        [5, 8], [1, 2, 8], [5, 6, 7],
        [1, 2, 6, 7], [2, 3, 4, 8], [2, 3, 4, 6, 7]]));

    auto emptyComplex = SimplicialComplex!()();
    assert(join(sc1, emptyComplex).facets.empty);
}
