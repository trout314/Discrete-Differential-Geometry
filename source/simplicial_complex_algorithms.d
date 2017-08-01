import unit_threaded : Name;
import simplicial_complex : SimplicialComplex, simplicialComplex, fVector;

import std.algorithm : all, canFind, chunkBy, find, joiner, map,
    setIntersection, sort, sum;
import std.exception : enforce;
import std.conv : to;
import std.range : array, empty, front, enumerate, iota, walkLength;
import utility : throwsWithMsg;

import simplex : simplex;

/*******************************************************************************
Returns a range containing this simplicial complex's connected components
(returned as simplicial complexes of the same vertex type.)
*/
auto connectedComponents(Vertex)(const ref SimplicialComplex!Vertex sc)
{
    static struct FacetRecord
    {
        Vertex[] facet;
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
    auto nVertFac = sc.facets!0.walkLength;
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

            foreach(f; toDo.front.facet.map!(v => sc.star(simplex(v))).joiner)
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
@Name("connectedComponents")
unittest
{
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
@Name("eulerCharacteristic")
unittest
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

    return sc.isSurfaceOfGenus(0);
}
///
@Name("is2Sphere")
unittest
{
    auto s1 = simplicialComplex([[1,2,3], [1,2,4], [1,3,4], [2,3,4]]);
    assert(s1.is2Sphere);

    auto disjoint2Spheres = SimplicialComplex!()(s1.facets
        ~ [[5,6,7], [5,6,8], [5,7,8], [6,7,8]]);
    assert(!disjoint2Spheres.is2Sphere);

    auto s2 = simplicialComplex(s1.facets
        ~ [[1,6,7], [1,6,8], [1,7,8], [6,7,8]]);
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
    return sc.isSurfaceOfGenus(1);
}
///
@Name("is2Torus")
unittest
{
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
        && sc.simplices!0.all!(v => sc.star(v).walkLength == 2)
        && sc.isConnected;
}
///
@Name("isCircle")
unittest
{
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
@Name("isConnected")
unittest
{
    auto disjointCircles = simplicialComplex([
        [1,2], [2,3], [1,3],
        [4,5], [5,6], [4,6]]);
        
    assert(!disjointCircles.isConnected);

    auto barbell = SimplicialComplex!()(disjointCircles.facets ~ [[3,4]]);
    assert(barbell.isConnected);

    auto emptyComplex = SimplicialComplex!()();
    assert(emptyComplex.isConnected);
}

/*******************************************************************************
Decide if a simplicial complex is pure of dimension `d`
*/
bool isPureOfDim(Vertex)(const ref SimplicialComplex!Vertex sc, int d)
{
    enforce(d >= 0, "expected a non-negative dimension but got " ~ d.to!string);
    return sc.facets(d).walkLength == sc.numFacets;
}
///
@Name("isPureOfDim")
unittest
{
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
    
    throwsWithMsg(sc.isPureOfDim(-2),
        "expected a non-negative dimension but got -2");
    
    auto emptyComplex = SimplicialComplex!()();
    assert(iota(16).all!(d => emptyComplex.isPureOfDim(d)));
}

/*******************************************************************************
Decide if a simplicial complex is homeomorphic to a surface of genus `g`
*/
bool isSurfaceOfGenus(Vertex)(const ref SimplicialComplex!Vertex sc, int g)
{
    alias SimpComp = SimplicialComplex!Vertex;

    return !sc.facets.empty
        && sc.isConnected
        && sc.isPureOfDim(2)
        && sc.simplices!0.all!(v => sc.link(v).to!SimpComp.isCircle)
        && sc.eulerCharacteristic == 2 - 2 * g;    
}
///
@Name("isSurfaceOfGenus")
unittest
{
    // http://page.math.tu-berlin.de/~lutz/stellar/manifolds_lex/manifolds_lex_d2_n10_o1_g2
    // Surface #514 in the genus 2 list 
    auto g2 = simplicialComplex([[1,2,3],[1,2,4],[1,3,5],[1,4,6],[1,5,6],
        [2,3,6],[2,4,7],[2,6,8],[2,7,9],[2,8,9],[3,5,8],[3,6,9],[3,7,8],
        [3,7,9],[4,5,7],[4,5,9],[4,6,10],[4,9,10],[5,6,9],[5,7,10],[5,8,10],
        [6,7,8],[6,7,10],[8,9,10]]);
    assert(g2.isSurfaceOfGenus(2));

    // http://page.math.tu-berlin.de/~lutz/stellar/manifolds_lex/manifolds_lex_d2_n12_o1_g6
    // Surface #30 in the genus 6 list
    auto g6 = simplicialComplex([[1,2,3],[1,2,4],[1,3,5],[1,4,6],[1,5,7],
        [1,6,8],[1,7,9],[1,8,10],[1,9,11],[1,10,12],[1,11,12],[2,3,6],[2,4,5],
        [2,5,9],[2,6,11],[2,7,8],[2,7,12],[2,8,11],[2,9,10],[2,10,12],[3,4,7],
        [3,4,9],[3,5,12],[3,6,10],[3,7,11],[3,8,11],[3,8,12],[3,9,10],[4,5,8],
        [4,6,12],[4,7,10],[4,8,9],[4,10,11],[4,11,12],[5,6,9],[5,6,12],[5,7,11],
        [5,8,10],[5,10,11],[6,7,8],[6,7,10],[6,9,11],[7,9,12],[8,9,12]]);
    assert(g6.isSurfaceOfGenus(6));

    auto emptyComplex = SimplicialComplex!()();
    assert(iota(0, 10).all!(g => !emptyComplex.isSurfaceOfGenus(g)));
}

/***************************************************************************
Returns the join of two simplicial complexes
*/
auto join(Vertex)(const SimplicialComplex!Vertex sc1, 
                  const SimplicialComplex!Vertex sc2)
{
    Vertex[][] result;
    foreach(f1; sc1.facets)
    {
        foreach(f2; sc2.facets)
        {
            assert(setIntersection(f1, f2).empty);
            result ~= (f1 ~ f2).sort().array;
        }
    }

    return SimplicialComplex!Vertex(result);
}
///
@Name("join")
unittest
{
    auto sc1 = simplicialComplex([[1,2], [2,3,4], [5]]);
    auto sc2 = simplicialComplex([[6,7], [8]]);

    assert(join(sc1, sc2).facets == [
        [5, 8], [1, 2, 8], [5, 6, 7],
        [1, 2, 6, 7], [2, 3, 4, 8], [2, 3, 4, 6, 7]]);

    auto emptyComplex = SimplicialComplex!()();
    assert(join(sc1, emptyComplex).facets.empty);
}
