import std.algorithm : all, canFind, countUntil, equal, filter, findAdjacent,
    remove, sort;
import std.array : array;
import std.conv : parse, to;
import std.datetime : Clock, DateTime;
import std.exception : assertThrown;
import std.getopt : getopt;
import std.range : back, front, iota, popBack;
import std.stdio : File, lines, writeln;
import std.string : strip;

/******************************************************************************/

/// Type used for vertex labels. Required to be an unsigned integer type.
alias Vertex = size_t;

/// Type for storing index of a facet record. Must be an unsigned integer type.
alias FacetIndex = size_t;

/// Type for storing hinge degrees. Must be an unsigned integer type.
alias Degree = size_t;

/// Type representing simplices during processing.
alias Simplex = size_t[];

/***
* Dimension for manifolds used in the program. Must be at least 2, and at most
* maxDimesnion.
*/
enum dimManifold = 3;

///
enum maxDimension = 16; // Can go higher, but unit tests will take longer
enum minDimension = 2; // Easier to not support (boring!) 1-manifolds
static assert(dimManifold <= maxDimension);
static assert(dimManifold >= minDimension);

/******************************************************************************
* Data structure for the top-dimensional simplices (called facets) in the
* manifold. Three things make up a FacetRecord. We of course must store the
* the vertices of the facet. In addition, to avoid iterating over every facet
* in the manifold, we store indices (in the list of all facets in the manifold)
* for the "neighboring" facets. These are the facets in the manifold which
* share a codimension-1 face (called a ridge) with this simplex. The neighbors
* are listed so that the i-th one shares the ridge missing the i-th vertex.
* Finally, for quick access we also store the degree of each codimension-2
* face (called a hinge) in this facet. In accordance with our ordering rule
* for neighbors, the degrees of these hinges are listed so that the indices
* the two vertices not part of the hinge appear in dictionary order.
*/
struct FacetRecord(size_t dimension = dimManifold)
{
    static assert(dimension >= 2);
    static assert(dimension <= maxDimension);
    ///
    Vertex[dimension + 1] facet;
    ///
    FacetIndex[dimension + 1] neighbors;
    ///
    Degree[dimension * (dimension + 1) / 2] degrees;
}

/******************************************************************************
* Data structure for holding an array of facet records and updating them in
* a coherent way. A manifold of the given dimension may be constructed from
* an array of facets, each given as an array of vertex labels.
*/
struct Manifold(size_t dim = dimManifold)
{
    static assert(dim >= 2);
    static assert(dim <= maxDimension);

    enum dimension = dim;
    FacetRecord!dimension[] facetRecords;
    Vertex[] freeVertexLabels;

    this(const size_t[][] facets)
    {
        foreach (const facet; facets)
        {
            assert(facet.length == dimension + 1,
                    "encountered facet with wrong number of vertices in construction of manifold");

            assert(!facet.isDegenerate,
                    "encountered degenerate facet (with repeated vertex labels) in construction of manifold");

            FacetRecord!dimension newFacetRecord;
            newFacetRecord.facet = facet;
            facetRecords ~= newFacetRecord;

            // TO DO: Check that vertex labels are exactly 0,1, .. #vertices
            // if they are not, relabel vertices to make this true
        }

        // TO DO: When dimension = 2, 3 check for manifold-ness
        // TO DO: (TON OF WORK!) Check manifold-ness for n = 4

        facetRecords.setAllNeighbors;
        assert(facetRecords.checkNeighbors2,
                "failed to correctly set neighbor data in construction of manifold");

        facetRecords.setAllDegrees;
        assert(facetRecords.checkDegrees,
                "failed to correctly set degree data in construction of manifold");
    }
}
///
unittest
{

    // The octahedral surface with eight triangular facets
    auto octahedron = Manifold!2([[0, 1, 2], [0, 2, 3], [0, 3, 4], [0, 4, 1],
            [5, 1, 2], [5, 2, 3], [5, 3, 4], [5, 4, 1]]);
    assert(octahedron.facetRecords[].all!(r => r.degrees[].all!(d => d == 4)));

    // The standard 4-sphere (given by boundary of 5-simplex.)
    auto sphere4 = Manifold!4(standardSphere(4));
    import std.algorithm : all;

    assert(sphere4.facetRecords[].all!(r => r.degrees[].all!(d => d == 3)));

    assertThrown!Error(Manifold!2([[1, 2, 3], [5, 5, 2]]));
    assertThrown!Error(Manifold!2([[1, 2, 3], [1, 2, 3, 5]]));
}

/******************************************************************************
* Prints in a nice way the list of facet records of a given manifold.
*/
void print(size_t dimension)(ref const Manifold!dimension manifold)
{
    foreach (index, record; manifold.facetRecords)
    {
        writeln("#", index, " facet ", record.facet, " neighbors ",
                record.neighbors, " degrees ", record.degrees);
    }
}

/******************************************************************************
* Returns true if and only if possibleFace is a face of simplex. That is, the
* function returns true if and only if each vertex label of possibleFace
* appears in simplex.
*/
bool hasFace(const Simplex simplex, const Simplex possibleFace)
{
    return possibleFace.all!(vertex => simplex.canFind(vertex));
}
///
unittest
{

    assert(hasFace([3, 2, 1], [2, 1, 3]));
    assert(hasFace([1, 3, 2, 5], [2, 1]));
    assert([1, 3, 2].hasFace([3]));
    assert(hasFace([1, 2, 3], []));
    assert(hasFace([], []));
    assert(![1, 3, 2].hasFace([1, 4]));
    assert(!hasFace([1, 3, 2], [5]));
    assert(!hasFace([], [5]));
}

/******************************************************************************
* Returns an array of codimension-1 faces (ridges) of the given simplex.
* The ridges are ordered so that the index of the missing vertex in a ridge
* goes in increasing order.
*/
Simplex[] ridges(in Simplex simplex)
{
    if (simplex.length <= 1)
        return []; // No ridges in dimension zero

    Simplex[] output;
    foreach (index; 0 .. simplex.length)
    {
        output ~= simplex.dup.remove(index);
    }
    return output;
}
///
unittest
{

    Simplex vertex = [3];
    assert(vertex.ridges == []);

    Simplex edge = [1, 8];
    assert(edge.ridges == [[8], [1]]);

    assert([5, 2, 3].ridges == [[2, 3], [5, 3], [5, 2]]);
    assert([1, 2, 3, 4].ridges == [[2, 3, 4], [1, 3, 4], [1, 2, 4], [1, 2, 3]]);
}

/******************************************************************************
* Returns list of the hinges that are faces of a given simplex. The hinges are
* ordered so that the indices of the two vertices not part of the hinge appear
* in dictionary order. Note that this is the same order as the hinge-degrees
* stored in a facetRecord.
*/
Simplex[] hinges(in Simplex simplex)
{
    if (simplex.length <= 2)
        return []; // No hinges if dimensions <= 1

    Simplex[] hinges;
    foreach (firstIndex; 0 .. simplex.length)
    {
        foreach (secondIndex; firstIndex .. simplex.length - 1)
        {
            hinges ~= simplex.dup.remove(firstIndex).remove(secondIndex);
        }
    }
    return hinges;
}
///
unittest
{

    assert([2].hinges == []);
    assert([1, 2].hinges == []);

    Simplex triangle = [4, 1, 2];
    assert(triangle.hinges == [[2], [1], [4]]);

    Simplex tetrahedron = [1, 2, 3, 4];
    assert(tetrahedron.hinges == [[3, 4], [2, 4], [2, 3], [1, 4], [1, 3], [1, 2]]);
}

/******************************************************************************
* Returns the "standard" simplex with vertex labels [0, 1, 2, ... , dimension]
*/
Simplex standardSimplex(size_t dimension = dimManifold)
{
    // The + 1 below is there to allow for construction of a standard sphere
    // with dimension = maxDimension from the boundary of one-higher dimensional
    // simplex.
    assert(dimension <= maxDimension + 1);
    return iota(dimension + 1).array;
}
///
unittest
{

    assert(standardSimplex(0) == [0]);
    assert(standardSimplex(3) == [0, 1, 2, 3]);
    assertThrown!Error(standardSimplex(maxDimension + 2));
}

/*****************************************************************************
* Returns facets in the "standard" sphere of given dimension. This is the
* boundary of the standard simplex of one higher dimension. The facets of this
* manifold are given by the ridges of this higher dimensional simplex
*/
Simplex[] standardSphere(size_t dimension = dimManifold)
{
    assert(dimension <= maxDimension);
    return standardSimplex(dimension + 1).ridges;
}
///
unittest
{

    assert(standardSphere(0) == [[1], [0]]);
    assert(standardSphere(1) == [[1, 2], [0, 2], [0, 1]]);
    assert(standardSphere(2) == [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]]);
    assert(standardSphere(maxDimension).length == maxDimension + 2);
    assertThrown!Error(standardSphere(maxDimension + 1));
}

/+*****************************************************************************
+ Sets the neighbor indices for a list of facet records (containing correct
+ vertex label information.)  This task is O(M*N^2) wherer N in the number of
+ facets and M is the number of ridges per facet.
+/
void setAllNeighbors(size_t dimension)(FacetRecord!dimension[] facetRecords)
{
    static assert(dimension >= 2);
    static assert(dimension <= maxDimension);

    foreach (thisIndex, ref thisRecord; facetRecords)
    {
        immutable thisFacet = thisRecord.facet;
        foreach (neighborIndex, commonRidge; thisFacet.ridges)
        {
            foreach (neighborRecordIndex, neighborRecord; facetRecords)
            {
                immutable neighborFacet = neighborRecord.facet;
                if (neighborFacet.hasFace(commonRidge) && (thisIndex != neighborRecordIndex))
                {
                    thisRecord.neighbors[neighborIndex] = neighborRecordIndex;
                }
            }
        }
    }
}

/+*****************************************************************************
+ Checks that the neighbor's in a given set of facet records match up as they
+ should. Returns true if they do and false if the don't.
+/
bool checkNeighbors(size_t dim)(const FacetRecord!dim[] facetRecords)
{
    static assert(dim >= 2);
    static assert(dim <= maxDimension);

    /+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   + TO DO: (Bug report?)
   + When the foreach loop below is written: foreach(record, facetRecords) 
   + something rotten happens here! We get instances of incorrectly formed
   + neighbor information. This occurs for dmd 2.066.1 but not for
   + gdc 4.9.1, and not for ldc 0.14.0. In all cases the alternate neighbor-
   + checking function checkNeighbors2 deems the neighbor info correct.
   +/

    // Check that neighboring facets share the proper common ridge
    foreach (record; facetRecords)
    {
        immutable thisFacet = record.facet;
        foreach (index, commonRidge; thisFacet.ridges)
        {
            immutable neighborIndex = record.neighbors[index];
            immutable neighborFacet = facetRecords[neighborIndex].facet;
            if ((!thisFacet.hasFace(commonRidge)) || (!neighborFacet.hasFace(commonRidge)))
            {
                return false;
            }
        }
    }
    return true;
}

/+*****************************************************************************
+ Checks that the neighbor's in a given set of facet records match up as they
+ should. Returns true if they do and false if they don't. This gives an
+ alternate check to checkNeighbors. This was written because of the
+ previous (fixed?) issue with that function. (See TO DO in that function's
+ main loop.)
+/
bool checkNeighbors2(size_t dim)(const FacetRecord!dim[] facetRecords)
{
    static assert(dim >= 2);
    static assert(dim <= maxDimension);

    // Buid an associative array of the indexes of all facet records having as a
    // face each ridge seen in the manifold
    FacetIndex[][immutable Simplex] recordIndicesByRidge;
    foreach (index, record; facetRecords)
    {
        foreach (ridge; record.facet.ridges)
        {
            recordIndicesByRidge[sort(ridge).array.idup] ~= index;
        }
    }

    // Check that things are as they should be. Nothing wrong at start.
    bool correct = true;

    // From double counting argument, should have 2f_{n-1} == (n+1)f_n
    if (2 * recordIndicesByRidge.length != (dim + 1) * facetRecords.length)
    {
        correct = false;
    }

    // Check that the listed facet indices are correct
    foreach (ridge, indices; recordIndicesByRidge)
    {
        assert(ridge.length == dim);
        assert(indices.length == 2);

        immutable firstRecord = facetRecords[indices.front];
        immutable firstFacet = firstRecord.facet;
        immutable secondRecord = facetRecords[indices.back];
        immutable secondFacet = secondRecord.facet;

        // Intersection of first two facets should be exactly the ridge
        if (!commonFace(firstFacet, secondFacet).isEqual(ridge))
        {
            correct = false;
        }

        // Each facet should contain the ridge's vertices plus one more      
        if ((firstFacet.faceOpposite(ridge).length != 1)
                || (secondFacet.faceOpposite(ridge).length != 1))
        {
            correct = false;
        }

        // Check facets under consideration are proper neighbors of eachother
        auto firstMissingVertex = firstFacet.faceOpposite(ridge).front;
        auto firstNeighbor = firstFacet.neighborIndex(firstMissingVertex);
        auto secondMissingVertex = secondFacet.faceOpposite(ridge).front;
        auto secondNeighbor = secondFacet.neighborIndex(secondMissingVertex);

        if ((firstRecord.neighbors[firstNeighbor] != indices.back)
                || (secondRecord.neighbors[secondNeighbor] != indices.front))
        {
            correct = false;
        }
    }

    // Check against another neighbor checker, which seems to have weird
    // behavior sometimes for DMD. See comments in main loop of checkNeighbors
    version (unittest)
    {
        if (facetRecords.checkNeighbors != correct)
        {
            writeln("Found discrepancy between neighbor checkers");
            writeln("facetRecords.length = ", facetRecords.length);
            writeln("facetRecords[0].facet.length = ", facetRecords[0].facet.length);
            writeln("facetRecords.checkNeighbors2 = ", correct);
            writeln("facetRecords.checkNeighbors  = ", facetRecords.checkNeighbors);
            writeln("facetRecords = ", facetRecords);
            assert(false);
        }
    }

    return correct;
}

unittest
{

    FacetRecord!2[] facetRecords;
    foreach (facet; standardSphere(2))
    {
        FacetRecord!2 newFacetRecord;
        newFacetRecord.facet = facet;
        facetRecords ~= newFacetRecord;
    }

    assert(!facetRecords.checkNeighbors2);

    foreach (ref record; facetRecords)
    {
        // Just so happens neighbor indices are same as vertex labels
        // in this situation (manifold made from standard sphere)
        record.neighbors = record.facet;
    }
    assert(facetRecords.checkNeighbors2);

    auto hyper = loadManifold!3("data/manifold_sampler_unittest_dodecahedral.dat");
    foreach (thisIndex, ref record; hyper.facetRecords)
    {
        foreach (ref neighborIndex; record.neighbors)
        {
            auto oldNeighborIndex = neighborIndex;
            while ((neighborIndex == thisIndex) || (neighborIndex == oldNeighborIndex))
            {
                neighborIndex = (neighborIndex + 1) % hyper.facetRecords.length;
            }
            assert(neighborIndex != oldNeighborIndex);
            assert(neighborIndex != thisIndex);

            assert(!hyper.facetRecords.checkNeighbors2);
            neighborIndex = oldNeighborIndex;
            assert(hyper.facetRecords.checkNeighbors2);
        }
    }

}

/+*****************************************************************************
+ Sets all of the degree counts to their proper values. Assumes that the given
+ array of facet records has valid neighbor information.
+/
void setAllDegrees(size_t dim)(ref FacetRecord!dim[] facetRecords)
{
    static assert(dim >= 2);
    static assert(dim <= maxDimension);
    assert(facetRecords.checkNeighbors2);
    foreach (baseRecordIndex, ref baseRecord; facetRecords)
    {
        auto baseFacet = baseRecord.facet;
        foreach (hingeIndex, ref degree; baseRecord.degrees)
        {
            if (degree >= 3)
                continue; // Skip already processed hinges
            assert(degree == 0);

            // Get indices of the missing vertices in hinge of interest
            immutable missingVertexIndices = indexPairInPosition!dim(hingeIndex);

            // Get the vertex labels in hinge of interest
            auto hinge = baseFacet.hinges[hingeIndex];

            // Save initial link vertex label, so we know when to stop
            auto startingVertex = baseFacet[missingVertexIndices.front];

            // Set inital values for data tracked (as we go around hinge's star)
            auto nextVertex = baseFacet[missingVertexIndices.back];
            auto currentVertex = startingVertex;
            auto currentFacet = baseFacet;
            auto currentRecord = baseRecord;

            auto starRecordIndices = [baseRecordIndex];
            auto hingeIndices = [hingeIndex];

            // Now, we go around the star of the given hinge
            auto degreeCount = 1;
            while (nextVertex != startingVertex)
            {
                // Get and store index of next neighbor around hinge
                immutable nextNeighbor = currentFacet.neighborIndex(currentVertex);
                starRecordIndices ~= currentRecord.neighbors[nextNeighbor];

                // Retrieve information about next facet
                auto nextRecord = facetRecords[starRecordIndices.back];
                auto nextFacet = nextRecord.facet;

                // Store index of hinge in next facet
                hingeIndices ~= nextFacet.hingeIndex!dim(hinge);

                // Update stuff we track
                currentVertex = nextVertex;
                currentFacet = nextFacet;
                currentRecord = nextRecord;
                nextVertex = nextFacet.faceOpposite(hinge ~ [nextVertex]).front;

                // Increment the number of facets seen in the star
                ++degreeCount;
            }

            // Update the degree counts for facetRecords seen
            foreach (i, recordIndex; starRecordIndices)
            {
                facetRecords[recordIndex].degrees[hingeIndices[i]] = degreeCount;
            }
        }
    }
}

unittest
{

    auto sphere = Manifold!2(standardSphere(2));
    foreach (record; sphere.facetRecords)
    {
        foreach (degree; record.degrees)
        {
            assert(degree == 3);
        }
    }

    auto sphere6 = Manifold!6(standardSphere(6));
    foreach (record; sphere6.facetRecords)
    {
        foreach (degree; record.degrees)
        {
            assert(degree == 3);
        }
    }

    enum dim = maxDimension;
    auto sphereMax = Manifold!dim(standardSphere(dim));
    foreach (record; sphereMax.facetRecords)
    {
        foreach (degree; record.degrees)
        {
            assert(degree == 3);
        }
    }

    auto surface = Manifold!2([[0, 1, 2], [0, 2, 3], [0, 3, 4], [0, 4, 5], [0,
            5, 6], [0, 6, 1], [1, 2, 6], [2, 3, 5], [2, 5, 6], [3, 4, 5]]);

    auto expectedFacets = [[0, 1, 2], [0, 2, 3], [0, 3, 4], [
        0, 4, 5
    ], [0, 5, 6], [0, 6, 1], [1, 2, 6], [2, 3, 5], [2, 5, 6], [3, 4, 5]];

    auto expectedNeighbors = [[6, 1, 5], [7, 2, 0], [9, 3, 1], [
        9, 4, 2
    ], [8, 5, 3], [6, 0, 4], [8, 5, 0], [9, 8, 1], [4, 6, 7], [3, 7, 2]];

    auto expectedDegrees = [[5, 3, 6], [4, 5, 6], [3, 4, 6], [
        5, 3, 6
    ], [4, 5, 6], [3, 4, 6], [4, 5, 3], [5, 4, 5], [4, 5, 5], [5, 3, 4]];

    foreach (index, record; surface.facetRecords)
    {
        with (record)
        {
            assert(facet == expectedFacets[index]);
            assert(neighbors == expectedNeighbors[index]);
            assert(degrees == expectedDegrees[index]);
        }
    }
}

/+*****************************************************************************
+ Checks that the degree counts in a given set of facet records are as they
+ should be. Returns true if they do and fale if they don't. This involves
+ counting the number of facets in which each hinge appears as a face, and then
+ cheking that each corresponding degree value is correct in the
+ given array of facet records.
+/
bool checkDegrees(size_t dim)(const FacetRecord!dim[] facetRecords)
{
    static assert(dim >= 2);
    static assert(dim <= maxDimension);

    size_t[immutable size_t[]] hingeDegrees;

    // Build associative array containing number of times each hinge is seen
    // in a facet.
    foreach (record; facetRecords)
    {
        foreach (hinge; record.facet.hinges)
        {
            if (sort(hinge).array.idup in hingeDegrees)
            {
                ++hingeDegrees[sort(hinge).array.idup];
            }
            else
            {
                hingeDegrees[sort(hinge).array.idup] = 1;
            }
        }
    }

    // Check that all the hinges are in our associative array and their
    // corresponding degrees are correct in each facet record
    foreach (record; facetRecords)
    {
        foreach (index, degree; record.degrees)
        {
            auto hinge = sort(record.facet.hinges[index]).array.idup;
            if ((hinge !in hingeDegrees) || (hingeDegrees[hinge] != degree))
            {
                return false;
            }
        }
    }
    return true;
}

unittest
{

    enum dimension = 8;
    FacetRecord!dimension[] facetRecords;
    foreach (facet; standardSphere(dimension))
    {
        FacetRecord!dimension newFacetRecord;
        newFacetRecord.facet = facet;
        facetRecords ~= newFacetRecord;
    }

    assert(!facetRecords.checkDegrees);

    foreach (ref record; facetRecords)
    {
        foreach (ref degree; record.degrees)
        {
            degree = 3; // In standard sphere all hinges have degree three
        }
    }
    assert(facetRecords.checkDegrees);
}

/+*****************************************************************************
+ Returns the pair of indices [i, j] with 0 <= i < j <= dimension*(dimension+1)
+ which appears in given position in dictionary order.
+/
size_t[] indexPairInPosition(size_t dimension)(size_t position)
{
    // TO DO: Replace this with look-up table version computed at compile-time?

    static assert(dimension >= 2);
    static assert(dimension <= maxDimension);

    assert(position < dimension * (dimension + 1) / 2);

    auto remaining = position;
    auto toRemove = dimension;
    size_t i = 0;
    while (remaining >= toRemove)
    {
        remaining -= toRemove;
        --toRemove;
        ++i;
    }
    assert(i <= maxDimension);

    immutable j = i + remaining + 1;
    assert(j <= maxDimension);

    return [i, j];
}

unittest
{

    // dimension two
    auto expectedIndexPairs2 = [[0, 1], [0, 2], [1, 2]];
    foreach (position, indexPair; expectedIndexPairs2)
    {
        assert(position.indexPairInPosition!2 == indexPair);
    }

    // dimension three
    auto expectedIndexPairs3 = [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]];
    foreach (position, indexPair; expectedIndexPairs3)
    {
        assert(position.indexPairInPosition!3 == indexPair);
    }

    enum dim = maxDimension;
    assert(indexPairInPosition!dim(dim * (dim + 1) / 2 - 1) == [dim - 1, dim]);
    assert(indexPairInPosition!dim(dim * (dim + 1) / 2 - 2) == [dim - 2, dim]);

    assertThrown!Error(indexPairInPosition!2(3));
    assertThrown!Error(indexPairInPosition!5(15));
}

unittest
{

    // dimension four
    auto expectedIndexPairs4 = [[0, 1], [0, 2], [0, 3], [0, 4], [
        1, 2
    ], [1, 3], [1, 4], [2, 3], [2, 4], [3, 4]];
    foreach (position, indexPair; expectedIndexPairs4)
    {
        assert(position.indexPairInPosition!4 == indexPair);
    }
}

/+*****************************************************************************
+ Returns the position (in dictionary order) of a given pair of indices [i, j]
+ among all such pairs of indices with 0 <= i < j <= dimension*(dimension+1).
+/
size_t positionOfIndexPair(size_t dimension)(const size_t[] indexPair)
{
    assert(dimension >= 2);
    assert(dimension <= maxDimension);

    assert(indexPair.length == 2);
    assert(indexPair.front < indexPair.back);
    immutable i = indexPair.front;
    immutable j = indexPair.back;

    foreach (index; 0 .. dimension * (dimension + 1) / 2)
    {
        if (indexPairInPosition!dimension(index) == [i, j])
        {
            return index;
        }
    }
    assert(false, "asked for position of index-pair with invalid indices");
}

unittest
{

    // dimension two
    size_t[][] indexPairs2 = [[0, 1], [0, 2], [1, 2]];
    foreach (position, indexPair; indexPairs2)
    {
        assert(positionOfIndexPair!2(indexPair) == position);
    }

    // dimension three
    size_t[][] indexPairs3 = [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]];
    foreach (position, indexPair; indexPairs3)
    {
        assert(positionOfIndexPair!3(indexPair) == position);
    }

    // positionOfIndexPair and indexPairInPosition are inverse functions

    enum dim = 5;
    foreach (position; 0 .. dim * (dim + 1) / 2)
    {
        auto indexPair = indexPairInPosition!dim(position);
        assert(positionOfIndexPair!dim(indexPair) == position);
    }

    foreach (i; 0 .. dim + 1)
    {
        foreach (j; i + 1 .. dim + 1)
        {
            auto position = positionOfIndexPair!dim([i, j]);
            assert(indexPairInPosition!dim(position) == [i, j]);
        }
    }

    assertThrown!Error(positionOfIndexPair!2([0, 1, 2]));
    assertThrown!Error(positionOfIndexPair!2([4, 5]));
}

/******************************************************************************
* Function which returns the hinge index of given hinge in a simplex
*/
size_t hingeIndex(size_t dimension)(const Simplex facet, const Simplex hinge)
{
    static assert(dimension >= 2);
    static assert(dimension <= maxDimension);
    assert(facet.length == dimension + 1);
    assert(hinge.length == dimension - 1);
    assert(facet.hasFace(hinge));

    // Find index of first vertex of facet missing from the hinge
    immutable i = facet.countUntil!(vertex => !hinge.canFind(vertex));
    assert(i >= 0);
    assert(i < dimension);

    // Find index of second vertex of facet missing from the hinge
    immutable j = facet[i + 1 .. $].countUntil!(vertex => !hinge.canFind(vertex)) + i + 1;
    assert(j > i);
    assert(j <= dimension);

    return positionOfIndexPair!dimension([i, j]);
}
///
unittest
{

    assert([1, 3, 5].hingeIndex!2([5]) == 0);
    assert([1, 3, 5].hingeIndex!2([3]) == 1);
    assert([1, 3, 5].hingeIndex!2([1]) == 2);

    assert([2, 1, 3, 0].hingeIndex!3([0, 3]) == 0);
    assert([2, 1, 3, 0].hingeIndex!3([0, 1]) == 1);
    assert([2, 1, 3, 0].hingeIndex!3([3, 1]) == 2);
    assert([2, 1, 3, 0].hingeIndex!3([2, 0]) == 3);
    assert([2, 1, 3, 0].hingeIndex!3([2, 3]) == 4);
    assert([2, 1, 3, 0].hingeIndex!3([2, 1]) == 5);

    assertThrown!Error([1, 3, 5].hingeIndex!2([3, 5, 6]));
    assertThrown!Error([1, 3, 5].hingeIndex!2([7]));
}

/******************************************************************************
* Returns the neighbor index for the neighbor opposite a given vertex. This is
* equivalent to the position of the given vertex in this simplex's list of
* vertices.
*/
size_t neighborIndex(const Simplex simplex, Vertex vertex)
{
    assert(simplex.countUntil(vertex) != -1, "tried to find index of a vertex not in given simplex");
    return simplex.countUntil(vertex);
}
///
unittest
{

    assert([3, 2, 7].neighborIndex(3) == 0);
    assert([3, 2, 7].neighborIndex(2) == 1);
    assert([3, 2, 7].neighborIndex(7) == 2);
    assert([5, 4, 1].neighborIndex(5) == 0);
    assert([11].neighborIndex(11) == 0);
    assertThrown!Error([1, 3, 5].neighborIndex(7));
}

/******************************************************************************
* Returns whether firstSimplex and secondSimplex are the same or not. Note that
* this is not the same as firstSimplex == secondSimplex since that requires the
* vertices in firstSimplex and secondSimplex to appear in the same order.
*/
bool isEqual(const Simplex firstSimplex, const Simplex secondSimplex)
{
    return sort(firstSimplex.dup) == sort(secondSimplex.dup);
}
///
unittest
{

    assert(isEqual([1, 3, 2], [3, 2, 1]));
    assert(!isEqual([5, 7], [5, 7, 1]));
    assert([1, 4, 6, 3].isEqual([3, 4, 1, 6]));
    assert(isEqual([], []));
    assert(!isEqual([], [6]));
}

/******************************************************************************
* Returns whether a simplex is degenerate. That is, whether or not there are
* repeated vertex labels in the simplex.
*/
bool isDegenerate(const Simplex simplex)
{
    return sort(simplex.dup).findAdjacent.length != 0;
}
///
unittest
{

    assert([1, 1, 3, 7, 6].isDegenerate);
    assert([7, 7].isDegenerate);
    assert(!standardSimplex(maxDimension).isDegenerate);
    assert(!isDegenerate([13]));
    assert(!isDegenerate([]));
}

/******************************************************************************
* Returns the opposite face of a given face in a simplex.
*/
Simplex faceOpposite(in Simplex simplex, in Simplex face)
{
    assert(simplex.hasFace(face));
    return simplex.filter!(vertex => !face.canFind(vertex)).array.dup;
}
///
unittest
{

    assert([4, 1, 0, 7].faceOpposite([0]).isEqual([7, 4, 1]));
    assert([1, 3, 2].faceOpposite([1, 3]).isEqual([2]));
    assert([1, 2].faceOpposite([1]) != [1]);
    assert([3].faceOpposite([3]) == []);
    assert([3, 1, 4].faceOpposite([]).isEqual([3, 1, 4]));
    assert([].faceOpposite([]) == []);
    assertThrown!Error([2, 3, 5].faceOpposite([7]));
    assertThrown!Error([2, 3, 5].faceOpposite([3, 5, 1]));
}

/******************************************************************************
* Returns the opposite face of a given face in a simplex.
*/
Simplex commonFace(in Simplex first, in Simplex second)
{
    return first.filter!(vertex => second.canFind(vertex)).array.dup;
}
///
unittest
{

    assert(commonFace([3, 5, 0], [0, 1, 5]).isEqual([5, 0]));
    assert(commonFace([3, 1, 7, 5], [7, 2, 4]) == [7]);
    assert(commonFace([3, 1, 7], [0, 2, 4, 5]) == []);
    assert(commonFace([], [0, 2]) == []);
    assert(commonFace([1], [1]) == [1]);
    assert(commonFace([1, 2, 3], []) == []);
    assert(commonFace([], []) == []);
}

/******************************************************************************
* Returns a manifold loaded from the file specified by fileName. If fileName
* is the empty string, the returned manifold is the standard sphere.
*/
Manifold!dim loadManifold(size_t dim = dimManifold)(string fileName)
{
    auto manifoldFile = File(fileName, "r"); // Open file in read-only mode

    string facets;
    foreach (string line; manifoldFile.lines)
    {
        // We allow comments starting with '#'
        if (line.front != '#')
        {
            facets ~= line.strip;
        }
    }

    try
    {
        auto loaded = Manifold!dim(facets.parse!(Simplex[]));
        return loaded;
    }
    catch (Exception ex)
    {
        ex.msg ~= "\n\n\t ERROR: encountered malformed facet list in initial manifold triangulation file: "
            ~ fileName ~ "\n";
        throw ex;
    }
}
///
unittest
{

    auto manifold = loadManifold!2("data/manifold_sampler_unittest_load.dat");
    auto expectedManifold = Manifold!2([[0, 1, 2], [0, 2, 3], [0, 3, 4], [0, 4,
            5], [0, 5, 6], [0, 6, 1], [1, 2, 6], [2, 3, 5], [2, 5, 6], [3, 4, 5]]);
    assert(manifold == expectedManifold);

    assertThrown(loadManifold!2("data/manifold_sampler_unittest_bad_load.dat"));
}

/******************************************************************************
* Saves a manifold to file specified by fileName.
*/
void saveTo(size_t dimension)(Manifold!dimension manifold, string fileName)
{
    auto saveFile = File(fileName, "w"); // Open in write-only mode

    saveFile.writeln("# created ", Clock.currTime.to!DateTime);
    saveFile.write("[");
    foreach (record; manifold.facetRecords[0 .. $ - 1])
    {
        saveFile.write(record.facet, ",");
    }
    saveFile.writeln(manifold.facetRecords.back.facet, "]");
}
///
unittest
{

    auto fileName = "data/manifold_sampler_unittest_save.dat";
    auto sphere = Manifold!4(standardSphere(4));
    sphere.saveTo(fileName);
    Manifold!4 loadedSphere = loadManifold!4(fileName);
    assert(loadedSphere == sphere);
}

/+*****************************************************************************
+ Applies an action to each facetRecord whose facet is in the star of a given
+ simplex. NOTE: This is done by iterating over *all* the facet records in the
+ manifold. Do not make any assumptions about the order in which the action is
+ applied to the facetRecords.
+/
void doInStarByIteration(size_t dim)(Manifold!dim manifold, in Simplex simplex,
        void delegate(FacetRecord!dim) action)
{
    foreach (record; manifold.facetRecords)
    {
        if (record.facet.hasFace(simplex))
        {
            action(record);
        }
    }
}

unittest
{

    auto sphere = Manifold!3(standardSphere(3));
    int count = 0;
    sphere.doInStarByIteration!3([], (record) { ++count; });
    assert(count == 5);

    auto octahedron = Manifold!2([[0, 1, 2], [0, 2, 3], [0, 3, 4], [0, 4, 1],
            [5, 1, 2], [5, 2, 3], [5, 3, 4], [5, 4, 1]]);

    foreach (vertex; 0 .. 6)
    {
        count = 0;
        octahedron.doInStarByIteration!2([vertex], (record) { ++count; });
        assert(count == 4);
    }

    auto surface = Manifold!2([[0, 1, 2], [0, 2, 3], [0, 3, 4], [0, 4, 5], [0,
            5, 6], [0, 6, 1], [1, 2, 6], [2, 3, 5], [2, 5, 6], [3, 4, 5]]);

    count = 0;
    surface.doInStarByIteration!2([0], (record) { ++count; });
    assert(count == 6);

    count = 0;
    surface.doInStarByIteration!2([1], (record) { ++count; });
    assert(count == 3);

    count = 0;
    surface.doInStarByIteration!2([2], (record) { ++count; });
    assert(count == 5);
}

void doInStar(size_t dim)(Manifold!dim manifold, const Simplex simplex,
        void delegate(FacetRecord!dim) action, FacetIndex index)
{
    static assert(dim >= minDimension);
    static assert(dim <= maxDimension);
    assert(index < manifold.facetRecords.length);
    assert(manifold.facetRecords[index].facet.hasFace(simplex),
            "doInStar called with simplex not the face of facet at index given");

    FacetIndex[] indicesSeen = [index];
    FacetRecord!dim[] recordsToCheck = [manifold.facetRecords[index]];
    FacetRecord!dim[] recordsFound;

    while (recordsToCheck.length > 0)
    {
        // Get record to process, removing it from those left to check, and
        // adding it to those found
        auto record = recordsToCheck.back;
        recordsToCheck.popBack;
        recordsFound ~= record;

        // Now consider which neighboring facet records must be processed
        foreach (ridgeIndex, ridge; record.facet.ridges)
        {
            immutable neighborIndex = record.neighbors[ridgeIndex];
            if (ridge.hasFace(simplex) && !indicesSeen.canFind(neighborIndex))
            {
                indicesSeen ~= neighborIndex;
                recordsToCheck ~= manifold.facetRecords[neighborIndex];
            }
        }
    }

    foreach (record; recordsFound)
    {
        action(record);
    }
}

unittest
{

    auto sphere = Manifold!2(standardSphere(2));
    int count = 0;

    // Count the facets in the star of vertex 1 starting with facet #0)
    // Here x is the record under examination (which we don't use here.)
    sphere.doInStar!2([1], (r) { ++count; }, 0);
    assert(count == 3);

    auto octahedron = Manifold!2([[0, 1, 2], [0, 2, 3], [0, 3, 4], [0, 4, 1],
            [5, 1, 2], [5, 2, 3], [5, 3, 4], [5, 4, 1]]);

    Simplex[] starFacets;
    Simplex[] starFacets2; // For cross-checking with by-iteration algorithm

    octahedron.doInStar!2([0, 1], (r) { starFacets ~= r.facet.dup; }, 0);
    octahedron.doInStarByIteration!2([0, 1], (r) { starFacets2 ~= r.facet.dup; });
    assert(starFacets == [[0, 1, 2], [0, 4, 1]]);
    assert(equal(sort(starFacets2), sort(starFacets))); // Cross-check

    starFacets = [];
    starFacets2 = [];
    octahedron.doInStar!2([5], (r) { starFacets ~= r.facet.dup; }, 7);
    octahedron.doInStarByIteration!2([5], (r) { starFacets2 ~= r.facet.dup; });
    assert(starFacets == [[5, 4, 1], [5, 3, 4], [5, 2, 3], [5, 1, 2]]);
    assert(equal(sort(starFacets2), sort(starFacets)));

    auto surface = Manifold!2([[0, 1, 2], [0, 2, 3], [0, 3, 4], [0, 4, 5], [0,
            5, 6], [0, 6, 1], [1, 2, 6], [2, 3, 5], [2, 5, 6], [3, 4, 5]]);

    starFacets = [];
    starFacets2 = [];
    surface.doInStar!2([2, 5, 6], (r) { starFacets ~= r.facet.dup; }, 8);
    surface.doInStarByIteration!2([2, 5, 6], (r) { starFacets2 ~= r.facet.dup; });
    assert(starFacets == [[2, 5, 6]]);
    assert(starFacets2 == starFacets); // No need to sort here

    starFacets = [];
    starFacets2 = [];
    surface.doInStar!2([6], (r) { starFacets ~= r.facet.dup; }, 8);
    surface.doInStarByIteration!2([6], (r) { starFacets2 ~= r.facet.dup; });
    assert(starFacets == [[2, 5, 6], [1, 2, 6], [0, 6, 1], [0, 5, 6]]);
    assert(equal(sort(starFacets2), sort(starFacets)));

    auto hyperbolicDodecahedral = loadManifold!3("data/manifold_sampler_unittest_dodecahedral.dat");

    size_t vertexDegree;
    hyperbolicDodecahedral.doInStar!3([19], (r) { ++vertexDegree; }, 168);
    assert(vertexDegree == 34);

    vertexDegree = 0;
    hyperbolicDodecahedral.doInStarByIteration!3([19], (r) { ++vertexDegree; });
    assert(vertexDegree == 34);
}

bool anyInStar(size_t dim)(Manifold!dim manifold, const Simplex simplex,
        bool delegate(FacetRecord!dim) predicate, FacetIndex index)
{
    static assert(dim >= minDimension);
    static assert(dim <= maxDimension);
    assert(index < manifold.facetRecords.length);
    assert(manifold.facetRecords[index].facet.hasFace(simplex),
            "anyInStar called with simplex not the face of facet at index given");

    FacetIndex[] indicesSeen = [index];
    FacetRecord!dim[] recordsToCheck = [manifold.facetRecords[index]];

    while (recordsToCheck.length > 0)
    {
        // Get record to process, check if it satisfies the predicate...
        auto record = recordsToCheck.back;
        if (predicate(record))
        {
            return true;
        }
        //... if not, add remove it from those left to check
        recordsToCheck.popBack;

        // Now consider which neighboring facet records must be processed
        foreach (ridgeIndex, ridge; record.facet.ridges)
        {
            immutable neighborIndex = record.neighbors[ridgeIndex];
            if (ridge.hasFace(simplex) && !indicesSeen.canFind(neighborIndex))
            {
                indicesSeen ~= neighborIndex;
                recordsToCheck ~= manifold.facetRecords[neighborIndex];
            }
        }
    }
    return false;
}

unittest
{

    auto surface = Manifold!2([[0, 1, 2], [0, 2, 3], [0, 3, 4], [0, 4, 5], [0,
            5, 6], [0, 6, 1], [1, 2, 6], [2, 3, 5], [2, 5, 6], [3, 4, 5]]);

    assert(surface.anyInStar!2([0, 1], (r) => r.facet.hasFace([6, 1]), 5) == true);
    assert(surface.anyInStar!2([0, 5, 6], (r) => r.degrees[0] == 5, 4) == false);

    size_t nEvals = 0;
    assert(surface.anyInStar!2([6], (r) { ++nEvals; return r.facet.hasFace([6]); }, 8) == true);
    // anyInStar should return immediately after first record for which
    // predicate(record) is true. Here, all facets have vertex 6 as face so
    // only one evaluation of predicate should be made.
    assert(nEvals == 1);

    nEvals = 0;
    assert(surface.anyInStar!2([0], (r) {
            ++nEvals;
            return r.facet.hasFace([5, 3]);
        }, 3) == false);
    // none of the six facets satisfy predicate but all six should have been
    // checked.
    assert(nEvals == 6);
}

Simplex[] getStar(size_t dim)(Manifold!dim manifold, in Simplex simplex, FacetIndex index)
{
    Simplex[] facets;
    manifold.doInStar!dim(simplex, (r) { facets ~= r.facet.dup; }, index);
    return facets;
}

unittest
{

    auto surface = Manifold!2([[0, 1, 2], [0, 2, 3], [0, 3, 4], [0, 4, 5], [0,
            5, 6], [0, 6, 1], [1, 2, 6], [2, 3, 5], [2, 5, 6], [3, 4, 5]]);

    assert(surface.getStar!2([0, 1], 5) == [[0, 6, 1], [0, 1, 2]]);
    assert(surface.getStar!2([0, 5, 6], 4) == [[0, 5, 6]]);
    assert(surface.getStar!2([6], 8) == [[2, 5, 6], [1, 2, 6], [0, 6, 1], [0, 5, 6]]);
}

/******************************************************************************
* Returns a list of all the faces of a given simples. They are returned sorted
* by dimension (lowest to highest).
*/
Simplex[] allFaces(in Simplex simplex)
{
    Simplex[] faces;
    foreach (bitChoice; 1 .. 2 ^^ simplex.length)
    {
        Simplex face;
        foreach (index, vertex; simplex)
        {
            if ((1 << index) & bitChoice)
            {
                face ~= vertex;
            }
        }
        faces ~= face.dup;
    }
    faces.sort!((a, b) => a.length < b.length);
    return faces;
}

unittest
{
    assert([].allFaces == []);
    assert([9].allFaces == [[9]]);
    assert([7, 5].allFaces == [[7], [5], [7, 5]]);
    assert([1, 2, 3].allFaces == [[1], [2], [3], [1, 2], [1, 3], [2, 3], [1, 2, 3]]);
}

/******************************************************************************
* Returns the f-vector of the manifold. The i-th value in the returned array
* gives the number of i-dimensional simplices in the manifold.
*/
size_t[] fVector(size_t dimension)(in Manifold!dimension manifold)
{
    assert(dimension <= maxDimension);
    bool[immutable Simplex] simplicesSeen; // only keys matter, values unused
    size_t[dimension + 1] fVec;

    foreach (record; manifold.facetRecords)
    {
        foreach (simplex; record.facet.allFaces)
        {
            if (sort(simplex).array !in simplicesSeen)
            {
                simplicesSeen[sort(simplex).array.idup] = true;
                ++fVec[simplex.length - 1];
            }
        }
    }
    return fVec.dup;
}
///
unittest
{

    auto octahedron = Manifold!2([[0, 1, 2], [0, 2, 3], [0, 3, 4], [0, 4, 1],
            [5, 1, 2], [5, 2, 3], [5, 3, 4], [5, 4, 1]]);

    assert(octahedron.fVector == [6, 12, 8]);

    auto sphere7 = Manifold!7(standardSphere(7));

    // fVector should be (9 choose 1), (9 choose 2), ... , (9 choose 8)
    assert(sphere7.fVector == [9, 9 * 8 / 2, 9 * 8 * 7 / (3 * 2),
            9 * 8 * 7 * 6 / (4 * 3 * 2), 9 * 8 * 7 * 6 / (4 * 3 * 2),
            (9 * 8 * 7) / (3 * 2), (9 * 8) / 2, 9]);

    auto hyperbolicDodecahedral = loadManifold!3("data/manifold_sampler_unittest_dodecahedral.dat");

    // TO DO: Get reference for this...
    assert(hyperbolicDodecahedral.fVector == [21, 190, 338, 169]);

}

/******************************************************************************
* Returns...
*/
size_t[][] getPotentialPachnerMoves(size_t dim)(const ref FacetRecord!dim facetRecord)
{
    return [];
}
///
unittest
{

    //auto testRecord = FacetRecord!3([1, 2, 3, 4], [0, 0, 0, 0], [3, 3, 3, 3, 3, 3]);
    //testRecord.getPotentialPachnerMoves!3.writeln;
}

// TO DO: Function doing a pachner move at given facetRecord and face
// TO DO: Function for removing a facetRecord from manifold
// TO DO: Function for doing special moves like 4 -> 4 in dim=3. Other dims?
// TO DO: Logic for choosing random pachner move and accept/reject
