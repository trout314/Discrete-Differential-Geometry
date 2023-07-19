/// TO DO: Module description
module simplicial_complex;

import std.algorithm, std.conv, std.datetime, std.exception, std.random,
    std.range, std.stdio, std.typecons, std.traits;
import unit_threaded;
import algorithms, utility;

alias isIRof = isInputRangeOf;
alias isIRofIRof = isInputRangeOfInputRangeOf;

/// Basic Functionality
@Name("doc tests") pure @safe unittest
{
    // create an empty simplicial complex with vertices of default type `int`
    SimplicialComplex!() sc;
    assert(sc.numFacets == 0);
    assert(sc.facets.empty);

    // get the type of vertex used
    static assert(is(sc.Vertex == int));

    /* while you can specify the vertex type used as a template parameter,
    we will stick to `int` for now */
    static assert(is(SimplicialComplex!int == typeof(sc)));
    
    // insert some facets
    sc.insertFacet([1,2,3]);
    sc.insertFacet([2,3,4]);
    sc.insertFacet([4,5]);
    sc.insertFacet([5,6]);

    /* can also create simplicial complexes directly from an array of arrays of
    vertices which specify the facets */
    auto sc2 = SimplicialComplex!()([[4,5], [5,6], [1,2,3], [2,3,4]]);

    /* can compare simplicial complexes for equality. NOTE: this allocates new
    memory to copy all facets into, then sorts each array of facets. */
    assert(sc == sc2);

    // a helper function template that returns a newly constructed complex
    auto sc3 = simplicialComplex([[5,6], [2,3,4], [1,2,3], [4,5]]);
    assert(sc3 == sc);

    // simplicial complexes are value types
    auto sc4 = sc3;
    sc4.insertFacet([6,7]);
    assert(sc3 == sc);
    assert(sc4 != sc);

    /* get the vertices in the facets, returned in order of increasing dimension
    and dictionary order within a dimension */
    sc.facets.shouldBeSameSetAs([[4,5], [5,6], [1,2,3], [2,3,4]]);

    // Can get a sorted array of all the facet simplices of a given dimension
    assert(sc.facets(0).empty);
    sc.facets(1).shouldBeSameSetAs([[4,5], [5,6]]);
    sc.facets(2).shouldBeSameSetAs([[1,2,3], [2,3,4]]);
    assert(sc.facets(3).empty);

    // Can get a sorted array of all the simplices of a given dimension
    sc.simplices(0).shouldBeSameSetAs([[1], [2], [3], [4], [5], [6]]);
    sc.simplices(1).shouldBeSameSetAs([[1,2], [1,3], [2,3], [2,4], [3,4], 
        [4,5], [5,6]]);
    sc.simplices(2).shouldBeSameSetAs(sc.facets(2));
    assert(sc.simplices(3).empty);

    // get the f-vector of the simplicial complex
    assert(sc.fVector == [6,7,2]);

    // check for the presence of simplices
    assert(sc.contains([4,5]));
    assert(sc.contains([2,3,4]));
    assert(sc.contains([6]));
    assert(!sc.contains([1,2,5]));
    assert(!sc.contains([4,6]));

    // get the star of a simplex an array of facets
    sc.star([5]).shouldBeSameSetAs([[4,5], [5,6]]);
    sc.star([4]).shouldBeSameSetAs([[4,5], [2,3,4]]);
    sc.star([2,3]).shouldBeSameSetAs([[1,2,3], [2,3,4]]);

    /* get link of a simplex as list of facets. NOTE: map!array is needed here
    due to library bug. See https://github.com/gedaiu/fluent-asserts/issues/85 */
    sc.link([5]).map!array.shouldBeSameSetAs([[4], [6]]);
    sc.link([4]).map!array.shouldBeSameSetAs([[5], [2,3]]);
    sc.link([2,3]).map!array.shouldBeSameSetAs([[1], [4]]);

    // can print out a simplicial complex
    assert(sc.toString == "[[4, 5], [5, 6], [1, 2, 3], [2, 3, 4]]");

    // can construct from a range of ranges of vertices
    auto sc5 = SimplicialComplex!()(iota(3).map!(i => iota(3*i, 3*i + 3)));
    sc5.facets.shouldBeSameSetAs([[0, 1, 2], [3, 4, 5], [6, 7, 8]]);
}

///
@Name("toHash (@safe)") @safe unittest
{
    auto sc1 = simplicialComplex([[4,5], [1,2,3], [2,3,4], [5,6]]);
    auto sc2 = simplicialComplex([[4,5], [5,6], [2,3,4], [1,2,3]]);  
    assert(sc1.toHash == sc2.toHash);   // Note: this isn't pure

    sc1.insertFacet([6,7]);
    // Note: this assert *could* fail, but useful check anyway!
    assert(sc1.toHash != sc2.toHash);
} 

///
@Name("works as AA KeyType") pure @safe unittest
{
    int[SimplicialComplex!()] aa;
    auto sc1 = simplicialComplex([[4,5]]);
    auto sc2 = simplicialComplex([[4,5], [5,6]]);
    auto sc3 = simplicialComplex([[4,5], [5,6], [1,2,3]]);

    aa[sc1] = 1;
    aa[sc2] = 2;

    assert(sc1 in aa);
    assert(sc2 in aa);
    assert(sc3 !in aa);

    aa[sc3] = 3;
    assert(sc3 in aa);    
}

/// Some restrictions
@Name("doc tests (errors)") pure @system unittest
{
    auto sc = simplicialComplex([[4,5], [5,6], [1,2,3], [2,3,4]]);

    // cannot insert a simplex if it is already the face of an existing facet
    sc.insertFacet([4,5]).throwsWithMsg(
        "expected a simplex not already in the simplicial complex");

    sc.insertFacet([2]).throwsWithMsg(
        "expected a simplex not already in the simplicial complex");

    // cannot remove a facet that isn't part of the simplicial complex
    simplicialComplex([[3,4]]).removeFacet([2,3]).throwsWithMsg(
        "tried to remove a facet not in the simplicial complex");

    // not allowed to ask for link of a simplex not in the simplicial complex
    simplicialComplex([[1,3,4], [2,3,4]]).link([1,2]).throwsWithMsg(
        "expected a simplex in the simplicial complex");

    sc.facets(-1).throwsWithMsg("expected a non-negative dimension");
}

@Name("additional tests") pure @safe unittest
{
    alias sComp = simplicialComplex;

    // ---------------- TEST A DEFAULT INITIALIZED COMPLEX ---------------------
    auto sc = SimplicialComplex!()();
    static assert(is(sc.Vertex == int)); // Default vertex type is int

    sc.facets.empty.shouldEqual(true);
    sc.numFacets.shouldEqual(0);
  
    // ---------------------------- ADD SOME FACETS ----------------------------
    sc.insertFacet([1, 2]);
    sc.insertFacet([2, 3]);
    sc.insertFacet([3, 4, 5]);
    
    sc.shouldEqual(sComp([[1,2], [2,3], [3,4,5]]));

    // -------------------------- TEST FUNCTIONALITY ---------------------------
    sc.facets.walkLength.shouldEqual(3);
    sc.numFacets.shouldEqual(3);

    sc.facets(0).empty.shouldEqual(true);
    sc.facets(1).shouldBeSameSetAs([[1, 2], [2, 3]]);
    sc.facets(2).shouldBeSameSetAs([[3, 4, 5]]);

    sc.simplices(0).shouldBeSameSetAs([[1], [2], [3], [4], [5]]);
    sc.simplices(1).shouldBeSameSetAs([[1, 2], [2, 3], [3, 4], [3, 5], [4, 5]]);
    sc.simplices(2).shouldBeSameSetAs([[3, 4, 5]]);

    foreach (d; iota(3, 16))
    {
        sc.facets(d).empty.shouldEqual(true);
        sc.simplices(d).empty.shouldEqual(true);
    }
}

@Name("insertFacet") pure @safe unittest
{
    auto sc = SimplicialComplex!()();
    sc.insertFacet([1,2]);
    sc.facets.shouldBeSameSetAs([[1,2]]);
    sc.insertFacet([2,3,4]);
    sc.facets.shouldBeSameSetAs([[1,2], [2,3,4]]);
    sc.insertFacet([1,5,6]);
    sc.facets.shouldBeSameSetAs([[1,2], [1,5,6], [2,3,4]]);
    sc.insertFacet([1,5,6,7]);
    sc.facets.shouldBeSameSetAs([[1,2], [2,3,4], [1,5,6,7]]);
}

@Name("insertFacet/removeFacet") pure @safe unittest
{
    auto sc = simplicialComplex([[1,2], [1,3], [1,4], [4,5,6]]);

    sc.facets.shouldBeSameSetAs([[1,2], [1,3], [1,4], [4,5,6]]);
    sc.removeFacet([1,3]);
    sc.facets.shouldBeSameSetAs([[1,2], [1,4], [4,5,6]]);
    sc.removeFacet([1,2]);
    sc.facets.shouldBeSameSetAs([[1,4], [4,5,6]]);
    sc.removeFacet([4,5,6]);
    sc.facets.shouldBeSameSetAs([[1,4]]);
    sc.removeFacet([1,4]);
    assert(sc.facets.empty);

    auto sc2 = simplicialComplex([[1,2,3], [1,2,4], [1,3,4], [2,3,4]]);
    sc2.removeFacet([1,2,3]);
    sc2.facets.shouldBeSameSetAs([[1,2,4], [1,3,4], [2,3,4]]);
    sc2.insertFacet([1,2,5]);
    sc2.insertFacet([1,3,5]);
    sc2.insertFacet([2,3,5]);
    sc2.facets.shouldBeSameSetAs([[1,2,4], [1,2,5], [1,3,4], [1,3,5], [2,3,4],
        [2,3,5]]);

    SimplicialComplex!int sc3;
    sc3.insertFacet([1,2,3]);
    assert(sc3.facets.equal([[1,2,3]]));
    sc3.removeFacet([1,2,3]);
    assert(sc3.facets.empty);
}

@Name("insertFacet!(No.checkForFacetFaces)") pure @safe unittest
{
    SimplicialComplex!() sc;
    sc.insertFacet!(No.checkForFacetFaces)([1,2]);
    sc.insertFacet!(No.checkForFacetFaces)([3,4]);
    sc.facets.shouldBeSameSetAs([[1,2],[3,4]]);

    sc.insertFacet!(No.checkForFacetFaces)([2,3,5]);
    sc.facets.shouldBeSameSetAs([[1,2],[3,4],[2,3,5]]);
}

/*******************************************************************************
Checks if a range or array of vertices represents a valid simplex or not.
Throws error if anything is wrong.
*/
void assertValidSimplex(V)(V vertices_, int dim) if (isInputRange!V)
{
    assert(vertices_.walkLength == dim + 1, "wrong number of vertices");
    static if (is(ElementType!V == class))
    {
        assert(vertices_.all!(v => v !is null),
            "null class references not allowed");
    }
    assert(vertices_.isSorted, "vertices must occur in increasing order");
    assert(vertices_.findAdjacent.walkLength == 0, "vertices must be distinct");
}


/*******************************************************************************
A simplicial complex type whose vertices are of type `Vertex`.
*/
struct SimplicialComplex(Vertex_ = int, int maxFacetDimension = 16)
{
private:
    /* We store arrays containing all the vertices in the facets of a given 
    dimension. These arrays are flattened, so the array storing the 2-dim
    facets has the form [v1, v2, v3, w1, w2, w3, ... ] where [v1, v2, v3] give
    the first simplex, [w1, w2, w3] the second, etc.
    */
    SmallMap!(int, Vertex[]) facetVertices;

    alias NSimplex = StackArray!(Vertex_, maxFacetDimension + 1);
    static NSimplex toNSimp(R)(R range)
    if (isInputRange!R && is(ElementType!R : Vertex_))
    {
        NSimplex toReturn;
        range.each!(r => toReturn ~= r);
        return toReturn;
    }

    // indexOfFacet[f] = i means the vertices of f begins at facetVertices[i]
    size_t[NSimplex] indexOfFacet;
public:
    /***************************************************************************
    The type of the vertices in this simplicial complex.
    */
    alias Vertex = Vertex_;

    /***************************************************************************
    Construct a simplicial complex from an input range of input ranges with
    element type implicitly convertible to Vertex
    */
    this(R)(R facets) if (isIRofIRof!(R, const(Vertex)))
    {
        facets.each!(f => this.insertFacet(f));
    }

    // Postblit makes sure copy doesn't share data
    this(this) pure @safe
    {
        facetVertices = facetVertices.dup;
        facetVertices.byKey.each!(k => facetVertices[k] = facetVertices[k].dup);
        indexOfFacet = indexOfFacet.dup;
    }

    /***************************************************************************
    Inserts a facet (given as an input range of vertices) into the simplicial
    complex.
    */
    void insertFacet(Flag!"checkForFacetFaces" checkForFacetFaces = Yes.checkForFacetFaces, S)(S simplex)
    if (isIRof!(S, const(Vertex)))
    {
        NSimplex simplex_ = toNSimp(simplex);
        int dim = simplex_.length.to!int - 1;
        assert(dim >= 0, simplex_.length.to!string);
        simplex_[].assertValidSimplex(dim);
        assert(!this.contains(simplex),
            "expected a simplex not already in the simplicial complex");

        /* We must remove any existing facets which are faces of the inserted
        facet. Also, we need independent copies of the facets to remove. */

        static if (checkForFacetFaces)
        {
            simplex_[].subsets
                .filter!(s => containsFacet(s))
                .each!(s => removeFacet(s));
        }
        else
        {
            assert(!simplex_[].subsets.any!(f => this.containsFacet(f)),
                    "expected a simplex without any faces that are facets");
        }
        
        if (dim in facetVertices)
        {
            facetVertices[dim] ~= simplex_;
            assert(simplex_ !in indexOfFacet);                       
            indexOfFacet[simplex_] = facetVertices[dim].length - dim - 1;
        }
        else
        {
            facetVertices.insert(dim, simplex_[].dup);
            indexOfFacet[simplex_] = 0;
        }
    }

    /***************************************************************************
    Removes an existing facet of the simplicial complex a facet (given as an 
    input range of vertices)
    */
    void removeFacet(S)(S simplex) if (isIRof!(S, const(Vertex)))
    {
        NSimplex simplex_ = toNSimp(simplex);

        int dim = simplex_[].length.to!int - 1;
        assert(dim >= 0);
        simplex_[].assertValidSimplex(dim);

        assert(this.contains(simplex),
            "tried to remove a facet not in the simplicial complex");
            
        auto indx = indexOfFacet[simplex_];

        // copy last facet into removed ones spot        
        copy(facetVertices[dim][$ - dim - 1  .. $],
             facetVertices[dim][indx .. indx + dim + 1]);

        indexOfFacet[toNSimp(facetVertices[dim][$ - dim - 1  .. $])] = indx;

        facetVertices[dim] = facetVertices[dim][0 .. $ - dim - 1];
        indexOfFacet.remove(simplex_);        
    }

    /***************************************************************************
    Returns the facets of the simplicial complex of a particular dimension as an 
    array of arrays of vertices.
    */
    auto facets(int dim) const
    {
        assert(dim >= 0, "expected a non-negative dimension");

        if (dim !in facetVertices)
        {
            typeof(facetVertices[dim]) empty;
            return empty.chunks(dim + 1);
        }
        else
        {
            return facetVertices[dim].chunks(dim + 1);
        }
    }

    /***************************************************************************
    Returns the facets of the simplicial complex. These are simplicies that 
    are not the face of another simplex. They are returned in increasing order 
    of dimension but in an implementation defined order within dimension.
    */
    auto facets() const
    {
       return FacetRange!Vertex(this.facetVertices);
    }

    /***************************************************************************
    Returns the number of facets
    */
    auto numFacets() const
    {
        return this.facets.walkLength;
    }

    /***************************************************************************
    Returns the star of the given simplex a forward range of arrays of vertices
    giving the facets.
    */
    auto link(S)(S simplex) const if (isIRof!(S, const(Vertex)))
    {
        assert(this.contains(simplex), "expected a simplex in the simplicial complex");
        return LinkRange!Vertex(StarRange!Vertex(simplex, this.facets));
    }

    /***************************************************************************
    Returns the star of the given simplex a forward range of arrays of vertices
    giving the facets.
    */ 
    auto star(S)(S simplex) const if (isIRof!(S, const(Vertex)))
    {
        return StarRange!Vertex(simplex, this.facets);
    }

    /***************************************************************************
    Returns true if the given simplex is in this simplicial complex and false
    otherwise.
    */
    bool contains(S)(S simplex) const if (isIRof!(S, const(Vertex)))
    {
        return this.facets.any!(f => simplex.isSubsetOf(f));
    }

    /***************************************************************************
    Returns true if the given simplex is a facet in this simplicial complex and
    false otherwise.
    */
    bool containsFacet(F)(F facet) const if (isIRof!(F, const(Vertex)))
    {
        return cast(bool) (toNSimp(facet) in indexOfFacet);
    }

    /***************************************************************************
    Returns a range containing a randomly chosen facet of dimension `dim`
    */
    const(Vertex)[] randomFacetOfDim()(int dim) const
    {
        assert(dim in this.facetVertices);
        size_t nVerts = (dim + 1).to!size_t;
        auto indx = uniform(0, facetVertices[dim].length / nVerts);
        return facetVertices[dim][indx * nVerts .. (indx + 1) * nVerts];
    }

    /***************************************************************************
    Get the simplices of dimension `dim` as lists of vertices.
    */
    const(Vertex)[][] simplices()(int dim) const pure nothrow @safe
    {
        // TO DO: Reduce gc presure here. (Can't make @nogc I think.)

        assert(dim >= 0, "dimension must be non-negative, but got "
            ~ dim.to!string);

        const(Vertex)[][] simplicesSeen;
        auto dims = facetVertices.byKey;
        foreach (d; dims.filter!(d => d >= dim))
        {
            foreach (f; this.facets(d))
            {
                simplicesSeen ~= f.subsetsOfSize(dim + 1).map!array.array;
            }
        }
        simplicesSeen = simplicesSeen.sort.uniq.array;
        return simplicesSeen;
    }

    bool opEquals()(auto ref const(SimplicialComplex!Vertex) sc) const
    {
        auto dims = sc.facetVertices.byKey;
        auto thisDims = this.facetVertices.byKey;
        if (dims.walkLength != thisDims.walkLength)
        {
            return false;
        }
        foreach (dim, thisDim; zip(dims, thisDims))
        {
            if (dim != thisDim)
            {
                return false;
            }
            const(Vertex)[][] f = sc.facets(dim).map!array.array.sort.array;
            const(Vertex)[][] thisF = this.facets(dim).map!array.array.sort.array;
            if (f != thisF)
            {
                return false;
            }
        }

        return true;
    }

    size_t toHash() const nothrow @safe
    {
        auto f = facets.map!array.array.sort.array;
        return () @trusted {return typeid(f).getHash(&f);} ();
    }

    /*******************************************************************************
    Returns a nice looking representation of the simplicial complex as a string.
    */
    string toString() const pure @safe
    {
        return this.facets.to!string;
    }

    /*******************************************************************************
    Returns a string containing more detailed, implementation specific info on the
    simplicial complex
    */
    string toDetailedString()() const
    {
        string output;
        output ~= "facets   : " ~ facets.to!string ~ "\n";
        output ~= "dims     : " ~ facetVertices.byKey.to!string ~ "\n";
        foreach (dim; facetVertices.byKey)
        {
            output ~= "vertex list (dim " ~ dim.to!string ~ "): "
                ~ facetVertices[dim].to!string ~ "\n";
        }
        output ~= indexOfFacet.to!string; 

        return output;
    }

    void toString(W)(ref W writer) const if (isOutputRange!(W, char)) 
    {
        put(writer, facets.toString);
    }
}

///
@Name("randomFacetOfDim") @safe unittest
{
    auto sc = simplicialComplex([[1,2,3]]);
    assert(sc.randomFacetOfDim(2) == [1,2,3]);

    sc.insertFacet([4,5,6]);
    sc.insertFacet([7]);

    int count;
    foreach (i; 0 .. 100)
    {
        auto f = sc.randomFacetOfDim(2);  
        if (f == [1,2,3])
        {
            count++;
        }
        else
        {
            assert(f == [4,5,6]);
        }
    }
    assert(count > 35);
    assert(count < 65);
}

///
@Name("containsFacet") pure @safe unittest
{
    auto sc = simplicialComplex([[1], [2,3], [3,4,5], [3,4,6]]);
    assert(sc.containsFacet([1]));
    assert(sc.containsFacet([2,3]));
    assert(sc.containsFacet([3,4,5]));
    assert(sc.containsFacet([3,4,6]));

    assert(!sc.containsFacet([2]));
    assert(!sc.containsFacet([4]));
    assert(!sc.containsFacet([3,6]));
    assert(!sc.containsFacet([1,7]));
    assert(!sc.containsFacet([1,6]));

    assert(simplicialComplex([[1,2]]).containsFacet([1,2]));
    assert(!simplicialComplex([[1,2]]).containsFacet([2]));
}

@Name("opEquals (pure @safe nothrow)") pure @safe unittest
{
    auto s1 = simplicialComplex([[1,2], [2,3,4]]);
    auto s2 = simplicialComplex([[1,3], [2,3,4]]);
    auto s3 = simplicialComplex([[2,3,4]]);
    auto s4 = simplicialComplex([[2,3,4], [1,2]]);
    auto s5 = simplicialComplex([[1,2], [2,3,4,5]]);

    () pure @safe nothrow {
        assert(s1 != s2);
        assert(s2 != s3);
        assert(s1 != s3);
        assert(s1 != s5);
        assert(s1 == s4);    
    }();
}

/*******************************************************************************
Get the f-vector of the simplicial complex. The returned array lists the
number of simplices in each dimension.
*/
size_t[] fVector(Vertex)(const ref SimplicialComplex!Vertex sc)
{   
    immutable maxDim = sc.facetVertices.byKey.maxElement.to!int;
    return iota(maxDim + 1).map!(dim => sc.simplices(dim).walkLength).array; 
}
///
@Name("fVector (pure @safe)") pure @safe unittest
{
    SimplicialComplex!() sc;
    sc.insertFacet([1,2]);
    assert(sc.fVector == [2,1]);        
}

/*******************************************************************************
Helper function template returning a newly constructed simplicial complex from
an array of facets (given as arrays of vertices.)
*/
auto simplicialComplex(F)(F initialFacets)
if (isInputRange!F && isInputRange!(ElementType!F))
{
    return SimplicialComplex!(Unqual!(ElementType!(ElementType!F)))(initialFacets);
}
///
@Name("simplicialComplex (pure @safe)") pure @safe unittest
{
    auto sc = simplicialComplex([[1,2], [2,3], [3,4,5], [6,7,8]]);
    sc.facets.shouldBeSameSetAs([[1,2], [2,3], [3,4,5], [6,7,8]]);

    assert(sc == simplicialComplex(sc.facets.array));

    const(int)[][] noFacets;
    assert(simplicialComplex(noFacets).facets.empty);
}

///
@Name("facets(dim) (pure nothrow @nogc @safe)") pure @safe unittest
{
    auto sc = simplicialComplex([[1,2], [2,4,5]]);
    int[2] edge = [1,2];
    int[3] triangle = [2,4,5];

    () pure nothrow @nogc @safe {
        auto facs1 = sc.facets(1);
        auto saved1 = facs1.save;

        auto facs2 = sc.facets(2);
        auto saved2 = facs2.save;

        assert(facs1.front == edge[]);
        facs1.popFront;
        assert(facs1.empty);    
        assert(!saved1.empty);    
        
        assert(facs2.front == triangle[]);
        facs2.popFront;
        assert(facs2.empty);    
        assert(!saved2.empty);    
    }();
}

///
@Name("facets (pure nothrow @nogc @safe)") pure @safe unittest
{
    auto sc = simplicialComplex([[1,2], [2,4,5]]);
    int[2] edge = [1,2];
    int[3] triangle = [2,4,5];

    static assert(isForwardRange!(typeof(sc.facets())));

    () pure nothrow @nogc @safe {
        auto facs = sc.facets;
        auto saved = facs.save;
        assert(facs.front == edge[]);
        facs.popFront;
        assert(facs.front == triangle[]);
        facs.popFront;
        assert(facs.empty);    
        assert(!saved.empty);    
    }();
}

///
@Name("star(range) (pure nothrow @nogc @safe)") pure @safe unittest
{
    auto sc = simplicialComplex([[1,2], [2,4,5]]);
    int[2] edge = [1,2];

    () pure nothrow @nogc @safe {
        auto starRange = sc.star(edge[]);
        auto saved = starRange.save;

        assert(starRange.front == edge[]);

        // Can still access the captured simplex
        assert(starRange.centerSimplex == edge[]);

        starRange.popFront;
        assert(starRange.empty);
        assert(!saved.empty);
    }();
}

///
@Name("link(range) (pure nothrow @nogc @safe)") pure @safe unittest
{
    auto sc = simplicialComplex([[1,2], [2,4,5]]);
    int[1] vertex = [1];

    () pure nothrow @nogc @safe {
        auto linkRange = sc.link(vertex[]);
        auto saved = linkRange.save;

        assert(linkRange.front.front == 2);

        linkRange.popFront;
        assert(linkRange.empty);
        assert(!saved.empty);
    }();  
}

///
@Name("contains (pure nothrow @nogc @safe)") pure @safe unittest
{
    auto sc = simplicialComplex([[1,2], [2,4,5]]);
    int[2] edge = [1,2];
    int[3] triangle = [2,4,5];
    int[4] tetrahedron = [1,4,5,7];

    () pure nothrow @nogc @safe {
        assert(sc.contains(edge[]));
        assert(sc.contains(triangle[]));
        assert(!sc.contains(tetrahedron[]));
    }();  
}

///
@Name("removeFacet (pure @safe)") pure @safe unittest
{
    auto sc = simplicialComplex([[1,2], [2,4,5]]);
    int[2] edge = [1,2];

    () pure @safe {
        sc.removeFacet(edge[]);
    }();  
}

///
@Name("numFacets (pure nothrow @nogc @safe)") pure @safe unittest
{
    auto sc = simplicialComplex([[1,2], [2,4,5]]);
    () pure nothrow @nogc @safe {
        assert(sc.numFacets == 2);
    }();  
}

private struct FacetRange(Vertex_ = int)
{
    const(SmallMap!(int, Vertex_[])) facetVertices;
    typeof(SmallMap!(int, Vertex_[])().byKey()) facetDims;
    const(Vertex_)[] vertices;

    this()(const(SmallMap!(int, Vertex_[])) facetVertices_)
    {
        facetVertices = facetVertices_;
        facetDims = facetVertices.byKey;
        if (!facetDims.empty)
        {
            vertices = facetVertices[facetDims.front];
        }
    }

    @property const(Vertex_)[] front() /* const */ pure nothrow @nogc @safe
    {
        assert(!this.empty);
        return vertices[0 .. facetDims.front + 1];
    }

    @property bool empty() const pure nothrow @nogc @safe
    {
        return vertices.empty;
    }

    void popFront() pure nothrow @nogc @safe
    {
        assert(!this.empty);
        vertices = vertices[facetDims.front + 1 .. $];
        if (vertices.empty)
        {
            facetDims.popFront;
            if (!facetDims.empty)
            {
                vertices = facetVertices[facetDims.front];
            }
        }
    }

    FacetRange!Vertex_ save() const pure nothrow @nogc @safe
    {
        return FacetRange!Vertex_(this.facetVertices);
    }
}

private struct StarRange(Vertex_ = int)
{
private:
    alias SimpComp = SimplicialComplex!Vertex_;
    alias Facets = typeof(SimpComp().facets()); 
    Facets facetsLeft;
public:
    const(Vertex_)[] centerSimplex;

    this()(const(Vertex_)[] centerSimplex_, const(Facets) facets_)
    {
        centerSimplex = centerSimplex_;
        facetsLeft = facets_.save.find!(f => centerSimplex.isSubsetOf(f));
    }

    // TO DO: facetsLeft.front isn't const, so this cant be const. FIX?
    @property const(Vertex_)[] front() /* const */ pure nothrow @nogc @safe
    {
        assert(!this.empty);
        return facetsLeft.front;
    }

    @property bool empty() const pure nothrow @nogc @safe
    {
        return facetsLeft.empty;
    }

    void popFront() pure nothrow @nogc @safe
    {
        assert(!this.empty);
        facetsLeft.popFront;
        while ((!facetsLeft.empty) && (!centerSimplex.isSubsetOf(facetsLeft.front)))
        {
            facetsLeft.popFront;
        }
    }

    StarRange!Vertex_ save() const pure nothrow @nogc @safe
    {
        return StarRange!Vertex_(this.centerSimplex, this.facetsLeft);
    }
}

@Name("StarRange (pure nothrow @nogc @safe)") pure @safe unittest
{
    auto sc = simplicialComplex([[1,2], [2,3,4], [2,3,5]]);
    int[] vertex = [3];
    int[] t1 = [2,3,4];
    int[] t2 = [2,3,5];

    int[][] answer1 = [t1, t2];
    int[][] answer2 = [t2, t1];

    static assert(isForwardRange!(StarRange!int));

    () pure nothrow @nogc @safe {
        auto linkRange = StarRange!int(vertex[], sc.facets);
        auto saved = linkRange.save;

        assert(linkRange.equal(answer1) || linkRange.equal(answer2));

        linkRange.popFront;
        assert(linkRange.front.equal(t1) || linkRange.front.equal(t2));

        linkRange.popFront;
        assert(linkRange.empty);
        assert(!saved.empty);
        assert(saved.equal(answer1) || saved.equal(answer2));
    }();
}

private struct LinkRange(Vertex_ = int)
{
private:
    StarRange!Vertex_ facetsLeft;
public:
    this()(const(StarRange!Vertex_) starFacets)
    {
        facetsLeft = starFacets.save;
    }

    // TO DO: facetsLeft.front isn't const, so this cant be const. FIX?
    @property auto front() /* const */ pure nothrow @nogc @safe
    {
        assert(!this.empty);
        return facetsLeft.front.setDifference(facetsLeft.centerSimplex);
    }

    @property bool empty() const pure nothrow @nogc @safe
    {
        return facetsLeft.empty;
    }

    void popFront() pure nothrow @nogc @safe
    {
        assert(!this.empty);
        facetsLeft.popFront;
    }

    LinkRange!Vertex_ save() const pure nothrow @nogc @safe
    {
        return LinkRange!Vertex_(this.facetsLeft);
    }
}

@Name("LinkRange (pure nothrow @nogc @safe)") pure @safe unittest
{
    auto sc = simplicialComplex([[1,2], [2,3,4], [2,3,5]]);
    int[] vertex = [3];
    int[][] ans1 = [[2, 4], [2, 5]];
    int[][] ans2 = [[2, 5], [2, 4]];

    static assert(isForwardRange!(LinkRange!int));

    () pure nothrow @nogc @safe {
        auto lnk = LinkRange!int(sc.star(vertex));
        auto saved = lnk.save;
        assert(lnk.equal!equal(ans1) || lnk.equal!equal(ans2));
        lnk.popFront;
        lnk.popFront;
        assert(lnk.empty);
        assert(!saved.empty);
        assert(saved.walkLength == 2);
        assert(saved.equal!equal(ans1) || saved.equal!equal(ans2));
    }();
}

/// Exception thrown when loading a simplicial complex from a file fails
class BadSimpCompLoad : Exception
{
    this(string msg, string file = __FILE__, size_t line = __LINE__) {
        super(msg, file, line);
    }
}

/******************************************************************************
Returns a simplicial complex loaded from the file specified by fileName.
*/
SimplicialComplex!Vertex loadSimplicialComplex(Vertex = int)(string fileName)
{
    auto facetString = File(fileName, "r").byLineCopy
        .filter!(line => line.front != '#').joiner.array;

    // accidentally have '#' or '╔' at end of the manifold facet line on some
    // files. Only include the part before any of these characters
    facetString = facetString.findSplit('#'.only)[0];
    facetString = facetString.findSplit('╔'.only)[0];

    Vertex[][] facets;
    try
    {
        facets = facetString.to!(Vertex[][]);
    }
    catch (Exception ex)
    {
        throw new BadSimpCompLoad("malformed facet list in file "
            ~ fileName ~ "\n    " ~ ex.msg);
    }
    facets.each!sort;
    return SimplicialComplex!Vertex(facets);
}

///
@Name("loadSimplicialComplex") @system unittest
{
    auto sc = loadSimplicialComplex(
        "data/manifold_sampler_unittest_load.dat");
    auto expected = simplicialComplex([[0, 1, 2], [0, 2, 3], [0, 3, 4], [0, 4, 5],
        [0, 5, 6], [0, 1, 6], [1, 2, 6], [2, 3, 5], [2, 5, 6], [3, 4, 5]]);
    assert(sc == expected);

    assertThrown(loadSimplicialComplex(
        "data/manifold_sampler_unittest_bad_load.dat"));
}

/******************************************************************************
* Saves a simplicial complex to a file specified by fileName.
*/
void saveTo(Vertex)(SimplicialComplex!Vertex sc, string fileName)
{
    auto saveFile = File(fileName, "w"); // Open in write-only mode
    saveFile.writeln("# created ", Clock.currTime.to!DateTime);
    saveFile.write(sc);
}

///
@Name("saveTo") unittest
{
    auto fileName = "data/simplicial_complex_unittest_save.dat";
    auto sc = simplicialComplex([[1,2],[3],[2,4,5], [2,4,6]]);
    sc.saveTo(fileName);
    auto loaded = loadSimplicialComplex(fileName);
    assert(loaded == sc);
}