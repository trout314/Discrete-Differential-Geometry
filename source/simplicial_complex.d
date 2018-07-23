import algorithms : connectedComponents, eulerCharacteristic,
    isCircle, isConnected,  isPureOfDim, join;
import fluent.asserts : should;
import std.algorithm : all, any, canFind, chunkBy, copy, countUntil, each,
    equal, filter, find, findAdjacent, isSorted, joiner, map, maxElement,
    setDifference, sort, sum, uniq;
import std.conv : to;
import std.range : array, chunks, dropExactly, ElementType, empty, enumerate, 
    front, iota, isForwardRange, isInputRange, isOutputRange, popFront, put, refRange, save, walkLength, zip;
import unit_threaded : Name;
import utility : isSubsetOf, SmallMap, StackArray, staticIota, subsets,
    subsetsOfSize, throwsWithMsg;

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
    assert(sc.facets.equal([[4,5], [5,6], [1,2,3], [2,3,4]]));
    assert(sc.numFacets == 4);

    // Can get a sorted array of all the facet simplices of a given dimension
    assert(sc.facets(0).empty);
    assert(sc.facets(1).equal([[4,5], [5,6]]));
    assert(sc.facets(2).equal([[1,2,3], [2,3,4]]));
    assert(sc.facets(3).empty);

    // Can get a sorted array of all the simplices of a given dimension
    assert(sc.simplices(0).equal([[1], [2], [3], [4], [5], [6]]));
    assert(sc.simplices(1).equal([[1,2], [1,3], [2,3], [2,4], [3,4], 
        [4,5], [5,6]]));
    assert(sc.simplices(2).equal(sc.facets(2)));
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
    assert(sc.star([5]).equal([[4,5], [5,6]]));

    assert(sc.star([4]).equal([[4,5], [2,3,4]]));
    assert(sc.star([2,3]).equal([[1,2,3], [2,3,4]]));

    // get link of a simplex as list of facets
    assert(sc.link([5]).equal!equal([[4], [6]]));
    assert(sc.link([4]).equal!equal([[5], [2,3]]));
    assert(sc.link([2,3]).equal!equal([[1], [4]]));

    // can print out a simplicial complex
    assert(sc.toString == "[[4, 5], [5, 6], [1, 2, 3], [2, 3, 4]]");

    // can construct from a range of ranges of vertices
    auto sc5 = SimplicialComplex!()(iota(3).map!(i => iota(3*i, 3*i + 3)));
    assert(sc5.facets.equal([[0, 1, 2], [3, 4, 5], [6, 7, 8]]));
}

///
@Name("toHash (@safe nothrow)") @safe unittest
{
    auto sc1 = simplicialComplex([[4,5], [1,2,3], [2,3,4], [5,6]]);
    auto sc2 = simplicialComplex([[4,5], [5,6], [2,3,4], [1,2,3]]);
    
    () nothrow @safe {
        assert(sc1.toHash == sc2.toHash);
    }();

    sc1.insertFacet([6,7]);

    // Note: this assert *could* fail, but useful check anyway!
    () nothrow @safe {
        assert(sc1.toHash != sc2.toHash);
    }();
}

///
@Name("works as AA KeyType") @safe unittest
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

@Name("additional tests") @safe unittest
{
    alias sComp = simplicialComplex;

    // ---------------- TEST A DEFAULT INITIALIZED COMPLEX ---------------------
    auto sc = SimplicialComplex!()();
    static assert(is(sc.Vertex == int)); // Default vertex type is int

    sc.facets.empty.should.equal(true);
    sc.numFacets.should.equal(0);
  
    // ---------------------------- ADD SOME FACETS ----------------------------
    sc.insertFacet([1, 2]);
    sc.insertFacet([2, 3]);
    sc.insertFacet([3, 4, 5]);
    
    sc.should.equal(sComp([[1,2], [2,3], [3,4,5]]));

    // -------------------------- TEST FUNCTIONALITY ---------------------------
    sc.facets.walkLength.should.equal(3);
    sc.numFacets.should.equal(3);

    sc.facets(0).empty.should.equal(true);
    sc.facets(1).should.containOnly([[1, 2], [2, 3]].to!(const(int)[][]));
    sc.facets(2).should.containOnly([[3, 4, 5]].to!(const(int)[][]));

    sc.simplices(0).should.containOnly([[1], [2], [3], [4], [5]]);
    sc.simplices(1).should.containOnly([[1, 2], [2, 3], [3, 4], [3, 5], [4, 5]]);
    sc.simplices(2).should.containOnly([[3, 4, 5]]);

    foreach(d; iota(3, 16))
    {
        sc.facets(d).empty.should.equal(true);
        sc.simplices(d).empty.should.equal(true);
    }
}

@Name("insertFacet") @safe unittest
{
    auto sc = SimplicialComplex!()();
    sc.insertFacet([1,2]);
    assert(sc.facets.equal([[1,2]]));
    sc.insertFacet([2,3,4]);
    assert(sc.facets.equal([[1,2], [2,3,4]]));
    sc.insertFacet([1,5,6]);
    sc.facets.array.should.containOnly([[1,2], [1,5,6], [2,3,4]]);
    sc.insertFacet([1,5,6,7]);
    assert(sc.facets.equal([[1,2], [2,3,4], [1,5,6,7]]));
}

@Name("insertFacet/removeFacet") @safe unittest
{
    auto sc = simplicialComplex([[1,2], [1,3], [1,4], [4,5,6]]);

    assert(sc.facets.equal([[1,2], [1,3], [1,4], [4,5,6]]));
    sc.removeFacet([1,3]);
    assert(sc.facets.equal([[1,2], [1,4], [4,5,6]]));
    sc.removeFacet([1,2]);
    assert(sc.facets.equal([[1,4], [4,5,6]]));
    sc.removeFacet([4,5,6]);
    assert(sc.facets.equal([[1,4]]));
    sc.removeFacet([1,4]);
    assert(sc.facets.empty);

    auto sc2 = simplicialComplex([[1,2,3], [1,2,4], [1,3,4], [2,3,4]]);
    sc2.removeFacet([1,2,3]);
    assert(sc2.facets.equal([[1,2,4], [1,3,4], [2,3,4]]));
    sc2.insertFacet([1,2,5]);
    sc2.insertFacet([1,3,5]);
    sc2.insertFacet([2,3,5]);
    sc2.facets.should.containOnly([[1,2,4], [1,2,5], [1,3,4], [1,3,5], [2,3,4],
        [2,3,5]]);

    SimplicialComplex!int sc3;
    sc3.insertFacet([1,2,3]);
    assert(sc3.facets.equal([[1,2,3]]));
    sc3.removeFacet([1,2,3]);
    assert(sc3.facets.empty);
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
struct SimplicialComplex(Vertex_ = int)
{
private:
    /* We store arrays containing all the vertices in the facets of a given 
    dimension. These arrays are flattened, so the array storing the 2-dim
    facets has the form [v1, v2, v3, w1, w2, w3, ... ] where [v1, v2, v3] give
    the first simplex, [w1, w2, w3] the second, etc.
    */
    SmallMap!(int, Vertex[]) facetVertices;

    /* indexOfFacet[f] = i means the vertices of f begins at facetVertices[i]
    */
    size_t[Vertex[]] indexOfFacet;
public:
    /***************************************************************************
    The type of the vertices in this simplicial complex.
    */
    alias Vertex = Vertex_;

    /***************************************************************************
    Construct a simplicial complex from a range of ranges with element type
    implicitly convertible to Vertex
    */
    this(R)(R facets) if (isInputRange!R && isInputRange!(ElementType!R) 
        && is(ElementType!(ElementType!R) : Vertex))
    {
        facets.each!(f => this.insertFacet(f));
    }

    // Postblit makes sure copy doesn't share data
    this(this) pure @safe
    {
        facetVertices = facetVertices.dup;
        facetVertices.byKey.each!(k => facetVertices[k] = facetVertices[k].dup);
    }

    /***************************************************************************
    Inserts a facet (given as an input range of vertices) into the simplicial
    complex.
    */
    void insertFacet(S)(S simplex) if (isInputRange!S && is(ElementType!S : Vertex))
    {
        // TO DO: Improve this!
        StackArray!(Vertex, 16) simplex_;
        simplex.each!(v => simplex_ ~= v);

        int dim = simplex_[].walkLength.to!int - 1;
        assert(dim >= 0);
        simplex_[].assertValidSimplex(dim);

        assert(!this.contains(simplex),
            "expected a simplex not already in the simplicial complex");

        /* We must remove any existing facets which are faces of the inserted
        facet. Also, we need independent copies of the facets to remove.
        TO DO: Create rawInsertFacet that skips this? NOTE: This is what
        is responsible for function allocating a closure!  */
        auto toRemove = simplex_[].subsets.map!array.filter!(
            s => this.facets.canFind(s)).array;
        toRemove.each!(s => this.removeFacet(s));

        if (dim in facetVertices)
        {
            facetVertices[dim] ~= simplex_;
        }
        else
        {
            facetVertices.insert(dim, simplex_);
        }
    }

    /***************************************************************************
    Removes an existing facet of the simplicial complex a facet (given as an 
    input range of vertices)
    */
    void removeFacet(S)(S simplex) if (isInputRange!S && is(ElementType!S : Vertex))
    {
        assert(this.contains(simplex),
            "tried to remove a facet not in the simplicial complex");
        
        assert(simplex.walkLength <= int.max);
        auto dim = cast(int) simplex.walkLength - 1;
        simplex.assertValidSimplex(dim);
     
        auto indx = facetVertices[dim].chunks(dim + 1).countUntil(simplex);
        
        copy(facetVertices[dim][(dim + 1) * (indx + 1)  .. $],
             facetVertices[dim][(dim + 1) * indx .. $ - (dim + 1)]);
        
        facetVertices[dim] = facetVertices[dim][0 .. $ - (dim + 1)];
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
    auto link(S)(S simplex) const if (isInputRange!S && is(ElementType!S : Vertex))
    {
        assert(this.contains(simplex), "expected a simplex in the simplicial complex");
        return LinkRange!Vertex(StarRange!Vertex(simplex, this.facets));
    }

    /***************************************************************************
    Returns the star of the given simplex a forward range of arrays of vertices
    giving the facets.
    */ 
    auto star(S)(S simplex) const if (isInputRange!S && is(ElementType!S : Vertex))
    {
        return StarRange!Vertex(simplex, this.facets);
    }

    /***************************************************************************
    Returns true if the given simplex is in this simplicial complex and false
    otherwise.
    */
    bool contains(S)(S simplex) const if (isInputRange!S && is(ElementType!S : Vertex))
    {
        return this.facets.any!(f => simplex.isSubsetOf(f));
    }

    /***************************************************************************
    Get the simplices of dimension `dim` as lists of vertices.
    */
    Vertex[][] simplices(int dim) const pure nothrow @safe
    {
        // TO DO: Reduce gc presure here. (Can't make @nogc I think.)

        assert(dim >= 0, "dimension must be non-negative, but got "
            ~ dim.to!string);

        Vertex[][] simplicesSeen;
        auto dims = facetVertices.byKey;
        foreach(d; dims.filter!(d => d >= dim))
        {
            foreach(f; this.facets(d))
            {
                simplicesSeen ~= f.subsetsOfSize(dim + 1)
                    .map!(s => s.array.dup).array;
            }
        }
        simplicesSeen = simplicesSeen.sort.uniq.array;
        return simplicesSeen;
    }

    bool opEquals()(auto ref const SimplicialComplex!Vertex sc) const
    {
        auto dims = sc.facetVertices.byKey;
        auto thisDims = this.facetVertices.byKey;
        if(dims.walkLength != thisDims.walkLength)
        {
            return false;
        }
        foreach(dim, thisDim; zip(dims, thisDims))
        {
            if(dim != thisDim)
            {
                return false;
            }
            const(Vertex)[][] f = sc.facets(dim).map!array.array.sort.array;
            const(Vertex)[][] thisF = this.facets(dim).map!array.array.sort.array;
            if(f != thisF)
            {
                return false;
            }
        }

        return true;
    }

    size_t toHash() const @safe nothrow
    {
        auto f = facets.map!array.array.sort.array;
        return () @trusted {return typeid(f).getHash(&f);} ();
    }

    /*******************************************************************************
    Returns a nice looking representation of the simplicial complex as a string.
    */
    string toString() const
    {
        return this.facets.to!string;
    }

    /*******************************************************************************
    Returns a string containing more detailed, implementation specific info on the
    simplicial complex
    */
    string toDetailedString() const
    {
        string output;
        output ~= "facets   : " ~ facets.to!string ~ "\n";
        output ~= "dims     : " ~ facetVertices.byKey.to!string ~ "\n";
        foreach(dim; facetVertices.byKey)
        {
            output ~= "vertex list (dim " ~ dim.to!string ~ "): "
                ~ facetVertices[dim].to!string ~ "\n";
        }
        
        return output;
    }

    void toString(W)(ref W writer) const if (isOutputRange!(W, char)) 
    {
        put(writer, facets.toString);
    }
}

@Name("opEquals (pure @safe)") pure @safe unittest
{
    auto s1 = simplicialComplex([[1,2], [2,3,4]]);
    auto s2 = simplicialComplex([[1,3], [2,3,4]]);
    auto s3 = simplicialComplex([[2,3,4]]);
    auto s4 = simplicialComplex([[2,3,4], [1,2]]);
    auto s5 = simplicialComplex([[1,2], [2,3,4,5]]);

    () pure @safe {
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
SimplicialComplex!Vertex simplicialComplex(Vertex)(const(Vertex)[][] initialFacets)
{
    return SimplicialComplex!Vertex(initialFacets);
}
///
@Name("simplicialComplex") pure @safe unittest
{
    auto sc = simplicialComplex([[1,2], [2,3], [3,4,5], [6,7,8]]);
    assert(sc.facets.equal([[1,2], [2,3], [3,4,5], [6,7,8]]));

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
@Name("removeFacet (pure nothrow @nogc @safe)") pure @safe unittest
{
    auto sc = simplicialComplex([[1,2], [2,4,5]]);
    int[2] edge = [1,2];

    () pure nothrow @nogc @safe {
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
        if(!facetDims.empty)
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
        if(vertices.empty)
        {
            facetDims.popFront;
            if(!facetDims.empty)
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
        while((!facetsLeft.empty) && (!centerSimplex.isSubsetOf(facetsLeft.front)))
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