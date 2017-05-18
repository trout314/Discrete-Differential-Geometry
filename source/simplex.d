import std.conv : to;
import std.container : make;
import std.traits : CommonType, isArray, isInstanceOf, isImplicitlyConvertible;
import std.algorithm : all, canFind, copy, filter, findAdjacent, isSorted,
    joiner, map, sort;
import std.range : chain, empty, front, popFront, walkLength, isInputRange, ElementType;
import std.array : array;

version (unittest)
{
    import utility : throwsWithMsg;
    import std.stdio : writeln;
    import std.exception : assertThrown;
}

/*******************************************************************************
* Represents a non-degenerate simplex represented as set of vertices of
* user-specified type Vertex. Elements of type Vertex must be comparable with
* less-than "<" and also comparable for equality with "==".
*/
struct Simplex(int dim, Vertex = int)
{
    import utility : isEqualityComparable, isLessThanComparable;

    static assert(isLessThanComparable!Vertex,
            "Vertex type must support opCmp. See " ~
            "https://dlang.org/operatoroverloading.html#eqcmp");
    static assert(isEqualityComparable!Vertex,
            "Vertex type must support opEquals. See " ~
            "https://dlang.org/operatoroverloading.html#eqcmp");

    /***************************************************************************
   * a simplex can be constructed from a range with element type Vertex. There
   * must be exactly dim + 1 vertices in the range, they must occur in order,
   * and there can be no repeated vertices.
   */
    this(R)(R vertices_)
    {
        static assert(isInputRange!R || isArray!R,
            "Simplex must be constructed from a range or array");

        import utility : isConstructible;
        static assert(isConstructible!(ElementType!R, Vertex),
            "The vertex type (" ~ Vertex.stringof ~ ") must be constructible "
            ~ "from the element type (" ~ ElementType!R.stringof ~ ") of the "
            ~ "range or array used to construct a simplex");
        
        assert(vertices_.walkLength == dim + 1,
            "tried to create a simplex of dimension " ~ dim.stringof
            ~ " using the wrong number of vertices: " ~ vertices_.to!string);
        
        verts_[] = vertices_.map!(to!Vertex).array[];

        // We need to treat different types of vertice differently
        import std.traits : isPointer;
        static if (isPointer!Vertex)
        {
            // Probably don't want to sort by pointer values
            // TO DO: Make this optional. Sorting by address would be ok with
            // a custom allocator!
            alias sortingCriterion = (v1, v2) => *v1 < *v2;

            // Also don't want to check for equality by pointer values
            alias comparisonCriterion = (v1, v2) => *v1 == *v2;
        }
        else
        {
            // For non-pointer types we do the standard things
            alias sortingCriterion = (v1, v2) => v1 < v2;
            alias comparisonCriterion = (v1, v2) => v1 == v2;
        }

        static if (is(Vertex == class))
        {
            // Make sure to check for any null values
            assert(vertices.all!(v => v !is null),
                "class references used as vertices cannot be null");
        }

        static if (isPointer!Vertex)
        {
            assert(vertices.all!(v => v !is null),
                "pointers used as vertices cannot be null");
        }

        assert(verts_[].isSorted!sortingCriterion,
                "tried to create a simplex: " ~ this.toString ~ " with "
                ~"unsorted vertices");

        assert(verts_[].findAdjacent!comparisonCriterion.length == 0,
                "tried to create a simplex " ~ this.toString ~ " containing a "
                ~ "repeated vertex");
    }

    /***************************************************************************
    * a simplex can be copy constructed from any other simplex type S if S has 
    * the same dimension and has a vertex type (S.VertexType) implicitly
    * convertible to this simplex type's vertex type (Vertex).
    */
    this(S)(ref S simplexToCopy)
            if (isInstanceOf!(Simplex, S)
                && isImplicitlyConvertible!(S.VertexType, this.VertexType))
    {
        static assert(S.dimension == this.dimension);
        this.verts_ = simplexToCopy.verts_;
    }

    void opAssign(S)(S rhs)
            if (isInstanceOf!(Simplex, S)
                && isImplicitlyConvertible!(S.VertexType, this.VertexType))
    {
        static assert(S.dimension == this.dimension);
        this.verts_ = rhs.verts_;
    }

    /// dimension of the simplex. one less than the number of vertices.
    enum dimension = dim;

    /// the type of vertices used in the simplex
    alias VertexType = Vertex;

    /// provides read-only access to the vertices by returning a slice
    // TO DO: Think about the safety of this! Will "return ref" improve it?
    const(Vertex)[] vertices() const
    {
        return verts_[];
    }

    /// nice looking string representation
    string toString() const
    {
        return "[" ~ verts_[].map!(to!string).joiner(",").to!string ~ "]";
    }

private:
    Vertex[dim + 1] verts_;
}
/// Some basic examples
unittest
{
    // create a simplex of dimension 1 with vertices [1,2]
    auto s1 = Simplex!1([1, 2]);

    // Need the proper number of vertices
    assertThrown!Error(Simplex!1([7]));
    throwsWithMsg(
        Simplex!1([1,2,3]),
        "tried to create a simplex of dimension 1 using the wrong number of vertices: [1, 2, 3]"
    );


    // The type used to store the vertices is accessible through the
    // symbol VertexType. This type is set to an int by default
    static assert(is(s1.VertexType == int));

    // The dimension of a simplex is also easily accesible
    static assert(s1.dimension == 1);

    // The user can specify the vertex type used
    auto s2 = Simplex!(1, ubyte)([ubyte(3), ubyte(5)]);
    static assert(is(s2.VertexType == ubyte));

    // Simplices are value types, so s4 becomes an independent copy of s1
    auto s4 = s1;

    // Modifying s1 therefore has no effect on the copy s4.
    s1 = Simplex!1([5, 10]);
    assert(s1 != s4);
    assert(s4 == Simplex!1([1, 2]));

    // Of course, constant or immutable simplices can't be modified
    const s5 = s1;
    immutable s6 = s1;
    static assert(!__traits(compiles, s5 = s4));
    static assert(!__traits(compiles, s6 = s4));

    immutable(int)[] v1 = [1,2];
    immutable(int)[] v2 = [3,4];
    auto s = simplex(v1, v2);
}

///
unittest
{

    // User specified vertex types must be comparable using < and ==
    // many types work already.
    static assert(__traits(compiles, Simplex!(5, string)));
    static assert(__traits(compiles, Simplex!(5, int[2]*)));
    static assert(__traits(compiles, Simplex!(5, double[])));

    // User defined structs and classes will not work automatically
    struct A {}
    static assert(!__traits(compiles, Simplex!(3, A)));

    class B {}
    static assert(!__traits(compiles, Simplex!(5, B)));

    // To make a struct work we need to define opCmp and opEquals
    struct C
    {
        int label;
        int opCmp()(auto ref const C rhs) const
        {
            if (this.label < rhs.label)
            {
                return -1;
            }
            else if (this.label > rhs.label)
            {
                return 1;
            }
            else
            {
                return 0;
            }
        }
        // Note: the compiler-generated opEquals is OK here.
        // See $(LINK https://dlang.org/spec/operatoroverloading.html#compare).
    }

    auto testC = simplex(C(2), C(3));
    static assert(is(testC.VertexType == C));
    assertThrown!Error(simplex(C(1), C(0)));
    assertThrown!Error(simplex(C(2), C(2)));

    // We can use pointers to stucts and everything works
    // (the dereferencing is done internally)
    auto dynamicTestC = simplex(new C(2), new C(3));
    static assert(is(dynamicTestC.VertexType == C*));

    // null pointer vertices are not allowed
    throwsWithMsg(simplex(null, new C(0)),
        "pointers used as vertices cannot be null");

    // the usual things are disallowed too...
    assertThrown!Error(simplex(new C(1), new C(0)));
    assertThrown!Error(simplex(new C(2), new C(2)));

    // To make a class work we must define opCmp and opEquals as well
    class D
    {
        int label;

        this(int label_)
        {
            label = label_;
        }

        override int opCmp(Object rhs) const
        {
            if (this.label < rhs.to!D.label)
            {
                return -1;
            }
            else if (this.label > rhs.to!D.label)
            {
                return 1;
            }
            else
            {
                return 0;
            }
        }

        override bool opEquals(Object rhs) const
        {
            return this.label == rhs.to!D.label;
        }

        override string toString() const
        {
            import std.conv : to;
            return label.to!string;
        }
    }

    auto testD = simplex(new D(2), new D(3));
    
    static assert(is(testD.VertexType == D));
    assertThrown!Error(simplex(new D(1), new D(0)));
    assertThrown!Error(simplex(new D(2), new D(2)));

    // null class references are not allowed as vertices
    throwsWithMsg(simplex(new D(1), new D(3), null ),
        "class references used as vertices cannot be null");   
}

/******************************************************************************
* This helper template returns (by value) a newly constructed simplex from
* a list of vertices given as the arguments of simplex(...).
*/
auto simplex(size_t numVertices, V)(V[numVertices] vertices...)
{
    return Simplex!(numVertices - 1, V)(vertices[]);
}

///
unittest
{
    auto s1 = simplex(1, 2, 4, 7);
    static assert(is(s1.VertexType == int));

    auto s2 = simplex(1L, 2L, 3L);
    static assert(is(s2.VertexType == long));

    immutable i = 3, j = 5;
    auto s = simplex(i, j);
    static assert(is(typeof(s) == Simplex!(1, immutable int)));

   auto s4 = simplex("alice", "bob", "carol");
   static assert(is(s4.VertexType == string));

    assertThrown!Error(simplex("bob", "alice", "carol"));
    assertThrown!Error(simplex(1, 3, 3, 5));

    static assert(!__traits(compiles, s4.hasFace(simplex(1,2,4,7))));
}

/******************************************************************************
* Returns true if simplex has possibleFace as a face. Otherwise returns false
*/
bool hasFace(S, F)(const ref S simplex, const ref F possibleFace)
        if (isInstanceOf!(Simplex, S) && isInstanceOf!(Simplex, F))
{
    static assert(!is(CommonType!(F.VertexType, S.VertexType) == void),
            "hasFace called on simplices with vertices of incompatible types: "
            ~ "simplex.VertexType = " ~ S.VertexType.stringof ~ " and "
            ~ "possibleFace.VertexType = " ~ F.VertexType.stringof);

    alias CommonVertex = CommonType!(F.VertexType, S.VertexType);

    // Check using phobos algorithms
    bool phobosAnswer; // Can't put this in assert block. Why? no {} stripping?
    version (assert)
    {
        auto possFace = possibleFace.vertices.map!(to!CommonVertex);
        auto simpToCheck = simplex.vertices.map!(to!CommonVertex);
        phobosAnswer = possFace.all!(v => simpToCheck.canFind(v));
    }

    static if (S.dimension < F.dimension)
    {
        assert(phobosAnswer == false);
        return false;
    }
    else
    {
        auto simp = simplex.vertices;
        auto face = possibleFace.vertices;
        while ((simp.length > 0) && (face.length > 0))
        {
            if (face.front < simp.front)
            {
                assert(phobosAnswer == false);
                return false;
            }
            else if (face.front == simp.front)
            {
                face.popFront;
            }
            simp.popFront;
        }
        assert(phobosAnswer == (face.length == 0));
        return face.length == 0;
    }
}

/******************************************************************************
* Returns ...
*/
auto oppositeFace(S, F)(S simplex, F face)
        if (isInstanceOf!(Simplex, S) && isInstanceOf!(Simplex, F))
{
    static assert(F.dimension < S.dimension);
    enum coDimension = S.dimension - F.dimension;
    assert(simplex.hasFace(face));

    S.VertexType[coDimension] verticesInAnswer;
    size_t index = 0;
    foreach (vertex; simplex.vertices)
    {
        if (!face.vertices.canFind(vertex))
        {
            verticesInAnswer[index] = vertex;
            ++index;
        }
    }

    return Simplex!(coDimension - 1, S.VertexType)(verticesInAnswer[]);
}

/******************************************************************************
* Returns ...
*/
auto oppositeVertices(S, F)(const ref S simplex, const ref F face)
        if (isInstanceOf!(Simplex, S) && isInstanceOf!(Simplex, F))
{
}

///
unittest
{
    alias s = simplex;

    assert(s(1,2).oppositeFace(s(2)) == s(1));
    assert(s(1,2).oppositeFace(s(1)) == s(2));
    assert(s(1,2,4,7).oppositeFace(s(2)) == s(1,4,7));
    assert(s(1,2,4,7).oppositeFace(s(1,4,7)) == s(2));
    assert(s(1,2,4,7).oppositeFace(s(4,7)) == s(1,2));
    assert(s(1,2,4,7).oppositeFace(s(1,2)) == s(4,7));

    static assert(!__traits(compiles, s(2).oppositeFace(s(1))));
    assertThrown!Error(s(1,2,4,7).oppositeFace(s(1,2,3)));
}

/******************************************************************************
* Returns... 
*/
auto commonVertices(S)(const ref S simplex) if (isInstanceOf!(Simplex, S))
{
}

/******************************************************************************
* Returns a newly allocated dynamic array containing all faces of dimension...
*/
auto facesOfDim(size_t dim, S)(const ref S simplex)
if (isInstanceOf!(Simplex, S))
{
    return 4;
}

/******************************************************************************
* Returns a newly allocated dynamic array containing arrays of vertices. These
* are the vertices in all the faces of the given simplex (simplex). They are
* returned sorted by dimension (lowest to highest) and in dictionary order
* within each dimension.
*/
auto faces(S)(const ref S simplex) if (isInstanceOf!(Simplex, S))
{
    S.VertexType[][] faces;
    foreach (bitChoice; 1 .. 2 ^^ simplex.vertices.length)
    {
        S.VertexType[] face;
        foreach (index, vertex; simplex.vertices)
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
///
unittest
{
    auto s0 = simplex(9);
    auto s1 = simplex(5, 7);
    auto s2 = simplex(1, 2, 3);

    assert(s0.faces == [[9]]);
    assert(s1.faces == [[5], [7], [5, 7]]);
    assert(s2.faces == [[1], [2], [3], [1, 2], [1, 3], [2, 3], [1, 2, 3]]);
}

///
unittest
{

    auto s1 = simplex(1, 2, 3, 5);
    auto s2 = simplex(2, 3, 5);
    auto s3 = simplex(1, 2, 3);
    auto s4 = simplex(3, 5);
    auto s5 = simplex(3);
    auto s6 = simplex(3, 7);

    assert(s1.dimension == 3);
    assert(s5.dimension == 0);

    assert(s1.hasFace(s1));
    assert(hasFace(s1, s2));
    assert(s1.hasFace(s3));
    assert(s1.hasFace(s4));
    assert(s1.hasFace(s5));
    assert(s2.hasFace(s4));
    assert(s4.hasFace(s5));
    assert(!s2.hasFace(s1));
    assert(!s2.hasFace(s3));
    assert(!s3.hasFace(s2));
    assert(!s1.hasFace(s6));
}

///
unittest
{
    auto s1 = Simplex!(2, ubyte)([ubyte(1), ubyte(2), ubyte(3)]);
    auto s2 = Simplex!(0, long)([257L]);
    auto s3 = Simplex!(1, ulong)([2UL, 3UL]);

    assert(!s1.hasFace(s2));
    assert(!s2.hasFace(s1));
    assert(s1.hasFace(s3));
}

/******************************************************************************
* Returns ...
*/
auto hinges(S)(const ref S simplex) if (isInstanceOf!(Simplex, S))
{
}

/******************************************************************************
* Returns ...
*/
auto ridges(S)(ref S simplex) if (isInstanceOf!(Simplex, S))
{
    struct RidgeRange
    {
        S.VertexType[] vertices;
        size_t missing; // index of missing vertex plus one

        this(ref S simplex)
        {
            vertices = simplex.verts_[];
            missing = simplex.verts_.length;
        }

        @property bool empty() const
        {
            return missing == 0;
        }

        @property auto front()
        {
            return Simplex!(S.dimension - 1)(chain(vertices[0 .. missing - 1],
                    vertices[missing .. $]).array);
        }

        void popFront()
        {
            assert(!empty, "tried to popFront on empty range of ridges");
            --missing;
        }

        @property size_t length() const
        {
            return missing;
        }
    }

    return RidgeRange(simplex);
}