import std.conv : to;
import std.container : make;
import std.traits : CommonType, isArray, isInstanceOf, isImplicitlyConvertible;
import std.algorithm : all, canFind, copy, filter, findAdjacent, isSorted,
    joiner, map, sort;
import std.range : chain, empty, front, popFront, walkLength, isInputRange, ElementType;
import std.array : array;

version (unittest)
{
    import std.stdio : writeln;
    import std.exception : assertThrown;
}

/*******************************************************************************
* Represents a non-degenerate simplex whose dimension (dim) is known at
* compile time. The simplex is represented as an array of vertices of
* user-specified type Vertex. Elements of type Vertex must be comparable with
* the less-than operator "<" and also comparible for equality with "==". We
* reserve the use of Vertex.init for vertices that have yet to be specified
* by the user.
*/
struct Simplex(size_t dim, Vertex = int)
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
   * and there can be no repeated vert!Rices.
   */
    this(R)(R vertices_)
    {
        static assert(isInputRange!R || isArray!R,
            "Simplex must be constructed from a range or array");

        static assert(is(ElementType!R : Vertex), "The element type ("
            ~ ElementType!R.stringof ~ ") of the range or array used to "
            ~ "constuct a Simplex must be implicitly convertible to the "
            ~ "vertex type (" ~ Vertex.stringof ~ ")");

        assert(vertices_.walkLength == dim + 1,
                "tried to create a simplex of dimension " ~ dim
                ~ " with too few vertices " ~ vertices_.to!string);

        copy(vertices_[].map!(to!Vertex), verts_[]);

        assert(vertices.isSorted,
                "tried to create a simplex " ~ this.toString ~ " with unsorted "
                ~ "vertices");
        assert(vertices.findAdjacent.length == 0,
                "tried to create a simplex " ~ this.toString ~ " containing a "
                ~ "repeated vertex");
        assert(!vertices.canFind(VertexType.init),
                "tried to create a simplex " ~ this.toString ~ " using a vertex"
                ~ " with value " ~ to!string(Vertex.init) ~ " which is reserved"
                ~ " for un-initialized vertices");
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
        return "[" ~ vertices.map!(to!string).joiner(",").to!string ~ "]";
    }

private:
    Vertex[dim + 1] verts_;
}
/// Some basic examples of creating and using the Simplex struct
unittest
{

    // create a simplex of dimension 1 with vertices [1,2]
    auto s1 = Simplex!1([1, 2]);
    pragma(msg, typeof(s1).stringof);

    // Need the proper number of vertices
    assertThrown!Error(Simplex!1([1, 2, 3]));
    assertThrown!Error(Simplex!1([7]));

    // The type used to store the vertices is accessible through the
    // symbol VertexType. This type is set to an int by default
    static assert(is(s1.VertexType == int));

    // The dimension of a simplex is also easily accesible
    static assert(s1.dimension == 1);

    // The user can specify the vertex type used
    auto s2 = Simplex!(1, ubyte)([ubyte(1), ubyte(2)]);
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
}

///
unittest
{

    // User specified vertex types must be comparable using < and ==
    // many types work already.
    static assert(__traits(compiles, Simplex!(5, string)));
    //    static assert(__traits(compiles, Simplex!(5, int[2]*)));
    auto s = Simplex!(5, int[2]*)();
    static assert(__traits(compiles, Simplex!(5, double[])));

    // User defined struct and class types will not work automatically
    struct A
    {
    }

    class B
    {
    }

    static assert(!__traits(compiles, Simplex!(5, A)));
    static assert(!__traits(compiles, Simplex!(3, B)));

    // Less than comparison using < is implemented by defining opCmp
    // TO DO: Explanation of why we need _two_ opCmp functions here?
    struct C
    {
        int opCmp(ref const C c) const
        {
            return 0;
        }

        int opCmp(const C c) const
        {
            return 0;
        }
    }
    //   static assert(__traits(compiles, Simplex!(5, C)));

    Simplex!(5, C) testC;

    // Equality comparison using == is implemented by defining opEquals, though
    // by default a struct defines opEquals
    struct D
    {
        int opCmp(ref const D d) const
        {
            return 0;
        }

        int opCmp(const D d) const
        {
            return 0;
        }

        @disable bool opEquals()(auto ref const D d) const;
    }

    static assert(!__traits(compiles, Simplex!(5, D)));

    class E
    {
        override int opCmp(Object e) const
        {
            return 0;
        }

        override bool opEquals(Object e) const
        {
            return false;
        }
    }

    auto x = new E;
    auto y = new E;
    //writeln(x < y, y < x, x == y);

    // TO DO: Fix this! Probably don't want to compare pointers via
    // the class opCmp right? comparison by actual pointers?
    //    auto z = Simplex!(5, E);
    //    static assert(__traits(compiles, Simplex!(5, E)));

//    static assert(simplex(1, 5, 9).toString == "[1,5,9]");
}

/******************************************************************************
* This helper template returns (by value) a newly constructed simplex from
* a list of vertices given as the arguments of simplex(...).
*/
auto simplex(V, size_t numVertices)(V[numVertices] vertices...)
{
    return Simplex!(numVertices - 1, V)(vertices);
}
///
unittest
{

//    auto s1 = simplex(1, 2, 4, 7);
//    static assert(is(s1.VertexType == int));

//    auto s2 = simplex(1L, 2L, 3L);
//    static assert(is(s2.VertexType == long));

    //    immutable i = 3, j = 5;
    //    auto s(1,2,4,7) = simplex(i, j);
    //    static assert(is(typeof(s(1,2,4,7)) == Simplex!(1, immutable int)));

//    auto s4 = simplex("alice", "bob", "carol");
//    static assert(is(s4.VertexType == string));

    // assertThrown!Error(simplex("bob", "alice", "carol"));
    // assertThrown!Error(simplex(0, 1, 2));
    // assertThrown!Error(simplex(1, 3, 3, 5));
    // static assert(!__traits(compiles, s4.hasFace(s(1,2,4,7))));
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

    return Simplex!(coDimension - 1, S.VertexType)(verticesInAnswer);
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

    // assert(s(1,2).oppositeFace(s(2)) == s(1));
    // assert(s(1,2).oppositeFace(s(1)) == s(2));
    // assert(s(1,2,4,7).oppositeFace(s(2)) == s(1,4,7));
    // assert(s(1,2,4,7).oppositeFace(s(1,4,7)) == s(2));
    // assert(s(1,2,4,7).oppositeFace(s(4,7)) == s(1,2));
    // assert(s(1,2,4,7).oppositeFace(s(1,2)) == s(4,7));

    // static assert(!__traits(compiles, s(2).oppositeFace(s(1))));
    // assertThrown!Error(s(1,2,4,7).oppositeFace(s(1,2,3)));
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

    // auto s0 = simplex(9);
    // auto s1 = simplex(5, 7);
    // auto s2 = simplex(1, 2, 3);

    // assert(s0.faces == [[9]]);
    // assert(s1.faces == [[5], [7], [5, 7]]);
    // assert(s2.faces == [[1], [2], [3], [1, 2], [1, 3], [2, 3], [1, 2, 3]]);
}

///
unittest
{

    // auto s1 = simplex(1, 2, 3, 5);
    // auto s2 = simplex(2, 3, 5);
    // auto s3 = simplex(1, 2, 3);
    // auto s4 = simplex(3, 5);
    // auto s5 = simplex(3);
    // auto s6 = simplex(3, 7);

    // assert(s1.dimension == 3);
    // assert(s5.dimension == 0);

    // assert(s1.hasFace(s1));
    // assert(hasFace(s1, s2));
    // assert(s1.hasFace(s3));
    // assert(s1.hasFace(s4));
    // assert(s1.hasFace(s5));
    // assert(s2.hasFace(s4));
    // assert(s4.hasFace(s5));
    // assert(!s2.hasFace(s1));
    // assert(!s2.hasFace(s3));
    // assert(!s3.hasFace(s2));
    // assert(!s1.hasFace(s6));
}

///
unittest
{

    // auto s1 = Simplex!(2, ubyte)([1, 2, 3]);
    // auto s2 = Simplex!(0, long)([257]);
    // auto s3 = Simplex!(1, ulong)([2, 3]);

    // assert(!s1.hasFace(s2));
    // assert(!s2.hasFace(s1));
    // assert(s1.hasFace(s3));
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