
import std.algorithm : all, canFind, copy, equal, filter, findAdjacent, isSorted, 
    joiner, map, setDifference, sort;
import std.conv : to;
import std.exception : assertThrown, enforce;
import std.format : format;
import std.meta : staticMap;
import std.range : array, chain, ElementType, empty, enumerate, front, iota, 
    isInputRange, popFront, walkLength;
import std.traits : CommonType, isArray, isImplicitlyConvertible, isInstanceOf, 
    isPointer, PointerTarget;
import utility : binarySequences, isConstructible, isEqualityComparable, 
    isLessThanComparable, isPrintable, isSubsetOf, subsetsOfSize, throwsWithMsg;

/*******************************************************************************
Represents a non-degenerate simplex represented as set of vertices of user
specified type Vertex. Elements of type Vertex must be comparable with less-than
"<" and also comparable for equality with "==". */
struct Simplex(int dim, Vertex = int)
{
    static assert(dim >= 0, "dimension cannot be negative");
    static assert(!isPointer!Vertex, "Vertex type cannot be a pointer type");
    static assert(!is(Vertex == union), "Vertex type cannot be a union type");
    static assert(isEqualityComparable!Vertex, "Vertex type must be equality "
        ~ "comparable");
    static assert(isLessThanComparable!Vertex, "Vertex type must be less-than "
        ~ "comparable");
    static assert(isPrintable!Vertex, "Vertex type must be printable");

    /***************************************************************************
    A simplex can be constructed from a range with element type Vertex. There
    must be exactly dim + 1 vertices in the range, they must occur in order,
    and there can be no repeated vertices. */
    this(R)(R vertices_) if (!isInstanceOf!(.Simplex, R))
    {
        static assert(isInputRange!R || isArray!R,
            "simplex must be constructed from a range or array");

        static assert(isConstructible!(ElementType!R, Vertex),
            "the vertex type: " ~ Vertex.stringof ~ " must be constructible "
            ~ "from the element type: " ~ ElementType!R.stringof ~ " of the "
            ~ "range or array used to construct a simplex");
        
        enforce(vertices_.walkLength == dim + 1,
            "a dimension " ~ dim.to!string ~ " simplex needs "
            ~ (dim+1).to!string ~ " vertices, but got vertices: " 
            ~ vertices_.to!string);
        
        copy(vertices_.map!(to!Vertex), verts_[]);

        static if (is(Vertex == class))
        {
            enforce(vertices.all!(v => v !is null), "null class references "
                ~ "cannot be vertices, but got vertices: " ~ this.toString);
        }

        enforce(verts_[].isSorted, "vertices must occur in increasing order, "
            ~ "but got vertices: " ~ this.toString);

        enforce(verts_[].findAdjacent.length == 0, "vertices must be distinct, "
            ~ "but got vertices: " ~ this.toString);
    }

    /***************************************************************************
    A simplex can be copy constructed from any other simplex type S if S has the
    same dimension and has a vertex type (S.VertexType) implicitly convertible 
    to this simplex type's vertex type (Vertex). TO DO: Make sure I understand
    why we need .Simplex in if-constraint instead of just Simplex */
    this(S)(auto ref const S simplexToCopy) if (isInstanceOf!(.Simplex, S))
    {
        static assert(S.dimension == this.dimension, "wrong dimension");   
        static assert(isImplicitlyConvertible!(S.VertexType, this.VertexType),
            "vertex types not implicitly convertible");
        copy(simplexToCopy.verts_[].map!(to!Vertex), verts_[]);
    }

    /***************************************************************************
    A simplex can be assigned to another simplex if they have the same dimension
    and the vertex type of vertex `rhs` is implicitly convertible to the vertex 
    type of the lhs vertex (here, `this`). */
    void opAssign(S)(auto ref const S rhs) if (isInstanceOf!(.Simplex, S))
    {
        static assert(isImplicitlyConvertible!(S.VertexType, this.VertexType));
        static assert(S.dimension == this.dimension);
        copy(rhs.verts_[].map!(to!Vertex), verts_[]);
    }

    /***************************************************************************
    The dimension of the simplex. Always one less than the number of vertices.*/
    enum dimension = dim;

    /***************************************************************************
    The type of vertices used in the simplex. */
    alias VertexType = Vertex;

    /***************************************************************************
    Provides read-only access to the vertices by returning a slice. TO DO: Think
    about DIP 1000 and how it relates to this. */
    const(Vertex)[] vertices() const pure nothrow @nogc @safe
    {
        return verts_[];
    }

    /***************************************************************************
    Returns a nice looking string representation. */
    string toString() const
    {
        return "[" ~ verts_[].map!(to!string).joiner(",").to!string ~ "]";
    }

private:
    Vertex[dim + 1] verts_;
}

/// Some basic examples
pure @safe unittest
{
    // Create a simplex of dimension 1 with vertices [1,2]
    auto s1 = Simplex!1([1, 2]);

    // Need the proper number of vertices
    throwsWithMsg(Simplex!1([1,2,3]),
        "a dimension 1 simplex needs 2 vertices, but got vertices: [1, 2, 3]");

    // Vertices must be in increasing order
    throwsWithMsg(Simplex!3([5,2,3,4]),
        "vertices must occur in increasing order, but got vertices: [5,2,3,4]");

    // Vertices must be distinct
    throwsWithMsg(Simplex!2([1,2,2]),
        "vertices must be distinct, but got vertices: [1,2,2]");

    /* The type used to store the vertices is accessible through the symbol
    VertexType. This type defaults to an int. */
    static assert(is(s1.VertexType == int));

    // The dimension of a simplex is also easily accesible
    static assert(s1.dimension == 1);

    // The user can specify the vertex type used
    auto sByte = Simplex!(2, byte)([3, 5, 7]);
    static assert(is(sByte.VertexType == byte));

    // Simplices are value types, so s2 becomes an independent copy of s1
    auto s2 = s1;

    // Modifying s1 therefore has no effect on the copy s2.
    s1 = Simplex!1([5, 10]);
    assert(s1 != s2);
    assert(s2 == Simplex!1([1, 2]));
    
    // Can use ranges to construct a simplex
    auto s3 = iota(5).to!(Simplex!4);
    assert(s3.vertices.array == [0,1,2,3,4]);
    static assert(is(s3.VertexType == int));

    auto s4 = Simplex!2(iota(3));
    assert(s4.vertices.array == [0,1,2]);
    assert(s3.hasFace(s4));

    auto s5 = s3.vertices.filter!(v => v > 2).to!(Simplex!1);
    assert(s5.vertices.array == [3,4]);

    /* Can explicitly construct one simplex from another if they have the same
    dimension and the vertex type of the origin simplex can be implicitly
    converted to the vertex type of the desitination simplex. */
    alias SI = Simplex!(1, int);
    alias SB = Simplex!(1, byte);

    // OK, byte implicitly converts to int
    auto s6 = SI(SB([1, 2])); 
    assert(s6 == SI([1,2]));

    // NOT OK, int will not convert to byte
    static assert(!__traits(compiles, SB(SI([1,2]))));
    
    // Here is an alternate way. TO DO: Why terrible errors here on wrong dim?
    auto s7 = SB([1, 2]).to!SI;
    assert(s7 == SI([1,2]));

    // Rules for assigning one simplex to another are same as for construction.
    auto s8 = SB([2,9]);
    auto s9 = SI([1,2]);
    s9 = s8;    // OK, same dimension, and byte implicitly convertes int

    // NOT OK, int doesn't implicitly convert to byte
    static assert(!__traits(compiles, s8 = s9));
}   

// additional tests
pure @safe unittest
{
    // Of course, constant or immutable simplices can't be modified
    const sCon = simplex(1,2,3);
    immutable sImm = simplex(4,5,6); 
    auto sMut = simplex(7,8,9);

    sMut = sCon; // OK
    sMut = sImm; // OK
    static assert(!__traits(compiles, sCon = sMut));
    static assert(!__traits(compiles, sImm = sMut));
    static assert(!__traits(compiles, sCon = sImm));
    static assert(!__traits(compiles, sImm = sCon));

    // Check that conversions work as expected
    auto sByte = simplex(byte(10), byte(11), byte(12));
    static assert(is(sByte.VertexType == byte));

    Simplex!(2, byte) sByte3 = sMut.vertices.to!(Simplex!(2, byte));

    auto sStr = simplex("1", "2", "3");
    auto sMut2 = simplex(13, 14, 15);
    sMut = sStr.vertices.to!(int[]).to!(Simplex!2);
    sMut2 = Simplex!(2, int)(sByte);

    assertThrown(Simplex!1([7]));
}

/// Shorter way to construct simplices
pure @safe unittest
{
    /* For convenience we have a function template "simplex" which can be used
    to costruct Simplex types. */
    auto s = simplex(1,2,3);
    static assert(is(typeof(s) == Simplex!(2, int)));
    assert(s.vertices.array == [1,2,3]);
}

/// Suitable Vertex types
unittest
{
    /* A user specified vertex types must be comparable using < and == and also
    be printable with writeln. Many types work already. */

    /* Arrays of comparable types like doubles work. The ordering used is the
    "dictionary order" on the arrays. */
    auto s1 = simplex([1.0, -3.4], [1.0, 1.5], [2.3, -1.6]);
    static assert(is(s1.VertexType == double[]));

    /* This would allows us to represent points in Euclidean plan for example.
    Although we would likely want to use static arrays. */
    double[2] v1 = [1.0, -3.4],
              v2 = [1.0,  1.5],
              v3 = [2.3, -1.6];
    auto s2 = simplex(v1, v2, v3);
    static assert(is(s2.VertexType == double[2]));

    // As arrays of comparable elements, strings work too.
    auto s3 = simplex("alice", "bob", "dinesh");

    // Pointer and union types not allowed
    static assert(!__traits(compiles, Simplex!(5, int*)));
    union U {}
    static assert(!__traits(compiles, Simplex!(2, U)));

    // User defined structs and classes will not work automatically
    struct A {}
    static assert(!__traits(compiles, Simplex!(3, A)));
    class B {}
    static assert(!__traits(compiles, Simplex!(5, B)));
    
    /* For a struct to work as a vertex type we must define appropriate opCmp
    and opEquals methods.
    
    Here is an example that works. Note that the compiler-generated opEquals
    which does member-by-member comparison is OK here. The rule is that opEquals
    and opCmp must be consistent with eachother. That is, for any c1 and c2 of 
    type C we must have: (c1 <= c2) && (c2 <= c1) if and only if (c1 == c2).

    See $(LINK https://dlang.org/spec/operatoroverloading.html#compare).*/
    static struct C
    {
        int label;

        // Here, opCmp just compares labels and returns the appropriate result
        int opCmp()(auto ref const C rhs) const
        {
            return (label < rhs.label) ? -1 : (label > rhs.label);
        }

        /* The compiler-generated toString would be fine, although the error 
        messages and unittests will be easier to read with the following. */ 
        string toString() const { return label.to!string; }
    }

    auto s4 = simplex(C(2), C(3));
    static assert(is(s4.VertexType == C));
    assertThrown(simplex(C(1), C(0)));
    assertThrown(simplex(C(2), C(2)));

    /* To make a class work as a vertex type we must define suitable overrides
    for opCmp, opEquals and toString. */ 
    class D
    {
        int label;

        this(int label_) { label = label_;}

        override int opCmp(Object rhs) const
        {
            return (label < rhs.to!D.label) ? -1 : (label > rhs.to!D.label);
        }

        override bool opEquals(Object rhs) const
        {
            return label == rhs.to!D.label;
        }

        override string toString() const { return label.to!string; }
    }

    auto s5 = simplex(new D(2), new D(3));
    static assert(is(s5.VertexType == D));
    simplex(new D(1), new D(0)).assertThrown;
    simplex(new D(2), new D(2)).assertThrown;

    // null class references are not allowed as vertices
    simplex(new D(3), null).throwsWithMsg(
        "null class references cannot be vertices, but got vertices: [3,null]");
}

/******************************************************************************
This helper template returns (by value) a newly constructed simplex from a list 
of vertices given as the arguments of simplex(...).
*/
auto simplex(size_t numVertices, V)(V[numVertices] vertices...)
{
    return Simplex!(numVertices - 1, V)(vertices[]);
}

///
pure @safe unittest
{
    auto s1 = simplex(1, 2, 4, 7);
    static assert(is(s1.VertexType == int));

    auto s2 = simplex(1L, 2L, 3L);
    static assert(is(s2.VertexType == long));

    immutable i = 3, j = 5;
//    auto s = simplex(i, j);
//    static assert(is(typeof(s) == Simplex!(1, immutable int)));

    auto s4 = simplex("alice", "bob", "carol");
    static assert(is(s4.VertexType == string));

    simplex("bob", "alice", "carol").assertThrown;
    simplex(1, 3, 3, 5).assertThrown;

    static assert(!__traits(compiles, s4.hasFace(simplex(1,2,4,7))));
}

/*******************************************************************************
Returns true if simplex has possibleFace as a face. Otherwise, returns false. */
bool hasFace(S, F)(auto ref const S simplex, auto ref const F possibleFace)
        if (isInstanceOf!(Simplex, S) && isInstanceOf!(Simplex, F))
{
    static assert(!is(CommonType!(F.VertexType, S.VertexType) == void),
            "hasFace called on simplices with vertices of incompatible types: "
            ~ "simplex.VertexType = " ~ S.VertexType.stringof ~ " and "
            ~ "possibleFace.VertexType = " ~ F.VertexType.stringof);
    alias CommonVertex = CommonType!(F.VertexType, S.VertexType);

    static if (S.dimension < F.dimension)
    {
        return false;
    }
    else
    {
        auto possibleFaceVertices = 
            possibleFace.vertices.map!(v => v.to!CommonVertex);
        auto simplexVertices = simplex.vertices.map!(v => v.to!CommonVertex);
        return possibleFaceVertices.isSubsetOf(simplexVertices);
    }
}

/******************************************************************************
Returns the face opposite `face` in `simplex`. */
auto oppositeFace(S, F)(S simplex, F face)
        if (isInstanceOf!(Simplex, S) && isInstanceOf!(Simplex, F))
{
    static assert(F.dimension < S.dimension);
    enforce(simplex.hasFace(face), "oppositeFace expected a face of " ~ 
        simplex.toString ~ " but got " ~ face.toString);

    return Simplex!(S.dimension - F.dimension - 1, S.VertexType)(
        setDifference(simplex.vertices, face.vertices));
}

///
pure @safe unittest
{
    alias s = simplex; // For clarity

    assert(s(1,2).oppositeFace(s(2)) == s(1));
    assert(s(1,2).oppositeFace(s(1)) == s(2));
    assert(s(1,2,4,7).oppositeFace(s(2)) == s(1,4,7));
    assert(s(1,2,4,7).oppositeFace(s(1,4,7)) == s(2));
    assert(s(1,2,4,7).oppositeFace(s(4,7)) == s(1,2));
    assert(s(1,2,4,7).oppositeFace(s(1,2)) == s(4,7));

    static assert(s(1,2,4,7).oppositeFace(s(1,4)) == s(2,7));

    static assert(!__traits(compiles, s(2).oppositeFace(s(1))));
    s(1,2,4,7).oppositeFace(s(1,2,3)).throwsWithMsg(
        "oppositeFace expected a face of [1,2,4,7] but got [1,2,3]");
}

/******************************************************************************
Returns array containing all faces of simplex `s` with dimension `dim`. Faces
given in dictionary order. */
auto facesOfDim(int dim, S)(auto ref const S s) if (isInstanceOf!(Simplex, S))
{
    static assert(dim <= s.dimension, "faces");

    alias SimplexType = Simplex!(dim, s.VertexType);
    return s.vertices.subsetsOfSize(dim + 1).map!(verts => SimplexType(verts));
}
///
pure @safe unittest
{
    alias s = simplex;   // For clarity

    assert(s(1,2,3).facesOfDim!0.array == [s(1), s(2), s(3)]);
    assert(s(1,2,3).facesOfDim!1.array == [s(1,2), s(1,3), s(2,3)]);
    assert(s(1,2,3).facesOfDim!2.array == [s(1,2,3)]);

    static assert(s(1,2,3).facesOfDim!0.array == [s(1), s(2), s(3)]);
    static assert(s(1,2,3).facesOfDim!1.array == [s(1,2), s(1,3), s(2,3)]);
    static assert(s(1,2,3).facesOfDim!2.array == [s(1,2,3)]);
}

/******************************************************************************
Returns a newly allocated dynamic array containing arrays of vertices. These are
the vertices in all the faces of the given simplex `s`. They are returned 
sorted by dimension (lowest to highest) and in dictionary order within each
dimension. */
auto faces(S)(auto ref const S s) if (isInstanceOf!(Simplex, S))
{
    alias makeArgs = () => iota(s.dimension + 1).map!(dim => 
        "s.facesOfDim!" ~ dim.to!string ~ ".map!(f => f.vertices.dup)")
        .joiner(", ").to!string;

    mixin("return chain(" ~ makeArgs() ~ ");");
}
///
pure @safe unittest
{
    assert(simplex(9).faces.array == [[9]]);
    assert(simplex(5,7).faces.array == [[5], [7], [5, 7]]);
    assert(simplex(1,2,3).faces.array == [[1], [2], [3], [1, 2], [1, 3], [2, 3],
        [1, 2, 3]]);
}

pure @safe unittest // NOTE: could have @nogc here if exceptions were @nogc
{
    auto pt = simplex(9);

    // TO DO: Can we make this @nogc?
    // assert(simplex(9).faces.front.front == 9);
}

///
pure @safe unittest
{
    alias s = simplex; // For clarity

    assert(s(1, 2, 3, 5).dimension == 3);
    assert(s(3).dimension == 0);

    assert(s(1, 2, 3, 5).hasFace(s(1, 2, 3, 5)));
    assert(s(1, 2, 3, 5).hasFace(s(2, 3, 5)));
    assert(s(1, 2, 3, 5).hasFace(s(1, 2, 3)));
    assert(s(1, 2, 3, 5).hasFace(s(3, 5)));
    assert(s(1, 2, 3, 5).hasFace(s(3)));
    assert(s(2, 3, 5).hasFace(s(3, 5)));
    assert(s(3, 5).hasFace(s(3)));
    assert(!s(2, 3, 5).hasFace(s(1, 2, 3, 5)));
    assert(!s(2, 3, 5).hasFace(s(1, 2, 3)));
    assert(!s(1, 2, 3).hasFace(s(2, 3, 5)));
    assert(!s(1, 2, 3, 5).hasFace(s(3,7)));
}

///
pure @safe unittest
{
    auto s1 = Simplex!(2, ubyte)([ubyte(1), ubyte(2), ubyte(3)]);
    auto s2 = Simplex!(0, long)([257L]);
    auto s3 = Simplex!(1, ulong)([2UL, 3UL]);

    assert(!s1.hasFace(s2));
    assert(!s2.hasFace(s1));
    assert(s1.hasFace(s3));
}

/******************************************************************************
Returns the codimension-2 faces (hinges) */
auto hinges(S)(auto ref const S simplex) if (isInstanceOf!(Simplex, S))
{
    static assert(simplex.dimension >= 2);
    return simplex.facesOfDim!(simplex.dimension - 2);
}
///
pure @safe unittest
{
    alias s = simplex;

    assert(s(1,2,3).hinges.array == [s(1), s(2), s(3)]);
    assert(s(1,3,5,7).hinges.array
        == [s(1,3), s(1,5), s(1,7), s(3,5), s(3,7), s(5,7)]);

    static assert(s(1,2,3).hinges.array == [s(1), s(2), s(3)]);
    static assert(s(1,3,5,7).hinges.array 
        == [s(1,3), s(1,5), s(1,7), s(3,5), s(3,7), s(5,7)]);
}

/******************************************************************************
Returns the codimension-1 faces (ridges) */
auto ridges(S)(auto ref S simplex) if (isInstanceOf!(Simplex, S))
{
    static assert(simplex.dimension >= 1);
    return simplex.facesOfDim!(simplex.dimension - 1);
}
///
pure @safe unittest
{
    alias s = simplex;

    assert(s(1, 2, 3).ridges.array == [s(1, 2), s(1, 3), s(2, 3)]);
    assert(s(2, 3, 5, 7).ridges.array == s(2, 3, 5, 7).facesOfDim!2.array);

    static assert(s(1,2 ).ridges.array == [s(1), s(2)]);
    static assert(s(2, 3, 5).ridges.array == s(2, 3, 5).facesOfDim!1.array);
}

auto join(S1, S2)(auto ref const S1 s1, auto ref const S2 s2)
    if (isInstanceOf!(Simplex, S1) && isInstanceOf!(Simplex, S2))
{
    auto commonVertices = chain(s1.vertices, s2.vertices).array.dup.sort();
    enforce(commonVertices.findAdjacent.empty, "join expected two simplices "
        ~ "without any common vertices, but got: " ~ s1.toString ~ " and " 
        ~ s2.toString);

    return Simplex!(s1.dimension + s2.dimension + 1)(commonVertices);
}
///
pure @safe unittest
{
    alias s = simplex;

    assert(join(s(1, 3), s(2, 4)) == s(1,2,3,4));
    assert(join(s(1, 3), s(2)) == s(1,2,3));

    throwsWithMsg(join(s(1, 2), s(1, 3)), "join expected two simplices without "
        ~"any common vertices, but got: [1,2] and [1,3]");
}