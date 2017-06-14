import std.algorithm : all, canFind, filter, find, map, sort;
import std.conv : to;
import std.exception : enforce;
import std.meta : AliasSeq, allSatisfy, anySatisfy;
import std.range : array, chain, enumerate, front, only, repeat, take, walkLength;
import std.traits : lvalueOf, rvalueOf;

/*******************************************************************************
Checks if items of type T can be compared with the less-than operation. Note 
that this means all other comparison operations are also valid, per dlang rules 
about operator overloading. Note also that the comparison must handle both
rvalues and lvalues. See 
$(LINK https://dlang.org/operatoroverloading.html#eqcmp)
*/
template isLessThanComparable(T)
{
    static if (is(T == class))
    {
        enum isLessThanComparable = __traits(isOverrideFunction, T.opCmp);
    }
    else
    {
        enum isLessThanComparable = is(typeof({bool b = T.init < T.init; }));
    }
}

/// Basic type tests
@safe pure nothrow @nogc unittest
{
    alias ComparableTypes = AliasSeq!(bool, byte, ubyte, short, ushort, int, 
        uint, long, ulong, float, double, real, ifloat, idouble, ireal, char,
        wchar, dchar, int*, void*);
    static assert(allSatisfy!(isLessThanComparable, ComparableTypes));
    
    alias NonComparableTypes = AliasSeq!(cfloat, cdouble, creal, void);
    static assert(!anySatisfy!(isLessThanComparable, NonComparableTypes));
}

/// Struct tests
@safe pure nothrow @nogc unittest
{
    struct A {}
    static assert(!isLessThanComparable!A);

    /* The "auto ref" template here allows us to have a single function that
    accepts both rvalues and lvalues. lvalues are taken by reference an rvalues
    are taken by copy. See
    https://dlang.org/spec/template.html#auto-ref-parameters. */
    struct B
    {
        int opCmp()(auto ref const B rhs) const
        {
            return 0; // b1 <= b2 and b2 <= b1 for all objects b1, b2
        }
    }
    static assert(isLessThanComparable!B);

    /* Making opCmp a normal non-template function taking the struct to compare
    by copy works but may result in unneccesary copying. */
    struct C
    {
        int opCmp(const C rhs) const
        {
            return 0;
        }
    }
    static assert(isLessThanComparable!C);

    /* Making opCmp a normal non-template function taking the strut to compare
    by reference does not work, since this opCmp cannot accept rvalues. */
    struct D
    {
        int opCmp(const ref D rhs) const
        {
            return 0;
        }
    }
    static assert(!isLessThanComparable!D);

    // Check that rvalues and lvalues work for types B, C
    import std.meta : AliasSeq;
    foreach(T; AliasSeq!(B, C))
    {
        T lvalue;
        assert(lvalue <= lvalue);
        assert(T.init <= T.init);
    }

    // Check that lvalues work but rvalues don't for type D
    {
        D lvalue;
        assert(lvalue <= lvalue);
        static assert(!__traits(compiles, D.init <= D.init));
    }

    /* TO DO: Find out why lvalueOf and rvalueOf don't seem to work
    correctly on A, B, C, D above. Bug? */
    foreach(T; AliasSeq!(A, B, C, D))
    {
        static assert(!__traits(compiles, lvalueOf!T <= lvalueOf!T));
        static assert(!__traits(compiles, rvalueOf!T <= rvalueOf!T));
    }
}

/// Class tests    
@safe pure nothrow @nogc unittest
{
    class A {}
    static assert(!isLessThanComparable!A);

    class B
    {
        override int opCmp(Object rhs) const
        {
            return 0;
        }
    }
    static assert(isLessThanComparable!B);
}

/*******************************************************************************
Checks if items of type T can be converted to a string. TO DO: Clean this up
and make sure it works as advertised!
*/
template isPrintable(T)
{
    static if(is(Vertex == class))
    {
        enum isPrintable = __traits(isOverrideFunction, T.toString);
    }
    else static if (is(Vertex == struct))
    {
        enum isPrintable = __traits(compiles, T.init.to!string);
    }
    else
    {
        enum isPrintable = __traits(compiles, T.init.to!string);
    }
}

/*******************************************************************************
Checks if items of type T can be compared for equality with ==. See
$(LINK https:dlang.org/operatoroverloading.html#eqcmp)
*/
template isEqualityComparable(T)
{
    static if (is(T == class))
    {
        enum isEqualityComparable = __traits(isOverrideFunction, T.opEquals);
    }
    else
    {
        enum isEqualityComparable = is(typeof({
            bool b = T.init == T.init;   
        }));
    }
}

/// Basic type tests
pure nothrow @nogc @safe unittest
{
    import std.meta : AliasSeq, allSatisfy, anySatisfy;

    alias ComparableTypes = AliasSeq!(bool, byte, ubyte, short, ushort, int, 
        uint, long, ulong, float, double, real, ifloat, idouble, ireal, char,
        wchar, dchar, cfloat, cdouble, creal,  int*, void*);
    static assert(allSatisfy!(isEqualityComparable, ComparableTypes));    
}

/// Struct Tests
pure nothrow @nogc @safe unittest
{
    struct A 
    {
        @disable bool opEquals()(auto ref const A rhs);
    }

    static assert(!isEqualityComparable!A);
 
    struct B
    {
        bool opEquals()(auto ref const B rhs) const
        {
            return true;
        }
    }
    
    static assert(isEqualityComparable!B);

    struct C
    {
        bool opEquals(const C rhs) const
        {
            return true;
        }
    }
    static assert(isEqualityComparable!C);

    struct D
    {
        bool opEquals(ref const D rhs) const
        {
            return true;
        }
    }
    static assert(!isEqualityComparable!D);
    
    // Check that rvalues and lvalues work for types B, C
    import std.meta : AliasSeq;
    foreach(T; AliasSeq!(B, C))
    {
        T lvalue;
        assert(lvalue == lvalue);
        assert(T.init == T.init);
    }

    // Check that lvalues work but rvalues don't for type D
    {
        D lvalue;
        assert(lvalue == lvalue);
        static assert(!__traits(compiles, D.init == D.init));
    }

    /* TO DO: Find out why lvalueOf and rvalueOf don't seem to work
    correctly on A, B, C, D above. Bug? */
    foreach(T; AliasSeq!(A, B, C, D))
    {
        import std.traits : lvalueOf, rvalueOf;
        static assert(!__traits(compiles, lvalueOf!T == lvalueOf!T));
        static assert(!__traits(compiles, rvalueOf!T == rvalueOf!T));
    }
}

/// Class tests
pure nothrow @nogc @safe unittest
{
    class A {}
    static assert(!isEqualityComparable!A);
 
    class B
    {
        override bool opEquals(Object rhs) const
        {
            return true;
        }
    }
    static assert(isEqualityComparable!B);
}


/*******************************************************************************
Checks if an instance of type S can be constructed from an instance of type T
TO DO: Didn't see any phobos trait for this, ask about it on forums.
TO DO: Make sure this actually works as advertised!
*/
enum bool isConstructible(From, To) = is(From : To) ||
is(typeof({
    auto t = To(From());
}));

pure nothrow @nogc @safe unittest
{
    // static assert(!isConstructible!(int, ubyte));
    // static assert(!isConstructible!(int, byte));
    // static assert(!isConstructible!(int, ushort));
    // static assert(!isConstructible!(int, short));
    // static assert(!isConstructible!(int, string));

    static assert(isConstructible!(ubyte, int));
    static assert(isConstructible!(byte, int));
    static assert(isConstructible!(ushort, int));
    static assert(isConstructible!(short, int));
    static assert(isConstructible!(string, string));
//    static assert(!isConstructible!(string, int));

    struct A {}
    struct B {}
    static assert(!isConstructible!(A, B));
    static assert(isConstructible!(A, A));

    struct C
    {
        this(A a) {}
    }
//    static assert(isConstructible!(A, C));
    auto c = C(A());

}

/*******************************************************************************
Asserts that `expression` throws a `ThrownType` which must be
a subtype of `Throwable`. All `Throwable`s are caught and will not escape
throwsWithMsg. If a type other than `ThrownType` is thrown, or if the message in
the thrown exception is not `msg` this funtion throws an AssertError. 

Params:
    ThrownType = The `Throwable` type to test for. Default is `Error`
    expression = The expression to test.
    msg        = Optional message to output on test failure.
    file       = The file where the error occurred.
                 Defaults to `__FILE__`.
    line       = The line where the error occurred.
                 Defaults to `__LINE__`.

    Throws:
        `AssertError` if the given `Throwable` with message `msg` is
        not thrown.
*/  
void throwsWithMsg(ThrownType : Throwable = Exception, E)
                  (lazy E expression,
                  string msg = null,
                  string file = __FILE__,
                  size_t line = __LINE__)
{
    try
    {
        expression();
    }
    catch(Exception exception)
    {
        auto thrown = cast(ThrownType) exception;
        if(thrown is null)
        {
            throw new Error("throwsWithMsg failed with wrong throwable "
                ~ "type.\n  Actual type  : " ~ typeid(exception).to!string 
                ~ "\n  Expected type: " ~ ThrownType.stringof, file, line);
        }
        else if (thrown.msg != msg)
        {
            throw new Error("throwsWithMsg failed with wrong message.\n"
                ~ "  Actual message  : " ~ thrown.msg ~ "\n  Expected message: "
                ~ msg, file, line);
        }
        else
        {
            return;
        }
    }
}

/*******************************************************************************
Returns all sequences of length `length` which contain `numOnes` ones and all 
other elements zero. Sequences given in increasing dictionary order. (This is 
the same as saying the corresponding binary numbers are in increasing order.)
*/
auto binarySequences(ulong length, ulong numOnes)
{
    assert(length >= 0);
    assert(numOnes >= 0);
    assert(length >= numOnes);

    if(numOnes == length)
    {
        return [1.repeat(numOnes).array];
    }

    if(numOnes == 0)
    {
        return [0.repeat(length).array];
    }

    return chain(
        binarySequences(length - 1, numOnes    ).map!(seq => [0] ~ seq),
        binarySequences(length - 1, numOnes - 1).map!(seq => [1] ~ seq)
    ).array;
}
///
@safe pure unittest
{
    assert(binarySequences(3, 0) == [[0,0,0]]);
    assert(binarySequences(3, 1) == [[0,0,1],[0,1,0],[1,0,0]]);
    assert(binarySequences(3, 2) == [[0,1,1],[1,0,1],[1,1,0]]);
    assert(binarySequences(3, 3) == [[1,1,1]]);

    assert(binarySequences(4, 2) == [[0,0,1,1],[0,1,0,1],[0,1,1,0],
                                     [1,0,0,1],[1,0,1,0],[1,1,0,0]]);
                                     
    assert(binarySequences(0, 0) == [[]]);

    /* Careful! The number of such sequences is "length choose numOnes" which
    can grow super-exponentially (depending on how numOnes is chosen). It is 
    easy to exhaust memory. For example, there are (50 choose 25) length 50
    sequences with 25 ones and 25 zeros, more than 126 trillion sequences! */

    // So the following would probably exhaust memory!
    // auto seqs = binarySequences(50, 25); 

    // Listing length 50 sequences with only a small numbers of ones is fine.
    assert(binarySequences(50, 0).length == 1);
    assert(binarySequences(50, 1).length == 50);
    assert(binarySequences(50, 2).length == (50 * 49) / 2);
}

/*******************************************************************************
An set class which uses an array to hold the underlying data. Should be very
fast for small data sets. (TO DO: Benchmarks!)
*/
struct SmallMap(KeyType, ValueType)
{
    private static struct Record
    {
        KeyType key;
        ValueType value;
    }

    /// insert a (key, value) pair into the map
    void insert(KeyType key, ValueType value)
    {
        enforce(key !in this, "key already present");
        data ~= Record(key, value);
        data.sort!((r1, r2) => r1.key < r2.key);
    }

    /// We support the (key in smallMap) syntax 
    bool opBinaryRight(string op : "in")(KeyType key)
    {
        return this.keys.canFind(key);
    }

    /// Get a lazy range returning the keys in increasing order sorted
    auto keys()
    {
        return data.map!(r => r.key);
    }

    /// Get a lazy range returning the values
    auto values()
    {
        return data.map!(r => r.value);
    }

    /// Provide index operator axes like this: smallMap[key]
    ref ValueType opIndex(KeyType key_)
    {
        auto found = data.find!(r => r.key == key_);

        if (found.length > 0)
        {
            return found.front.value;
        }
        else
        {
            this.insert(key_, ValueType.init);
            return data.find!(r => r.key == key_).front.value;
        }
    }

private:
    Record[] data;
}
///
pure @safe unittest
{
    SmallMap!(int, string) sm;
    sm.insert(5, "hello");
    assert(5 in sm);
    assert(2 !in sm);

    sm.insert(2, "goodbye");
    assert(2 in sm);

    sm[2] = "bubba";
    assert(sm[2] == "bubba");

    sm.insert(5, "nope").throwsWithMsg("key already present");
}

/*******************************************************************************
Returns true if set A is contained in set B
*/
auto isSubsetOf(A, B)(A setA, B setB)
{
    return setA.all!(element => setB.canFind(element));
}
///
unittest
{
    assert([1,3].isSubsetOf([1,3,4]));
}

/*******************************************************************************
Returns all the sub-ranges of length `sizeOfSubset` from the ordered range `set` 
*/
auto subsetsOfSize(R)(R set, int subsetSize)
{
    auto setSize = set.walkLength;

    assert(setSize >= 0);
    assert(subsetSize >= 0);
    assert(subsetSize <= setSize);

    auto elementChoices = binarySequences(setSize, setSize - subsetSize);
    return elementChoices.map!(
        choice => enumerate(choice)
                  .filter!(c => c.value == 0)
                  .map!(c => set[c.index])
                  .array
    ).array;

}
///
unittest
{
    import std.stdio : writeln;
    [1,2,3,4,5].subsetsOfSize(2).writeln;
}
