module utility;

import core.bitop : popcnt;
import std.algorithm : all, any, canFind, cartesianProduct, copy, each, equal, filter, find, findSplit, fold, joiner, map, merge,
    sort, startsWith, sum, uniq;  
import std.conv : to;
import std.exception : enforce, assertThrown;
import std.format : format;
import std.meta : AliasSeq, allSatisfy, anySatisfy;
import std.range : array, back, popBack, chain, cycle, drop, ElementType, empty, enumerate, front,
    iota, isForwardRange, isInputRange, join, popFront, repeat, retro, save, take,
    walkLength;
import std.stdio : File, writeln;
import std.string : strip;
import std.traits : lvalueOf, rvalueOf, hasAliasing, Unqual;


import unit_threaded : Name, shouldBeSameSetAs, shouldEqual;

import std.array : staticArray;

//dfmt off

// I wish we could compute these, but acos dosn't work at compile time. These
// come from wolfram alpha input: Table[N[2 Pi/ArcCos[1/k],30],{k, 2, 16}]
static immutable real[17] flatDegreeInDim = [
    real.nan,   // Not applicable in dimension zero!
    2.00000000000000000000000000000,
    6.00000000000000000000000000000,
    5.10429931211954041318017937918,
    4.76679212271567770016331231084,
    4.58814743301323607389345550261,
    4.47728161419381561316532718870,
    4.40168886795573189523776294354,
    4.34681580829256787810763853238,
    4.30515772121519709317292314208,
    4.27244785078781511448809296727,
    4.24607958792781091933915226732,
    4.22436998865935854222871451330,
    4.20618365430421015310353357902,
    4.19072666439811610044839625288,
    4.17742710626470998673691504053,
    4.16586250565979517934736897387];

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
        enum isLessThanComparable = is(typeof({ bool b = T.init < T.init; }));
    }
}

/// Basic type tests
@Name("isLessThanComparable (basic types)")
@safe pure nothrow @nogc unittest
{
    alias ComparableTypes = AliasSeq!(bool, byte, ubyte, short, ushort, int,
            uint, long, ulong, float, double, real,
            char, wchar, dchar, int*, void*);
    static assert(allSatisfy!(isLessThanComparable, ComparableTypes));

    alias NonComparableTypes = AliasSeq!(void);
    static assert(!anySatisfy!(isLessThanComparable, NonComparableTypes));
}

/// Struct tests
@Name("isLessThanComparable (structs)")
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

    foreach (T; AliasSeq!(B, C))
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
    foreach (T; AliasSeq!(A, B, C, D))
    {
        static assert(!__traits(compiles, lvalueOf!T <= lvalueOf!T));
        static assert(!__traits(compiles, rvalueOf!T <= rvalueOf!T));
    }
}

/// Class tests    
@Name("isLessThanComparable (classes)")
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
    static if (is(Vertex == class))
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
        enum isEqualityComparable = is(typeof({ bool b = T.init == T.init; }));
    }
}

/// Basic type tests
@Name("isEqualityComparable (basic types)")
pure nothrow @nogc @safe unittest
{
    import std.meta : AliasSeq, allSatisfy, anySatisfy;

    alias ComparableTypes = AliasSeq!(bool, byte, ubyte, short, ushort, int,
            uint, long, ulong, float, double, real,
            char, wchar, dchar, int*, void*);
    static assert(allSatisfy!(isEqualityComparable, ComparableTypes));
}

/// Struct Tests
@Name("isEqualityComparable (structs)")
pure nothrow @nogc @safe unittest
{
    struct A
    {
        @disable bool opEquals()(auto ref const A rhs) const;
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

    foreach (T; AliasSeq!(B, C))
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
    foreach (T; AliasSeq!(A, B, C, D))
    {
        import std.traits : lvalueOf, rvalueOf;

        static assert(!__traits(compiles, lvalueOf!T == lvalueOf!T));
        static assert(!__traits(compiles, rvalueOf!T == rvalueOf!T));
    }
}

/// Class tests
@Name("isEqualityComparable (classes)") pure nothrow @nogc @safe unittest
{
    class A
    {
    }

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
enum bool isConstructible(From, To) = is(From : To) || is(typeof({
            auto t = To(From());
        }));

@Name("isConstructible")
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
    static assert(!isConstructible!(string, int));

    struct A {}

    struct B {}

    static assert(!isConstructible!(A, B));
    static assert(isConstructible!(A, A));

    static struct C
    {
        this(A a) {}
    }

    /* Note, if we declare C as a non-static struct this assert fail. TO DO: Why
    does this happen? Something to do with these structs being inside a unittest
    and needing context pointers? */
    static assert(isConstructible!(A, C));
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
void throwsWithMsg(ThrownType : Throwable = Error, E)(lazy E expression,
        string msg = null, string file = __FILE__, size_t line = __LINE__)
{
    try
    {
        expression();
    }
    catch (Throwable exception)
    {
        auto thrown = cast(ThrownType) exception;
        if (thrown is null)
        {
            throw new Error("throwsWithMsg failed with wrong throwable type.\n  Actual type  : " 
                ~ typeid(exception).to!string ~ "\n  Expected type: " ~ ThrownType.stringof, file, line);
        }
        else if (thrown.msg != msg)
        {
            throw new Error("throwsWithMsg failed with wrong message.\n  Actual message  : " 
                ~ thrown.msg ~ "\n  Expected message: " ~ msg, file, line);
        }
        else
        {
            return;
        }
    }
    throw new Error("throwsWithMsg failed because expression did not throw", file, line);
}

///
@Name("throwsWithMsg") pure @system unittest
{
    static void throwException()
    {
        enforce(false, "kapow!");
    }
    throwException().throwsWithMsg!Exception("kapow!");

    static void throwError()
    {
        assert(false, "boom!");
    }

    throwError().throwsWithMsg("boom!");

    // throwsWithMsg throws an Error if the message is wrong
    throwError().throwsWithMsg("kaboom").throwsWithMsg(
        "throwsWithMsg failed with wrong message.\n"
        ~ "  Actual message  : boom!\n"
        ~ "  Expected message: kaboom");

    // throwsWithMsg throws an Error if the wrong Throwable type is given
    throwError().throwsWithMsg!Exception("boom!").throwsWithMsg(        
        "throwsWithMsg failed with wrong throwable type.\n"
        ~ "  Actual type  : core.exception.AssertError\n"
        ~ "  Expected type: Exception");

    static void doesNotThrow(){}

    // throwsWithMsg throws an Error if the expression does not throw
    doesNotThrow().throwsWithMsg("nope").throwsWithMsg(
        "throwsWithMsg failed because expression did not throw");
}

/*******************************************************************************
Returns all sequences of length `length` which contain `numOnes` ones and all 
other elements zero. Sequences given in increasing dictionary order. (This is 
the same as saying the corresponding binary numbers are in increasing order.)
*/
auto binarySequences(size_t length, size_t numOnes)
{
    assert(numOnes <= length, "the number of ones must be at most the length");

    if (numOnes == length)
    {
        return [1.repeat(numOnes).array];
    }

    if (numOnes == 0)
    {
        return [0.repeat(length).array];
    }

    return chain(binarySequences(length - 1, numOnes).map!(seq => [0] ~ seq),
            binarySequences(length - 1, numOnes - 1).map!(seq => [1] ~ seq)).array;
}
///
@Name("binarySequences") @safe pure nothrow unittest
{
    assert(binarySequences(3, 0) == [[0, 0, 0]]);
    assert(binarySequences(3, 1) == [[0, 0, 1], [0, 1, 0], [1, 0, 0]]);
    assert(binarySequences(3, 2) == [[0, 1, 1], [1, 0, 1], [1, 1, 0]]);
    assert(binarySequences(3, 3) == [[1, 1, 1]]);

    assert(binarySequences(4, 2) == [[0, 0, 1, 1], [0, 1, 0, 1], [0, 1, 1, 0],
            [1, 0, 0, 1], [1, 0, 1, 0], [1, 1, 0, 0]]);

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

///
@Name("binarySequences (errors)") pure @system unittest
{
    binarySequences(5,6).throwsWithMsg("the number of ones must be at most the length");
}

/*******************************************************************************
An set class which uses an array to hold the underlying data. Should be very
fast for small data sets. (TO DO: Benchmarks!)
*/
struct SmallMap(KeyType, ValueType)
{
private:
    Record[] data;
public:
    private static struct Record
    {
        KeyType key;
        ValueType value;
    }

    /// insert a (key, value) pair into the map
    void insert(const KeyType key, /* const */ ValueType value)
    {
        assert(key !in this, "key already present");
        data ~= Record(key, value);
        data.sort!((r1, r2) => r1.key < r2.key);
    }

    /// Get a new copy of this map
    auto dup()
    {
        return SmallMap!(KeyType, ValueType)(data.dup);
    }

    /// We support the (key in smallMap) syntax 
    bool opBinaryRight(string op : "in")(const KeyType key) const
    {
        return this.byKey.canFind(key);
    }

    /// Returns a range containing the keys in increasing order.
    auto byKey() const
    {
        return data.map!(r => r.key);
    }

    /// Returns a range containing the values in increasing order of their corresponding keys.
    auto byValue() const
    {
        return data.map!(r => r.value);
    }

    /// Provide index operator axes like this: smallMap[key]
    ref inout(ValueType) opIndex(const KeyType key_) inout
    {
        auto found = data.find!(r => r.key == key_);
        assert(found.length > 0, "SmallMap access error");
        return found.front.value;
    }
}
///
@Name("SmallMap") pure @safe unittest
{
    SmallMap!(int, string) sm;
    sm.insert(5, "hello");
    assert(5 in sm);
    assert(2 !in sm);

    sm.insert(2, "goodbye");
    assert(2 in sm);

    sm[2] = "Wayne";
    assert(sm[2] == "Wayne");

    assert(sm.byKey.equal([2, 5]));
    assert(sm.byValue.equal(["Wayne", "hello"]));

    // WARNING: A SmallMap instance behaves like a slice!

    // 1) If no insertions are done, we get reference-type behavior for accessing values
    auto sm1 = sm;
    sm1[5] = "all aboard!";
    assert(sm1.byKey.equal([2, 5]));
    assert(sm1.byValue.equal(["Wayne", "all aboard!"]));
    assert(sm.byKey.equal([2, 5]));
    assert(sm.byValue.equal(["Wayne", "all aboard!"]));

    // 2) After an insertion, two SmallMap instances become independent.
    sm1.insert(7,"good luck!");
    assert(sm1.byKey.equal([2, 5, 7]));
    assert(sm1.byValue.equal(["Wayne", "all aboard!", "good luck!"]));
    assert(sm.byKey.equal([2, 5]));
    assert(sm.byValue.equal(["Wayne", "all aboard!"]));
    sm1[2] = "Barb";
    assert(sm1.byValue.equal(["Barb", "all aboard!", "good luck!"]));
    assert(sm.byValue.equal(["Wayne", "all aboard!"]));
}

///
@Name("SmallMap (errors)") pure @system unittest
{
    SmallMap!(int, string) sm;
    sm.insert(5, "hello");
    sm.insert(5, "nope").throwsWithMsg("key already present");
}

///
@Name("SmallMap.keys (pure nothrow @nogc @safe)") pure @safe unittest
{
    SmallMap!(int, string) sm;
    sm.insert(1, "hello");
    () pure nothrow @nogc @safe{ assert(sm.byKey.front == 1); }();
}

///
@Name("SmallMap.values (pure nothrow @nogc @safe)") pure @safe unittest
{
    SmallMap!(int, string) sm;
    sm.insert(1, "hello");
    () pure nothrow @nogc @safe{ assert(sm.byValue.front == "hello"); }();
}

///
@Name("SmallMap.opIndex (pure nothrow @nogc @safe)") pure @safe unittest
{
    SmallMap!(int, string) sm;
    sm.insert(1, "hello");
    sm.insert(2, "goodbye");

    () pure nothrow @nogc @safe{ assert(sm[1] == "hello"); }();
}

/*******************************************************************************
Returns true if set A is contained in set B
*/
bool isSubsetOf(A, B)(A setA, B setB)
        if (isInputRange!A && isInputRange!B && is(ElementType!A : ElementType!B))
{
    return setA.all!(element => setB.canFind(element));
}
///
@Name("isSubsetOf") pure nothrow @safe unittest
{
    assert([1, 3].isSubsetOf([1, 3, 4]));
    assert(![2, 3].isSubsetOf([1, 3, 4]));

    int[] empty;
    assert(empty.isSubsetOf([6, 8]));
    assert(empty.isSubsetOf(empty));

    /* Note that using an empty string literal in the above example will not
    work since [] has type void[] */
    static assert(!__traits(compiles, [].isSubsetOf([6, 8])));
}
///
@Name("isSubsetOf (pure nothrow @nogc @safe)") pure nothrow @safe unittest
{
    assert([1, 2, 3].isSubsetOf([1, 2, 3, 4]));
    assert(2.iota.isSubsetOf(8.iota));

    () pure nothrow @nogc @safe{ assert(iota(3, 7).isSubsetOf(iota(10))); }();
}

/*******************************************************************************
Returns true if the uint `x` has a 1 at position `pos` and 0 otherwise (where
position 0 is the lesat significant bit and position 31 the highest.)
*/
bool hasOneAtBit(uint x, size_t pos) pure nothrow @nogc @safe
{
    assert(pos < uint.sizeof * 8, "bad bit position");
    return ((1 << pos) & x) > 0;
}
///
@Name("hasOneAtBit") pure nothrow @nogc @safe unittest
{
    assert(1.hasOneAtBit(0));
    assert(iota(1, 32).all!(pos => !1.hasOneAtBit(pos)));

    assert(iota(32).all!(pos => uint.max.hasOneAtBit(pos)));

    assert(!12.hasOneAtBit(0));
    assert(!12.hasOneAtBit(1));
    assert(12.hasOneAtBit(2));
    assert(12.hasOneAtBit(3));
    assert(iota(4, 32).all!(pos => !12.hasOneAtBit(pos)));
}

///
@Name("hasOneAtBit (errors)") pure @system unittest
{
    // An error is thrown if the bit position is outside the uint
    7.hasOneAtBit(32).throwsWithMsg("bad bit position");
}

private auto subsetFromUint(R)(R set, uint whichToKeep) if (isInputRange!R)
{
    assert(whichToKeep < (1 << set.walkLength),
            "1 bits found in positions not corresponsing to elements in set");

    static struct SubsetFromUintRange
    {
        R set_;
        uint whichToKeep_;

        void advance()()
        {
            while (!this.empty && !whichToKeep_.hasOneAtBit(0))
            {
                assert(!set_.empty);
                set_.popFront;
                whichToKeep_ >>= 1;
            }
        }

        this(R s, uint w)
        {
            set_ = s;
            whichToKeep_ = w;
            this.advance;
        }

        auto front()
        {
            assert(!this.empty,
                "called popFront on an empty SubsetFromUintRange");
            assert(!set_.empty);
            return set_.front;
        }

        auto popFront()
        {
            assert(!this.empty,
                "called popFront on an empty SubsetFromUintRange");
            assert(!set_.empty);
            set_.popFront;
            whichToKeep_ >>= 1;
            this.advance;
        }

        auto empty()
        {
            return whichToKeep_ == 0;
        }

        static if (isForwardRange!R)
        {
            auto save()
            {
                return SubsetFromUintRange(set_, whichToKeep_);
            }
        }
    }

    return SubsetFromUintRange(set, whichToKeep);
}
///
@Name("subsetsFromUint") pure nothrow @safe unittest
{
    assert([3, 4, 5].subsetFromUint(0).empty);
    assert([3, 4, 5].subsetFromUint(1).array == [3]);
    assert([3, 4, 5].subsetFromUint(2).array == [4]);
    assert([3, 4, 5].subsetFromUint(3).array == [3, 4]);
    assert([3, 4, 5].subsetFromUint(4).array == [5]);
    assert([3, 4, 5].subsetFromUint(5).array == [3, 5]);
    assert([3, 4, 5].subsetFromUint(6).array == [4, 5]);
    assert([3, 4, 5].subsetFromUint(7).array == [3, 4, 5]);
}

///
@Name("subsetsFromUint (errors)") pure nothrow @system unittest
{
    [3, 4, 5].subsetFromUint(1 << 3).throwsWithMsg(
            "1 bits found in positions not corresponsing to elements in set");
}

///
@Name("subsetsFromUint (pure nothrow @nogc @safe)") pure nothrow @nogc @safe unittest
{
    auto r = iota(1, 5).subsetFromUint(7);
    assert(r.front == 1);
    r.popFront;
    assert(r.front == 2);
    r.popFront;
    assert(r.front == 3);
    r.popFront;
    assert(r.empty);
}

/*******************************************************************************
*/
auto subsetsOfSize(R)(R set, int subsetSize) if (isInputRange!R)
{
    immutable setSize = set.walkLength;

    assert(subsetSize >= 0, "subset size must be non-negative");
    assert(subsetSize <= setSize, "subset size must be at most the size of the set");
    assert(setSize <= 31, "subset size must be at most 31");

    static struct SubsetsOfSizeRange
    {
    private:
        /* We will think of the bits in whichToKeep as indexed starting with
        zero for the least significant bit. Highest bit is empty flag. */
        uint whichToKeep;
        R set_; // Underlying set from which to draw elements
    public:
        auto front()
        {
            assert(!this.empty);
            return set_.subsetFromUint(whichToKeep);
        }

        auto popFront()
        {
            assert(!this.empty);

            immutable len = set_.walkLength;
            immutable subLen = whichToKeep.popcnt;
            immutable uint lastToKeep = ((1 << subLen) - 1) << (len - subLen);

            // Check for impending emptiness
            if (whichToKeep == lastToKeep)
            {
                whichToKeep = 1 << 31;
                return;
            }

            auto currentPos = len - 1;
            if (whichToKeep.hasOneAtBit(currentPos))
            {
                // Skip over additional ones
                while (whichToKeep.hasOneAtBit(currentPos))
                {
                    --currentPos;
                }
                assert(currentPos > 0 && currentPos < len);
                assert(!whichToKeep.hasOneAtBit(currentPos));

                immutable numOnesSeen = len - currentPos - 1;

                // Find another one to move
                while (!whichToKeep.hasOneAtBit(currentPos))
                {
                    --currentPos;
                }
                assert(currentPos >= 0 && currentPos < len);
                assert(whichToKeep.hasOneAtBit(currentPos));

                whichToKeep &= ~(1 << currentPos);

                iota(currentPos + 1, currentPos + 2 + numOnesSeen)
                    .each!(pos => whichToKeep |= (1 << pos));

                iota(currentPos + 2 + numOnesSeen, len)
                    .each!(pos => whichToKeep &= ~(1 << pos));
            }
            else // Zero at current position
            {
                assert(!whichToKeep.hasOneAtBit(currentPos));

                // Find a one to move
                while (!whichToKeep.hasOneAtBit(currentPos))
                {
                    --currentPos;
                }
                assert(currentPos >= 0 && currentPos < len);
                assert(whichToKeep.hasOneAtBit(currentPos));

                whichToKeep &= ~(1 << currentPos);
                whichToKeep |= (1 << (currentPos + 1));
            }
        }

        auto empty()
        {
            return whichToKeep.hasOneAtBit(31);
        }

        static if (isForwardRange!R)
        {
            auto save()
            {
                return SubsetsOfSizeRange(whichToKeep, set_.save);
            }
        }
    }

    static assert(isInputRange!SubsetsOfSizeRange);

    return SubsetsOfSizeRange((1 << subsetSize) - 1, set);
}
///
@Name("subsetsOfSize") pure @safe unittest
{

    [1, 2, 3, 4].subsetsOfSize(1).map!array.shouldBeSameSetAs([[1], [2], [3], [4]]);

    [1, 2, 3, 4].subsetsOfSize(2).map!array.shouldBeSameSetAs([[1, 2], [1, 3],
            [1, 4], [2, 3], [2, 4], [3, 4]]);

    [1, 2, 3, 4].subsetsOfSize(3).map!array.shouldBeSameSetAs([[1, 2, 3], [1, 2, 4],
        [1, 3, 4], [2, 3, 4]]);

    [1, 2, 3, 4].subsetsOfSize(4).map!array.shouldBeSameSetAs([[1, 2, 3, 4]]);

    [1].subsetsOfSize(1).map!array.shouldBeSameSetAs([[1]]);

    0.iota.subsetsOfSize(0).map!array.shouldBeSameSetAs([[]]);
    assert(iota(3).subsetsOfSize(0).front.empty);

    3.iota.subsetsOfSize(1).map!array.shouldBeSameSetAs([[0], [1], [2]]);
    3.iota.subsetsOfSize(2).map!array.shouldBeSameSetAs([[0, 1], [0, 2], [1, 2]]);
    3.iota.subsetsOfSize(3).map!array.shouldBeSameSetAs([[0, 1, 2]]);

    iota(31).subsetsOfSize(1).walkLength.shouldEqual(31);
    iota(10).subsetsOfSize(2).walkLength.shouldEqual(45);
}

///
@Name("subsetsOfSize (errors)") pure nothrow @system unittest
{
    [1,2].subsetsOfSize(-1).throwsWithMsg("subset size must be non-negative");
    3.iota.subsetsOfSize(4).throwsWithMsg("subset size must be at most the size of the set");
    40.iota.subsetsOfSize(32).throwsWithMsg("subset size must be at most 31");
}

///
@Name("subsetsOfSize (pure nothrow @nogc @safe)") pure nothrow @safe unittest
{
    int[4] set = [1, 2, 3, 4];
    int[2][6] expectedSubsets = [[1, 2], [1, 3], [1, 4], [2, 3], [2, 4], [3, 4]];
    int[2][6] actualSubsets;

    () pure nothrow @nogc @safe{
        auto subsetsRange = set[].subsetsOfSize(2);
        auto saved = subsetsRange.save;

        foreach (indx; 0 .. 6)
        {
            copy(subsetsRange.front, actualSubsets[indx][]);
            subsetsRange.popFront;
        }
        assert(subsetsRange.empty);
        assert(!saved.empty);

        // Don't want to test the specific order, so sort
        actualSubsets[].sort();

        foreach (indx; 0 .. 6)
        {
            assert(actualSubsets[indx] == expectedSubsets[indx]);
        }

        // Some additional tests
        assert(31.iota.subsetsOfSize(1).walkLength == 31);
        assert(31.iota.subsetsOfSize(31).walkLength == 1);
        assert(10.iota.subsetsOfSize(2).walkLength == 45);
        assert(1.iota.subsetsOfSize(1).walkLength == 1);
    }();
}

/******************************************************************************
Returns a lazy range that lists the subsets of the given input range.
*/
auto subsets(R)(R set) if (isInputRange!R)
{
    assert(set.walkLength <= 31);

    static struct SubsetsRange
    {
    private:
        R set_;
        uint whichToKeep_;
    public:
        auto front()
        {
            assert(!empty);
            return subsetFromUint(set_, whichToKeep_);
        }

        auto popFront()
        {
            assert(!empty);
            ++whichToKeep_;
        }

        auto empty()
        {
            return whichToKeep_ == (1 << set_.walkLength);
        }

        static if (isForwardRange!R)
        {
            auto save()
            {
                return SubsetsRange(set_, whichToKeep_);
            }
        }
    }

    return SubsetsRange(set, 1);
}
///
@Name("subsets") pure @safe unittest
{
    [1, 2, 3].subsets.map!array.shouldBeSameSetAs([[1], [2], [3], [1, 2],
            [1, 3], [2, 3], [1, 2, 3]]);

    iota(1, 5).subsets.map!array.shouldBeSameSetAs([[1], [2], [3], [4], [1, 2],
            [1, 3], [1, 4], [2, 3], [2, 4], [3, 4], [1, 2, 3], [1, 2, 4],
            [1, 3, 4], [2, 3, 4], [1, 2, 3, 4]]);

    int[] emptySet;
    assert(emptySet.subsets.empty);

    assert(iota(6).subsets.walkLength == 63);
}
///
@Name("subsets (pure nothrow @nogc @safe)") pure nothrow @nogc @safe unittest
{
    auto subsetRange = 3.iota.subsets;
    auto saved = subsetRange.save;
    assert(subsetRange.front.front == 0);
    subsetRange.popFront;
    assert(subsetRange.front.front == 1);
    assert(saved.front.front == 0);
    assert(!subsetRange.empty);

    // TO DO: Improve this test. Ugly...
}

@Name("subsets bug") pure @safe unittest
{
    auto set = [1,2,3];
   
    foreach(s; set.subsets)
    {
        assert(s.walkLength == s.array.walkLength);        
    }
}

/*******************************************************************************
Gives AliasSeq containing integers `begin`, `begin + 1`, ... ,`end - 1`
*/
template staticIota(int begin, int end)
{
    static if (begin + 1 >= end)
    {
        static if (begin >= end)
        {
            alias staticIota = AliasSeq!();
        }
        else
        {
            alias staticIota = AliasSeq!(begin);
        }
    }
    else
    {
        enum middle = begin + (end - begin) / 2;
        alias staticIota = AliasSeq!(staticIota!(begin, middle), staticIota!(middle, end));
    }
}
///
@Name("staticIota") pure nothrow @nogc @safe unittest
{
    alias seq = staticIota!(3, 6);
    static assert(seq[0] == 3);
    static assert(seq[1] == 4);
    static assert(seq[2] == 5);

    alias seq2 = staticIota!(1, 5);
    alias func = (a, b, c, d) => a + b + c + d;
    static assert(func(seq2) == 10);
}


/***************************************************************************
TO DO: docs...
*/
auto productUnion(Vertex = int, R1, R2)(R1 range1, R2 range2)
if (isInputRangeOfInputRangeOf!(R1, const(Vertex)) 
    && isInputRangeOfInputRangeOf!(R2, const(Vertex))) 
{
    return cartesianProduct(range1, range2).map!(p => merge(p[0], p[1]));
}
///
@Name("productUnion") pure @safe unittest
{
    auto f1 = [[1,2],[3]];
    auto f2 = [[4,5,6], [7,8], [9]];
    productUnion(f1, f2).map!array.shouldBeSameSetAs(
        [[1, 2, 4, 5, 6], [1, 2, 7, 8], [1, 2, 9],
         [3, 4, 5, 6], [3, 7, 8], [3, 9]]);
    
    // TO DO: More tests...
}

/*******************************************************************************
auto capture(Range, Data...)(Range range, Data data) if (isInputRange!Range)
{
    alias dataDeclarations = 
        (len) => len.iota.map!(indx => "Data[%s] d%s;".format(indx, indx)).joiner("\n").array;
    alias dataArgs =
        (len) => len.iota.map!(indx => "d%s".format(indx)).joiner(",").array;

    static struct CaptureRange
    {
    private:
        Range range_;
        mixin(dataDeclarations(data.length));
    public:
        auto front()
        {
            static struct CaptureFront
            {
                typeof(range_.front) result;
                mixin(dataDeclarations(data.length));
                alias result this;
            }

            mixin("return CaptureFront(range_.front, " ~ dataArgs(data.length) ~ ");");
        }

        auto popFront()
        {
            range_.popFront;
        }

        auto empty()
        {
            return range_.empty;
        }

        static if (isForwardRange!Range)
        {
            auto save()
            {
                mixin("return CaptureRange(range_, " ~ dataArgs(data.length) ~ ");");
            }
        }

        // TO DO: Support random access ranges, and bi-directional ranges
    }

    return CaptureRange(range, data);
}
///
@Name("capture (pure nothrow @nogc @safe)") pure nothrow @nogc @safe unittest
{
    auto r = 10.iota.capture(1, "hello");
    assert(r.front == 0);
    auto s = r.save;
    r.popFront;
    assert(r.front == 1);
    assert(s.front == 0);

    // Access captured data from each element of the range
    assert(r.front.d0 == 1);
    assert(r.front.d1 == "hello");

    // TO DO: More tests! This function will probably get lots of use!
}
*/

/*******************************************************************************
*/
struct StackArray(T, size_t maxLength)
{
private:
    T[maxLength] data;
    size_t currentLength;
public:
    /// Returns the currently used length
    auto length() const
    {
        return currentLength;
    }

    /// Sets length, additional elements are .init initialized
    auto length(size_t newLength)
    {
        assert(newLength <= maxLength, "new length must be at most max length");

        if (newLength < currentLength)
        {
            data[newLength .. currentLength] = T.init;
        }

        currentLength = newLength;
    }

    @disable void opBinary(string op : "~")(T){}

    @disable void opBinary(string op : "~")(T[]){}

    /// Append an element
    T opOpAssign(string op : "~")(T itemToAppend)
    {
        assert(currentLength < maxLength, "StackArray at maximum length. cannot append");
        length = currentLength + 1;
        return data[currentLength - 1] = itemToAppend;
    }

    /// Concatenate with a slice
    T[] opOpAssign(string op : "~")(T[] sliceToConcatenate)
    {
        assert(currentLength + sliceToConcatenate.length <= maxLength,
            "slice to concatenate is too long");
        auto start = currentLength;
        currentLength += sliceToConcatenate.length;
        return data[start .. currentLength] = sliceToConcatenate[];
    }

    /// Remove all elements
    void clear()
    {
        length = 0;
    }

    /// Get a slice of the stored data
    inout(T)[] opSlice() inout
    {
        return data[0 .. currentLength];
    }

    /// Access stored element at indx
    ref inout(T) opIndex(size_t indx) inout
    {
        return data[indx];
    }

    alias opSlice this;
}

///
StackArray!(T, maxLen) toStackArray(T, size_t maxLen, R)(R range)
if (isForwardRange!R && (is(ElementType!R : T)))
{
    assert(range.walkLength <= maxLen);
    StackArray!(T, maxLen) sa;
    range.each!(r => sa ~= r);
    return sa;
}

///
@Name("StackArray") pure @safe unittest
{
    int[2] toAppend = [4, 5];

    () pure nothrow @nogc @safe {
        auto sa = StackArray!(int, 5)();
        assert(sa.length == 0);

        sa ~= 3;
        assert(sa.length == 1);
        
        sa ~= toAppend[];
        assert(sa.length == 3);
        assert(sa[0] == 3);
        assert(sa[1] == 4);
        assert(sa[2] == 5);

        sa.length = 5;
        assert(sa.length == 5);
        assert(sa[0] == 3);
        assert(sa[1] == 4);
        assert(sa[2] == 5);
        assert(sa[3] == 0);
        assert(sa[4] == 0);

        sa.clear;
        assert(sa.length == 0);
        assert(sa[].empty);
        assert(sa[].walkLength == 0);
        assert(sa[].length == 0);

        sa.length = 3.iota.length;
        copy(3.iota, sa[]);
        assert(sa[0] == 0);
        assert(sa[1] == 1);
        assert(sa[2] == 2);        
    } ();
}

///
@Name("StackArray (errors)") pure @system unittest
{
    StackArray!(int, 3) sa;
    sa ~= 1;   
    sa ~= [2,3];    // No more space now!   

    (sa ~= [4, 5]).throwsWithMsg("slice to concatenate is too long");
    (sa ~= 4).throwsWithMsg("StackArray at maximum length. cannot append");

    sa.clear;
    sa ~= [1,2];
    (sa ~= [3,4]).throwsWithMsg("slice to concatenate is too long");
}

/******************************************************************************
Returns the factorial of the input. Supports inputs up to 20.
*/
ulong factorial(ulong n) pure nothrow @nogc @safe
{
    assert(n <= 20, "factorial only accepts arguments up to 20");   
    return n == 0 ? 1 : iota(1UL, n + 1).fold!((a, b) => a * b)(1UL);
}
///
@Name("factorial") pure nothrow @nogc @safe unittest
{
    // Might as well test all the input values
    
    assert(factorial(0) == 1UL);    
    assert(factorial(1) == 1UL);
    assert(factorial(2) == 1UL * 2);
    assert(factorial(3) == 1UL * 2 * 3);
    assert(factorial(4) == 1UL * 2 * 3 * 4);
    assert(factorial(5) == 1UL * 2 * 3 * 4 * 5);
    assert(factorial(6) == 1UL * 2 * 3 * 4 * 5 * 6);
    assert(factorial(7) == 1UL * 2 * 3 * 4 * 5 * 6 * 7);
    assert(factorial(8) == 1UL * 2 * 3 * 4 * 5 * 6 * 7 * 8);
    assert(factorial(9) == 1UL * 2 * 3 * 4 * 5 * 6 * 7 * 8 * 9);
    assert(factorial(10) == 1 * 2 * 3 * 4 * 5 * 6 * 7 * 8 * 9 * 10);
    assert(factorial(11) == 1UL * 2 * 3 * 4 * 5 * 6 * 7 * 8 * 9 * 10 * 11);
    assert(factorial(12) == 1UL * 2 * 3 * 4 * 5 * 6 * 7 * 8 * 9 * 10 * 11 
        * 12);
    assert(factorial(13) == 1UL * 2 * 3 * 4 * 5 * 6 * 7 * 8 * 9 * 10 * 11 
        * 12 * 13);
    assert(factorial(14) == 1UL * 2 * 3 * 4 * 5 * 6 * 7 * 8 * 9 * 10 * 11 
        * 12 * 13 * 14);
    assert(factorial(15) == 1UL * 2 * 3 * 4 * 5 * 6 * 7 * 8 * 9 * 10 * 11 
        * 12 * 13 * 14 * 15);
    assert(factorial(16) == 1UL * 2 * 3 * 4 * 5 * 6 * 7 * 8 * 9 * 10 * 11 
        * 12 * 13 * 14 * 15 * 16);
    assert(factorial(17) == 1UL * 2 * 3 * 4 * 5 * 6 * 7 * 8 * 9 * 10 * 11
        * 12 * 13 * 14 * 15 * 16 * 17);
    assert(factorial(18) == 1UL * 2 * 3 * 4 * 5 * 6 * 7 * 8 * 9 * 10 * 11
        * 12 * 13 * 14 * 15 * 16 * 17 * 18);
    assert(factorial(19) == 1UL * 2 * 3 * 4 * 5 * 6 * 7 * 8 * 9 * 10 * 11
        * 12 * 13 * 14 * 15 * 16 * 17 * 18 * 19);
    assert(factorial(20) == 1UL * 2 * 3 * 4 * 5 * 6 * 7 * 8 * 9 * 10 * 11
        * 12 * 13 * 14 * 15 * 16 * 17 * 18 * 19 * 20);
}

///
@Name("factorial (errors)") pure @system unittest
{
    factorial(21).throwsWithMsg("factorial only accepts arguments up to 20");
}

ulong binomial(ulong n, ulong k) pure nothrow @nogc @safe
{
    assert(n <= 20, "binomial only accepts arguments up to 20");
    assert(k <= n, "bad binomial input");
    return (factorial(n) / factorial(k)) / factorial(n - k);
}
///
@Name("binomial") pure nothrow @nogc @safe unittest
{
    // Test exhaustively up to n=5
    assert(binomial(0, 0) == 1);
   
    assert(binomial(1, 0) == 1);
    assert(binomial(1, 1) == 1);

    assert(binomial(2, 0) == 1);
    assert(binomial(2, 1) == 2);
    assert(binomial(2, 2) == 1);

    assert(binomial(3, 0) == 1);
    assert(binomial(3, 1) == 3);
    assert(binomial(3, 2) == 3);
    assert(binomial(3, 3) == 1);

    assert(binomial(4, 0) == 1);
    assert(binomial(4, 1) == 4);
    assert(binomial(4, 2) == 6);
    assert(binomial(4, 3) == 4);    
    assert(binomial(4, 4) == 1);    

    assert(binomial(5, 0) == 1);
    assert(binomial(5, 1) == 5);
    assert(binomial(5, 2) == 10);
    assert(binomial(5, 3) == 10);    
    assert(binomial(5, 4) == 5);
    assert(binomial(5, 5) == 1);

    // A few larger tests, answers verified with Mathematica
    assert(binomial(20, 10) == 184_756);
    assert(binomial(17, 7) == 19_448);
}


/******************************************************************************
template evaluating to true if T is an input range with element type
implicitly convertible to E. TO DO: maybe a more general facility for this?
*/
template isInputRangeOf(T, E)
{
    enum isEmptySlice = is(T == typeof([]));
    enum isInputRangeOf = isEmptySlice || (isInputRange!T && is(ElementType!T : E));
}
///
@Name("isInputRangeOf") pure @safe unittest
{
    alias AI = int[];
    alias ACI = const(int)[];
    static assert(isInputRangeOf!(AI, const(int)));
    static assert(isInputRangeOf!(ACI, const(int)));
    static assert(isInputRangeOf!(AI, int));
    static assert(isInputRangeOf!(ACI, int));

    alias API = int*[];
    alias ACPI = const(int*)[];
    static assert(isInputRangeOf!(API, const(int*)));
    static assert(isInputRangeOf!(ACPI, const(int*)));
    static assert(isInputRangeOf!(API, int*));
    static assert(!isInputRangeOf!(ACPI, int*));

    alias R = typeof(5.iota);
    static assert(isInputRangeOf!(R, int));
    static assert(isInputRangeOf!(R, uint));
    static assert(!isInputRangeOf!(R, string));
}

/******************************************************************************
template evaluating to true if T is an input range of input ranges of types
implicitly convertible to E. TO DO: maybe a more general facility for this?
*/
template isInputRangeOfInputRangeOf(T, E)
{
    enum isInputRangeOfInputRangeOf = isInputRange!T && isInputRange!(ElementType!T)
        && is(ElementType!(ElementType!T) : E);
}
///
@Name("isInputRangeOfInputRangeOf") unittest
{
    alias AAI = int[][];
    alias AACI = const(int)[][];
    static assert(isInputRangeOfInputRangeOf!(AAI, const(int)));
    static assert(isInputRangeOfInputRangeOf!(AACI, const(int)));
    static assert(isInputRangeOfInputRangeOf!(AAI, int));
    static assert(isInputRangeOfInputRangeOf!(AACI, int));

    alias AAPI = int*[][];
    alias AACPI = const(int*)[][];
    static assert(isInputRangeOfInputRangeOf!(AAPI, const(int*)));
    static assert(isInputRangeOfInputRangeOf!(AACPI, const(int*)));
    static assert(isInputRangeOfInputRangeOf!(AAPI, int*));
    static assert(!isInputRangeOfInputRangeOf!(AACPI, int*));

    alias RoR = typeof(5.iota.map!(k => k.iota));
    static assert(isInputRangeOfInputRangeOf!(RoR, int));
    static assert(isInputRangeOfInputRangeOf!(RoR, uint));
    static assert(!isInputRangeOfInputRangeOf!(RoR, string));
}

void swapPop(S, T)(ref S[] unorderedArray, T index)
{
    assert(index < unorderedArray.length);
    unorderedArray[index] = unorderedArray[$-1];
    unorderedArray = unorderedArray[0..$-1];
}

string prettyTime(T)(T timer)
{
    import std.algorithm : findSplit;
    import std.string : replace;
    return timer.peek.to!string.findSplit(",")[0].findSplit("and")[0].replace("hnsecs", "ns");
}

auto replaceEmptyLiteral(R, T)(T input)
{
    static if(is(T == typeof([])))
    {
        // This function is only designed to fix inputs
        // that are the empty literal [] 
        assert(input.empty);
        return R[].init;
    }
    else
    {
        return input;
    }
}
///
@Name("replaceEmptyLiteral") pure nothrow @nogc @safe unittest
{
    static assert(is(typeof(replaceEmptyLiteral!int([])) == int[]));
    int[] x;
    assert(replaceEmptyLiteral!int([]) == x);
}

void dump(alias variable)()
{
  import std.stdio : writefln;
  writefln("%s = %s",
           variable.stringof,
           variable);
}

auto parseParameterFile(string[][] parametersUsed)(string parameterFileName)
{
    auto numberedLines = File(parameterFileName, "r").byLineCopy.array.enumerate(1);
    string[string] valueStrings;

    foreach(lineNum, line; numberedLines)
    {
        if(line.strip.empty || line.startsWith("#")) continue;
        auto lineParts = line.findSplit("=");
        auto paramName = lineParts[0].strip;
        auto separator = lineParts[1];
        auto paramValueString = lineParts[2].strip;
        
        auto location = "line %s of parameter file: %s".format(lineNum, parameterFileName);
        assert(separator == "=", "no '=' found on " ~ location);
        assert(!paramName.empty, "no parameter name found on " ~ location);
        assert(!paramValueString.empty, "no parameter value found on " ~ location);
        assert(paramName !in valueStrings, "repeated parameter found on " ~ location);
        bool isValidParam = parametersUsed.any!(p => paramName == p[1]);
        assert(isValidParam, "unknown parameter found on " ~ location);

        valueStrings[paramName] = paramValueString;
    }

    enum declarations = parametersUsed.map!(typeAndNamePair => "%s    %s;\n"
        .format(typeAndNamePair[0], typeAndNamePair[1]));
    mixin("struct Parameters {\n" ~ declarations.join ~ "}\n");
    Parameters parameters;

    static foreach(typeAndNamePair; parametersUsed)
    {{
        auto paramType = typeAndNamePair[0];
        auto paramName = typeAndNamePair[1];
        assert(paramName in valueStrings,
            "parameter '" ~ paramName ~ "' not found in parameter file: " ~ parameterFileName);

        auto paramValueString = valueStrings[paramName];
        mixin("parameters.%s = paramValueString.to!%s;".format(paramName, paramType));        
    }}

    return parameters;
}

string toPrettyString(T)(T value)
{
    enum fields = __traits(allMembers, T);
    string prettyString;
    static foreach (fieldName; fields)
    {
        prettyString ~= fieldName ~ " = ";
        mixin("prettyString ~= value.%s.to!string;".format(fieldName));
        prettyString ~= "\n";
    }
    return prettyString;
}

// TO DO: Unittests for above function