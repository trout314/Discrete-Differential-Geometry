import core.bitop : popcnt;
import fluent.asserts;
import std.algorithm : all, canFind, copy, each, filter, find, joiner,
    map, sort, sum;
import std.conv : to;
import std.format : format;
import std.meta : AliasSeq, allSatisfy, anySatisfy;
import std.range : array, chain, drop, ElementType, empty, enumerate, front,
    iota, isForwardRange, isInputRange, popFront, repeat, save, take, walkLength;
import std.traits : lvalueOf, rvalueOf;
import unit_threaded : Name;

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
@Name("isLessThanComparable (basic types)")
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
@Name("isEqualityComparable (basic types)")
pure nothrow @nogc @safe unittest
{
    import std.meta : AliasSeq, allSatisfy, anySatisfy;

    alias ComparableTypes = AliasSeq!(bool, byte, ubyte, short, ushort, int, 
        uint, long, ulong, float, double, real, ifloat, idouble, ireal, char,
        wchar, dchar, cfloat, cdouble, creal,  int*, void*);
    static assert(allSatisfy!(isEqualityComparable, ComparableTypes));    
}

/// Struct Tests
@Name("isEqualityComparable (structs)")
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
@Name("isEqualityComparable (classes)")
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
    catch(Throwable exception)
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
    throw new Error("throwsWithMsg failed because expression did not throw",
        file, line);
}

// TO DO: Tests for throwsWithMsg

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
@Name("binarySequences") @safe pure unittest
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
        return this.keys.canFind(key);
    }

    /// Get a lazy range returning the keys in increasing order sorted
    auto keys() const
    {
        return data.map!(r => r.key);
    }

    /// Get a lazy range returning the values
    auto values() const
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
private:
    Record[] data;
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

    sm[2] = "bubba";
    assert(sm[2] == "bubba");


    assert(sm.keys.array == [2, 5]);
    assert(sm.values.array == ["bubba", "hello"]);
}

///
@Name("SmallMap (errors)") pure @system unittest
{
    SmallMap!(int, string) sm;
    sm.insert(5, "hello");
    sm.insert(5, "nope").throwsWithMsg!Error("key already present");
}

///
@Name("SmallMap.keys (pure nothrow @nogc @safe)") pure @safe unittest
{    
    SmallMap!(int, string) sm;
    sm.insert(1, "hello");
    sm.insert(2, "goodbye");

    () pure nothrow @nogc @safe {
        assert(sm.keys.front == 1);
    }();
}

///
@Name("SmallMap.values (pure nothrow @nogc @safe)") pure @safe unittest
{    
    SmallMap!(int, string) sm;
    sm.insert(1, "hello");
    sm.insert(2, "goodbye");

    () pure nothrow @nogc @safe {
        assert(sm.values.front == "hello");
    }();
}

///
@Name("SmallMap.opIndex (pure nothrow @nogc @safe)")  pure @safe unittest
{    
    SmallMap!(int, string) sm;
    sm.insert(1, "hello");
    sm.insert(2, "goodbye");

    () pure nothrow @nogc @safe {
        assert(sm[1] == "hello");
    }();
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
    assert([1,3].isSubsetOf([1,3,4]));
    assert(![2,3].isSubsetOf([1,3,4]));

    int[] empty;
    assert(empty.isSubsetOf([6,8]));
    assert(empty.isSubsetOf(empty));

    /* Note that using an empty string literal in the above example will not
    work since [] has type void[] */
    static assert(!__traits(compiles, [].isSubsetOf([6, 8])));
}
///
@Name("isSubsetOf (pure nothrow @nogc @safe)") pure nothrow @safe unittest
{
    assert([1,2,3].isSubsetOf([1,2,3,4]));
    assert(2.iota.isSubsetOf(8.iota));

    () pure nothrow @nogc @safe {
        assert(iota(3,7).isSubsetOf(iota(10)));
    }();
}


/*******************************************************************************
Returns true if the uint `x` has a 1 at position `pos` and 0 otherwise (where
position 0 is the lesat significant bit and position 31 the highest.)
*/
bool hasOneAtBit(uint x, size_t pos) pure nothrow @nogc @safe
{
    assert(pos < uint.sizeof * 8, "bad bit position");
    return  ((1 << pos) & x) > 0;
}
///
@Name("hasOneAtBit") pure nothrow @nogc @safe unittest
{
    assert(1.hasOneAtBit(0));
    assert(iota(1,32).all!(pos => !1.hasOneAtBit(pos)));

    assert(iota(32).all!(pos => uint.max.hasOneAtBit(pos)));

    assert(!12.hasOneAtBit(0));
    assert(!12.hasOneAtBit(1));
    assert(12.hasOneAtBit(2));
    assert(12.hasOneAtBit(3));
    assert(iota(4,32).all!(pos => !12.hasOneAtBit(pos)));
}

///
@Name("hasOneAtBit (errors)") pure @system unittest
{
    // An error is thrown if the bit position is outside the uint
    7.hasOneAtBit(32).throwsWithMsg!Error("bad bit position");
}

auto subsetFromUint(R)(R set, uint whichToKeep) if (isInputRange!R)
{
    assert(whichToKeep < (1 << set.walkLength),
        "1 bits found in positions not corresponsing to elements in set");

    static struct SubsetFromUintRange
    {
        uint whichToKeep_;
        R set_;

        auto front()
        {
            assert(!set_.empty);
            while(!whichToKeep_.hasOneAtBit(0))
            {
                set_.popFront;
                whichToKeep_ >>= 1;
            }

            assert(!set_.empty);
            return set_.front;
        }

        auto popFront()
        {
            assert(!set_.empty, 
                "tried to popFront an empty SubsetFromUintRange");

            set_.popFront;
            whichToKeep_ >>= 1;       

            while(!whichToKeep_.hasOneAtBit(0) && !empty)
            {
                set_.popFront;
                whichToKeep_ >>= 1;
            }
            assert(whichToKeep_.hasOneAtBit(0) || empty);
        }

        auto empty()
        {
            return whichToKeep_ == 0;
        }

        static if(isForwardRange!R)
        {
            auto save()
            {
                return SubsetFromUintRange(whichToKeep_, set_);
            }
        }
    }

    return SubsetFromUintRange(whichToKeep, set);
}
///
@Name("subsetsFromUint") pure nothrow @safe unittest
{
    assert([3,4,5].subsetFromUint(0).empty);
    assert([3,4,5].subsetFromUint(1).array == [3]);
    assert([3,4,5].subsetFromUint(2).array == [4]);
    assert([3,4,5].subsetFromUint(3).array == [3,4]);
    assert([3,4,5].subsetFromUint(4).array == [5]);
    assert([3,4,5].subsetFromUint(5).array == [3,5]);
    assert([3,4,5].subsetFromUint(6).array == [4,5]);
    assert([3,4,5].subsetFromUint(7).array == [3,4,5]);
}

///
@Name("subsetsFromUint (errors)") pure nothrow @system unittest
{
    [3,4,5].subsetFromUint(1 << 3).throwsWithMsg!Error(
        "1 bits found in positions not corresponsing to elements in set");
}

///
@Name("subsetsFromUint (pure nothrow @nogc @safe)") pure nothrow @nogc @safe
unittest
{
    auto r = iota(1,5).subsetFromUint(7);
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

    assert(subsetSize > 0, "subset size must be positive");
    assert(subsetSize <= setSize, "subset size must be at most the size of the set");
    assert(setSize <= 31, "subset size must be at most 31");

    static struct SubsetsOfSizeRange
    {
    private:
        /* We will think of the bits in whichToKeep as indexed starting with
        zero for the least significant bit. Highest bit is empty flag. */
        uint whichToKeep;
        R set_;            // Underlying set from which to draw elements
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
            if(whichToKeep == lastToKeep)
            {
                whichToKeep = 1 << 31;
                return;
            }

            auto currentPos = len - 1;
            if(whichToKeep.hasOneAtBit(currentPos))
            {
                // Skip over additional ones
                while(whichToKeep.hasOneAtBit(currentPos))
                {
                    --currentPos;
                }
                assert(currentPos > 0 && currentPos < len);
                assert(!whichToKeep.hasOneAtBit(currentPos));

                immutable numOnesSeen = len - currentPos - 1;

                // Find another one to move
                while(!whichToKeep.hasOneAtBit(currentPos))
                {
                    --currentPos;
                }
                assert(currentPos >= 0 && currentPos < len);
                assert(whichToKeep.hasOneAtBit(currentPos));

                whichToKeep &= ~(1 << currentPos);

                iota(currentPos + 1, currentPos + 2 + numOnesSeen).
                    each!(pos => whichToKeep |= (1 << pos));
                
                iota(currentPos + 2 + numOnesSeen, len).
                    each!(pos => whichToKeep &= ~(1 << pos));  
            }
            else // Zero at current position
            {
                assert(!whichToKeep.hasOneAtBit(currentPos));

                // Find a one to move
                while(!whichToKeep.hasOneAtBit(currentPos))
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

        static if(isForwardRange!R)
        {
            auto save()
            {
                return SubsetsOfSizeRange(whichToKeep, set_.save);
            }
        }
    }

    return SubsetsOfSizeRange((1 << subsetSize) - 1, set);
}
///
@Name("subsetsOfSize") @system unittest 
{
    [1,2,3,4].subsetsOfSize(1).map!array.should.containOnly(
        [[1], [2], [3], [4]]);

    [1,2,3,4].subsetsOfSize(2).map!array.should.containOnly(
        [[1,2], [1,3], [1,4], [2,3], [2,4], [3,4]]);

    [1,2,3,4].subsetsOfSize(3).map!array.should.containOnly(
        [[1, 2, 3], [1, 2, 4], [1, 3, 4], [2, 3, 4]]);
        
    [1,2,3,4].subsetsOfSize(4).map!array.should.containOnly(
        [[1, 2, 3, 4]]);

    [1].subsetsOfSize(1).map!array.should.containOnly([[1]]);

    3.iota.subsetsOfSize(1).map!array.should.containOnly([[0], [1], [2]]);
    3.iota.subsetsOfSize(2).map!array.should.containOnly([[0,1], [0,2], [1,2]]);
    3.iota.subsetsOfSize(3).map!array.should.containOnly([[0,1,2]]);

    iota(31).subsetsOfSize(1).walkLength.should.equal(31);
    iota(10).subsetsOfSize(2).walkLength.should.equal(45);
}

///
@Name("subsetsOfSize (pure nothrow @nogc @safe)") pure nothrow @safe unittest
{
    int[4] set = [1,2,3,4];
    int[2][6] expectedSubsets = [[1,2], [1,3], [1,4], [2,3], [2,4], [3,4]];
    int[2][6] actualSubsets;

    () pure nothrow @nogc @safe {
        auto subsetsRange = set[].subsetsOfSize(2);
        auto saved = subsetsRange.save;

        foreach(indx; 0 .. 6)
        {
            copy(subsetsRange.front, actualSubsets[indx][]);
            subsetsRange.popFront;
        }
        assert(subsetsRange.empty);
        assert(!saved.empty);

        // Don't want to test the specific order, so sort
        actualSubsets[].sort();

        foreach(indx; 0 .. 6)
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

auto subsets(R)(R set) if (isInputRange!R)
{
    assert(set.walkLength <= 31);

    static struct SubsetsRange
    {
    private:
        R set_;
        uint whichToKeep;
        bool empty_;
    public:
        auto front()
        {
            assert(!empty);
            return subsetFromUint(set_, whichToKeep);
        }

        auto popFront()
        {
            assert(!empty);
            ++whichToKeep;
        }

        auto empty()
        {
            return whichToKeep == 1 << set_.walkLength;
        }

        static if(isForwardRange!R)
        {
            auto save()
            {
                return SubsetsRange(set_, whichToKeep);
            }
        }
    }

    return SubsetsRange(set, 1);
}
///
@Name("subsets") @system unittest
{
    [1,2,3].subsets.map!array.should.containOnly([
        [1], [2], [3],
        [1, 2], [1, 3], [2, 3],
        [1, 2, 3]
    ]);

    iota(1,5).subsets.map!array.should.containOnly([
        [1], [2], [3], [4],
        [1,2], [1,3], [1,4], [2,3], [2,4], [3,4],
        [1,2,3], [1,2,4], [1,3,4], [2,3,4],
        [1,2,3,4]
    ]);

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
        alias staticIota = AliasSeq!(
            staticIota!(begin, middle), staticIota!(middle, end));
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
    alias func = (a,b,c,d) => a+b+c+d;
    static assert(func(seq2) == 10);
}

/*******************************************************************************
*/
auto capture(Range, Data...)(Range range, Data data) if (isInputRange!Range)
{
    alias dataDeclarations = (len) => len.iota.map!(
        indx => "Data[%s] d%s;".format(indx, indx)).joiner("\n").array;

    alias dataArgs = (len) => len.iota.map!(
        indx => "d%s".format(indx)).joiner(",").array;

    static struct CaptureRange
    {
        Range range_;
        mixin(dataDeclarations(data.length));

        auto front()
        {
            static struct CaptureFront
            {
                typeof(range_.front) result;
                alias result this;
                mixin(dataDeclarations(data.length));
            }

            mixin("return CaptureFront(range_.front, " 
                ~ dataArgs(data.length) ~ ");");
        }

        auto popFront()
        {
            range_.popFront;
        }

        auto empty()
        {
            return range_.empty;
        }

        static if(isForwardRange!Range)
        {
            auto save()
            {
                mixin("return CaptureRange(range_, " 
                    ~ dataArgs(data.length) ~ ");");
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

/*******************************************************************************
*/
struct StackArray(T, size_t maxLength)
{
private:
    T[maxLength] data;
    size_t currentLength;
public:
    auto length() const
    {
        return currentLength;
    }

    auto length(size_t newLength)
    {
        assert(newLength <= maxLength, "new length must be at most max length");

        if(newLength < currentLength)
        {
            data[newLength .. currentLength] = T.init;
        }

        currentLength = newLength;      
    }

    @disable void opBinary(string op : "~")(T){}
    @disable void opBinary(string op : "~")(T[]){}

    T opOpAssign(string op : "~")(T itemToConcatenate)
    {
        ++currentLength;
        length = currentLength;
        return data[currentLength - 1] = itemToConcatenate;
    }

    T[] opOpAssign(string op : "~")(T[] sliceToConcatenate)
    {
        auto start = currentLength;
        currentLength += sliceToConcatenate.length;
        return data[start .. currentLength] = sliceToConcatenate[];
    }

    void clear()
    {
        length = 0;
    }

    inout(T)[] opSlice() inout
    {
        return data[0 .. currentLength];
    }

    ref inout(T) opIndex(size_t indx) inout
    {
        return data[indx];
    }

    alias opSlice this;
}
///

@Name("StackArray") @nogc unittest
{
    int[2] toAppend = [4,5];

    auto a = StackArray!(int, 5)();
    a ~= 3;
    assert(a.length == 1);
    a ~= toAppend[];
    assert(a.length == 3);

    assert(a[0] == 3);
    assert(a[1] == 4);
    assert(a[2] == 5);
    a.length = 5;
}