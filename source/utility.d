import std.algorithm : all, canFind, copy, each, equal, filter, find, joiner,
    map, sort, sum;
import std.conv : to;
import std.exception : enforce;
import std.meta : AliasSeq, allSatisfy, anySatisfy;
import std.range : array, chain, drop, ElementType, empty, enumerate, front,
    iota, save, isForwardRange, only, popFront, repeat, take, walkLength;
import std.traits : lvalueOf, rvalueOf;
import core.bitop : popcnt;
import unit_threaded;
import fluent.asserts;

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
    throw new Error("throwsWithMsg failed because expression did not throw",
        file, line);
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
@Name("binarySequences")
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
    void insert(const KeyType key, /* const */ ValueType value)
    {
        enforce(key !in this, "key already present");
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
        enforce(found.length > 0, "SmallMap access error");
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

    sm.insert(5, "nope").throwsWithMsg("key already present");

    assert(sm.keys.array == [2, 5]);
    assert(sm.values.array == ["bubba", "hello"]);
}

/*******************************************************************************
Returns true if set A is contained in set B
*/
auto isSubsetOf(A, B)(A setA, B setB)
if (isForwardRange!A && isForwardRange!B && is(ElementType!A : ElementType!B))
{
    return setA.all!(element => setB.canFind(element));        
}
///
@Name("isSubsetOf")
unittest
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
/*******************************************************************************
Returns true if the uint `x` has a 1 at position `pos` and 0 otherwise (where
position 0 is the lesat significant bit and position 31 the highest.)
*/
auto hasOneAtBit(uint x, size_t pos)
{
    return  (1 << pos) & x;
}

auto subsetFromUint(R)(R set, uint whichToKeep) if (isForwardRange!R)
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

        auto save()
        {
            return SubsetFromUintRange(whichToKeep_, set_);
        }
    }

    return SubsetFromUintRange(whichToKeep, set);
}
///
@Name("subsetsFromUint") unittest
{
    assert([3,4,5].subsetFromUint(0).empty);
    assert([3,4,5].subsetFromUint(1).array == [3]);
    assert([3,4,5].subsetFromUint(2).array == [4]);
    assert([3,4,5].subsetFromUint(3).array == [3,4]);
    assert([3,4,5].subsetFromUint(4).array == [5]);
    assert([3,4,5].subsetFromUint(5).array == [3,5]);
    assert([3,4,5].subsetFromUint(6).array == [4,5]);
    assert([3,4,5].subsetFromUint(7).array == [3,4,5]);

    // This will throw errors
    // assert([3,4,5].subsetFromUint(1 << 31).array == [3,4,5]);

    // TO DO: More tests
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
auto subsetsOfSize(R)(R set, int subsetSize) if (isForwardRange!R)
{
    auto setSize = set.walkLength;

    assert(subsetSize > 0);
    assert(subsetSize <= setSize);
    assert(setSize <= 32);

    static struct SubsetsOfSizeRange
    {
        /* We will think of the bits in whichToKeep as indexed starting with
        zero for the least significant bit. Highest bit is empty flag. */
        uint whichToKeep;
        R set_;            // Underlying set from which to draw elements
        bool empty_;

        auto front()
        {
            assert(!this.empty);
            return set_.subsetFromUint(whichToKeep);
        }

        auto popFront() @nogc
        {
            assert(!this.empty);

            immutable len = set_.walkLength;
            immutable subLen = whichToKeep.popcnt;
            immutable uint lastToKeep = ((1 << subLen) - 1) << (len - subLen);

            // Check for impending emptiness
            if(whichToKeep == lastToKeep)
            {
                empty_ = true;
                return;
            }

            auto currentPos = len - 1;
            if(whichToKeep.hasOneAtBit(currentPos))
            {
                // Skip over additional ones
                while(whichToKeep.hasOneAtBit(currentPos) && (currentPos >= 0))
                {
                    --currentPos;
                }
                assert(currentPos > 0);
                assert(currentPos < len);
                assert(!whichToKeep.hasOneAtBit(currentPos));

                immutable numOnesSeen = len - currentPos - 1;

                // Find another one to move
                while(!whichToKeep.hasOneAtBit(currentPos) && (currentPos >= 0))
                {
                    --currentPos;
                }
                assert(currentPos >= 0);
                assert(currentPos < len);
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
                assert(currentPos >= 0);
                assert(currentPos < len);
                assert(whichToKeep.hasOneAtBit(currentPos));

                whichToKeep &= ~(1 << currentPos);
                whichToKeep |= (1 << (currentPos + 1));
            }
        }

        auto empty()
        {
            return empty_;
        }

        auto save()
        {
            return SubsetsOfSizeRange(whichToKeep, set_, empty_);
        }
    }

    uint initToKeep = (1 << subsetSize) - 1;
    return SubsetsOfSizeRange(initToKeep, set);
}
///
@Name("subsetsOfSize") unittest 
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
}

///
@Name("subsetsOfSize @nogc") @nogc unittest
{
    // TO DO: Improve this test! Ugly...

    int[4] s = [1,2,3,4];

    auto r = s[].subsetsOfSize(2);
    assert(r.front.front == 1);
    assert(!r.empty);
    r.popFront;
    r.popFront;
    r.popFront;
    assert(r.front.front == 2);
    auto r1 = r.save;
    r.popFront;
    r.popFront;
    assert(r.front.front == 3);
    assert(r1.front.front == 2);
}

auto subsets(R)(R set) if (isForwardRange!R)
{
    static struct SubsetsRange
    {
        R set_;
        uint whichToKeep;
        bool empty_;

        auto front()
        {
            assert(!empty);
            return subsetFromUint(set_, whichToKeep);
        }

        auto popFront()
        {
            assert(!empty);
            if(whichToKeep == (1 << set_.walkLength) - 1)
            {
                empty_ = true;
            }
            ++whichToKeep;
        }

        auto empty()
        {
            return empty_;
        }

        auto save()
        {
            return SubsetsRange(set_, whichToKeep, empty_);
        }
    }

    if(set.empty)
    {
        return SubsetsRange(set, 0, true);
    }
    else
    {
        return SubsetsRange(set, 1, false);
    }    
}
///
@Name("subsets") unittest
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
}
///
@Name("subsets pure nothrow @nogc @safe") pure nothrow @nogc @safe unittest
{
    auto r = 3.iota.subsets;
    auto s = r.save;
    assert(r.front.front == 0);    
    r.popFront;
    assert(r.front.front == 1);    
    assert(s.front.front == 0);
    assert(!r.empty);
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
@Name("staticIota")
unittest
{
    alias seq = staticIota!(3, 6);
    static assert(seq[0] == 3);
    static assert(seq[1] == 4);
    static assert(seq[2] == 5);

    alias seq2 = staticIota!(1, 5);
    alias func = (a,b,c,d) => a+b+c+d;
    static assert(func(seq2) == 10);
}