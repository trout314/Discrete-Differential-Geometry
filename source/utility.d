import std.conv : to;
import std.meta : AliasSeq, allSatisfy, anySatisfy;
import std.traits : lvalueOf, rvalueOf;

version(unittest)
{
    import std.stdio : writeln;
}

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
unittest
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
import std.conv : to;
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
void throwsWithMsg(ThrownType : Throwable = Error, E)
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
}