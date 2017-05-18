version(unittest)
{
    import std.stdio : writeln;
}


/*******************************************************************************
* Checks if items of type T can be compared with the less-than operation. Note 
* that this means all other comparison operations are also valid, per dlang
* rules about operator overloading. See
* $(LINK https://dlang.org/operatoroverloading.html#eqcmp)
*/
template isLessThanComparable(T)
{
    static if (is(T == class))
    {
        enum isLessThanComparable = __traits(isOverrideFunction, T.opCmp);
    }
    else
    {
        enum isLessThanComparable = is(typeof({
            bool b = T.init < T.init;   
        }));
    }
}

/// Basic type tests
@safe pure nothrow @nogc unittest
{
    import std.meta : AliasSeq, allSatisfy, anySatisfy;

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

    struct B
    {
        /* Note that the "auto ref" template here allows us to have a single
           function that accepts both rvalues and lvalues. See
           https://dlang.org/spec/template.html#auto-ref-parameters. */
        int opCmp()(auto ref const B rhs) const
        {
            return 0;
        }
    }
    static assert(isLessThanComparable!B);
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
* Checks if items of type T can be compared for equality with ==. See
* $(LINK https:dlang.org/operatoroverloading.html#eqcmp)
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
}

/// Class tests
pure nothrow @nogc @safe unittest
{
    class A {}
    // static assert(isEqualityComparable!A);
 
    // class B
    // {
    //     bool opEquals()(auto ref const B rhs) const
    //     {
    //         return true;
    //     }
    // }
    // static assert(isEqualityComparable!B);
}


/*******************************************************************************
* Checks if an instance of type S can be constructed from an instance of type T
* TO DO: Didn't see any phobos trait for this, ask about it on forums.
*/
enum bool isConstructible(From, To) = is(typeof({
    // Test to see if we can construct a To from a From. Note that we can't
    // use To(From.init) since VRP (value range propagation) will allow us to
    // construct a byte from an int.
    import std.traits : lvalueOf;
    To t = To(lvalueOf!From);
    // TO DO: Why do we need the special case of string?
})) || is(From == To);

pure nothrow @nogc @safe unittest
{
    static assert(!isConstructible!(int, ubyte));
}

/*******************************************************************************
Asserts that the given expression throws the given type of $(D Throwable).
The $(D Throwable) is caught and does not escape assertThrown. However,
any other $(D Throwable)s $(I will) escape, and if no $(D Throwable)
of the given type is thrown, then an $(D AssertError) is thrown. Also, throws

Params:
    ThrownType = The Throwable type to test for. Default is $(D AssertError)
    expression = The expression to test.
    msg        = Optional message to output on test failure.
    file       = The file where the error occurred.
                 Defaults to $(D __FILE__).
    line       = The line where the error occurred.
                 Defaults to $(D __LINE__).

    Throws:
        $(D AssertError) if the given $(D Throwable) with message $(D msg) is
        not thrown.
*/  
void throwsWithMsg(ThrownType : Throwable = Error, E)
                  (lazy E expression,
                  string msg = null,
                  string file = __FILE__,
                  size_t line = __LINE__)
{
    import std.conv : to;
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
                ~ "type.\n  Actual  : " ~ typeid(exception).to!string ~ "\n  "
                ~ "Expected: " ~ ThrownType.stringof, file, line);
        }
        else if (thrown.msg != msg)
        {
            throw new Error("throwsWithMsg failed with wrong message.\n"
                ~ "  Actual  : " ~ thrown.msg ~ "\n  Expected: " ~ msg, file,
                line);
        }
        else
        {
            return;
        }
    }

}