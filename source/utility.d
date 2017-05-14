/*******************************************************************************
* Checks if items of type T can be compared with the less-than operation. Note 
* that this means all other comparison operations are also valid, per dlang
* rules about operator overloading. See
* $(LINK https://dlang.org/operatoroverloading.html#eqcmp)
*/ 
enum bool isLessThanComparable(T) = is(typeof({
    bool b = T.init < T.init;   
}));

///
@safe pure nothrow @nogc unittest
{
    struct NotComparable;
    static assert(!isLessThanComparable!NotComparable);
 
    struct Comparable
    {
        /* Note that the "auto ref" template here allows us to have a single
           function that accepts both rvalues and lvalues. See
           https://dlang.org/spec/template.html#auto-ref-parameters. */
        int opCmp()(auto ref const Comparable rhs) const
        {
            return 0;
        }
    }
    static assert(isLessThanComparable!Comparable);
}

/*******************************************************************************
* Checks if items of type T can be compared for equality with ==. See
* $(LINK https:dlang.org/operatoroverloading.html#eqcmp)
*/
enum bool isEqualityComparable(T) = is(typeof(
{
    bool b = T.init == T.init;
}));

///
pure nothrow @nogc @safe unittest
{
    struct NotComparable;
    static assert(!isEqualityComparable!NotComparable);
 
    struct Comparable
    {
        bool opEquals()(auto ref const Comparable rhs) const
        {
            return true;
        }
    }
    static assert(isEqualityComparable!Comparable);
}
