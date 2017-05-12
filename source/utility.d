///
enum bool isLessThanComparable(T) = is(typeof({
    bool b = T.init < T.init;
}));

///
unittest
{
    struct NotComparable;
    static assert(!isLessThanComparable!NotComparable);
 
    struct Comparable
    {
        // opCmp is used to define "<", "<=", etc. For reference, see
        // https://dlang.org/operatoroverloading.html#eqcmp
        // Note that the "auto ref" template here allows us to have a single
        // function that accepts both rvalues and lvalues. See
        // https://dlang.org/spec/template.html#auto-ref-parameters
        int opCmp()(auto ref const Comparable rhs) const
        {
            return 0;
        }
    }
    static assert(isLessThanComparable!Comparable);
}

///
enum bool isEqualityComparable(T) = is(typeof(
{
    bool b = T.init == T.init;
}));

///
unittest
{

}
