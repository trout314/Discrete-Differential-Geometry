import std.rational; // Note: not actually a standard package yet

auto simplexRoots(int dimension)
{
    assert(dimension >= 0, "simplex dimension must be non-negative");

    import std.algorithm : map;
    import std.range : iota;
    import factoring : squareFreePart;

    return iota(1, dimension + 1).map!(k => squareFreePart(k * (k + 1) / 2));
}

unittest
{
    foreach (radicand; simplexRoots(50))
    {
        assert(radicand > 0);
        assert(squarePart(radicand) == 1);
        assert(squareFreePart(radicand) == radicand);
    }

    import std.range : array;

    assert(simplexRoots(4).array == [1, 3, 6, 10]);
}

auto visibilityCoefsRange(int dim, int simplexPointIndx)
{

    import std.range : chain, iota, only;
    import std.algorithm : map;
    import factoring : sqrtSquarePart;

}

Rational!int[] visibilityCoefs(int dim, int simplexPointIndx)
in
{
    assert(dim > 0);
    assert(simplexPointIndx >= 0);
    assert(simplexPointIndx <= dim);
}
out (result)
{
    // TO DO: Put in some checks...
}
body
{
    import factoring : sqrtSquarePart;

    Rational!int[] result;
    foreach (basisIndx; 1 .. dim + 1)
    {
        auto n = (basisIndx * (basisIndx + 1)) / 2;
        if (basisIndx < simplexPointIndx)
            result ~= sqrtSquarePart(n) * Rational!int(1, 2 * n);
        else if (basisIndx == simplexPointIndx)
            result ~= sqrtSquarePart(n) * Rational!int(1, basisIndx);
        else
            result ~= Rational!int(0, 1);
    }
    return result;
}

RationalExtensionVector!dim[] simplexVecs(int dim)()
{
    RationalExtensionVector!dim[] result;
    foreach (j; 0 .. dim + 1)
    {
        result ~= RationalExtensionVector!dim(visibilityCoefs(dim, j));
    }
    return result;
}

struct RationalExtensionVector(int dim)
{
    string toString()
    {
        import std.conv : to;
        import std.format : format;

        string result = "[";
        foreach (i, coef; rationalCoefs)
        {

            // Take care of initial rational part
            if (coef == 0)
            {
                result ~= "0";
            }
            else if ((roots[i] == 1) && (coef == 1))
            {
                result ~= "1";
            }
            else if (roots[i] == 1)
            {
                result ~= to!string(coef);
            }
            else
            {
                result ~= "(" ~ to!string(coef) ~ ")";
            }

            // Append root portion if needed
            if ((roots[i] > 1) && (coef != 1) && (coef != 0))
            {
                result ~= x"E2 88 9A"; // square root symbol UTF-8 (hex)
                result ~= to!string(roots[i]);
            }

            // Append comma and space if needed
            if (i + 1 < rationalCoefs.length)
                result ~= ", ";
        }
        result ~= "]";
        return result;
    }

    this(Rational!int[] coefs)
    in
    {
        assert(coefs.length == dim);
    }
    body
    {
        rationalCoefs = coefs.dup;
    }

    import std.range : array;

    static immutable int[dim] roots = simplexRoots(dim).array;
private:
    Rational!int[dim] rationalCoefs;
}
