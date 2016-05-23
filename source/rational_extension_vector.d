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
    import factoring : squarePart, squareFreePart;

    foreach (radicand; simplexRoots(50))
    {
        assert(radicand > 0);
        assert(squarePart(radicand) == 1);
        assert(squareFreePart(radicand) == radicand);
    }

    import std.range : array;

    assert(simplexRoots(4).array == [1, 3, 6, 10]);
}

auto simplexCoefs(int dim, int simplexPointIndex)
{

    import std.range : iota;
    import std.algorithm : map;
    import factoring : sqrtSquarePart;
    import std.rational : rational;

    return iota(1, dim + 1).map!((basisIndex) {
        auto n = basisIndex * (basisIndex + 1) / 2;
        auto m = sqrtSquarePart(n);
        if (basisIndex < simplexPointIndex)
        {
            return m * rational(1, 2 * n);
        }
        else if (basisIndex == simplexPointIndex)
        {
            return m * rational(1, basisIndex);
        }
        else
        {
            return rational(0, 1);
        }
    });
}

unittest
{
    // TO DO...
}

auto simplexVecs(int dim)()
{
    import std.range : iota, array;
    import std.algorithm : map;

    return iota(0, dim + 1).map!(k => RationalExtensionVector!dim(simplexCoefs(dim, k)));
}

unittest
{
    // TO DO...
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
            {
                result ~= ", ";
            }
        }
        result ~= "]";
        return result;
    }

    // TO DO: Add some constraints to this template...
    this(T)(T coefs)
    {
        import std.range : walkLength;
        import std.algorithm : copy;

        assert(coefs.walkLength == dim);
        copy(coefs, rationalCoefs[]);
    }

    import std.rational : Rational;
    import std.range : array;

    static immutable int[dim] roots = simplexRoots(dim).array;
private:
    Rational!int[dim] rationalCoefs;
}
