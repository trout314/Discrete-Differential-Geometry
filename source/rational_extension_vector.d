auto simplexPoints(int dim)()
{
    import std.range : iota, array;
    import std.algorithm : map;

    return iota(0, dim + 1).map!(
        k => RationalExtensionVector!dim(simplexCoefs(dim)[k]));
}

unittest
{
    // Here we check that the 

    import std.rational : rational;
    alias Point = RationalExtensionVector!3;

    auto x0 = Point([rational(0,1), rational(0,1), rational(0,1)]);
    auto x1 = Point([rational(1,1), rational(0,1), rational(0,1)]);
    auto x2 = Point([rational(1,2), rational(1,2), rational(0,1)]);
    auto x3 = Point([rational(1,2), rational(1,6), rational(1,3)]);

    import std.algorithm : equal;
    assert(simplexPoints!3.equal([x0,x1,x2,x3]));

    import std.traits : hasFunctionAttributes;
    static assert(hasFunctionAttributes!(simplexPoints!3,
        "pure", "nothrow", "@nogc", "@safe"));
}

struct RationalExtensionVector(int dim)
{
    string toString()
    {
        import std.conv : to;
        import std.format : format;

        string result = "(";
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
        result ~= ")";
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

    static this()
    {
        import std.range : enumerate;

        foreach (index, root; enumerate(simplexRoots(dim)))
        {
            roots[index] = root;
        }
    }

private:
    static immutable int[dim] roots;

    import std.rational : Rational;

    Rational!int[dim] rationalCoefs;
}

private:

auto simplexRoots(int dim)
{
    assert(dim >= 0);

    import std.algorithm : map;
    import std.range : iota;
    import factoring : squareFreePart;

    return iota(1, dim + 1).map!(k => squareFreePart(k * (k + 1) / 2));
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
    assert(simplexRoots(10).array == [1, 3, 6, 10, 15, 21, 7, 1, 5, 55]);
}

auto simplexCoefs(int dim)
{
    import std.range : iota;
    import std.algorithm : map;
    import std.rational : rational;
    import factoring : sqrtSquarePart;

    auto computeCoef(int pointIndex, int basisIndex)
    {
        auto n = basisIndex * (basisIndex + 1) / 2;
        auto m = sqrtSquarePart(n);
        if (basisIndex < pointIndex)
        {
            return m * rational(1, 2 * n);
        }
        else if (basisIndex == pointIndex)
        {
            return m * rational(1, basisIndex);
        }
        else
        {
            return rational(0, 1);
        }
    }

    return iota(0, dim + 1).map!(pointIndex => iota(1, dim + 1)
            .map!(basisIndex => computeCoef(pointIndex, basisIndex)));

}

unittest
{
    // TO DO...
}