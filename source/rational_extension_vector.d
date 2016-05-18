import std.rational; // Note: not actually a standard package yet
import factoring;

int[] visibilityRoots(int dim)
in
{
    assert(dim > 0);
}
out (result)
{
    assert(result.length == dim);
    foreach (d; result)
    {
        assert(d > 0);
        assert(squarePart(d) == 1);
        assert(squareFreePart(d) == d);
    }
}
body
{
    int[] result;
    foreach (k; 1 .. dim + 1)
        result ~= squareFreePart((k * (k + 1)) / 2);
    return result;
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
            if (i < rationalCoefs.length - 1)
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

    static immutable int[dim] roots = visibilityRoots(dim);
private:
    Rational!int[dim] rationalCoefs;
}
