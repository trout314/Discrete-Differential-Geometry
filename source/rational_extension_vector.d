import std.algorithm : all, copy, equal, map;
import std.conv : to;
import std.format : format;
import std.range : array, enumerate, iota, isForwardRange, walkLength;
import std.traits : hasFunctionAttributes, ReturnType;

import std.rational : rational, Rational;
import factoring : squarePart, sqrtSquarePart, squareFreePart;

/*******************************************************************************
Returns a forward range that ... TO DO

Params:
    dim = 

Returns:
    Forward range of `REVector`s giving the  ... 
*/
auto simplexPoints(int dim)() pure nothrow @nogc @safe
{
    alias Vec = REVector!dim;
    return iota(0, dim + 1).map!(k => Vec(simplexCoefs(dim)[k]));
}

///
pure nothrow @safe unittest
{
    /* We check that the vertices of the 3-simplex are as they should be:
    v0 = (0,         0,       0)
    v1 = (1,         0,       0)
    v2 = (1/2, (1/2)√3,       0)
    v3 = (1/2, (1/6)√3, (1/3)√6)
    */

    alias Vec = reVector;
    alias r = rational;
    auto pts = [Vec(  r(0),   r(0),   r(0)),
                Vec(  r(1),   r(0),   r(0)),
                Vec(r(1,2), r(1,2),   r(0)),
                Vec(r(1,2), r(1,6), r(1,3))];

    assert(simplexPoints!3.equal(pts));
    assert(simplexPoints!3.all!(v => v.roots == [1,3,6]));

    static assert(isForwardRange!(ReturnType!(simplexPoints!3)));
}

struct REVector(int dim)
{
    ///
    string toString()
    {
        string result = "(";
        foreach (i, coef; rationalCoefs)
        {
            // Take care of initial rational part
            if (coef == 0)
            {
                result ~= "0";
            }
            else if ((roots_[i] == 1) && (coef == 1))
            {
                result ~= "1";
            }
            else if (roots_[i] == 1)
            {
                result ~= to!string(coef);
            }
            else
            {
                result ~= "(" ~ to!string(coef) ~ ")";
            }

            // Append root portion if needed
            if ((roots_[i] > 1) && (coef != 1) && (coef != 0))
            {
                result ~= x"E2 88 9A"; // square root symbol UTF-8 (hex)
                result ~= to!string(roots_[i]);
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

    // TO DO: Add some constraints to this template?
    this(T)(T coefs)
    {
        assert(coefs.walkLength == dim);
        copy(coefs, rationalCoefs[]);
    }

    static this()
    {
        foreach (index, root; enumerate(simplexRoots(dim)))
        {
            roots_[index] = root;
        }
    }

    @property auto roots() const
    {
        return roots_;
    }

private:
    static immutable int[dim] roots_;
    Rational!int[dim] rationalCoefs;
}

auto reVector(T, size_t dim)(Rational!T[dim] coefficients...)
{
    alias r = Rational!T;
    return REVector!dim(coefficients[]);
}

auto dotProduct(size_t dim)(REVector!dim v0, REVector!dim v1)
{

} 


private:

auto simplexRoots(int dim)
{
    assert(dim >= 0);
    return iota(1, dim + 1).map!(k => squareFreePart(k * (k + 1) / 2));
}

unittest
{
    assert(simplexRoots(10).array == [1, 3, 6, 10, 15, 21, 7, 1, 5, 55]);
    foreach (radicand; simplexRoots(50))
    {
        assert(radicand > 0);
        assert(squarePart(radicand) == 1);
        assert(squareFreePart(radicand) == radicand);
    }

    static assert(hasFunctionAttributes!(simplexRoots,
        "pure", "nothrow", "@nogc", "@safe"));   
}

auto simplexCoefs(int dim)
{
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
