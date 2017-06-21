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
    alias Vec = REVector!(
        iota(1, dim + 1).map!(k => squareFreePart(k * (k + 1) / 2)).array);
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

    alias r = rational;
    enum roots = [1,3,6];
    alias vec = reVector!roots;
    auto pts = [vec(  r(0),   r(0),   r(0)),
                vec(  r(1),   r(0),   r(0)),
                vec(r(1,2), r(1,2),   r(0)),
                vec(r(1,2), r(1,6), r(1,3))];
    
    assert(simplexPoints!3.equal(pts));
    assert(simplexPoints!3.all!(v => v.roots == [1,3,6]));
    static assert(isForwardRange!(ReturnType!(simplexPoints!3)));
}

struct REVector(int[] theRoots)
{
    

    // TO DO: Add some constraints to this template?
    this(T)(T coefs)
    {
        assert(coefs.walkLength == dim);
        copy(coefs, rationalCoefs[]);
    }

    static this()
    {
        roots_[] = theRoots[];
    }

    auto roots() const
    {
        return roots_[];
    }
    
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
            else if ((roots_[i] == 1) && (coef != 1))
            {
                result ~= to!string(coef);
            }
            else if ((roots_[i] != 1) && (coef == 1))
            {
                result ~= x"E2 88 9A"; // square root symbol UTF-8 (hex)
                result ~= to!string(roots_[i]);
            }
            else if ((roots_[i] != 1) && (coef != 1))
            {
                result ~= "(" ~ to!string(coef) ~ ")";            
                result ~= x"E2 88 9A";
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
private:
    static immutable dim = theRoots.walkLength;
    static immutable int[dim] roots_;
    Rational!int[dim] rationalCoefs;
}

template reVector(int[] roots)
{
    auto reVector(T, size_t dim)(Rational!T[dim] coefficients...)
    {
        static assert(coefficients.length == roots.walkLength);
        return REVector!roots(coefficients[]);
    }
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

unittest
{
    alias r = rational;

    assert(REVector!([1,3,6])([r(1, 1), r(3, 7), r(11, 2)]).toString ==
        "(1, (3/7)√3, (11/2)√6)");

    assert(REVector!([1,3,6])([r(1, 2), r(0, 1), r(1, 1)]).toString ==
        "(1/2, 0, √6)");

    auto v = REVector!([1,3,6])([r(1,2), r(1,3), r(0,1)]);


    import std.stdio : writeln;
    v.writeln;
    v.roots.writeln;

}