import std.algorithm : all, copy, equal, map, sum;
import std.conv : to;
import std.format : format;
import std.meta : allSatisfy;
import std.range : array, drop, enumerate, iota, isForwardRange, walkLength,
    zip;
import std.traits : hasFunctionAttributes, ReturnType, isInstanceOf;

import utility : staticIota, subsetsOfSize;
import rational : rational, Rational;
import factoring : squarePart, sqrtSquarePart, squareFreePart;

/*******************************************************************************
Returns a forward range which lists the points in the regular unit-edge simplex
of dimension `dim`. The points are returned as REVectors. The first point is
always the origin (0,0,...,0) and from then the non-zero entries appear on the
left. For example, in dimension three the points would be: (0,0,0), (1,0,0),
(1/2,(1/2)√3,0), (1/2,(1/6)√3,(1/3)√6). 

Params:
    dim = dimension of the simplex

Returns:
    Forward range of `REVector`s listing points in the simplex. 
*/
auto simplexPoints(int dim)() pure nothrow @nogc @safe
{
    alias Vec = REVector!(simplexRoots(dim).array);
    return iota(0, dim + 1).map!(k => Vec(simplexCoefs(dim)[k]));
}

///
unittest
{
    /* We check that the vertices of the 3-simplex are as they should be:
    v0 = (0,         0,       0)
    v1 = (1,         0,       0)
    v2 = (1/2, (1/2)√3,       0)
    v3 = (1/2, (1/6)√3, (1/3)√6)
    */

    alias r = rational;
    alias vec = reVector!(1, 3, 6);
    auto pts = [vec(  r(0),   r(0),   r(0)),
                vec(  r(1),   r(0),   r(0)),
                vec(r(1,2), r(1,2),   r(0)),
                vec(r(1,2), r(1,6), r(1,3))];
    
    assert(simplexPoints!3.equal(pts));
    assert(simplexPoints!3.all!(v => v.roots == [1,3,6]));
    static assert(isForwardRange!(ReturnType!(simplexPoints!3)));
 
    foreach(pair; simplexPoints!4.subsetsOfSize(2))
    {
        assert(distanceSquared(pair[0], pair[1]) == 1);
        assert(distanceSquared(pair[1], pair[0]) == 1);       
    }
}


/*******************************************************************************
A "rational-extension vector" that stores vectors of the form:
(r1√d1, r2√d2, ... , rN√dN) where r1, ..., rN are rational numbers and
d1, ... , dN are integers. The rational coefficients may vary at runtime, but 
the integers under the roots are part of the type.
*/
struct REVector(int[] roots_) if (roots_.all!(r => r>0))
{
    ///
    alias roots = roots_;

    ///
    static immutable dimension = roots_.length;

    ///
    typeof(this) opUnary(string op)() if (op == "-" || op == "+")
    {
        static if(op == "-")
        {
            return typeof(this)(this.rationalCoefs[].map!(c => -c));
        }
        else
        {
            return this;
        }
    }

    ///
    typeof(this) opBinary(string s)(typeof(this) rhs) if (s == "-" || s == "+")
    {
        mixin(q{
            return zip(this.rationalCoefs[], rhs.rationalCoefs[])
                .map!(pair => pair[0]} ~ s ~ q{pair[1])
                .to!(typeof(this));
        });
    }

    ///
    typeof(this) opBinaryRight(string op)(Rational!int scalar) if (op == "*")
    {
        return this.rationalCoefs[].map!(c => scalar * c).to!(typeof(this));
    }

    ///
    this(T)(T coefficients) if (isForwardRange!T)
    {
        assert(coefficients.walkLength == dimension);
        copy(coefficients, rationalCoefs[]);
    }

    // square root symbol UTF-8 (hex)
    auto rootSym = x"E2 88 9A";

    ///
    string toString()
    {
        string result = "(";
        foreach (indx, coef, root; zip(rationalCoefs[], roots[]).enumerate)
        {
            // Take care of initial rational part
            if (coef == 0)
            {
                result ~= "0";
            }
            else if ((root == 1) && (coef == 1))
            {
                result ~= "1";
            }
            else if ((root == 1) && (coef != 1))
            {
                result ~= coef.to!string;
            }
            else if ((root != 1) && (coef == 1))
            {
                result ~= rootSym ~ root.to!string;
            }
            else if ((root != 1) && (coef != 1))
            {
                result ~= "(" ~ coef.to!string ~ ")" ~ rootSym ~ root.to!string;
            }

            // Append comma and space if needed
            if (indx + 1 < dimension)
            {
                result ~= ", ";
            }
        }
        result ~= ")";
        return result;
    }
private:
    Rational!int[dimension] rationalCoefs;
}

///
unittest
{
    auto pts = simplexPoints!3.array;
    alias vec = reVector!(1,3,6);
    alias r = rational;
    
    assert(-pts[3] == vec(-r(1,2), -r(1,6), -r(1,3)));
    assert(+pts[3] == pts[3]);

    assert(pts[3] - pts[2] == vec(r(0), -r(1,3), r(1,3)));
    assert(pts[3] + pts[2] == vec(r(1),  r(2,3), r(1,3)));
    assert(  r(6) * pts[3] == vec(r(3),    r(1),   r(2)));

    // Check distributivity
    assert(r(2,3) * (pts[2] + pts[3]) == r(2,3) * pts[2] + r(2,3) * pts[3]);
}

/// Convenience function for creating REVectors
template reVector(R...)
{
    enum int[] roots = [R];
    static assert(roots.length > 1, "must have at least one root");

    auto reVector(T, size_t dim)(Rational!T[dim] coefficients...)
    {
        static assert(coefficients.length == roots.length, "expected "
            ~ roots.length.to!string ~ " coefficients but recieved "
            ~ coefficients.length.to!string);
        return REVector!roots(coefficients[]);
    }
}

///
unittest
{
    alias r = rational;

    auto v = reVector!(1,3,6)(r(1), r(2), r(3));
    assert(dotProduct(v,v) == 67);
}

Rational!int dotProduct(int[] roots)(REVector!roots v, REVector!roots w)
{
    static assert(v.roots.equal(w.roots));
    static assert(v.rationalCoefs.length == w.rationalCoefs.length);
    static assert(v.roots.length == v.rationalCoefs.length);
    
    auto total = rational(0);
    foreach(indx; 0 .. roots.length)
    {
        total += roots[indx] * v.rationalCoefs[indx] * w.rationalCoefs[indx];
    }
    return total;
}

unittest
{
    alias vec = reVector!(1,3);
    alias r = rational;

    auto v = vec(r(1,3), r(2,5));
    auto w = vec(r(1,2), r(2,3));

    assert(dotProduct(v, w) == r(29,30));
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

    assert(reVector!(1,3,6)(r(1, 1), r(3, 7), r(11, 2)).toString ==
        "(1, (3/7)√3, (11/2)√6)");

    assert(reVector!(1,3,6)(r(1, 2), r(0, 1), r(1, 1)).toString ==
        "(1/2, 0, √6)");
}

///
Rational!int distanceSquared(int[] roots)(REVector!roots v, REVector!roots w)
{
    return dotProduct(v - w, v - w);
}

///
unittest
{
    alias r = rational;
    alias vec = reVector!(1,2,3);
    assert(distanceSquared(
        vec(r(1), r(2, 5), r(11, 2)),    
        vec(r(0), r(1, 5), r(11, 2))) == r(27,25));
}

// Additional unittests
unittest
{
    foreach(dim; staticIota!(1, 8))
    {
         foreach(pair; simplexPoints!dim.subsetsOfSize(2))
        {
            assert(distanceSquared(pair[0], pair[1]) == 1);
            assert(distanceSquared(pair[1], pair[0]) == 1);       
        }
    }
}