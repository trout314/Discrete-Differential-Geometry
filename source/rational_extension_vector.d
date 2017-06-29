import factoring : sqrtSquarePart, squareFreePart, squarePart;
import rational : rational, Rational;
import std.algorithm : all, copy, equal, map, sum;
import std.bigint : BigInt;
import std.conv : to;
import std.format : format;
import std.meta : allSatisfy;
import std.range : array, drop, ElementType, enumerate, iota, isForwardRange,
    walkLength, zip;
import std.traits : hasFunctionAttributes, isInstanceOf, ReturnType;

import utility : staticIota, subsetsOfSize;

import std.stdio : writeln;

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
auto simplexPoints(int dim)()
{
    alias Vec = REVector!(simplexRoots(dim).array);
    return iota(0, dim + 1).map!(k => Vec(simplexCoefs(dim)[k]));
}

///
unittest // TO DO: BigInt won't let me put @safe here!
{
    /* We check that the vertices of the 3-simplex are as they should be:
    v0 = (0,         0,       0)
    v1 = (1,         0,       0)
    v2 = (1/2, (1/2)√3,       0)
    v3 = (1/2, (1/6)√3, (1/3)√6)
    */
    assert(simplexPoints!3.map!(pt => pt.toString).equal([
        "(0, 0, 0)",
        "(1, 0, 0)",
        "(1/2, (1/2)√3, 0)",
        "(1/2, (1/6)√3, (1/3)√6)"]));

    alias vec = reVector!(1, 3, 6);
    alias r = (a, b) => rational(BigInt(a), BigInt(b));

    auto pts = [vec(r(0,1), r(0,1), r(0,1)),
                vec(r(1,1), r(0,1), r(0,1)),
                vec(r(1,2), r(1,2), r(0,1)),
                vec(r(1,2), r(1,6), r(1,3))];
    
    assert(simplexPoints!3.equal(pts));
    assert(simplexPoints!3.all!(v => v.roots == [1,3,6]));
    static assert(isForwardRange!(ReturnType!(simplexPoints!3)));

    foreach(pair; simplexPoints!3.subsetsOfSize(2))
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
struct REVector(int[] roots_, RationalType_ = Rational!BigInt) 
if (roots_.all!(r => r>0))
{
    static assert(isInstanceOf!(Rational, RationalType));
    alias RationalType = RationalType_;
    alias IntType = RationalType_.IntType;

    ///
    static immutable roots = roots_;

    ///
    static immutable dimension = roots_.length;

    ///
    const(RationalType)[] coefficients() const
    {
        return coefs[];
    }

    ///
    auto opUnary(string op)() if (op == "-" || op == "+")
    {
        static if(op == "-")
        {
            return typeof(this)(coefs[].map!(c => -c));
        }
        else
        {
            return this;
        }
    }

    ///
    auto opBinary(string op)(REVector!(roots, RationalType) rhs) 
    if (op == "-" || op == "+")
    {
        mixin(q{
            return REVector!(roots, RationalType)(zip(coefs[], rhs.coefs[])
                .map!(pair => pair[0]} ~ op ~ q{pair[1]));
        });
    }

    ///
    auto opBinaryRight(string op)(RationalType scalar) if (op == "*")
    {
        return typeof(this)(coefs[].map!(c => scalar * c));
    }

    ///
    this(R)(R coefficients_) if (isForwardRange!R)
    {
        assert(coefficients_.walkLength == dimension);
        copy(coefficients_, coefs[]);
    }

    ///
    string toString() const
    {
        // square root symbol UTF-8 (hex)
        auto rootSym = x"E2 88 9A";

        string result = "(";
        foreach (indx, coef, root; zip(coefs[], roots[]).enumerate)
        {
            // Take care of initial rational part
            if (coef == RationalType(IntType(0)))
            {
                result ~= "0";
            }
            else if (root == 1 && coef == 1)
            {
                result ~= "1";
            }
            else if (root == 1 && coef != 1)
            {
                result ~= coef.to!string;
            }
            else if (root != 1 && coef == 1)
            {
                result ~= rootSym ~ root.to!string;
            }
            else if (root != 1 && coef != 1)
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
    
    int opCmp()(const typeof(this) rhs) const
    {
        return (coefs < rhs.coefs) ? -1 : (coefs > rhs.coefs);
    }
private:
    RationalType[dimension] coefs;
}

///
unittest // TO DO: BigInt won't allow @safe.
{
    alias vec = reVector!(1, 3, 6);
    alias r = (a, b) => rational(BigInt(a), BigInt(b));

    auto v = simplexPoints!3.array[3];
    auto w = simplexPoints!3.array[2];

    // Negate a vector (and for completeness unary "+")
    assert(-v == vec(-r(1,2), -r(1,6), -r(1,3)));
    assert(+v == v);

    // Addition and subtration of vectors
    assert(v + w == vec(r(1,1),  r(2,3), r(1,3)));
    assert(v - w == vec(r(0,1), -r(1,3), r(1,3)));

    // Scalar multiplication
    assert(r(6,1) * v == vec(r(3,1),  r(1,1), r(2,1)));

    // Check distributivity
    assert(r(2,3) * (v + w) == r(2,3) * v + r(2,3) * w);
    
    // Can be compared using dictionary order on coords
    assert(w > v);
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
        return REVector!(roots, Rational!T)(coefficients[]);
    }
}

///
unittest
{
    alias r = rational;
    auto v = reVector!(3,5,11)(r(1), r(2, 3), r(6));
    static assert(is(typeof(v.coefs[0]) == Rational!int));
    assert(v.toString == "(√3, (2/3)√5, (6/1)√11)");
}

///
pure unittest
{
    alias r = rational;
    auto v = reVector!(1,3,6)(r(1), r(2), r(3));
    assert(dotProduct(v,v) == 67);
}

RationalType dotProduct(int[] roots, RationalType)(
    REVector!(roots, RationalType) v, REVector!(roots, RationalType) w)
{
    static assert(v.roots.equal(w.roots));
    static assert(v.coefs.length == w.coefs.length);
    static assert(v.roots.length == v.coefs.length);
    
    RationalType total;
    foreach(indx; 0 .. roots.length)
    {
        // NOTE: work-around BigInt not having this(int) (WHY?!)
        RationalType r;
        r = roots[indx];
        total += r * v.coefs[indx] * w.coefs[indx];
    }
    return total;
}

pure @safe unittest
{
    alias vec = reVector!(1,3);
    alias r = rational;

    auto v = vec(r(1,3), r(2,5));
    auto w = vec(r(1,2), r(2,3));

    assert(dotProduct(v, w) == r(29,30));
}

///
RationalType distanceSquared(int[] roots, RationalType)(
    REVector!(roots, RationalType) v, REVector!(roots, RationalType) w)
{
    return dotProduct(v - w, v - w);
}

///
pure @safe unittest
{
    alias r = rational;
    alias vec = reVector!(1,2,3);
    assert(distanceSquared(
        vec(r(1), r(2, 5), r(11, 2)),    
        vec(r(0), r(1, 5), r(11, 2))) == r(27,25));
}

// Additional tests. NOTE: BigInt won't let us put @safe here
pure unittest
{
    foreach(dim; staticIota!(1, 8))
    {
        foreach(pair; simplexPoints!dim.subsetsOfSize(2))
        {
            assert(distanceSquared(pair[0], pair[0]) == 0);
            assert(distanceSquared(pair[1], pair[1]) == 0);       
            assert(distanceSquared(pair[0], pair[1]) == 1);
            assert(distanceSquared(pair[1], pair[0]) == 1);       
        }
    }
}

//------------------------------------------------------------------------------
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

pure @safe unittest
{
    alias r = rational;

    assert(reVector!(1,3,6)(r(1, 1), r(3, 7), r(11, 2)).toString ==
        "(1, (3/7)√3, (11/2)√6)");

    assert(reVector!(1,3,6)(r(1, 2), r(0, 1), r(1, 1)).toString ==
        "(1/2, 0, √6)");
}

