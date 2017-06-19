import std.algorithm : equal, filter, group, map, reduce, setDifference, sort;
import std.range : array, chain, isForwardRange, stride;
import std.traits : ReturnType;
import std.typecons : staticIota, tuple;

/*******************************************************************************
Returns a forward range that computes the prime factors of `num`` in increasing 
order.

Params:
    num = a positive integer

Returns:
    Forward range of `int`s giving the prime factors of `num` in increasing
    order.
*/
PrimeFactorsRange primeFactors(int num) pure nothrow @nogc @safe
{
    assert(num > 0);
    return PrimeFactorsRange(num);
}

///
pure nothrow @safe unittest
{
    assert(primeFactors(1).array == []);
    assert(primeFactors(2 * 3 * 5).array == [2, 3, 5]);
    assert(primeFactors(2 * 2 * 2 * 3 * 3 * 5).array == [2, 2, 2, 3, 3, 5]);
    assert(primeFactors(1223).array == [1223]);

    // Everything works at compile-time too
    static assert(primeFactors(24).array == [2, 2, 2, 3]);
    static assert(isForwardRange!(ReturnType!primeFactors));
}

/*******************************************************************************
Returns a forward range containing (in increasing order) those primes which
appear to an odd power in the prime-factorization of `num`. The product of these
primes is `squareFreePart(num)` and is the largest divisor of `num` not 
divisible by any square. This is also what's left under the the root after
"simplifying" the square root of `num` by moving any squares outside the root.

Params:
    num = a positive integer

Returns:
    Forward range of `int`s giving (in increasing order) those prime factors of
    `num` which appear to an odd power

See_Also:
    `squareFreePart`
*/
auto squareFreePrimeFactors(int num) pure nothrow @nogc @safe
{
    assert(num > 0);
    return num
        .primeFactors
        .group  // group into tuples of form (prime, power)
        .filter!(r => r[1] % 2 == 1)    // only odd powers
        .map!(r => r[0]);               // want the primes
}

///
pure nothrow @safe unittest
{
    assert(squareFreePrimeFactors(1).array == []);
    assert(squareFreePrimeFactors(11 * 11).array == []);
    assert(squareFreePrimeFactors(2 * 2 * 5 * 5 * 5 * 7 * 11).array == [5, 7, 11]);
    static assert(squareFreePrimeFactors(24).array == [2, 3]);

    static assert(isForwardRange!(ReturnType!squareFreePrimeFactors));
}

/*******************************************************************************
Returns a forward range containing (in increasing order) the prime factors of
the largest square divisor of `num`. The product of these primes is
`squarePart(num)`.

Params:
    num = a positive integer

Returns:
    Forward range of `int`s containing (in increasing order) the prime factors 
    of the largest square divisor of `num`

See_Also:
    `squarePart`
*/
auto squarePrimeFactors(int num) pure nothrow @nogc @safe
{
    assert(num > 0);  
    return num.primeFactors.setDifference(num.squareFreePrimeFactors);
}

///
pure nothrow @safe unittest
{
    assert(squarePrimeFactors(1).array == []);
    assert(squarePrimeFactors(11 * 11).array == [11, 11]);
    assert(squarePrimeFactors(2 * 2 * 5 * 5 * 5 * 7).array == [2, 2, 5, 5]);
    static assert(squarePrimeFactors(24).array == [2, 2]);

    static assert(isForwardRange!(ReturnType!squarePrimeFactors));
}

/*******************************************************************************
Returns a forward range containing (in increasing order) the prime factors of 
the square root of the largest square divisor of `num`. The product of these 
primes is `sqrtSquarePart(num)`.

Params:
    num = a positive integer

Returns:
    Forward range of `int` s containing (in increasing order) the prime factors 
    of the square root of the largest square divisor of `num`

See_Also:
    `sqrtSquarePart`
*/
auto sqrtSquarePrimeFactors(int num) pure nothrow @nogc @safe
{
    assert(num > 0);
    return num.squarePrimeFactors.stride(2);
}

///
pure nothrow @safe unittest
{
    assert(sqrtSquarePrimeFactors(1).array == []);
    assert(sqrtSquarePrimeFactors(11 * 11).array == [11]);
    assert(sqrtSquarePrimeFactors(2 * 2 * 5 * 5 * 5 * 7 * 11).array == [2, 5]);

    static assert(sqrtSquarePrimeFactors(24).array == [2]);

    static assert(isForwardRange!(ReturnType!sqrtSquarePrimeFactors));
}

/*******************************************************************************
Returns the largest square-free divisor of `num`.

Params:
    num = a positive integer

Returns:
    Largest square-free divisor of `num`

See_Also:
    `squareFreePrimeFactors`
*/
int squareFreePart(int num) pure nothrow @nogc @safe
{
assert(num > 0);
    return 1.reduce!((a, b) => a * b)(num.squareFreePrimeFactors);
}

///
pure nothrow @nogc @safe unittest
{
    assert(squareFreePart(1) == 1);
    assert(squareFreePart(11 * 11) == 1);
    assert(squareFreePart(2 * 2 * 5 * 5 * 5 * 7 * 11) == 5 * 7 * 11);
    static assert(squareFreePart(24) == 2 * 3);
}

/*******************************************************************************
Returns the largest square divisor of `num`.

Params:
    num = a positive integer

Returns:
    Largest square divisor of `num`

See_Also:
    `squarePrimeFactors`
*/
int squarePart(int num) pure nothrow @nogc @safe
{
    assert(num > 0);
    return 1.reduce!((a, b) => a * b)(num.squarePrimeFactors);
}

///
pure nothrow @nogc @safe unittest
{
    assert(squarePart(1) == 1);
    assert(squarePart(11 * 11) == 11 * 11);
    assert(squarePart(2 * 2 * 5 * 5 * 5 * 7 * 11) == 2 * 2 * 5 * 5);
    static assert(squarePart(24) == 2 * 2);
}

/*******************************************************************************
Returns the square root of the largest square divisor of `num`.

Params:
    num = a positive integer

Returns:
    Square root of the largest square divisor of `num`

See_Also:
    `sqrtSquarePrimeFactors`
*/
int sqrtSquarePart(int num) pure nothrow @nogc @safe
{
    assert(num > 0);
    return 1.reduce!((a, b) => a * b)(num.sqrtSquarePrimeFactors);
}

///
pure nothrow @nogc @safe unittest
{
    assert(sqrtSquarePart(1) == 1);
    assert(sqrtSquarePart(11 * 11) == 11);
    assert(sqrtSquarePart(2 * 2 * 5 * 5 * 5 * 7 * 11) == 2 * 5);
    static assert(sqrtSquarePart(24) == 2);
}

// Some additional tests
pure nothrow @safe unittest
{
    static assert(primeFactors(1).array == []);
    static assert(squareFreePrimeFactors(1).array == []);
    static assert(squarePrimeFactors(1).array == []);
    static assert(sqrtSquarePrimeFactors(1).array == []);

    enum testFactorLists = tuple([2], [3], [7919], [2, 7], [3, 11], [5,
            7529], [2, 2], [13, 13], [7529, 7529], [2, 3, 5, 7, 11, 13, 17],
            [41, 41, 53, 67], [2, 3, 5, 11, 101, 101], [2, 2, 2, 11, 17], [2,
            2, 2, 11, 11, 17], [2, 11, 4259, 4259], [2, 2, 2, 2, 2, 2, 2, 2],
            [5, 5, 5, 5, 5, 5, 5], [3, 3, 3, 3, 7, 7, 7]);

    foreach (factors; testFactorLists)
    {
        // get the number n whose prime factorization is the given factor list 
        immutable n = 1.reduce!((a, b) => a * b)(factors);

        assert(n.primeFactors.equal(factors));
        static assert(n.primeFactors.equal(factors));

        auto allPrimes = chain(n.squareFreePrimeFactors, n.squarePrimeFactors);
        assert(sort(allPrimes.array).equal(factors));
        assert(n == n.squareFreePart * n.squarePart);
        assert(n == n.squareFreePart * n.sqrtSquarePart^^2);

        enum allPrimes_ = chain(n.squareFreePrimeFactors, n.squarePrimeFactors);
        static assert(sort(allPrimes_.array).equal(factors));
        static assert(n == n.squareFreePart * n.squarePart);
        static assert(n == n.squareFreePart * n.sqrtSquarePart^^2);
    }

    immutable maxRuntimeTest = 10_000;
    foreach (immutable n; 1 .. maxRuntimeTest + 1)
    {
        assert(1.reduce!((a,b) => a*b)(n.primeFactors) == n);
        auto allPrimes = chain(n.squareFreePrimeFactors, n.squarePrimeFactors);
        assert(sort(allPrimes.array).equal(n.primeFactors));
        assert(n == n.squareFreePart * n.squarePart);
        assert(n == n.squareFreePart * n.sqrtSquarePart^^2);
    }

    enum maxCompiletimeTest = 100;
    foreach (immutable n; staticIota!(1, maxCompiletimeTest))
    {
        static assert(1.reduce!((a,b) => a*b)(n.primeFactors) == n);
        enum allPrimes = chain(n.squareFreePrimeFactors, n.squarePrimeFactors);        
        static assert(sort(allPrimes.array).equal(n.primeFactors));
        static assert(n == n.squareFreePart * n.squarePart);
        static assert(n == n.squareFreePart * n.sqrtSquarePart^^2);
    }
}

/*******************************************************************************
Returns the smallest (non-trivial if possible) factor of `num`. If `num == 1` 
then the function returns `1` otherwise it returns the smallest factor of `num` 
greater than one. Note that if `num > 1` this factor is necessarily prime.

Params:
    num = a positive integer

Returns:
    `int` which is the smallest prime factor of `num` if such exists, otherwise 1.
*/
int lowestFactor(int num) pure nothrow @nogc @safe
{
    assert(num > 0);

    if (num <= 3)
    {
        return num;
    }

    int possibleFactor = 2;
    while (num % possibleFactor != 0)
    {
        if (possibleFactor == 2)
        {
            possibleFactor = 3;
        }
        else
        {
            possibleFactor += 2;
        }
    }

    assert(num % possibleFactor == 0);   
    return possibleFactor;
}

///
pure nothrow @nogc @safe unittest
{
    assert(lowestFactor(1) == 1);
    assert(lowestFactor(2) == 2);
    assert(lowestFactor(3) == 3);
    assert(lowestFactor(4) == 2);
    assert(lowestFactor(11) == 11);
    assert(lowestFactor(2 * 2 * 5 * 11) == 2);
    assert(lowestFactor(13 * 13) == 13);

    static assert(lowestFactor(2*3) == 2);
}

private:

struct PrimeFactorsRange
{
    int number;
    int factor;

    this(int num) pure nothrow @nogc @safe
    {
        assert(num > 0);

        number = num;
        factor = 1;
    }

    @property int front() pure nothrow @nogc @safe
    {
        assert(!empty); // Must have prime factor remaining
        if (factor == 1)
        {
            factor = number.lowestFactor;
        }
        return factor;
    }

    @property bool empty() const pure nothrow @nogc @safe
    {
        return number == 1;
    }

    void popFront() pure nothrow @nogc @safe
    {
        assert(!empty); // Must have prime factor remaining
        if (factor == 1)
        {
            factor = number.lowestFactor;
        }
        number /= factor;
        factor = number.lowestFactor;
    }

    @property PrimeFactorsRange save() pure nothrow @nogc @safe const
    {
        return this;
    }
}

pure nothrow @safe unittest
{
    auto pf = primeFactors(6);
    auto pfCopy = pf.save;

    pf.popFront;
    assert(pfCopy.array == [2, 3]);
    assert(pf.array == [3]);
}