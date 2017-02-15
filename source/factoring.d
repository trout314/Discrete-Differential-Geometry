version (unittest)
{
    import std.algorithm : equal, reduce, sort;
    import std.range : array, isForwardRange, chain;
    import std.traits : ReturnType, hasFunctionAttributes;
    import std.typecons : staticIota, tuple;
}

/++
Returns a forward range that enumerates the prime factors of num in increasing order.

Params:
    num =   a positive integer

Returns:
    a forward range of the prime factors of num in increasing order.
+/
auto primeFactors(int num)
{
    assert(num > 0);
    return PrimeFactorsRange(num);
}

///
unittest
{
    static assert(isForwardRange!(ReturnType!primeFactors));
    static assert(hasFunctionAttributes!(primeFactors, "pure", "nothrow", "@nogc", "@safe"));

    assert(primeFactors(1).array == []);
    assert(primeFactors(2 * 3 * 5).array == [2, 3, 5]);
    assert(primeFactors(2 * 2 * 2 * 3 * 3 * 5).array == [2, 2, 2, 3, 3, 5]);
    assert(primeFactors(1223).array == [1223]);

    // Everything works at compile-time too
    static assert(primeFactors(24).array == [2, 2, 2, 3]);
}

/++
Returns a range containing (in increasing order) those primes which appear to an odd power in the prime-factorization of `num`. The product of these primes is `squareFreePart(num)` and is the largest divisor of `num` not divisible by any square. This is also what's left under the the root after "simplifying" the square root of `num` by moving any squares outside the root.

Params:
    num = a positive integer

Returns:
    a range of `int`s containing (in increasing order) those prime factors of `num` which appear to an odd power.

See_Also:
    `squareFreePart`
+/
auto squareFreePrimeFactors(int num)
{
    assert(num > 0);

    import std.algorithm : filter, group, map;

    return num
        .primeFactors
        .group  // group into tuples of form (prime, power)
        .filter!(r => r[1] % 2 == 1)    // only odd powers
        .map!(r => r[0]);               // want the primes
}

///
unittest
{
    assert(squareFreePrimeFactors(1).array == []);
    assert(squareFreePrimeFactors(11 * 11).array == []);
    assert(squareFreePrimeFactors(2 * 2 * 5 * 5 * 5 * 7 * 11).array == [5, 7, 11]);
    static assert(squareFreePrimeFactors(24).array == [2, 3]);

    static assert(hasFunctionAttributes!(squareFreePrimeFactors, "pure", "nothrow", "@nogc", "@safe"));
}

/++
Returns a range containing (in increasing order) the prime factors of the largest square divisor of `num`. The product of these primes is `squarePart(num)`.

Params:
    num = a positive integer

Returns:
    a range of `int`s containing (in increasing order) the prime factors of the largest square divisor of `num`

See_Also:
    `squarePart`
+/
auto squarePrimeFactors(int num)
{
    assert(num > 0);

    import std.algorithm : setDifference;
    
    return num.primeFactors.setDifference(num.squareFreePrimeFactors);
}

///
unittest
{
    assert(squarePrimeFactors(1).array == []);
    assert(squarePrimeFactors(11 * 11).array == [11, 11]);
    assert(squarePrimeFactors(2 * 2 * 5 * 5 * 5 * 7 * 11).array == [2, 2, 5, 5]);
    static assert(squarePrimeFactors(24).array == [2, 2]);

    static assert(hasFunctionAttributes!(squarePrimeFactors, "pure", "nothrow", "@nogc", "@safe"));
}

/++
Returns a range containing (in increasing order) the prime factors of the square root of the largest square divisor of `num`. The product of these primes is `sqrtSquarePart(num)`.

Params:
    num =   a positive integer

Returns:
    a range of `int`s containing (in increasing order) the prime factors of the square root of the largest square divisor of `num`

See_Also:
    `sqrtSquarePart`
+/
auto sqrtSquarePrimeFactors(int num)
{
    assert(num > 0);
    import std.range : stride;

    return num.squarePrimeFactors.stride(2);
}

///
unittest
{
    assert(sqrtSquarePrimeFactors(1).array == []);
    assert(sqrtSquarePrimeFactors(11 * 11).array == [11]);
    assert(sqrtSquarePrimeFactors(2 * 2 * 5 * 5 * 5 * 7 * 11).array == [2, 5]);

    static assert(sqrtSquarePrimeFactors(24).array == [2]);

    static assert(hasFunctionAttributes!(sqrtSquarePrimeFactors, "pure", "nothrow", "@nogc", "@safe"));
}

/++
Returns the largest square-free divisor of `num`.

Params:
    num =   a positive integer

Returns:
    largest square-free divisor of `num`

See_Also:
    `squareFreePrimeFactors`
+/
auto squareFreePart(int num)
{
    assert(num > 0);

    import std.algorithm : reduce;

    return 1.reduce!((a, b) => a * b)(num.squareFreePrimeFactors);
}

///
unittest
{
    assert(squareFreePart(1) == 1);
    assert(squareFreePart(11 * 11) == 1);
    assert(squareFreePart(2 * 2 * 5 * 5 * 5 * 7 * 11) == 5 * 7 * 11);
    static assert(squareFreePart(24) == 2 * 3);

    static assert(hasFunctionAttributes!(squareFreePart, "pure", "nothrow", "@nogc", "@safe"));    
}

/++
Returns the largest square divisor of `num`.

Params:
    num =   a positive integer

Returns:
    largest square divisor of `num`

See_Also:
    `squarePrimeFactors`
+/
auto squarePart(int num)
{
    assert(num > 0);

    import std.algorithm : reduce;

    return 1.reduce!((a, b) => a * b)(num.squarePrimeFactors);
}

///
unittest
{
    assert(squarePart(1) == 1);
    assert(squarePart(11 * 11) == 11 * 11);
    assert(squarePart(2 * 2 * 5 * 5 * 5 * 7 * 11) == 2 * 2 * 5 * 5);
    static assert(squarePart(24) == 2 * 2);

    static assert(hasFunctionAttributes!(squarePart, "pure", "nothrow", "@nogc", "@safe"));    
}

/++
Returns the square root of the largest square divisor of `num`.

Params:
    num =   a positive integer

Returns:
    square root of the largest square divisor of `num`

See_Also:
    `sqrtSquarePrimeFactors`
+/
auto sqrtSquarePart(int num)
{
    assert(num > 0);

    import std.algorithm : reduce;

    return 1.reduce!((a, b) => a * b)(num.sqrtSquarePrimeFactors);
}

///
unittest
{
    assert(sqrtSquarePart(1) == 1);
    assert(sqrtSquarePart(11 * 11) == 11);
    assert(sqrtSquarePart(2 * 2 * 5 * 5 * 5 * 7 * 11) == 2 * 5);
    static assert(sqrtSquarePart(24) == 2);

    static assert(hasFunctionAttributes!(sqrtSquarePart, "pure", "nothrow", "@nogc", "@safe"));    
}

// Some additional tests
unittest
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

        auto allPrimeFactors = sort(chain(n.squareFreePrimeFactors, n.squarePrimeFactors).array);
        assert(allPrimeFactors.equal(factors));
        assert(n == n.squarePart * n.squareFreePart);
        assert(n == n.sqrtSquarePart * n.sqrtSquarePart * n.squareFreePart);

        enum allPrimeFactors_ = sort(chain(n.squareFreePrimeFactors, n.squarePrimeFactors).array);
        static assert(allPrimeFactors_.equal(factors));
        static assert(n == n.squarePart * n.squareFreePart);
        static assert(n == n.sqrtSquarePart * n.sqrtSquarePart * n.squareFreePart);
    }

    immutable maxRuntimeTest = 10_000;
    foreach (immutable n; 1 .. maxRuntimeTest + 1)
    {
        assert(1.reduce!((a,b) => a*b)(n.primeFactors) == n);
        auto allPrimeFactors = sort(chain(n.squareFreePrimeFactors, n.squarePrimeFactors).array);
        assert(allPrimeFactors.equal(n.primeFactors));
        assert(n == n.squarePart * n.squareFreePart);
        assert(n == n.sqrtSquarePart * n.sqrtSquarePart * n.squareFreePart);
    }

    enum maxCompiletimeTest = 100;
    foreach (immutable n; staticIota!(1, maxCompiletimeTest))
    {
        static assert(1.reduce!((a,b) => a*b)(n.primeFactors) == n);
        enum allPrimeFactors = sort(chain(n.squareFreePrimeFactors, n.squarePrimeFactors).array);        
        static assert(allPrimeFactors.equal(n.primeFactors));
        static assert(n == n.squarePart * n.squareFreePart);
        static assert(n == n.sqrtSquarePart * n.sqrtSquarePart * n.squareFreePart);
    }
}

/++
Returns the lowest (non-trivial if possible) factor of `num`. If `num == 1` then the function returns `1` otherwise it returns the smallest factor of `num` greater than one.

Params:
    num =   a positive integer

Returns:
    `int` which is lowest non-trivial factor of `num` if such exists, otherwise 1 if `num` is 1.
+/
auto lowestFactor(int num)
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
unittest
{
    assert(lowestFactor(1) == 1);
    assert(lowestFactor(2) == 2);
    assert(lowestFactor(3) == 3);
    assert(lowestFactor(4) == 2);
    assert(lowestFactor(11) == 11);
    assert(lowestFactor(2 * 2 * 5 * 11) == 2);
    assert(lowestFactor(13 * 13) == 13);

    static assert(lowestFactor(2*3) == 2);

    static assert(hasFunctionAttributes!(lowestFactor, "pure", "nothrow", "@nogc", "@safe"));    
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