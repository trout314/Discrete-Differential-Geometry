/++
Returns a newly allocated dynamic array containing the prime factors of `num`. This array has the same content as `primeFactorsRange(num).array`.

Params:
    num =   a positive integer

Returns:
    a dynamic array of ints containing the prime factors of `num`.

See_Also:
    `primeFactorsRange`
+/
int[] primeFactors(int num) pure nothrow @safe
{
    assert(num > 0);

    int n = num;
    int[] answer = [];

    if (n < 2)
    {
        return answer;
    }

    int factor = 2;
    while (factor < n)
    {
        if (n % factor == 0)
        {
            answer ~= factor;
            n /= factor;
        }
        factor = n.nextFactorFrom(factor);
    }
    answer ~= factor;
    return answer;
}

///
unittest
{
    assert(primeFactors(1) == []);
    assert(primeFactors(30) == [2, 3, 5]);
    assert(primeFactors(360) == [2, 2, 2, 3, 3, 5]);
    assert(primeFactors(1223) == [1223]);

    // This works, but it is slow. TO DO: Figure out why.
    // assert(primeFactorsRange(30).array == [2, 3, 5]);
}

/++
Returns a newly allocated dynamic array containing (in increasing order) those primes which appear to an odd power in the prime-factorization of `num`. The product of these primes is `squareFreePart(num)` and is the largest divisor of `num` not divisible by any square. This is also what's left under the the root after "simplifying" the square root of `num` by moving any squares outside the root.

Params:
    num =   a positive integer

Returns:
    a dynamic array of ints containing (in increasing order) those prime factors of `num` which appear to an odd power.

See_Also:
    `squareFreePart`
+/
int[] squareFreePrimeFactors(int num) pure nothrow @safe
{
    assert(num > 0);

    import std.algorithm : uniq, count, filter;
    import std.range : array;

    auto factors = num.primeFactors;
    return factors.uniq.filter!(f => factors.count(f) % 2 == 1).array;
}

///
unittest
{
    assert(squareFreePrimeFactors(1) == []);
    assert(squareFreePrimeFactors(11 * 11) == []);
    assert(squareFreePrimeFactors(2 * 2 * 5 * 5 * 5 * 7 * 11) == [5, 7, 11]);
}

/++
Returns a newly allocated dynamic array containing (in increasing order) the prime factors of the largest square divisor of `num`. The product of the primes in this array will be `squarePart(num)`.

Params:
    num =   a positive integer

Returns:
    a dynamic array of ints containing (in increasing order) the prime factors of the largest square divisor of `num`

See_Also:
    `squarePart`
+/
int[] squarePrimeFactors(int num) pure nothrow @safe
{
    assert(num > 0);

    import std.algorithm : setDifference;
    import std.range : array;

    return num.primeFactors.setDifference(num.squareFreePrimeFactors).array;
}

///
unittest
{
    assert(squarePrimeFactors(1) == []);
    assert(squarePrimeFactors(11 * 11) == [11, 11]);
    assert(squarePrimeFactors(2 * 2 * 5 * 5 * 5 * 7 * 11) == [2, 2, 5, 5]);
}

/++
Returns a newly allocated dynamic array containing (in increasing order) the prime factors of the square root of the largest square divisor of `num`. The product of the primes in this array will be `sqrtSquarePart(num)`.

Params:
    num =   a positive integer

Returns:
    a dynamic array of ints containing (in increasing order) the prime factors of the square root of the largest square divisor of `num`

See_Also:
    `sqrtSquarePart`
+/
int[] sqrtSquarePrimeFactors(int num) pure nothrow @safe
{
    assert(num > 0);
    import std.range : array, stride;

    return num.squarePrimeFactors.stride(2).array;
}

///
unittest
{
    assert(sqrtSquarePrimeFactors(1) == []);
    assert(sqrtSquarePrimeFactors(11 * 11) == [11]);
    assert(sqrtSquarePrimeFactors(2 * 2 * 5 * 5 * 5 * 7 * 11) == [2, 5]);
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
int squareFreePart(int num) pure nothrow @safe
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
int squarePart(int num) pure nothrow @safe
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
int sqrtSquarePart(int num) pure nothrow @safe
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
}

/++
Returns an input range that enumerates the prime factors of num. This range has the same content as the array returned by `primeFactors(num)`.

Params:
    num =   a positive integer

Returns:
    an input range containing the prime factors of num.

See_Also:
    `primeFactors`
+/
auto primeFactorsRange(int num) pure nothrow @nogc @safe
{
    assert(num > 0);

    return PrimeFactorsRange(num);
}

///
unittest
{
    static assert(isForwardRange!(ReturnType!primeFactorsRange));

    // assert(primeFactorsRange(20).array == [2, 2, 5]);
    assert(primeFactors(20) == [2, 2, 5]);

    // assert(primeFactorsRange(1).empty);

    // This should work, but eats all memory then crashes. 
    // static assert(primeFactorsRange(6).array == [2,3]);
}

// Some additional tests
unittest
{
    static assert(primeFactors(1) == []);
    static assert(squareFreePrimeFactors(1) == []);
    static assert(squarePrimeFactors(1) == []);
    static assert(sqrtSquarePrimeFactors(1) == []);

    enum testFactorLists = tuple([2], [3], [7919], [2, 7], [3, 11], [5,
            7529], [2, 2], [13, 13], [7529, 7529], [2, 3, 5, 7, 11, 13, 17],
            [41, 41, 53, 67], [2, 3, 5, 11, 101, 101], [2, 2, 2, 11, 17], [2,
            2, 2, 11, 11, 17], [2, 11, 4259, 4259], [2, 2, 2, 2, 2, 2, 2, 2],
            [5, 5, 5, 5, 5, 5, 5], [3, 3, 3, 3, 7, 7, 7]);

    foreach (factors; testFactorLists)
    {
        immutable number = 1.reduce!((a, b) => a * b)(factors);
        assert(number.primeFactors.equal(factors));
        static assert(number.primeFactors.equal(factors));

        // TO DO: Why is this so slow?    
        // assert(number.primeFactorsRange.equal(factors));

        // TO DO: Why is this slow and use so much memory?        
        // static assert(number.primeFactorsRange.equal(factors));
    }

    immutable maxRuntimeTest = 10_000;
    foreach (n; 1 .. maxRuntimeTest + 1)
    {
        assert(merge(n.squareFreePrimeFactors, n.squarePrimeFactors).equal(n.primeFactors));
        assert(n == n.squarePart * n.squareFreePart);
        assert(n == n.sqrtSquarePart * n.sqrtSquarePart * n.squareFreePart);
    }

    enum maxCompiletimeTest = 100;
    foreach (n; staticIota!(1, maxCompiletimeTest))
    {
        static assert(merge(n.squareFreePrimeFactors, n.squarePrimeFactors).equal(n.primeFactors));
        static assert(n == n.squarePart * n.squareFreePart);
        static assert(n == n.sqrtSquarePart * n.sqrtSquarePart * n.squareFreePart);
    }
}

private:

/++
Returns an input range that enumerates the prime factors of num. This range has the same content as the array returned by `primeFactors(num)`.

Params:
    num =   a positive integer

Returns:
    an input range containing the prime factors of num.

See_Also:
    `primeFactors`
+/
int nextFactorFrom(int number, int factor) pure nothrow @nogc @safe
{
    assert(number > 0);
    assert(factor > 0);

    int nextFactor = factor;

    if (factor == 1)
    {
        nextFactor = 2;
    }

    while (number % nextFactor != 0)
    {
        if (nextFactor == 2)
        {
            nextFactor = 3;
        }
        else
        {
            nextFactor += 2;
        }
    }

    assert(number % nextFactor == 0);
    return nextFactor;
}

struct PrimeFactorsRange
{
    int remainingPart;
    int currentPart;

    this(int num) pure nothrow @nogc @safe
    {
        assert(num > 0);

        remainingPart = num;
        currentPart = 1;
        popFront(); // Set currentPart to first prime factor
    }

    @property int front() const pure nothrow @nogc @safe
    {
        assert(!empty); // Must have prime factor remaining
        return currentPart;
    }

    @property bool empty() const pure nothrow @nogc @safe
    {
        return remainingPart == 1;
    }

    void popFront() pure nothrow @nogc @safe
    {
        remainingPart /= currentPart;
        currentPart = remainingPart.nextFactorFrom(currentPart);
    }

    @property PrimeFactorsRange save() const {
        return this;
    }
}

version (unittest)
{
    import std.algorithm : equal, reduce, sort, merge;
    import std.range : array, isForwardRange;
    import std.traits : ReturnType;
    import std.typecons : staticIota, tuple;
}
