version(unittest)
{
    import std.range.primitives : isInputRange;
    import std.range : array;
    import std.traits : ReturnType;
    import std.typecons : tuple;
    import std.algorithm : reduce, equal;
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
auto primeFactorsRange(int num) pure nothrow @nogc
{
    assert(num > 0);

    return PrimeFactorsRange(num);
}

///
unittest
{
    static assert(isInputRange!(ReturnType!primeFactorsRange));

    assert(primeFactorsRange(20).array == [2, 2, 5]);
    assert(primeFactors(20) == [2, 2, 5]);

    assert(primeFactorsRange(1).empty);

    // This should work, but eats all memory then crashes. 
    // static assert(primeFactorsRange(1).empty);
}

/++
Returns a newly allocated dynamic array containing the prime factors of `num`. This array had the same content as `primeFactorsRange(num).array`.

Params:
    num =   a positive integer

Returns:
    a dynamic array of ints containing the prime factors of `num`.

See_Also:
    `primeFactorsRange`
+/
public int[] primeFactors(int num) pure nothrow
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

unittest
{
    static assert(primeFactors(1) == []);

    enum testFactorLists = tuple([2], [3], [7919], [2, 7], [3, 11],[5, 7529],
            [2, 2], [13, 13], [7529, 7529], [2, 3, 5, 7, 11, 13, 17],
            [41, 41, 53, 67], [2, 3, 5, 11, 101, 101], [2, 2, 2, 11, 17], [2,
            2, 2, 11, 11, 17], [2, 11, 4259, 4259], [2, 2, 2, 2, 2, 2, 2, 2],
            [5, 5, 5, 5, 5, 5, 5], [3, 3, 3, 3, 7, 7, 7]);

    foreach (factors; testFactorLists)
    {
        immutable number = 1.reduce!product(factors);
        assert(number.primeFactors.equal(factors));
        static assert(number.primeFactors.equal(factors));

        // TO DO: Why is this so slow?    
        // assert(number.primeFactorsRange.equal(factors));

        // TO DO: Why is this slow and use so much memory?        
        // static assert(number.primeFactorsRange.equal(factors));
    }
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
int[] squareFreePrimeFactors(int num)
{
    assert(num > 0);

    import std.algorithm : uniq, count, filter;
    import std.range : array;

    auto factors = num.primeFactors;
    alias oddPower = f => factors.count(f) % 2 == 1;
    return factors.uniq.filter!oddPower.array;
}

///
unittest
{
    assert(squareFreePrimeFactors(1) == []);
    static assert(squareFreePrimeFactors(1) == []);
}

int[] squarePrimeFactors(int num)
{
    assert(num > 0);

    import std.algorithm : setDifference;
    import std.range : array;

    return num.primeFactors.setDifference(num.squareFreePrimeFactors).array;
}

unittest
{
    assert(squarePrimeFactors(1) == []);
    static assert(squarePrimeFactors(1) == []);

}

int[] sqrtSquarePrimeFactors(int num)
{
    assert(num > 0);
    import std.range : array, stride;

    return num.squarePrimeFactors.stride(2).array;
}

unittest
{
    assert(sqrtSquarePrimeFactors(1) == []);
    static assert(sqrtSquarePrimeFactors(1) == []);
}

int squareFreePart(int num)
{
    assert(num > 0);

    import std.algorithm : reduce;

    return 1.reduce!product(num.squareFreePrimeFactors);
}

int squarePart(int num)
{
    assert(num > 0);

    import std.algorithm : reduce;

    return 1.reduce!product(num.squarePrimeFactors);
}

int sqrtSquarePart(int num)
{
    assert(num > 0);

    import std.algorithm : reduce;

    return 1.reduce!product(num.sqrtSquarePrimeFactors);
}

// Test all functions together on range of inputs (runtime checks only)
unittest
{
    import std.algorithm : sort, reduce, equal;
    import std.range : chain;
    import std.typecons : staticIota;

    immutable maxRuntimeTest = 10_000;
    foreach (n; 1 .. maxRuntimeTest + 1)
    {
        assert(sort(chain(n.squareFreePrimeFactors, n.squarePrimeFactors)).equal(n.primeFactors));
        assert(n == n.squarePart * n.squareFreePart);
        assert(n == n.sqrtSquarePart * n.sqrtSquarePart * n.squareFreePart);
    }

    enum maxCompiletimeTest = 100;
    foreach (n; staticIota!(1, maxCompiletimeTest))
    {
        static assert(sort(chain(n.squareFreePrimeFactors,
            n.squarePrimeFactors)).equal(n.primeFactors));
        static assert(n == n.squarePart * n.squareFreePart);
        static assert(n == n.sqrtSquarePart * n.sqrtSquarePart * n.squareFreePart);
    }
}

private:

alias product = (a, b) => a * b;

int nextFactorFrom(int number, int factor) pure nothrow @nogc
{
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
    return nextFactor;
}

struct PrimeFactorsRange
{
    int remainingPart;
    int currentPart;

    this(int num) pure nothrow @nogc
    {
        assert(num > 0);

        remainingPart = num;
        currentPart = 1;
        popFront(); // Set currentPart to first prime factor
    }

    @property int front() const pure nothrow @nogc
    {
        assert(!empty); // Must have prime factor remaining
        return currentPart;
    }

    @property bool empty() const pure nothrow @nogc
    {
        return remainingPart == 1;
    }

    void popFront() pure nothrow @nogc
    {
        remainingPart /= currentPart;
        currentPart = remainingPart.nextFactorFrom(currentPart);
    }
}
