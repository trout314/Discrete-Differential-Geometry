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

alias product = (a,b) => a*b;

struct PrimeFactorsRange
{
    int remainingPart;
    int currentPart;

    this(int n) pure nothrow @nogc
    {
        remainingPart = n;
        currentPart = 1;
        popFront();
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

auto primeFactorsRange(int num) pure nothrow @nogc
{
    return PrimeFactorsRange(num);
}

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

unittest
{
    assert(primeFactors(1) == []);
    static assert(primeFactors(1) == []);
    
    assert(squareFreePrimeFactors(1) == []);
    static assert(squareFreePrimeFactors(1) == []);

    assert(squarePrimeFactors(1) == []);
    static assert(squarePrimeFactors(1) == []);

    assert(sqrtSquarePrimeFactors(1) == []);
    static assert(sqrtSquarePrimeFactors(1) == []);


    import std.typecons : tuple;
    enum testFactorLists = tuple(
        [2], [3], [7919],[2, 7], [3, 11], [5, 7529], [2, 2], [13, 13], [7529, 7529],
        [2, 3, 5, 7, 11, 13, 17], [41, 41, 53, 67], [2, 3, 5, 11, 101, 101],[2, 2, 2, 11, 17],
        [2, 2, 2, 11, 11, 17],[2, 11, 4259, 4259], [2, 2, 2, 2, 2, 2, 2, 2],
        [5, 5, 5, 5, 5, 5, 5], [3, 3, 3, 3, 7, 7, 7]);

    foreach(factors; testFactorLists)
    {
        import std.algorithm : reduce, equal;

        immutable number = 1.reduce!product(factors);
        assert(number.primeFactors.equal(factors));
        static assert(number.primeFactors.equal(factors));

        // TO DO: Why is this so slow?    
        // assert(number.primeFactorsRange.equal(factors));

        // TO DO: Why is this slow and use so much memory?        
        // static assert(number.primeFactorsRange.equal(factors));
    }
}

int[] squareFreePrimeFactors(int num)
{
    assert(num > 0);

    import std.algorithm : uniq, count, filter;
    import std.range : array;

    auto factors = num.primeFactors;
    alias oddPower = f => factors.count(f) % 2 == 1;
    return factors.uniq.filter!oddPower.array;
}

int[] squarePrimeFactors(int num)
{
    assert(num > 0);

    import std.algorithm : setDifference;
    import std.range : array;
    return num.primeFactors.setDifference(num.squareFreePrimeFactors).array;
}

int[] sqrtSquarePrimeFactors(int num)
{
    assert(num > 0);

    import std.range : array, stride;
    return num.squarePrimeFactors.stride(2).array;
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

    immutable max_test = 10_000;
    foreach (n; 1 .. max_test + 1)
    {
        assert(sort(chain(n.squareFreePrimeFactors, n.squarePrimeFactors))
            .equal(n.primeFactors));
        assert(n == n.squarePart * n.squareFreePart);
        assert(n == n.sqrtSquarePart * n.sqrtSquarePart * n.squareFreePart);
    }
}
