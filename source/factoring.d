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


version(unittest)
{
    import std.typecons : tuple;
    enum factorLists = tuple(
        [2], [3], [7919],[2, 7], [3, 11], [5, 7529], [2, 2], [13, 13], [7529, 7529],
        [2, 3, 5, 7, 11, 13, 17], [41, 41, 53, 67], [2, 3, 5, 11, 101, 101],[2, 2, 2, 11, 17],
        [2, 2, 2, 11, 11, 17],[2, 11, 4259, 4259], [2, 2, 2, 2, 2, 2, 2, 2],
        [5, 5, 5, 5, 5, 5, 5], [3, 3, 3, 3, 7, 7, 7]);
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

    foreach(factors; factorLists)
    {
        import std.algorithm : reduce, equal, uniq, count, filter, sort, setDifference;
        import std.range : stride;

        alias oddPower = (factor) => (factors.count(factor) % 2 == 1);
        alias product = (a,b) => a*b;
        immutable number = reduce!product(1, factors);

        assert(number.primeFactors.equal(factors));
        static assert(number.primeFactors.equal(factors));

        assert(number.squareFreePrimeFactors.equal(
            sort(factors).uniq.filter!oddPower));
        static assert(number.squareFreePrimeFactors.equal(
            sort(factors).uniq.filter!oddPower));

        assert(number.squarePrimeFactors.equal(
            sort(factors).setDifference(number.squareFreePrimeFactors)));
        static assert(number.squarePrimeFactors.equal(
            sort(factors).setDifference(number.squareFreePrimeFactors)));

        assert(number.sqrtSquarePrimeFactors.equal(
            number.squarePrimeFactors.stride(2)));
        static assert(number.sqrtSquarePrimeFactors.equal(
            number.squarePrimeFactors.stride(2)));

        assert(number.squareFreePart 
            == reduce!product(1, number.squareFreePrimeFactors));
        static assert(number.squareFreePart 
            == reduce!product(1, number.squareFreePrimeFactors));

        assert(number.squarePart 
            == reduce!product(1, number.squarePrimeFactors));
        static assert(number.squarePart 
            == reduce!product(1, number.squarePrimeFactors));

        assert(number.sqrtSquarePart 
            == reduce!product(1, number.sqrtSquarePrimeFactors));
        static assert(number.sqrtSquarePart 
            == reduce!product(1, number.sqrtSquarePrimeFactors));

        // TO DO: Why is this so slow?        
        // assert(number.primeFactorsRange.equal(factors));

        // TO DO: Why is this slow and use so much memory?        
        // static assert(number.primeFactorsRange.equal(factors));
    }
}

int[] squareFreePrimeFactors(int num)
{
    import std.algorithm : uniq, count;

    assert(num > 0);
    int[] factors = num.primeFactors;
    int[] answer;
    foreach (prime; factors.uniq)
    {
        if (count(factors, prime) % 2 == 1)
        {
            answer ~= prime;
        }
    }
    return answer;
}

int[] squarePrimeFactors(int num)
{
    import std.algorithm : uniq, count;

    assert(num > 0);
    int[] factors = num.primeFactors;
    int[] answer;
    foreach (prime; factors.uniq)
    {
        foreach (i; 0 .. 2 * (factors.count(prime) / 2))
        {
            answer ~= prime;
        }
    }
    return answer;
}

int[] sqrtSquarePrimeFactors(int num)
{
    import std.algorithm : uniq, count;

    assert(num > 0);
    int[] factors = num.primeFactors;
    int[] answer;
    foreach (prime; factors.uniq)
    {
        foreach (i; 0 .. factors.count(prime) / 2)
        {
            answer ~= prime;
        }
    }
    return answer;
}

int squareFreePart(int num)
{
    assert(num > 0);
    int answer = 1;
    foreach (prime; num.squareFreePrimeFactors)
    {
        answer *= prime;
    }
    return answer;
}

int squarePart(int num)
{
    assert(num > 0);
    int answer = 1;
    foreach (prime; num.squarePrimeFactors)
    {
        answer *= prime;
    }
    return answer;
}

int sqrtSquarePart(int num)
{
    assert(num > 0);
    int answer = 1;
    foreach (prime; num.sqrtSquarePrimeFactors)
    {
        answer *= prime;
    }
    return answer;
}

// Test all functions together on range of inputs (runtime checks only)
unittest
{
    import std.algorithm : sort, reduce, equal;
    import std.range : chain;

    alias product = (a,b) => a*b;

    immutable max_test = 10_000;
    foreach (n; 1 .. max_test + 1)
    {
        assert(sort(chain(n.squareFreePrimeFactors, n.squarePrimeFactors))
            .equal(n.primeFactors));

        assert(reduce!product(1, n.squareFreePrimeFactors) == n.squareFreePart);
        assert(reduce!product(1, n.squarePrimeFactors) == n.squarePart);
        assert(reduce!product(1, n.sqrtSquarePrimeFactors) == n.sqrtSquarePart);
        assert(n == n.squarePart * n.squareFreePart);
        assert(n == n.sqrtSquarePart * n.sqrtSquarePart * n.squareFreePart);
    }
}
