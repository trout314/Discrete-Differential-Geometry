import std.algorithm;

public int[] primeFactors(int num) pure nothrow
{
    assert(num > 0);
    int n = num;
    int[] answer = [];

    if (n < 2)
    {
        return answer;
    }

    int divisor = 2;
    while (divisor < n)
    {
        if (n % divisor == 0)
        {
            answer ~= divisor;
            n /= divisor;
        }
        else
        {
            if (divisor == 2)
            {
                divisor = 3;
            }
            else
            {
                divisor += 2;
            }
        }
    }
    answer ~= divisor;
    return answer;
}

unittest
{
    static assert(primeFactors(1) == []);
    static assert(primeFactors(2) == [2]);
    static assert(primeFactors(3) == [3]);
    static assert(primeFactors(7919) == [7919]);

    static assert(primeFactors(2 * 7) == [2, 7]);
    static assert(primeFactors(3 * 11) == [3, 11]);
    static assert(primeFactors(5 * 7529) == [5, 7529]);
    static assert(primeFactors(2 * 2) == [2, 2]);
    static assert(primeFactors(13 * 13) == [13, 13]);
    static assert(primeFactors(7529 * 7529) == [7529, 7529]);

    static assert(primeFactors(2 * 3 * 5 * 7 * 11 * 13 * 17) == [2, 3, 5, 7, 11, 13, 17]);
    static assert(primeFactors(41 * 41 * 53 * 67) == [41, 41, 53, 67]);
    static assert(primeFactors(2 * 3 * 5 * 11 * 101 * 101) == [2, 3, 5, 11, 101, 101]);
    static assert(primeFactors(2 * 2 * 2 * 11 * 17) == [2, 2, 2, 11, 17]);
    static assert(primeFactors(2 * 2 * 2 * 11 * 11 * 17) == [2, 2, 2, 11, 11, 17]);
    static assert(primeFactors(2 * 11 * 4259 * 4259) == [2, 11, 4259, 4259]);
    static assert(primeFactors(2 * 2 * 2 * 2 * 2 * 2 * 2 * 2) == [2, 2, 2, 2, 2, 2, 2, 2]);
    static assert(primeFactors(5 * 5 * 5 * 5 * 5 * 5 * 5) == [5, 5, 5, 5, 5, 5, 5]);
    static assert(primeFactors(3 * 3 * 3 * 3 * 7 * 7 * 7) == [3, 3, 3, 3, 7, 7, 7]);

    assert(primeFactors(1) == []);
    assert(primeFactors(2) == [2]);
    assert(primeFactors(3) == [3]);
    assert(primeFactors(7919) == [7919]);

    assert(primeFactors(2 * 7) == [2, 7]);
    assert(primeFactors(3 * 11) == [3, 11]);
    assert(primeFactors(5 * 7529) == [5, 7529]);
    assert(primeFactors(2 * 2) == [2, 2]);
    assert(primeFactors(13 * 13) == [13, 13]);
    assert(primeFactors(7529 * 7529) == [7529, 7529]);

    assert(primeFactors(2 * 3 * 5 * 7 * 11 * 13 * 17) == [2, 3, 5, 7, 11, 13, 17]);
    assert(primeFactors(41 * 41 * 53 * 67) == [41, 41, 53, 67]);
    assert(primeFactors(2 * 3 * 5 * 11 * 101 * 101) == [2, 3, 5, 11, 101, 101]);
    assert(primeFactors(2 * 2 * 2 * 11 * 17) == [2, 2, 2, 11, 17]);
    assert(primeFactors(2 * 2 * 2 * 11 * 11 * 17) == [2, 2, 2, 11, 11, 17]);
    assert(primeFactors(2 * 11 * 4259 * 4259) == [2, 11, 4259, 4259]);
    assert(primeFactors(2 * 2 * 2 * 2 * 2 * 2 * 2 * 2) == [2, 2, 2, 2, 2, 2, 2, 2]);
    assert(primeFactors(5 * 5 * 5 * 5 * 5 * 5 * 5) == [5, 5, 5, 5, 5, 5, 5]);
    assert(primeFactors(3 * 3 * 3 * 3 * 7 * 7 * 7) == [3, 3, 3, 3, 7, 7, 7]);
}

int[] squareFreePrimeFactors(in int num)
{
    assert(num > 0);
    int[] factors = primeFactors(num);
    auto uniqueFactors = uniq(factors);
    int[] answer;
    foreach (prime; uniqueFactors)
    {
        if (count(factors, prime) % 2 == 1)
        {
            answer ~= prime;
        }
    }
    return answer;
}

unittest
{
    static assert(squareFreePrimeFactors(1) == []);
    static assert(squareFreePrimeFactors(2) == [2]);
    static assert(squareFreePrimeFactors(3) == [3]);
    static assert(squareFreePrimeFactors(7919) == [7919]);

    static assert(squareFreePrimeFactors(2 * 7) == [2, 7]);
    static assert(squareFreePrimeFactors(3 * 11) == [3, 11]);
    static assert(squareFreePrimeFactors(5 * 7529) == [5, 7529]);
    static assert(squareFreePrimeFactors(2 * 2) == []);
    static assert(squareFreePrimeFactors(13 * 13) == []);
    static assert(squareFreePrimeFactors(7529 * 7529) == []);

    static assert(squareFreePrimeFactors(2 * 3 * 5 * 7 * 11 * 13 * 17) == [2, 3, 5, 7, 11, 13, 17]);
    static assert(squareFreePrimeFactors(41 * 41 * 53 * 67) == [53, 67]);
    static assert(squareFreePrimeFactors(2 * 3 * 5 * 11 * 101 * 101) == [2, 3, 5, 11]);
    static assert(squareFreePrimeFactors(2 * 2 * 2 * 11 * 17) == [2, 11, 17]);
    static assert(squareFreePrimeFactors(2 * 2 * 2 * 11 * 11 * 17) == [2, 17]);
    static assert(squareFreePrimeFactors(2 * 11 * 4259 * 4259) == [2, 11]);
    static assert(squareFreePrimeFactors(2 * 2 * 2 * 2 * 2 * 2 * 2 * 2) == []);
    static assert(squareFreePrimeFactors(5 * 5 * 5 * 5 * 5 * 5 * 5) == [5]);
    static assert(squareFreePrimeFactors(3 * 3 * 3 * 3 * 7 * 7 * 7) == [7]);

    assert(squareFreePrimeFactors(1) == []);
    assert(squareFreePrimeFactors(2) == [2]);
    assert(squareFreePrimeFactors(3) == [3]);
    assert(squareFreePrimeFactors(7919) == [7919]);

    assert(squareFreePrimeFactors(2 * 7) == [2, 7]);
    assert(squareFreePrimeFactors(3 * 11) == [3, 11]);
    assert(squareFreePrimeFactors(5 * 7529) == [5, 7529]);
    assert(squareFreePrimeFactors(2 * 2) == []);
    assert(squareFreePrimeFactors(13 * 13) == []);
    assert(squareFreePrimeFactors(7529 * 7529) == []);

    assert(squareFreePrimeFactors(2 * 3 * 5 * 7 * 11 * 13 * 17) == [2, 3, 5, 7, 11, 13, 17]);
    assert(squareFreePrimeFactors(41 * 41 * 53 * 67) == [53, 67]);
    assert(squareFreePrimeFactors(2 * 3 * 5 * 11 * 101 * 101) == [2, 3, 5, 11]);
    assert(squareFreePrimeFactors(2 * 2 * 2 * 11 * 17) == [2, 11, 17]);
    assert(squareFreePrimeFactors(2 * 2 * 2 * 11 * 11 * 17) == [2, 17]);
    assert(squareFreePrimeFactors(2 * 11 * 4259 * 4259) == [2, 11]);
    assert(squareFreePrimeFactors(2 * 2 * 2 * 2 * 2 * 2 * 2 * 2) == []);
    assert(squareFreePrimeFactors(5 * 5 * 5 * 5 * 5 * 5 * 5) == [5]);
    assert(squareFreePrimeFactors(3 * 3 * 3 * 3 * 7 * 7 * 7) == [7]);
}

int[] squarePrimeFactors(in int num)
{
    assert(num > 0);
    int[] factors = primeFactors(num);
    int[] answer;
    foreach (prime; uniq(factors))
    {
        foreach (i; 0 .. 2 * (count(factors, prime) / 2))
        {
            answer ~= prime;
        }
    }
    return answer;
}

unittest
{
    assert(squarePrimeFactors(1) == []);
    assert(squarePrimeFactors(2) == []);
    assert(squarePrimeFactors(3) == []);
    assert(squarePrimeFactors(7919) == []);

    assert(squarePrimeFactors(2 * 7) == []);
    assert(squarePrimeFactors(3 * 11) == []);
    assert(squarePrimeFactors(5 * 7529) == []);
    assert(squarePrimeFactors(2 * 2) == [2, 2]);
    assert(squarePrimeFactors(13 * 13) == [13, 13]);
    assert(squarePrimeFactors(7529 * 7529) == [7529, 7529]);

    assert(squarePrimeFactors(2 * 3 * 5 * 7 * 11 * 13 * 17) == []);
    assert(squarePrimeFactors(41 * 41 * 53 * 67) == [41, 41]);
    assert(squarePrimeFactors(2 * 3 * 5 * 11 * 101 * 101) == [101, 101]);
    assert(squarePrimeFactors(2 * 2 * 2 * 11 * 17) == [2, 2]);
    assert(squarePrimeFactors(2 * 2 * 2 * 11 * 11 * 17) == [2, 2, 11, 11]);
    assert(squarePrimeFactors(2 * 11 * 4259 * 4259) == [4259, 4259]);
    assert(squarePrimeFactors(2 * 2 * 2 * 2 * 2 * 2 * 2 * 2) == [2, 2, 2, 2, 2, 2, 2, 2]);
    assert(squarePrimeFactors(5 * 5 * 5 * 5 * 5 * 5 * 5) == [5, 5, 5, 5, 5, 5]);
    assert(squarePrimeFactors(3 * 3 * 3 * 3 * 7 * 7 * 7) == [3, 3, 3, 3, 7, 7]);

    static assert(squarePrimeFactors(1) == []);
    static assert(squarePrimeFactors(2) == []);
    static assert(squarePrimeFactors(3) == []);
    static assert(squarePrimeFactors(7919) == []);

    static assert(squarePrimeFactors(2 * 7) == []);
    static assert(squarePrimeFactors(3 * 11) == []);
    static assert(squarePrimeFactors(5 * 7529) == []);
    static assert(squarePrimeFactors(2 * 2) == [2, 2]);
    static assert(squarePrimeFactors(13 * 13) == [13, 13]);
    static assert(squarePrimeFactors(7529 * 7529) == [7529, 7529]);

    static assert(squarePrimeFactors(2 * 3 * 5 * 7 * 11 * 13 * 17) == []);
    static assert(squarePrimeFactors(41 * 41 * 53 * 67) == [41, 41]);
    static assert(squarePrimeFactors(2 * 3 * 5 * 11 * 101 * 101) == [101, 101]);
    static assert(squarePrimeFactors(2 * 2 * 2 * 11 * 17) == [2, 2]);
    static assert(squarePrimeFactors(2 * 2 * 2 * 11 * 11 * 17) == [2, 2, 11, 11]);
    static assert(squarePrimeFactors(2 * 11 * 4259 * 4259) == [4259, 4259]);
    static assert(squarePrimeFactors(2 * 2 * 2 * 2 * 2 * 2 * 2 * 2) == [2, 2, 2, 2, 2, 2, 2, 2]);
    static assert(squarePrimeFactors(5 * 5 * 5 * 5 * 5 * 5 * 5) == [5, 5, 5, 5, 5, 5]);
    static assert(squarePrimeFactors(3 * 3 * 3 * 3 * 7 * 7 * 7) == [3, 3, 3, 3, 7, 7]);
}

int[] sqrtSquarePrimeFactors(in int num)
{
    assert(num > 0);
    int[] factors = primeFactors(num);
    int[] answer;
    foreach (prime; uniq(factors))
    {
        foreach (i; 0 .. count(factors, prime) / 2)
        {
            answer ~= prime;
        }
    }
    return answer;
}

unittest
{
    assert(sqrtSquarePrimeFactors(1) == []);
    assert(sqrtSquarePrimeFactors(2) == []);
    assert(sqrtSquarePrimeFactors(3) == []);
    assert(sqrtSquarePrimeFactors(7919) == []);

    assert(sqrtSquarePrimeFactors(2 * 7) == []);
    assert(sqrtSquarePrimeFactors(3 * 11) == []);
    assert(sqrtSquarePrimeFactors(5 * 7529) == []);
    assert(sqrtSquarePrimeFactors(2 * 2) == [2]);
    assert(sqrtSquarePrimeFactors(13 * 13) == [13]);
    assert(sqrtSquarePrimeFactors(7529 * 7529) == [7529]);

    assert(sqrtSquarePrimeFactors(2 * 3 * 5 * 7 * 11 * 13 * 17) == []);
    assert(sqrtSquarePrimeFactors(41 * 41 * 53 * 67) == [41]);
    assert(sqrtSquarePrimeFactors(2 * 3 * 5 * 11 * 101 * 101) == [101]);
    assert(sqrtSquarePrimeFactors(2 * 2 * 2 * 11 * 17) == [2]);
    assert(sqrtSquarePrimeFactors(2 * 2 * 2 * 11 * 11 * 17) == [2, 11]);
    assert(sqrtSquarePrimeFactors(2 * 11 * 4259 * 4259) == [4259]);
    assert(sqrtSquarePrimeFactors(2 * 2 * 2 * 2 * 2 * 2 * 2 * 2) == [2, 2, 2, 2]);
    assert(sqrtSquarePrimeFactors(5 * 5 * 5 * 5 * 5 * 5 * 5) == [5, 5, 5]);
    assert(sqrtSquarePrimeFactors(3 * 3 * 3 * 3 * 7 * 7 * 7) == [3, 3, 7]);

    static assert(sqrtSquarePrimeFactors(1) == []);
    static assert(sqrtSquarePrimeFactors(2) == []);
    static assert(sqrtSquarePrimeFactors(3) == []);
    static assert(sqrtSquarePrimeFactors(7919) == []);

    static assert(sqrtSquarePrimeFactors(2 * 7) == []);
    static assert(sqrtSquarePrimeFactors(3 * 11) == []);
    static assert(sqrtSquarePrimeFactors(5 * 7529) == []);
    static assert(sqrtSquarePrimeFactors(2 * 2) == [2]);
    static assert(sqrtSquarePrimeFactors(13 * 13) == [13]);
    static assert(sqrtSquarePrimeFactors(7529 * 7529) == [7529]);

    static assert(sqrtSquarePrimeFactors(2 * 3 * 5 * 7 * 11 * 13 * 17) == []);
    static assert(sqrtSquarePrimeFactors(41 * 41 * 53 * 67) == [41]);
    static assert(sqrtSquarePrimeFactors(2 * 3 * 5 * 11 * 101 * 101) == [101]);
    static assert(sqrtSquarePrimeFactors(2 * 2 * 2 * 11 * 17) == [2]);
    static assert(sqrtSquarePrimeFactors(2 * 2 * 2 * 11 * 11 * 17) == [2, 11]);
    static assert(sqrtSquarePrimeFactors(2 * 11 * 4259 * 4259) == [4259]);
    static assert(sqrtSquarePrimeFactors(2 * 2 * 2 * 2 * 2 * 2 * 2 * 2) == [2, 2, 2, 2]);
    static assert(sqrtSquarePrimeFactors(5 * 5 * 5 * 5 * 5 * 5 * 5) == [5, 5, 5]);
    static assert(sqrtSquarePrimeFactors(3 * 3 * 3 * 3 * 7 * 7 * 7) == [3, 3, 7]);
}

int squareFreePart(in int num)
{
    assert(num > 0);
    int answer = 1;
    foreach (prime; squareFreePrimeFactors(num))
    {
        answer *= prime;
    }
    return answer;
}

unittest
{
    assert(squareFreePart(1) == 1);
    assert(squareFreePart(2) == 2);
    assert(squareFreePart(3) == 3);
    assert(squareFreePart(7919) == 7919);

    assert(squareFreePart(2 * 7) == 2 * 7);
    assert(squareFreePart(3 * 11) == 3 * 11);
    assert(squareFreePart(5 * 7529) == 5 * 7529);
    assert(squareFreePart(2 * 2) == 1);
    assert(squareFreePart(13 * 13) == 1);
    assert(squareFreePart(7529 * 7529) == 1);

    assert(squareFreePart(2 * 3 * 5 * 7 * 11 * 13 * 17) == 2 * 3 * 5 * 7 * 11 * 13 * 17);
    assert(squareFreePart(41 * 41 * 53 * 67) == 53 * 67);
    assert(squareFreePart(2 * 3 * 5 * 11 * 101 * 101) == 2 * 3 * 5 * 11);
    assert(squareFreePart(2 * 2 * 2 * 11 * 17) == 2 * 11 * 17);
    assert(squareFreePart(2 * 2 * 2 * 11 * 11 * 17) == 2 * 17);
    assert(squareFreePart(2 * 11 * 4259 * 4259) == 2 * 11);
    assert(squareFreePart(2 * 2 * 2 * 2 * 2 * 2 * 2 * 2) == 1);
    assert(squareFreePart(5 * 5 * 5 * 5 * 5 * 5 * 5) == 5);
    assert(squareFreePart(3 * 3 * 3 * 3 * 7 * 7 * 7) == 7);

    static assert(squareFreePart(1) == 1);
    static assert(squareFreePart(2) == 2);
    static assert(squareFreePart(3) == 3);
    static assert(squareFreePart(7919) == 7919);

    static assert(squareFreePart(2 * 7) == 2 * 7);
    static assert(squareFreePart(3 * 11) == 3 * 11);
    static assert(squareFreePart(5 * 7529) == 5 * 7529);
    static assert(squareFreePart(2 * 2) == 1);
    static assert(squareFreePart(13 * 13) == 1);
    static assert(squareFreePart(7529 * 7529) == 1);

    static assert(squareFreePart(2 * 3 * 5 * 7 * 11 * 13 * 17) == 2 * 3 * 5 * 7 * 11 * 13 * 17);
    static assert(squareFreePart(41 * 41 * 53 * 67) == 53 * 67);
    static assert(squareFreePart(2 * 3 * 5 * 11 * 101 * 101) == 2 * 3 * 5 * 11);
    static assert(squareFreePart(2 * 2 * 2 * 11 * 17) == 2 * 11 * 17);
    static assert(squareFreePart(2 * 2 * 2 * 11 * 11 * 17) == 2 * 17);
    static assert(squareFreePart(2 * 11 * 4259 * 4259) == 2 * 11);
    static assert(squareFreePart(2 * 2 * 2 * 2 * 2 * 2 * 2 * 2) == 1);
    static assert(squareFreePart(5 * 5 * 5 * 5 * 5 * 5 * 5) == 5);
    static assert(squareFreePart(3 * 3 * 3 * 3 * 7 * 7 * 7) == 7);

}

int squarePart(in int num)
{
    assert(num > 0);
    int answer = 1;
    foreach (prime; squarePrimeFactors(num))
    {
        answer *= prime;
    }
    return answer;
}

unittest
{
    assert(squarePart(1) == 1);
    assert(squarePart(2) == 1);
    assert(squarePart(3) == 1);
    assert(squarePart(7919) == 1);

    assert(squarePart(2 * 7) == 1);
    assert(squarePart(3 * 11) == 1);
    assert(squarePart(5 * 7529) == 1);
    assert(squarePart(2 * 2) == 2 * 2);
    assert(squarePart(13 * 13) == 13 * 13);
    assert(squarePart(7529 * 7529) == 7529 * 7529);

    assert(squarePart(2 * 3 * 5 * 7 * 11 * 13 * 17) == 1);
    assert(squarePart(41 * 41 * 53 * 67) == 41 * 41);
    assert(squarePart(2 * 3 * 5 * 11 * 101 * 101) == 101 * 101);
    assert(squarePart(2 * 2 * 2 * 11 * 17) == 2 * 2);
    assert(squarePart(2 * 2 * 2 * 11 * 11 * 17) == 2 * 2 * 11 * 11);
    assert(squarePart(2 * 11 * 4259 * 4259) == 4259 * 4259);
    assert(squarePart(2 * 2 * 2 * 2 * 2 * 2 * 2 * 2) == 2 * 2 * 2 * 2 * 2 * 2 * 2 * 2);
    assert(squarePart(5 * 5 * 5 * 5 * 5 * 5 * 5) == 5 * 5 * 5 * 5 * 5 * 5);
    assert(squarePart(3 * 3 * 3 * 3 * 7 * 7 * 7) == 3 * 3 * 3 * 3 * 7 * 7);

    static assert(squarePart(1) == 1);
    static assert(squarePart(2) == 1);
    static assert(squarePart(3) == 1);
    static assert(squarePart(7919) == 1);

    static assert(squarePart(2 * 7) == 1);
    static assert(squarePart(3 * 11) == 1);
    static assert(squarePart(5 * 7529) == 1);
    static assert(squarePart(2 * 2) == 2 * 2);
    static assert(squarePart(13 * 13) == 13 * 13);
    static assert(squarePart(7529 * 7529) == 7529 * 7529);

    static assert(squarePart(2 * 3 * 5 * 7 * 11 * 13 * 17) == 1);
    static assert(squarePart(41 * 41 * 53 * 67) == 41 * 41);
    static assert(squarePart(2 * 3 * 5 * 11 * 101 * 101) == 101 * 101);
    static assert(squarePart(2 * 2 * 2 * 11 * 17) == 2 * 2);
    static assert(squarePart(2 * 2 * 2 * 11 * 11 * 17) == 2 * 2 * 11 * 11);
    static assert(squarePart(2 * 11 * 4259 * 4259) == 4259 * 4259);
    static assert(squarePart(2 * 2 * 2 * 2 * 2 * 2 * 2 * 2) == 2 * 2 * 2 * 2 * 2 * 2 * 2 * 2);
    static assert(squarePart(5 * 5 * 5 * 5 * 5 * 5 * 5) == 5 * 5 * 5 * 5 * 5 * 5);
    static assert(squarePart(3 * 3 * 3 * 3 * 7 * 7 * 7) == 3 * 3 * 3 * 3 * 7 * 7);
}

int sqrtSquarePart(in int num)
{
    assert(num > 0);
    int answer = 1;
    foreach (prime; sqrtSquarePrimeFactors(num))
    {
        answer *= prime;
    }
    return answer;
}

unittest
{
    assert(sqrtSquarePart(1) == 1);
    assert(sqrtSquarePart(2) == 1);
    assert(sqrtSquarePart(3) == 1);
    assert(sqrtSquarePart(7919) == 1);

    assert(sqrtSquarePart(2 * 7) == 1);
    assert(sqrtSquarePart(3 * 11) == 1);
    assert(sqrtSquarePart(5 * 7529) == 1);
    assert(sqrtSquarePart(2 * 2) == 2);
    assert(sqrtSquarePart(13 * 13) == 13);
    assert(sqrtSquarePart(7529 * 7529) == 7529);

    assert(sqrtSquarePart(2 * 3 * 5 * 7 * 11 * 13 * 17) == 1);
    assert(sqrtSquarePart(41 * 41 * 53 * 67) == 41);
    assert(sqrtSquarePart(2 * 3 * 5 * 11 * 101 * 101) == 101);
    assert(sqrtSquarePart(2 * 2 * 2 * 11 * 17) == 2);
    assert(sqrtSquarePart(2 * 2 * 2 * 11 * 11 * 17) == 2 * 11);
    assert(sqrtSquarePart(2 * 11 * 4259 * 4259) == 4259);
    assert(sqrtSquarePart(2 * 2 * 2 * 2 * 2 * 2 * 2 * 2) == 2 * 2 * 2 * 2);
    assert(sqrtSquarePart(5 * 5 * 5 * 5 * 5 * 5 * 5) == 5 * 5 * 5);
    assert(sqrtSquarePart(3 * 3 * 3 * 3 * 7 * 7 * 7) == 3 * 3 * 7);

    static assert(sqrtSquarePart(1) == 1);
    static assert(sqrtSquarePart(2) == 1);
    static assert(sqrtSquarePart(3) == 1);
    static assert(sqrtSquarePart(7919) == 1);

    static assert(sqrtSquarePart(2 * 7) == 1);
    static assert(sqrtSquarePart(3 * 11) == 1);
    static assert(sqrtSquarePart(5 * 7529) == 1);
    static assert(sqrtSquarePart(2 * 2) == 2);
    static assert(sqrtSquarePart(13 * 13) == 13);
    static assert(sqrtSquarePart(7529 * 7529) == 7529);

    static assert(sqrtSquarePart(2 * 3 * 5 * 7 * 11 * 13 * 17) == 1);
    static assert(sqrtSquarePart(41 * 41 * 53 * 67) == 41);
    static assert(sqrtSquarePart(2 * 3 * 5 * 11 * 101 * 101) == 101);
    static assert(sqrtSquarePart(2 * 2 * 2 * 11 * 17) == 2);
    static assert(sqrtSquarePart(2 * 2 * 2 * 11 * 11 * 17) == 2 * 11);
    static assert(sqrtSquarePart(2 * 11 * 4259 * 4259) == 4259);
    static assert(sqrtSquarePart(2 * 2 * 2 * 2 * 2 * 2 * 2 * 2) == 2 * 2 * 2 * 2);
    static assert(sqrtSquarePart(5 * 5 * 5 * 5 * 5 * 5 * 5) == 5 * 5 * 5);
    static assert(sqrtSquarePart(3 * 3 * 3 * 3 * 7 * 7 * 7) == 3 * 3 * 7);
}

// Test all functions together on range of inputs (runtime checks only)
unittest
{
    immutable max_test = 10_000;
    foreach (n; 1 .. max_test + 1)
    {
        auto all = squareFreePrimeFactors(n) ~ squarePrimeFactors(n);
        sort(all);
        assert(all == primeFactors(n));

        int productSqF = 1;
        foreach (prime; squareFreePrimeFactors(n))
        {
            productSqF *= prime;
        }
        assert(productSqF == squareFreePart(n));

        int productSq = 1;
        foreach (prime; squarePrimeFactors(n))
        {
            productSq *= prime;
        }
        assert(productSq == squarePart(n));

        int productSqrtSq = 1;
        foreach (prime; sqrtSquarePrimeFactors(n))
        {
            productSqrtSq *= prime;
        }
        assert(productSqrtSq == sqrtSquarePart(n));

        assert(n == squarePart(n) * squareFreePart(n));
        assert(n == sqrtSquarePart(n) * sqrtSquarePart(n) * squareFreePart(n));
    }
}
