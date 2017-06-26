import rational : rational, Rational;
import std.bigint : BigInt;
import std.stdio : writefln, writeln;

enum maxTerm = 30;

Rational!(BigInt) getTerm(int termNumber) {
    auto addFactor = 2 * termNumber - 1;

    if(termNumber == maxTerm) {
        return rational(BigInt(addFactor));
    }

    auto termNumberSquared = BigInt(termNumber * termNumber);
    auto continued = termNumberSquared / getTerm(termNumber + 1);

    continued += addFactor;
    return continued;
}

void main() {

    auto pi = rational(BigInt(4)) / getTerm(1);

    // Display the result in rational form.
    writeln(pi);

    // Display the decimal equivalent, which is accurate to 18 decimal places.
    writefln("%.18f", cast(real) pi);
}