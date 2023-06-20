module app_2;

version(unittest) {} else {

void main()
{
    import std.stdio : writeln;
    ("hello from " ~ __FILE__).writeln;
}
}