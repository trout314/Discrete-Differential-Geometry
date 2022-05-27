module applications.app_1;

version(unittest) {} else {

void main()
{
    import std.stdio : writeln;
    ("hello from " ~ __FILE__).writeln;
}

}