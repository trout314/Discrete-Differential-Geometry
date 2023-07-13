module applications.app_2;

version(unittest) {} else {
void main()
{
    import std.array : array;
    import std.range : enumerate;
    import std.stdio : writeln;

    auto data = ["A", "B", "C"];
    auto enumeratedData = data.enumerate;
    auto enumeratedDataArray = data.enumerate.array;
 
    foreach(index, value; enumeratedData)
    {
        writeln(index, "   ", value);
    }

    foreach(index, value; enumeratedDataArray)
    {
        writeln(index, "   ", value);
    }
}
}