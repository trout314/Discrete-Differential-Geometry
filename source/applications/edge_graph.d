module applications.edge_graph;

import std.datetime, std.format, std.getopt, std.range, std.stdio;

import simplicial_complex, utility;

version (unittest) {} else {
int main(string[] args)
{
    string simpCompFilename;
    string edgeGraphFilename;
    auto helpInformation = getopt(
        args,
        std.getopt.config.required, "simplicialComplex|s", &simpCompFilename,
        "edgeGraph|e", &edgeGraphFilename);

    if (helpInformation.helpWanted)
    {
        // TO DO: some actual help info here!
        defaultGetoptPrinter("TO DO: Help info",
            helpInformation.options);
        return 0; // exit with return code zero (success)
    }

    auto simpComp = loadSimplicialComplex(simpCompFilename);
    if (edgeGraphFilename.empty)
    {
        edgeGraphFilename = simpCompFilename.split(".")[0 .. $-1].join(".") ~ ".edge_graph";
    }

    simpComp.saveEdgeGraphTo(edgeGraphFilename);
    auto edgeGraphFile = File(edgeGraphFilename, "a");
    edgeGraphFile.writeln(
        "# Edge graph created from file: %s on %s".format(simpCompFilename, Clock.currTime));

    return 0;
}}