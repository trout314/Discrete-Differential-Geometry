module applications.edge_graph;

import std.getopt;

import simplicial_complex, utility;

version (unittest) {} else {
int main(string[] args)
{
    string surfaceFileName;
    auto helpInformation = getopt(args, std.getopt.config.required,
        "simplicialComplex|s", &surfaceFileName);

    if (helpInformation.helpWanted)
    {
        // TO DO: some actual help info here!
        defaultGetoptPrinter("TO DO: Help info",
            helpInformation.options);
        return 0; // exit with return code zero (success)
    }

    auto surface = loadSimplicialComplex(surfaceFileName);

    foreach(vertex; surface.simplices(0))
    {
        dump!vertex;
    }




    dump!surface;

    return 0;
}}