module applications.dual_graph;

import std.datetime, std.format, std.getopt, std.range, std.stdio;

import manifold, utility;

version (unittest) {} else {
int main(string[] args)
{
    enum maxDim = 8;

    string manifoldFilename;
    string dualGraphFilename;
    int dim;
    auto helpInformation = getopt(args,
        std.getopt.config.required, "manifold|m", &manifoldFilename,
        std.getopt.config.required, "dimension|d", &dim,
        "dualGraph|g", &dualGraphFilename);

    if (helpInformation.helpWanted)
    {
        // TO DO: some actual help info here!
        defaultGetoptPrinter("TO DO: Help info",
            helpInformation.options);
        return 0; // exit with return code zero (success)
    }

    if (dualGraphFilename.empty)
    {
        dualGraphFilename = manifoldFilename.split(".")[0] ~ ".dual_graph";
    }

    static foreach(d; 2 .. maxDim)
    {{
        if (dim == d)
        {
            auto mfd = loadManifold!d(manifoldFilename);
            mfd.saveDualGraphTo(dualGraphFilename);
        }
    }}

    assert((dim <= maxDim) && (dim >= 2), "unsupported dimension for manifold");

    auto dualGraphFile = File(dualGraphFilename, "a");
    dualGraphFile.writeln(
        "# Dual graph created from file: %s on %s".format(manifoldFilename, Clock.currTime));

    return 0;
}}