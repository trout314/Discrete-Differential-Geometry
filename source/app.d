version(unittest) {}
else
{
    // --------------- GET EDGE GRAPH FROM MANIFOLD ------------------------
    // void main(string[] args)
    // {
    //     enum dim = 3;

    //     import std.getopt : getopt, defaultGetoptPrinter;
    //     import std.stdio : write, writeln, writefln;
    //     import std.algorithm : findSplitBefore, findSplit;
    //     import std.range : empty;
    //     import manifold : loadManifold, saveEdgeGraphTo;
    //     import std.datetime.stopwatch : StopWatch;
    //     import std.format : format;

    //     // TO DO: add option to compute dual graph, filter nodes by
    //     // by various manifold properties etc

    //     string[] manifoldFiles;
    //     string[] outputGraphFiles;
    //     auto helpInformation = getopt(args,
    //         "manifoldFile|m", &manifoldFiles,
    //         "outputGraphFile|o", &outputGraphFiles);

    //     if (helpInformation.helpWanted)
    //     {
    //         // TO DO: some actual help info here!
    //         defaultGetoptPrinter("TO DO: Help info",
    //         helpInformation.options);
    //     }

    //     foreach(indx, mfdFile; manifoldFiles)
    //     {
    //         string outputFile;
    //         if(outputGraphFiles.empty)
    //         {
    //             outputFile = mfdFile.findSplitBefore(".mfd")[0]
    //                 ~ ".edge_graph";
    //         }
    //         else
    //         {
    //             outputFile = outputGraphFiles[indx];
    //         }

    //         StopWatch timer;
    //         timer.start;
    //         "Loading manifold file: %s".writefln(mfdFile);
    //         "... ".write;
    //         auto mfd = loadManifold!dim(mfdFile);
    //         "done. Took %s.".format(timer.peek)
    //             .findSplit(",")[0].findSplit("and")[0].writeln;
            
    //         timer.reset;

    //         "Creating edge graph file: %s".writefln(outputFile);
    //         "... ".write;
    //         mfd.saveEdgeGraphTo(outputFile);
    //         "done. Took %s.".format(timer.peek)
    //             .findSplit(",")[0].findSplit("and")[0].writeln;
    //     }
    // }
    
    //------------- CHECK MANIFOLD FILE FOR PROBLEMS -----------------------
    // void main(string[] args)
    // {
    //     enum dim = 3;

    //     import std.algorithm : each, findSplit;
    //     import std.getopt : getopt, defaultGetoptPrinter;
    //     import std.stdio : write, writeln, writefln;

    //     import manifold : loadManifold, findProblems;
    //     import std.datetime.stopwatch : StopWatch;
    //     import std.format : format;

    //     // TO DO: add option to compute dual graph, filter nodes by
    //     // by various manifold properties etc

    //     string[] manifoldFiles;
    //     auto helpInformation = getopt(args,
    //         "manifoldFile|m", &manifoldFiles);
    //     if (helpInformation.helpWanted)
    //     {
    //         // TO DO: some actual help info here!
    //         defaultGetoptPrinter("TO DO: Help info",
    //         helpInformation.options);
    //     }

    //     foreach(fileName; manifoldFiles)
    //     {
    //         auto mfd = loadManifold!dim(fileName);
    //         "checking input manifold file: %s".writefln(fileName);
    //         "... ".write;
    //         StopWatch timer;
    //         timer.start;
    //         mfd.findProblems.each!writeln;
    //         "done. Took %s.".format(timer.peek)
    //             .findSplit(",")[0].findSplit("and")[0].writeln;
    //     }
    // }

    //------------- SAMPLE MANIFOLDS BY OBJECTIVE  -----------------------
    void main(string[] args)
    {
        enum dim = 3;

        import std.algorithm : each, findSplit;
        import std.getopt : getopt, defaultGetoptPrinter;
        import std.stdio : write, writeln, writefln, File;

        import manifold : loadManifold, Manifold, findProblems, saveTo, saveEdgeGraphTo, standardSphere;
        import sampler : sample, Sampler;
        import std.datetime.stopwatch : StopWatch;
        import std.format : format;
        import std.range : empty;

        // TO DO: add option to compute dual graph, filter nodes by
        // by various manifold properties etc

        string mfdFile;
        auto helpInformation = getopt(args,
            "manifoldFile|m", &mfdFile);
        if (helpInformation.helpWanted)
        {
            // TO DO: some actual help info here!
            defaultGetoptPrinter("TO DO: Help info",
            helpInformation.options);
        }

        StopWatch timer;
        timer.start;
        Manifold!dim mfd;
        if (!mfdFile.empty)
        {
            "Loading manifold file: %s".writefln(mfdFile);
            "... ".write;
            mfd = loadManifold!dim(mfdFile);
            "done. Took %s.".format(timer.peek)
                .findSplit(",")[0].findSplit("and")[0].writeln;
            timer.reset;
        }
        else
        {
            // TO DO: User chosen options from library of triangs
            mfd = standardSphere!dim;
        }

        auto s = Sampler!(int, dim)(mfd);
        s.sample;

        auto start = "8k_test";

        auto fileName = start ~ ".mfd";
        s.manifold.saveTo(fileName);
        auto saveFile = File(fileName, "a");
        saveFile.writeln;
        s.report(saveFile);
        s.manifold.saveEdgeGraphTo(start ~ ".edge_graph");

        foreach(d; 0 .. dim - 1)
        {
            "dim %s historic averages:".writefln(d);
            "  mean deg. = %s".writefln(
                s.historicMeanDeg(d));
            "  std. dev. = %s".writefln(
                s.historicStdDevDeg(d));
        }

    }
}
