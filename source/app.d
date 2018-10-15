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
        import sampler : Parameters, sample, Sampler;
        import std.datetime.stopwatch : StopWatch;
        import std.format : format;
        import std.range : empty;
        import std.uuid : randomUUID;
        import utility : flatDegreeInDim;

        Parameters params;

        with (params)
        {
            numFacetsTarget = 64_000;
            hingeDegreeTarget = flatDegreeInDim[3];
            
            numFacetsCoef = 0.01;
            numHingesCoef = 0.01;
            hingeDegreeVarCoef = 0.05;
            cd3DegVarCoef = hingeDegreeVarCoef / (22.8 / 5.1);
            maxSweeps = 100;

            // Time increment (in units of sweeps) used for finer-graned intervals
            dt = 0.1;
            dtPerHistory = 5;
            // dtPerFileReport = 20;
            // dtPerSave = 1000;

            triesPerStdoutReport = 50_000;

            useHingeMoves = true;
            disableGC = true;
            triesPerCollect = 500;
        }

        writeln(randomUUID());

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

        auto filePrefix = mfdFile.findSplit(".mfd")[0];
        params.saveFilePrefix = filePrefix ~ "_sample";

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
        s.setParameters(params);
        s.sample;

        s.manifold.saveTo(filePrefix ~ "_final.mfd");
        auto saveFile = File(filePrefix ~ "_final.mfd", "a");
        s.report(saveFile);
        s.manifold.saveEdgeGraphTo(filePrefix ~ "_final.edge_graph");

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
