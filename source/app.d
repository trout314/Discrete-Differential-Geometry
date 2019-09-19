version(unittest) {} else {

import std.stdio : write, writeln, writefln;
void main(string[] args)
{
    import std.getopt : getopt;
    string task;
    getopt(args, "task", &task);

    if(task == "sampler")
    {
        //------------- SAMPLE MANIFOLDS BY OBJECTIVE  -----------------------
        import manifold : loadManifold, Manifold, findProblems, saveTo, saveEdgeGraphTo;
        import manifold_examples : standardSphere;
        import sampler : Parameters, sample, Sampler;
        import std.algorithm : each, findSplit;
        import std.conv : to;
        import std.getopt : getopt, defaultGetoptPrinter;
        import std.stdio : write, writeln, writefln, File;
        import std.datetime.stopwatch : StopWatch;
        import std.format : format;
        import std.range : empty;
        import std.uuid : randomUUID;
        import utility : flatDegreeInDim;

        enum dim = 3;

        Parameters params;
        with (params)
        {
            saveFilePrefix = "test";

            numFacetsTarget = 100;
            hingeDegreeTarget = flatDegreeInDim[3];
            
            numFacetsCoef = 0.01;
            numHingesCoef = 0.0;
            hingeDegreeVarCoef = 0.0;
            cd3DegVarCoef = 0.0;
            maxSweeps = 10;

            // Time increment (in units of sweeps) used for finer-graned intervals
            dt = 0.1;
            dtPerHistory = 50;
            dtPerFileReport = 50;
            dtPerSave = 250;

            triesPerStdoutReport = 200;

            useHingeMoves = true;
            disableGC = true;
            triesPerCollect = 50;
        }

        auto runID = randomUUID();

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

        params.saveFilePrefix = mfdFile.findSplit(".mfd")[0];
        if (params.saveFilePrefix.empty)
        {
            params.saveFilePrefix = runID.to!string[0..8];
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
        s.setParameters(params);
        s.sample;
    }
    else if (task=="get_edge_graph")
    {
        //--------------- GET EDGE GRAPH FROM MANIFOLD ------------------------
        enum dim = 3;

        import std.getopt : getopt, defaultGetoptPrinter;
        import std.stdio : write, writeln, writefln;
        import std.algorithm : findSplitBefore, findSplit;
        import std.range : empty;
        import manifold : loadManifold, saveEdgeGraphTo;
        import std.datetime.stopwatch : StopWatch;
        import std.format : format;

        // TO DO: add option to compute dual graph, filter nodes by
        // by various manifold properties etc

        string[] manifoldFiles;
        string[] outputGraphFiles;
        auto helpInformation = getopt(args,
            "manifoldFile|m", &manifoldFiles,
            "outputGraphFile|o", &outputGraphFiles);

        if (helpInformation.helpWanted)
        {
            // TO DO: some actual help info here!
            defaultGetoptPrinter("TO DO: Help info",
            helpInformation.options);
        }

        foreach(indx, mfdFile; manifoldFiles)
        {
            string outputFile;
            if(outputGraphFiles.empty)
            {
                outputFile = mfdFile.findSplitBefore(".mfd")[0]
                    ~ ".edge_graph";
            }
            else
            {
                outputFile = outputGraphFiles[indx];
            }

            StopWatch timer;
            timer.start;
            "Loading manifold file: %s".writefln(mfdFile);
            "... ".write;
            auto mfd = loadManifold!dim(mfdFile);
            "done. Took %s.".format(timer.peek)
                .findSplit(",")[0].findSplit("and")[0].writeln;
            
            timer.reset;

            "Creating edge graph file: %s".writefln(outputFile);
            "... ".write;
            mfd.saveEdgeGraphTo(outputFile);
            "done. Took %s.".format(timer.peek)
                .findSplit(",")[0].findSplit("and")[0].writeln;
        }
    }
    else if (task=="problem_check")
    {
        //------------- CHECK MANIFOLD FILE FOR PROBLEMS -----------------------
        enum dim = 3;

        import std.algorithm : each, findSplit;
        import std.getopt : getopt, defaultGetoptPrinter;
        import std.stdio : write, writeln, writefln;

        import manifold : loadManifold, findProblems;
        import std.datetime.stopwatch : StopWatch;
        import std.format : format;

        // TO DO: add option to compute dual graph, filter nodes by
        // by various manifold properties etc

        string[] manifoldFiles;
        auto helpInformation = getopt(args,
            "manifoldFile|m", &manifoldFiles);
        if (helpInformation.helpWanted)
        {
            // TO DO: some actual help info here!
            defaultGetoptPrinter("TO DO: Help info",
            helpInformation.options);
        }

        foreach(fileName; manifoldFiles)
        {
            auto mfd = loadManifold!dim(fileName);
            "checking input manifold file: %s".writefln(fileName);
            "... ".write;
            StopWatch timer;
            timer.start;
            mfd.findProblems.each!writeln;
            "done. Took %s.".format(timer.peek)
                .findSplit(",")[0].findSplit("and")[0].writeln;
        }
    }
    else if (task=="test_pachner_moves")
    {    
        import std.stdio : writeln;
        import moves : Move;
        import manifold : Manifold, computePachnerMoves, findCoCenter, doPachner, computeMBPMoves;
        import manifold_examples : standardSphere;
        import std.range : enumerate, iota;

        // trigonal bi-pyramid.
        auto m = standardSphere!2;

        auto report = () {
            writeln("   current manifold: ", m);
            "   moves:".writeln;
            foreach(indx, move; m.moves.enumerate)
            {
                write("      ", indx, ": ", move);
                if (m.contains(move.coCenter))
                {
                    writeln(", blocked");
                }
                else
                {
                    writeln;
                }
            }
            foreach(cen, indx; m.indxOfCenter)
            {
                writeln("      ", cen, ": ", indx);            
            }
            "   indicesOfCoCenter:".writeln;
            foreach(coCen, indices; m.indicesOfCoCenter)
            {
                writeln("      ", coCen, ": ", indices);            
            }
            writeln("   num valid moves: ", m.numValidMoves);

            import unit_threaded : shouldBeSameSetAs;
            import std.algorithm : each, sort;
            import std.range : array;

            m.moves.array.sort.each!writeln;
            "-------".writeln;
            m.computeMBPMoves.array.sort.each!writeln;
            m.moves.shouldBeSameSetAs(m.computeMBPMoves.dup.sort);
        };

        "starting out with standard 2-sphere... ".writeln;
        report();

        "m.doPachner([1,2,3], [4]);".writeln;
        m.doPachner([1,2,3], [4]);
        report();

        "m.doPachner([1,2], [0,4]);".writeln;
        m.doPachner([1,2], [0,4]);
        report();

        "m.doPachner([0,4], [1,2])".writeln;
        m.doPachner([0,4], [1,2]);
        report();    

        "m.doPachner([4], [1,2,3]);".writeln;
        m.doPachner([4], [1,2,3]);
        report();

        // "m.doPachner([1,2,3], [4]);".writeln;
        // m.doPachner([1,2,3], [4]);
        // report();

        // "m.doPachner([1,2,3], [5]);".writeln;
        // m.doPachner([1,2,4], [5]);
        // report();

    }
    else
    {
        import std.stdio : writeln;
        writeln("unrecognized task: " ~ task);
    }
}}