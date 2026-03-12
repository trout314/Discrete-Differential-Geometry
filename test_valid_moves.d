/// Standalone test for incremental valid move tracking.
/// Build with: dub build --config=test_valid_moves
module test_valid_moves;

version (unittest) {} else
{

import std.algorithm, std.array, std.conv, std.random, std.stdio, std.format;
import manifold, manifold_moves, manifold_examples;

enum dim = 3;

void main()
{
    writeln("=== Incremental Valid Move Count Test ===");

    version (TrackValidMoves)
    {
        // Test dim=3 sphere with stellar subdivisions and random moves
        alias BM = BistellarMove!(dim, int);
        auto mfd = Manifold!dim([[0,1,2,3],[0,1,2,4],[0,1,3,4],[0,2,3,4],[1,2,3,4]]);
        auto rng = Mt19937(42);
        int nextV = 5;

        assert(mfd.validMoveCount == mfd.countValidBistellarMoves,
            "initial: incremental=%d full=%d".format(
                mfd.validMoveCount, mfd.countValidBistellarMoves));
        writefln("  initial: %d valid moves (verified)", mfd.validMoveCount);

        int nMoves = 0;
        foreach (_; 0 .. 100)
        {
            // Stellar subdivision to grow
            auto facet = mfd.randomFacetOfDim(dim).array;
            mfd.doMove(BM(facet, [nextV]));
            nextV++;
            nMoves++;

            auto inc = mfd.validMoveCount;
            auto full = mfd.countValidBistellarMoves;
            assert(inc == full,
                "after move %d (stellar): incremental=%d full=%d".format(nMoves, inc, full));

            // Random valid move
            auto moves = mfd.allBistellarMoves;
            if (!moves.empty)
            {
                auto idx = uniform(0, cast(int) moves.length, rng);
                mfd.doMove(moves[idx]);
                nMoves++;

                inc = mfd.validMoveCount;
                full = mfd.countValidBistellarMoves;
                assert(inc == full,
                    "after move %d (bistellar): incremental=%d full=%d".format(nMoves, inc, full));
            }

            if (_ % 10 == 0)
                writefln("  move %d: %d tets, %d valid moves (verified)",
                    nMoves, mfd.fVector[dim], mfd.validMoveCount);
        }

        writefln("  PASSED: %d moves, all incremental counts matched full recount", nMoves);
        writefln("  final: %d tets, %d valid moves", mfd.fVector[dim], mfd.validMoveCount);
    }
    else
    {
        writeln("  TrackValidMoves not enabled. Skipping.");
    }
}

}
