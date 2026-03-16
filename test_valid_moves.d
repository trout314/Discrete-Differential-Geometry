/// Test for incremental valid move tracking.
/// Verifies that validMoveCount stays in sync with countValidBistellarMoves
/// after every move (stellar subdivisions and random bistellar moves).
module test_valid_moves;

version (unittest) {} else
{

import std.algorithm, std.array, std.conv, std.format, std.random, std.stdio;
import manifold, manifold_moves, manifold_examples;

void main()
{
    version (TrackValidMoves)
    {
        testDim2();
        testDim3();
        writeln("PASSED: all incremental valid move count tests");
    }
    else
    {
        writeln("TrackValidMoves not enabled — skipping");
    }
}

version (TrackValidMoves):

void testDim2()
{
    alias BM = BistellarMove!(2, int);
    // Start with octahedron
    auto mfd = Manifold!2([
        [0,1,2], [0,2,3], [0,3,4], [0,1,4],
        [1,2,5], [2,3,5], [3,4,5], [1,4,5],
    ]);

    assert(mfd.validMoveCount == mfd.countValidBistellarMoves);

    auto rng = Mt19937(42);
    int nextV = 6;

    foreach (_; 0 .. 50)
    {
        auto moves = mfd.allBistellarMoves;
        if (moves.empty) break;

        auto idx = uniform(0, cast(int) moves.length, rng);
        mfd.doMove(moves[idx]);

        if (uniform(0, 3, rng) == 0)
        {
            auto facet = mfd.facets.front.array;
            mfd.doMove(BM(facet, [nextV]));
            nextV++;
        }

        assert(mfd.validMoveCount == mfd.countValidBistellarMoves,
            "dim=2: incremental=%d full=%d".format(
                mfd.validMoveCount, mfd.countValidBistellarMoves));
    }
}

void testDim3()
{
    alias BM = BistellarMove!(3, int);
    auto mfd = Manifold!3([
        [0,1,2,3], [0,1,2,4], [0,1,3,4], [0,2,3,4], [1,2,3,4],
    ]);

    assert(mfd.validMoveCount == mfd.countValidBistellarMoves);

    auto rng = Mt19937(123);
    int nextV = 5;

    foreach (_; 0 .. 100)
    {
        // Stellar subdivision to grow
        auto facet = mfd.facets.front.array;
        mfd.doMove(BM(facet, [nextV]));
        nextV++;

        assert(mfd.validMoveCount == mfd.countValidBistellarMoves,
            "dim=3: incremental=%d full=%d after stellar subdiv".format(
                mfd.validMoveCount, mfd.countValidBistellarMoves));

        // Random valid move
        auto moves = mfd.allBistellarMoves;
        if (!moves.empty)
        {
            auto idx = uniform(0, cast(int) moves.length, rng);
            mfd.doMove(moves[idx]);

            assert(mfd.validMoveCount == mfd.countValidBistellarMoves,
                "dim=3: incremental=%d full=%d after bistellar move".format(
                    mfd.validMoveCount, mfd.countValidBistellarMoves));
        }
    }
}

}
