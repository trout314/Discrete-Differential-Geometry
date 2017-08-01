import std.algorithm : any, filter, find, map;
import std.range : array, enumerate, front, walkLength;
import manifold_small : SmallManifold;
import unit_threaded : Name;

import utility : subsetsOfSize;

/*******************************************************************************
Returns true if the given manifold is orientable and false otherwise.
*/
bool isOrientable(Vertex, int dim)(SmallManifold!(dim, Vertex) manifold)
{
    /* We must choose a compatible orientation for each facet. Since the facets
    already come equipped with an ordering for the vertices, we need only
    indicate whether this orientation is the one we want (GivenOrder) or not
    (OppositeOrder). We let NotSet indicate that the facet has not yet been
    processed. */
    enum Orientation
    {
        NotSet,
        GivenOrder,
        OppositeOrder
    }

    static struct FacetRecord
    {
        Vertex[] facet;
        Orientation label;
        bool done;
    }

    // An empty manifold is orientable
    if(manifold.numFacets == 0)
    {
        return true;
    }

    auto records = manifold.facets.map!(f => FacetRecord(f)).array;

    // We (arbitrarily) use the given vertex order to orient the first simplex
    assert(records.walkLength > 0);
    records.front.label = Orientation.GivenOrder;

    while(records.any!(r => !r.done))
    {
        /* If we're not done there must be some oriented facet that hasn't had
        its neighbors orientations set (or checked OK if already set) */
        auto toDo = records.find!(r => !r.done && r.label != Orientation.NotSet);
        assert(toDo.walkLength > 0);

        foreach(i, ridge; toDo.front.facet.subsetsOfSize(dim).enumerate)
        {
            auto oppFacet = manifold.star(ridge).filter!(f => f != toDo.front.facet);
            assert(!oppFacet.empty);

            auto j = oppFacet.front.subsetsOfSize(dim).enumerate
                .find!(p => p.value == ridge).front.index;

            Orientation oppFacetLabel;
            if(i % 2 != j % 2)
            {
                oppFacetLabel = toDo.front.label;
            }
            else
            {
                if(toDo.front.label == Orientation.GivenOrder)
                {
                    oppFacetLabel = Orientation.OppositeOrder;
                }
                else
                {
                    oppFacetLabel = Orientation.GivenOrder;
                }
            }

            auto oppRecord = records.find!(r => r.facet == oppFacet.front);
            if(oppRecord.front.label != Orientation.NotSet)
            {
                if (oppRecord.front.label != oppFacetLabel)
                {
                    return false;
                }
            }
            else
            {
                oppRecord.front.label = oppFacetLabel;
            }           
        }
        toDo.front.done = true;
    }
    return true;
}
///
@Name("isOrientable") unittest
{
    assert(SmallManifold!2().isOrientable);

    // http://page.math.tu-berlin.de/~lutz/stellar/manifolds_lex/manifolds_lex_d2_n10_o0_g5
    // Surface #4941 on non-orientable genus 5 list
    auto g5 = SmallManifold!2([[1, 2, 3], [1, 2, 4], [1, 3, 5], [1, 4, 6], [1,
            5, 7], [1, 6, 7], [2, 3, 6], [2, 4, 8], [2, 5, 7], [2, 5, 9], [2,
            6, 10], [2, 7, 8], [2, 9, 10], [3, 4, 6], [3, 4, 10], [3, 5, 8],
            [3, 7, 8], [3, 7, 10], [4, 5, 9], [4, 5, 10], [4, 8, 9], [5, 8,
            10], [6, 7, 9], [6, 8, 9], [6, 8, 10], [7, 9, 10]]);
    assert(!g5.isOrientable);

    // http://page.math.tu-berlin.de/~lutz/stellar/manifolds_lex/manifolds_lex_d2_n9_o1_g0
    // Surface #15 on orientable genus 0 list
    auto g0 = SmallManifold!2([[1, 2, 3], [1, 2, 4], [1, 3, 4], [2, 3, 5], [2,
            4, 5], [3, 4, 6], [3, 5, 7], [3, 6, 7], [4, 5, 8], [4, 6, 7], [4,
            7, 9], [4, 8, 9], [5, 7, 8], [7, 8, 9]]);
    assert(g0.isOrientable);

    // http://page.math.tu-berlin.de/~lutz/stellar/library_of_triangulations/poincare
    auto poincare = SmallManifold!3([[1, 2, 4, 9], [1, 2, 4, 15], [1, 2, 6,
            14], [1, 2, 6, 15], [1, 2, 9, 14], [1, 3, 4, 12], [1, 3, 4, 15],
            [1, 3, 7, 10], [1, 3, 7, 12], [1, 3, 10, 15], [1, 4, 9, 12], [1, 5,
            6, 13], [1, 5, 6, 14], [1, 5, 8, 11], [1, 5, 8, 13], [1, 5, 11,
            14], [1, 6, 13, 15], [1, 7, 8, 10], [1, 7, 8, 11], [1, 7, 11, 12],
            [1, 8, 10, 13], [1, 9, 11, 12], [1, 9, 11, 14], [1, 10, 13, 15],
            [2, 3, 5, 10], [2, 3, 5, 11], [2, 3, 7, 10], [2, 3, 7, 13], [2, 3,
            11, 13], [2, 4, 9, 13], [2, 4, 11, 13], [2, 4, 11, 15], [2, 5, 8,
            11], [2, 5, 8, 12], [2, 5, 10, 12], [2, 6, 10, 12], [2, 6, 10, 14],
            [2, 6, 12, 15], [2, 7, 9, 13], [2, 7, 9, 14], [2, 7, 10, 14], [2,
            8, 11, 15], [2, 8, 12, 15], [3, 4, 5, 14], [3, 4, 5, 15], [3, 4,
            12, 14], [3, 5, 10, 15], [3, 5, 11, 14], [3, 7, 12, 13], [3, 11,
            13, 14], [3, 12, 13, 14], [4, 5, 6, 7], [4, 5, 6, 14], [4, 5, 7,
            15], [4, 6, 7, 11], [4, 6, 10, 11], [4, 6, 10, 14], [4, 7, 11, 15],
            [4, 8, 9, 12], [4, 8, 9, 13], [4, 8, 10, 13], [4, 8, 10, 14], [4,
            8, 12, 14], [4, 10, 11, 13], [5, 6, 7, 13], [5, 7, 9, 13], [5, 7,
            9, 15], [5, 8, 9, 12], [5, 8, 9, 13], [5, 9, 10, 12], [5, 9, 10,
            15], [6, 7, 11, 12], [6, 7, 12, 13], [6, 10, 11, 12], [6, 12, 13,
            15], [7, 8, 10, 14], [7, 8, 11, 15], [7, 8, 14, 15], [7, 9, 14,
            15], [8, 12, 14, 15], [9, 10, 11, 12], [9, 10, 11, 16], [9, 10, 15,
            16], [9, 11, 14, 16], [9, 14, 15, 16], [10, 11, 13, 16], [10, 13,
            15, 16], [11, 13, 14, 16], [12, 13, 14, 15], [13, 14, 15, 16]]);
    assert(poincare.isOrientable);

    /* This should be an 5-sphere, hence orientable. See:
    http://page.math.tu-berlin.de/~lutz/stellar/5_manifolds
    http://page.math.tu-berlin.de/~lutz/stellar/5_manifolds.type */
    auto manifold_5_10_3_1 = SmallManifold!5([[1, 2, 3, 4, 5, 6], [1, 2, 3, 4,
            5, 10], [1, 2, 3, 4, 6, 7], [1, 2, 3, 4, 7, 8], [1, 2, 3, 4, 8, 9],
            [1, 2, 3, 4, 9, 10], [1, 2, 3, 5, 6, 10], [1, 2, 3, 6, 7, 10], [1,
            2, 3, 7, 8, 10], [1, 2, 3, 8, 9, 10], [1, 2, 4, 5, 6, 7], [1, 2, 4,
            5, 7, 8], [1, 2, 4, 5, 8, 9], [1, 2, 4, 5, 9, 10], [1, 2, 5, 6, 7,
            8], [1, 2, 5, 6, 8, 9], [1, 2, 5, 6, 9, 10], [1, 2, 6, 7, 8, 9],
            [1, 2, 6, 7, 9, 10], [1, 2, 7, 8, 9, 10], [1, 3, 4, 5, 6, 10], [1,
            3, 4, 6, 7, 10], [1, 3, 4, 7, 8, 10], [1, 3, 4, 8, 9, 10], [1, 4,
            5, 6, 7, 10], [1, 4, 5, 7, 8, 10], [1, 4, 5, 8, 9, 10], [1, 5, 6,
            7, 8, 10], [1, 5, 6, 8, 9, 10], [1, 6, 7, 8, 9, 10], [2, 3, 4, 5,
            6, 7], [2, 3, 4, 5, 7, 8], [2, 3, 4, 5, 8, 9], [2, 3, 4, 5, 9, 10],
            [2, 3, 5, 6, 7, 8], [2, 3, 5, 6, 8, 9], [2, 3, 5, 6, 9, 10], [2, 3,
            6, 7, 8, 9], [2, 3, 6, 7, 9, 10], [2, 3, 7, 8, 9, 10], [3, 4, 5, 6,
            7, 8], [3, 4, 5, 6, 8, 9], [3, 4, 5, 6, 9, 10], [3, 4, 6, 7, 8, 9],
            [3, 4, 6, 7, 9, 10], [3, 4, 7, 8, 9, 10], [4, 5, 6, 7, 8, 9], [4,
            5, 6, 7, 9, 10], [4, 5, 7, 8, 9, 10], [5, 6, 7, 8, 9, 10]]);
    assert(manifold_5_10_3_1.isOrientable);

    // http://page.math.tu-berlin.de/~lutz/stellar/manifolds_lex/manifolds_lex_d2_n9_o0_g1
    // surface #9 on the genus 1 non-orientable list
    auto g1 = SmallManifold!2([[1, 2, 3], [1, 2, 4], [1, 3, 4], [2, 3, 5], [2,
            4, 5], [3, 4, 6], [3, 5, 7], [3, 6, 8], [3, 7, 8], [4, 5, 8], [4,
            6, 9], [4, 7, 8], [4, 7, 9], [5, 6, 8], [5, 6, 9], [5, 7, 9]]);
    assert(!g1.isOrientable);

    // http://page.math.tu-berlin.de/~lutz/stellar/RP3
    auto rp3 = SmallManifold!3([[1, 2, 3, 7], [1, 2, 3, 11], [1, 2, 6, 9], [1,
            2, 6, 11], [1, 2, 7, 9], [1, 3, 5, 10], [1, 3, 5, 11], [1, 3, 7,
            10], [1, 4, 7, 9], [1, 4, 7, 10], [1, 4, 8, 9], [1, 4, 8, 10], [1,
            5, 6, 8], [1, 5, 6, 11], [1, 5, 8, 10], [1, 6, 8, 9], [2, 3, 4, 8],
            [2, 3, 4, 11], [2, 3, 7, 8], [2, 4, 6, 10], [2, 4, 6, 11], [2, 4,
            8, 10], [2, 5, 7, 8], [2, 5, 7, 9], [2, 5, 8, 10], [2, 5, 9, 10],
            [2, 6, 9, 10], [3, 4, 5, 9], [3, 4, 5, 11], [3, 4, 8, 9], [3, 5, 9,
            10], [3, 6, 7, 8], [3, 6, 7, 10], [3, 6, 8, 9], [3, 6, 9, 10], [4,
            5, 6, 7], [4, 5, 6, 11], [4, 5, 7, 9], [4, 6, 7, 10], [5, 6, 7, 8]]);
    assert(rp3.isOrientable);

    auto rp4 = SmallManifold!4([[1, 2, 4, 5, 11], [1, 2, 4, 5, 14], [1, 2, 4, 11,
        13], [1, 2, 4, 13, 14], [1, 2, 5, 11, 15], [1, 2, 5, 14, 15], [1, 2,
        7, 8, 13], [1, 2, 7, 8, 15], [1, 2, 7, 13, 14], [1, 2, 7, 14, 15], [1,
        2, 8, 11, 13], [1, 2, 8, 11, 15], [1, 3, 4, 10, 13], [1, 3, 4, 10, 16],
        [1, 3, 4, 11, 13], [1, 3, 4, 11, 16], [1, 3, 8, 9, 11], [1, 3, 8, 9,
        12], [1, 3, 8, 11, 13], [1, 3, 8, 12, 13], [1, 3, 9, 11, 16], [1, 3, 9,
        12, 16], [1, 3, 10, 12, 13], [1, 3, 10, 12, 16], [1, 4, 5, 11, 16], [1,
        4, 5, 14, 16], [1, 4, 10, 13, 14], [1, 4, 10, 14, 16], [1, 5, 6, 9,
        15], [1, 5, 6, 9, 16], [1, 5, 6, 14, 15], [1, 5, 6, 14, 16], [1, 5,
        9, 11, 15], [1, 5, 9, 11, 16], [1, 6, 7, 10, 12], [1, 6, 7, 10, 14],
        [1, 6, 7, 12, 15], [1, 6, 7, 14, 15], [1, 6, 9, 12, 15], [1, 6, 9, 12,
        16], [1, 6, 10, 12, 16], [1, 6, 10, 14, 16], [1, 7, 8, 12, 13], [1, 7,
        8, 12, 15], [1, 7, 10, 12, 13], [1, 7, 10, 13, 14], [1, 8, 9, 11, 15],
        [1, 8, 9, 12, 15], [2, 3, 5, 10, 12], [2, 3, 5, 10, 15], [2, 3, 5, 12,
        14], [2, 3, 5, 14, 15], [2, 3, 7, 9, 14], [2, 3, 7, 9, 16], [2, 3, 7,
        14, 15], [2, 3, 7, 15, 16], [2, 3, 9, 12, 14], [2, 3, 9, 12, 16], [2,
        3, 10, 12, 16], [2, 3, 10, 15, 16], [2, 4, 5, 11, 12], [2, 4, 5, 12,
        14], [2, 4, 6, 9, 12], [2, 4, 6, 9, 13], [2, 4, 6, 11, 12], [2, 4, 6,
        11, 13], [2, 4, 9, 12, 14], [2, 4, 9, 13, 14], [2, 5, 10, 11, 12], [2,
        5, 10, 11, 15], [2, 6, 8, 10, 11], [2, 6, 8, 10, 16], [2, 6, 8, 11,
        13], [2, 6, 8, 13, 16], [2, 6, 9, 12, 16], [2, 6, 9, 13, 16], [2, 6,
        10, 11, 12], [2, 6, 10, 12, 16], [2, 7, 8, 13, 16], [2, 7, 8, 15, 16],
        [2, 7, 9, 13, 14], [2, 7, 9, 13, 16], [2, 8, 10, 11, 15], [2, 8, 10,
        15, 16], [3, 4, 6, 7, 11], [3, 4, 6, 7, 15], [3, 4, 6, 11, 13], [3, 4,
        6, 13, 15], [3, 4, 7, 11, 16], [3, 4, 7, 15, 16], [3, 4, 10, 13, 15],
        [3, 4, 10, 15, 16], [3, 5, 6, 8, 13], [3, 5, 6, 8, 14], [3, 5, 6, 13,
        15], [3, 5, 6, 14, 15], [3, 5, 8, 12, 13], [3, 5, 8, 12, 14], [3, 5,
        10, 12, 13], [3, 5, 10, 13, 15], [3, 6, 7, 11, 14], [3, 6, 7, 14, 15],
        [3, 6, 8, 11, 13], [3, 6, 8, 11, 14], [3, 7, 9, 11, 14], [3, 7, 9, 11,
        16], [3, 8, 9, 11, 14], [3, 8, 9, 12, 14], [4, 5, 7, 8, 12], [4, 5, 7,
        8, 16], [4, 5, 7, 11, 12], [4, 5, 7, 11, 16], [4, 5, 8, 12, 14], [4, 5,
        8, 14, 16], [4, 6, 7, 11, 12], [4, 6, 7, 12, 15], [4, 6, 9, 12, 15],
        [4, 6, 9, 13, 15], [4, 7, 8, 12, 15], [4, 7, 8, 15, 16], [4, 8, 9, 10,
        14], [4, 8, 9, 10, 15], [4, 8, 9, 12, 14], [4, 8, 9, 12, 15], [4, 8,
        10, 14, 16], [4, 8, 10, 15, 16], [4, 9, 10, 13, 14], [4, 9, 10, 13,
        15], [5, 6, 8, 13, 16], [5, 6, 8, 14, 16], [5, 6, 9, 13, 15], [5, 6,
        9, 13, 16], [5, 7, 8, 12, 13], [5, 7, 8, 13, 16], [5, 7, 9, 10, 11],
        [5, 7, 9, 10, 13], [5, 7, 9, 11, 16], [5, 7, 9, 13, 16], [5, 7, 10, 11,
        12], [5, 7, 10, 12, 13], [5, 9, 10, 11, 15], [5, 9, 10, 13, 15], [6, 7,
        10, 11, 12], [6, 7, 10, 11, 14], [6, 8, 10, 11, 14], [6, 8, 10, 14,
        16], [7, 9, 10, 11, 14], [7, 9, 10, 13, 14], [8, 9, 10, 11, 14], [8, 9, 
        10, 11, 15]]);
    assert(!rp4.isOrientable);
}
