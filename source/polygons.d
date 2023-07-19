/// Contains functionality related to polygons and their triangulations
module polygons;

import std.algorithm, std.range;
import unit_threaded;

/******************************************************************************
Returns the set of symmetries of the a regular polygon with n sides.
*/
auto nGonSymmetries()(ulong n)
{
    return chain(nGonRotations(n), nGonReflections(n));
}
///
@Name("nGonSymmetries") pure @safe unittest
{
    3.nGonSymmetries.shouldBeSameSetAs([
        [0,1,2], [1,2,0], [2,0,1],  // rotations
        [2,1,0], [0,2,1], [1,0,2]   // reflections
    ]);

    4.nGonSymmetries.shouldBeSameSetAs([
        [0,1,2,3], [1,2,3,0], [2,3,0,1], [3,0,1,2], // rotations
        [3,2,1,0], [0,3,2,1], [1,0,3,2], [2,1,0,3]  // reflections
    ]);
}

/******************************************************************************
Returns the set of rotations of the a regular polygon with n sides.
*/
auto nGonRotations(ulong n)
{
    return n.iota.map!(k => n.iota.cycle(k).take(n).array);
}
///
@Name("nGonRotations") pure @safe unittest
{
    3.nGonRotations.shouldBeSameSetAs([
        [0,1,2], [1,2,0], [2,0,1]
    ]);
    4.nGonRotations.shouldBeSameSetAs([
        [0,1,2,3], [1,2,3,0], [2,3,0,1], [3,0,1,2]
    ]);    
}

/******************************************************************************
Returns the set of reflections of the n-gon.
*/
auto nGonReflections(ulong n)
{
    return n.iota.map!(k => n.iota.cycle(k).take(n).retro.array);
}
///
@Name("nGonReflections") pure @safe unittest
{
    3.nGonReflections.shouldBeSameSetAs([
        [2,1,0], [0,2,1], [1,0,2]
    ]);

    4.nGonReflections.shouldBeSameSetAs([
        [3,2,1,0], [0,3,2,1], [1,0,3,2], [2,1,0,3]
    ]);        
}

/******************************************************************************
Returns the set distinct triangulations of the labelled n-gon [0,1, ... ,n-1]
up to the action of the isometry group of the n-gon
*/
const(ulong)[][][] nGonTriangReps()(ulong n)
{
    if(n == 3)
    {
        return [[[0,1,2]]];
    }

    if(n == 4)
    {
        return [[[0,1,2],[0,2,3]]];
    }

    if(n == 5)
    {
        return [[[0,1,2],[0,2,3],[0,3,4]]];
    }

    if(n == 6)
    {
        return [
            [[0,1,2],[0,2,3],[0,3,4],[0,4,5]], 
            [[0,1,2],[0,2,3],[0,3,5],[3,4,5]],
            [[0,1,5],[1,2,3],[1,3,5],[3,4,5]]
        ];
    }

    if(n == 7)
    {
        return [
            [[0,1,2],[0,2,3],[0,3,4],[0,4,5],[0,5,6]], 
            [[0,1,3],[0,3,4],[0,4,5],[0,5,6],[1,2,3]],
            [[0,1,3],[0,3,4],[0,4,6],[1,2,3],[4,5,6]],
            [[0,1,2],[0,2,4],[0,4,5],[0,5,6],[2,3,4]]
        ];
    }

    // TO DO: Make function to generate these
    assert(0, "unsupported nubmer of sides");
}

@Name("nGonTriangReps") pure @safe unittest
{
    foreach(n; 3 .. 8)
    {
        foreach(triang; nGonTriangReps(n))
        {
            import simplicial_complex : simplicialComplex;
            import algorithms : isPureOfDim;
            import std.stdio : writeln;

            auto sc = simplicialComplex(triang);

            // Test that this is n-gon triangulation
            assert(sc.isPureOfDim(2));
            assert(sc.simplices(1).all!(
                edge => sc.star(edge).walkLength <= 2));

            auto boundaryEdges = sc.simplices(1).filter!(
                edge => sc.star(edge).walkLength == 1);
            
            auto ulongeriorEdges = sc.simplices(1).filter!(
                edge => sc.star(edge).walkLength == 2);

            boundaryEdges.walkLength.shouldEqual(n);
            ulongeriorEdges.walkLength.shouldEqual(n - 3);
            sc.simplices(2).walkLength.shouldEqual(n - 2);
        
            auto nGonEdges = n.iota.map!(
                k => (k < n - 1) ? [k, k + 1] : [0 , n - 1]);
            boundaryEdges.shouldBeSameSetAs(nGonEdges);

            sc.simplices(0).map!front.shouldBeSameSetAs(n.iota);
        }
    }
}

ulong[][][] nGonTriangs()(ulong n)
{
    ulong[][][] ans;
    foreach(rep; nGonTriangReps(n))
    {
        foreach(sym_; nGonSymmetries(n))
        {
            alias sym = f => f.map!(v => sym_[v]).array.sort.array;
            ans ~= rep.map!sym.array.sort.array;
        }
    }
    return ans.sort.uniq.array;
}
///
@Name("nGonTriangs") pure @safe unittest
{
    // import std.stdio : writeln;
    // 7.nGonTriangs.each!writeln;
    // 7.nGonTriangs.length.writeln;
    assert(3.nGonTriangs.walkLength == 1);
    assert(4.nGonTriangs.walkLength == 2);
    assert(5.nGonTriangs.walkLength == 5);
    assert(6.nGonTriangs.walkLength == 14);
    assert(7.nGonTriangs.walkLength == 42);

    // TO DO: More tests!
}


auto numNgonTriangs(size_t degree)
{
    static immutable result = [1,2,5,14,42];
    return result[degree-3];
}

unittest
{
    assert(numNgonTriangs(4) == 2);
    assert(numNgonTriangs(7) == 42);
}
