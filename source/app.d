// version(unittest) {}
// else
// {
//     void main()
//     {
//         // import sampler : sample;
//         // import gperftools_d.profiler;
//         // ProfilerStart();
//         // Test comment
//         // sample;
//         // ProfilerStop();
//         // import std.stdio : writeln;
//         // import std.algorithm : equal, map;
//         // import simplicial_complex : simplicialComplex;

//         // auto sc2 = simplicialComplex([[4,5], [5,6], [1,2,3], [2,3,4]]);
//         // sc2.toDetailedString.writeln;
//         // writeln("star([4]) = ", sc2.star([4]));
//         // writeln("link([5]) = ", sc2.link([5]));

//         // auto lnk = sc2.link([5]);
//         // auto lnk2 = lnk.save;
//         // assert(lnk.front.equal([4]));
//         // lnk.popFront;
//         // assert(lnk.front.equal([6]));
//         // lnk.popFront;
//         // assert(lnk.empty);


//         import std.algorithm : equal, map;
//         import fluent.asserts : should;
//         import std.range : iota, array;

//         auto ror = iota(1,4).map!iota;
//         assert(ror.equal!equal([[0],[0,1],[0,1,2]]));

//         // This slowly uses up all available memory:
// //        ror.should.equal(([[0],[0,1],[0,1,2]]));
//         ror.should.containOnly(([[0],[0,1],[0,1,2]]));

//         // However, this works:
//         ror.map!array.should.containOnly([[0],[0,1],[0,1,2]]);
//     }
// }


void main()
{
   import std.algorithm : equal, map;
   import fluent.asserts : should;
   import std.range : iota, array;

   auto ror = iota(1,4).map!iota;
   assert(ror.equal!equal([[0],[0,1],[0,1,2]]));

   // These slowly use up all available memory:
   ror.should.equal(([[0],[0,1],[0,1,2]]));
   // ror.should.containOnly(([[0],[0,1],[0,1,2]]));

   // However, these work:
   ror.map!array.should.containOnly([[0],[0,1],[0,1,2]]);
   ror.map!array.should.equal([[0],[0,1],[0,1,2]]);
}