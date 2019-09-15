module moves;

import std.conv : to;
import utility : isInputRangeOf, isInputRangeOfInputRangeOf;
import std.range : walkLength;
import std.algorithm : copy;
import unit_threaded : Name, shouldBeSameSetAs, shouldEqual, shouldBeTrue, shouldBeFalse;

alias isIRof = isInputRangeOf;
alias isIRofIRof = isInputRangeOfInputRangeOf;

struct Move(int dim, Vertex = int, int maxHingeDeg=7)
{
public:   
    alias center = getCenter;
    alias coCenter = getCoCenter;

    bool isHinge()() const
    {
        return lenCenter + lenCoCenter > dim + 2;
    }

    bool isPachner()() const
    {
        return lenCenter + lenCoCenter == dim + 2;
    }

    int triangIndx()() const
    {
        assert(isHinge, "no triangulation index for a pachner move");
        return triangIndx_;
    }

    real weight()() const
    {
        if(isHinge)
        {
            return 1.0 / lenCoCenter;
        }
        else
        {
            return 1.0 / (dim + 2 - lenCenter);
        }
    }

    string toString()() const
    {
        string move_name;
        if (isPachner)
        {
            return "pachner move, center=" ~ center.to!string 
                ~ " coCenter=" ~ coCenter.to!string;
        }
        else
        {
            return "hinge move, center=" ~ center.to!string 
                ~ " coCenter=" ~ coCenter.to!string
                ~ " triangIndx=" ~ triangIndx.to!string;
        }
    }

    this(R1, R2)(R1 center_, R2 coCenter_, int hingeMoveTriangIndx = -1)
    if (isIRof!(R1, Vertex) && isIRof!(R2, Vertex))
    {
        lenCenter = center_.walkLength;
        lenCoCenter = coCenter_.walkLength;
        triangIndx_ = hingeMoveTriangIndx;

        assert(lenCenter + lenCoCenter <= vertices.length, "too many vertices for a move");
        assert(lenCenter + lenCoCenter >= dim + 2, "not enough vertices for a move");

        if(lenCenter + lenCoCenter > dim + 2)
        {
            assert(lenCenter == dim - 1, "center is not a hinge");
            assert(lenCoCenter >= 4);
            assert(triangIndx_ >= 0, "must set a valid trianulation index");
        }

        copy(center_, center);
        copy(coCenter_, coCenter);
    }
private:
    Vertex[dim + maxHingeDeg - 1] vertices;
    size_t lenCenter;
    size_t lenCoCenter;

    // indicates which triangulation of a disk to use
    // (where the boundary of this disk is the coCenter)
    int triangIndx_;
    
    inout(Vertex)[] getCenter()() inout
    {
        return vertices[0..lenCenter];
    }

    inout(Vertex)[] getCoCenter()() inout
    {
        return vertices[lenCenter..lenCenter + lenCoCenter];
    }
}

///
@Name("Moves doc tests") pure @safe unittest
{
    auto hm = Move!3([1,2], [4, 5, 6, 7], 42);
    hm.center.shouldBeSameSetAs([1,2]);
    hm.coCenter.shouldBeSameSetAs([4,5,6,7]);
    hm.triangIndx.shouldEqual(42);
    hm.isHinge.shouldBeTrue;
    hm.isPachner.shouldBeFalse;

    auto pm = Move!4([1,2,3], [4,5,6]);
    pm.center.shouldBeSameSetAs([1,2,3]);
    pm.coCenter.shouldBeSameSetAs([4,5,6]);
    pm.isHinge.shouldBeFalse;
    pm.isPachner.shouldBeTrue;

}

    
        
