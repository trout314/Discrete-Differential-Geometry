/**
 * Open-addressing hash map with linear probing.
 *
 * Drop-in replacement for D's built-in associative arrays (V[K]).
 * Uses flat arrays for cache-friendly access — no pointer chasing,
 * no GC-allocated nodes per entry.
 *
 * Ported from the pachner_data_structure project.
 */
module hashmap;

import std.traits : isIntegral;
import utility : StackArray;

/**
 * Open-addressing hash map with inline key/value storage.
 *
 * Invariants:
 *   - Load factor kept below 0.75 (resize at 3/4 capacity)
 *   - Capacity is always a power of 2 (mask-based indexing)
 *   - Tombstones used for deletion to preserve probe chains
 */
struct HashMap(K, V)
{
    private enum ubyte EMPTY = 0;
    private enum ubyte OCCUPIED = 1;

    private K[] _keys;
    private V[] _values;
    private ubyte[] _ctrl;
    private size_t _length;
    private size_t _mask; // capacity - 1

    /// Return the number of entries.
    size_t length() const pure nothrow @nogc @safe { return _length; }

    /// Lookup: returns pointer to value if found, null otherwise.
    inout(V)* opBinaryRight(string op : "in")(K key) inout
    {
        if (_ctrl.length == 0) return null;
        auto i = locate(key);
        if (i == size_t.max) return null;
        return &(cast(inout(V)[])_values)[i];
    }

    /// Read by key (asserts key exists).
    ref inout(V) opIndex(K key) inout
    {
        auto p = key in this;
        assert(p !is null, "Key not found in HashMap");
        return *p;
    }

    /// Insert or update.
    void opIndexAssign(V value, K key)
    {
        ensureCapacity();
        auto i = locateOrInsert(key);
        _values[i] = value;
    }

    /// Prefix increment/decrement (e.g. ++map[key], --map[key]).
    /// Auto-inserts with default value if key is absent.
    ref V opIndexUnary(string op)(K key) if (op == "++" || op == "--")
    {
        ensureCapacity();
        auto i = locateOrInsert(key);
        mixin("return " ~ op ~ "_values[i];");
    }

    /// Compound assignment (e.g. map[key] += 5).
    /// Auto-inserts with default value if key is absent.
    ref V opIndexOpAssign(string op)(V rhs, K key)
    {
        ensureCapacity();
        auto i = locateOrInsert(key);
        mixin("return _values[i] " ~ op ~ "= rhs;");
    }

    /// Remove a key using backward-shift deletion (Knuth's Algorithm R).
    /// Avoids tombstones entirely, so the table never needs rehashing
    /// just to clear accumulated tombstones.
    bool remove(K key)
    {
        if (_ctrl.length == 0) return false;
        auto i = locate(key);
        if (i == size_t.max) return false;

        // Backward-shift: scan forward from the deleted slot.  For each
        // occupied slot j, compute its ideal slot k.  If the entry at j
        // was displaced past the gap at i (i.e. its probe distance from k
        // is >= the distance from i to j), shift it back into the gap.
        auto j = i;
        while (true)
        {
            j = (j + 1) & _mask;
            if (_ctrl[j] != OCCUPIED)
                break;

            auto k = computeHash(_keys[j]) & _mask;
            // Move if entry at j was probed past the gap at i.
            // Works correctly with wrapping for power-of-2 tables.
            if (((j - k) & _mask) >= ((j - i) & _mask))
            {
                _keys[i] = _keys[j];
                _values[i] = _values[j];
                _ctrl[i] = OCCUPIED;
                i = j;
            }
        }
        _ctrl[i] = EMPTY;
        _keys[i] = K.init;
        _values[i] = V.init;
        _length--;
        return true;
    }

    /// Clear all entries and deallocate.
    void clear()
    {
        _keys = null;
        _values = null;
        _ctrl = null;
        _length = 0;
        _mask = 0;
    }

    /// Return all keys as an array.
    K[] keys() const
    {
        auto result = new K[](_length);
        size_t j = 0;
        foreach (i; 0 .. _ctrl.length)
            if (_ctrl[i] == OCCUPIED)
                result[j++] = _keys[i];
        return result;
    }

    /// Return a shallow copy.
    HashMap dup() const
    {
        HashMap copy;
        copy._keys = _keys.dup;
        copy._values = _values.dup;
        copy._ctrl = _ctrl.dup;
        copy._length = _length;
        copy._mask = _mask;
        return copy;
    }

    /// Iterate over key-value pairs.
    auto byKeyValue() const
    {
        return KVRange!(const(K), const(V))(_keys, _values, _ctrl);
    }

    /// ditto
    auto byKeyValue()
    {
        return KVRange!(K, V)(_keys, _values, _ctrl);
    }

    /// Iterate over keys.
    auto byKey() const
    {
        return KeyRange!(const(K))(_keys, _ctrl);
    }

    /// ditto
    auto byKey()
    {
        return KeyRange!(K)(_keys, _ctrl);
    }

    // -----------------------------------------------------------------------
    // Internals
    // -----------------------------------------------------------------------

    /// Hash function. Fibonacci hashing for integers, FNV-1a for arrays/structs.
    private static size_t computeHash(K key) pure nothrow @nogc @trusted
    {
        static if (isIntegral!K)
        {
            // Fibonacci hashing
            return cast(size_t)(cast(ulong) key * 11_400_714_819_323_198_485UL);
        }
        else static if (__traits(isStaticArray, K))
        {
            return fnvHash(key);
        }
        else static if (isStackArray!K)
        {
            // Hash only the used portion of the StackArray
            return fnvHash(key[]);
        }
        else
        {
            return .hashOf(key);
        }
    }

    /// FNV-1a hash over elements of a range/array.
    private static size_t fnvHash(R)(R data) pure nothrow @nogc @trusted
    {
        size_t h = 14_695_981_039_346_656_037UL;
        foreach (e; data)
        {
            auto bytes = (cast(const(ubyte)*)&e)[0 .. typeof(e).sizeof];
            foreach (b; bytes)
            {
                h ^= b;
                h *= 1_099_511_628_211UL;
            }
        }
        return h;
    }

    /// Find the index of an existing key, or size_t.max if not found.
    private size_t locate(K key) const pure nothrow @nogc @safe
    {
        if (_ctrl.length == 0) return size_t.max;
        auto h = computeHash(key) & _mask;
        while (true)
        {
            if (_ctrl[h] == EMPTY) return size_t.max;
            if (_ctrl[h] == OCCUPIED && _keys[h] == key) return h;
            h = (h + 1) & _mask;
        }
    }

    /// Find slot for key (existing or first EMPTY for insert).
    private size_t locateOrInsert(K key)
    {
        auto h = computeHash(key) & _mask;
        while (true)
        {
            if (_ctrl[h] == EMPTY)
            {
                _ctrl[h] = OCCUPIED;
                _keys[h] = key;
                _length++;
                return h;
            }
            if (_keys[h] == key)
            {
                return h;
            }
            h = (h + 1) & _mask;
        }
    }

    /// Grow the table if load factor exceeds 75%.
    /// With backward-shift deletion (no tombstones), this only triggers
    /// when the actual number of entries exceeds the threshold.
    private void ensureCapacity()
    {
        if (_ctrl.length == 0)
        {
            allocate(16);
            return;
        }
        if (_length * 4 >= _ctrl.length * 3)
            grow();
    }

    private void allocate(size_t cap)
    {
        _keys = new K[cap];
        _values = new V[cap];
        _ctrl = new ubyte[cap];
        _ctrl[] = EMPTY;
        _mask = cap - 1;
    }

    private void grow()
    {
        auto oldKeys = _keys;
        auto oldValues = _values;
        auto oldCtrl = _ctrl;

        auto newCap = _ctrl.length * 2;
        allocate(newCap);
        _length = 0;

        foreach (i; 0 .. oldCtrl.length)
        {
            if (oldCtrl[i] == OCCUPIED)
            {
                auto slot = locateOrInsert(oldKeys[i]);
                _values[slot] = oldValues[i];
            }
        }
    }

    /// Key-value range
    private static struct KVRange(KT, VT)
    {
        KT[] _keys;
        VT[] _values;
        const(ubyte)[] _ctrl;
        size_t _pos;

        this(KT[] k, VT[] v, const(ubyte)[] c)
        {
            _keys = k;
            _values = v;
            _ctrl = c;
            _pos = 0;
            advance();
        }

        private void advance()
        {
            while (_pos < _ctrl.length && _ctrl[_pos] != OCCUPIED)
                _pos++;
        }

        bool empty() const pure nothrow @nogc @safe
        {
            return _pos >= _ctrl.length;
        }

        auto front()
        {
            struct KV { KT key; VT value; }
            return KV(_keys[_pos], _values[_pos]);
        }

        void popFront()
        {
            _pos++;
            advance();
        }
    }

    /// Key-only range
    private static struct KeyRange(KT)
    {
        KT[] _keys;
        const(ubyte)[] _ctrl;
        size_t _pos;

        this(KT[] k, const(ubyte)[] c)
        {
            _keys = k;
            _ctrl = c;
            _pos = 0;
            advance();
        }

        private void advance()
        {
            while (_pos < _ctrl.length && _ctrl[_pos] != OCCUPIED)
                _pos++;
        }

        bool empty() const pure nothrow @nogc @safe
        {
            return _pos >= _ctrl.length;
        }

        KT front() const
        {
            return _keys[_pos];
        }

        void popFront()
        {
            _pos++;
            advance();
        }
    }
}

/// Detect StackArray types
private template isStackArray(T)
{
    static if (is(T == StackArray!(U, n), U, size_t n))
        enum isStackArray = true;
    else
        enum isStackArray = false;
}

// ---------------------------------------------------------------------------
// Unit tests
// ---------------------------------------------------------------------------

unittest
{
    // Basic int -> int
    HashMap!(int, int) map;
    map[42] = 100;
    assert(map.length == 1);
    assert(map[42] == 100);
    assert(42 in map);
    assert(99 !in map);

    map[42] = 200;
    assert(map[42] == 200);
    assert(map.length == 1);

    map.remove(42);
    assert(42 !in map);
    assert(map.length == 0);
}

unittest
{
    // StackArray keys (matches degreeMap usage)
    alias SA = StackArray!(int, 4);
    HashMap!(SA, size_t) map;

    SA key1;
    key1 ~= 1; key1 ~= 2; key1 ~= 3;

    SA key2;
    key2 ~= 1; key2 ~= 2; key2 ~= 4;

    map[key1] = 10;
    map[key2] = 20;
    assert(map.length == 2);
    assert(map[key1] == 10);
    assert(map[key2] == 20);

    // Increment operator
    map[key1] += 5;
    assert(map[key1] == 15);

    // dup
    auto copy = map.dup;
    assert(copy.length == 2);
    assert(copy[key1] == 15);
    copy[key1] = 999;
    assert(map[key1] == 15); // original unchanged
}

unittest
{
    // Static array keys (matches ridgeLinks usage)
    HashMap!(int[3], int) map;
    int[3] k1 = [1, 2, 3];
    int[3] k2 = [4, 5, 6];

    map[k1] = 100;
    map[k2] = 200;
    assert(map.length == 2);
    assert(map[k1] == 100);

    map.remove(k1);
    assert(k1 !in map);
    assert(map.length == 1);
}

unittest
{
    // Stress test: many insertions trigger resizing
    HashMap!(int, int) map;
    foreach (i; 0 .. 1000)
        map[i] = i * 2;

    assert(map.length == 1000);
    foreach (i; 0 .. 1000)
        assert(map[i] == i * 2);

    // Remove half
    foreach (i; 0 .. 500)
        map.remove(i);
    assert(map.length == 500);

    // Remaining entries still accessible
    foreach (i; 500 .. 1000)
        assert(map[i] == i * 2);
}

unittest
{
    // byKeyValue iteration
    HashMap!(int, int) map;
    map[1] = 10;
    map[2] = 20;
    map[3] = 30;

    int sum = 0;
    foreach (kv; map.byKeyValue())
        sum += kv.value;
    assert(sum == 60);
}

unittest
{
    // Increment on absent key (default-initializes to 0)
    HashMap!(int, size_t) map;
    ++map[5];
    assert(map[5] == 1);
    ++map[5];
    assert(map[5] == 2);
}

unittest
{
    // Decrement operator on absent and present keys
    HashMap!(int, int) map;
    --map[1];
    assert(map[1] == -1);
    --map[1];
    assert(map[1] == -2);
    assert(map.length == 1);
}

unittest
{
    // Compound assignment operators: -=, *=
    HashMap!(int, int) map;

    // -= on absent key (starts from V.init == 0)
    map[1] -= 3;
    assert(map[1] == -3);

    // *= on present key
    map[2] = 5;
    map[2] *= 4;
    assert(map[2] == 20);

    // += on present key
    map[2] += 10;
    assert(map[2] == 30);
}

unittest
{
    // 'in' operator returns usable pointer to value
    HashMap!(int, int) map;
    map[7] = 42;

    auto p = 7 in map;
    assert(p !is null);
    assert(*p == 42);

    // Mutate through pointer
    *p = 99;
    assert(map[7] == 99);

    // Absent key returns null
    auto q = 8 in map;
    assert(q is null);
}

unittest
{
    // 'in' on empty map returns null
    HashMap!(int, int) map;
    assert(42 !in map);
}

unittest
{
    // remove on empty map returns false
    HashMap!(int, int) map;
    assert(!map.remove(1));
}

unittest
{
    // remove on non-existent key returns false
    HashMap!(int, int) map;
    map[1] = 10;
    assert(!map.remove(2));
    assert(map.length == 1);
}

unittest
{
    // clear deallocates and resets
    HashMap!(int, int) map;
    map[1] = 10;
    map[2] = 20;
    assert(map.length == 2);

    map.clear();
    assert(map.length == 0);
    assert(1 !in map);
    assert(2 !in map);

    // Can reuse after clear
    map[3] = 30;
    assert(map.length == 1);
    assert(map[3] == 30);
}

unittest
{
    // keys() returns all inserted keys
    HashMap!(int, int) map;
    map[5] = 50;
    map[3] = 30;
    map[8] = 80;

    auto k = map.keys();
    assert(k.length == 3);

    // Sort for deterministic comparison
    import std.algorithm : sort;
    k.sort();
    assert(k == [3, 5, 8]);
}

unittest
{
    // keys() on empty map returns empty array
    HashMap!(int, int) map;
    assert(map.keys().length == 0);
}

unittest
{
    // keys() excludes removed entries
    HashMap!(int, int) map;
    map[1] = 10;
    map[2] = 20;
    map[3] = 30;
    map.remove(2);

    auto k = map.keys();
    assert(k.length == 2);

    import std.algorithm : sort;
    k.sort();
    assert(k == [1, 3]);
}

unittest
{
    // byKey iteration
    HashMap!(int, int) map;
    map[10] = 100;
    map[20] = 200;
    map[30] = 300;

    int sum = 0;
    int count = 0;
    foreach (k; map.byKey())
    {
        sum += k;
        count++;
    }
    assert(count == 3);
    assert(sum == 60);
}

unittest
{
    // byKey on empty map yields nothing
    HashMap!(int, int) map;
    foreach (k; map.byKey())
        assert(false, "should not iterate");
}

unittest
{
    // byKeyValue on empty map yields nothing
    HashMap!(int, int) map;
    foreach (kv; map.byKeyValue())
        assert(false, "should not iterate");
}

unittest
{
    // byKeyValue skips tombstones
    HashMap!(int, int) map;
    map[1] = 10;
    map[2] = 20;
    map[3] = 30;
    map.remove(2);

    int keySum = 0;
    int valSum = 0;
    foreach (kv; map.byKeyValue())
    {
        keySum += kv.key;
        valSum += kv.value;
    }
    assert(keySum == 4);  // 1 + 3
    assert(valSum == 40); // 10 + 30
}

unittest
{
    // Tombstone reuse: re-insert a removed key
    HashMap!(int, int) map;
    map[1] = 10;
    map.remove(1);
    assert(1 !in map);
    assert(map.length == 0);

    map[1] = 20;
    assert(map.length == 1);
    assert(map[1] == 20);
}

unittest
{
    // Insert-remove-insert cycle with multiple keys (tombstone handling)
    HashMap!(int, int) map;
    foreach (i; 0 .. 50)
        map[i] = i;
    foreach (i; 0 .. 50)
        map.remove(i);
    assert(map.length == 0);

    // Reinsert — triggers rehash since tombstones fill the table
    foreach (i; 0 .. 50)
        map[i] = i * 3;
    assert(map.length == 50);
    foreach (i; 0 .. 50)
        assert(map[i] == i * 3);
}

unittest
{
    // Update does not change length
    HashMap!(int, int) map;
    map[1] = 10;
    map[1] = 20;
    map[1] = 30;
    assert(map.length == 1);
    assert(map[1] == 30);
}

unittest
{
    // dup produces independent copy — mutations don't cross
    HashMap!(int, int) map;
    map[1] = 10;
    map[2] = 20;

    auto copy = map.dup;
    copy[3] = 30;
    copy.remove(1);

    assert(map.length == 2);
    assert(map[1] == 10);
    assert(3 !in map);

    assert(copy.length == 2);
    assert(1 !in copy);
    assert(copy[3] == 30);
}

unittest
{
    // String keys
    HashMap!(string, int) map;
    map["hello"] = 1;
    map["world"] = 2;
    assert(map.length == 2);
    assert(map["hello"] == 1);
    assert(map["world"] == 2);

    map.remove("hello");
    assert("hello" !in map);
    assert(map.length == 1);
}

unittest
{
    // Const map: lookup and iteration work
    HashMap!(int, int) map;
    map[1] = 10;
    map[2] = 20;

    const cmap = map;
    assert(cmap.length == 2);
    assert(1 in cmap);
    assert(cmap[1] == 10);

    int sum = 0;
    foreach (kv; cmap.byKeyValue())
        sum += kv.value;
    assert(sum == 30);

    int keySum = 0;
    foreach (k; cmap.byKey())
        keySum += k;
    assert(keySum == 3);
}

unittest
{
    // Negative integer keys
    HashMap!(int, int) map;
    map[-1] = 10;
    map[-100] = 20;
    map[0] = 30;
    assert(map.length == 3);
    assert(map[-1] == 10);
    assert(map[-100] == 20);
    assert(map[0] == 30);
}

unittest
{
    // Large resize: ensure correctness across multiple grow() calls
    HashMap!(int, int) map;
    enum N = 10_000;
    foreach (i; 0 .. N)
        map[i] = i;
    assert(map.length == N);

    // Spot-check values
    assert(map[0] == 0);
    assert(map[N / 2] == N / 2);
    assert(map[N - 1] == N - 1);

    // Remove all
    foreach (i; 0 .. N)
        map.remove(i);
    assert(map.length == 0);
}
