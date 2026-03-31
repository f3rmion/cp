const std = @import("std");
const Managed = std.math.big.int.Managed;
const Limb = std.math.big.Limb;
const Log2Limb = std.math.big.Log2Limb;

const limb_bits = @typeInfo(Limb).int.bits;

/// Modular exponentiation: result = base^exp mod m
/// Uses right-to-left binary square-and-multiply
pub fn powMod(
    result: *Managed,
    base: *const Managed,
    exp: *const Managed,
    m: *const Managed,
) !void {
    const allocator = result.allocator;

    // Temporaries
    var q = try Managed.init(allocator);
    defer q.deinit();
    var rem = try Managed.init(allocator);
    defer rem.deinit();
    var tmp = try Managed.init(allocator);
    defer tmp.deinit();

    // b = base % m
    var b = try Managed.init(allocator);
    defer b.deinit();
    try q.divFloor(&b, base, m);

    // result = 1
    try result.set(1);
    const bit_count = exp.bitCountAbs();
    for (0..bit_count) |i| {
        // Check if bit i of exp is set
        const limb_index = i / limb_bits;
        const bit_pos: Log2Limb = @truncate(i);
        const is_set = (exp.limbs[limb_index] >> bit_pos) & 1 != 0;

        if (is_set) {
            // result = (result * b) % m
            try tmp.mul(result, &b);
            try q.divFloor(result, &tmp, m);
        }

        // b = (b * b) % m
        try tmp.mul(&b, &b);
        try q.divFloor(&b, &tmp, m);
    }
}

test "powMod basic" {
    const ally = std.testing.allocator;

    // 3^13 mod 7 = 3
    var base = try Managed.initSet(ally, 3);
    defer base.deinit();
    var exp = try Managed.initSet(ally, 13);
    defer exp.deinit();
    var m = try Managed.initSet(ally, 7);
    defer m.deinit();
    var result = try Managed.init(ally);
    defer result.deinit();

    try powMod(&result, &base, &exp, &m);

    const val = try result.toInt(u64);
    try std.testing.expectEqual(@as(u64, 3), val);
}

/// Returns true if n is probably prime
/// Uses deterministic witness sufficient for number uo to 2^64,
/// and probabilistic witnesses for larger numbers.
pub fn isPrime(n: *const Managed) !bool {
    const allocator = n.allocator;

    // Handle small cases
    if (!n.isPositive()) return false;
    const cmp_2 = n.toConst().orderAgainstScalar(2);
    if (cmp_2 == .lt) return false; // n < 2
    if (cmp_2 == .eq) return true; // n == 2

    // Write n-1 = 2^r * d
    var n_minus_1 = try Managed.init(allocator);
    defer n_minus_1.deinit();
    var one = try Managed.initSet(allocator, 1);
    defer one.deinit();
    try n_minus_1.sub(n, &one);

    // Find r and d: n-1 = 2^r * d where d is odd
    var d = try n_minus_1.clone();
    defer d.deinit();
    var r: usize = 0;
    while (!d.isOdd()) {
        try d.shiftRight(&d, 1);
        r += 1;
    }

    // Deterministic witnesses for small numbers, probabilistic for large
    const small_witnesses = [_]u64{ 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37 };

    var a = try Managed.init(allocator);
    defer a.deinit();
    var x = try Managed.init(allocator);
    defer x.deinit();
    var n_minus_2 = try Managed.init(allocator);
    defer n_minus_2.deinit();
    var two = try Managed.initSet(allocator, 2);
    defer two.deinit();
    try n_minus_2.sub(n, &two);

    for (small_witnesses) |w| {
        // Skip witness if a >= n-1
        try a.set(w);
        if (a.order(n_minus_1) != .lt) continue;

        // x = a^d mod n
        try powMod(&x, &a, &d, n);

        if (x.toConst().orderAgainstScalar(1) == .eq) continue;
        if (x.eql(n_minus_1)) continue;

        var composite = true;
        for (0..r -| 1) |_| {
            // x = x^2 mod n
            try powMod(&x, &x, &two, n);
            if (x.eql(n_minus_1)) {
                composite = false;
                break;
            }
        }

        if (composite) return false;
    }

    return true;
}

test "isPrime" {
    const ally = std.testing.allocator;

    var p = try Managed.initSet(ally, 104729);
    defer p.deinit();
    try std.testing.expect(try isPrime(&p));

    var c = try Managed.initSet(ally, 104730);
    defer c.deinit();
    try std.testing.expect(!try isPrime(&c));

    // Mersenne prime 2^61 - 1
    var m61 = try Managed.initSet(ally, 2305843009213693951);
    defer m61.deinit();
    try std.testing.expect(try isPrime(&m61));
}
