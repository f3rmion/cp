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

/// Legendre symbol: returns 1, -1, or 0
pub fn legendreSymbol(a: *const Managed, p: *const Managed) !i2 {
    const allocator = a.allocator;

    // Compute a^((p-1)/2) mod p
    var one = try Managed.initSet(allocator, 1);
    defer one.deinit();

    var p_minus_1 = try Managed.init(allocator);
    defer p_minus_1.deinit();
    try p_minus_1.sub(p, &one);

    var exp = try Managed.init(allocator);
    defer exp.deinit();
    try exp.shiftRight(&p_minus_1, 1); // (p-1)/2

    var result = try Managed.init(allocator);
    defer result.deinit();
    try powMod(&result, a, &exp, p);

    if (result.eqlZero()) return 0;
    if (result.toConst().orderAgainstScalar(1) == .eq) return 1;
    return -1;
}

test "legendreSymbol" {
    const ally = std.testing.allocator;

    var a = try Managed.initSet(ally, 4);
    defer a.deinit();
    var p = try Managed.initSet(ally, 7);
    defer p.deinit();
    try std.testing.expectEqual(@as(i2, 1), try legendreSymbol(&a, &p));

    try a.set(3);
    try std.testing.expectEqual(@as(i2, -1), try legendreSymbol(&a, &p));
}

/// Modular square root via Tonelli-Shanks. Returns x where x^2 == a (mod p)
/// Caller must verify a is a QR first (legendreSymbol == 1).
pub fn sqrtMod(result: *Managed, a: *const Managed, p: *const Managed) !void {
    const allocator = result.allocator;

    var one = try Managed.initSet(allocator, 1);
    defer one.deinit();
    var two = try Managed.initSet(allocator, 2);
    defer two.deinit();

    // Factor out powers of 2: p-1 == 2^s * q
    var p_minus_1 = try Managed.init(allocator);
    defer p_minus_1.deinit();
    try p_minus_1.sub(p, &one);

    var q = try p_minus_1.clone();
    defer q.deinit();
    var s: usize = 0;
    while (!q.isOdd()) {
        try q.shiftRight(&q, 1);
        s += 1;
    }

    // Find a quadratic non-residue z
    var z = try Managed.initSet(allocator, 2);
    defer z.deinit();
    while (try legendreSymbol(&z, p) != -1) {
        try z.add(&z, &one);
    }

    var m: usize = s;
    var c = try Managed.init(allocator);
    defer c.deinit();
    try powMod(&c, &z, &q, p); // c = z^q mod p

    var t = try Managed.init(allocator);
    defer t.deinit();
    try powMod(&t, a, &q, p); // t = a^q mod p

    // r = a^((q+1)/2) mod p
    var q_plus_1 = try Managed.init(allocator);
    defer q_plus_1.deinit();
    try q_plus_1.add(&q, &one);
    var exp = try Managed.init(allocator);
    defer exp.deinit();
    try exp.shiftRight(&q_plus_1, 1);
    try powMod(result, a, &exp, p);

    // Temporaries for the loop
    var tmp = try Managed.init(allocator);
    defer tmp.deinit();
    var qr = try Managed.init(allocator);
    defer qr.deinit();
    var b = try Managed.init(allocator);
    defer b.deinit();
    var pow_exp = try Managed.init(allocator);
    defer pow_exp.deinit();

    while (true) {
        if (t.toConst().orderAgainstScalar(1) == .eq) break;

        // Find least i such that t^(2^i) == 1 mod p
        var i: usize = 1;
        try tmp.copy(t.toConst());
        while (i < m) {
            try powMod(&tmp, &tmp, &tmp, p);
            if (tmp.toConst().orderAgainstScalar(1) == .eq) break;
            i += 1;
        }

        // b = c^(2^(m-i-1)) mod p
        try pow_exp.set(1);
        try pow_exp.shiftLeft(&pow_exp, m - i - 1);
        try powMod(&b, &c, &two, p);

        m = i;

        // c = b^2 mod p
        try powMod(&c, &b, &two, p);

        // t = t * c mod p
        try tmp.mul(&t, &c);
        try qr.divFloor(&t, &tmp, p);

        // result = result * b mod p
        try tmp.mul(result, &b);
        try qr.divFloor(result, &tmp, p);
    }
}

test "sqrtMod" {
    const ally = std.testing.allocator;

    // sqrt(4) mod 7 should be 2 or 5
    var a = try Managed.initSet(ally, 4);
    defer a.deinit();
    var p = try Managed.initSet(ally, 7);
    defer p.deinit();
    var result = try Managed.init(ally);
    defer result.deinit();
    try sqrtMod(&result, &a, &p);

    // Verify: result^2 mod p == a
    var two = try Managed.initSet(ally, 2);
    defer two.deinit();
    var check = try Managed.init(ally);
    defer check.deinit();
    try powMod(&check, &result, &two, &p);
    try std.testing.expect(check.toConst().orderAgainstScalar(4) == .eq);
}

/// Find z such that z^k == 1 mod r and z != 1
/// r must be prime and k must divide r-1.
pub fn findKthRootOfUnity(result: *Managed, k: u64, r: *const Managed) !void {
    const allocator = result.allocator;

    var one = try Managed.initSet(allocator, 1);
    defer one.deinit();
    var r_minus_1 = try Managed.init(allocator);
    defer r_minus_1.deinit();
    try r_minus_1.sub(r, &one);

    // exp = (r-1) / k
    var k_big = try Managed.initSet(allocator, k);
    defer k_big.deinit();
    var exp = try Managed.init(allocator);
    defer exp.deinit();
    var rem = try Managed.init(allocator);
    defer rem.deinit();
    try exp.divFloor(&rem, &r_minus_1, &k_big);

    // Try g = 2, 3, 4, ... until g^((r-1)/k) != 1 mod r
    var g = try Managed.initSet(allocator, 2);
    defer g.deinit();

    while (true) {
        try powMod(result, &g, &exp, r);
        if (result.toConst().orderAgainstScalar(1) != .eq) return;
        try g.add(&g, &one);
    }
}

test "findKthRootOfUnity" {
    const ally = std.testing.allocator;

    // r = 13, k = 3. 3 | 12 so a cube root or unity exists.
    var r = try Managed.initSet(ally, 13);
    defer r.deinit();
    var z = try Managed.init(ally);
    defer z.deinit();
    try findKthRootOfUnity(&z, 3, &r);

    // Verify z^3 == 1 mod 13
    var three = try Managed.initSet(ally, 3);
    defer three.deinit();
    var check = try Managed.init(ally);
    defer check.deinit();
    try powMod(&check, &z, &three, &r);
    try std.testing.expect(check.toConst().orderAgainstScalar(1) == .eq);

    // And z != 1
    try std.testing.expect(z.toConst().orderAgainstScalar(1) != .eq);
}
