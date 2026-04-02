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
            try powMod(&tmp, &tmp, &two, p);
            if (tmp.toConst().orderAgainstScalar(1) == .eq) break;
            i += 1;
        }

        // b = c^(2^(m-i-1)) mod p
        try pow_exp.set(1);
        try pow_exp.shiftLeft(&pow_exp, m - i - 1);
        try powMod(&b, &c, &pow_exp, p);

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

/// Modular inverse: result = a^(-1) mod m
/// Uses Fermat's little theorem: a^(-1) = a^(m-2) mod m (m must be prime)
pub fn modInverse(result: *Managed, a: *const Managed, m: *const Managed) !void {
    const allocator = result.allocator;
    var two = try Managed.initSet(allocator, 2);
    defer two.deinit();
    var exp = try Managed.init(allocator);
    defer exp.deinit();
    try exp.sub(m, &two); // m-2
    try powMod(result, a, &exp, m);
}

/// Result of the Cocks-Pinch method
pub const CocksPinchResult = struct {
    p: Managed,
    r: Managed,
    t: Managed,
    D: u64,
    y: Managed,

    pub fn deinit(self: *CocksPinchResult) void {
        self.p.deinit();
        self.r.deinit();
        self.t.deinit();
        self.y.deinit();
    }
};

/// Cocks-Pinch method
/// Find a prime p and trace t such that:
///  - r is prime, k | r-1
///  - t == z+1 mod r where z is a k-th root of unity
///  - p = (t^2 + D*y^2) / 4 is prime
///  - r | p + 1 - t
pub fn cocksPinch(allocator: std.mem.Allocator, k: u64, r: *const Managed) !CocksPinchResult {
    var one = try Managed.initSet(allocator, 1);
    defer one.deinit();
    var two = try Managed.initSet(allocator, 2);
    defer two.deinit();
    var four = try Managed.initSet(allocator, 4);
    defer four.deinit();

    var q = try Managed.init(allocator);
    defer q.deinit();
    var rem = try Managed.init(allocator);
    defer rem.deinit();
    var tmp = try Managed.init(allocator);
    defer tmp.deinit();

    // Step 1: find k-th root of unity z mod r
    var z = try Managed.init(allocator);
    defer z.deinit();
    try findKthRootOfUnity(&z, k, r);

    // Step 2: t = z + 1 (mod r) - just take the representation in [0, r)
    var t = try Managed.init(allocator);
    try t.add(&z, &one);
    // Reduce mod r (in case z+1 == r)
    try q.divFloor(&rem, &t, r);
    try t.copy(rem.toConst());

    // Step 3: try small CM discriminants D = 1, 2, 3, ...
    // Need: -D is a QR mod r, and p = (t^2 + D*y^2)/4 is prime
    //
    // From the CM equation and r | p+1-t:
    //  4p = t^2 + D*y^2
    //  4(t-1) === t^2 + D*y^2 (mod r) [since p == t-1 mod r]
    //  D*y^2 === -(t-2)^2 (mod r)
    //  y === (t-2) / sqrt(-D) (mod r)

    var t_minus_2 = try Managed.init(allocator);
    defer t_minus_2.deinit();
    try t_minus_2.sub(&t, &two);

    var neg_D = try Managed.init(allocator);
    defer neg_D.deinit();
    var D_big = try Managed.init(allocator);
    defer D_big.deinit();
    var sqrt_neg_D = try Managed.init(allocator);
    defer sqrt_neg_D.deinit();
    var inv = try Managed.init(allocator);
    defer inv.deinit();
    var y = try Managed.init(allocator);
    var p = try Managed.init(allocator);

    var D: u64 = 1;
    while (D < 1000) : (D += 1) {
        // Compute -D mod r
        try D_big.set(D);
        try tmp.sub(r, &D_big); // tmp = r - D
        try q.divFloor(&neg_D, &tmp, r); // neg_D = (-D) mod r

        // Check if -D is a QR mod r
        if (try legendreSymbol(&neg_D, r) != 1) continue;

        // sqrt(-D) mod r
        try sqrtMod(&sqrt_neg_D, &neg_D, r);

        // y = (t - 2) * sqrt(-D)^(-1) mod r
        try modInverse(&inv, &sqrt_neg_D, r);
        try tmp.mul(&t_minus_2, &inv);
        try q.divFloor(&y, &tmp, r);

        // p = (t^2 + D * y^2) / 4
        try tmp.mul(&t, &t); // t^2
        try p.mul(&y, &y); // y^2
        try p.mul(&p, &D_big); // D * y^2
        try p.add(&p, &tmp); // t^2 + D*y^2

        // Check divisible by 4
        try q.divFloor(&rem, &p, &four);
        if (!rem.eqlZero()) continue;
        try p.copy(q.toConst()); // p = (t^2 + D*y^2) / 4

        // Check p is prime
        if (p.toConst().orderAgainstScalar(2) == .lt) continue;
        if (try isPrime(&p)) {
            return CocksPinchResult{
                .p = p,
                .r = try r.toConst().toManaged(allocator),
                .t = t,
                .D = D,
                .y = y,
            };
        }
    }

    // Cleanup
    p.deinit();
    y.deinit();
    t.deinit();
    return error.NoCurveFound;
}

test "cocksPinch" {
    const ally = std.testing.allocator;

    // Use a small prime r where k=6 divides r-1
    // r = 13: r-1 = 12, 6 | 12
    var r = try Managed.initSet(ally, 13);
    defer r.deinit();

    var result = try cocksPinch(ally, 6, &r);
    defer result.deinit();

    // Verify r | p + 1 - t
    var tmp = try Managed.init(ally);
    defer tmp.deinit();
    var one = try Managed.initSet(ally, 1);
    defer one.deinit();
    var check = try Managed.init(ally);
    defer check.deinit();
    var q = try Managed.init(ally);
    defer q.deinit();
    var rem = try Managed.init(ally);
    defer rem.deinit();

    try tmp.add(&result.p, &one); // p+1
    try check.sub(&tmp, &result.t); // p+1-t
    try q.divFloor(&rem, &check, &r);
    try std.testing.expect(rem.eqlZero());

    // Verify p is prime
    try std.testing.expect(try isPrime(&result.p));
}
