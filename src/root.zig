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
