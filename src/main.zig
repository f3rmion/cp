const std = @import("std");
const cp = @import("cocks_pinch");
const Managed = std.math.big.int.Managed;

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const ally = gpa.allocator();

    // secp256k1 group order
    var r = try Managed.init(ally);
    defer r.deinit();
    try r.setString(16, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141");

    const r_str = try r.toString(ally, 10, .lower);
    defer ally.free(r_str);
    std.debug.print("r = {s}\n", .{r_str});

    // Target k=6: F_{p^6} ≈ 2^3072 → ~128-bit security
    const k_values = [_]u64{6};

    var result = try cp.cocksPinchCM(ally, &r, &k_values);
    defer result.deinit();

    // Print parameters
    const p_str = try result.p.toString(ally, 10, .lower);
    defer ally.free(p_str);
    const t_str = try result.t.toString(ally, 10, .lower);
    defer ally.free(t_str);
    const y_str = try result.y.toString(ally, 10, .lower);
    defer ally.free(y_str);

    std.debug.print("=== Pairing-friendly curve found ===\n", .{});
    std.debug.print("k = {}\n", .{result.k});
    std.debug.print("D = {}\n", .{result.D});
    std.debug.print("j = {}\n", .{result.j});
    std.debug.print("p = {s}\n", .{p_str});
    std.debug.print("t = {s}\n", .{t_str});
    std.debug.print("y = {s}\n", .{y_str});
    std.debug.print("p bits = {}\n", .{result.p.bitCountAbs()});

    // Construct curve E: y^2 = x^3 + ax + b
    var a = try Managed.init(ally);
    defer a.deinit();
    var b = try Managed.init(ally);
    defer b.deinit();
    try cp.constructCurveFromJ(&a, &b, result.j, &result.p);

    const a_str = try a.toString(ally, 10, .lower);
    defer ally.free(a_str);
    const b_str = try b.toString(ally, 10, .lower);
    defer ally.free(b_str);
    std.debug.print("\nE: y^2 = x^3 + {s}*x + {s}\n", .{ a_str, b_str });

    // Verify r | p+1-t
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

    try tmp.add(&result.p, &one);
    try check.sub(&tmp, &result.t);
    try q.divFloor(&rem, &check, &result.r);

    std.debug.print("\n=== Verification ===\n", .{});
    std.debug.print("r | p+1-t: {}\n", .{rem.eqlZero()});
    std.debug.print("p is prime: {}\n", .{try cp.isPrime(&result.p)});
}
