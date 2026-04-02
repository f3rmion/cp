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

    // Construct curve E: y^2 = x^3 + ax + b
    var a = try Managed.init(ally);
    defer a.deinit();
    var b = try Managed.init(ally);
    defer b.deinit();
    try cp.constructCurveFromJ(&a, &b, result.j, &result.p);

    // Compute group order n = p + 1 - t and cofactor h = n / r
    var one = try Managed.initSet(ally, 1);
    defer one.deinit();
    var n = try Managed.init(ally);
    defer n.deinit();
    var tmp = try Managed.init(ally);
    defer tmp.deinit();
    try tmp.add(&result.p, &one);
    try n.sub(&tmp, &result.t);

    var h = try Managed.init(ally);
    defer h.deinit();
    var rem = try Managed.init(ally);
    defer rem.deinit();
    try h.divFloor(&rem, &n, &result.r);

    // String conversions (decimal)
    const a_str = try a.toString(ally, 10, .lower);
    defer ally.free(a_str);
    const b_str = try b.toString(ally, 10, .lower);
    defer ally.free(b_str);
    const n_str = try n.toString(ally, 10, .lower);
    defer ally.free(n_str);
    const h_str = try h.toString(ally, 10, .lower);
    defer ally.free(h_str);

    // String conversions (hex)
    const p_hex = try result.p.toString(ally, 16, .lower);
    defer ally.free(p_hex);
    const r_hex = try result.r.toString(ally, 16, .lower);
    defer ally.free(r_hex);
    const a_hex = try a.toString(ally, 16, .lower);
    defer ally.free(a_hex);
    const b_hex = try b.toString(ally, 16, .lower);
    defer ally.free(b_hex);
    const n_hex = try n.toString(ally, 16, .lower);
    defer ally.free(n_hex);
    const h_hex = try h.toString(ally, 16, .lower);
    defer ally.free(h_hex);
    const t_hex = try result.t.toString(ally, 16, .lower);
    defer ally.free(t_hex);

    // Dump
    std.debug.print("========================================\n", .{});
    std.debug.print("  Pairing-Friendly Elliptic Curve\n", .{});
    std.debug.print("========================================\n\n", .{});

    std.debug.print("--- Cocks-Pinch parameters ---\n", .{});
    std.debug.print("  embedding degree k = {}\n", .{result.k});
    std.debug.print("  CM discriminant  D = {}\n", .{result.D});
    std.debug.print("  j-invariant      j = {}\n", .{result.j});
    std.debug.print("  Frobenius trace  t = {s}\n", .{t_str});
    std.debug.print("                     = 0x{s}\n", .{t_hex});

    std.debug.print("\n--- Curve: E(F_p): y^2 = x^3 + a*x + b ---\n", .{});

    std.debug.print("\n  p  = {s}\n", .{p_str});
    std.debug.print("     = 0x{s}\n", .{p_hex});
    std.debug.print("     ({} bits)\n", .{result.p.bitCountAbs()});

    std.debug.print("\n  a  = {s}\n", .{a_str});
    std.debug.print("     = 0x{s}\n", .{a_hex});

    std.debug.print("\n  b  = {s}\n", .{b_str});
    std.debug.print("     = 0x{s}\n", .{b_hex});

    std.debug.print("\n--- Group structure ---\n", .{});

    std.debug.print("\n  #E(F_p) = p + 1 - t\n", .{});
    std.debug.print("  n  = {s}\n", .{n_str});
    std.debug.print("     = 0x{s}\n", .{n_hex});

    std.debug.print("\n  subgroup order (secp256k1)\n", .{});
    std.debug.print("  r  = {s}\n", .{r_str});
    std.debug.print("     = 0x{s}\n", .{r_hex});
    std.debug.print("     ({} bits)\n", .{result.r.bitCountAbs()});

    std.debug.print("\n  cofactor h = n / r\n", .{});
    std.debug.print("  h  = {s}\n", .{h_str});
    std.debug.print("     = 0x{s}\n", .{h_hex});

    std.debug.print("\n--- Security ---\n", .{});
    std.debug.print("  embedding degree k = {}\n", .{result.k});
    std.debug.print("  F_p^k size  ≈ 2^{}\n", .{result.p.bitCountAbs() * result.k});
    std.debug.print("  rho security ≈ {}-bit (ECDLP)\n", .{result.r.bitCountAbs() / 2});
    std.debug.print("  MOV security ≈ NFS on {}-bit field\n", .{result.p.bitCountAbs() * result.k});

    std.debug.print("\n--- Verification ---\n", .{});
    std.debug.print("  r | #E(F_p):  {}\n", .{rem.eqlZero()});
    std.debug.print("  p is prime:   {}\n", .{try cp.isPrime(&result.p)});
    std.debug.print("========================================\n", .{});
}
