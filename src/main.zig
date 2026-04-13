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

    const p_bits = result.p.bitCountAbs();
    const r_bits = result.r.bitCountAbs();
    const verified_r_divides = rem.eqlZero();
    const verified_p_prime = try cp.isPrime(&result.p);

    // === Build curve.txt content ===
    const txt = try std.fmt.allocPrint(ally,
        \\========================================
        \\  Pairing-Friendly Elliptic Curve
        \\  Cocks-Pinch with CM construction
        \\========================================
        \\
        \\--- Cocks-Pinch parameters ---
        \\  embedding degree k = {[k]}
        \\  CM discriminant  D = {[D]}
        \\  j-invariant      j = {[j]}
        \\  Frobenius trace  t = {[t_dec]s}
        \\                     = 0x{[t_hex]s}
        \\
        \\--- Curve: E(F_p): y^2 = x^3 + a*x + b ---
        \\
        \\  p  = {[p_dec]s}
        \\     = 0x{[p_hex]s}
        \\     ({[p_bits]} bits)
        \\
        \\  a  = {[a_dec]s}
        \\     = 0x{[a_hex]s}
        \\
        \\  b  = {[b_dec]s}
        \\     = 0x{[b_hex]s}
        \\
        \\--- Group structure ---
        \\
        \\  #E(F_p) = p + 1 - t
        \\  n  = {[n_dec]s}
        \\     = 0x{[n_hex]s}
        \\
        \\  subgroup order (secp256k1)
        \\  r  = {[r_dec]s}
        \\     = 0x{[r_hex]s}
        \\     ({[r_bits]} bits)
        \\
        \\  cofactor h = n / r
        \\  h  = {[h_dec]s}
        \\     = 0x{[h_hex]s}
        \\
        \\--- Security ---
        \\  embedding degree k = {[k2]}
        \\  F_p^k size    ≈ 2^{[ext_bits]}
        \\  rho security  ≈ {[rho_bits]}-bit (ECDLP)
        \\  MOV security  ≈ NFS on {[mov_bits]}-bit field
        \\
        \\--- Verification ---
        \\  r | #E(F_p):  {[v_r]}
        \\  p is prime:   {[v_p]}
        \\========================================
        \\
    , .{
        .k = result.k,
        .D = result.D,
        .j = result.j,
        .t_dec = t_str,
        .t_hex = t_hex,
        .p_dec = p_str,
        .p_hex = p_hex,
        .p_bits = p_bits,
        .a_dec = a_str,
        .a_hex = a_hex,
        .b_dec = b_str,
        .b_hex = b_hex,
        .n_dec = n_str,
        .n_hex = n_hex,
        .r_dec = r_str,
        .r_hex = r_hex,
        .r_bits = r_bits,
        .h_dec = h_str,
        .h_hex = h_hex,
        .k2 = result.k,
        .ext_bits = p_bits * result.k,
        .rho_bits = r_bits / 2,
        .mov_bits = p_bits * result.k,
        .v_r = verified_r_divides,
        .v_p = verified_p_prime,
    });
    defer ally.free(txt);

    // === Build curve.json content ===
    const json = try std.fmt.allocPrint(ally,
        \\{{
        \\  "curve": {{
        \\    "equation": "y^2 = x^3 + a*x + b",
        \\    "p": {{
        \\      "decimal": "{[p_dec]s}",
        \\      "hex": "0x{[p_hex]s}",
        \\      "bits": {[p_bits]}
        \\    }},
        \\    "a": {{
        \\      "decimal": "{[a_dec]s}",
        \\      "hex": "0x{[a_hex]s}"
        \\    }},
        \\    "b": {{
        \\      "decimal": "{[b_dec]s}",
        \\      "hex": "0x{[b_hex]s}"
        \\    }}
        \\  }},
        \\  "group": {{
        \\    "order": {{
        \\      "decimal": "{[n_dec]s}",
        \\      "hex": "0x{[n_hex]s}"
        \\    }},
        \\    "subgroup_order": {{
        \\      "decimal": "{[r_dec]s}",
        \\      "hex": "0x{[r_hex]s}",
        \\      "bits": {[r_bits]},
        \\      "source": "secp256k1"
        \\    }},
        \\    "cofactor": {{
        \\      "decimal": "{[h_dec]s}",
        \\      "hex": "0x{[h_hex]s}"
        \\    }}
        \\  }},
        \\  "pairing": {{
        \\    "embedding_degree": {[k]},
        \\    "cm_discriminant": {[D]},
        \\    "j_invariant": {[j]},
        \\    "frobenius_trace": {{
        \\      "decimal": "{[t_dec]s}",
        \\      "hex": "0x{[t_hex]s}"
        \\    }}
        \\  }},
        \\  "security": {{
        \\    "extension_field_bits": {[ext_bits]},
        \\    "rho_security_bits": {[rho_bits]},
        \\    "mov_field_bits": {[mov_bits]}
        \\  }},
        \\  "verification": {{
        \\    "r_divides_group_order": {[v_r]},
        \\    "p_is_prime": {[v_p]}
        \\  }}
        \\}}
        \\
    , .{
        .p_dec = p_str,
        .p_hex = p_hex,
        .p_bits = p_bits,
        .a_dec = a_str,
        .a_hex = a_hex,
        .b_dec = b_str,
        .b_hex = b_hex,
        .n_dec = n_str,
        .n_hex = n_hex,
        .r_dec = r_str,
        .r_hex = r_hex,
        .r_bits = r_bits,
        .h_dec = h_str,
        .h_hex = h_hex,
        .k = result.k,
        .D = result.D,
        .j = result.j,
        .t_dec = t_str,
        .t_hex = t_hex,
        .ext_bits = p_bits * result.k,
        .rho_bits = r_bits / 2,
        .mov_bits = p_bits * result.k,
        .v_r = verified_r_divides,
        .v_p = verified_p_prime,
    });
    defer ally.free(json);

    // === Write files ===
    {
        const file = try std.fs.cwd().createFile("curve.txt", .{});
        defer file.close();
        try file.writeAll(txt);
    }
    {
        const file = try std.fs.cwd().createFile("curve.json", .{});
        defer file.close();
        try file.writeAll(json);
    }

    // Also print to stderr
    std.debug.print("{s}", .{txt});
    std.debug.print("\nWrote curve.txt and curve.json\n", .{});
}
