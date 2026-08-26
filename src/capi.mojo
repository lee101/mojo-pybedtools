"""Compiled sweep-line kernels for BED's half-open interval semantics."""

comptime IPtr = UnsafePointer[Int64, AnyOrigin[mut=True]]


def valid_overlap(
    ac: Int64, ast: Int64, aen: Int64, astr: Int64,
    bc: Int64, bst: Int64, ben: Int64, bstr: Int64,
    min_a: Float64, min_b: Float64, reciprocal: Int, either: Int, strand_mode: Int,
) -> Bool:
    if ac != bc or aen <= bst or ben <= ast:
        return False
    if strand_mode == 1 and (astr == 0 or bstr == 0 or astr != bstr):
        return False
    if strand_mode == 2 and (astr == 0 or bstr == 0 or astr == bstr):
        return False
    var lo = ast
    if bst > lo:
        lo = bst
    var hi = aen
    if ben < hi:
        hi = ben
    var overlap = hi - lo
    var pass_a = Float64(overlap) >= min_a * Float64(aen - ast)
    var pass_b = Float64(overlap) >= min_b * Float64(ben - bst)
    if reciprocal != 0:
        pass_b = Float64(overlap) >= min_a * Float64(ben - bst)
    if either != 0:
        return pass_a or pass_b
    return pass_a and pass_b


@export("mpbt_intersect_count")
def mpbt_intersect_count(
    achrom: Int, astart: Int, aend: Int, astrand: Int,
    bchrom: Int, bstart: Int, bend: Int, bstrand: Int,
    na: Int, nb: Int, min_a: Float64, min_b: Float64,
    reciprocal: Int, either: Int, strand_mode: Int,
) abi("C") -> Int:
    var ac = IPtr(unsafe_from_address=achrom)
    var ast = IPtr(unsafe_from_address=astart)
    var aen = IPtr(unsafe_from_address=aend)
    var ass = IPtr(unsafe_from_address=astrand)
    var bc = IPtr(unsafe_from_address=bchrom)
    var bst = IPtr(unsafe_from_address=bstart)
    var ben = IPtr(unsafe_from_address=bend)
    var bss = IPtr(unsafe_from_address=bstrand)
    var total = 0
    var first = 0
    for i in range(na):
        while first < nb and (bc[first] < ac[i] or (bc[first] == ac[i] and ben[first] <= ast[i])):
            first += 1
        for j in range(first, nb):
            if bc[j] != ac[i] or bst[j] >= aen[i]:
                break
            if valid_overlap(ac[i], ast[i], aen[i], ass[i], bc[j], bst[j], ben[j], bss[j], min_a, min_b, reciprocal, either, strand_mode):
                total += 1
    return total


@export("mpbt_intersect_pairs")
def mpbt_intersect_pairs(
    achrom: Int, astart: Int, aend: Int, astrand: Int,
    bchrom: Int, bstart: Int, bend: Int, bstrand: Int,
    na: Int, nb: Int, min_a: Float64, min_b: Float64,
    reciprocal: Int, either: Int, strand_mode: Int, left: Int, right: Int,
) abi("C") -> Int:
    var ac = IPtr(unsafe_from_address=achrom)
    var ast = IPtr(unsafe_from_address=astart)
    var aen = IPtr(unsafe_from_address=aend)
    var ass = IPtr(unsafe_from_address=astrand)
    var bc = IPtr(unsafe_from_address=bchrom)
    var bst = IPtr(unsafe_from_address=bstart)
    var ben = IPtr(unsafe_from_address=bend)
    var bss = IPtr(unsafe_from_address=bstrand)
    var ol = IPtr(unsafe_from_address=left)
    var orr = IPtr(unsafe_from_address=right)
    var total = 0
    var first = 0
    for i in range(na):
        while first < nb and (bc[first] < ac[i] or (bc[first] == ac[i] and ben[first] <= ast[i])):
            first += 1
        for j in range(first, nb):
            if bc[j] != ac[i] or bst[j] >= aen[i]:
                break
            if valid_overlap(ac[i], ast[i], aen[i], ass[i], bc[j], bst[j], ben[j], bss[j], min_a, min_b, reciprocal, either, strand_mode):
                ol[total] = Int64(i)
                orr[total] = Int64(j)
                total += 1
    return total


@export("mpbt_intersect_counts")
def mpbt_intersect_counts(
    achrom: Int, astart: Int, aend: Int, astrand: Int,
    bchrom: Int, bstart: Int, bend: Int, bstrand: Int,
    na: Int, nb: Int, min_a: Float64, min_b: Float64,
    reciprocal: Int, either: Int, strand_mode: Int, counts: Int,
) abi("C"):
    var ac = IPtr(unsafe_from_address=achrom)
    var ast = IPtr(unsafe_from_address=astart)
    var aen = IPtr(unsafe_from_address=aend)
    var ass = IPtr(unsafe_from_address=astrand)
    var bc = IPtr(unsafe_from_address=bchrom)
    var bst = IPtr(unsafe_from_address=bstart)
    var ben = IPtr(unsafe_from_address=bend)
    var bss = IPtr(unsafe_from_address=bstrand)
    var oc = IPtr(unsafe_from_address=counts)
    var first = 0
    for i in range(na):
        var count = Int64(0)
        while first < nb and (bc[first] < ac[i] or (bc[first] == ac[i] and ben[first] <= ast[i])):
            first += 1
        for j in range(first, nb):
            if bc[j] != ac[i] or bst[j] >= aen[i]:
                break
            if valid_overlap(ac[i], ast[i], aen[i], ass[i], bc[j], bst[j], ben[j], bss[j], min_a, min_b, reciprocal, either, strand_mode):
                count += 1
        oc[i] = count


@export("mpbt_coverage")
def mpbt_coverage(
    achrom: Int, astart: Int, aend: Int, bchrom: Int, bstart: Int, bend: Int,
    aorder: Int, na: Int, nb: Int, bases: Int, counts: Int,
) abi("C"):
    var ac = IPtr(unsafe_from_address=achrom)
    var ast = IPtr(unsafe_from_address=astart)
    var aen = IPtr(unsafe_from_address=aend)
    var bc = IPtr(unsafe_from_address=bchrom)
    var bst = IPtr(unsafe_from_address=bstart)
    var ben = IPtr(unsafe_from_address=bend)
    var order = IPtr(unsafe_from_address=aorder)
    var ob = IPtr(unsafe_from_address=bases)
    var oc = IPtr(unsafe_from_address=counts)
    var first = 0
    for i in range(na):
        var sum = Int64(0)
        var hits = Int64(0)
        var covered_until = ast[i]
        while first < nb and (bc[first] < ac[i] or (bc[first] == ac[i] and ben[first] <= ast[i])):
            first += 1
        for j in range(first, nb):
            if bc[j] != ac[i] or bst[j] >= aen[i]:
                break
            if bc[j] == ac[i] and ben[j] > ast[i]:
                var lo = ast[i]
                if bst[j] > lo:
                    lo = bst[j]
                var hi = aen[i]
                if ben[j] < hi:
                    hi = ben[j]
                if hi > lo:
                    hits += 1
                    if hi > covered_until:
                        var uncovered_start = lo
                        if uncovered_start < covered_until:
                            uncovered_start = covered_until
                        sum += hi - uncovered_start
                        covered_until = hi
        ob[order[i]] = sum
        oc[order[i]] = hits


@export("mpbt_merge")
def mpbt_merge(chrom: Int, start: Int, end: Int, n: Int, distance: Int, out_chrom: Int, out_start: Int, out_end: Int) abi("C") -> Int:
    var c = IPtr(unsafe_from_address=chrom)
    var st = IPtr(unsafe_from_address=start)
    var en = IPtr(unsafe_from_address=end)
    var oc = IPtr(unsafe_from_address=out_chrom)
    var os = IPtr(unsafe_from_address=out_start)
    var oe = IPtr(unsafe_from_address=out_end)
    if n == 0:
        return 0
    var k = 0
    var current_c = c[0]
    var current_s = st[0]
    var current_e = en[0]
    for i in range(1, n):
        if c[i] == current_c and st[i] <= current_e + Int64(distance):
            if en[i] > current_e:
                current_e = en[i]
        else:
            oc[k] = current_c
            os[k] = current_s
            oe[k] = current_e
            k += 1
            current_c = c[i]
            current_s = st[i]
            current_e = en[i]
    oc[k] = current_c
    os[k] = current_s
    oe[k] = current_e
    return k + 1
