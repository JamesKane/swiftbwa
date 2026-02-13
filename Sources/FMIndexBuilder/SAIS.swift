/// SA-IS (Suffix Array by Induced Sorting) — O(n) suffix array construction.
///
/// Port of sais-lite by Yuta Mori (MIT license).
/// Copyright (c) 2008-2010 Yuta Mori. All Rights Reserved.
/// Modified (C) 2020 Intel Corporation, Heng Li.
///
/// Uses Int64 index type to support genomes exceeding 2^31.
/// Internally works on Int64 text to handle recursive subproblems
/// where "characters" are lexicographic names (> 255).
struct SAIS {

    /// Build a suffix array for `text[0..<n]` with alphabet `[0, alphabetSize)`.
    /// Returns the suffix array as `[Int64]` of length `n`.
    static func buildSuffixArray(
        text: UnsafePointer<UInt8>,
        n: Int64,
        alphabetSize: Int64
    ) -> [Int64] {
        guard n > 1 else {
            if n == 1 { return [0] }
            return []
        }

        // Convert UInt8 text to Int64 for uniform internal processing
        let t = UnsafeMutablePointer<Int64>.allocate(capacity: Int(n))
        defer { t.deallocate() }
        for i in 0..<Int(n) { t[i] = Int64(text[i]) }

        var sa = [Int64](repeating: 0, count: Int(n))
        sa.withUnsafeMutableBufferPointer { saBuf -> Void in
            suffixSort(
                T: t,
                SA: saBuf.baseAddress!,
                fs: 0,
                n: n,
                k: alphabetSize,
                isbwt: false
            )
        }
        return sa
    }

    // MARK: - Private SAIS Implementation

    private static func getCounts(
        _ T: UnsafePointer<Int64>, _ C: UnsafeMutablePointer<Int64>,
        _ n: Int64, _ k: Int64
    ) {
        for i in 0..<Int(k) { C[i] = 0 }
        for i in 0..<Int(n) { C[Int(T[i])] += 1 }
    }

    private static func getBuckets(
        _ C: UnsafePointer<Int64>, _ B: UnsafeMutablePointer<Int64>,
        _ k: Int64, _ end: Bool
    ) {
        var sum: Int64 = 0
        if end {
            for i in 0..<Int(k) { sum += C[i]; B[i] = sum }
        } else {
            for i in 0..<Int(k) { sum += C[i]; B[i] = sum - C[i] }
        }
    }

    private static func lmsSSort1(
        _ T: UnsafePointer<Int64>, _ SA: UnsafeMutablePointer<Int64>,
        _ C: UnsafeMutablePointer<Int64>, _ B: UnsafeMutablePointer<Int64>,
        _ n: Int64, _ k: Int64, _ recount: Bool
    ) {
        if recount { getCounts(T, C, n, k) }
        getBuckets(C, B, k, false)

        var j = n - 1
        var c1 = T[Int(j)]
        var b = B[Int(c1)]
        j -= 1
        SA[Int(b)] = (T[Int(j)] < c1) ? ~j : j
        b += 1

        for i in 0..<Int(n) {
            j = SA[i]
            if j > 0 {
                let c0 = T[Int(j)]
                if c0 != c1 { B[Int(c1)] = b; b = B[Int(c0)]; c1 = c0 }
                j -= 1
                SA[Int(b)] = (T[Int(j)] < c1) ? ~j : j
                b += 1
                SA[i] = 0
            } else if j < 0 {
                SA[i] = ~j
            }
        }

        if recount { getCounts(T, C, n, k) }
        getBuckets(C, B, k, true)

        c1 = 0
        b = B[0]
        for i in stride(from: Int(n) - 1, through: 0, by: -1) {
            j = SA[i]
            if j > 0 {
                let c0 = T[Int(j)]
                if c0 != c1 { B[Int(c1)] = b; b = B[Int(c0)]; c1 = c0 }
                j -= 1
                b -= 1
                SA[Int(b)] = (T[Int(j)] > c1) ? ~(j + 1) : j
                SA[i] = 0
            }
        }
    }

    private static func lmsPostproc1(
        _ T: UnsafePointer<Int64>, _ SA: UnsafeMutablePointer<Int64>,
        _ n: Int64, _ m: Int64
    ) -> Int64 {
        var i: Int64 = 0
        var p: Int64

        // Compact sorted substrings into first m items
        while true {
            p = SA[Int(i)]
            if p >= 0 { break }
            SA[Int(i)] = ~p
            i += 1
        }
        if i < m {
            var j2 = i
            i += 1
            while true {
                p = SA[Int(i)]
                if p < 0 {
                    SA[Int(j2)] = ~p
                    SA[Int(i)] = 0
                    j2 += 1
                    if j2 == m { break }
                }
                i += 1
            }
        }

        // Store lengths of all substrings
        i = n - 1
        var j = n - 1
        var c0 = T[Int(n - 1)]
        var c1: Int64
        repeat { c1 = c0; i -= 1; guard i >= 0 else { break }; c0 = T[Int(i)] } while c0 >= c1
        while i >= 0 {
            repeat { c1 = c0; i -= 1; guard i >= 0 else { break }; c0 = T[Int(i)] } while c0 <= c1
            if i >= 0 {
                SA[Int(m + ((i + 1) >> 1))] = j - i
                j = i + 1
                repeat { c1 = c0; i -= 1; guard i >= 0 else { break }; c0 = T[Int(i)] } while c0 >= c1
            }
        }

        // Find lexicographic names
        var name: Int64 = 0
        var q: Int64 = n
        var qlen: Int64 = 0
        for idx in 0..<Int(m) {
            p = SA[idx]
            let plen = SA[Int(m + (p >> 1))]
            var diff = true
            if plen == qlen && (q + plen) < n {
                var jj: Int64 = 0
                while jj < plen && T[Int(p + jj)] == T[Int(q + jj)] { jj += 1 }
                if jj == plen { diff = false }
            }
            if diff { name += 1; q = p; qlen = plen }
            SA[Int(m + (p >> 1))] = name
        }

        return name
    }

    private static func induceSA(
        _ T: UnsafePointer<Int64>, _ SA: UnsafeMutablePointer<Int64>,
        _ C: UnsafeMutablePointer<Int64>, _ B: UnsafeMutablePointer<Int64>,
        _ n: Int64, _ k: Int64, _ recount: Bool
    ) {
        if recount { getCounts(T, C, n, k) }
        getBuckets(C, B, k, false)

        var j = n - 1
        var c1 = T[Int(j)]
        var b = B[Int(c1)]
        SA[Int(b)] = (j > 0 && T[Int(j - 1)] < c1) ? ~j : j
        b += 1

        for i in 0..<Int(n) {
            j = SA[i]
            SA[i] = ~j
            if j > 0 {
                j -= 1
                let c0 = T[Int(j)]
                if c0 != c1 { B[Int(c1)] = b; b = B[Int(c0)]; c1 = c0 }
                SA[Int(b)] = (j > 0 && T[Int(j - 1)] < c1) ? ~j : j
                b += 1
            }
        }

        if recount { getCounts(T, C, n, k) }
        getBuckets(C, B, k, true)

        c1 = 0
        b = B[0]
        for i in stride(from: Int(n) - 1, through: 0, by: -1) {
            j = SA[i]
            if j > 0 {
                j -= 1
                let c0 = T[Int(j)]
                if c0 != c1 { B[Int(c1)] = b; b = B[Int(c0)]; c1 = c0 }
                b -= 1
                SA[Int(b)] = (j == 0 || T[Int(j - 1)] > c1) ? ~j : j
            } else {
                SA[i] = ~j
            }
        }
    }

    private static func stage1Sort(
        _ T: UnsafePointer<Int64>, _ SA: UnsafeMutablePointer<Int64>,
        _ C: UnsafeMutablePointer<Int64>, _ B: UnsafeMutablePointer<Int64>,
        _ n: Int64, _ k: Int64, _ flags: UInt32
    ) -> (m: Int64, name: Int64) {
        getCounts(T, C, n, k)
        getBuckets(C, B, k, true)

        for i in 0..<Int(n) { SA[i] = 0 }

        var bIdx = n - 1
        var i = n - 1
        var j = n
        var m: Int64 = 0
        var c0 = T[Int(n - 1)]
        var c1: Int64

        repeat { c1 = c0; i -= 1; guard i >= 0 else { break }; c0 = T[Int(i)] } while c0 >= c1

        while i >= 0 {
            repeat { c1 = c0; i -= 1; guard i >= 0 else { break }; c0 = T[Int(i)] } while c0 <= c1
            if i >= 0 {
                SA[Int(bIdx)] = j
                B[Int(c1)] -= 1
                bIdx = B[Int(c1)]
                j = i
                m += 1
                repeat { c1 = c0; i -= 1; guard i >= 0 else { break }; c0 = T[Int(i)] } while c0 >= c1
            }
        }
        SA[Int(n - 1)] = 0

        let name: Int64
        if m > 1 {
            lmsSSort1(T, SA, C, B, n, k, (flags & (4 | 64)) != 0)
            name = lmsPostproc1(T, SA, n, m)
        } else if m == 1 {
            SA[Int(bIdx)] = j + 1
            name = 1
        } else {
            name = 0
        }

        return (m, name)
    }

    private static func stage3Sort(
        _ T: UnsafePointer<Int64>, _ SA: UnsafeMutablePointer<Int64>,
        _ C: UnsafeMutablePointer<Int64>, _ B: UnsafeMutablePointer<Int64>,
        _ n: Int64, _ m: Int64, _ k: Int64, _ flags: UInt32
    ) {
        if (flags & 8) != 0 { getCounts(T, C, n, k) }

        if m > 1 {
            getBuckets(C, B, k, true)
            var i2 = m - 1
            var j2 = n
            var p = SA[Int(m - 1)]
            var c1 = T[Int(p)]
            var c0: Int64

            while true {
                c0 = c1
                let q = B[Int(c0)]
                while q < j2 { j2 -= 1; SA[Int(j2)] = 0 }
                while true {
                    j2 -= 1
                    SA[Int(j2)] = p
                    i2 -= 1
                    if i2 < 0 { break }
                    p = SA[Int(i2)]
                    c1 = T[Int(p)]
                    if c1 != c0 { break }
                }
                if i2 < 0 { break }
            }
            while j2 > 0 { j2 -= 1; SA[Int(j2)] = 0 }
        }

        induceSA(T, SA, C, B, n, k, (flags & (4 | 64)) != 0)
    }

    @discardableResult
    private static func suffixSort(
        T: UnsafePointer<Int64>,
        SA: UnsafeMutablePointer<Int64>,
        fs: Int64,
        n: Int64,
        k: Int64,
        isbwt: Bool
    ) -> Int64 {
        var C: UnsafeMutablePointer<Int64>?
        var B: UnsafeMutablePointer<Int64>?
        var flags: UInt32 = 0

        if k <= 256 {
            C = .allocate(capacity: Int(k))
            if k <= fs {
                B = SA + Int(n + fs - k)
                flags = 1
            } else {
                B = .allocate(capacity: Int(k))
                flags = 3
            }
        } else if k <= fs {
            C = SA + Int(n + fs - k)
            if k <= (fs - k) {
                B = C! - Int(k)
                flags = 0
            } else {
                B = .allocate(capacity: Int(k))
                flags = 2
            }
        } else {
            C = .allocate(capacity: Int(k))
            B = C
            flags = 4 | 8
        }

        let result = stage1Sort(T, SA, C!, B!, n, k, flags)
        let m = result.m
        let name = result.name

        if m < 0 {
            if flags & (1 | 4) != 0 { C?.deallocate() }
            if flags & 2 != 0 { B?.deallocate() }
            return -2
        }

        // Stage 2: recursively sort reduced problem if names not unique
        if name < m {
            if flags & 4 != 0 { C?.deallocate() }
            if flags & 2 != 0 { B?.deallocate() }
            var newfs = (n + fs) - (m * 2)
            if (flags & (1 | 4 | 64)) == 0 {
                if (k + name) <= newfs { newfs -= k }
                else { flags |= 8 }
            }

            let RA = SA + Int(m + newfs)
            var j2 = m - 1
            for i2 in stride(from: Int(m + (n >> 1) - 1), through: Int(m), by: -1) {
                if SA[i2] != 0 { RA[Int(j2)] = SA[i2] - 1; j2 -= 1 }
            }

            // Recursive call: RA is Int64* text with alphabet [0, name)
            suffixSort(T: UnsafePointer(RA), SA: SA, fs: newfs, n: m, k: name, isbwt: false)

            // Re-scan LMS positions from original T
            var i3 = n - 1
            j2 = m - 1
            var c0_3 = T[Int(n - 1)]
            var c1_3: Int64
            repeat { c1_3 = c0_3; i3 -= 1; guard i3 >= 0 else { break }; c0_3 = T[Int(i3)] } while c0_3 >= c1_3
            while i3 >= 0 {
                repeat { c1_3 = c0_3; i3 -= 1; guard i3 >= 0 else { break }; c0_3 = T[Int(i3)] } while c0_3 <= c1_3
                if i3 >= 0 {
                    RA[Int(j2)] = i3 + 1
                    j2 -= 1
                    repeat { c1_3 = c0_3; i3 -= 1; guard i3 >= 0 else { break }; c0_3 = T[Int(i3)] } while c0_3 >= c1_3
                }
            }
            for i4 in 0..<Int(m) { SA[i4] = RA[Int(SA[i4])] }

            if flags & 4 != 0 {
                C = .allocate(capacity: Int(k))
                B = C
            }
            if flags & 2 != 0 {
                B = .allocate(capacity: Int(k))
            }
        }

        // Stage 3: induce the result for the original problem
        stage3Sort(T, SA, C!, B!, n, m, k, flags)

        if flags & (1 | 4) != 0 { C?.deallocate() }
        if flags & 2 != 0 { B?.deallocate() }

        return 0
    }
}
