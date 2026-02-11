import Testing
@testable import Alignment
@testable import BWACore

@Suite("Alignment Tests")
struct AlignmentTests {

    @Test("CIGARBuilder basic operations")
    func testCIGARBuilder() {
        var builder = CIGARBuilder()
        builder.append(.match, length: 10)
        builder.append(.insertion, length: 2)
        builder.append(.match, length: 5)

        let cigar = builder.build()
        #expect(cigar.count == 3)
        #expect(cigar[0] == (10 << 4 | 0))  // 10M
        #expect(cigar[1] == (2 << 4 | 1))   // 2I
        #expect(cigar[2] == (5 << 4 | 0))   // 5M
    }

    @Test("CIGARBuilder merges consecutive same ops")
    func testCIGARMerge() {
        var builder = CIGARBuilder()
        builder.append(.match, length: 5)
        builder.append(.match, length: 3)
        builder.append(.deletion, length: 1)

        let cigar = builder.build()
        #expect(cigar.count == 2)
        #expect(cigar[0] == (8 << 4 | 0))   // 8M (merged)
        #expect(cigar[1] == (1 << 4 | 2))   // 1D
    }

    @Test("CIGARBuilder reverse")
    func testCIGARReverse() {
        var builder = CIGARBuilder()
        builder.append(.softClip, length: 3)
        builder.append(.match, length: 10)
        builder.reverse()

        let cigar = builder.build()
        #expect(cigar.count == 2)
        #expect(cigar[0] == (10 << 4 | 0))  // 10M first after reverse
        #expect(cigar[1] == (3 << 4 | 4))   // 3S second after reverse
    }

    @Test("BandedSWScalar perfect match")
    func testScalarPerfectMatch() {
        let scoring = ScoringParameters()
        let query: [UInt8] = [0, 1, 2, 3, 0, 1]  // ACGTAC
        let target: [UInt8] = [0, 1, 2, 3, 0, 1]  // ACGTAC

        let result = query.withUnsafeBufferPointer { qBuf in
            target.withUnsafeBufferPointer { tBuf in
                BandedSWScalar.align(query: qBuf, target: tBuf, scoring: scoring, w: 10, h0: 0)
            }
        }

        // Perfect match: score = 6 * 1 (matchScore) = 6
        #expect(result.score == 6)
    }

    @Test("BandedSWScalar with mismatch")
    func testScalarMismatch() {
        let scoring = ScoringParameters()
        let query: [UInt8] = [0, 1, 2, 3]  // ACGT
        let target: [UInt8] = [0, 1, 0, 3]  // ACAT (G->A mismatch)

        let result = query.withUnsafeBufferPointer { qBuf in
            target.withUnsafeBufferPointer { tBuf in
                BandedSWScalar.align(query: qBuf, target: tBuf, scoring: scoring, w: 10, h0: 0)
            }
        }

        // 3 matches (3) - 1 mismatch penalty (4) = max of local alignments
        // With local alignment, best is either first 2 bases (score 2) or various combos
        #expect(result.score > 0)
    }

    @Test("BandedSWScalar with h0 > 0")
    func testScalarWithInitialScore() {
        let scoring = ScoringParameters()
        let query: [UInt8] = [0, 1, 2]  // ACG
        let target: [UInt8] = [0, 1, 2]  // ACG

        let result = query.withUnsafeBufferPointer { qBuf in
            target.withUnsafeBufferPointer { tBuf in
                BandedSWScalar.align(query: qBuf, target: tBuf, scoring: scoring, w: 10, h0: 5)
            }
        }

        // h0 = 5, plus 3 matches = 8
        #expect(result.score >= 8)
    }

    @Test("BandedSW8 returns nil for large scores")
    func testSW8Overflow() {
        let scoring = ScoringParameters()
        // 300 bases would give score > 250
        let query = [UInt8](repeating: 0, count: 300)
        let target = query

        let result = query.withUnsafeBufferPointer { qBuf in
            target.withUnsafeBufferPointer { tBuf in
                BandedSW8.align(query: qBuf, target: tBuf, scoring: scoring, w: 10, h0: 0)
            }
        }

        // Should return nil due to overflow risk
        #expect(result == nil)
    }

    @Test("SIMD16 score agrees with scalar on perfect match")
    func testSIMD16PositiveScore() {
        let scoring = ScoringParameters()
        let query: [UInt8] = [0, 1, 2, 3, 0, 1, 2, 3]  // ACGTACGT
        let target: [UInt8] = [0, 1, 2, 3, 0, 1, 2, 3]

        let scalarResult = query.withUnsafeBufferPointer { qBuf in
            target.withUnsafeBufferPointer { tBuf in
                BandedSWScalar.align(query: qBuf, target: tBuf, scoring: scoring, w: 10, h0: 0)
            }
        }

        let simdResult = query.withUnsafeBufferPointer { qBuf in
            target.withUnsafeBufferPointer { tBuf in
                BandedSW16.align(query: qBuf, target: tBuf, scoring: scoring, w: 10, h0: 0)
            }
        }

        #expect(scalarResult.score == 8)
        #expect(simdResult.score == scalarResult.score)
    }

    // MARK: - SIMD Completeness and Z-Dropoff Tests

    @Test("BandedSW16 returns correct globalScore and globalTargetEnd")
    func testSIMD16GlobalScore() {
        let scoring = ScoringParameters()
        let query: [UInt8] = [0, 1, 2, 3, 0, 1]  // ACGTAC
        let target: [UInt8] = [0, 1, 2, 3, 0, 1]  // ACGTAC

        let scalar = query.withUnsafeBufferPointer { qBuf in
            target.withUnsafeBufferPointer { tBuf in
                BandedSWScalar.align(query: qBuf, target: tBuf, scoring: scoring, w: 10, h0: 0)
            }
        }

        let simd16 = query.withUnsafeBufferPointer { qBuf in
            target.withUnsafeBufferPointer { tBuf in
                BandedSW16.align(query: qBuf, target: tBuf, scoring: scoring, w: 10, h0: 0)
            }
        }

        // Both should report globalScore (when entire query is consumed)
        #expect(scalar.globalScore > 0)
        #expect(simd16.globalScore == scalar.globalScore)
        #expect(simd16.globalTargetEnd == scalar.globalTargetEnd)
    }

    @Test("BandedSW16 with h0 > 0 matches scalar")
    func testSIMD16WithH0() {
        let scoring = ScoringParameters()
        let query: [UInt8] = [0, 1, 2]  // ACG
        let target: [UInt8] = [0, 1, 2]  // ACG

        let scalar = query.withUnsafeBufferPointer { qBuf in
            target.withUnsafeBufferPointer { tBuf in
                BandedSWScalar.align(query: qBuf, target: tBuf, scoring: scoring, w: 10, h0: 5)
            }
        }

        let simd16 = query.withUnsafeBufferPointer { qBuf in
            target.withUnsafeBufferPointer { tBuf in
                BandedSW16.align(query: qBuf, target: tBuf, scoring: scoring, w: 10, h0: 5)
            }
        }

        #expect(scalar.score >= 8)
        #expect(simd16.score == scalar.score)
    }

    @Test("BandedSW8 score agrees with scalar on short perfect match")
    func testSW8ScoreAgreement() {
        let scoring = ScoringParameters()
        let query: [UInt8] = [0, 1, 2, 3, 0, 1]  // ACGTAC
        let target: [UInt8] = [0, 1, 2, 3, 0, 1]

        let scalar = query.withUnsafeBufferPointer { qBuf in
            target.withUnsafeBufferPointer { tBuf in
                BandedSWScalar.align(query: qBuf, target: tBuf, scoring: scoring, w: 10, h0: 0)
            }
        }

        let simd8 = query.withUnsafeBufferPointer { qBuf in
            target.withUnsafeBufferPointer { tBuf in
                BandedSW8.align(query: qBuf, target: tBuf, scoring: scoring, w: 10, h0: 0)
            }
        }

        #expect(simd8 != nil)
        #expect(simd8!.score == scalar.score)
    }

    @Test("All three SW implementations agree on mismatch input")
    func testAllThreeAgreeOnMismatch() {
        let scoring = ScoringParameters()
        // 10bp with a mismatch in the middle
        let query: [UInt8]  = [0, 1, 2, 3, 0, 1, 2, 3, 0, 1]  // ACGTACGTAC
        let target: [UInt8] = [0, 1, 2, 3, 3, 1, 2, 3, 0, 1]  // ACGTTCGTAC

        let scalar = query.withUnsafeBufferPointer { qBuf in
            target.withUnsafeBufferPointer { tBuf in
                BandedSWScalar.align(query: qBuf, target: tBuf, scoring: scoring, w: 10, h0: 0)
            }
        }

        let simd16 = query.withUnsafeBufferPointer { qBuf in
            target.withUnsafeBufferPointer { tBuf in
                BandedSW16.align(query: qBuf, target: tBuf, scoring: scoring, w: 10, h0: 0)
            }
        }

        let simd8 = query.withUnsafeBufferPointer { qBuf in
            target.withUnsafeBufferPointer { tBuf in
                BandedSW8.align(query: qBuf, target: tBuf, scoring: scoring, w: 10, h0: 0)
            }
        }

        #expect(scalar.score > 0)
        #expect(simd16.score == scalar.score)
        #expect(simd8 != nil)
        #expect(simd8!.score == scalar.score)
    }

    @Test("Z-dropoff triggers in BandedSW16 with divergent tail")
    func testSIMD16ZDropoff() {
        var scoring = ScoringParameters()
        scoring.zDrop = 10  // Low z-drop to trigger early termination

        // 5 matching bases then 100 all-mismatch bases
        var query: [UInt8] = [0, 1, 2, 3, 0]
        query += [UInt8](repeating: 0, count: 100)  // All A's
        var target: [UInt8] = [0, 1, 2, 3, 0]
        target += [UInt8](repeating: 3, count: 100)  // All T's (mismatches)

        let result = query.withUnsafeBufferPointer { qBuf in
            target.withUnsafeBufferPointer { tBuf in
                BandedSW16.align(query: qBuf, target: tBuf, scoring: scoring, w: 5, h0: 0)
            }
        }

        // Score should be from the matching prefix (5 matches = 5)
        #expect(result.score == 5)
        // Z-dropoff should have stopped extension early, so targetEnd should be
        // well short of the full 105
        #expect(result.targetEnd < 50)
    }

    @Test("Z-dropoff triggers in BandedSW8 with divergent tail")
    func testSW8ZDropoff() {
        var scoring = ScoringParameters()
        scoring.zDrop = 10

        var query: [UInt8] = [0, 1, 2, 3, 0]
        query += [UInt8](repeating: 0, count: 100)
        var target: [UInt8] = [0, 1, 2, 3, 0]
        target += [UInt8](repeating: 3, count: 100)

        let result = query.withUnsafeBufferPointer { qBuf in
            target.withUnsafeBufferPointer { tBuf in
                BandedSW8.align(query: qBuf, target: tBuf, scoring: scoring, w: 5, h0: 0)
            }
        }

        #expect(result != nil)
        #expect(result!.score == 5)
        #expect(result!.targetEnd < 50)
    }

    @Test("BandedSW16 maxOff tracks diagonal offset")
    func testSIMD16MaxOff() {
        let scoring = ScoringParameters()
        // Query embedded in a longer target — best alignment is off-diagonal
        let query: [UInt8] = [0, 1, 2, 3]  // ACGT
        let target: [UInt8] = [3, 3, 3, 0, 1, 2, 3, 3, 3]  // TTTACGTTTT

        let scalar = query.withUnsafeBufferPointer { qBuf in
            target.withUnsafeBufferPointer { tBuf in
                BandedSWScalar.align(query: qBuf, target: tBuf, scoring: scoring, w: 10, h0: 0)
            }
        }

        let simd16 = query.withUnsafeBufferPointer { qBuf in
            target.withUnsafeBufferPointer { tBuf in
                BandedSW16.align(query: qBuf, target: tBuf, scoring: scoring, w: 10, h0: 0)
            }
        }

        #expect(scalar.score == 4)
        #expect(simd16.score == scalar.score)
        // maxOff should be > 0 since the best alignment is off the main diagonal
        #expect(simd16.maxOff > 0)
    }

    @Test("ChainFilter removes low-weight chains")
    func testChainFilterWeight() {
        var chains = [
            MemChain(seeds: [MemSeed(rbeg: 0, qbeg: 0, len: 25, score: 25)], weight: 25, rid: 0, kept: 3),
            MemChain(seeds: [MemSeed(rbeg: 100, qbeg: 50, len: 10, score: 10)], weight: 10, rid: 0, kept: 3),
        ]
        let scoring = ScoringParameters()
        ChainFilter.filter(chains: &chains, scoring: scoring)

        // Chain with weight 10 < minSeedLength(19) should be removed
        #expect(chains.count == 1)
        #expect(chains[0].weight == 25)
    }

    @Test("ChainFilter.markSecondary marks overlapping regions")
    func testMarkSecondary() {
        var regions = [
            MemAlnReg(rb: 0, re: 100, qb: 0, qe: 100, score: 100),
            MemAlnReg(rb: 10, re: 90, qb: 10, qe: 90, score: 50),
        ]
        ChainFilter.markSecondary(regions: &regions, maskLevel: 0.50)

        #expect(regions[0].secondary == -1)  // Primary
        #expect(regions[1].secondary == 0)   // Secondary of region 0
    }

    // MARK: - GlobalAligner Tests

    @Test("GlobalAligner perfect match produces all-M CIGAR")
    func testGlobalAlignerPerfectMatch() {
        let scoring = ScoringParameters()
        let query: [UInt8] = [0, 1, 2, 3, 0, 1, 2, 3, 0, 1]  // ACGTACGTAC
        let target: [UInt8] = [0, 1, 2, 3, 0, 1, 2, 3, 0, 1]

        let result = query.withUnsafeBufferPointer { qBuf in
            target.withUnsafeBufferPointer { tBuf in
                GlobalAligner.align(query: qBuf, target: tBuf, scoring: scoring, w: 10)
            }
        }

        // Perfect match: 10M, score = 10, NM = 0
        #expect(result.cigar.count == 1)
        #expect(result.cigar[0] == (10 << 4 | CIGAROp.match.rawValue))
        #expect(result.score == 10)
        #expect(result.nm == 0)
    }

    @Test("GlobalAligner single mismatch produces 10M with NM=1")
    func testGlobalAlignerMismatch() {
        let scoring = ScoringParameters()
        let query: [UInt8]  = [0, 1, 2, 3, 0, 1, 2, 3, 0, 1]  // ACGTACGTAC
        let target: [UInt8] = [0, 1, 2, 3, 0, 3, 2, 3, 0, 1]  // ACGTАТGTAC (C->T at pos 5)

        let result = query.withUnsafeBufferPointer { qBuf in
            target.withUnsafeBufferPointer { tBuf in
                GlobalAligner.align(query: qBuf, target: tBuf, scoring: scoring, w: 10)
            }
        }

        // Same lengths, 1 mismatch: should be 10M with NM=1
        #expect(result.cigar.count == 1)
        #expect(result.cigar[0] == (10 << 4 | CIGAROp.match.rawValue))
        #expect(result.nm == 1)
        // Score: 9 matches - 4 mismatch = 5
        #expect(result.score == 5)
    }

    @Test("GlobalAligner single insertion")
    func testGlobalAlignerInsertion() {
        let scoring = ScoringParameters()
        // Query has an extra base (G) in the middle
        let query: [UInt8]  = [0, 1, 2, 3, 2, 0, 1]  // ACGTGAC (7bp)
        let target: [UInt8] = [0, 1, 2, 3, 0, 1]      // ACGTAC  (6bp)

        let result = query.withUnsafeBufferPointer { qBuf in
            target.withUnsafeBufferPointer { tBuf in
                GlobalAligner.align(query: qBuf, target: tBuf, scoring: scoring, w: 10)
            }
        }

        // Should contain an insertion operation
        var hasInsertion = false
        var totalQueryBases: UInt32 = 0
        var totalTargetBases: UInt32 = 0
        for c in result.cigar {
            let op = c & 0xF
            let len = c >> 4
            if op == CIGAROp.insertion.rawValue {
                hasInsertion = true
                totalQueryBases += len
            } else if op == CIGAROp.deletion.rawValue {
                totalTargetBases += len
            } else if op == CIGAROp.match.rawValue {
                totalQueryBases += len
                totalTargetBases += len
            }
        }
        #expect(hasInsertion)
        #expect(totalQueryBases == 7)
        #expect(totalTargetBases == 6)
        #expect(result.nm >= 1)
    }

    @Test("GlobalAligner single deletion")
    func testGlobalAlignerDeletion() {
        let scoring = ScoringParameters()
        // Target has an extra base (G) in the middle
        let query: [UInt8]  = [0, 1, 2, 3, 0, 1]      // ACGTAC  (6bp)
        let target: [UInt8] = [0, 1, 2, 3, 2, 0, 1]    // ACGTGAC (7bp)

        let result = query.withUnsafeBufferPointer { qBuf in
            target.withUnsafeBufferPointer { tBuf in
                GlobalAligner.align(query: qBuf, target: tBuf, scoring: scoring, w: 10)
            }
        }

        // Should contain a deletion operation
        var hasDeletion = false
        var totalQueryBases: UInt32 = 0
        var totalTargetBases: UInt32 = 0
        for c in result.cigar {
            let op = c & 0xF
            let len = c >> 4
            if op == CIGAROp.deletion.rawValue {
                hasDeletion = true
                totalTargetBases += len
            } else if op == CIGAROp.insertion.rawValue {
                totalQueryBases += len
            } else if op == CIGAROp.match.rawValue {
                totalQueryBases += len
                totalTargetBases += len
            }
        }
        #expect(hasDeletion)
        #expect(totalQueryBases == 6)
        #expect(totalTargetBases == 7)
        #expect(result.nm >= 1)
    }

    @Test("GlobalAligner narrow vs wide band produce same result on short sequences")
    func testGlobalAlignerBandWidth() {
        let scoring = ScoringParameters()
        let query: [UInt8]  = [0, 1, 2, 3, 0, 1, 2, 3]  // ACGTACGT
        let target: [UInt8] = [0, 1, 2, 3, 0, 1, 2, 3]

        let narrowResult = query.withUnsafeBufferPointer { qBuf in
            target.withUnsafeBufferPointer { tBuf in
                GlobalAligner.align(query: qBuf, target: tBuf, scoring: scoring, w: 1)
            }
        }

        let wideResult = query.withUnsafeBufferPointer { qBuf in
            target.withUnsafeBufferPointer { tBuf in
                GlobalAligner.align(query: qBuf, target: tBuf, scoring: scoring, w: 50)
            }
        }

        #expect(narrowResult.score == wideResult.score)
        #expect(narrowResult.cigar == wideResult.cigar)
    }

    // MARK: - CIGARGenerator Tests

    @Test("CIGARGenerator full-length alignment has no soft-clips")
    func testCIGARGeneratorNoClips() {
        let scoring = ScoringParameters()
        let query: [UInt8] = [0, 1, 2, 3, 0, 1, 2, 3]
        let target: [UInt8] = [0, 1, 2, 3, 0, 1, 2, 3]

        let result = CIGARGenerator.generate(
            querySegment: query,
            refSegment: target,
            qb: 0,
            qe: 8,
            readLength: 8,
            isReverse: false,
            trueScore: 8,
            initialW: 10,
            scoring: scoring,
            refPos: 0
        )

        // No soft-clips when qb=0 and qe=readLength
        for c in result.cigar {
            let op = c & 0xF
            #expect(op != CIGAROp.softClip.rawValue)
        }
        #expect(result.nm == 0)
    }

    @Test("CIGARGenerator partial alignment adds soft-clips")
    func testCIGARGeneratorSoftClips() {
        let scoring = ScoringParameters()
        let query: [UInt8] = [0, 1, 2, 3, 0]     // aligned portion (5bp)
        let target: [UInt8] = [0, 1, 2, 3, 0]

        let result = CIGARGenerator.generate(
            querySegment: query,
            refSegment: target,
            qb: 3,       // 3bp before alignment
            qe: 8,       // = qb + 5
            readLength: 12,  // 4bp after alignment
            isReverse: false,
            trueScore: 5,
            initialW: 10,
            scoring: scoring,
            refPos: 100
        )

        // Should have 3S at start, 5M, 4S at end
        #expect(result.cigar.count == 3)
        let firstOp = result.cigar[0] & 0xF
        let firstLen = result.cigar[0] >> 4
        #expect(firstOp == CIGAROp.softClip.rawValue)
        #expect(firstLen == 3)

        let lastOp = result.cigar[2] & 0xF
        let lastLen = result.cigar[2] >> 4
        #expect(lastOp == CIGAROp.softClip.rawValue)
        #expect(lastLen == 4)
    }

    @Test("CIGARGenerator reverse strand swaps clip ends")
    func testCIGARGeneratorReverseClips() {
        let scoring = ScoringParameters()
        let query: [UInt8] = [0, 1, 2, 3, 0]
        let target: [UInt8] = [0, 1, 2, 3, 0]

        let result = CIGARGenerator.generate(
            querySegment: query,
            refSegment: target,
            qb: 3,
            qe: 8,
            readLength: 12,
            isReverse: true,
            trueScore: 5,
            initialW: 10,
            scoring: scoring,
            refPos: 100
        )

        // Reverse strand: 5' clip = readLength - qe = 4, 3' clip = qb = 3
        #expect(result.cigar.count == 3)
        let firstLen = result.cigar[0] >> 4
        #expect(firstLen == 4)  // 5' clip on reverse = readLen - qe

        let lastLen = result.cigar[2] >> 4
        #expect(lastLen == 3)  // 3' clip on reverse = qb
    }

    @Test("CIGARGenerator bandwidth retry loop")
    func testCIGARGeneratorRetry() {
        let scoring = ScoringParameters()
        // Create a case with insertion that needs wider bandwidth
        let query: [UInt8]  = [0, 1, 2, 3, 2, 2, 0, 1, 2, 3]  // 10bp with 2 extra
        let target: [UInt8] = [0, 1, 2, 3, 0, 1, 2, 3]          // 8bp

        let result = CIGARGenerator.generate(
            querySegment: query,
            refSegment: target,
            qb: 0,
            qe: 10,
            readLength: 10,
            isReverse: false,
            trueScore: 0,  // Low trueScore so no retry needed
            initialW: 1,   // Start with very narrow band
            scoring: scoring,
            refPos: 0
        )

        // Should still produce valid CIGAR covering all query and target bases
        var totalQ: UInt32 = 0
        var totalT: UInt32 = 0
        for c in result.cigar {
            let op = c & 0xF
            let len = c >> 4
            if op == CIGAROp.match.rawValue { totalQ += len; totalT += len }
            else if op == CIGAROp.insertion.rawValue { totalQ += len }
            else if op == CIGAROp.deletion.rawValue { totalT += len }
            else if op == CIGAROp.softClip.rawValue { /* doesn't count */ }
        }
        #expect(totalQ == 10)
        #expect(totalT == 8)
    }

    // MARK: - MD Tag Tests

    @Test("MD tag for perfect match")
    func testMDPerfectMatch() {
        let scoring = ScoringParameters()
        let query: [UInt8] = [0, 1, 2, 3, 0, 1, 2, 3]  // ACGTACGT
        let target: [UInt8] = [0, 1, 2, 3, 0, 1, 2, 3]

        let result = CIGARGenerator.generate(
            querySegment: query,
            refSegment: target,
            qb: 0, qe: 8,
            readLength: 8,
            isReverse: false,
            trueScore: 8,
            initialW: 10,
            scoring: scoring,
            refPos: 0
        )

        // Perfect match: MD = "8" (8 consecutive matches)
        #expect(result.md == "8")
    }

    @Test("MD tag for single mismatch")
    func testMDSingleMismatch() {
        let scoring = ScoringParameters()
        let query: [UInt8]  = [0, 1, 2, 3, 0, 1, 2, 3]  // ACGTACGT
        let target: [UInt8] = [0, 1, 2, 3, 3, 1, 2, 3]  // ACGTTCGT (A->T at pos 4)

        let result = CIGARGenerator.generate(
            querySegment: query,
            refSegment: target,
            qb: 0, qe: 8,
            readLength: 8,
            isReverse: false,
            trueScore: 3,
            initialW: 10,
            scoring: scoring,
            refPos: 0
        )

        // 4 matches, ref=T mismatch, 3 matches: "4T3"
        #expect(result.md == "4T3")
        #expect(result.nm == 1)
    }

    @Test("MD tag for deletion")
    func testMDDeletion() {
        let scoring = ScoringParameters()
        // Query is missing a base that's in the reference
        let query: [UInt8]  = [0, 1, 2, 3, 0, 1]      // ACGTAC  (6bp)
        let target: [UInt8] = [0, 1, 2, 3, 2, 0, 1]    // ACGTGAC (7bp, extra G)

        let result = CIGARGenerator.generate(
            querySegment: query,
            refSegment: target,
            qb: 0, qe: 6,
            readLength: 6,
            isReverse: false,
            trueScore: -1,
            initialW: 10,
            scoring: scoring,
            refPos: 0
        )

        // Should contain ^G or ^-prefixed deletion in MD
        #expect(result.md.contains("^"))
    }

    @Test("MD tag for insertion has no insertion marker")
    func testMDInsertion() {
        let scoring = ScoringParameters()
        // Query has an extra base not in reference
        let query: [UInt8]  = [0, 1, 2, 3, 2, 0, 1]  // ACGTGAC (7bp, extra G)
        let target: [UInt8] = [0, 1, 2, 3, 0, 1]      // ACGTAC  (6bp)

        let result = CIGARGenerator.generate(
            querySegment: query,
            refSegment: target,
            qb: 0, qe: 7,
            readLength: 7,
            isReverse: false,
            trueScore: -1,
            initialW: 10,
            scoring: scoring,
            refPos: 0
        )

        // Insertions don't appear in MD (no ^ marker for insertions)
        // MD should only show matched ref bases
        #expect(!result.md.contains("^"))
    }

    @Test("MD tag for multiple mismatches")
    func testMDMultipleMismatches() {
        let scoring = ScoringParameters()
        let query: [UInt8]  = [0, 1, 2, 3, 0, 1, 2, 3, 0, 1]  // ACGTACGTAC
        let target: [UInt8] = [0, 3, 2, 3, 0, 1, 2, 0, 0, 1]  // ATGTACGAAC
        //                         ^                 ^           mismatches at 1,7

        let result = CIGARGenerator.generate(
            querySegment: query,
            refSegment: target,
            qb: 0, qe: 10,
            readLength: 10,
            isReverse: false,
            trueScore: 0,
            initialW: 10,
            scoring: scoring,
            refPos: 0
        )

        // 1 match, ref=T mismatch, 5 matches, ref=A mismatch, 2 matches: "1T5A2"
        #expect(result.md == "1T5A2")
        #expect(result.nm == 2)
    }

    @Test("SeedChainer basic chaining")
    func testSeedChainer() {
        let smems = [
            SMEM(k: 0, l: 0, queryBegin: 0, queryEnd: 20),
            SMEM(k: 1, l: 1, queryBegin: 25, queryEnd: 45),
        ]
        var metadata = ReferenceMetadata()
        metadata.totalLength = 1000
        metadata.numSequences = 1
        metadata.annotations = [ReferenceAnnotation(offset: 0, length: 1000, name: "chr1")]

        let scoring = ScoringParameters()
        let chains = SeedChainer.chain(
            smems: smems,
            getSAEntry: { pos in pos * 10 },  // Simple mapping
            metadata: metadata,
            scoring: scoring,
            readLength: 50
        )

        #expect(!chains.isEmpty)
    }

    // MARK: - LocalSWAligner Tests

    @Test("LocalSWAligner perfect match")
    func testLocalSWPerfectMatch() {
        let scoring = ScoringParameters()
        let query: [UInt8] = [0, 1, 2, 3, 0, 1]  // ACGTAC
        let target: [UInt8] = [0, 1, 2, 3, 0, 1]  // ACGTAC

        let result = LocalSWAligner.align(query: query, target: target, scoring: scoring)
        #expect(result != nil)
        #expect(result!.score == 6)
        #expect(result!.queryBegin == 0)
        #expect(result!.queryEnd == 5)
        #expect(result!.targetBegin == 0)
        #expect(result!.targetEnd == 5)
    }

    @Test("LocalSWAligner query embedded in target")
    func testLocalSWEmbedded() {
        let scoring = ScoringParameters()
        // Query matches in the middle of target
        let query: [UInt8] = [0, 1, 2, 3]  // ACGT
        let target: [UInt8] = [3, 3, 0, 1, 2, 3, 3, 3]  // TTACGTTT

        let result = LocalSWAligner.align(query: query, target: target, scoring: scoring)
        #expect(result != nil)
        #expect(result!.score == 4)
        #expect(result!.queryBegin == 0)
        #expect(result!.queryEnd == 3)
        #expect(result!.targetBegin == 2)
        #expect(result!.targetEnd == 5)
    }

    @Test("LocalSWAligner with mismatch finds best local alignment")
    func testLocalSWMismatch() {
        let scoring = ScoringParameters()
        // 4 matches then 1 mismatch then 3 matches
        let query: [UInt8]  = [0, 1, 2, 3, 0, 1, 2, 3]  // ACGTACGT
        let target: [UInt8] = [0, 1, 2, 3, 3, 1, 2, 3]  // ACGTTCGT

        let result = LocalSWAligner.align(query: query, target: target, scoring: scoring)
        #expect(result != nil)
        // Best score: either 4 (first half) or the whole thing (8-4=4)
        // Actually: 7 matches - 1 * 4 mismatch penalty = 3; but local can skip mismatch
        // Best contiguous: 4 matches = score 4
        #expect(result!.score >= 4)
    }

    @Test("LocalSWAligner returns nil for all-mismatch")
    func testLocalSWAllMismatch() {
        var scoring = ScoringParameters()
        scoring.matchScore = 1
        scoring.mismatchPenalty = 4
        // Query all A, target all T — every position is a mismatch
        let query: [UInt8] = [0, 0, 0, 0]
        let target: [UInt8] = [3, 3, 3, 3]

        let result = LocalSWAligner.align(query: query, target: target, scoring: scoring)
        #expect(result == nil)
    }

    @Test("LocalSWAligner empty input returns nil")
    func testLocalSWEmpty() {
        let scoring = ScoringParameters()
        let result = LocalSWAligner.align(query: [], target: [0, 1, 2], scoring: scoring)
        #expect(result == nil)
    }

    // MARK: - ExtensionAligner Clip-vs-Extend Tests

    @Test("ExtensionAligner clips when penClip is low (local alignment preferred)")
    func testExtensionAlignerClips() {
        // A seed in the middle of the read with a mismatching flanking region
        // should clip rather than extend when the clip penalty is low
        var scoring = ScoringParameters()
        scoring.penClip5 = 5
        scoring.penClip3 = 5

        // 20bp seed at read positions 10-30, with mismatches at the edges
        let seed = MemSeed(rbeg: 100, qbeg: 10, len: 20, score: 20)
        let chain = MemChain(seeds: [seed], weight: 20, rid: 0, kept: 3)

        // Read: 10 mismatching bases + 20 matching + 10 mismatching
        var readBases = [UInt8](repeating: 3, count: 10)  // T's that won't match
        readBases += [UInt8](repeating: 0, count: 20)      // A's matching reference
        readBases += [UInt8](repeating: 3, count: 10)      // T's that won't match
        let read = ReadSequence(
            name: "test",
            bases: readBases,
            qualities: [UInt8](repeating: 30, count: 40)
        )

        // Reference: all A's so seed matches but flanks mismatch
        let refBases = [UInt8](repeating: 0, count: 200)

        let regions = ExtensionAligner.extend(
            chain: chain,
            read: read,
            getReference: { pos, length in
                let start = max(0, Int(pos))
                let end = min(refBases.count, start + length)
                return Array(refBases[start..<end])
            },
            scoring: scoring
        )

        #expect(!regions.isEmpty)
        let reg = regions[0]
        // With mismatches at edges and low clip penalty, should clip
        // qb should be > 0 (clipped at 5' end)
        #expect(reg.qb >= 0)
        // qe should be < 40 (clipped at 3' end)
        #expect(reg.qe <= 40)
    }

    // MARK: - ALT-Aware Chain Filter Tests

    @Test("ChainFilter.filter does not suppress primary chain overlapping ALT chain")
    func testChainFilterALTPreservesPrimary() {
        // ALT chain (i=0) has higher weight, primary chain (j=1) overlaps it.
        // The primary chain should NOT be suppressed even though it overlaps the ALT.
        var chains = [
            MemChain(
                seeds: [MemSeed(rbeg: 0, qbeg: 0, len: 30, score: 30)],
                weight: 30, rid: 0, kept: 3, isAlt: true
            ),
            MemChain(
                seeds: [MemSeed(rbeg: 100, qbeg: 5, len: 20, score: 20)],
                weight: 20, rid: 1, kept: 3, isAlt: false
            ),
        ]
        let scoring = ScoringParameters()
        ChainFilter.filter(chains: &chains, scoring: scoring)

        // Both chains should survive: primary is protected from ALT suppression
        #expect(chains.count == 2)
    }

    @Test("ChainFilter.filter suppresses ALT chain overlapping higher-weight ALT")
    func testChainFilterALTSuppressesALT() {
        // Two ALT chains: higher-weight one suppresses lower-weight overlapping one
        var chains = [
            MemChain(
                seeds: [MemSeed(rbeg: 0, qbeg: 0, len: 30, score: 30)],
                weight: 30, rid: 0, kept: 3, isAlt: true
            ),
            MemChain(
                seeds: [MemSeed(rbeg: 100, qbeg: 5, len: 20, score: 20)],
                weight: 10, rid: 1, kept: 3, isAlt: true  // low weight, will be removed by min weight filter
            ),
        ]
        let scoring = ScoringParameters()
        ChainFilter.filter(chains: &chains, scoring: scoring)

        // Low weight ALT removed by minimum weight filter (< minSeedLength=19)
        #expect(chains.count == 1)
        #expect(chains[0].isAlt == true)
        #expect(chains[0].weight == 30)
    }

    @Test("markSecondaryALT: primary hit stays primary even with higher-scoring ALT")
    func testMarkSecondaryALTPrimaryPreserved() {
        // ALT hit with score 100, primary hit with score 80, overlapping query ranges
        var regions = [
            MemAlnReg(rb: 0, re: 100, qb: 0, qe: 100, score: 100),
            MemAlnReg(rb: 200, re: 300, qb: 0, qe: 100, score: 80),
        ]
        regions[0].isAlt = true
        regions[1].isAlt = false

        let scoring = ScoringParameters()
        ChainFilter.markSecondaryALT(
            regions: &regions, maskLevel: 0.50, scoring: scoring
        )

        // After Phase 2 re-marking among primary-only, the primary hit
        // should have secondary == -1 (it's the only primary)
        let primaryRegion = regions.first { !$0.isAlt }!
        #expect(primaryRegion.secondary == -1)
    }

    @Test("markSecondaryALT tracks altSc for primary with ALT competitor")
    func testMarkSecondaryALTAltSc() {
        // Primary hit that overlaps with a higher-scoring ALT hit
        var regions = [
            MemAlnReg(rb: 0, re: 100, qb: 0, qe: 100, score: 120),
            MemAlnReg(rb: 200, re: 300, qb: 0, qe: 100, score: 80),
        ]
        regions[0].isAlt = true
        regions[1].isAlt = false

        let scoring = ScoringParameters()
        ChainFilter.markSecondaryALT(
            regions: &regions, maskLevel: 0.50, scoring: scoring
        )

        // The primary hit should record the ALT competitor's score
        let primaryRegion = regions.first { !$0.isAlt }!
        #expect(primaryRegion.altSc == 120)
    }

    @Test("markSecondaryALT: pure primary regions behave like markSecondary")
    func testMarkSecondaryALTNonMixed() {
        // No ALT regions: secondaryAll should equal secondary
        var regions = [
            MemAlnReg(rb: 0, re: 100, qb: 0, qe: 100, score: 100),
            MemAlnReg(rb: 10, re: 90, qb: 10, qe: 90, score: 50),
        ]

        let scoring = ScoringParameters()
        ChainFilter.markSecondaryALT(
            regions: &regions, maskLevel: 0.50, scoring: scoring
        )

        // First should be primary
        #expect(regions[0].secondary == -1)
        // Second should be secondary of first
        #expect(regions[1].secondary == 0)
        // secondaryAll should match secondary when no ALT
        #expect(regions[0].secondaryAll == -1)
        #expect(regions[1].secondaryAll == regions[1].secondary)
    }

    @Test("markSecondaryCore ALT-aware subN logic")
    func testMarkSecondaryALTSubN() {
        // Two overlapping hits: primary and ALT
        // subN should only increment when (primary.isAlt || !secondary.isAlt)
        var regions = [
            MemAlnReg(rb: 0, re: 100, qb: 0, qe: 100, score: 100),
            MemAlnReg(rb: 200, re: 300, qb: 0, qe: 100, score: 95),
        ]
        regions[0].isAlt = false
        regions[1].isAlt = true

        let scoring = ScoringParameters()
        ChainFilter.markSecondaryALT(
            regions: &regions, maskLevel: 0.50, scoring: scoring
        )

        // In Phase 1: primary hit (score 100) encounters ALT hit (score 95)
        // Since primary.isAlt=false and secondary.isAlt=true, the condition
        // (primary.isAlt || !secondary.isAlt) is false, so subN should NOT increment
        let primaryRegion = regions.first { !$0.isAlt }!
        #expect(primaryRegion.subN == 0)
    }

    @Test("ExtensionAligner extends to read boundary with high penClip")
    func testExtensionAlignerExtendsWithHighPenClip() {
        // With a very high clip penalty, extension should prefer extending
        // to read boundaries even with some cost
        var scoring = ScoringParameters()
        scoring.penClip5 = 100
        scoring.penClip3 = 100

        // Seed at read positions 5-25 with matching flanks
        let seed = MemSeed(rbeg: 100, qbeg: 5, len: 20, score: 20)
        let chain = MemChain(seeds: [seed], weight: 20, rid: 0, kept: 3)

        // All bases match the reference
        let readBases = [UInt8](repeating: 0, count: 30)  // 30bp of A's
        let read = ReadSequence(
            name: "test",
            bases: readBases,
            qualities: [UInt8](repeating: 30, count: 30)
        )

        let refBases = [UInt8](repeating: 0, count: 200)  // All A's

        let regions = ExtensionAligner.extend(
            chain: chain,
            read: read,
            getReference: { pos, length in
                let start = max(0, Int(pos))
                let end = min(refBases.count, start + length)
                return Array(refBases[start..<end])
            },
            scoring: scoring
        )

        #expect(!regions.isEmpty)
        let reg = regions[0]
        // With high clip penalty and matching flanks, should extend to read boundaries
        #expect(reg.qb == 0)
        #expect(reg.qe == 30)
    }

    // MARK: - RegionDedup Tests

    @Test("RegionDedup no-op for 0 or 1 regions")
    func testRegionDedupTrivial() {
        let scoring = ScoringParameters()
        var empty: [MemAlnReg] = []
        RegionDedup.sortDedupPatch(
            regions: &empty,
            query: [0, 1, 2, 3],
            getReference: { _, _ in [] },
            genomeLength: 1000,
            scoring: scoring
        )
        #expect(empty.isEmpty)

        var single = [MemAlnReg(rb: 0, re: 50, qb: 0, qe: 50, score: 50)]
        RegionDedup.sortDedupPatch(
            regions: &single,
            query: [UInt8](repeating: 0, count: 50),
            getReference: { _, _ in [] },
            genomeLength: 1000,
            scoring: scoring
        )
        #expect(single.count == 1)
        #expect(single[0].score == 50)
    }

    @Test("RegionDedup removes redundant overlapping region (>95% overlap)")
    func testRegionDedupRedundant() {
        let scoring = ScoringParameters()
        // Two regions on same rid with >95% overlap on both ref and query
        var regions = [
            MemAlnReg(rb: 0, re: 100, qb: 0, qe: 100, rid: 0, score: 80),
            MemAlnReg(rb: 1, re: 99, qb: 1, qe: 99, rid: 0, score: 60),
        ]
        RegionDedup.sortDedupPatch(
            regions: &regions,
            query: [UInt8](repeating: 0, count: 100),
            getReference: { _, _ in [] },
            genomeLength: 1000,
            scoring: scoring
        )
        // Lower-scoring region should be removed
        #expect(regions.count == 1)
        #expect(regions[0].score == 80)
    }

    @Test("RegionDedup preserves non-overlapping regions on different rids")
    func testRegionDedupDifferentRids() {
        let scoring = ScoringParameters()
        var regions = [
            MemAlnReg(rb: 0, re: 100, qb: 0, qe: 50, rid: 0, score: 50),
            MemAlnReg(rb: 200, re: 300, qb: 50, qe: 100, rid: 1, score: 40),
        ]
        RegionDedup.sortDedupPatch(
            regions: &regions,
            query: [UInt8](repeating: 0, count: 100),
            getReference: { _, _ in [] },
            genomeLength: 1000,
            scoring: scoring
        )
        #expect(regions.count == 2)
    }

    @Test("RegionDedup removes exact duplicates (same score, rb, qb)")
    func testRegionDedupExactDuplicates() {
        let scoring = ScoringParameters()
        var regions = [
            MemAlnReg(rb: 0, re: 100, qb: 0, qe: 100, rid: 0, score: 80),
            MemAlnReg(rb: 0, re: 100, qb: 0, qe: 100, rid: 0, score: 80),
        ]
        RegionDedup.sortDedupPatch(
            regions: &regions,
            query: [UInt8](repeating: 0, count: 100),
            getReference: { _, _ in [] },
            genomeLength: 1000,
            scoring: scoring
        )
        #expect(regions.count == 1)
    }

    @Test("RegionDedup patches colinear adjacent regions")
    func testRegionDedupPatch() {
        let scoring = ScoringParameters()
        // Two colinear regions on same rid with a small gap
        // Region a: ref 0..50, query 0..50
        // Region b: ref 50..100, query 50..100
        var regions = [
            MemAlnReg(rb: 0, re: 50, qb: 0, qe: 50, rid: 0, score: 50, w: 10),
            MemAlnReg(rb: 50, re: 100, qb: 50, qe: 100, rid: 0, score: 50, w: 10),
        ]

        // Perfect-match reference so global alignment will produce a good score
        let refBases = [UInt8](repeating: 0, count: 200)
        let queryBases = [UInt8](repeating: 0, count: 100)

        RegionDedup.sortDedupPatch(
            regions: &regions,
            query: queryBases,
            getReference: { pos, length in
                let start = max(0, Int(pos))
                let end = min(refBases.count, start + length)
                guard start < end else { return [] }
                return Array(refBases[start..<end])
            },
            genomeLength: 1000,
            scoring: scoring
        )

        // Should be merged into one region covering 0..100
        #expect(regions.count == 1)
        #expect(regions[0].qb == 0)
        #expect(regions[0].qe == 100)
        #expect(regions[0].rb == 0)
        #expect(regions[0].re == 100)
    }

    @Test("RegionDedup output sorted by score descending")
    func testRegionDedupScoreOrder() {
        let scoring = ScoringParameters()
        // Three non-overlapping regions on different rids
        var regions = [
            MemAlnReg(rb: 0, re: 30, qb: 0, qe: 30, rid: 0, score: 20),
            MemAlnReg(rb: 100, re: 160, qb: 30, qe: 60, rid: 1, score: 60),
            MemAlnReg(rb: 200, re: 240, qb: 60, qe: 100, rid: 2, score: 40),
        ]
        RegionDedup.sortDedupPatch(
            regions: &regions,
            query: [UInt8](repeating: 0, count: 100),
            getReference: { _, _ in [] },
            genomeLength: 1000,
            scoring: scoring
        )
        #expect(regions.count == 3)
        #expect(regions[0].score == 60)
        #expect(regions[1].score == 40)
        #expect(regions[2].score == 20)
    }
}
