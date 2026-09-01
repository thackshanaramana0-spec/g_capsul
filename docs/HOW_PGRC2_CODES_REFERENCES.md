# How PgRC2 codes references, and what it means for us

Read from their source, not inferred. `matching/SimplePgMatcher.cpp`.

## The mechanism

Encoding, lines 104-133:

    for each match:
        copy the literal gap [pos, match.posDestText) into destPg
        destPg[nPos++] = MATCH_MARK                    // one '%' INSIDE the sequence
        writeValue<uint32_t>(pgMapOffDest, match.posSrcText)
        writeUIntByteFrugal(pgMapLenDest, match.length - minMatchLength)
        pos = match.endPosDestText()

Decoding, lines 327-340:

    while ((markPos = destPg.find(MATCH_MARK, posDest)) != npos) {
        resPg.append(destPg, posDest, markPos - posDest);   // literal up to the mark
        posDest = markPos + 1;
        read offset, read length, append srcPg[off, off+len)
    }

**The destination is never stored.** It is *where the mark sits* in the literal
stream: the decoder is already at that position when it reads the mark. They
persist only offsets and lengths. We persist destinations explicitly as
mem_dstgap.

Their streams, measured by instrumenting writeCompressedCollectiveParallel
(CompressionJob carries a label), on S. acidocaldarius:

    Good sequence mapping - offsets   raw 14,136   coded  9,504
    lengths                           raw  5,727   coded  4,744
    Bad sequence mapping - offsets    raw  2,428   coded  1,629
    lengths                           raw    608   coded    571

Offsets are 4 B each, so raw 16,564 B = **4,141 matches**.

## Per match, like for like

| | ours | PgRC2 |
|---|---|---|
| source | 21.8 bits | 21.5 bits |
| length | 8.8 bits | 10.3 bits |
| destination | 4.5 bits | **0.0 bits** |
| total | 35.2 bits | 31.8 bits |
| matches | **21,000** | **4,141** |

Their entire per-match advantage is the destination, and it is exactly the
4.5 bits the MATCH_MARK removes. Our length coding is already cheaper than
theirs and our source coding is within 1.5%.

## Adopting the MATCH_MARK would make us WORSE

Placing 21,000 marks among 2,500,501 literal symbols costs at minimum
log2(C(L+N,N)) = **21,903 B**. We pay **11,929 B** for the same information --
**0.54x the bound** -- because our destination gaps are clustered, not uniform,
and a delta-coded gap stream exploits that while a positional marker cannot.

Their marks are also not free. The cost is absorbed into the 595,838 B sequence
stream as a fifth symbol and cannot be read off separately, so "0.0 bits" is
where their accounting puts it, not what it costs.

**Conclusion: the in-band match marker is the wrong thing to copy.** It is
elegant and it is cheaper *for them*, but our explicit stream is already below
what the technique can achieve.

## The real difference is match COUNT, and ours is already optimal

21,000 against 4,141, a factor of 5.07. Swept MINMEM on S. acidocaldarius with
the whole archive as the objective:

| MINMEM | archive | matches | vs PgRC2 |
|---|---|---|---|
| 20 | 3,145,674 | 21,215 | -0.95% |
| **24** | **3,140,692** | **21,000** | **-0.79%** |
| 28 | 3,143,181 | 19,871 | -0.87% |
| 32 | 3,145,672 | 18,928 | -0.95% |
| 45 (their p45) | 3,155,833 | 16,841 | -1.28% |
| 64 | 3,181,112 | 13,814 | -2.09% |
| 90 | 3,220,352 | 10,973 | -3.35% |

24 is the optimum, and it is a genuine interior optimum -- worse on both sides.
Adopting their minimum match length of 45 costs us 15,141 B. Our extra matches
pay for themselves; they are not waste.

Their 4,141 matches come from a differently-shaped pseudogenome (the hq/lq/n
three-way split, each self-matched separately), not from a better reference
encoding. That is a construction difference, not a coding one.

## What this closes

The reference gap was the largest single line in the stream comparison
(+68,547 B) and it looked like the obvious target. It is not:

- Their coding technique would cost us more than what we already do.
- Our per-match cost is within 11% of theirs and cheaper on lengths.
- Our match count is at a verified interior optimum.

The remaining honest statement is that the two tools reach a similar total by
different decompositions: we spend more on references and far less on mismatch
symbols (124,280 against 208,234, -40%). On this dataset their decomposition
wins by 24,662 B, 0.79%.
