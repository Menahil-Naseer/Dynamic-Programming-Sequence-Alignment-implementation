
import random
import time


# 1. DYNAMIC PROGRAMMING (Needleman-Wunsch)

def build_dp_table(seq1: str, seq2: str,
                   match_score: int = 1,
                   mismatch_penalty: int = -1,
                   gap_penalty: int = -2) -> list[list[int]]:
    


    m, n = len(seq1), len(seq2)
    dp = [[0] * (n + 1) for _ in range(m + 1)]

    # ── base cases ──────────────────────────────────────────
    for i in range(m + 1):
        dp[i][0] = i * gap_penalty
    for j in range(n + 1):
        dp[0][j] = j * gap_penalty


    for i in range(1, m + 1):
        for j in range(1, n + 1):
            diag = match_score if seq1[i-1] == seq2[j-1] else mismatch_penalty
            dp[i][j] = max(
                dp[i-1][j-1] + diag,          # match / mismatch
                dp[i-1][j]   + gap_penalty,    # deletion
                dp[i][j-1]   + gap_penalty,    # insertion
            )

    return dp


def traceback(dp: list[list[int]],
              seq1: str, seq2: str,
              match_score: int = 1,
              mismatch_penalty: int = -1,
              gap_penalty: int = -2) -> tuple[str, str, str]:
   


   
    a1, a2, mid = [], [], []
    i, j = len(seq1), len(seq2)

    while i > 0 or j > 0:
        if i > 0 and j > 0:
            diag = match_score if seq1[i-1] == seq2[j-1] else mismatch_penalty
            if dp[i][j] == dp[i-1][j-1] + diag:          # came from diagonal
                a1.append(seq1[i-1]);  a2.append(seq2[j-1])
                mid.append('|' if seq1[i-1] == seq2[j-1] else ':')
                i -= 1;  j -= 1;  continue

        if i > 0 and dp[i][j] == dp[i-1][j] + gap_penalty:  # came from above
            a1.append(seq1[i-1]);  a2.append('-');  mid.append(' ')
            i -= 1
        else:                                                  # came from left
            a1.append('-');  a2.append(seq2[j-1]);  mid.append(' ')
            j -= 1

    aligned1 = ''.join(reversed(a1))
    aligned2 = ''.join(reversed(a2))
    midline  = ''.join(reversed(mid))
    return aligned1, aligned2, midline


def align(seq1: str, seq2: str,
          match_score: int = 1,
          mismatch_penalty: int = -1,
          gap_penalty: int = -2) -> dict:
    """Top-level DP alignment pipeline."""
    dp = build_dp_table(seq1, seq2, match_score, mismatch_penalty, gap_penalty)
    a1, a2, mid = traceback(dp, seq1, seq2, match_score, mismatch_penalty, gap_penalty)
    return {
        "seq1": seq1, "seq2": seq2,
        "aligned1": a1, "aligned2": a2, "midline": mid,
        "score": dp[len(seq1)][len(seq2)],
        "dp_table": dp,
    }


# 2. GREEDY APPROACH

def greedy_align(seq1: str, seq2: str,
                 match_score: int = 1,
                 mismatch_penalty: int = -1,
                 gap_penalty: int = -2) -> dict:
  


    a1, a2, mid = [], [], []
    i = j = 0
    score = 0

    while i < len(seq1) and j < len(seq2):
        if seq1[i] == seq2[j]:                             # match
            a1.append(seq1[i]); a2.append(seq2[j]); mid.append('|')
            score += match_score;  i += 1;  j += 1
        elif i+1 < len(seq1) and seq1[i+1] == seq2[j]:    # skip seq1[i]
            a1.append(seq1[i]); a2.append('-'); mid.append(' ')
            score += gap_penalty;  i += 1
        elif j+1 < len(seq2) and seq1[i] == seq2[j+1]:    # skip seq2[j]
            a1.append('-'); a2.append(seq2[j]); mid.append(' ')
            score += gap_penalty;  j += 1
        else:                                               # mismatch
            a1.append(seq1[i]); a2.append(seq2[j]); mid.append(':')
            score += mismatch_penalty;  i += 1;  j += 1

    while i < len(seq1):
        a1.append(seq1[i]); a2.append('-'); mid.append(' ')
        score += gap_penalty;  i += 1
    while j < len(seq2):
        a1.append('-'); a2.append(seq2[j]); mid.append(' ')
        score += gap_penalty;  j += 1

    return {"aligned1": ''.join(a1), "aligned2": ''.join(a2),
            "midline": ''.join(mid), "score": score}


# 3. DIVIDE & CONQUER (Hirschberg)

def _score_row(seq1: str, seq2: str,
               match_score: int, mismatch_penalty: int,
               gap_penalty: int) -> list[int]:
    """Return only the last DP row — O(n) space."""
    n = len(seq2)
    prev = [j * gap_penalty for j in range(n + 1)]
    for ch in seq1:
        curr = [0] * (n + 1)
        curr[0] = prev[0] + gap_penalty
        for j in range(1, n + 1):
            diag = match_score if ch == seq2[j-1] else mismatch_penalty
            curr[j] = max(prev[j-1] + diag,
                          prev[j]   + gap_penalty,
                          curr[j-1] + gap_penalty)
        prev = curr
    return prev


def hirschberg(seq1: str, seq2: str,
               match_score: int = 1,
               mismatch_penalty: int = -1,
               gap_penalty: int = -2) -> tuple[str, str]:
 

    m, n = len(seq1), len(seq2)

    if m == 0:  return '-' * n, seq2
    if n == 0:  return seq1, '-' * m
    if m == 1:
        r = align(seq1, seq2, match_score, mismatch_penalty, gap_penalty)
        return r["aligned1"], r["aligned2"]

    mid    = m // 2
    top    = _score_row(seq1[:mid], seq2,
                        match_score, mismatch_penalty, gap_penalty)
    bottom = _score_row(seq1[mid:][::-1], seq2[::-1],
                        match_score, mismatch_penalty, gap_penalty)

    scores = [t + b for t, b in zip(top, reversed(bottom))]
    split  = scores.index(max(scores))

    l1, l2 = hirschberg(seq1[:mid],  seq2[:split],
                        match_score, mismatch_penalty, gap_penalty)
    r1, r2 = hirschberg(seq1[mid:],  seq2[split:],
                        match_score, mismatch_penalty, gap_penalty)
    return l1 + r1, l2 + r2


# 4. DISPLAY HELPERS

def print_dp_table(dp: list[list[int]], seq1: str, seq2: str) -> None:
    """Pretty-print the DP matrix with row/column labels."""
    W = 5
    print("  " + " " * W + "  —  " + "".join(f"{c:>{W}}" for c in seq2))
    for i, row in enumerate(dp):
        label = "—" if i == 0 else seq1[i-1]
        print(f"  {label:>{W}}" + "".join(f"{v:>{W}}" for v in row))


def print_alignment(result: dict, width: int = 60) -> None:
    """Print the alignment in readable blocks."""
    a1, mid, a2 = result["aligned1"], result["midline"], result["aligned2"]
    score = result["score"]
    length = len(mid)
    matches    = mid.count('|')
    mismatches = mid.count(':')
    gaps       = mid.count(' ')

    print("\n" + "─" * 62)
    print("  Optimal Alignment")
    print("─" * 62)
    for s in range(0, length, width):
        e = s + width
        mid_vis = mid[s:e].replace('|','│').replace(':','·').replace(' ','·')
        print(f"  Seq1:  {a1[s:e]}")
        print(f"         {mid_vis}")
        print(f"  Seq2:  {a2[s:e]}\n")

    pct = lambda n: f"{100*n/length:.1f}%"
    print(f"  Score      : {score}")
    print(f"  Length     : {length}")
    print(f"  Matches    : {matches}  ({pct(matches)})")
    print(f"  Mismatches : {mismatches}  ({pct(mismatches)})")
    print(f"  Gaps       : {gaps}  ({pct(gaps)})")
    print("─" * 62 + "\n")


# 5. EXPERIMENTAL ANALYSIS

def random_dna(length: int, seed: int | None = None) -> str:
    if seed is not None:
        random.seed(seed)
    return ''.join(random.choices("ACGT", k=length))


def benchmark(lengths: list[int],
              match_score: int = 1,
              mismatch_penalty: int = -1,
              gap_penalty: int = -2) -> None:
    """Time DP alignment on sequences of increasing length."""
    print("\n" + "─" * 48)
    print("  Experimental Analysis — Runtime vs Length")
    print("─" * 48)
    print(f"  {'Length':>8}  {'Time (s)':>10}  {'Score':>8}  {'n²':>10}")
    print("  " + "-" * 44)
    for L in lengths:
        s1 = random_dna(L, seed=42)
        s2 = random_dna(L, seed=99)
        t0 = time.perf_counter()
        result = align(s1, s2, match_score, mismatch_penalty, gap_penalty)
        elapsed = time.perf_counter() - t0
        print(f"  {L:>8}  {elapsed:>10.4f}  {result['score']:>8}  {L*L:>10,}")
    print("─" * 48 + "\n")


# 6. DEMO ALL MODULES

if __name__ == "__main__":

    MATCH, MISMATCH, GAP = 1, -1, -2
    SEP = "═" * 62

    # ── Example 1: Core DP + traceback ────────────────────────
    print(SEP)
    print("  EXAMPLE 1 — Core DP Alignment (AGCTGAC vs GCTAC)")
    print(SEP)

    seq1, seq2 = "AGCTGAC", "GCTAC"
    result = align(seq1, seq2, MATCH, MISMATCH, GAP)

    print(f"\n  Seq1 : {seq1}")
    print(f"  Seq2 : {seq2}")
    print(f"  Scoring : match={MATCH}, mismatch={MISMATCH}, gap={GAP}\n")

    print("  DP Table:")
    print_dp_table(result["dp_table"], seq1, seq2)
    print_alignment(result)

    # Example 2: Greedy counterexample
    print(SEP)
    print("  EXAMPLE 2 — Greedy Counterexample (ACTTGAC vs CTTAC)")
    print(SEP)

    s1, s2 = "ACTTGAC", "CTTAC"
    dp_res = align(s1, s2, MATCH, MISMATCH, GAP)
    gr_res = greedy_align(s1, s2, MATCH, MISMATCH, GAP)

    print(f"\n  Seq1 : {s1}   Seq2 : {s2}\n")
    print(f"  DP alignment  (score = {dp_res['score']}):")
    print(f"    {dp_res['aligned1']}")
    print(f"    {dp_res['midline'].replace('|','│').replace(':','·')}")
    print(f"    {dp_res['aligned2']}\n")

    print(f"  Greedy alignment  (score = {gr_res['score']}):")
    print(f"    {gr_res['aligned1']}")
    print(f"    {gr_res['midline'].replace('|','│').replace(':','·')}")
    print(f"    {gr_res['aligned2']}\n")

    diff = dp_res['score'] - gr_res['score']
    if diff > 0:
        print(f"  ✔ DP outperforms Greedy by {diff} point(s).\n")
    else:
        print(f"  ~ Both methods score equally on this example.\n")

    # Example 3: Hirschberg D&C
    print(SEP)
    print("  EXAMPLE 3 — Hirschberg Divide & Conquer")
    print(SEP)

    s1, s2 = "GCATGCU", "GATTACA"
    dp_r      = align(s1, s2, MATCH, MISMATCH, GAP)
    h1, h2    = hirschberg(s1, s2, MATCH, MISMATCH, GAP)

    print(f"\n  NW  : {dp_r['aligned1']}   score = {dp_r['score']}")
    print(f"        {dp_r['aligned2']}")
    print(f"\n  D&C : {h1}")
    print(f"        {h2}\n")

    # Example 4: Effect of gap penalty
    print(SEP)
    print("  EXAMPLE 4 — Effect of Gap Penalty")
    print(SEP)

    s1 = random_dna(20, seed=7)
    s2 = random_dna(20, seed=13)
    print(f"\n  Seq1: {s1}")
    print(f"  Seq2: {s2}\n")
    print(f"  {'Gap':>5}  {'Score':>6}  {'Gaps':>5}  {'Alignment (S1)'}")
    print(f"  {'-'*5}  {'-'*6}  {'-'*5}  {'-'*25}")
    for gp in [-1, -2, -4, -8]:
        r = align(s1, s2, MATCH, MISMATCH, gp)
        print(f"  {gp:>5}  {r['score']:>6}  {r['midline'].count(' '):>5}  {r['aligned1'][:30]}")

    # Benchmark 
    benchmark([10, 50, 100, 200, 500, 1000])
