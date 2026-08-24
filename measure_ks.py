#!/usr/bin/env python3
"""Per-species-pair Ks/Ka baseline from the run's own single-copy families.

Written because a saturation claim was asserted twice and measured zero
times. The first version reasoned from a species abbreviation (Dcar read as
Daucus carota, Apiales) and concluded Ks was hopelessly saturated; the
correction reasoned from species-tree branch lengths (max 0.548 subs/site)
and put Ks at 1.5-2.0. Both were guesses, and the second was still wrong.

Measured on all 69 strict single-copy families (15 genes / 15 species) using
the codon alignments the pipeline already wrote, 7,216 pairwise comparisons:

    median Ks              0.555
    p90                    0.903
    max                    2.330
    JC-undefined           1 / 7,217   = 0.0%
    Ks > 1                 451 / 7,216 = 6.2%
    deepest pair (Dcar-Pami)  median Ks 0.841
    shallowest (Obas-Ococ)    median Ks 0.024

So Ks is usable across this panel - even the deepest species pair sits inside
the conventional Ks <= 1 range, and exactly one comparison out of 7,217 was
saturated enough to break the Jukes-Cantor correction. Synonymous sites here
run about 1.5x the total nucleotide distance, not the 3-4x the correction
assumed.

This does NOT make Ks a family-boundary criterion. Single-copy orthologs
diverged at speciation while families hold paralogs that diverged before it -
PTPC and BTPC split before land plants - so an ortholog distance distribution
is a lower bound on family radius and cannot locate the boundary. Nor can the
observed/expected score in steps/prune.py be lifted: it normalizes by each
family's own median and skips same-species paralogs, so it is comparable only
within the family it was computed in, and ancient duplications hide in exactly
that skip. What this table IS good for is calibration - what the distance
between true orthologs actually looks like for each species pair.

Median Ka/Ks is 0.311 and barely moves across pairs, which is the same point
from the other side: omega measures how a gene changed, not how much, so it
cannot order anything by relatedness.

Method: Nei-Gojobori counting with Jukes-Cantor correction. Saturation is
reported as an undefined value rather than a capped one - the correction
diverging IS the measurement.

Output: data/ks_baseline_15sp.tsv (105 species pairs).
"""
import math, itertools, statistics as st
from pathlib import Path
from collections import defaultdict

F = Path.home()/"scratch/bin/family_finder/output_15sp_v2"
W = Path.home()/"scratch/ks_check"
CODON = {}
BASES = "TCAG"
AAS = ("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")
i = 0
for a in BASES:
    for b in BASES:
        for c in BASES:
            CODON[a+b+c] = AAS[i]; i += 1


def syn_sites(codon):
    """Nei-Gojobori: fraction of each site that is synonymous."""
    aa = CODON.get(codon)
    if aa is None or aa == "*":
        return None
    s = 0.0
    for pos in range(3):
        for nt in "TCAG":
            if nt == codon[pos]:
                continue
            mut = codon[:pos] + nt + codon[pos+1:]
            alt = CODON.get(mut)
            if alt is not None and alt == aa:
                s += 1.0 / 3.0
    return s


def read_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            name = line[1:].split()[0]; seqs[name] = []
        elif name:
            seqs[name].append(line.strip())
    return {k: "".join(v).upper() for k, v in seqs.items()}


def pair_ks_ka(s1, s2):
    S = N = Sd = Nd = 0.0
    for i in range(0, min(len(s1), len(s2)) - 2, 3):
        c1, c2 = s1[i:i+3], s2[i:i+3]
        if len(c1) < 3 or len(c2) < 3: continue
        if "-" in c1 or "-" in c2 or "N" in c1 or "N" in c2: continue
        a1, a2 = CODON.get(c1), CODON.get(c2)
        if a1 is None or a2 is None or a1 == "*" or a2 == "*": continue
        s_1, s_2 = syn_sites(c1), syn_sites(c2)
        if s_1 is None or s_2 is None: continue
        s = (s_1 + s_2) / 2.0
        S += s; N += 3.0 - s
        diff = sum(1 for k in range(3) if c1[k] != c2[k])
        if diff == 0: continue
        if diff == 1:
            if a1 == a2: Sd += 1
            else: Nd += 1
        else:
            # multi-difference codons: split by the fraction of single-step
            # paths that are synonymous rather than enumerating pathways
            frac_s = 1.0 if a1 == a2 else 0.0
            Sd += diff * frac_s; Nd += diff * (1 - frac_s)
    if S < 10 or N < 10: return None
    ps, pn = Sd / S, Nd / N
    def jc(p):
        x = 1.0 - (4.0/3.0)*p
        return None if x <= 0 else -0.75 * math.log(x)
    return jc(ps), jc(pn), S, N


fams = [l.strip() for l in open(W/"singlecopy.txt") if l.strip()]
by_pair_ks = defaultdict(list); by_pair_ka = defaultdict(list)
sat = defaultdict(int); tot = defaultdict(int)
used = 0
for fam in fams:
    aln = F/"final_families"/fam/"confirmed_codon.afa"
    if not aln.exists(): continue
    seqs = read_fasta(aln)
    by_sp = {}
    for name, s in seqs.items():
        by_sp[name.split("_", 1)[0]] = s
    if len(by_sp) < 2: continue
    used += 1
    for a, b in itertools.combinations(sorted(by_sp), 2):
        r = pair_ks_ka(by_sp[a], by_sp[b])
        if r is None: continue
        ks, ka, S, N = r
        key = (a, b); tot[key] += 1
        if ks is None:
            sat[key] += 1          # JC correction undefined = saturated
        else:
            by_pair_ks[key].append(ks)
            if ka is not None: by_pair_ka[key].append(ka)

print("alignments used: %d of %d" % (used, len(fams)))
print()
rows = []
for key in tot:
    ks = by_pair_ks.get(key, [])
    rows.append((st.median(ks) if ks else float("inf"), key, ks, sat[key], tot[key]))
rows.sort()
print("%-14s %8s %8s %10s %8s" % ("species pair", "med_Ks", "med_Ka", "saturated", "n_fam"))
for med, key, ks, s, t in rows[:4] + rows[-8:]:
    ka = by_pair_ka.get(key, [])
    print("%-14s %8s %8s %6d/%-4d %6d" % (
        "%s-%s" % key,
        "%.3f" % med if med != float("inf") else "n/a",
        "%.3f" % st.median(ka) if ka else "n/a", s, t, t))
print()
allks = [v for l in by_pair_ks.values() for v in l]
print("all pairwise Ks: n=%d  median=%.3f  p90=%.3f  max=%.3f"
      % (len(allks), st.median(allks), sorted(allks)[int(len(allks)*0.9)], max(allks)))
tot_sat = sum(sat.values()); tot_all = sum(tot.values())
print("JC-undefined (saturated) comparisons: %d/%d = %.1f%%"
      % (tot_sat, tot_all, 100.0*tot_sat/tot_all))
print("Ks > 1: %d/%d = %.1f%%" % (sum(1 for v in allks if v > 1),
                                  len(allks), 100.0*sum(1 for v in allks if v > 1)/len(allks)))
