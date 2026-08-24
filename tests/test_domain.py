"""Domain architecture as a scalable source of NEGATIVES (issue #47).

The outgroup problem is that a defensible negative - something known to lie
outside the family - has to come from somewhere, and PEPC was the lucky case:
the literature had already split plant-type from bacterial-type PEPC, and
anchors.tsv already labelled ATH_AT1G68750.1 bacterial. The other 23,743
families have no such gift.

Pfam does, for every family at once. Two proteins built from different domain
architectures are different families, and that is curated rather than inferred
from the very graph whose boundaries are in question - which is what sank
drawing outgroups from below-cut vote edges.

The whole difficulty is telling a real architectural difference from this
repository's signature failure: truncated and mis-annotated models lose
domains, and a lost domain looks exactly like a different architecture. So
nothing here compares two genes. A family's signature is built only from
domains carried by a MAJORITY of its members, and two families are called
incompatible only when a domain is common in one and near-absent in the other.
A short model cannot move a majority on its own.
"""
import sys
import types
from pathlib import Path

sys.modules.setdefault("ete4", types.ModuleType("ete4"))
sys.modules["ete4"].Tree = object
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from steps.domain import (architecture_signature, incompatible_families,
                          parse_pfam_domtblout, pick_domain_outgroup)


def _domtbl(path, rows):
    """rows: (target_pfam_acc, query_gene, ali_from, ali_to, i_evalue)"""
    lines = ["# hmmscan domtblout"]
    for acc, gene, start, end, ev in rows:
        f = ["Dom", acc, "-", gene, "-", "300", "1e-40", "100.0", "0", "1",
             "1", str(ev), str(ev), "80.0", "0.0", "1", "50",
             str(start), str(end)]
        lines.append(" ".join(f))
    lines.append("# [ok]")
    path.write_text("\n".join(lines) + "\n")


# ---------------------------------------------------------------------------
# parsing
# ---------------------------------------------------------------------------

def test_parse_orders_domains_along_the_protein(tmp_path):
    # Arrange: hmmscan emits per profile, so rows arrive out of positional
    # order; architecture is meaningless unless they are sorted by position
    p = tmp_path / "d.domtblout"
    _domtbl(p, [("PF00016.1", "Sp1_a", 200, 260, 1e-30),
                ("PF00311.1", "Sp1_a", 10, 90, 1e-50)])

    # Act
    got = parse_pfam_domtblout(p)

    # Assert: accession version stripped - PF00311.12 and PF00311.13 are the
    # same domain and a version bump must not look like a new architecture
    assert got == {"Sp1_a": ["PF00311", "PF00016"]}


def test_parse_applies_the_evalue_gate(tmp_path):
    # Arrange
    p = tmp_path / "d.domtblout"
    _domtbl(p, [("PF00311.1", "Sp1_a", 10, 90, 1e-50),
                ("PF09999.1", "Sp1_a", 100, 140, 1.0)])

    # Act
    got = parse_pfam_domtblout(p, max_evalue=1e-5)

    # Assert
    assert got == {"Sp1_a": ["PF00311"]}


def test_parse_keeps_repeats_but_not_overlaps(tmp_path):
    # Arrange: a genuine tandem repeat is architecture; two profiles hitting
    # the same span are one region, and the weaker must not inflate it
    p = tmp_path / "d.domtblout"
    _domtbl(p, [("PF00400.1", "Sp1_a", 10, 50, 1e-40),
                ("PF00400.1", "Sp1_a", 60, 100, 1e-40),
                ("PF00401.1", "Sp1_a", 15, 45, 1e-10)])

    # Act
    got = parse_pfam_domtblout(p)

    # Assert
    assert got == {"Sp1_a": ["PF00400", "PF00400"]}


# ---------------------------------------------------------------------------
# family signature
# ---------------------------------------------------------------------------

def test_signature_keeps_only_domains_a_majority_carries():
    # Arrange: three of four members carry PF00311; one truncated model also
    # shows a spurious PF09999
    arch = {"g1": ["PF00311"], "g2": ["PF00311"], "g3": ["PF00311"],
            "g4": ["PF09999"]}

    # Act
    sig = architecture_signature(["g1", "g2", "g3", "g4"], arch)

    # Assert
    assert sig == {"PF00311"}


def test_signature_of_a_family_whose_members_have_no_domains_is_empty():
    # Arrange / Act / Assert: no signature means no claim, not a merge
    assert architecture_signature(["g1", "g2"], {}) == set()


def test_a_lone_truncated_member_cannot_move_the_signature():
    # Arrange: the failure this repository keeps hitting - a short model drops
    # its domains, and the loss reads as a different architecture
    arch = {f"g{i}": ["PF00311", "PF00016"] for i in range(9)}
    arch["g_trunc"] = []

    # Act
    sig = architecture_signature(sorted(arch), arch)

    # Assert
    assert sig == {"PF00311", "PF00016"}


# ---------------------------------------------------------------------------
# incompatibility -> negative anchors
# ---------------------------------------------------------------------------

def test_families_with_disjoint_signatures_are_incompatible():
    # Arrange
    sigs = {"A": {"PF00311"}, "B": {"PF00400"}}

    # Act / Assert
    assert incompatible_families("A", sigs) == ["B"]


def test_a_shared_domain_is_not_enough_to_be_compatible_or_not():
    # Arrange: sharing a domain says related, not same family; the rule only
    # asserts the NEGATIVE, so overlap simply yields no claim
    sigs = {"A": {"PF00311", "PF00016"}, "B": {"PF00311"}}

    # Act / Assert
    assert incompatible_families("A", sigs) == []


def test_an_empty_signature_never_makes_anything_incompatible():
    # Arrange: absence of domain evidence must not become evidence of absence
    sigs = {"A": {"PF00311"}, "B": set()}

    # Act / Assert
    assert incompatible_families("A", sigs) == []
    assert incompatible_families("B", sigs) == []


def test_domain_outgroup_prefers_the_nearest_incompatible_family():
    # Arrange: an outgroup should be just outside, not arbitrarily far - a
    # distant monophyletic panel merges everything it brackets
    sigs = {"A": {"PF00311"}, "N1": {"PF00400"}, "N2": {"PF00500"}}
    proximity = {"N1": 0.4, "N2": 0.1}

    # Act
    got = pick_domain_outgroup(["A"], sigs, proximity, n=1)

    # Assert
    assert got == ["N1"]


def test_domain_outgroup_never_returns_a_cluster_member():
    # Arrange
    sigs = {"A": {"PF00311"}, "B": {"PF00400"}}

    # Act / Assert: B is incompatible with A but both are in the cluster
    assert pick_domain_outgroup(["A", "B"], sigs, {"B": 0.9}, n=2) == []
