"""Configuration for the family_finder pipeline."""

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path


@dataclass
class Config:
    # Tool paths
    orthofinder_bin: str = "orthofinder"
    mafft_bin: str = "mafft"
    fasttree_bin: str = "FastTree"
    iqtree_bin: str = "iqtree"
    codeml_bin: str = "codeml"
    pal2nal_bin: str = "pal2nal.pl"
    hmmbuild_bin: str = "hmmbuild"
    hmmsearch_bin: str = "hmmsearch"
    hmmpress_bin: str = "hmmpress"

    # HMMER rescue parameters
    hmmer_rescue: bool = False
    hmmer_evalue: float = 1e-5

    # Per-round profile assignment (issue #13; steps/profile_assign.py)
    profile_assign_per_round: bool = False  # off by default: behavior change is opt-in
    profile_min_coverage: float = 0.5        # min fraction of HMM covered by domains
    profile_min_query_coverage: float = 0.4  # min fraction of query covered by domains
    profile_margin_bits: float = 10.0        # best-vs-second bit-score margin (absolute)
    profile_margin_frac: float = 0.05        # ... or this fraction of best bits, if larger
    profile_reassign_margin_nbits: float = 0.15  # length-normalized bits margin for moves
    merge_min_reciprocal: float = 0.6        # min reciprocal cross-hit fraction for merges

    # Pipeline parameters
    max_rounds: int = 10
    min_orthogroup_size: int = 4  # floor to START align/tree work (cost gate)
    min_family_size: int = 2      # floor to EMIT a family after pruning
    convergence_threshold: int = 5
    convergence_no_new_families: int = 2

    # Pruning parameters
    distance_ratio_threshold: float = 5.0
    # Pruning criterion (issue #14):
    #   "relative" (default) — terminal-branch-subtracted internal-path ratio,
    #     normalized by the family-wide median ratio. Prune iff
    #     S_norm > prune_relative_threshold AND S_raw > prune_score_floor.
    #   "absolute" — legacy behavior: median(d_obs/d_exp) compared against
    #     distance_ratio_threshold (kept for reproducibility of old runs).
    prune_criterion: str = "relative"
    prune_relative_threshold: float = 3.0
    prune_score_floor: float = 2.0
    min_species_for_pruning: int = 3
    treeshrink_quantile: float = 0.05  # TreeShrink quantile (lower = stricter)

    # OrthoFinder parameters
    orthofinder_threads: int = 8
    # MCL inflation (-I). OrthoFinder v3 default is 1.2 (v2.x was 1.5); passing
    # it explicitly makes the value visible in logs and sweepable via config.
    orthofinder_inflation: float = 1.2
    # Search program (-S). Empty = OrthoFinder default (diamond --more-sensitive).
    # Options include: diamond, diamond_ultra_sens, blast, mmseqs.
    orthofinder_search_program: str = ""
    # Stop after orthogroup inference (-og). Skips per-OG gene trees, species
    # tree, and ortholog inference that this pipeline never reads — a large
    # wall-clock saving. Caveat: also skips Comparative_Genomics_Statistics/
    # (per-species stats). NOTE: -og is accepted by OrthoFinder v3.1.3 but
    # absent from its -h output — do not probe support by parsing -h.
    orthofinder_stop_after_orthogroups: bool = True
    # Reuse a previous WorkingDirectory of cached DIAMOND results (-b).
    # Only valid for re-running the SAME sequence set (e.g. inflation sweeps
    # on round 1) — later rounds have different pools and cannot reuse.
    orthofinder_reuse_blast_dir: str = ""
    orthofinder_extra_args: str = ""

    # MAFFT parameters
    mafft_strategy: str = "auto"

    # Tree builder: "fasttree" or "iqtree"
    tree_builder: str = "fasttree"

    # Parallelism
    n_workers: int = 8

    # codeml
    run_codeml: bool = False
    codeml_models: list = field(default_factory=lambda: ["M0", "M1a", "M2a", "M7", "M8"])

    # DeepLoc (issue #16: retargeting / neofunctionalization module)
    deeploc_bin: str = "deeploc2"
    deeploc_model: str = "Fast"      # Fast | Accurate
    deeploc_device: str = ""         # e.g. "cuda"; empty = DeepLoc default
    deeploc_min_prob: float = 0.5    # class prob below this = uninformative leaf

    # Pseudogene detection
    pseudogene_detection: bool = True
    pseudogene_truncation_threshold: float = 0.5  # flag if gene < 50% of family median
    pseudogene_species_filter: str = ""           # restrict to one species (e.g. "Ococ")

    # Structural tier + annotation layers (issues #17, #18, #20).
    # The *_cmd fields are command TEMPLATES ("{fasta}"/"{input}"/"{outdir}"
    # placeholders, shlex-split after formatting) because the cluster tooling
    # varies (esm-extract vs custom scripts, separate conda envs). Empty
    # template = that step is unavailable and its wrapper raises.
    esm_embed_cmd: str = ""    # e.g. "esm-extract esm2_t33_650M_UR50D {fasta} {outdir} --include mean"
    esmfold_cmd: str = ""      # e.g. "esm-fold -i {fasta} -o {outdir}"
    foldseek_bin: str = "foldseek"
    plddt_flag_below: float = 50.0    # mean pLDDT below this -> pseudogene-evidence flag
    tier3_min_tmscore: float = 0.5    # min alntmscore for a tier-3 structural assignment
    tier3_margin_tm: float = 0.05     # best-vs-second alntmscore margin (ties -> ambiguous)
    embed_tiebreak_min_delta: float = 0.02  # mean-cosine margin for the #13 ambiguous tie-break
    glm_score_cmd: str = ""    # plant gLM adapter (PlantCaduceus/AgroNT/Evo 2), issue #18
    ecforest_cmd: str = ""     # ESM-ECForest wrapper (own conda env), issue #20

    # Gene ID format: species extracted from prefix before this delimiter
    species_delimiter: str = "_"

    @classmethod
    def from_json(cls, path: str) -> "Config":
        import logging
        with open(path) as f:
            data = json.load(f)
        unknown = set(data.keys()) - set(cls.__dataclass_fields__.keys())
        if unknown:
            logging.getLogger("family_finder").warning(
                f"Unknown config keys in {path} (ignored): {unknown}"
            )
        return cls(**{k: v for k, v in data.items() if k in cls.__dataclass_fields__})

    def to_json(self, path: str):
        with open(path, "w") as f:
            json.dump(asdict(self), f, indent=2)
