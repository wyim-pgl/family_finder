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

    # Pipeline parameters
    max_rounds: int = 10
    min_orthogroup_size: int = 4
    convergence_threshold: int = 5
    convergence_no_new_families: int = 2

    # Pruning parameters
    distance_ratio_threshold: float = 5.0
    min_species_for_pruning: int = 3
    treeshrink_quantile: float = 0.05  # TreeShrink quantile (lower = stricter)

    # OrthoFinder parameters
    orthofinder_threads: int = 8
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
