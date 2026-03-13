"""Deterministic Golden Gate adapter builder for binder export."""

from __future__ import annotations

import hashlib
import re
import secrets
from dataclasses import dataclass
from typing import Optional, Tuple, Union

DNA_RE = re.compile(r"^[ACGTacgt]+$")


@dataclass(frozen=True)
class GoldenGateSeqBuilder:
    """
    Standalone builder for:
      FULL = LEFT + INSERT + RIGHT

    Default adapter pattern:
      LEFT  = [barcode_len nt deterministic handle] + left_suffix
      RIGHT = right_const
    """

    left_suffix: str = "ggtctcGCCTCA"
    right_const: str = "GGACAGgagaccatcgatcg"
    barcode_len: int = 20
    barcode_prefix: str = "GC"

    enforce_constraints: bool = False
    gc_min: float = 0.40
    gc_max: float = 0.60
    max_homopolymer: int = 4
    banned_motifs: Tuple[str, ...] = ("GGTCTC", "GAGACC")
    max_tries: int = 10000

    def build(
        self,
        insert_seq: str,
        seed: Optional[Union[int, str]] = None,
        adapter: Optional[str] = None,
    ) -> dict:
        insert = self._validate_insert(insert_seq)

        used_seed: Union[int, str]
        if seed is None:
            used_seed = secrets.randbits(64)
        else:
            used_seed = seed

        if adapter is not None:
            left, right = self._split_adapter(adapter)
            full = f"{left}{insert}{right}"
            return {
                "full_seq": full,
                "seed": used_seed,
                "barcode": None,
                "left": left,
                "right": right,
                "adapter_seqs": f"{left}/{right}",
            }

        barcode = self._make_barcode(used_seed, self.barcode_len)
        left = barcode + self.left_suffix
        right = self.right_const
        full = f"{left}{insert}{right}"
        return {
            "full_seq": full,
            "seed": used_seed,
            "barcode": barcode,
            "left": left,
            "right": right,
            "adapter_seqs": f"{left}/{right}",
            "barcode_gc": self._gc_fraction(barcode),
            "barcode_max_homopolymer": self._max_homopolymer_run(barcode),
        }

    def _validate_insert(self, seq: str) -> str:
        seq = seq.strip()
        if not DNA_RE.match(seq):
            raise ValueError("insert_seq must contain only A/C/G/T (no gaps, no Ns).")
        return seq.upper()

    def _split_adapter(self, adapter_seqs: str) -> Tuple[str, str]:
        if "/" not in adapter_seqs:
            raise ValueError("adapter must contain '/' (insertion point), e.g. LEFT/RIGHT")
        left, right = adapter_seqs.split("/", 1)
        return left, right

    def _seed_to_bytes(self, seed: Union[int, str]) -> bytes:
        return str(seed).encode("utf-8") if isinstance(seed, int) else seed.encode("utf-8")

    def _deterministic_dna(self, seed: Union[int, str], length: int, counter: int) -> str:
        bases = "ACGT"
        seed_bytes = self._seed_to_bytes(seed)
        out = []
        block = 0
        while len(out) < length:
            digest = hashlib.sha256(
                seed_bytes
                + b"|barcode|"
                + str(counter).encode("utf-8")
                + b"|"
                + str(block).encode("utf-8")
            ).digest()
            block += 1
            for b in digest:
                out.append(bases[b & 0b11])
                if len(out) >= length:
                    break
        return "".join(out)

    def _make_barcode(self, seed: Union[int, str], length: int) -> str:
        prefix = str(self.barcode_prefix or "").upper()
        if len(prefix) > length:
            raise ValueError(
                f"barcode_len={length} is shorter than barcode_prefix='{prefix}' ({len(prefix)} nt)."
            )

        def _candidate(counter: int) -> str:
            tail_len = length - len(prefix)
            if tail_len <= 0:
                return prefix
            tail = self._deterministic_dna(seed, length=tail_len, counter=counter)
            return f"{prefix}{tail}"

        if not self.enforce_constraints:
            return _candidate(counter=0)
        for counter in range(self.max_tries):
            barcode = _candidate(counter=counter)
            if self._passes_constraints(barcode):
                return barcode
        raise ValueError(
            f"Could not find a barcode meeting constraints within max_tries={self.max_tries}. "
            "Relax constraints (GC range/homopolymer/banned motifs) or increase max_tries."
        )

    def _passes_constraints(self, barcode: str) -> bool:
        seq = barcode.upper()
        if any(motif.upper() in seq for motif in self.banned_motifs):
            return False
        gcf = self._gc_fraction(seq)
        if not (self.gc_min <= gcf <= self.gc_max):
            return False
        if self._max_homopolymer_run(seq) > self.max_homopolymer:
            return False
        return True

    def _gc_fraction(self, seq: str) -> float:
        seq = seq.upper()
        if not seq:
            return 0.0
        gc = sum(1 for c in seq if c in ("G", "C"))
        return gc / len(seq)

    def _max_homopolymer_run(self, seq: str) -> int:
        seq = seq.upper()
        best = 1
        cur = 1
        for idx in range(1, len(seq)):
            if seq[idx] == seq[idx - 1]:
                cur += 1
                if cur > best:
                    best = cur
            else:
                cur = 1
        return best
