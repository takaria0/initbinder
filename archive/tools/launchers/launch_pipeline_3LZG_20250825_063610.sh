#!/bin/bash
set -euo pipefail
echo "[LAUNCH] Conserved Stem Helix@A"
jid_rfd=$(sbatch /data/homezvol1/inagakit/Projects/initbinder/tools/rfa_rfdiff/submit_rfa-diff_3LZG_Conserved_Stem_Helix_hsA.sh | awk '{print $4}')
jid_mpnn=$(sbatch --dependency=afterok:${jid_rfd} /data/homezvol1/inagakit/Projects/initbinder/tools/rfa_mpnn/submit_rfa-mpnn_3LZG_Conserved_Stem_Helix_hsA.sh | awk '{print $4}')
jid_af3s1=$(sbatch --dependency=afterok:${jid_mpnn} /data/homezvol1/inagakit/Projects/initbinder/tools/rfa_af3/submit_rfa-af3seed_3LZG_Conserved_Stem_Helix.sh | awk '{print $4}')
DESIGNS_PER_TASK=100 jid_af3s2=$(sbatch --dependency=afterok:${jid_af3s1} /data/homezvol1/inagakit/Projects/initbinder/tools/rfa_af3/submit_rfa-af3batch_3LZG_Conserved_Stem_Helix.sh | awk '{print $4}')
echo "[LAUNCH] Conserved Stem Helix@B"
jid_rfd=$(sbatch /data/homezvol1/inagakit/Projects/initbinder/tools/rfa_rfdiff/submit_rfa-diff_3LZG_Conserved_Stem_Helix_hsB.sh | awk '{print $4}')
jid_mpnn=$(sbatch --dependency=afterok:${jid_rfd} /data/homezvol1/inagakit/Projects/initbinder/tools/rfa_mpnn/submit_rfa-mpnn_3LZG_Conserved_Stem_Helix_hsB.sh | awk '{print $4}')
jid_af3s1=$(sbatch --dependency=afterok:${jid_mpnn} /data/homezvol1/inagakit/Projects/initbinder/tools/rfa_af3/submit_rfa-af3seed_3LZG_Conserved_Stem_Helix.sh | awk '{print $4}')
DESIGNS_PER_TASK=100 jid_af3s2=$(sbatch --dependency=afterok:${jid_af3s1} /data/homezvol1/inagakit/Projects/initbinder/tools/rfa_af3/submit_rfa-af3batch_3LZG_Conserved_Stem_Helix.sh | awk '{print $4}')
echo "[LAUNCH] Conserved Stem Helix@C"
jid_rfd=$(sbatch /data/homezvol1/inagakit/Projects/initbinder/tools/rfa_rfdiff/submit_rfa-diff_3LZG_Conserved_Stem_Helix_hsC.sh | awk '{print $4}')
jid_mpnn=$(sbatch --dependency=afterok:${jid_rfd} /data/homezvol1/inagakit/Projects/initbinder/tools/rfa_mpnn/submit_rfa-mpnn_3LZG_Conserved_Stem_Helix_hsC.sh | awk '{print $4}')
jid_af3s1=$(sbatch --dependency=afterok:${jid_mpnn} /data/homezvol1/inagakit/Projects/initbinder/tools/rfa_af3/submit_rfa-af3seed_3LZG_Conserved_Stem_Helix.sh | awk '{print $4}')
DESIGNS_PER_TASK=100 jid_af3s2=$(sbatch --dependency=afterok:${jid_af3s1} /data/homezvol1/inagakit/Projects/initbinder/tools/rfa_af3/submit_rfa-af3batch_3LZG_Conserved_Stem_Helix.sh | awk '{print $4}')
