#!/bin/bash
set -euo pipefail
echo "[LAUNCH] Receptor Binding Motif Core@A"
jid_rfd=$(sbatch /data/homezvol1/inagakit/Projects/initbinder/tools/rfa_rfdiff/submit_rfa-diff_6M17_Receptor_Binding_Motif_Core_hsA.sh | awk '{print $4}')
jid_mpnn=$(sbatch --dependency=afterok:${jid_rfd} /data/homezvol1/inagakit/Projects/initbinder/tools/rfa_mpnn/submit_rfa-mpnn_6M17_Receptor_Binding_Motif_Core_hsA.sh | awk '{print $4}')
jid_af3s1=$(sbatch --dependency=afterok:${jid_mpnn} /data/homezvol1/inagakit/Projects/initbinder/tools/rfa_af3/submit_rfa-af3seed_6M17_Receptor_Binding_Motif_Core.sh | awk '{print $4}')
DESIGNS_PER_TASK=100 jid_af3s2=$(sbatch --dependency=afterok:${jid_af3s1} /data/homezvol1/inagakit/Projects/initbinder/tools/rfa_af3/submit_rfa-af3batch_6M17_Receptor_Binding_Motif_Core.sh | awk '{print $4}')
echo "[LAUNCH] Receptor Binding Motif Core@B"
jid_rfd=$(sbatch /data/homezvol1/inagakit/Projects/initbinder/tools/rfa_rfdiff/submit_rfa-diff_6M17_Receptor_Binding_Motif_Core_hsB.sh | awk '{print $4}')
jid_mpnn=$(sbatch --dependency=afterok:${jid_rfd} /data/homezvol1/inagakit/Projects/initbinder/tools/rfa_mpnn/submit_rfa-mpnn_6M17_Receptor_Binding_Motif_Core_hsB.sh | awk '{print $4}')
jid_af3s1=$(sbatch --dependency=afterok:${jid_mpnn} /data/homezvol1/inagakit/Projects/initbinder/tools/rfa_af3/submit_rfa-af3seed_6M17_Receptor_Binding_Motif_Core.sh | awk '{print $4}')
DESIGNS_PER_TASK=100 jid_af3s2=$(sbatch --dependency=afterok:${jid_af3s1} /data/homezvol1/inagakit/Projects/initbinder/tools/rfa_af3/submit_rfa-af3batch_6M17_Receptor_Binding_Motif_Core.sh | awk '{print $4}')
echo "[LAUNCH] Receptor Binding Motif Core@C"
jid_rfd=$(sbatch /data/homezvol1/inagakit/Projects/initbinder/tools/rfa_rfdiff/submit_rfa-diff_6M17_Receptor_Binding_Motif_Core_hsC.sh | awk '{print $4}')
jid_mpnn=$(sbatch --dependency=afterok:${jid_rfd} /data/homezvol1/inagakit/Projects/initbinder/tools/rfa_mpnn/submit_rfa-mpnn_6M17_Receptor_Binding_Motif_Core_hsC.sh | awk '{print $4}')
jid_af3s1=$(sbatch --dependency=afterok:${jid_mpnn} /data/homezvol1/inagakit/Projects/initbinder/tools/rfa_af3/submit_rfa-af3seed_6M17_Receptor_Binding_Motif_Core.sh | awk '{print $4}')
DESIGNS_PER_TASK=100 jid_af3s2=$(sbatch --dependency=afterok:${jid_af3s1} /data/homezvol1/inagakit/Projects/initbinder/tools/rfa_af3/submit_rfa-af3batch_6M17_Receptor_Binding_Motif_Core.sh | awk '{print $4}')
echo "[LAUNCH] RBM Flank and Crest@A"
jid_rfd=$(sbatch /data/homezvol1/inagakit/Projects/initbinder/tools/rfa_rfdiff/submit_rfa-diff_6M17_RBM_Flank_and_Crest_hsA.sh | awk '{print $4}')
jid_mpnn=$(sbatch --dependency=afterok:${jid_rfd} /data/homezvol1/inagakit/Projects/initbinder/tools/rfa_mpnn/submit_rfa-mpnn_6M17_RBM_Flank_and_Crest_hsA.sh | awk '{print $4}')
jid_af3s1=$(sbatch --dependency=afterok:${jid_mpnn} /data/homezvol1/inagakit/Projects/initbinder/tools/rfa_af3/submit_rfa-af3seed_6M17_RBM_Flank_and_Crest.sh | awk '{print $4}')
DESIGNS_PER_TASK=100 jid_af3s2=$(sbatch --dependency=afterok:${jid_af3s1} /data/homezvol1/inagakit/Projects/initbinder/tools/rfa_af3/submit_rfa-af3batch_6M17_RBM_Flank_and_Crest.sh | awk '{print $4}')
echo "[LAUNCH] RBM Flank and Crest@B"
jid_rfd=$(sbatch /data/homezvol1/inagakit/Projects/initbinder/tools/rfa_rfdiff/submit_rfa-diff_6M17_RBM_Flank_and_Crest_hsB.sh | awk '{print $4}')
jid_mpnn=$(sbatch --dependency=afterok:${jid_rfd} /data/homezvol1/inagakit/Projects/initbinder/tools/rfa_mpnn/submit_rfa-mpnn_6M17_RBM_Flank_and_Crest_hsB.sh | awk '{print $4}')
jid_af3s1=$(sbatch --dependency=afterok:${jid_mpnn} /data/homezvol1/inagakit/Projects/initbinder/tools/rfa_af3/submit_rfa-af3seed_6M17_RBM_Flank_and_Crest.sh | awk '{print $4}')
DESIGNS_PER_TASK=100 jid_af3s2=$(sbatch --dependency=afterok:${jid_af3s1} /data/homezvol1/inagakit/Projects/initbinder/tools/rfa_af3/submit_rfa-af3batch_6M17_RBM_Flank_and_Crest.sh | awk '{print $4}')
echo "[LAUNCH] RBM Flank and Crest@C"
jid_rfd=$(sbatch /data/homezvol1/inagakit/Projects/initbinder/tools/rfa_rfdiff/submit_rfa-diff_6M17_RBM_Flank_and_Crest_hsC.sh | awk '{print $4}')
jid_mpnn=$(sbatch --dependency=afterok:${jid_rfd} /data/homezvol1/inagakit/Projects/initbinder/tools/rfa_mpnn/submit_rfa-mpnn_6M17_RBM_Flank_and_Crest_hsC.sh | awk '{print $4}')
jid_af3s1=$(sbatch --dependency=afterok:${jid_mpnn} /data/homezvol1/inagakit/Projects/initbinder/tools/rfa_af3/submit_rfa-af3seed_6M17_RBM_Flank_and_Crest.sh | awk '{print $4}')
DESIGNS_PER_TASK=100 jid_af3s2=$(sbatch --dependency=afterok:${jid_af3s1} /data/homezvol1/inagakit/Projects/initbinder/tools/rfa_af3/submit_rfa-af3batch_6M17_RBM_Flank_and_Crest.sh | awk '{print $4}')
