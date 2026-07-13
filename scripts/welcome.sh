#!/usr/bin/env bash
# welcome.sh – printed at every interactive session start (Day 2 / course container)

cat << 'EOF'

╔══════════════════════════════════════════════════════════════════════════════╗
║                 THE EPI-CODE — Florence Training School 2026                 ║
║      "Epitranscriptomics and Epigenomics in Plants and Microorganisms"       ║
╠══════════════════════════════════════════════════════════════════════════════╣
║  Course Environment (Day 2 – WGBS)  │  University of Messina, ChiBioFarAm    ║
║  Instructors: D. Giosa & L. Giuffrè                                          ║
╠══════════════════════════════════════════════════════════════════════════════╣
║  QUICK REFERENCE                                                             ║
║  ──────────────────────────────────────────────────────────────────────────  ║
║  Day 1 – lncRNA / mRNA  →  runs on your HOST, in its own image:              ║
║      docker pull leogiuffre/lncrna-mnps-workshop:1.0                         ║
║      results land in   ./results/day1_lncRNA                                 ║
║                                                                              ║
║  Day 2 – WGBS methylation  →  runs HERE, in this container:                  ║
║      cd $RESULTS_DIR/day2_methylation                                        ║
║                                                                              ║
║  Data directory   →  $DATA_DIR                                               ║
║  Reference        →  $REF_DIR                                                ║
║  Scripts          →  $SCRIPTS_DIR                                            ║
║  Exercises        →  /course/exercises                                       ║
╠══════════════════════════════════════════════════════════════════════════════╣
║  Day 2 setup (once):  bash $SCRIPTS_DIR/00_setup_reference.sh                ║
║  Day 2 data  (once):  bash $DATA_DIR/download_data.sh                        ║
║  Day 2 pipeline:      bash $SCRIPTS_DIR/run_day2.sh                          ║
╚══════════════════════════════════════════════════════════════════════════════╝

EOF

echo "Working directory: $(pwd)"
echo "Available disk space: $(df -h /course | awk 'NR==2{print $4}') free"
echo ""
