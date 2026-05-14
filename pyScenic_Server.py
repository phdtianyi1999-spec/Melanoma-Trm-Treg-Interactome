#!/usr/bin/env python
# -*- coding: utf-8 -*-

# The top two lines are standard headers for Linux scripts, specifying the interpreter and encoding

import os
import sys
import subprocess
import multiprocessing

# ================= Core Configuration Area (Must be modified according to your server) =================

# --- 1. Working Directory Setup ---
# Linux path example: "/home/username/project/scenic_analysis/"
# Please navigate to your data folder in the server terminal, enter the `pwd` command to get the absolute path, and paste it here
work_dir = "/home/zhangtianyi/project/scenic_data/"  # <--- Please modify to your actual path!!!

# --- 2. Performance Settings ---
# Linux servers usually have many cores, can be automatically detected using multiprocessing.cpu_count()
try:
    total_cores = multiprocessing.cpu_count()
    # Leave 2-4 cores for the system, use the rest to run GRN
    GRN_WORKERS = max(1, total_cores - 4) 
    # Other steps don't need that many cores, 12-20 is enough
    OTHER_WORKERS = 12 
except:
    GRN_WORKERS = 12
    OTHER_WORKERS = 8

# --- 3. File Path Settings ---
# Please ensure these files have been uploaded to the work_dir directory above
INPUT_LOOM = "input.loom"
TF_FILE = "allTFs_hg38.txt"
DB_FILES = [
    "hg38_refseq-r80_10kb_up_and_down_tss.mc9nr.feather",
    "hg38_refseq-r80_500bp_up_and_100bp_down_tss.mc9nr.feather"
]
OUTPUT_LOOM = "output_scenic.loom"

# ======================================================================

def run_pipeline():
    # 1. Change working directory
    if os.path.exists(work_dir):
        os.chdir(work_dir)
        print(f"✅ Switched working directory to: {os.getcwd()}")
    else:
        print(f"❌ Directory does not exist: {work_dir}")
        print("Please modify the 'work_dir' variable in the script to the correct Linux path!")
        sys.exit(1)

    # 2. Print hardware information
    print(f"🚀 Detected server cores: {multiprocessing.cpu_count()}")
    print(f"⚙️  GRN step will use: {GRN_WORKERS} cores")
    print(f"⚙️  Other steps will use: {OTHER_WORKERS} cores")

    # 3. Check if files exist
    missing = []
    files_to_check = [INPUT_LOOM, TF_FILE] + DB_FILES
    for f in files_to_check:
        if not os.path.exists(f): missing.append(f)

    if missing:
        print("\n❌ Missing files, please check the upload status:")
        for f in missing: print(f" - {f}")
        sys.exit(1)

    # 4. Define execution function
    def run_command(cmd_str, step_name):
        print(f"\n{'='*15} {step_name} Started {'='*15}")
        print(f"Executing command: {cmd_str}\n")
        try:
            # On Linux, shell=True allows using commands from the conda environment
            subprocess.run(cmd_str, shell=True, check=True)
            print(f"\n✅ {step_name} successfully completed!")
        except subprocess.CalledProcessError as e:
            print(f"\n❌ {step_name} failed (Exit code: {e.returncode})")
            sys.exit(1)

    # ================= Pipeline Execution =================

    # Step 1: GRN
    print(f"--- Starting Step 1: GRN ---")
    # Added --seed to ensure reproducible results, --sparse to save memory
    cmd_grn = f"pyscenic grn --num_workers {GRN_WORKERS} --output adj.tsv --method grnboost2 --seed 777 {INPUT_LOOM} {TF_FILE}"
    run_command(cmd_grn, "Step 1: GRN")

    # Step 2: CTX
    print(f"--- Starting Step 2: CTX ---")
    db_args = " ".join(DB_FILES)
    # Note: If your loom file does not have 'HGNC' symbols, you may need to add --no_pruning or other parameters; defaults are usually fine
    cmd_ctx = f"pyscenic ctx adj.tsv {db_args} --expression_mtx_fname {INPUT_LOOM} --output regulons.csv --num_workers {OTHER_WORKERS} --mask_dropouts"
    run_command(cmd_ctx, "Step 2: CTX")

    # Step 3: AUCell
    print(f"--- Starting Step 3: AUCell ---")
    cmd_auc = f"pyscenic aucell {INPUT_LOOM} regulons.csv --output {OUTPUT_LOOM} --num_workers {OTHER_WORKERS} --seed 777"
    run_command(cmd_auc, "Step 3: AUCell")

    print("\n🎉🎉🎉 Entire pipeline finished!")
    print(f"Results saved to: {os.path.join(work_dir, OUTPUT_LOOM)}")

if __name__ == "__main__":
    run_pipeline()