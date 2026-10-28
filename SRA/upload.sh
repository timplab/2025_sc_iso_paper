#!/usr/bin/env bash
# Upload the PacBio Kinnex BAMs to an SRA preload folder with Aspera.
# Get KEY and FOLDER from the SRA Submission Portal:
#   https://submit.ncbi.nlm.nih.gov/subs/sra/ -> "My Preload Folder" -> "Aspera command line upload instructions"
# Usage: ./upload.sh <aspera_key_file> <uploads/your_account_folder>
# Run inside tmux/screen: ~144 GB.
set -euo pipefail

KEY=$1
FOLDER=$2
BAMDIR=${BAMDIR:-/atium/Data/projects/sheridansinglecell/pacbio/261006_sherbear_pacbio_kinnex_bams}
SUBDIR=pacbio_kinnex
HERE=$(cd "$(dirname "$0")" && pwd)

# file list straight from the metadata sheet, so only submitted files go up
tail -n +2 "$HERE/pacbio_sra_metadata.tsv" | cut -f13 | sed "s|^|$BAMDIR/|" > "$HERE/upload_files.txt"

# -k1 resumes partial transfers if this has to be re-run
ascp -i "$KEY" -QT -l500m -k1 --file-list="$HERE/upload_files.txt" \
    "subasp@upload.ncbi.nlm.nih.gov:$FOLDER/$SUBDIR"
