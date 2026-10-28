#!/usr/bin/env bash
# Pre-flight check of the PacBio Kinnex BAMs listed in pacbio_sra_metadata.tsv:
#   - file exists and passes samtools quickcheck
#   - @RG SM tag matches the sample in library_ID (SC_<sample>_t<min>_PB)
#   - writes md5sums.txt (for verifying the transfer / SRA load)
# Usage: ./check_bams.sh [bam_dir] [samtools]
set -euo pipefail

BAMDIR=${1:-/atium/Data/projects/sheridansinglecell/pacbio/261006_sherbear_pacbio_kinnex_bams}
SAMTOOLS=${2:-samtools}
HERE=$(cd "$(dirname "$0")" && pwd)
META=$HERE/pacbio_sra_metadata.tsv

fail=0
: > "$HERE/md5sums.txt"
# columns: 3 = library_ID, 13 = filename
while IFS=$'\t' read -r lib fn; do
    f=$BAMDIR/$fn
    sample=$(echo "$lib" | cut -d_ -f2)
    if [[ ! -s $f ]]; then echo "MISSING  $fn"; fail=1; continue; fi
    # -u: unaligned BAMs have no @SQ lines
    if ! "$SAMTOOLS" quickcheck -u "$f"; then echo "CORRUPT  $fn"; fail=1; continue; fi
    sm=$("$SAMTOOLS" view -H "$f" | grep '^@RG' | grep -o 'SM:[^[:space:]]*' | sort -u | cut -d: -f2)
    if [[ $sm != "$sample" ]]; then echo "SM MISMATCH  $fn  SM=$sm  library_ID=$lib"; fail=1; continue; fi
    echo "OK  $fn  ($sm)"
    (cd "$BAMDIR" && md5sum "$fn") >> "$HERE/md5sums.txt"
done < <(tail -n +2 "$META" | cut -f3,13)
exit $fail
