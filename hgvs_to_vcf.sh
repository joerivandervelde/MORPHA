#!/usr/bin/env bash
set -euo pipefail

# Script assumes normal operation of e.g. genename2uniquehaplo.sh with all resources in place
#
# Transcript/gene matching strategy:
# 1. Transcript target: check MANE GFF (exact -> unversioned -> versioned regex), then full RefSeq GFF3 (with RefSeq->UCSC chr mapping).
# 2. Gene target (or transcript fallback to gene): resolve via MANE GFF using resolve_gene_mane_transcript_id.

if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <TRANSCRIPT_OR_GENE:c.variant> [more variants...]" >&2
    echo "Examples:" >&2
    echo "  $0 NM_001009944.3:c.11137dup" >&2
    echo "  $0 NM_000092.5:c.1321_1369+3del" >&2
    echo "  $0 NM_080632.2:c.263_263+5dup" >&2
    echo "  $0 NM_021023.6:c.614-4del" >&2
    echo "  $0 NM_000091.5:c.5010_*14del" >&2
    echo "  $0 NM_014714.4:c.c.3350T>C" >&2
    echo "  $0 APOL1:c.1024A>G" >&2
    exit 1
fi

DATADIR="data"
TMPDIR="tmp/hgvs_bash"
mkdir -p "${DATADIR}" "${TMPDIR}"

GRCh38REFFA="GRCh38_full_analysis_set_plus_decoy_hla.fa"
GRCh38REFFAINDEX="${GRCh38REFFA}.fai"
GRCh38MANEGFFGZ="MANE.GRCh38.v1.5.refseq_genomic.gff.gz"
GRCh38REFSEQGFFGZ="GCF_000001405.40_GRCh38.p14_genomic.gff.gz"
GRCh38ASSEMBLYREPORT="GCF_000001405.40_GRCh38.p14_assembly_report.txt"

FA_PATH="${DATADIR}/${GRCh38REFFA}"
FAI_PATH="${DATADIR}/${GRCh38REFFAINDEX}"
MANE_GFF_GZ_PATH="${DATADIR}/${GRCh38MANEGFFGZ}"
REFSEQ_GFF_GZ_PATH="${DATADIR}/${GRCh38REFSEQGFFGZ}"
ASSEMBLY_REPORT_PATH="${DATADIR}/${GRCh38ASSEMBLYREPORT}"

[ -s "${FA_PATH}" ] || wget -q -P "${DATADIR}" "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/${GRCh38REFFA}"
[ -s "${FAI_PATH}" ] || wget -q -P "${DATADIR}" "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/${GRCh38REFFAINDEX}"
[ -s "${MANE_GFF_GZ_PATH}" ] || wget -q -P "${DATADIR}" "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.5/${GRCh38MANEGFFGZ}"
[ -s "${REFSEQ_GFF_GZ_PATH}" ] || wget -q -c -O "${REFSEQ_GFF_GZ_PATH}" "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/${GRCh38REFSEQGFFGZ}"
[ -s "${ASSEMBLY_REPORT_PATH}" ] || wget -q -P "${DATADIR}" \
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/${GRCh38ASSEMBLYREPORT}"

fa_slice() {
    local chrom="$1" start="$2" end="$3"
    samtools faidx "${FA_PATH}" "${chrom}:${start}-${end}" | awk '!/^>/ { printf "%s", toupper($0) }'
}

find_gff_transcript_id() {
    local gff_path="$1"
    local query_tx="$2"
    zgrep -F "${query_tx}" "${gff_path}" | awk -F'\t' -v q="$query_tx" '
    $0 ~ /^#/ { next }
    ($3 == "mRNA" || $3 == "transcript") {
        n = split($9, a, ";")
        tx_id = ""; name_val = ""; id_val = ""
        for (i = 1; i <= n; i++) {
            if (a[i] ~ /^transcript_id=/) {
                tx_id = substr(a[i], 15)
            } else if (a[i] ~ /^Name=/) {
                name_val = substr(a[i], 6)
            } else if (a[i] ~ /^ID=/) {
                id_val = substr(a[i], 4)
            }
        }
        if (tx_id == q) { print tx_id; exit }
        if (name_val == q) { print name_val; exit }
        if (id_val == q) { print id_val; exit }

        if (tx_id ~ ("^" q "\\.[0-9]+$")) { print tx_id; exit }
        if (name_val ~ ("^" q "\\.[0-9]+$")) { print name_val; exit }
        if (id_val ~ ("^" q "\\.[0-9]+$")) { print id_val; exit }
    }'
}

# Always resolve gene symbols against the MANE GFF (which already uses UCSC chr names)
resolve_gene_mane_transcript_id() {
    local gene="$1"
    zgrep -F "${gene}" "${MANE_GFF_GZ_PATH}" | awk -F'\t' -v gene="$gene" '
    $0 ~ /^#/ { next }
    $3 == "gene" && ($9 ~ ("gene=" gene "(;|$)") || $9 ~ ("Name=" gene "(;|$)")) {
        n = split($9, a, ";")
        for (i = 1; i <= n; i++) {
            if (a[i] ~ /^ID=/) {
                gene_ids[substr(a[i], 4)] = 1
            }
        }
    }
    ($3 == "mRNA" || $3 == "transcript") {
        n = split($9, a, ";")
        parent_match = 0; is_mane = 0; tx_id = ""; name_val = ""; id_val = ""
        for (i = 1; i <= n; i++) {
            if (a[i] ~ /^Parent=/) {
                m = split(substr(a[i], 8), p, ",")
                for (j = 1; j <= m; j++) {
                    if (p[j] in gene_ids) parent_match = 1
                }
            } else if (a[i] ~ /^transcript_id=/) {
                tx_id = substr(a[i], 15)
            } else if (a[i] ~ /^Name=/) {
                name_val = substr(a[i], 6)
            } else if (a[i] ~ /^ID=/) {
                id_val = substr(a[i], 4)
            }
        }
        if ($9 ~ /MANE/) {
            is_mane = 1
        }
        if ((parent_match || $9 ~ ("gene=" gene "(;|$)")) && is_mane) {
            if (tx_id != "") { print tx_id; exit }
            if (name_val != "") { print name_val; exit }
            if (id_val != "") { print id_val; exit }
        }
    }'
}

# Returns "<resolved_tx>\t<gff_source>" where gff_source is "MANE" or "REFSEQ"
resolve_transcript_source() {
    local target="$1"
    if [[ "${target}" =~ ^(NM_|NR_|XM_|XR_|ENST) ]]; then
        local resolved unversioned

        # 1. Exact match in MANE GFF
        resolved=$(find_gff_transcript_id "${MANE_GFF_GZ_PATH}" "${target}")
        if [ -n "${resolved}" ]; then
            printf "%s\tMANE\n" "${resolved}"
            return 0
        fi

        # 2. Unversioned prefix & any version match in MANE GFF
        unversioned="${target%.*}"
        if [ "${unversioned}" != "${target}" ]; then
            resolved=$(find_gff_transcript_id "${MANE_GFF_GZ_PATH}" "${unversioned}")
            if [ -n "${resolved}" ]; then
                printf "%s\tMANE\n" "${resolved}"
                return 0
            fi
        fi

        # 3. Exact match in full RefSeq GFF3
        resolved=$(find_gff_transcript_id "${REFSEQ_GFF_GZ_PATH}" "${target}")
        if [ -n "${resolved}" ]; then
            printf "%s\tREFSEQ\n" "${resolved}"
            return 0
        fi

        # 4. Unversioned prefix & any version match in full RefSeq GFF3
        if [ "${unversioned}" != "${target}" ]; then
            resolved=$(find_gff_transcript_id "${REFSEQ_GFF_GZ_PATH}" "${unversioned}")
            if [ -n "${resolved}" ]; then
                printf "%s\tREFSEQ\n" "${resolved}"
                return 0
            fi
        fi

        # 5. Gene MANE fallback: look up gene symbol in MANE first, then full RefSeq, and resolve gene's MANE transcript
        local gene_sym
        gene_sym=$(zgrep -F "${unversioned}" "${MANE_GFF_GZ_PATH}" | awk -F'\t' -v q="$unversioned" '
        $0 ~ /^#/ { next }
        $0 ~ q {
            n = split($9, a, ";")
            for (i = 1; i <= n; i++) {
                if (a[i] ~ /^gene=/) {
                    print substr(a[i], 6)
                    exit
                }
            }
        }')
        if [ -z "${gene_sym}" ]; then
            gene_sym=$(zgrep -F "${unversioned}" "${REFSEQ_GFF_GZ_PATH}" | awk -F'\t' -v q="$unversioned" '
            $0 ~ /^#/ { next }
            $0 ~ q {
                n = split($9, a, ";")
                for (i = 1; i <= n; i++) {
                    if (a[i] ~ /^gene=/) {
                        print substr(a[i], 6)
                        exit
                    }
                }
            }')
        fi
        if [ -n "${gene_sym}" ]; then
            resolved=$(resolve_gene_mane_transcript_id "${gene_sym}")
            if [ -n "${resolved}" ]; then
                printf "%s\tMANE\n" "${resolved}"
                return 0
            fi
        fi
        return 1
    fi

    # Gene symbol input: resolve via MANE GFF
    local mane_tx
    mane_tx=$(resolve_gene_mane_transcript_id "${target}")
    if [ -n "${mane_tx}" ]; then
        printf "%s\tMANE\n" "${mane_tx}"
        return 0
    fi
    return 1
}

normalize_chrom_to_ucsc() {
    local use_refseq_map="$1"
    local rows="$2"
    if [ "${use_refseq_map}" = "1" ]; then
        awk -F'\t' -v OFS='\t' '
            FNR==NR {
                if ($0 !~ /^#/) {
                    gsub(/\r/, "", $7)
                    gsub(/\r/, "", $10)
                    if ($7 != "" && $7 != "na" && $10 != "" && $10 != "na") {
                        refseq_to_ucsc[$7] = $10
                    }
                }
                next
            }
            {
                gsub(/\r/, "", $1)
                if ($1 in refseq_to_ucsc) {
                    $1 = refseq_to_ucsc[$1]
                }
                print
            }
        ' "${ASSEMBLY_REPORT_PATH}" <(printf "%s\n" "${rows}")
    else
        printf "%s\n" "${rows}" | awk -F'\t' -v OFS='\t' '{ gsub(/\r/, "", $1); print }'
    fi
}

RAW_VCF="${TMPDIR}/raw.vcf"
{
    printf "##fileformat=VCFv4.2\n"
    awk -F'\t' '{ printf "##contig=<ID=%s,length=%s>\n", $1, $2 }' "${FAI_PATH}"
    printf "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
} > "${RAW_VCF}"

for hgvs in "$@"; do
    target="${hgvs%%:*}"
    change="${hgvs#*:}"

    # Fix accidental duplicate leading "c." (e.g., c.c.3350T>C)
    if [[ "${change}" =~ ^c\.c\. ]]; then
        change="c.${change#c.c.}"
    fi

    resolved_pair=$(resolve_transcript_source "${target}") || resolved_pair=""
    [ -n "${resolved_pair}" ] || { echo "Error: Could not resolve transcript for '${target}' in ${hgvs}" >&2; exit 1; }

    tx="${resolved_pair%%$'\t'*}"
    gff_src="${resolved_pair#*$'\t'}"
    out_id="${tx}:${change}"

    if [ "${gff_src}" = "MANE" ]; then
        active_gff="${MANE_GFF_GZ_PATH}"
        use_refseq_map="0"
    else
        active_gff="${REFSEQ_GFF_GZ_PATH}"
        use_refseq_map="1"
    fi

    raw_cds_rows=$(
        zgrep -F "${tx}" "${active_gff}" | awk -F'\t' -v tx="$tx" '
        BEGIN { OFS="\t" }
        $0 ~ /^#/ { next }
        ($3 == "mRNA" || $3 == "transcript") && ($9 ~ ("ID=" tx "(;|$)") || $9 ~ ("Name=" tx "(;|$)") || $9 ~ ("transcript_id=" tx "(;|$)")) {
            n = split($9, a, ";")
            for (i = 1; i <= n; i++) {
                if (a[i] ~ /^ID=/) {
                    tx_ids[substr(a[i], 4)] = 1
                }
            }
        }
        $3 == "CDS" {
            match_tx = ($9 ~ ("transcript_id=" tx "(;|$)"))
            if (!match_tx) {
                n = split($9, a, ";")
                for (i = 1; i <= n; i++) {
                    if (a[i] ~ /^Parent=/) {
                        m = split(substr(a[i], 8), p, ",")
                        for (j = 1; j <= m; j++) {
                            if (p[j] in tx_ids) match_tx = 1
                        }
                    }
                }
            }
            if (match_tx) print $1, $4, $5, $7
        }'
    )

    [ -n "${raw_cds_rows}" ] || { echo "Error: No CDS found for ${tx} (${hgvs})" >&2; exit 1; }

    cds_rows=$(normalize_chrom_to_ucsc "${use_refseq_map}" "${raw_cds_rows}")

    strand=$(printf "%s\n" "${cds_rows}" | awk 'NR==1 { print $4 }')
    if [ "${strand}" = "+" ]; then
        sorted_cds=$(printf "%s\n" "${cds_rows}" | sort -k2,2n)
    elif [ "${strand}" = "-" ]; then
        sorted_cds=$(printf "%s\n" "${cds_rows}" | sort -k3,3nr)
    else
        echo "Error: Unsupported strand ${strand}" >&2
        exit 1
    fi

    cds_total_len=$(printf "%s\n" "${sorted_cds}" | awk '{ sum += ($3 - $2 + 1) } END { print sum + 0 }')
    last_cds_row=$(printf "%s\n" "${sorted_cds}" | tail -n 1)

    c_pt_to_g() {
        local pt_expr="$1"
        local offset_str="${2:-0}"

        if [[ "${pt_expr}" =~ ^\*([0-9]+)$ ]]; then
            local utr_off="${BASH_REMATCH[1]}"
            printf "%s\n" "${last_cds_row}" | awk -v base_len="$cds_total_len" -v utr_off="$utr_off" -v ioff="$offset_str" -v strand="$strand" '
            {
                chrom = $1; s = $2 + 0; e = $3 + 0
                step = utr_off + ioff
                if (strand == "+") {
                    print chrom ":" (e + step)
                } else {
                    print chrom ":" (s - step)
                }
            }'
            return 0
        fi

        if [[ ! "${pt_expr}" =~ ^[0-9]+$ ]]; then
            return 1
        fi

        printf "%s\n" "${sorted_cds}" | awk -v c="$pt_expr" -v ioff="$offset_str" -v strand="$strand" '
        BEGIN { acc = 0; found = 0 }
        {
            chrom = $1; s = $2 + 0; e = $3 + 0; len = e - s + 1
            if (acc + len >= c) {
                off = c - acc - 1
                gpos = (strand == "+") ? (s + off + ioff) : (e - off - ioff)
                print chrom ":" gpos
                found = 1
                exit
            }
            acc += len
        }
        END { if (!found) exit 1 }'
    }

    # Parse coding or intronic SNV: c.123A>G, c.817-1G>T, c.1150+1G>A, c.*14A>G
    if [[ "${change}" =~ ^c\.(\*[0-9]+|[0-9]+)([+-][0-9]+)?([ACGT])\>([ACGT])$ ]]; then
        p1="${BASH_REMATCH[1]}"
        off1="${BASH_REMATCH[2]:-0}"
        ref_c="${BASH_REMATCH[3]}"
        alt_c="${BASH_REMATCH[4]}"

        g1=$(c_pt_to_g "$p1" "$off1") || { echo "Error: Coordinate out of range for ${hgvs}" >&2; exit 1; }
        chrom="${g1%%:*}"
        pos="${g1#*:}"

        if [ "${strand}" = "-" ]; then
            ref_g=$(printf "%s" "$ref_c" | tr 'ACGT' 'TGCA')
            alt_g=$(printf "%s" "$alt_c" | tr 'ACGT' 'TGCA')
        else
            ref_g="$ref_c"
            alt_g="$alt_c"
        fi

        fa_ref=$(fa_slice "$chrom" "$pos" "$pos")
        if [ -z "${fa_ref}" ]; then
            echo "Error: Empty reference base at ${chrom}:${pos} for ${hgvs}" >&2
            exit 1
        fi
        if [ "${fa_ref}" != "${ref_g}" ]; then
            echo "Error: Reference mismatch at ${chrom}:${pos} for ${hgvs} (expected ${ref_g} from ${tx} strand ${strand}, got ${fa_ref} in FASTA)" >&2
            exit 1
        fi

        printf "%s\t%s\t%s\t%s\t%s\t.\tPASS\t.\n" "$chrom" "$pos" "$out_id" "$fa_ref" "$alt_g" >> "${RAW_VCF}"

    elif [[ "${change}" =~ ^c\.(\*[0-9]+|[0-9]+)([+-][0-9]+)?(_(\*[0-9]+|[0-9]+)([+-][0-9]+)?)?dup\[?[ACGT]*\]?$ ]]; then
        p1="${BASH_REMATCH[1]}"
        off1="${BASH_REMATCH[2]:-0}"
        p2="${BASH_REMATCH[4]:-$p1}"
        off2="${BASH_REMATCH[5]:-$off1}"

        g1_str=$(c_pt_to_g "$p1" "$off1") || { echo "Error: Start coordinate out of range in ${hgvs}" >&2; exit 1; }
        g2_str=$(c_pt_to_g "$p2" "$off2") || { echo "Error: End coordinate out of range in ${hgvs}" >&2; exit 1; }
        chrom="${g1_str%%:*}"
        g1="${g1_str#*:}"
        g2="${g2_str#*:}"

        if (( g1 < g2 )); then min_g="$g1"; max_g="$g2"; else min_g="$g2"; max_g="$g1"; fi

        anchor_pos=$(( min_g - 1 ))
        anchor_base=$(fa_slice "$chrom" "$anchor_pos" "$anchor_pos")
        dup_seq=$(fa_slice "$chrom" "$min_g" "$max_g")
        if [ -z "${anchor_base}" ] || [ -z "${dup_seq}" ]; then
            echo "Error: Empty sequence slice at ${chrom}:${anchor_pos}-${max_g} for ${hgvs}" >&2
            exit 1
        fi
        printf "%s\t%s\t%s\t%s\t%s\t.\tPASS\t.\n" "$chrom" "$anchor_pos" "$out_id" "$anchor_base" "${anchor_base}${dup_seq}" >> "${RAW_VCF}"

    elif [[ "${change}" =~ ^c\.(\*[0-9]+|[0-9]+)([+-][0-9]+)?(_(\*[0-9]+|[0-9]+)([+-][0-9]+)?)?del\[?[ACGT]*\]?$ ]]; then
        p1="${BASH_REMATCH[1]}"
        off1="${BASH_REMATCH[2]:-0}"
        p2="${BASH_REMATCH[4]:-$p1}"
        off2="${BASH_REMATCH[5]:-$off1}"

        g1_str=$(c_pt_to_g "$p1" "$off1") || { echo "Error: Start coordinate out of range in ${hgvs}" >&2; exit 1; }
        g2_str=$(c_pt_to_g "$p2" "$off2") || { echo "Error: End coordinate out of range in ${hgvs}" >&2; exit 1; }
        chrom="${g1_str%%:*}"
        g1="${g1_str#*:}"
        g2="${g2_str#*:}"

        if (( g1 < g2 )); then min_g="$g1"; max_g="$g2"; else min_g="$g2"; max_g="$g1"; fi

        anchor_pos=$(( min_g - 1 ))
        anchor_base=$(fa_slice "$chrom" "$anchor_pos" "$anchor_pos")
        del_seq=$(fa_slice "$chrom" "$min_g" "$max_g")
        if [ -z "${anchor_base}" ] || [ -z "${del_seq}" ]; then
            echo "Error: Empty sequence slice at ${chrom}:${anchor_pos}-${max_g} for ${hgvs}" >&2
            exit 1
        fi
        printf "%s\t%s\t%s\t%s\t%s\t.\tPASS\t.\n" "$chrom" "$anchor_pos" "$out_id" "${anchor_base}${del_seq}" "$anchor_base" >> "${RAW_VCF}"

    elif [[ "${change}" =~ ^c\.(\*[0-9]+|[0-9]+)([+-][0-9]+)?_(\*[0-9]+|[0-9]+)([+-][0-9]+)?ins\[?([ACGT]+)\]?$ ]]; then
        p1="${BASH_REMATCH[1]}"
        off1="${BASH_REMATCH[2]:-0}"
        p2="${BASH_REMATCH[3]}"
        off2="${BASH_REMATCH[4]:-0}"
        ins_c="${BASH_REMATCH[5]}"

        g1_str=$(c_pt_to_g "$p1" "$off1") || { echo "Error: Start coordinate out of range in ${hgvs}" >&2; exit 1; }
        g2_str=$(c_pt_to_g "$p2" "$off2") || { echo "Error: End coordinate out of range in ${hgvs}" >&2; exit 1; }
        chrom="${g1_str%%:*}"
        g1="${g1_str#*:}"
        g2="${g2_str#*:}"

        anchor_pos=$(( g1 < g2 ? g1 : g2 ))
        anchor_base=$(fa_slice "$chrom" "$anchor_pos" "$anchor_pos")
        if [ "${strand}" = "-" ]; then
            ins_g=$(printf "%s" "$ins_c" | tr 'ACGT' 'TGCA' | rev)
        else
            ins_g="$ins_c"
        fi
        if [ -z "${anchor_base}" ]; then
            echo "Error: Empty anchor base at ${chrom}:${anchor_pos} for ${hgvs}" >&2
            exit 1
        fi
        printf "%s\t%s\t%s\t%s\t%s\t.\tPASS\t.\n" "$chrom" "$anchor_pos" "$out_id" "$anchor_base" "${anchor_base}${ins_g}" >> "${RAW_VCF}"

    else
        echo "Error: Unsupported HGVS change '${change}' in ${hgvs}" >&2
        exit 1
    fi
done

bcftools sort "${RAW_VCF}" | bcftools norm -f "${FA_PATH}" -c w --no-version -Ov | { grep -v '^#' || true; }
