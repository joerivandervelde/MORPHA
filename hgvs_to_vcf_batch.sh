#!/usr/bin/env bash
set -euo pipefail

# run and sort output:
# ./hgvs_to_vcf_batch.sh variants.csv 1 | bcftools sort -Ov -o sorted.vcf

if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <input.csv> [hgvs_column_1_based]" >&2
    echo "Example: $0 variants.csv 1" >&2
    exit 1
fi

CSV_FILE="$1"
COL="${2:-1}"
HGVS_SCRIPT="${HGVS_SCRIPT:-./hgvs_to_vcf.sh}"
FAI_PATH="${FAI_PATH:-data/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai}"
SAMPLE_NAME="${SAMPLE_NAME:-SAMPLE1}"

if [ ! -s "${CSV_FILE}" ]; then
    echo "Error: CSV file '${CSV_FILE}' not found or empty." >&2
    exit 1
fi

if [ ! -x "${HGVS_SCRIPT}" ]; then
    echo "Error: HGVS script '${HGVS_SCRIPT}' not found or not executable." >&2
    exit 1
fi

if [ ! -s "${FAI_PATH}" ]; then
    mkdir -p "$(dirname "${FAI_PATH}")"
    wget -q -P "$(dirname "${FAI_PATH}")" \
        "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai"
fi

ERR_TMP=$(mktemp)
trap 'rm -f "${ERR_TMP}"' EXIT

# Stream valid VCFv4.2 header with GT FORMAT and 1 dummy sample
printf "##fileformat=VCFv4.2\n"
awk -F'\t' '{ printf "##contig=<ID=%s,length=%s>\n", $1, $2 }' "${FAI_PATH}"
printf "##FILTER=<ID=FAIL,Description=\"HGVS to VCF translation failed\">\n"
printf "##INFO=<ID=HGVS_INPUT,Number=1,Type=String,Description=\"Original HGVS input string\">\n"
printf "##INFO=<ID=ERROR_MSG,Number=1,Type=String,Description=\"Error message from HGVS translation script\">\n"
printf "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
printf "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n" "${SAMPLE_NAME}"

vcf_info_escape() {
    printf "%s" "$1" \
        | tr '\r\n\t' '   ' \
        | sed 's/  */ /g; s/^ //; s/ $//' \
        | sed 's/%/%25/g; s/ /%20/g; s/;/%3B/g; s/,/%2C/g; s/=/%3D/g'
}

vcf_id_clean() {
    printf "%s" "$1" | tr -d '\r\n' | tr ' \t' '__'
}

while IFS= read -r raw_line || [ -n "${raw_line}" ]; do
    clean_line="${raw_line%$'\r'}"

    hgvs=$(printf "%s\n" "${clean_line}" | awk -F',' -v col="${COL}" '
    {
        val = $col
        gsub(/^[[:space:]"]+|[[:space:]"]+$/, "", val)
        print val
    }')

    display_var="${hgvs:-<EMPTY>}"
    vcf_line=""
    status=0
    : > "${ERR_TMP}"

    if [ -n "${hgvs}" ]; then
        vcf_line=$("${HGVS_SCRIPT}" "${hgvs}" 2>"${ERR_TMP}") || status=$?
    else
        status=1
        printf "Empty HGVS value in column %s" "${COL}" > "${ERR_TMP}"
    fi

    line_count=$(printf "%s" "${vcf_line}" | awk 'END { print NR }')
    if [ "${status}" -ne 0 ] || [ -z "${vcf_line}" ] || [ "${line_count}" -ne 1 ]; then
        err_one_line=$(tr '\r\n' '  ' < "${ERR_TMP}" | sed 's/  */ /g; s/^ //; s/ $//')
        if [ -n "${err_one_line}" ]; then
            full_err="Error [${display_var}] (exit ${status}): ${err_one_line}"
        else
            full_err="Error [${display_var}] (exit ${status}): HGVS translation failed"
        fi

        escaped_err=$(vcf_info_escape "${full_err}")
        escaped_input=$(vcf_info_escape "${display_var}")
        id_field=$(vcf_id_clean "${hgvs:-MISSING_HGVS}")

        printf "chr3\t169764547\t%s\tT\tC\t.\tFAIL\tHGVS_INPUT=%s;ERROR_MSG=%s\tGT\t1/1\n" \
            "${id_field}" "${escaped_input}" "${escaped_err}"
    else
        printf "%s\tGT\t1/1\n" "${vcf_line}"
    fi
done < "${CSV_FILE}"
