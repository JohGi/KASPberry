rule compute_mash_distances:
    conda:
        "../envs/mash.yaml"
    input:
        CLEAN_FASTAS
    output:
        matrix=MASH_MATRIX
    threads: 1
    log:
        stdout=LOG_DIR / "compute_mash_distances" / "compute_mash_distances.stdout",
        stderr=LOG_DIR / "compute_mash_distances" / "compute_mash_distances.stderr"
    benchmark:
        BENCHMARK_DIR / "compute_mash_distances.tsv"
    shell:
        r"""
        mkdir -p "{MASH_DISTANCES_DIR}" "$(dirname "{log.stdout}")"

        tmpdir="$(mktemp -d "{MASH_DISTANCES_DIR}/mash.XXXXXX")"
        trap 'rm -rf "$tmpdir"' EXIT

        mash sketch \
            -k 21 \
            -s 1000 \
            -o "$tmpdir/regions" \
            {input} \
            1> "{log.stdout}" \
            2> "{log.stderr}"

        mash dist \
            -t \
            "$tmpdir/regions.msh" \
            "$tmpdir/regions.msh" \
            > "{output.matrix}" \
            2>> "{log.stderr}"
        """
