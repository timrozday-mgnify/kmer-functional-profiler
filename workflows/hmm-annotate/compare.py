"""One chunk's passing Pfam domains against the release's own Pfam table (mgy_proteins_pfam).

Writes ``annotate.py check``'s columns, with the release's table as the truth, for ``tables.py
check-summary``: per family, (protein, family) pairs and domains (same envelope) that the table has,
that are found here too, and that are here only. Meant for the Pfam release the table was made with.
"""

import argparse

import duckdb


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--domains", required=True)
    p.add_argument("--table", required=True, help="mgy_proteins_pfam.parquet, a path or a URL")
    p.add_argument("--lo", type=int, required=True)
    p.add_argument("--hi", type=int, required=True)
    p.add_argument("--sample", type=int, default=1, help="plan's --sample")
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--out", required=True)
    a = p.parse_args()
    con = duckdb.connect()
    con.execute(f"SET threads = {a.threads}")
    con.execute(
        f"""
        CREATE TEMP TABLE here AS SELECT DISTINCT protein_id, CAST(model[3:] AS INTEGER) AS family,
            env_from, env_to FROM read_parquet('{a.domains}') WHERE passes;
        CREATE TEMP TABLE there AS SELECT DISTINCT protein_id, pfam_accession AS family,
            env_from, env_to FROM read_parquet('{a.table}')
            WHERE protein_id BETWEEN {a.lo} AND {a.hi - 1} AND protein_id % {a.sample} = 0;
        """
    )

    def counts(name: str, keys: str) -> str:
        return f"""
            SELECT family,
                count(*) FILTER (WHERE t) AS {name}_truth,
                count(*) FILTER (WHERE t AND h) AS {name}_found,
                count(*) FILTER (WHERE h AND NOT t) AS {name}_extra
            FROM (
                SELECT family, bool_or(t) AS t, bool_or(h) AS h FROM (
                    SELECT {keys}, true AS t, false AS h FROM there
                    UNION ALL SELECT {keys}, false, true FROM here
                ) GROUP BY {keys}
            ) GROUP BY family"""

    con.execute(
        f"""
        COPY (
            SELECT format('PF{{:05d}}', family) AS model, pairs_truth, pairs_found, pairs_extra,
                domains_truth, domains_found, domains_extra
            FROM ({counts("pairs", "family, protein_id")})
            JOIN ({counts("domains", "family, protein_id, env_from, env_to")}) USING (family)
            ORDER BY model
        ) TO '{a.out}' (FORMAT csv, DELIMITER '\t', HEADER)
        """
    )


if __name__ == "__main__":
    main()
