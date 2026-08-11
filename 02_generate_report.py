#!/usr/bin/env python3
"""
Build a clinical variant report from a ClinVar-annotated VCF
(output of 01_pipeline.sh, e.g. sample1.clinvar_annotated.vcf).

Usage:
    python 02_generate_report.py <annotated.vcf> <sample_name> <outdir>

Produces:
    <outdir>/<sample>_variant_report.md    (human-readable, plain text/markdown)
    <outdir>/<sample>_variant_report.html  (styled, browsable)
    <outdir>/<sample>_variants.tsv         (flat table, all variants)
"""

import sys
import gzip
import csv
from pathlib import Path
from collections import defaultdict

# ClinVar CLNSIG categories -> our report buckets
PATHOGENIC_TERMS = {"pathogenic", "likely_pathogenic", "pathogenic/likely_pathogenic"}
VUS_TERMS = {"uncertain_significance", "conflicting_interpretations_of_pathogenicity"}
BENIGN_TERMS = {"benign", "likely_benign", "benign/likely_benign"}
OTHER_TERMS = {"drug_response", "risk_factor", "association", "protective", "affects", "not_provided"}

EXPLANATIONS = {
    "pathogenic": "ClinVar submitters agree this variant causes or strongly predisposes to the "
                   "associated disease, based on clinical, functional, and/or segregation evidence.",
    "likely_pathogenic": "Evidence favors disease causation (>90% confidence per ACMG/AMP guidelines) "
                          "but does not yet meet the full threshold for 'Pathogenic'.",
    "uncertain_significance": "Available evidence is insufficient or conflicting to classify this "
                               "variant as disease-causing or benign (a 'VUS'). Not actionable on its own; "
                               "may be reclassified as more data accumulates.",
    "conflicting_interpretations_of_pathogenicity": "Different submitting laboratories have reached "
                               "different classifications for this variant; treat as uncertain pending review.",
    "likely_benign": "Evidence favors a benign (non-disease-causing) interpretation but does not fully "
                      "meet criteria for 'Benign'.",
    "benign": "This variant is common and/or has strong evidence showing it does not cause the "
              "associated disease.",
    "drug_response": "This variant is associated with a response to a specific medication "
                      "(pharmacogenomic relevance), not with disease risk directly.",
    "risk_factor": "This variant is associated with increased risk of a condition but is not "
                    "sufficient by itself to cause disease.",
    "association": "Statistical association with a trait/disease reported, typically from GWAS-type "
                    "evidence rather than direct causal proof.",
    "affects": "Evidence that the variant affects a non-disease phenotype/trait.",
    "protective": "Evidence suggests this variant reduces disease risk.",
    "not_provided": "A ClinVar record exists but no classification was submitted.",
}


def norm_clnsig(raw):
    """Normalize ClinVar CLNSIG string to lowercase, underscore-joined token(s)."""
    if not raw or raw == ".":
        return []
    # ClinVar can pack multiple/combined terms with , or |
    parts = raw.replace("|", ",").split(",")
    return [p.strip().lower() for p in parts if p.strip()]


def bucket_for(clnsig_terms):
    terms = set(clnsig_terms)
    if terms & PATHOGENIC_TERMS or any("pathogenic" in t and "likely" not in t for t in terms) or \
       any(t.startswith("likely_pathogenic") for t in terms):
        return "Pathogenic / Likely Pathogenic"
    if terms & BENIGN_TERMS or any(t.startswith("benign") or t.startswith("likely_benign") for t in terms):
        return "Benign / Likely Benign"
    if terms & VUS_TERMS or any("uncertain" in t or "conflicting" in t for t in terms):
        return "Uncertain Significance (VUS)"
    if terms & OTHER_TERMS:
        return "Other (drug response / risk factor / association)"
    if not terms:
        return "Not in ClinVar"
    return "Uncertain Significance (VUS)"


def parse_info(info_str):
    d = {}
    for kv in info_str.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            d[k] = v
        else:
            d[kv] = True
    return d


def open_vcf(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def main():
    if len(sys.argv) != 4:
        print(__doc__)
        sys.exit(1)

    vcf_path, sample, outdir = sys.argv[1], sys.argv[2], sys.argv[3]
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    variants = []
    with open_vcf(vcf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom, pos, vid, ref, alt, qual, flt, info = fields[:8]
            info_d = parse_info(info)

            clnsig_raw = info_d.get("CLNSIG", "")
            clndn_raw = info_d.get("CLNDN", "").replace("_", " ")
            gene = info_d.get("GENEINFO", "").split(":")[0] if info_d.get("GENEINFO") else "-"
            revstat = info_d.get("CLNREVSTAT", "").replace("_", " ")
            hgvs = info_d.get("CLNHGVS", "-")

            terms = norm_clnsig(clnsig_raw)
            bucket = bucket_for(terms)
            explanation_parts = [EXPLANATIONS.get(t, "") for t in terms if EXPLANATIONS.get(t)]
            explanation = " ".join(explanation_parts) if explanation_parts else \
                "No ClinVar classification available for this variant."

            variants.append({
                "chrom": chrom, "pos": pos, "id": vid, "ref": ref, "alt": alt,
                "gene": gene, "clnsig": clnsig_raw or "-", "clndn": clndn_raw or "-",
                "review_status": revstat or "-", "hgvs": hgvs,
                "bucket": bucket, "explanation": explanation,
            })

    # ---- write flat TSV ----
    tsv_path = outdir / f"{sample}_variants.tsv"
    with open(tsv_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(variants[0].keys()) if variants else
                            ["chrom","pos","id","ref","alt","gene","clnsig","clndn",
                             "review_status","hgvs","bucket","explanation"], delimiter="\t")
        w.writeheader()
        for v in variants:
            w.writerow(v)

    # ---- group by bucket ----
    grouped = defaultdict(list)
    for v in variants:
        grouped[v["bucket"]].append(v)

    bucket_order = [
        "Pathogenic / Likely Pathogenic",
        "Uncertain Significance (VUS)",
        "Other (drug response / risk factor / association)",
        "Benign / Likely Benign",
        "Not in ClinVar",
    ]

    # ---- Markdown report ----
    md_lines = [f"# Clinical Variant Report — {sample}", ""]
    md_lines.append(f"Total variants processed: **{len(variants)}**")
    for b in bucket_order:
        md_lines.append(f"- {b}: **{len(grouped.get(b, []))}**")
    md_lines.append("")

    for b in bucket_order:
        vs = grouped.get(b, [])
        if not vs:
            continue
        md_lines.append(f"## {b} ({len(vs)})")
        md_lines.append("")
        for v in vs:
            md_lines.append(f"### {v['gene']} — {v['chrom']}:{v['pos']} {v['ref']}>{v['alt']}")
            md_lines.append(f"- **HGVS**: {v['hgvs']}")
            md_lines.append(f"- **ClinVar significance**: {v['clnsig']}")
            md_lines.append(f"- **Associated condition(s)**: {v['clndn']}")
            md_lines.append(f"- **Review status**: {v['review_status']}")
            md_lines.append(f"- **Explanation**: {v['explanation']}")
            md_lines.append("")

    md_path = outdir / f"{sample}_variant_report.md"
    md_path.write_text("\n".join(md_lines))

    # ---- HTML report ----
    color = {
        "Pathogenic / Likely Pathogenic": "#c0392b",
        "Uncertain Significance (VUS)": "#e67e22",
        "Other (drug response / risk factor / association)": "#8e44ad",
        "Benign / Likely Benign": "#27ae60",
        "Not in ClinVar": "#7f8c8d",
    }
    html = [
        "<html><head><meta charset='utf-8'><title>Variant Report</title>",
        "<style>",
        "body{font-family:Arial,Helvetica,sans-serif;margin:40px;color:#222;}",
        "h1{border-bottom:3px solid #333;padding-bottom:8px;}",
        ".summary{background:#f4f4f4;padding:15px;border-radius:8px;margin-bottom:25px;}",
        ".card{border-left:6px solid #999;background:#fafafa;padding:12px 16px;margin:10px 0;border-radius:4px;}",
        ".tag{display:inline-block;color:white;padding:3px 10px;border-radius:12px;font-size:12px;margin-bottom:6px;}",
        "table{border-collapse:collapse;width:100%;margin-top:10px;}",
        "td,th{border:1px solid #ddd;padding:4px 8px;font-size:13px;text-align:left;}",
        "</style></head><body>",
        f"<h1>Clinical Variant Report — {sample}</h1>",
        "<div class='summary'><b>Summary</b><ul>",
    ]
    html.append(f"<li>Total variants processed: {len(variants)}</li>")
    for b in bucket_order:
        html.append(f"<li>{b}: {len(grouped.get(b, []))}</li>")
    html.append("</ul></div>")

    for b in bucket_order:
        vs = grouped.get(b, [])
        if not vs:
            continue
        html.append(f"<h2>{b} ({len(vs)})</h2>")
        for v in vs:
            html.append(f"<div class='card' style='border-color:{color[b]}'>")
            html.append(f"<span class='tag' style='background:{color[b]}'>{b}</span><br>")
            html.append(f"<b>{v['gene']}</b> — {v['chrom']}:{v['pos']} {v['ref']}&gt;{v['alt']}<br>")
            html.append("<table>")
            html.append(f"<tr><th>HGVS</th><td>{v['hgvs']}</td></tr>")
            html.append(f"<tr><th>ClinVar significance</th><td>{v['clnsig']}</td></tr>")
            html.append(f"<tr><th>Associated condition(s)</th><td>{v['clndn']}</td></tr>")
            html.append(f"<tr><th>Review status</th><td>{v['review_status']}</td></tr>")
            html.append("</table>")
            html.append(f"<p><i>{v['explanation']}</i></p>")
            html.append("</div>")

    html.append("</body></html>")
    html_path = outdir / f"{sample}_variant_report.html"
    html_path.write_text("\n".join(html))

    print(f"Wrote:\n  {tsv_path}\n  {md_path}\n  {html_path}")


if __name__ == "__main__":
    main()
