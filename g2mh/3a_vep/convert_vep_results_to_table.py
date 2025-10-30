import argparse
from cyvcf2 import VCF
import polars as pl
import sys
from collections import defaultdict
from itertools import count

csq_format = ("Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|ALLELE_NUM|DISTANCE|STRAND|FLAGS|MINIMISED|SYMBOL_SOURCE|HGNC_ID|CANONICAL|SOURCE|pLI_gene_value|MPC_score|LOEUF|am_class|am_pathogenicity|gnomADe_V4|gnomADe_V4_AF|gnomADe_V4_nfe|gnomADe_V4_nfe_AF_nfe"
).split("|")
csq_dict = {i: f for i, f in enumerate(csq_format)}

lof_consequences = {
    "transcript_ablation", "splice_acceptor_variant", "splice_donor_variant",
    "stop_gained", "frameshift_variant", "stop_lost", "start_lost"
}
consequences = lof_consequences.union({"missense_variant"})

def get_gt(sid, sample_idx, rec):
    """Return genotype string for sid or '' if missing/not in vcf."""
    if sid == "0" or sid not in sample_idx:
        return ""
    g = rec.genotypes[ sample_idx[sid] ]
    return f"{g[0]}/{g[1]}"

def process_annotation(csq_annotation: str) -> dict:
    parts = csq_annotation.split("|")
    return { csq_dict[i]: val for i, val in enumerate(parts) }

def main():
    p = argparse.ArgumentParser(
        description="Make per‑sample variant table with child/parent genotypes + VEP CSQ"
    )
    p.add_argument("--vcf_genot", required=True,
                   help="bgzip‑compressed & indexed VCF of genotypes")
    p.add_argument("--vcf_annot", required=True,
                   help="bgzip‑compressed & indexed VCF of VEP annotations")
    p.add_argument("--psam", required=True,
                   help="TSV (no header): family_id, sample_id, father_id, mother_id, sex, phen")
    p.add_argument("--regions", nargs="+", default=None,
                   help="optional regions (e.g. chr1:10000-20000)")
    p.add_argument("--regions_file", default=None,
                   help="BED-like file with chrom  start  end  [...]  (0-based)")
    p.add_argument("--out", default="variant_matrix.csv",
                   help="output CSV filename")
    args = p.parse_args()

    # Load pedigree
    ped = pl.read_csv(args.psam, separator="\t", has_header=True)
    ped.columns = ["family_id","sample_id","father_id","mother_id","sex","phen"]
    ped_map = {
        str(e["sample_id"]): {
            "family_id": str(e["family_id"]),
            "father":    str(e["father_id"]),
            "mother":    str(e["mother_id"]),
            "sex":       str(e["sex"]),
            "phen":      str(e["phen"])
        }
        for e in ped.to_dicts()
    }
    
    #print(ped_map)

    children_by_family = defaultdict(list)
    for sid, info in ped_map.items():
        if info["father"] != "0" or info["mother"] != "0":
            children_by_family[ info["family_id"] ].append(sid)


    # Open VCFs
    vcf       = VCF(args.vcf_genot, strict_gt=False)
    vcf_annot = VCF(args.vcf_annot, strict_gt=False)
    samples   = vcf.samples
    sample_idx = {s: i for i, s in enumerate(samples)}

    rows = []
    regions = args.regions or [None]
    # -------------------------------
    # decide which regions to use
    # -------------------------------
    if args.regions_file:                 # chunk file from earlier step
        regions = []
        with open(args.regions_file) as fh:
            for ln in fh:
                if ln.startswith("#") or not ln.strip():
                    continue
                chrom, start, end, *_ = ln.rstrip().split("\t")
                # BED is 0-based; cyvcf2 expects 1-based
                regions.append(f"{chrom}:{int(start)+1}-{end}")
    elif args.regions:
        regions = args.regions
    else:
        regions = [None]                  # whole genome

    seen = set()  # will hold (family_id, chrom, start, ref, alt)

    for region in regions:
        ann_iter = vcf_annot(region) if region else vcf_annot
        gen_iter = vcf(region)       if region else vcf

        print(region,file=sys.stderr)

        for rec_annot, rec in zip(ann_iter, gen_iter):
            if "*" in rec.ALT: continue
            csq_raw = rec_annot.INFO.get("CSQ")
            if not csq_raw: continue
            #if rec.QUAL < 20: continue
            #if rec.INFO.get('DP') < 8: continue

            for ann in csq_raw.split(","):
                cd = process_annotation(ann)
                cons = cd["Consequence"]
                if cons not in consequences: continue

                alt = rec.ALT[0]
                for sample, gt in zip(samples, rec.genotypes):
                    #if sample not in children_by_family: continue
                    if gt[0] + gt[1] > 0:
                        # gather pedigree info
                        pinfo   = ped_map.get(sample, {})
                        fam     = pinfo.get("family_id")
                        dad_id  = pinfo.get("father")
                        mom_id  = pinfo.get("mother")
                        # determine role
                        role = ("offspring"
                                if dad_id != "0" or mom_id != "0"
                                else "parent")
                        if fam == None: continue
                        #print(sample)

                        key = (fam, rec.CHROM, rec.start, rec.REF, alt)
                        #print(key)
                        if key in seen: continue
                        seen.add(key)

                        # emit one row *per child* in this family
                        for child in children_by_family.get(fam, []):
                            #print(rec.CHROM)
                            rows.append({
                                "family_id":                    fam,
                                "offspring":                    child,
                                "father":                       ped_map[child]["father"],
                                "mother":                       ped_map[child]["mother"],
                                "offspring_gt":                 get_gt(child,sample_idx,rec),
                                "father_gt":                    get_gt(ped_map[child]["father"],sample_idx,rec),
                                "mother_gt":                    get_gt(ped_map[child]["mother"],sample_idx,rec),
                                "sex":                          ped_map[child]["sex"],
                                "phen":                         ped_map[child]["phen"],
                                "chrom":                        rec.CHROM,
                                "start":                        rec.start,
                                "end":                          rec.end,
                                "ref":                          rec.REF,
                                "alt":                          alt,
                                "symbol":                       cd["SYMBOL"],
                                "consequence":                  cons,
                                "canonical":                    cd["CANONICAL"],
                                "MPC_score":                    cd["MPC_score"],
                                "pLI_gene_value":               cd["pLI_gene_value"],
                                "LOEUF":                        cd["LOEUF"],
                                "gnomADe_V4_AF":                cd["gnomADe_V4_AF"],
                                "gnomADe_V4_nfe_AF_nfe":        cd["gnomADe_V4_nfe_AF_nfe"],
                            })
    print(len(rows))
    print("Done processing... making df now.",file=sys.stderr)

    df = pl.DataFrame(rows).select([
        "chrom","start","end","ref","alt",
        "family_id","offspring","father","mother","offspring_gt","father_gt","mother_gt","sex","phen",
        "symbol","consequence","canonical","MPC_score","pLI_gene_value","LOEUF",
        "gnomADe_V4_AF","gnomADe_V4_nfe_AF_nfe"
    ])

    df = df.with_columns(
        pl.when(
            # both parents are hom‑ref and child is het
            (pl.col("father_gt") == "0/0")
          & (pl.col("mother_gt") == "0/0")
          & (pl.col("offspring_gt") == "0/1")
        )
        .then(1)
        .when(
            # at least one parent is het
            (pl.col("father_gt").is_in(["0/1", "1/0"]))
          | (pl.col("mother_gt").is_in(["0/1", "1/0"]))
        )
        .then(0)
        .otherwise(None)  # covers missing parents or other combos
        .alias("denovo")
    )

    df.write_csv(args.out)
    print(f"Wrote {len(rows)} rows to {args.out}")

if __name__ == "__main__":
    main()
