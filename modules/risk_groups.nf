process ANNOTATE_BLAST_RISK_GROUPS {
  tag "${sample_id}, ${query_class}"
  label "low"

  input:
  tuple val(sample_id), val(query_class), path(blast_tsv)
  path risk_group_lookup
  val taxonomy_dir

  output:
  tuple val(sample_id), val(query_class), path("${sample_id}.${query_class}.blast.with_risk_groups.tsv")

  script:
  """
  annotate_risk_groups.py blast \
      --input ${blast_tsv} \
      --lookup ${risk_group_lookup} \
      --taxonomy-dir ${taxonomy_dir} \
      --output ${sample_id}.${query_class}.blast.with_risk_groups.tsv
  """
}
