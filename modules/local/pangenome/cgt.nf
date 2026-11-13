// Ported from EBI-Metagenomics/genomes-catalogue-pipeline modules/cgt.nf.
// The genome completeness comes from the catalogue metadata, not from a CheckM2 table.
// CGT runs only for the species that PANGENOME_FIX corrected. So it also emits the corrected
// files again (emit: pangenome): one publish rule then writes all the files of a corrected
// species, and no folder is made for the other species.
process CGT {
    tag "$meta.id"
    label 'process_medium'

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/cgt:0.1.1--hab7d0fd_0'
        : 'quay.io/biocontainers/cgt:0.1.1--hab7d0fd_0'}"

    input:
    tuple val(meta), path(csv), path(rtab), path(fasta), path(renamed)
    path(metadata)

    output:
    tuple val(meta), path("gene_prevalence_corrected.txt.gz"), emit: cgt
    tuple val(meta), path(csv, includeInputs: true), path(rtab, includeInputs: true), path(fasta, includeInputs: true), path(renamed, includeInputs: true), emit: pangenome
    path "versions.yml"                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def VERSION = '0.1.1' // WARN: cgt_bacpop has no option to print the version. This is the container version
    """
    gzip -cdf ${rtab} > gene_presence_absence.Rtab

    # Completeness table (genome, completeness, contamination) for the genomes of this species only
    head -1 gene_presence_absence.Rtab | tr '\\t' '\\n' | tail -n +2 > rtab_genomes.txt
    awk -F'\\t' 'FNR==NR {g[\$1]=1; next} FNR==1 {print "genome\\tcompleteness\\tcontamination"; next} \$1 in g {print \$1 "\\t" \$7 "\\t" \$8}' \\
        rtab_genomes.txt ${metadata} > completeness.tsv

    cgt_bacpop \\
        ${args} \\
        --completeness-column 2 \\
        --output-file cgt.txt \\
        completeness.tsv \\
        gene_presence_absence.Rtab > cgt.log 2>&1

    # Keep the core and rare thresholds as comment lines, as in the catalogue.
    # cgt_bacpop writes the genes in a random order; sort them so the same input gives the same file.
    {
        grep -E "^(Core|Rare) threshold:" cgt.log | sed 's/^/# /'
        head -1 cgt.txt
        tail -n +2 cgt.txt | LC_ALL=C sort
    } | gzip -n > gene_prevalence_corrected.txt.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cgt: $VERSION
    END_VERSIONS
    """

    stub:
    def VERSION = '0.1.1'
    """
    echo "" | gzip -n > gene_prevalence_corrected.txt.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cgt: $VERSION
    END_VERSIONS
    """
}
