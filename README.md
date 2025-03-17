# Overview

This repo contains the source code and structure to comprehensively describe
the analysis. Data are downloaded automatically from SRA, and config files are
written to work with these data. Other genomes or assemblies should be able to
be used with appropriate edits.

This repo is based off of [lcdb-wf
v1.10.2](https://github.com/lcdb/lcdb-wf/releases/tag/v1.10.2), with
substantial modifications for RIL-seq. The complete set of
[Snakemake](https://snakemake.readthedocs.io/en/stable/) workflows,
configuration, supporting code, and downstream RMarkdown files are included.

## Requirements

- [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) must be installed
- this code repository must be downloaded and unpacked. Below, we refer to the
  unpacked directory as `$WORKDIR`. Unless otherwise specified, paths provided
  are assumed to be relative to `$WORKDIR`.

## Create conda environments

Create the conda environments in the top level of `$WORKDIR`:

```bash
conda env create -p env-rilseq --file env-rilseq.yml
conda env create -p env-r --file env-r.yml
```

## Run annotation workflow

**Optional.** The files created by this workflow are included in this repo,
under `workflows/annotation/data`. However, for other genomes or assemblies
this workflow would need to be run.

Optionally run the annotation workflow. This will create the annotation file
formatted for RILseq analysis. The annotation file will contain the CDS
features, 5UTR, 3UTR, AS (antisense) and IGR (intergenic regions) that are
needed for RILseq analysis. The `workflows/annotation/config/config.yaml`
configures this. From the `workflows/annotation` directory, activate the
environment `env-rilseq` and run the workflow. Either use a [snakemake
profile](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles)
appropriate for your cluster, or run on a single machine (specifying the number
of cores; here we use the placeholder `$CORES`)

```bash
# in the workflows/annotation directory
conda activate ../../env-rilseq  # activate environment
snakemake --cores $CORES  # run workflow
```

The main output is the RILseq-formatted annotation file,
`workflows/annotation/data/ecoli_default.patch1.gff`. This is already
configured to be used in other workflows. For example, for the
`workflows/rnaseq/config/config.yaml` and `workflows/rilseq/config/config.yaml`
files, the entry `references:ecoli:default:annotation:url` points to this
path by `file://../annotation/data/ecoli_default.patch1.gff`.

## Configure RIL-seq workflows

There are two RIL-seq workflows, `rilseq-2016` and `rilseq-2020`. 

If necessary, edit the sample tables in each of the dataset directories (the
`config/sampletable.tsv` files, for example
`workflows/rilseq-2016/config/sampletable.tsv`). These are currently configured
to pull down data automatically from SRA, but if running on different genomes
or assemblies these would need additional configuration (see [lcdb-wf
docs](https://lcdb.github.io/lcdb-wf/sampletable.html#rna-seq-sample-table) for
more).

Circos plots, generated by the rilseq workflow, need a configuration file that
matches the reference genome and annotation. Here we used the [default E. coli
file provided in the RILseq Github repo](https://github.com/asafpr/RILseq/blob/master/data/E_coli_K12/ver2/draw_links.conf). For
other organisms or assemblies, edit
`workflows/rilseq/RILseq-master/data/E_coli_K12/ver2/draw_links.conf` as
needed.

## Generate trimmed fastq files

For each RIL-seq workflow (`workflows/rilseq-2016` and
`workflows/rilseq-2020`), move to the corresponding directory and run the
`Snakefile-rnaseq` workflow. This will create cutadapt-trimmed FASTQ files, as
well as run extensive QC that is summarized in the respective
`data/rnaseq_aggregation/multiqc.html` file. However, subsequent steps below
only need the trimmed FASTQ files, so this workflow can be optionally run in
a truncated form to stop as soon as those required files are created:

```bash
# in workflows/rilseq-2016 or workflows/rilseq-2020:

conda activate ../../env-rilseq # activate the environment, if not already done

# optionally only run until trimmed fastq files:
snakemake --cores $CORES -s Snakefile-rnaseq --until cutadapt

# or run entire workflow to get all QC
snakemake --cores=$CORES -s Snakefile-rnaseq

```

The cutadapt-trimmed FASTQ files can then be found in
`workflows/$RILSEQ/data/rnaseq_samples/$SAMPLE/`.

## Run RIL-seq workflows

Run the `rilseq` workflow for each of `rilseq-2016` and `rilseq-2020`. This
uses the cutadapt-trimmed FASTQ files from the above workflow. It needs
a different environment activated (`env-rilseq`).

```bash
# in workflows/rilseq-2016 or workflows/rilseq-2020:

conda activate ../../env-rilseq
snakemake --cores $CORES -s Snakefile-rilseq
```

The main outputs of this workflow are the `all_fragments` and
`significant_fragments` files, stored in
`data/rnaseq_rilseq/{sample}/{sample}_all_fragments.txt` and
`data/rnaseq_rilseq/{sample}/{sample}_significant_fragments.txt`. These are
TSVs of single and chimeric reads, the later interpreted as interactions.
They need to be further processed to map the fragments to features. This
workflow also generates the Circos plots of interactions, stored in
`data/rnaseq_rilseq/{sample}/{sample}_circos_plot.svg`.


## Render RMarkdown files

Plots and tables summarizing the mapped fragments results are generated through
downstream RMarkdown files rendered into HTML. Computation time takes several
hours per condition. To support parallel execution, samples are processed in
subgroups and have been split into multiple subdirectories differing only by
the set of samples processed.

Each RMarkdown file in the `workflows/rilseq-*/downstream-*` directories should
be run in the `env-r` environment. Each directory has a `summary.Rmd` to be run
first, followed by `cds-shrink.Rmd`. Here, using
`workflows/rilseq-2016/downstream-log-noCL` as an example:

```bash
# from workflows/rilseq-2016/downstream-log-noCL:
conda activate ../../../env-r

Rscript -e "rmarkdown::render('summary.Rmd')"
Rscript -e "rmarkdown::render('cds-shrink.Rmd')"
```

The output will be an HTML file named after the respective RMarkdown file
(here, `workflows/rilseq-2016/downstream-log-noCL/summary.html` and
`workflows/rilseq-2016/downstream-log-noCL/cds-shrink.html`). These will also
create Excel files with the summarized RIL-seq results.
