from typing import List, Tuple, Union

from latch import workflow
from latch.resources.launch_plan import LaunchPlan
from latch.types import LatchDir, LatchFile

from .docs import METATAXANN_DOCS
from .kaiju import (
    kaiju2krona_task,
    kaiju2table_task,
    plot_krona_task,
    taxonomy_classification_task,
)
from .metassembly import megahit, metabat2, metaquast
from .prodigal import prodigal
from .types import ProdigalOutput, TaxonRank


@workflow
def tax_classification(
    read1: LatchFile,
    read2: LatchFile,
    kaiju_ref_db: LatchFile,
    kaiju_ref_nodes: LatchFile,
    kaiju_ref_names: LatchFile,
    sample_name: str,
    taxon_rank: TaxonRank,
) -> Tuple[LatchFile, LatchFile]:
    kaiju_out = taxonomy_classification_task(
        read1=read1,
        read2=read2,
        kaiju_ref_db=kaiju_ref_db,
        kaiju_ref_nodes=kaiju_ref_nodes,
        sample=sample_name,
    )
    kaiju2table_out = kaiju2table_task(
        kaiju_out=kaiju_out,
        sample=sample_name,
        kaiju_ref_nodes=kaiju_ref_nodes,
        kaiju_ref_names=kaiju_ref_names,
        taxon=taxon_rank,
    )
    kaiju2krona_out = kaiju2krona_task(
        kaiju_out=kaiju_out,
        sample=sample_name,
        kaiju_ref_nodes=kaiju_ref_nodes,
        kaiju_ref_names=kaiju_ref_names,
    )
    krona_plot = plot_krona_task(krona_txt=kaiju2krona_out, sample=sample_name)

    return kaiju2table_out, krona_plot


@workflow
def metassembly(
    read1: LatchFile,
    read2: LatchFile,
    sample_name: str,
    min_count: str,
    k_min: str,
    k_max: str,
    k_step: str,
    min_contig_len: str,
    prodigal_output_format: ProdigalOutput,
) -> Tuple[LatchDir, LatchDir, LatchDir]:

    assembly_dir = megahit(
        read_1=read1,
        read_2=read2,
        sample_name=sample_name,
        min_count=min_count,
        k_min=k_min,
        k_max=k_max,
        k_step=k_step,
        min_contig_len=min_contig_len,
    )
    metassembly_results = metaquast(assembly_dir=assembly_dir, sample_name=sample_name)
    binning_results = metabat2(assembly_dir=assembly_dir, sample_name=sample_name)
    ann = prodigal(
        assembly_dir=assembly_dir,
        sample_name=sample_name,
        output_format=prodigal_output_format,
    )

    return metassembly_results, binning_results, ann


@workflow(METATAXANN_DOCS)
def metataxann(
    read1: LatchFile,
    read2: LatchFile,
    kaiju_ref_db: LatchFile,
    kaiju_ref_nodes: LatchFile,
    kaiju_ref_names: LatchFile,
    sample_name: str = "metataxann_sample",
    taxon_rank: TaxonRank = TaxonRank.species,
    min_count: str = "2",
    k_min: str = "21",
    k_max: str = "141",
    k_step: str = "12",
    min_contig_len: str = "200",
    prodigal_output_format: ProdigalOutput = ProdigalOutput.gbk,
) -> List[Union[LatchFile, LatchDir]]:
    """Metagenomic assembly, binning and annotation

    MetaTaxAnn
    ----------

    MetaTaxAnn is a workflow for taxonomic classification, assembly, binning
    and annotation of long-read metagenomics datasets. It's composed of:

    - MEGAHIT for assembly
    - MetaQuast for assembly evaluation
    - MetaBAT2 for binning

    - Kaiju for taxonomic classification

    - Prodigal for protein-coding gene prediction on assembly data.

    In the next sections, you can read brief descriptions of
    all subworkflows contained within MetaTaxAnn.

    ---

    ## MetAssembly

    MetAssembly is a workflow for assembly of metagenomics data.
    It provides as end results both the assembled contigs, binned contigs
    evaluation reports of said assembly.

    MetAssembly is a workflow composed of:
    - [MEGAHIT](https://github.com/voutcn/megahit) for assembly of input reads
    - [Quast](https://github.com/ablab/quast), specifically MetaQuast, for assembly evaluation.
    - [MetaBAT2](https://bitbucket.org/berkeleylab/metabat/src/master/) for binning of assemblies

    ---

    ## Taxonomic classification with Kaiju

    Kaiju performs taxonomic classification of
    whole-genome sequencing metagenomics reads.
    Reads are assigned to taxa by using a reference database
    of protein sequences.
    Read more about it [here](https://github.com/bioinformatics-centre/kaiju)

    ---

    ## Prodigal

    Prodigal is a protein-coding gene predictor for prokaryotic genomes.
    Read more about it [here](https://github.com/hyattpd/Prodigal).

    ---

    ### References

    Li, D., Luo, R., Liu, C.M., Leung, C.M., Ting, H.F., Sadakane, K., Yamashita, H. and Lam, T.W., 2016. MEGAHIT v1.0: A Fast and Scalable Metagenome Assembler driven by Advanced Methodologies and Community Practices. Methods.

    Alla Mikheenko, Vladislav Saveliev, Alexey Gurevich,
    MetaQUAST: evaluation of metagenome assemblies,
    Bioinformatics (2016) 32 (7): 1088-1090. doi: 10.1093/bioinformatics/btv697

    Kang DD, Li F, Kirton E, Thomas A, Egan R, An H, Wang Z. 2019. MetaBAT 2: an
    adaptive binning algorithm for robust and efficient genome reconstruction
    from metagenome assemblies. PeerJ 7:e7359 https://doi.org/10.7717/peerj.7359

    Menzel, P., Ng, K. & Krogh, A. Fast and sensitive taxonomic classification for
    metagenomics with Kaiju. Nat Commun 7, 11257 (2016).
    https://doi.org/10.1038/ncomms11257

    Hyatt, D., Chen, GL., LoCascio, P.F. et al. Prodigal: prokaryotic gene recognition
    and translation initiation site identification.
    BMC Bioinformatics 11, 119 (2010). https://doi.org/10.1186/1471-2105-11-119
    """
    kaiju2table_out, krona_plot = tax_classification(
        read1=read1,
        read2=read2,
        kaiju_ref_db=kaiju_ref_db,
        kaiju_ref_nodes=kaiju_ref_nodes,
        kaiju_ref_names=kaiju_ref_names,
        sample_name=sample_name,
        taxon_rank=taxon_rank,
    )

    metassembly_results, binning_results, ann = metassembly(
        read1=read1,
        read2=read2,
        sample_name=sample_name,
        min_count=min_count,
        k_min=k_min,
        k_max=k_max,
        k_step=k_step,
        min_contig_len=min_contig_len,
        prodigal_output_format=prodigal_output_format,
    )

    return [kaiju2table_out, krona_plot, metassembly_results, binning_results, ann]


LaunchPlan(
    metataxann,  # workflow name
    "Example Metagenome (Crohn's disease gut microbiome)",  # name of test data
    {
        "read1": LatchFile("s3://latch-public/test-data/4318/SRR579292_1.fastq"),
        "read2": LatchFile("s3://latch-public/test-data/4318/SRR579292_2.fastq"),
        "kaiju_ref_db": LatchFile(
            "s3://latch-public/test-data/4318/kaiju_db_plasmids.fmi"
        ),
        "kaiju_ref_nodes": LatchFile("s3://latch-public/test-data/4318/nodes.dmp"),
        "kaiju_ref_names": LatchFile("s3://latch-public/test-data/4318/names.dmp"),
        "sample_name": "SRR579292",
        "taxon_rank": TaxonRank.species,
        "min_count": "2",
        "k_min": "21",
        "k_max": "141",
        "k_step": "12",
        "min_contig_len": "200",
        "prodigal_output_format": ProdigalOutput.gff,
    },
)
