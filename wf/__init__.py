from typing import List, Union

from latch import workflow
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


@workflow(METATAXANN_DOCS)
def metataxann(
    read1: LatchFile,
    read2: LatchFile,
    kaiju_ref_db: LatchFile,
    kaiju_ref_nodes: LatchFile,
    kaiju_ref_names: LatchFile,
    sample_name: str = "metataxann_sample",
    taxon_rank: str = "species",
    min_count: str = "2",
    k_min: str = "21",
    k_max: str = "141",
    k_step: str = "12",
    min_contig_len: str = "200",
    prodigal_output_format: str = "gbk",
) -> List[Union[LatchFile, LatchDir]]:
    """Metagenome assembly, binning and annotation

    MetaTaxAnn
    ----------

    # MetAssembly

    MetAssembly is a workflow for assembly of metagenomics data.
    It provides as end results both the assembled contigs as well as
    evaluation reports of said assembly.

    MetAssembly is a workflow composed of:
    - [MEGAHIT](https://github.com/voutcn/megahit) for assembly of input reads
    - [Quast](https://github.com/ablab/quast), specifically MetaQuast, for assembly evaluation.

    ---

    # Taxonomic classification with Kaiju

    Kaiju performs taxonomic classification of
    whole-genome sequencing metagenomics reads.
    Reads are assigned to taxa by using a reference database
    of protein sequences.
    Read more about it [here](https://github.com/bioinformatics-centre/kaiju)

    Kaiju paper: https://doi.org/10.1038/ncomms11257

    ---

    # Prodigal

    Prodigal is a protein-coding gene predictor for prokaryotic genomes.
    Read more about it [here](https://github.com/hyattpd/Prodigal).

    ---

    ### References

    Li, D., Luo, R., Liu, C.M., Leung, C.M., Ting, H.F., Sadakane, K., Yamashita, H. and Lam, T.W., 2016. MEGAHIT v1.0: A Fast and Scalable Metagenome Assembler driven by Advanced Methodologies and Community Practices. Methods.

    Alla Mikheenko, Vladislav Saveliev, Alexey Gurevich,
    MetaQUAST: evaluation of metagenome assemblies,
    Bioinformatics (2016) 32 (7): 1088-1090. doi: 10.1093/bioinformatics/btv697

    """
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

    annotation = prodigal(
        assembly_dir=assembly_dir,
        sample_name=sample_name,
        output_format=prodigal_output_format,
    )
    return [
        kaiju2table_out,
        krona_plot,
        metassembly_results,
        binning_results,
        annotation,
    ]
