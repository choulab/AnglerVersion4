import argparse
from datetime import datetime
import json
from logging import getLogger
from pathlib import Path
from tempfile import mkstemp

from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd

from angler.util import configure_logger, validate_directory

logger = getLogger(__name__)


def transcriptome_blast(query_fasta, organism="House mouse"):
    """
    BLAST a multi fasta of probe candidates against a user-selected transcriptome, returns \
       a single JSON contaning the BLAST results 
    """

    _, save_path = mkstemp(".json", text=True)
    logger.debug(f"Saving multi fasta blast results to to {save_path}")

    cline = NcbiblastnCommandline(
        query=query_fasta,
        task="blastn-short",
        db="refseq_rna",
        evalue=1000,
        gapopen=5,
        gapextend=2,
        out=save_path,
        outfmt=15,
        entrez_query='"{} [ORGN]"'.format(organism),
        remote=True,
        reward=1,
        penalty=-3,
        word_size=7,
        strand="minus",
    )

    cline()
    return save_path


def parse_probes(fpath):
    """
    A function to assemble probe candidates into a BLASTable format
    """

    db = pd.read_csv(fpath)
    records = []
    for kmer in db.to_dict(orient="records"):
        record_1 = SeqRecord(seq=Seq(kmer["HCR3-P1"]), id=kmer["name"] + "-P1")
        record_2 = SeqRecord(seq=Seq(kmer["HCR3-P2"]), id=kmer["name"] + "-P2")
        records.append(record_1)
        records.append(record_2)

    _, save_path = mkstemp(".fasta")
    logger.debug(f"Saving probe candidates to {save_path}")

    SeqIO.write(records, save_path, "fasta")

    return save_path


def parse_results(fpath):
    output_1 = []
    output_2 = []

    with open(fpath) as f_in:
        results = json.load(f_in)

        for probe in range(len(results["BlastOutput2"])):
            kmer = results["BlastOutput2"][probe]["report"]["results"]["search"][
                "query_title"
            ]
            hits = results["BlastOutput2"][probe]["report"]["results"]["search"]["hits"]
            align_lengths = [
                {hit["description"][0]["title"]: hit["hsps"][0]["align_len"]}
                for hit in hits
                if hit["hsps"][0]["qseq"] == hit["hsps"][0]["hseq"]
            ]

            if len(align_lengths) == 0:
                align_lengths = ["NA"]

            if "-P1" in kmer:
                output_1.append(align_lengths[:20])
            else:
                output_2.append(align_lengths[:20])

    return [output_1, output_2]


def run(organism: str, input_filepath: str, output_dir: str):
    """
    Get the BLAST probes

        Parameters
            organism (str) : the organism to analyze
            input_filepath : the path to the file to analyze
            output_dir (str): the absolute path to the directory where the results should be written

        Returns:
            (str) the path to the results file
    """

    logger.debug(f"running BLAST probes on files located at {input_filepath}")


    output_dir = validate_directory(output_dir)

    if not Path(input_filepath).exists:
        raise FileNotFoundError(f"{input_filepath} does not exist!")

    probe_candidate_path = parse_probes(input_filepath)
    blast_results_path = transcriptome_blast(probe_candidate_path, organism=organism)

    probes_db = pd.read_csv(input_filepath)
    probes_db["BLAST_1"] = parse_results(blast_results_path)[0]
    probes_db["BLAST_2"] = parse_results(blast_results_path)[1]

    blast_probes_save_path = (
        output_dir
        / f"blast-probe-results-{datetime.now().strftime('%Y-%m-%d-%H%M%S')}.csv"
    )

    probes_db.to_csv(blast_probes_save_path)

    logger.debug(f"Saving BLAST results to {blast_probes_save_path}")

    return blast_probes_save_path


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="AnglerFISH, BLAST Script. Version 1.0."
    )
    parser.add_argument(
        "-organism",
        type=str,
        default="Mus musculus",
        required=False,
        help="Enter the name of the organism whose transcriptome you wish to BLAST probes against. Ex. Mus musculus.",
    )

    parser.add_argument(
        "-log_path",
        type=str,
        default=Path(__file__).parent.parent.absolute() / "log",
        required=False,
        help="The directory where the log file should be placed.",
    )

    parser.add_argument(
        "-debug",
        type=bool,
        default=False,
        required=False,
        help="Whether to enable verbose logging.",
    )

    parser.add_argument(
        "-output_dir",
        type=str,
        required=True,
        help="The absolute path to the directory where the results should be written.",
    )

    parser.add_argument(
        "-input_filepath",
        type=str,
        required=True,
        help="The absolute path to the probe csv.",
    )

    args = parser.parse_args()

    configure_logger(args.log_path, args.debug)

    run(args.organism, args.input_filepath, args.output_dir)
