import argparse
from pathlib import Path
from tempfile import mkdtemp
from logging import getLogger

from angler.BLAST_probes import run as run_blast
from angler.target_kmers import run as run_kmers
from angler.util import configure_logger

logger = getLogger(__name__)


def run_angler(
    input_dir: str,
    output_dir: str,
    organism="Mus musculus",
    min_GC=35,
    max_GC=55,
    spacing=3,
    min_Tm=40,
    max_Tm=60,
    max_lev=3,
    formamide=30,
    hyb_temp=37,
    critical_junction_length=3,
    min_length=52,
    max_length=52,
    debug=False,
    log_dir=Path(__file__).parent.parent.absolute() / "log",
):
    """
    Run the target_kmers and then blast_probes functions, passing the output of one to the other

        Parameters
            input_dir (str) : the directory of the file to be analyzed.
            output_dir (str) : the directory where the results will be stored.
            organism (str) : the organism to analyze.
            min (int) : The minumum allowable %GC content for a probe candidate.
            max_GC (int) : The maximum allowable %GC content for a probe candidate.
            spacing (int) : The minumum spacing between probes in nucleotides.
            min_Tm (int) : The minumum allowable formamide-adjusted Tm for a probe candidate.
            max_Tm (int) : The maximum allowable formamide-adjusted Tm for a probe candidate.
            max_lev (int) : The allowable Levenshtein Distance between a probe and an isoform.
            formamide (int) : Hybridization buffer formamide content.
            hyb_temp (float) : Probe hybridization temperature in degrees Celsius.
            critical_junction_length (float) : The size of the optimized critical junction for HCR 3.0 Probes. Default is 3 nucleotides.
            min_length (int)
            max_length (int)

        Returns: str
    """

    configure_logger(log_dir, debug)

    kmers_output_dir = mkdtemp()

    logger.debug(f"Saving kmers results to temp directory {kmers_output_dir}")

    probe_path = run_kmers(
        input_dir=input_dir,
        output_dir=kmers_output_dir,
        min_GC=min_GC,
        max_GC=max_GC,
        spacing=spacing,
        min_Tm=min_Tm,
        max_Tm=max_Tm,
        max_lev=max_lev,
        formamide=formamide,
        hyb_temp=hyb_temp,
        critical_junction_length=critical_junction_length,
        min_length=min_length,
        max_length=max_length,
    )

    logger.debug("kmers complete.")
    logger.debug(f"Running BLAST probes and saving to {output_dir}")
    print("Running BLAST probes. This may take several minutes....")

    try:
        result_path = run_blast(
            organism=organism, input_filepath=probe_path, output_dir=output_dir
        )
    except Exception as e:
        logger.error(e)
        print("BLAST probe creation failed. See logs for details.")
        return

    print(
        f"""Probe genreation complete! The result file is located at {result_path}.
        Note that if you are running this command via Docker, the actual parent directory will be the same one you passed as --ouput_dir."""
    )

    return result_path


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="AnglerFISH. Version 2.0.1")
    parser.add_argument(
        "--critical_junction_length",
        type=int,
        dest="critical_junction_length",
        default=3,
        required=False,
        help="The size of the optimized critical junction for HCR 3.0 Probes. Default is 3 nucleotides.",
    )
    parser.add_argument(
        "--debug",
        dest="debug",
        action="store_true",
        help="Enable verbose logging.",
    )

    parser.add_argument(
        "--formamide",
        dest="formamide",
        type=int,
        default=30,
        required=False,
        help="Hybridization buffer formamide content. Default is 30 percent.",
    )
    parser.add_argument(
        "--hyb_temp",
        dest="hyb_temp",
        type=int,
        default=37,
        required=False,
        help="Probe hybridization temperature in degrees Celsius. Default is 37 Celsius.",
    )
    parser.add_argument(
        "--input_dir",
        dest="input_dir",
        type=str,
        required=True,
        help="Absolute path to directory containing input files.",
    )
    parser.add_argument(
        "--log_dir",
        dest="log_dir",
        type=str,
        required=True,
        help="Absolute path to directory where logs should be saved.",
    )
    parser.add_argument(
        "--max_GC",
        dest="max_GC",
        type=int,
        default=55,
        required=False,
        help="The maximum allowable GC content for a probe candidate. Default is 35 percent.",
    )
    parser.add_argument(
        "--max_lev",
        dest="max_lev",
        type=int,
        default=3,
        required=False,
        help="The allowable Levenshtein Distance between a probe and an isoform. Default is 3 nucleotides.",
    )

    parser.add_argument(
        "--max_Tm",
        dest="max_Tm",
        type=int,
        default=60,
        required=False,
        help="The maximum allowable formamide-adjusted Tm for a probe candidate. Default is 60 Celsius.",
    )
    parser.add_argument(
        "--min_GC",
        dest="min_GC",
        type=int,
        default=35,
        required=False,
        help="The minumum allowable GC content for a probe candidate. Default is 35 percent.",
    )
    parser.add_argument(
        "--min_Tm",
        dest="min_Tm",
        type=int,
        default=40,
        required=False,
        help="The minumum allowable formamide-adjusted Tm for a probe candidate. Default is 40 Celsius.",
    )
    parser.add_argument(
        "--organism",
        dest="organism",
        type=int,
        default=3,
        required=False,
        help="The organism to analyze. Default is Mus musculus.",
    )
    parser.add_argument(
        "--output_dir",
        dest="output_dir",
        type=str,
        required=True,
        help="Absolute path to directory where output files should be saved.",
    )
    parser.add_argument(
        "--spacing",
        dest="spacing",
        type=int,
        default=3,
        required=False,
        help="The minumum spacing between probes in nucleotides. Default is 3 nucleotides.",
    )

    args = parser.parse_args()

    run_angler(**vars(args))
