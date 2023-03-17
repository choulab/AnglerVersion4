import argparse
from datetime import datetime
from logging import getLogger, DEBUG
import os
from pathlib import Path
import random
from tempfile import mkstemp

from Bio import AlignIO, SeqIO
from Bio.Align.Applications import MuscleCommandline
from Bio.SeqUtils import GC, MeltingTemp as mt
import Levenshtein as lev
from nupack import *
import pandas as pd
from tqdm import tqdm

from angler.util import configure_logger, validate_directory

logger = getLogger(__name__)


def flatten(t):
    """
    Flattens a list of lists containing tuples of nucleotide letters
    """

    return [item for sublist in t for item in sublist]


def probeTm(seq, fmd=30):
    """
    Calculates the Tm of a given sequence.
    """

    tmval = float(("%0.2f" % mt.Tm_NN(seq, Na=975, Tris=75)))
    fcorrected = "%0.2f" % mt.chem_correction(tmval, fmd=fmd)
    return fcorrected


def calc_energy(seq, T=37, fmd=30, Na=0.975, Mg=0):
    """
    Calculates ∆G of hybridization of given sequence to its perfect complement using NUPACK
    params:
        seq (str) : the sequence to be analyzed
    returns:
        seq_energy (float) : ∆G of hybridization of sequence in kcal/mol
    """

    T = T + (0.62 * fmd)
    len_seq = len(seq)
    seq_struct = len_seq * "(" + "+" + len_seq * ")"
    model1 = Model(
        material="dna", ensemble="stacking", celsius=T, sodium=Na, magnesium=Mg
    )
    seq_energy = structure_energy(
        strands=[seq, revcom(seq)], structure=seq_struct, model=model1
    )
    return seq_energy


def revcom(seq):
    """
    Function to calculate reverse complement of a sequence
    """

    alt_map = {"ins": "0"}
    complement = {"A": "T", "C": "G", "G": "C", "T": "A"}

    for k, v in alt_map.items():
        seq = seq.replace(k, v)
    bases = list(seq)
    bases = reversed([complement.get(base, base) for base in bases])
    bases = "".join(bases)
    for k, v in alt_map.items():
        bases = bases.replace(v, k)
    return bases


def window(seq, name, n=3, temp=37, form=30):
    """
    Identifies k-mers of a given length. Returns a list of dictionaries annotated with k-mer and transcript coordinates.
    params:
        n : k-mer length
    """

    print("Calculating {}-mers for {}".format(n, name))
    output = []
    for count in tqdm(range(n, len(seq) + 1)):
        kmer = seq[count - n : count].upper()
        revcom_kmer = revcom(kmer)
        output.append(
            {
                "name": name,
                "start": count - n,
                "stop": count,
                "probe": revcom_kmer,
                "reference": kmer,
                "dG": calc_energy(revcom_kmer, T=temp, fmd=form),
                "Tm": probeTm(revcom_kmer, fmd=form),
                "length": len(revcom_kmer),
                "GC": GC(revcom_kmer),
                "GC_HCR3-P1": GC(revcom_kmer[:25]),
                "GC_HCR3-P2": GC(revcom_kmer[27:]),
                "length_HCR3-P1": len(revcom_kmer[:25]),
                "length_HCR3-P2": len(revcom_kmer[27:]),
                "Tm_HCR3-P1": probeTm(revcom_kmer[:25]),
                "Tm_HCR3-P2": probeTm(revcom_kmer[27:]),
            }
        )
    return output


def create_kmer(file, min, max, form=30, temp=37):
    """
    Identifies all k-mers across a sliding window given a FASTA file for the mRNA target of interest
    """

    fasta_sequences = SeqIO.parse(open(file), "fasta")
    kmer_db = []
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        counter = min
        while counter <= max:
            kmer_db.append(window(sequence, name, n=counter, temp=temp, form=form))
            counter += 1
    return flatten(kmer_db)


def concat_fasta(dir):
    """
    Concatenates each of the mRNA FASTA files into a single FASTA file prior to multi sequence alignment
    """
    _, save_path = mkstemp(".fasta")
    logger.debug(f"Saving concatenated fasta to {save_path}")
    records = []
    for file in os.listdir(dir):
        for rec in SeqIO.parse(open(dir / file), "fasta"):
            records.append(rec)
    SeqIO.write(records, save_path, "fasta")
    return save_path


def multi_alignment(fpath):
    """
    Utilizes the MUSCLE alignment tool to generate a multi sequence alignment of each mRNA sequence
    """

    _, muscle_log_path = mkstemp()
    _, save_path = mkstemp("fasta")

    logger.debug("Saving muscle logs to {muscle_log_path}")
    logger.debug("Saving multi alignment fast file to {save_path}")

    muscle_cline = MuscleCommandline(
        input=fpath,
        out=save_path,
        diags=True,
        maxiters=3,
        verbose=True,
        log=muscle_log_path,
    )
    muscle_cline()
    return save_path


def load_alignment(fpath):
    """
    Loads the alignment file to a biopython Aligh.IO object
    """

    alignment = AlignIO.read(fpath, "fasta")
    print("Alignment length %i" % alignment.get_alignment_length())
    return alignment


def lev_distance(str1, str2):
    """calculates the levenshtein distance between two strings"""

    return lev.distance(str1.lower(), str2.lower())


def partial_fuzz(str1, str2):
    """calculates the maximum levenshtein distance between a short and long string"""

    if len(str1) > len(str2):
        seq = str1
        kmer = str2
    else:
        seq = str2
        kmer = str1

    distances = []

    for count in range(len(kmer), len(seq) + 1):
        ref = seq[count - len(kmer) : count]
        distances.append(lev_distance(ref, kmer))

    return min(distances)


def check_isoforms(kmer_dict, multi_fasta_fpath):
    """
    Returns information about which isoforms each k-mer is specific to
    """

    print("Calculating Probe Alignments. This may take several minutes.")

    isoform_kmers = []
    for kmer in tqdm(kmer_dict):
        kmer["isoforms_P1"] = {}
        kmer["isoforms_P2"] = {}
        for fasta in SeqIO.parse(open(multi_fasta_fpath), "fasta"):
            alignment_value_1 = partial_fuzz(kmer["reference"][:25], str(fasta.seq))
            alignment_value_2 = partial_fuzz(kmer["reference"][27:], str(fasta.seq))

            kmer["isoforms_P1"][fasta.id] = alignment_value_1
            kmer["isoforms_P2"][fasta.id] = alignment_value_2

        isoform_kmers.append(kmer)
    return isoform_kmers


def SNVscore(inFile):
    """
    Takes a multi-sequence alignment file and determines SNV score for each nucleotide position
    """

    seqs = [str(record.seq) for record in SeqIO.parse(inFile, "fasta")]
    seq_len, num_isoforms = len(seqs[0]), len(seqs)
    print("{} Isoforms".format(num_isoforms))
    snv_score = (
        []
    )  # snv scores for each nucleotide position in the multi-sequence alignment file

    # Crawl through each nt position in given multi-sequence alignment file
    for nt in range(0, seq_len):
        x = []

        # list of nt sequences for each isoform at specified nt position
        for isoform in seqs:
            x.append(isoform[nt])

        # scores variation based on frequency of bases
        score = (
            1
            / (
                max(
                    x.count("A"), x.count("T"), x.count("C"), x.count("G"), x.count("-")
                )
                / num_isoforms
            )
            - 1
        )
        snv_score.append(score)

    return snv_score


def SNV_sum(pos, probe_len, snv_score):
    """
    Calculates the sum of the SNV scores for a given probe along a multi-sequence alignemnt
    """

    score = round(sum(snv_score[pos : pos + probe_len]) / probe_len, 2)
    return score


def map_kmer(kmer, reference):
    """
    Locates a k-mer along an aligned reference sequence inclusive of gaps and returns the
    start coordinate in the multi-sequence alignment file
    """

    reference = str(reference)

    for nt in range(0, len(reference)):
        remaining = reference[nt:]
        remaining_cleaned = reference[nt:].replace("-", "")

        if remaining_cleaned[: len(kmer)] == kmer and remaining[0] != "-":
            return nt


def filter_db(
    db,
    prohibited_seqs=["AAAA", "TTTT", "CCCC", "GGGG"],
    min_Tm=40,
    max_Tm=60,
    min_GC=35,
    max_GC=55,
    min_dG=-100,
    max_dG=100,
    lev_cutoff=3,
):
    """
    Eliminates k-mers from JSON database according to user-defined criterion.
    """

    filtered_db = []
    for kmer in tqdm(db):
        res = any(ele in kmer["reference"] for ele in prohibited_seqs)
        if res == False:
            if min_Tm <= float(kmer["Tm_HCR3-P1"]) <= max_Tm:
                if min_Tm <= float(kmer["Tm_HCR3-P2"]) <= max_Tm:
                    if min_GC <= float(kmer["GC_HCR3-P1"]) <= max_GC:
                        if min_GC <= float(kmer["GC_HCR3-P2"]) <= max_GC:
                            if min_dG <= float(kmer["dG"]) <= max_dG:
                                if kmer["MSA_position"] != None:
                                    if (
                                        sorted(list(kmer["isoforms_P1"].values()))[1]
                                        >= lev_cutoff
                                    ):
                                        if (
                                            sorted(list(kmer["isoforms_P2"].values()))[
                                                1
                                            ]
                                            >= lev_cutoff
                                        ):
                                            filtered_db.append(kmer)
    return filtered_db


def sort_isoform(db):
    """
    Given a flat database of k-mers, this function idexes them by parent isoform
    """

    isoforms = {}
    for kmer in db:
        if kmer["name"] not in isoforms.keys():
            isoforms[kmer["name"]] = [kmer]
        else:
            isoforms[kmer["name"]].append(kmer)
    for isoform in isoforms.keys():
        for idx, kmer in enumerate(isoforms[isoform]):
            kmer["name"] = "{}-{}".format(kmer["name"], idx)
    return isoforms


def remove_overlaps(kmer_list, min_spacing=3):
    """
    Compare overlapping target regions and pick the one with better SNV score
    """

    seed = random.randrange(0, len(kmer_list))
    champions = []
    champions.append(kmer_list[seed])

    for contender in kmer_list:
        # check whether contender overlaps with existing probes
        l = len(contender["reference"]) + min_spacing
        d = [abs(probe["start"] - contender["start"]) < l for probe in champions]
        if True in d:
            idxs_ = [i for i, x in enumerate(d) if x]
            if len(idxs_) == 1:  # if this does not overlap multiple existing probes
                # compare who is better
                idx_ = d.index(True)
                if contender["HCR3_SNV"] > champions[idx_]["HCR3_SNV"]:
                    champions[idx_] = contender
        else:
            champions.append(contender)
    return champions


def append_HCR(db, initiator):
    if initiator == "B2":
        db["HCR3-P1-B2"] = db["HCR3-P1"].astype(str) + "AAATCATCCAGTAAACCGCC"
        db["HCR3-P2-B2"] = "CCTCGTAAATCCTCATCAAA" + db["HCR3-P2"].astype(str)
    elif initiator == "B3":
        db["HCR3-P1-B3"] = db["HCR3-P2"].astype(str) + "TTCCACTCAACTTTAACCCG"
        db["HCR3-P2-B3"] = "GTCCCTGCCTCTATATCTTT" + db["HCR3-P2"].astype(str)
    else:
        raise ValueError("Invalid HCR Initiator")
    return db


def run(
    input_dir: str,
    output_dir: str,
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
    version=1.00,
):
    """
    Run the target_kmers main function

        Parameters
            input_dir (str) : the directory of the file to be analyzed
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
            version (float) : Program version.

        Returns:
            (str) : Path to the results CSV file containing probes.
    """

    logger.debug("AnglerFISH v.{}".format(version))

    for arg, value in sorted(vars(args).items()):
        logger.info("{} = {}".format(arg, value))

    input_dir = validate_directory(input_dir)
    output_dir = validate_directory(output_dir)

    mrna_fnames = os.listdir(input_dir)
    # create a database of probe candidates
    concat_fasta_path = concat_fasta(input_dir)
    multi_align_path = multi_alignment(concat_fasta_path)
    aligned_sequences = load_alignment(multi_align_path)
    isoforms_dict = {rec.id: rec.seq for rec in aligned_sequences}
    conservation_score = SNVscore(multi_align_path)

    # calculate relevant features for each probe candidate
    output_db = []

    for target in tqdm(mrna_fnames, desc="Overall Progress"):
        file = input_dir / target
        output = create_kmer(
            file, min_length, max_length, form=formamide, temp=hyb_temp
        )
        isoform_output = check_isoforms(output, concat_fasta_path)

        for kmer in tqdm(
            isoform_output, desc="Mapping probes to MSA. This may take several minutes."
        ):
            kmer["MSA_position"] = map_kmer(
                kmer["reference"], isoforms_dict[kmer["name"]]
            )

            if kmer["MSA_position"] != None:
                kmer["SNV_score"] = SNV_sum(
                    kmer["MSA_position"], kmer["length"], conservation_score
                )
                HCR3_probe_len = int(len(kmer["probe"]) / 2) - 1
                kmer["HCR3-P1"] = kmer["probe"][:HCR3_probe_len]
                kmer["HCR3-P2"] = kmer["probe"][HCR3_probe_len + 2 :]
                P1_jun_pos = (
                    kmer["MSA_position"] + HCR3_probe_len - critical_junction_length
                )
                P2_jun_pos = kmer["MSA_position"] + HCR3_probe_len + 2
                kmer["HCR3-P1_junction_SNV"] = SNV_sum(
                    P1_jun_pos, critical_junction_length, conservation_score
                )
                kmer["HCR3-P2-_junction_SNV"] = SNV_sum(
                    P2_jun_pos, critical_junction_length, conservation_score
                )
                kmer["HCR3_SNV"] = (
                    kmer["HCR3-P1_junction_SNV"] + kmer["HCR3-P2-_junction_SNV"]
                )
                output_db.append(kmer)

            else:
                print("ERROR: Unable to map probe candidate to MSA")

    # remove probe candidates that do not compply with user-specified criteria
    output_db = filter_db(
        output_db,
        min_Tm=min_Tm,
        max_Tm=max_Tm,
        lev_cutoff=max_lev,
        min_GC=min_GC,
        max_GC=max_GC,
    )
    isoforms = sort_isoform(output_db)

    # Instantiate an isoform-indexed dictionary of probe candidates
    iso_flat_dict = {}

    for i in isoforms.keys():
        non_overlap = remove_overlaps(isoforms[i], min_spacing=spacing)
        iso_flat_dict[i] = sorted(non_overlap, key=lambda d: d["start"])

    output_list = []

    for k in iso_flat_dict:
        for j in iso_flat_dict[k]:
            output_list.append(j)

    db_out = pd.DataFrame(output_list)
    logger.debug(db_out)
    db_HCR = append_HCR(db_out, initiator="B2")
    output_filepath = (
        output_dir
        / f"target-kmers-results-{datetime.now().strftime('%Y-%m-%d-%H%M%S')}.csv"
    )
    db_HCR.to_csv(output_filepath, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="AnglerFISH. Version 1.0.")
    parser.add_argument(
        "-input_dir",
        type=str,
        required=True,
        help="The absolute path to the directory containing the file(s) to be analyzed.",
    )
    parser.add_argument(
        "-output_dir",
        type=str,
        required=True,
        help="The absolute path to the directory where the results should be written.",
    )
    parser.add_argument(
        "-min_GC",
        type=int,
        default=35,
        required=False,
        help="The minumum allowable GC content for a probe candidate. Default is 35 percent.",
    )
    parser.add_argument(
        "-max_GC",
        type=int,
        default=55,
        required=False,
        help="The maximum allowable GC content for a probe candidate. Default is 35 percent.",
    )
    parser.add_argument(
        "-S",
        type=int,
        default=3,
        required=False,
        help="The minumum spacing between probes in nucleotides. Default is 3 nucleotides.",
    )
    parser.add_argument(
        "-L",
        type=int,
        default=3,
        required=False,
        help="The allowable Levenshtein Distance between a probe and an isoform. Default is 3 nucleotides.",
    )
    parser.add_argument(
        "-F",
        type=int,
        default=30,
        required=False,
        help="Hybridization buffer formamide content. Default is 30 percent.",
    )
    parser.add_argument(
        "-T",
        type=int,
        default=37,
        required=False,
        help="Probe hybridization temperature in degrees Celsius. Default is 37 Celsius.",
    )
    parser.add_argument(
        "-J",
        type=int,
        default=3,
        required=False,
        help="The size of the optimized critical junction for HCR 3.0 Probes. Default is 3 nucleotides.",
    )
    parser.add_argument(
        "-min_Tm",
        type=int,
        default=40,
        required=False,
        help="The minumum allowable formamide-adjusted Tm for a probe candidate. Default is 40 Celsius.",
    )
    parser.add_argument(
        "-max_Tm",
        type=int,
        default=60,
        required=True,
        help="The maximum allowable formamide-adjusted Tm for a probe candidate. Default is 60 Celsius.",
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

    args = parser.parse_args()

    configure_logger(
        args.log_path,
        args.debug,
    )

    run(
        args.input_dir,
        args.output_dir,
        args.min_GC,
        args.max_GC,
        args.S,
        args.min_Tm,
        args.max_Tm,
        args.L,
        args.F,
        args.T,
        args.J,
        52,
        52,
        1.00,
    )
