import os
from logging import getLogger
from pathlib import Path

# for now testing the command line tools only, not individual functions

logger = getLogger(__name__)


def test_blast_defaults(capsys, tmp_path):
    input_filepath = Path(__file__).parent.parent / "static" / "non_overlap_probes.csv"
    assert len(list(tmp_path.glob("*.csv"))) == 0
    exit_status = os.system(
        f"python -m angler.BLAST_probes -input_filepath {input_filepath} -output_dir {tmp_path}"
    )
    captured = capsys.readouterr()
    print(captured.out)
    print(captured.err)
    logger.error(captured.err)
    assert exit_status == 0
    assert len(list(tmp_path.glob("*.csv"))) == 1


def test_kmers_defaults(capsys, tmp_path):
    input_dir = Path(__file__).parent.parent / "static" / "mrna_fasta"
    assert len(list(tmp_path.glob("*.csv"))) == 0
    exit_status = os.system(
        f"python -m angler.target_kmers -max_Tm 60 -input_dir {input_dir} -output_dir {tmp_path}"
    )
    captured = capsys.readouterr()
    print(captured.out)
    print(captured.err)
    logger.error(captured.err)
    assert exit_status == 0
    assert len(list(tmp_path.glob("*.csv"))) == 1
