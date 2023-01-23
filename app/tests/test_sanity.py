import os
from logging import getLogger

# here we need to test both command line tools, with particular concern for dependency resolution

logger = getLogger(__name__)

def test_blast_defaults(capsys):
    exit_status = os.system('python -m src.BLAST_probes')
    captured = capsys.readouterr()
    print(captured.out)
    print(captured.err)
    logger.error(captured.err)
    assert exit_status == 0

def test_kmers_defaults(capsys):
    exit_status = os.system('python -m src.target_kmers -max_Tm 60')
    captured = capsys.readouterr()
    print(captured.out)
    print(captured.err)
    logger.error(captured.err)
    assert exit_status == 0