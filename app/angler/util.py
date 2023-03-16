from logging import basicConfig, DEBUG, INFO
from os.path import isdir, isabs, join
from tempfile import mkdtemp


def configure_logger(path: str = None, verbose=False):
    """
    Set the path for the log files and the log level
    """
    if path and not isdir(path):
        raise FileNotFoundError(f"{path} either does not exist or is a directory!")
    elif not isabs(path):
        raise ValueError(f"Log path must be absolute!")
    elif path:
        log_path = join(path, "angler.log")
    else:
        log_path = join(mkdtemp(), "angler.log")

    if verbose:
        print(f"Logs will be written to {log_path}")
        level = DEBUG
    else:
        level = INFO

    basicConfig(
        filename=log_path,
        format="%(asctime)s %(message)s",
        encoding="utf-8",
        level=level,
    )

    pass
