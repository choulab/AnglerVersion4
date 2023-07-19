from logging import basicConfig, DEBUG, INFO
from logging.handlers import TimedRotatingFileHandler
from os.path import isdir, isabs, join
from pathlib import Path
from tempfile import mkdtemp


def validate_directory(dir_path: str):
    """
    Validate that provided path is absolute and is a directory

        Parameters
            dir_path (str) : the path to a directory

        Returns:
            (bool) True

        Raises:
            Path
    """
    dir_path = Path(dir_path)

    if not dir_path.is_absolute():
        raise ValueError(f"{dir_path} is not absolute!")

    if not dir_path.is_dir():
        raise ValueError(f"{dir_path} does not exist or is not a directory")

    return dir_path


def configure_logger(path: str = None, verbose=False):
    """
    Set the path for the log files and the log level

        Parameters
            path (str | None): the directory where the the logs should be stored
            verbose (bool): whether to enable verbose logs (sets logging level to DEBUG)
    """
    if path and not isdir(path):
        raise FileNotFoundError(f"{path} either does not exist or is not a directory!")
    elif not isabs(path):
        raise ValueError(f"Log path must be absolute!")
    elif path:
        log_path = join(path, "angler.log")
    else:
        log_path = join(mkdtemp(), "angler.log")

    print(f"Logs will be written to {log_path}")

    if verbose:
        level = DEBUG
        print("**VERBOSE LOGGING ENABLED**")
    else:
        level = INFO

    handler = TimedRotatingFileHandler(log_path, when="D", interval=1, backupCount=10)

    basicConfig(
        handlers=[handler],
        format="%(asctime)s %(message)s",
        encoding="utf-8",
        level=level,
    )
