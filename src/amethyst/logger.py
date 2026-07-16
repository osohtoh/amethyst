import logging

logger = logging.getLogger("amethyst")
logger.setLevel(logging.DEBUG)

if not logger.handlers:
    _formatter = logging.Formatter(
        "%(asctime)s | %(levelname)-8s | %(name)s:%(funcName)s:%(lineno)d - %(message)s"
    )
    _file_handler = logging.FileHandler("amethyst.log")
    _file_handler.setFormatter(_formatter)
    logger.addHandler(_file_handler)
