from logging.config import dictConfig


def config_loggers(prefix: str = '', *args, **kwargs):
    logging_config = dict(
        version=1,
        formatters={
            "verbose": {"format": f"%(levelname)s %(asctime)s {prefix} %(message)s"}
        },
        handlers={
            "console": {"class": "logging.StreamHandler", "formatter": "verbose"}
        },
        root={"handlers": ["console"], "level": "INFO"},
    )
    dictConfig(logging_config)
