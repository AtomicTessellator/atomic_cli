from logging.config import dictConfig


def config_loggers(*args, **kwargs):
    logging_config = dict(
        version=1,
        formatters={"verbose": {"format": "%(levelname)s %(asctime)s %(message)s"}},
        handlers={
            "console": {"class": "logging.StreamHandler", "formatter": "verbose"}
        },
        root={"handlers": ["console"], "level": "INFO"},
    )

    dictConfig(logging_config)
