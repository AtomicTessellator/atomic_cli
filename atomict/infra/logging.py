from logging.config import dictConfig
import google.cloud.logging

def config_loggers(prefix: str = '', *args, **kwargs):
    client = google.cloud.logging.Client()
    logging_config = dict(
        version=1,
        formatters={
            "verbose": {"format": f"%(levelname)s %(asctime)s {prefix} %(message)s"}
        },
        handlers={
            "cloud_logging_handler": {
                "class": "google.cloud.logging.handlers.CloudLoggingHandler",
                "client": client,
            },
            "structured_log_handler": {
                "class": "google.cloud.logging.handlers.StructuredLogHandler"
            },
            "console": {
                "class": "logging.StreamHandler", "formatter": "verbose"
            }
        },
        root={"handlers": ["console"], "level": "INFO"},
        loggers={
            "cloud_logger": {"handlers": ["cloud_logging_handler"], "level": "INFO"},
            "structured_logger": {
                "handlers": ["structured_log_handler"],
                "level": "INFO",
            },
            "atomict.api": {
                "handlers": ["console"],
                "level": "INFO",
                "propagate": False
            }
        }
    )
    dictConfig(logging_config)
