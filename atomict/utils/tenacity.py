
import typing

from tenacity import _utils

if typing.TYPE_CHECKING:
    import logging

    from tenacity import RetryCallState


def before_log(
    logger: "logging.Logger", log_level: int
) -> typing.Callable[["RetryCallState"], None]:
    """
    Before call strategy that logs to some logger the attempt.

    This one is a copy of the tenacity.before_log function, but it doesnt log the first
    attempt because that generates a lot of noise.
    """

    def log_it(retry_state: "RetryCallState") -> None:
        if retry_state.fn is None:
            # NOTE(sileht): can't really happen, but we must please mypy
            fn_name = "<unknown>"
        else:
            fn_name = _utils.get_callback_name(retry_state.fn)

        if retry_state.attempt_number > 1:
            logger.log(
                log_level,
            f"Starting call to '{fn_name}', "
            f"this is the {_utils.to_ordinal(retry_state.attempt_number)} time calling it.",
        )

    return log_it
