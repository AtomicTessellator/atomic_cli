import logging
import os

TOPIC_CANTERA_SIMULATION = "procsim_cantera_simulate"
TOPIC_CATALYSIS_EXPLORE = "catalysis_explore"
TOPIC_CATALYSIS_SIMULATE = "catalysis_simulate"
TOPIC_COUNTERFACTUAL_SIMULATE = "counterfactual_simulate"
TOPIC_MOL_FINGERPRINT = "mol_fingerprint"
TOPIC_KNOWLEDGE_EXTRACT = "knowledge_extract"
TOPIC_RMG_RUN = "rmg_run"


def get_consumer(
    c_group: str,
    auto_commit: bool = True
):
    # This is guarded because it's pending removal
    from confluent_kafka import Consumer

    host = os.environ.get("AT_CONFLUENT_HOST")
    port = os.environ.get("AT_CONFLUENT_PORT")
    c_group = os.environ.get("AT_CONFLUENT_CONSUMER_GROUP")

    session_timeout = os.environ.get("AT_CONFLUENT_SESSION_TIMEOUT_MS")
    heartbeat_interval = os.environ.get("AT_CONFLUENT_HEARTBEAT_INTERVAL_MS")

    try:
        session_timeout = int(session_timeout)
        logging.info(
            f"AT_CONFLUENT_SESSION_TIMEOUT_MS configuration found: {session_timeout}"
        )
    except (ValueError, TypeError):
        logging.warning(
            f"Invalid session timeout: {session_timeout}, using default of 30000"
        )
        session_timeout = 30000

    try:
        heartbeat_interval = int(heartbeat_interval)
        logging.info(
            f"AT_CONFLUENT_HEARTBEAT_INTERVAL_MS configuration found: {heartbeat_interval}"
        )
    except (ValueError, TypeError):
        logging.warning(
            f"Invalid heartbeat interval: {heartbeat_interval}, using default of 10000"
        )
        heartbeat_interval = 10000

    c = Consumer(
        {
            "bootstrap.servers": f"{host}:{port}",
            "group.id": c_group,
            "auto.offset.reset": "earliest",
            "enable.auto.commit": auto_commit,
            "heartbeat.interval.ms": heartbeat_interval,
            "session.timeout.ms": session_timeout,
            "max.poll.interval.ms": session_timeout
        }
    )

    return c
