import os

from confluent_kafka import Consumer

TOPIC_CANTERA_SIMULATION = "procsim_cantera_simulate"
TOPIC_CATALYSIS_EXPLORE = "catalysis_explore"
TOPIC_CATALYSIS_SIMULATE = "catalysis_simulate"
TOPIC_COUNTERFACTUAL_SIMULATE = "counterfactual_simulate"
TOPIC_MOL_FINGERPRINT = "mol_fingerprint"

GROUP_ID_CANTERA_SIMULATION = "procsim-cantera-simulation"
GROUP_ID_CATALYSIS_EXPLORE = "catalysis-exploration"
GROUP_ID_MOL_FINGERPRINT = "mol-fingerprint"


def get_consumer():
    host = os.environ.get("AT_CONFLUENT_HOST")
    port = os.environ.get("AT_CONFLUENT_PORT")
    app_id = os.environ.get("AT_CONFLUENT_APP_ID")

    c = Consumer(
        {
            "bootstrap.servers": f"{host}:{port}",
            "group.id": app_id,
            "auto.offset.reset": "earliest",
        }
    )

    return c
