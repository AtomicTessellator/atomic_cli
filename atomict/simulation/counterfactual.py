from atomict.api import get, post


def get_simulation(simulation_id: str):
    return get(f"api/cf-run/{simulation_id}/")


def create_simulation(data: dict):
    return post("api/cf-run/", data=data)


def create_cf_sample_space_example(payload: dict, extra_headers={}):

    return post(
        "api/cf-sample-space-example/",
        payload=payload,
        extra_headers=extra_headers
    )
