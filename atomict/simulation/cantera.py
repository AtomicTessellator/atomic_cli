from atomict.api import get, post


def get_simulation(simulation_id: str):
    return get(f"api/cantera-simulation/{simulation_id}/")


def post_ct_reaction(ct_reaction: dict):
    return post("api/ct-reaction/", ct_reaction)


def post_ct_species(ct_species: dict):
    return post("api/ct-species/", ct_species)


def post_ct_observation(ct_observation: dict):
    return post("api/ct-observation/", ct_observation)


def post_ct_thermo(ct_thermo: dict):
    return post("api/ct-thermo/", ct_thermo)
