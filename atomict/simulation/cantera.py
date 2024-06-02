import json

from atomict.api import get, post


def get_simulation(simulation_id: str):
    return get(f"api/cantera-simulation/{simulation_id}/")


def post_reaction(ct_reaction: dict):

    for json_prop in ['input_data', 'reactants', 'products']:

        if isinstance(ct_reaction[json_prop], dict):
            ct_reaction[json_prop] = json.dumps(ct_reaction[json_prop])

    return post("api/ct-reaction/", ct_reaction)


def post_species(ct_species: dict):

    for json_prop in ["input_data", ]:

        if isinstance(ct_species[json_prop], dict):
            ct_species[json_prop] = json.dumps(ct_species[json_prop])

    return post("api/ct-species/", ct_species)


def post_observation(ct_observation: dict):
    return post("api/ct-observation/", ct_observation)


def post_thermo(ct_thermo: dict):
    return post("api/ct-thermo/", ct_thermo)
