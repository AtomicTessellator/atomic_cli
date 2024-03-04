import json

from atomict.api import get, post


def get_simulation(simulation_id: str):
    return get(f"api/cantera-simulation/{simulation_id}/")


def post_cantera_simulation(cantera_simulation: dict):
    cantera_simulation = {
        "code": {
            "name": "Combust 5k High O2",
            "chemical_compositions": "CH4:1, O2:3, AR:1",
            "temperature": 5000,
            "pressure": 101325,
        },
        "project": "f8e4311a-11d1-43c9-a1e9-07997bf531b1",
        "draft": False,
        "simulation_update_id": None,
        "simulation_fork_id": "fa124169-3a3d-43a5-9f66-01242fff4af6",
    }

    return post("api/cantera-simulation/", cantera_simulation)


def post_ct_reaction(ct_reaction: dict):

    for json_prop in ['input_data', 'reactants', 'products']:

        if isinstance(ct_reaction[json_prop], dict):
            ct_reaction[json_prop] = json.dumps(ct_reaction[json_prop])

    return post("api/ct-reaction/", ct_reaction)


def post_ct_species(ct_species: dict):

    for json_prop in ["input_data", ]:

        if isinstance(ct_species[json_prop], dict):
            ct_species[json_prop] = json.dumps(ct_species[json_prop])

    return post("api/ct-species/", ct_species)


def post_ct_observation(ct_observation: dict):
    return post("api/ct-observation/", ct_observation)


def post_ct_thermo(ct_thermo: dict):
    return post("api/ct-thermo/", ct_thermo)
