from atomict.api import get


def get_model(model_name: str):
    return get(f"api/ml-model/?short_desc={model_name}")
