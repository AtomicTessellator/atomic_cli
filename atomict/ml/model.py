from atomict.api import get


def get_model(model_name: str):
    return get(f"api/ml-model/?short_desc={model_name}")


def get_model_by_id(model_id: str):
    return get(f"api/ml-model/{model_id}/?include_weights=true")


def get_hot_models():
    return get("ml-models/hot/")
