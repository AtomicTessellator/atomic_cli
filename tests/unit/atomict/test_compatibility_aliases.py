from unittest.mock import patch

from atomict.simulation.ea import (
    create_ea_exploration_sample,
    get_soec_exploration,
    list_soec_explorations,
)
from atomict.simulation.fhi_aims import (
    get_fhiaims_simulation,
    list_fhiaims_simulations,
)
from atomict.simulation.mlrelax import get_mlrelaxation


@patch("atomict.simulation.mlrelax.get_mlrelax")
def test_mlrelaxation_alias_delegates(mock_get_mlrelax):
    mock_get_mlrelax.return_value = {"id": "ml-1"}

    result = get_mlrelaxation("ml-1", include_ht=True)

    assert result == {"id": "ml-1"}
    mock_get_mlrelax.assert_called_once_with("ml-1", include_ht=True)


@patch("atomict.simulation.fhi_aims.get_simulation")
def test_fhiaims_get_alias_delegates(mock_get_simulation):
    mock_get_simulation.return_value = {"id": "sim-1"}

    result = get_fhiaims_simulation("sim-1", include_ht=True)

    assert result == {"id": "sim-1"}
    mock_get_simulation.assert_called_once_with("sim-1", include_ht=True)


@patch("atomict.simulation.fhi_aims.list_simulations")
def test_fhiaims_list_alias_delegates(mock_list_simulations):
    mock_list_simulations.return_value = {"results": []}

    result = list_fhiaims_simulations(search="alpha")

    assert result == {"results": []}
    mock_list_simulations.assert_called_once_with(search="alpha")


@patch("atomict.simulation.ea.get_ea_exploration")
def test_soec_get_alias_delegates(mock_get_exploration):
    mock_get_exploration.return_value = {"id": "ea-1"}

    result = get_soec_exploration("ea-1", full=True)

    assert result == {"id": "ea-1"}
    mock_get_exploration.assert_called_once_with("ea-1", full=True)


@patch("atomict.simulation.ea.list_ea_explorations")
def test_soec_list_alias_delegates(mock_list_explorations):
    mock_list_explorations.return_value = {"results": []}

    result = list_soec_explorations(ordering="-created_at")

    assert result == {"results": []}
    mock_list_explorations.assert_called_once_with(ordering="-created_at")


@patch("atomict.simulation.ea.create_exploration_sample")
def test_ea_sample_create_alias_delegates(mock_create_sample):
    mock_create_sample.return_value = {"id": "sample-1"}

    result = create_ea_exploration_sample("ea-1", simulation_id="sim-1")

    assert result == {"id": "sample-1"}
    mock_create_sample.assert_called_once_with(
        exploration_id="ea-1",
        simulation_id="sim-1",
        mlrelax_id=None,
        strain=None,
        matrix=None,
    )
