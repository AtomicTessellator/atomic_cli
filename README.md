# Atomic Tessellator - CLI package

## Installation
```
pip install atomict
```


## CLI Usage
```
# Show available public datasets

    > at list_datasets

# Download a dataset

    > at download_dataset tox21
```

## SDK Usage
```
from atomict.datasets import get

df = get('tox21')
```

## Roadmap
- [x] CLI
- [x] Download Public Datasets
- [ ] Local dataset caching
- [ ] Authentication
- [ ] Reality Server API endpoints - Simulations

