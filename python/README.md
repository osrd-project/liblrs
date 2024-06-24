# LibLRS Python

Expose [liblrs](https://github.com/osrd-project/liblrs/) through [pyO3](https://pyo3.rs).

## Usage

Simply add `liblrs-python` to your dependencies.

```python
import liblrs_python as lrs

plm = lrs.Lrs(open("path_to_your_lrs_file", "rb").read())

# We build a dict mapping for each LRM id to its handle
lrms = {plm.get_lrm_scale_id(i): i for i in range(plm.lrm_len())}

# We find the handle for the Via Aurelia (https://en.wikipedia.org/wiki/Via_Aurelia)
via_aurelia_handle = [v for k,v in lrms.items() if k.startswith("Via Aurelia")][0]

# We define two measures meaning “100 passus after milestone 50”
a = lrs.LrmScaleMeasure("50", 100)
b = lrs.LrmScaleMeasure("60", 200)

# Get the coordinates between those two measures
coordinates = [[p.x, p.y] for p in plm.resolve_range(via_aurelia_handle, a, b)]
```

## Developpment

Create your virtualenv and install [maturin](https://www.maturin.rs):

```sh
python -m venv venv-liblrs
source venv-liblrs/bin/activate
pip install
maturin develop
```


https://gist.github.com/Tristramg/32a2fff35eb8e0bb065eed25d987a2c7

Bindings with PyO3 to liblrs

## Publishing

We publish the library to PYPI running:

```sh
docker run --rm -e MATURIN_PYPI_TOKEN=PYPI_key -w /io/python -v $(pwd)/..:/io ghcr.io/pyo3/maturin publish
```