New readme for repo VSW_forecasting, designed to avoid all the mess with the original one.

Need to organise scripts etc. so I know where everything is, and to allow for a lot of flexibility without losing track of what is actually being acheived.

First steps will be to have scripts for collecting OMNI data, and calculating outflow/PFSS fields. In all such things, timestamps are necessary, as are lists of the parameters used etc.

Installation:

Currently needs python 3.12 as 3.14 is too new for various of the dependencies. That seems to be fine for now, though

```
python3.12 -m venv .venv
source .venv/bin/activate.csh
pip install -e .
```

Should do all the dependencies and things.
