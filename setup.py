from setuptools import setup

setup(
    name="wind_forecast",
    version="0.0.1",
    packages=["wind_forecast"],
    install_requires=["numpy","matplotlib","cma","scipy","outflowpy","appdirs","h5py","pyhdf", "numba","httplib2","h5netcdf","dtaidistance","cdflib"]
)
