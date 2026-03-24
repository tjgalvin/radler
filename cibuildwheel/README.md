# Using `cibuildwheel`

The`cibuildwheel` tool is intended to be run on a CI server. However, for testing purposes, it is also possible to run `cibuildwheel` locally. First, make sure that you have `docker` installed on your system, and that you have the rights to start docker containers. Next, create a Python virtual environment and install `cibuildwheel` in it. Then go to the root directory of the project (i.e. the directory containing your `setup.py`, `setup.cfg`, and/or `pyproject.toml` file). Finally, assuming you're on Linux, enter the following command:

```
cibuildwheel --platform linux
```

For more information check the `cibuildwheel --help` output, or goto https://cibuildwheel.readthedocs.io/en/stable/
