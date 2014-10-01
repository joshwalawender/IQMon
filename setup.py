from setuptools import setup, find_packages
from IQMon import __version__
setup(
    name = "IQMon",
    version=__version__,
    author='Josh Walawender',
    packages = find_packages(),
)
