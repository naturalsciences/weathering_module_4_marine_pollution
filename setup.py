from setuptools import setup, find_packages
setup(
    name='weathering_mar_pol',
    version='0.1',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "matplotlib>=3.10.3",
        "numpy>=2.3.1"
    ],
)