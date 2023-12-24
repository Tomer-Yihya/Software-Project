from setuptools import setup, find_packages, Extension

setup(
    name='mykmeanssp',
    version="1.0.0",
    install_requires=["invoke"],
    packages=find_packages(),
    ext_modules=[
        Extension("mykmeanssp", ["kmeans.c"])
    ]
)
