#!/usr/bin/env python
from setuptools import setup, find_packages

setup(
    name="ptp",
    version="1.1",
    packages=find_packages(exclude=['contrib', 'docs', 'tests']),
    install_requires=[],
    tests_require=[],
    author="Jiajie Zhang",
    author_email="bestzhangjiajie@gmail.com",
    description="",
    keywords="ptp",
    url="",
    zip_safe=True,
    include_package_data=True,
    package_data={
        '': ['*.dat', '*.gctx'],
    }
)
