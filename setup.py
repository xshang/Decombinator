#!/usr/bin/env python
import os
import sys
import pdb
import pkg_resources

from ez_setup import use_setuptools
use_setuptools()

from platform import system, processor
from setuptools import setup

classifiers = ["Natural Language :: English",
               "Programming Language :: Python"]

entry_points = """
[console_scripts]
decombinator-analyze = decombinator.analyze:main
"""

install_requires = "toolshed python-Levenshtein biopython \
                    numpy matplotlib".split()

arch = "_".join([system(), processor()])

class InstallationError(Exception):
    pass

if __name__ == "__main__":
    setup(name="decombinator",
          version="0.1",
          description="handling UMIs in sequences",
          author="Chain et al.",
          author_email="",
          classifiers=classifiers,
          long_description=open('README.md').read(),
          install_requires=install_requires,
          zip_safe=False,
          # XXX: this should be based off of __file__ instead
          packages=['decombinator'],
          include_package_data = True,
          entry_points=entry_points
          )