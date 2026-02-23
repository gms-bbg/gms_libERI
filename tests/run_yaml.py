"""
Runs arbitrary YAML commands
"""

import os
import sys
import yaml
import pytest


def parse_yaml(filename, key):
    """
    Opens a YAML file.
    """

    with open(filename, "r") as infile:
        ret = yaml.load(infile)

    return ret[key]
