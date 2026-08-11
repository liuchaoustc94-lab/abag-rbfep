#!/usr/bin/env python
import pytest
import os
from pmx.model import Model


@pytest.fixture(scope='session')
def gf():
    """Gets file in tests/data with correct path.
    Function taken from mdtraj: https://github.com/mdtraj/mdtraj
    """
    test_dir = os.path.dirname(os.path.abspath(__file__))
    def _gf(fn):
        return '{}/data/{}'.format(test_dir, fn)
    return _gf


# ------
# Models
# ------
@pytest.fixture(scope='session')
def empty_model(gf):
    """Returns empty Model"""
    return Model()


@pytest.fixture(scope='session')
def protein_model(gf):
    """Returns Model of a peptide"""
    return Model(gf('peptide.pdb'))


@pytest.fixture(scope='session')
def dna_model(gf):
    """Returns Model of a DNA molecule"""
    return Model(gf('dna.pdb'))


@pytest.fixture(scope='session')
def rna_model(gf):
    """Returns Model of a RNA molecule"""
    return Model(gf('rna.pdb'))
