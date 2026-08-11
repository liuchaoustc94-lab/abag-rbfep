#!/usr/bin/env python
import pytest
from pmx.model import Model, double_box
from filecmp import cmp


def test_title(empty_model):
    assert empty_model.title == 'PMX MODEL'


@pytest.mark.parametrize("inputpdb, expected_type", [
    ('peptide.pdb', 'protein'),
    ('dna.pdb', 'dna'),
    ('rna.pdb', 'rna')
])
def test_model_moltype(gf, inputpdb, expected_type):
    m = Model(gf(inputpdb))
    assert m.moltype == expected_type


def test_protein_residue_names(protein_model):
    residues = [r.resname for r in protein_model.residues]
    true_residues = ['NASN', 'LEU', 'TYR', 'ILE', 'GLN', 'TRP', 'LEU', 'LYS',
                     'ASP', 'GLY', 'GLY', 'PRO',  'SER', 'SER', 'GLY', 'ARG',
                     'PRO', 'PRO', 'PRO', 'CSER']
    assert residues == true_residues


@pytest.mark.parametrize("r, d, longest_axis, target_file", [
    (2.5, 1.5, False, 'r2.5_d1.5.pdb'),
    (2.5, 1.5, True,  'r2.5_d1.5_longaxis.pdb'),
    (3.0, 2.0, False, 'r3_d2.pdb'),
    (3.0, 2.0, True,  'r3_d2_longaxis.pdb'),
])
def test_double_box(r, d, longest_axis, target_file, gf, tmpdir):
    f1 = gf('doublebox/file1.gro')
    f2 = gf('doublebox/file2.gro')
    ref = gf('doublebox/{}'.format(target_file))

    m1 = Model(f1, bPDBTER=True, renumber_residues=False)
    m2 = Model(f2, bPDBTER=True, renumber_residues=False)
    mout = double_box(m1, m2, r=2.5, d=1.5, bLongestAxis=False, verbose=False)
    outfile = str(tmpdir.mkdir("test_model").join("doublebox.pdb"))
    mout.write(outfile, bPDBTER=True)

    cmp(ref, outfile)
