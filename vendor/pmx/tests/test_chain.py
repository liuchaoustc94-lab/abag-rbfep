#!/usr/bin/env python
from pmx.chain import Chain
from filecmp import cmp


def test_create(gf, tmpdir):
    ref_c = gf('chain_create_FRTLKNCWQ.pdb')
    new_c = ''
    new_c = str(tmpdir.mkdir("test_chain").join("out.pdb"))
    c = Chain().create("FRTLKNCWQ")
    c.write(new_c)
    cmp(ref_c, new_c)
