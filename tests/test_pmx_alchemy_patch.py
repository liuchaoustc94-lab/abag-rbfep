from dataclasses import dataclass

import pmx.alchemy as alchemy


@dataclass
class DummyAtom:
    name: str
    x: list[float]
    id: int | None = None
    resnr: int | str | None = None
    atomtypeB: str | None = None
    qB: float | None = None
    mB: float | None = None
    molecule: object | None = None

    def dihedral(self, atom2, atom3, atom4):
        return 0.0


class DummyChain:
    def __init__(self):
        self._residues = {}

    def add(self, residue):
        self._residues[residue.id] = residue
        residue.chain = self

    def fetch_residue(self, idx):
        return self._residues.get(idx)


class DummyResidue:
    def __init__(self, atoms, real_resname="VAL", resname="VAL", hybrid=False, residue_id=1):
        self.atoms = atoms
        self.real_resname = real_resname
        self.resname = resname
        self._by_name = {atom.name: atom for atom in atoms}
        self.hybrid = hybrid
        self.id = residue_id
        self.chain = None
        for atom in self.atoms:
            atom.molecule = self

    def get_real_resname(self):
        return self.real_resname

    def fetch(self, name):
        atom = self._by_name.get(name)
        return [atom] if atom is not None else []

    def fetchm(self, names):
        return [self._by_name[name] for name in names if name in self._by_name]

    def __getitem__(self, name):
        return self._by_name[name]

    def is_hybrid(self):
        return self.hybrid


def test_set_conformation_skips_missing_old_residue_atoms(monkeypatch) -> None:
    monkeypatch.setitem(alchemy.library._aa_dihedrals, "VAL", [["N", "CA", "CB", "HB", 1, 1]])

    old_res = DummyResidue(
        [
            DummyAtom("N", [0.0, 0.0, 0.0]),
            DummyAtom("CA", [1.0, 0.0, 0.0]),
            DummyAtom("CB", [1.0, 1.0, 0.0]),
            DummyAtom("C", [2.0, 0.0, 0.0]),
            DummyAtom("O", [3.0, 0.0, 0.0]),
        ]
    )
    new_res = DummyResidue(
        [
            DummyAtom("N", [9.0, 9.0, 9.0]),
            DummyAtom("CA", [9.0, 9.0, 9.0]),
            DummyAtom("CB", [9.0, 9.0, 9.0]),
            DummyAtom("HB", [9.0, 9.0, 9.0]),
            DummyAtom("CG1", [9.0, 9.0, 9.0]),
            DummyAtom("C", [9.0, 9.0, 9.0]),
            DummyAtom("O", [9.0, 9.0, 9.0]),
            DummyAtom("H", [9.0, 9.0, 9.0]),
        ]
    )

    alchemy._set_conformation(old_res, new_res, {"CA-CB": ["HB", "CG1"]})

    assert new_res["N"].x == [0.0, 0.0, 0.0]
    assert new_res["CA"].x == [1.0, 0.0, 0.0]
    assert new_res["C"].x == [2.0, 0.0, 0.0]
    assert new_res["H"].x == [9.0, 9.0, 9.0]


def test_get_hybrid_residues_matches_atoms_by_name(monkeypatch) -> None:
    residue = DummyResidue(
        [
            DummyAtom("N", [0.0, 0.0, 0.0]),
            DummyAtom("CA", [1.0, 0.0, 0.0]),
        ],
        real_resname="Y2F",
        resname="Y2F",
        hybrid=True,
    )
    hybrid_residue = DummyResidue(
        [
            DummyAtom("N", [0.0, 0.0, 0.0], atomtypeB="N_B", qB=1.0, mB=14.0),
            DummyAtom("CA", [1.0, 0.0, 0.0], atomtypeB="CA_B", qB=2.0, mB=12.0),
            DummyAtom("H", [2.0, 0.0, 0.0], atomtypeB="H_B", qB=3.0, mB=1.0),
        ],
        real_resname="Y2F",
        resname="Y2F",
        hybrid=True,
    )

    monkeypatch.setattr(alchemy, "get_mtp_file", lambda *args, **kwargs: "dummy.mtp")
    monkeypatch.setattr(
        alchemy,
        "_get_hybrid_residue",
        lambda **kwargs: (hybrid_residue, [], [], [], {}),
    )

    class DummyModel:
        residues = [residue]

    rlist, _ = alchemy._get_hybrid_residues(DummyModel(), ff="dummy")
    assert rlist == [residue]
    assert residue["N"].atomtypeB == "N_B"
    assert residue["CA"].atomtypeB == "CA_B"


def test_resolve_dihedral_atoms_returns_none_when_referenced_atom_is_missing() -> None:
    chain = DummyChain()
    prev_residue = DummyResidue([DummyAtom("C", [0.0, 0.0, 0.0], id=1, resnr=1)], residue_id=1)
    residue = DummyResidue(
        [
            DummyAtom("N", [1.0, 0.0, 0.0], id=2, resnr=2),
            DummyAtom("CA", [2.0, 0.0, 0.0], id=3, resnr=2),
        ],
        residue_id=2,
    )
    next_residue = DummyResidue([DummyAtom("N", [3.0, 0.0, 0.0], id=4, resnr=3)], residue_id=3)
    for entry in (prev_residue, residue, next_residue):
        chain.add(entry)

    assert alchemy._resolve_dihedral_atoms(residue, ["-C", "N", "CA", "+CB"]) is None


def test_check_dih_ildn_opls_skips_templates_with_missing_atoms() -> None:
    chain = DummyChain()
    residue = DummyResidue(
        [
            DummyAtom("N", [0.0, 0.0, 0.0], id=10, resnr=2),
            DummyAtom("CA", [1.0, 0.0, 0.0], id=11, resnr=2),
            DummyAtom("CB", [1.0, 1.0, 0.0], id=12, resnr=2),
            DummyAtom("CG", [1.0, 1.0, 1.0], id=13, resnr=2),
        ],
        resname="VAL",
        residue_id=2,
    )
    chain.add(residue)
    rdic = {"VAL": [None, None, [], [["N", "CA", "CB", "HB", "torsion_1", "un"]]]}

    counter = alchemy._check_dih_ILDN_OPLS(
        topol=None,
        rlist=[residue],
        rdic=rdic,
        a1=residue["N"],
        a2=DummyAtom("CA", [0.0, 0.0, 0.0], id=11, resnr=2),
        a3=residue["CB"],
        a4=residue["CG"],
    )

    assert counter == 0


def test_check_dih_ildn_opls_matches_insertion_coded_residue_atoms_via_molecule() -> None:
    chain = DummyChain()
    residue = DummyResidue(
        [
            DummyAtom("N", [0.0, 0.0, 0.0], id=10, resnr="52A"),
            DummyAtom("CA", [1.0, 0.0, 0.0], id=11, resnr="52A"),
            DummyAtom("CB", [1.0, 1.0, 0.0], id=12, resnr="52A"),
            DummyAtom("CG", [1.0, 1.0, 1.0], id=13, resnr="52A"),
        ],
        resname="VAL",
        residue_id=1,
    )
    chain.add(residue)
    rdic = {"VAL": [None, None, [], [["N", "CA", "CB", "CG", "un", "un"]]]}

    counter = alchemy._check_dih_ILDN_OPLS(
        topol=None,
        rlist=[residue],
        rdic=rdic,
        a1=residue["N"],
        a2=residue["CA"],
        a3=residue["CB"],
        a4=residue["CG"],
    )

    assert counter == 1
