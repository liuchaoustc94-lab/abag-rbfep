from dataclasses import dataclass

from pmx.geometry import _common_named_atoms


@dataclass
class DummyAtom:
    name: str
    x: list[float]


def test_common_named_atoms_ignores_missing_backbone_hydrogens() -> None:
    atoms1 = [
        DummyAtom("N", [0.0, 0.0, 0.0]),
        DummyAtom("CA", [0.0, 0.0, 0.0]),
        DummyAtom("C", [0.0, 0.0, 0.0]),
        DummyAtom("O", [0.0, 0.0, 0.0]),
    ]
    atoms2 = [
        DummyAtom("N", [1.0, 1.0, 1.0]),
        DummyAtom("CA", [1.0, 1.0, 1.0]),
        DummyAtom("C", [1.0, 1.0, 1.0]),
        DummyAtom("H", [1.0, 1.0, 1.0]),
        DummyAtom("O", [1.0, 1.0, 1.0]),
        DummyAtom("HA", [1.0, 1.0, 1.0]),
    ]

    pairs = _common_named_atoms(atoms1, atoms2, [atom.name for atom in atoms2])
    assert [left.name for left, _ in pairs] == ["N", "CA", "C", "O"]
    assert [right.name for _, right in pairs] == ["N", "CA", "C", "O"]
