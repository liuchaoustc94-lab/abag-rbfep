#!/usr/bin/env python

from pmx.utils import natural_sort

def test_natural_sort():
    lst = ['2', '1c', '10b', '5c', '7y', '100', '1']
    out = ['1', '1c', '2', '5c', '7y', '10b', '100']
    assert natural_sort(lst) == out
