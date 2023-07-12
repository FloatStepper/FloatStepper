import os
import pytest
import oftest
import pandas as pd
import matplotlib.pyplot as plt
from oftest import run_case, clean_case




def test_completed(run_case):
    log = oftest.path_log()
    assert oftest.case_status(log) == 'completed' # checks if run completes


def test_load_data():
    dirname = oftest.base_dir()
    inertTensor = oftest.read_functionObject(os.path.join(dirname,"ma.dat"))
    t = inertTensor.values
    print(t)
    assert True

def test_clean(clean_case):
    assert True
