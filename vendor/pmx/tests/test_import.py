import pmx

def test_name():
    try:
        assert pmx.__name__ == 'pmx'
    except Exception as e:
        raise e
