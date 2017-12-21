import aotools.turbulence

def test_calcSeeing():
    res=aotools.turbulence.calcSeeing(0.135,800.,30)
    assert type(res)==type(float)
