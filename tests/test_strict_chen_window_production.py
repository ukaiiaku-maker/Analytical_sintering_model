import strict_chen_window_production as s

def test_tier_a_rules():
    assert s.tier(.05,.05,25,300)=='Tier_A'
    assert s.tier(.051,.05,25,300)!='Tier_A'

def test_tier_b_rules():
    assert s.tier(.10,.10,25,450)=='Tier_B'
    assert s.tier(.10,.10,24,450)=='Tier_C'

def test_relaxed_and_reject_rules():
    assert s.tier(.20,.10,100,600)=='Tier_C'
    assert s.tier(.01,.01,40,200,lower=False)=='reject'
    assert s.tier(.01,.01,0,200)=='reject'

def test_practical_cutoff_required():
    assert s.tier(.01,.01,40,200,practical=False)=='reject'
