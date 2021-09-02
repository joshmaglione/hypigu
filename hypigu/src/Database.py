#
#   Copyright 2020 Joshua Maglione 
#
#   Distributed under MIT License
#

from .Globals import __TIME as _time

# Given a list of posets L, a poset P, and an integer n, decide if L has a poset
# isomorphic to P.
def _check(L, P, n):
    if len(L) == 0 or len(L) <= n:
        return False, None
    if L[n].is_isomorphic(P):
        return True, n
    return _check(L, P, n+1) 


class IADatabase():

    def __init__(self):
        self.poset_list = []
        self.gen_func_list = []

    def __repr__(self):
        return "A database indexed by %s posets" % (len(self.poset_list))

    def has_poset(self, P):
        L = self.poset_list
        return _check(L, P, 0)

    def get_gen_func(self, P, style):
        isit, k = self.has_poset(P)
        if not isit:
            return None
        return self.gen_func_list[k][style]

    def save_gen_func(self, P, style, F):
        assert P.rank() > 2
        assert style in ['Igusa', 'skele']
        isit, k = self.has_poset(P)
        if not isit:
            # New poset, so we save it.
            gen_dict = {
                'Igusa' : None,
                'skele' : None
            }
            gen_dict[style] = F
            self.poset_list += [P]
            self.gen_func_list += [gen_dict]
        else:
            if self.gen_func_list[k][style] == None:
                GFL = self.gen_func_list
                GFL[k][style] = F
                self.gen_func_list = GFL


def _initialize_main_DB():
    from sage.all import DiGraph, Poset, var
    import hypigu.src.init_data as init_data
    q = var('q')
    t = var('t')
    Y = var('Y')
    T = var('T')
    DB = IADatabase()

    # A3 arrangement 
    A3 = Poset(DiGraph(init_data.A3_rels))
    A3_Igusa = q**-6*(1-q**-1)*(q**4*(6-5*q+q**2)-4*q**4*(2-q)*t-q**2*(3-7*q+2*q**2)*t**2+q**2*(2-7*q+3*q**2)*t**3-4*q*(1-2*q)*t**4-(1-5*q+6*q**2)*t**5)/((1-q**-1*t)**2*(1-q**-2*t**3)*(1-q**-3*t**6))
    A3_skele = ((1 + 6*Y + 11*Y**2 + 6*Y**3) + (11 + 37*Y + 37*Y**2 + 11*Y**3)*T + (6 + 11*Y + 6*Y**2 + Y**3)*T**2)/((1 - T)**3)
    DB.save_gen_func(A3, 'Igusa', A3_Igusa)
    DB.save_gen_func(A3, 'skele', A3_skele)

    # A4 arrangement 
    A4 = Poset(DiGraph(init_data.A4_rels))
    A4_skele = ((1 + 10*Y + 35*Y**2 + 50*Y**3 + 24*Y**4) + (47 + 260*Y + 505*Y**2 + 400*Y**3 + 108*Y**4)*T + (108 + 400*Y + 505*Y**2 + 260*Y**3 + 47*Y**4)*T**2 + (24 + 50*Y + 35*Y**2 + 10*Y**3 + Y**4)*T**3)/((1 - T)**4)
    A4_Igusa = init_data.A4_Igusa_n(q, t) / init_data.A4_Igusa_d(q, t)
    DB.save_gen_func(A4, 'skele', A4_skele)
    DB.save_gen_func(A4, 'Igusa', A4_Igusa)
    
    # A5 arrangement
    A5 = Poset(DiGraph(init_data.A5_rels))
    A5_skele = (Y**5*T**4 + 197*Y**5*T**3 + 15*Y**4*T**4 + 1268*Y**5*T**2 + 1546*Y**4*T**3 + 85*Y**3*T**4 + 1114*Y**5*T + 7172*Y**4*T**2 + 4670*Y**3*T**3 + 225*Y**2*T**4 + 120*Y**5 + 4493*Y**4*T + 15320*Y**3*T**2 + 6700*Y**2*T**3 + 274*Y*T**4 + 274*Y**4 + 6700*Y**3*T + 15320*Y**2*T**2 + 4493*Y*T**3 + 120*T**4 + 225*Y**3 + 4670*Y**2*T + 7172*Y*T**2 + 1114*T**3 + 85*Y**2 + 1546*Y*T + 1268*T**2 + 15*Y + 197*T + 1)/((1 - T)**5)
    A5_Igusa = init_data.A5_Igusa_n(q, t) / init_data.A5_Igusa_d(q, t)
    DB.save_gen_func(A5, 'skele', A5_skele)
    DB.save_gen_func(A5, 'Igusa', A5_Igusa)

    # B3 arrangement
    B3 = Poset(DiGraph(init_data.B3_rels))
    B3_skele = (Y**3*T**2 + 20*Y**3*T + 9*Y**2*T**2 + 15*Y**3 + 76*Y**2*T + 23*Y*T**2 + 23*Y**2 + 76*Y*T + 15*T**2 + 9*Y + 20*T + 1)/((1 - T)**3)
    B3_Igusa = init_data.B3_Igusa_n(q, t) / init_data.B3_Igusa_d(q, t)
    DB.save_gen_func(B3, 'skele', B3_skele)
    DB.save_gen_func(B3, 'Igusa', B3_Igusa)

    return DB

global internal_database
internal_database = _initialize_main_DB()
