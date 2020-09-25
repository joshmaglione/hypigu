#
#   Copyright 2020 Joshua Maglione 
#
#   Distributed under MIT License
#


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
    from Globals import __DEFAULT_p, __DEFAULT_t
    import init_data
    p = var(__DEFAULT_p)
    t = var(__DEFAULT_t)
    DB = IADatabase()

    # A3 arrangement 
    A3 = Poset(DiGraph(init_data.A3_rels))
    A3_Igusa = p**-6*(1-p**-1)*(p**4*(6-5*p+p**2)-4*p**4*(2-p)*t-p**2*(3-7*p+2*p**2)*t**2+p**2*(2-7*p+3*p**2)*t**3-4*p*(1-2*p)*t**4-(1-5*p+6*p**2)*t**5)/((1-p**-1*t)**2*(1-p**-2*t**3)*(1-p**-3*t**6))
    A3_skele = ((1 + 6*p + 11*p**2 + 6*p**3) + (11 + 37*p + 37*p**2 + 11*p**3)*t + (6 + 11*p + 6*p**2 + p**3)*t**2)/((1 - t)**3)
    DB.save_gen_func(A3, 'Igusa', A3_Igusa)
    DB.save_gen_func(A3, 'skele', A3_skele)

    # A4 arrangement 
    A4 = Poset(DiGraph(init_data.A4_rels))
    A4_skele = ((1 + 10*p + 35*p**2 + 50*p**3 + 24*p**4) + (47 + 260*p + 505*p**2 + 400*p**3 + 108*p**4)*t + (108 + 400*p + 505*p**2 + 260*p**3 + 47*p**4)*t**2 + (24 + 50*p + 35*p**2 + 10*p**3 + p**4)*t**3)/((1 - t)**4)
    A4_Igusa = init_data.A4_Igusa_n(p, t) / init_data.A4_Igusa_d(p, t)
    DB.save_gen_func(A4, 'skele', A4_skele)
    DB.save_gen_func(A4, 'Igusa', A4_Igusa)
    
    # A5 arrangement
    A5 = Poset(DiGraph(init_data.A5_rels))
    A5_skele = (p**5*t**4 + 197*p**5*t**3 + 15*p**4*t**4 + 1268*p**5*t**2 + 1546*p**4*t**3 + 85*p**3*t**4 + 1114*p**5*t + 7172*p**4*t**2 + 4670*p**3*t**3 + 225*p**2*t**4 + 120*p**5 + 4493*p**4*t + 15320*p**3*t**2 + 6700*p**2*t**3 + 274*p*t**4 + 274*p**4 + 6700*p**3*t + 15320*p**2*t**2 + 4493*p*t**3 + 120*t**4 + 225*p**3 + 4670*p**2*t + 7172*p*t**2 + 1114*t**3 + 85*p**2 + 1546*p*t + 1268*t**2 + 15*p + 197*t + 1)/((1 - t)**5)
    A5_Igusa = init_data.A5_Igusa_n(p, t) / init_data.A5_Igusa_d(p, t)
    DB.save_gen_func(A5, 'skele', A5_skele)
    DB.save_gen_func(A5, 'Igusa', A5_Igusa)

    # B3 arrangement
    B3 = Poset(DiGraph(init_data.B3_rels))
    B3_skele = (p**3*t**2 + 20*p**3*t + 9*p**2*t**2 + 15*p**3 + 76*p**2*t + 23*p*t**2 + 23*p**2 + 76*p*t + 15*t**2 + 9*p + 20*t + 1)/((1 - t)**3)
    B3_Igusa = init_data.B3_Igusa_n(p, t) / init_data.B3_Igusa_d(p, t)
    DB.save_gen_func(B3, 'skele', B3_skele)
    DB.save_gen_func(B3, 'Igusa', B3_Igusa)

    return DB

global internal_database
internal_database = _initialize_main_DB()
