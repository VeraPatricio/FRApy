

this should call lensed_distance_map


create basic paramter file

def create_parameter():
    """This creates a basic parameter list, with the paramaters of the disctance map.
    Other parameters should be added to fit models"""
    p = Parameters()
    p.add_many(
    #  (Name,               Value,  Vary,    Min,     Max)
       ("cx",                 28,   False,     27,     31.), 
       ("cy",                 22,   False,     21,     25.),
       ("ellip",             0.8,   False,   0.01,    0.99),
       ("theta",               0,  False,    -90,     90))

    return p
