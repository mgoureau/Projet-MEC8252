import numpy as np
import cantera as ct
import time

__all__ = ["calc_frac_mol", "format_X_cantera", "calc_thermoprop_mix"]

def calc_frac_mol(N):
    """
    Renvoie une liste contenant les fractions molaires du mélange
    
    N: liste des quantités de matière de chaque espèce dans le mélange
    """
    
    return N/np.sum(N)
    

def format_X_cantera(species, X):
    """
    Renvoie une string contenant les fractions molaires du mélange formatée pour cantera
    
    species : liste des strings des espèces dans le mélange
    N: liste des quantités de matière de chaque espèce
    """
    for i, sp in enumerate(species):
        if i==0:
            X_cant = "{}:{}".format(species[i], X[i])
        else:
            X_cant = X_cant + ", {}:{}".format(species[i], X[i])

    return X_cant

def calc_thermoprop_mix(gas, species, N, T, P):
    """
    Renvoie l'énergie interne (U), l'enthalpie (H), l'entropie (S) et l'enthalpie  libre (G) du mélange
    
    gas : objet cantera représentant le mélange étudié
    species : liste des strings des espèces dans le mélange
    N : liste des quantités de matière de chaque espèce dans le mélange
    T : température du mélange
    P : pression du mélange
    """
    
    # Calcul des fraction molaires
    X = calc_frac_mol(N)
    
    # Formattage pour cantera
    X_cant = format_X_cantera(species, X)
    
    # Initialisation des paramètres du mélange
    gas.TPX = T, P, X_cant
    
    # Calcul des variables thermochimiques (penser aux unités)
    N_tot = np.sum(N)
    U = gas.u*1e-3*N_tot
    H = gas.h*1e-3*N_tot
    S = gas.s*1e-3*N_tot
    G = gas.g*1e-3*N_tot
    
    return [U, H, S, G]