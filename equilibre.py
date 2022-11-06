import numpy as np
import cantera as ct
from math import *

from thermoprop import *
from solveur import *

__all__ = ["calc_LHS", "calc_RHS", "calc_comp_TP_cst", "calc_H_TP_cst"]

def calc_LHS(gas, species, T, phi):
    """
    Renvoie le LHS avec les constantes d'équilibre Kp associées à chaque réaction calculées 
    avec les enthalpies libres et les bilans de conservation des espèces chimiques 
    
    gas : objet cantera représentant le mélange étudié
    species : liste des strings des espèces dans le mélange
    T : température du mélange
    phi : ratio d'équivalence du mélange
    """
    R = ct.gas_constant*1e-3       # Constante universelle des gaz parfaits [J/mol/K]
    P0 = ct.one_atm                # Pression de référence [Pa]

    # Récupération des enthalpie libre de chacune des espèces dans l'équation en [J/mol]
    g = {}
    for sp in species:
        g[sp] = calc_thermoprop_mix(gas, [sp], [1], T, P0)[3]
    
    # Calcul du LHS 
    LHS = # à remplir
    
    return np.array(LHS)

def calc_RHS(X, species, P):
    """
    Renvoie le RHS avec les constantes d'équilibre Kp associées à chaque réaction calculées 
    avec la composition du mélange et les bilans de conservation des espèces chimiques
    
    N: liste des fractions molaires de chaque espèce dans le mélange
    species : liste des strings des espèces dans le mélange
    P : pression dans le mélange
    """
    P0 = ct.one_atm                # Pression de référence [Pa]

    # Calcul du RHS de l'équation d'équilibre chimique (composition chimique)
    Xd = {}
    for i_sp, sp in enumerate(species):
        Xd[sp] = X[i_sp]
    
    # Calcul du RHS
    RHS = # à remplir
    
    return np.array(RHS)

def calc_comp_TP_cst(X_est, gas, species, T, P, phi):
    """
    Renvoie la liste des fractions molaires à T,P constant à l'équilibre du mélange
    
    X_est : liste des fractions molaires des espèces du mélange
    gas : objet cantera représentant le mélange étudié
    T : température du mélange
    species : liste des strings des espèces dans le mélange
    P : pression dans le mélange
    phi : ratio d'équivalence du mélange
    """
    # Paramètres pour le solveur
    epsi_X = 1e-8                     # Tolérance relative
    max_iter = 1000                       # Nombre maximal d'itération pour le solveur
    
    # Calcul du LHS de l'équation d'équilibre chimique (avec l'enthalpie libre)
    LHS = calc_LHS(gas, species, T, phi)
    
    # Appel du solveur pour trouver la composition chimique à T, P constants
    def fun_eq(X):
        RHS = calc_RHS(X, species, P)
        return RHS-LHS
    X_eq = solveur_Newton_Raphson(fun_eq, X_est, np.linalg.norm(LHS,2)**epsi_X, max_iter)
    
    return X_eq

def calc_H_TP_cst(T, X_est, gas, species, P, phi):
    """
    Renvoie l'enthalpie des produits après calcul de l'équilibre des espèces à la température T 
    
    X_est : liste des fractions molaires des espèces du mélange
    gas : objet cantera représentant le mélange étudié
    T : température du mélange
    species : liste des strings des espèces dans le mélange
    P : pression dans le mélange
    phi : ratio d'équivalence du mélange
    """
    X_eq = calc_comp_TP_cst(X_est, gas, species, T, P, phi)
    
    Xd = {}
    for i_sp, sp in enumerate(species):
        Xd[sp] = X_eq[i_sp]
    
    # Quantité de matière totale des produits pour 1 mole de carburant
    N_tot = # à remplir

    # Liste des quantités de matière des produits
    N_P = # à remplir
    
    # Enthalpie des produits dépendante de T
    H_P = calc_thermoprop_mix(gas, species, N_P, T, P)[1]
    
    return H_P