import numpy as np

__all__ = ["jacobien", "solveur_Newton_Raphson"]

def jacobien(fun, F0, X):
    """
    Renvoie le jacobien/la dérivée numérique de la fonction vectorielle/scalaire fun
    
    fun : fonction dont le jacobien/la dérivée est à évalué(e)
    F0 : évaluation de la fonction fun en X
    X : point d'évaluation du jacobien/de la dérivée
    """
    if type(X)!=np.ndarray:
        dFdX=0
        eps1=1e-5
        eps2=1e-6
        
        # Perturbation du scalaire X
        X_perturbed = X*(1+eps1)+eps2
        dX = X_perturbed-X
        
        # Evaluation de la fonction au point X perturbe
        F_perturbed=fun(X_perturbed)
        dF = F_perturbed-F0
        
        dFdX=dF/dX
    else:
        dFdX=np.zeros((len(X),len(X)))
        eps1=1e-8
        eps2=1e-10
    
        for i in range(len(X)):
            # Perturbation du vecteur X
            X_perturbed = np.copy(X)
            X_perturbed[i] = X[i]*(1+eps1)+eps2
            dX = X_perturbed-X

            # Evaluation de la fonction au point X perturbe
            F_perturbed=fun(X_perturbed)
            dF = F_perturbed-F0

            # Approximation de la colonne du Jacobien 
            dFdX[:,i]=dF/dX[i]
        
    return dFdX

def solveur_Newton_Raphson(fun, xsol, tol, maxiter):
    """
    Renvoie le zéro de la fonction vectorielle/scalaire fun
    
    fun : fonction dont on doit trouvé le zéro
    xsol : estimation initiale de la solution puis solution correspondant au zéro de fun
    tol : tolérance relative du solveur
    maxiter : nombre maximal d'itération du solveur
    """
    i=0
    flag_iter=1
    verbose=0
    
    residu = 1e12
    F=fun(xsol)
    
    while(residu>tol and flag_iter):
        i=i+1
        if type(xsol)!=np.ndarray:
            # Evaluation du jacobien de fun
            dFdx=jacobien(fun, F, xsol)
            
            # Iteration vers la solution
            xsol = xsol-F/dFdx
            
            # Evaluation de la fonction fun pour l'itération d'après et le résidu
            F=fun(xsol)
            residu=abs(F)
            
            # Détail sur le solveur
            if verbose>0:
                print("Itération T : {}".format(i))
                print("Résidu T : {}".format(residu))
        else:
            # Evaluation du jacobien de fun
            dFdX=jacobien(fun, F, xsol)
            
            # Transformation en structure colonne pour le calcul
            xsol=np.transpose(xsol)
            F=np.transpose(F)

            # Itération vers la solution
            xsol = xsol-np.dot(np.linalg.inv(dFdX),F)
            
            # Borne physique pour aider à la convergence
            for i_x in range(len(xsol)):
                xsol[i_x]=max(1e-15,min(1,xsol[i_x]))
            
            # Evaluation de la fonction fun pour l'itération d'après et le résidu
            F=fun(xsol)
            #print("{} : {}".format(i,F))
            residu=abs(np.sum(np.sqrt(F**2)))
            
            # Détail sur le solveur
            if verbose>1:
                print("Itération chimie : {}".format(i))
                print("Résidu chimie : {}".format(residu))
        
        if i>maxiter:
            flag_iter=0;
            print("Nombre maximum d'itération atteint")
    
    if type(xsol)!=np.ndarray:
        if verbose>0:
            print("Solution à la convergence : {}K\n".format(xsol))
        return xsol
    else:
        if verbose>0:
            print("Nombre d'itération pour la chimie : {}".format(i))
        return np.transpose(xsol)