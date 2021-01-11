## Importations, paramètres et fonctions

import numpy as np
from math import *
import matplotlib.pyplot as plt
#Paramètres géométriques
h1 = 3
h3 = 1.5
I1 = 7
I2 = 24
I3 = 75
I4 = 15
w1 = 2
w2 = 3
w3 = 10
w4 = 4
D0e = 13.65
D0i = 10.5
D1e = 5.33
D1i = 3.9
L = 110
Lmax = 20
B = 22.8
CI0 = D0e/2+D0i/2
CI1 = D1e/2+D1i/2
r = 4

#Paramètres de discrétisation
dZ = 0.5
eps = 0.5

#Paramètres du fluide/de la machine
eta = 1.00E-3
Q = 0.15E3
Q1 = 0.1E3
Q2 = 0.4E3
m = 3
rho = 0.9E-3
N = 200
e = 10.32



def De(Z):
    #Diamètre extérieur
    return D0e + Z*(D1e - D0e)/L

def Di(Z):
    #Diamètre intérieur
    return D0i + Z*(D1i - D0i)/L

def theta(Z):
    #angle de filet
    return atan(B/(pi*De(Z)))

def W(Z):
    #largeur du chenal
    return B*cos(theta(Z)) -e

def H(Z):
    #Profondeur de chenal
    return (De(Z)-Di(Z))/2

def V(Z):
    #Vitesse relative du fourreau Vf
    return (pi*N/60)*De(Z)*cos(theta(Z))

def dz(Z):
    #Incrément dz obtenu pour un incrément dZ en Z
    return dZ/sin(theta(Z))

def Sl(Z):
    #Surface libre
    C = CI0 + (CI1-CI0)*Z/L
    phi = acos(C/De(Z))
    alpha = pi/2 - 2*phi
    Sf = (pi - phi) * (De(Z)**2)/2 + (C*De(Z)*sin(phi))/2
    Sv = 2 *(phi*C**2 - C*De(Z)*sin(phi)/2) + alpha * ((De(Z)**2)/4 +(C-De(Z)/2)**2)
    return Sf, Sv, Sf - 2*Sv


def dp(Z, q=Q):
    #correspond à (dp(Z)/dz)*dz(Z) (delta de pression en Z)
    dpdz = 12*eta*(V(Z)*W(Z)*H(Z)/2 - q/3)/(W(Z)*H(Z)**3)
    return dpdz*dz(Z)


def fr(Z, q=Q):
    #Fonction à annuler : deltapcanal - deltapchenal(Z), quand on est nul on est à Z correspondant à la surface libre
    if Z < L:
        return fr(Z+dZ, q) - dp(Z, q)
    else :
        return deltap_canal - dp(L, q)

def h(q):
    #Fonction à minimiser : (L-Lrempli) - Z, quand on est nul le débit q donne le bon remplissage.
    P = deltap_canal - dp(L, q)
    c=0
    a, b = Lmax, L
    m = (a+b)/2
    while abs(a-b)>eps/40 and c<2000 :
        fa,fb,fm = fr(a,q), fr(b, q), fr(m, q)
        if fa*fb > 0:
            print("pbm")
        if fa*fm > 0:
            a = m
        else :
            b = m
        m = (a+b)/2
    return L-m-L_rempli, m, fm

## Affichage géométrie
DZ = np.arange(0, 110, dZ)
DDe = [De(i) for i in DZ]
DDi = [Di(i) for i in DZ]
Dtheta = [theta(i) for i in DZ]
DW = [W(i) for i in DZ]
DH = [H(i) for i in DZ]
DV = [V(i) for i in DZ]
Ddz = [dz(i) for i in DZ]
DSl = [Sl(i) for i in DZ]

plt.figure()
plt.subplot(2,3,1)
plt.plot(DZ,DDe, label = "Diamètre extérieur (mm)")
plt.plot(DZ,DDi, label ="Diamètre intérieur (mm)")
plt.legend()
plt.xlabel("Z(mm)")
plt.title("Evolution diamètres extérieur et intérieur")
plt.grid()
plt.subplot(2,3,2)
plt.plot(DZ,DW, label = "Largeur du chenal W (mm)")
plt.plot(DZ,DH, label ="Profondeur du chenal H (mm)")
plt.legend()
plt.xlabel("Z(mm)")
plt.title("Evolution largeur et profondeur du chenal")
plt.grid()
plt.subplot(2,3,3)
plt.plot(DZ,Dtheta, label ="Theta (rad)")
plt.legend()
plt.xlabel("Z(mm)")
plt.title("Evolution angle de filet")
plt.grid()
plt.subplot(2,3,4)
plt.plot(DZ,DV, label ="Vf (mm/s)")
plt.legend()
plt.xlabel("Z(mm)")
plt.title("Evolution vitesse relative du fourreau")
plt.grid()
plt.subplot(2,3,5)
plt.plot(DZ,DV, label ="Sl (mm²)")
plt.legend()
plt.xlabel("Z(mm)")
plt.title("Evolution surface libre")
plt.grid()
plt.subplot(2,3,6)
plt.plot(DZ,DV, label ="dz (mmm)")
plt.legend()
plt.xlabel("Z(mm)")
plt.title("Evolution incrément dz")
plt.grid()
plt.show()

## Détermination du débit

#Détermination de L_rempli
V_rempli = m/rho - ((w1*(2*I1) + w2*I2)*h1 + 5*h3*(w2 + (w3-w2)/2) + 5*(h1-h3)*w2/2 + (w3 - w2)*5*(h1-h3)/2 + h3*(w3*I3 - 2*pi*r**2 + (w4+(w3-w4)/2)*5 + w4*I4))
L_rempli, Vcalc, Z = L, 0, 0

while Z < L :
    Vcalc += Sl(Z)[2]*dZ
    Z+=dZ
while Vcalc > V_rempli :
    Vcalc -= Sl(L-L_rempli)[2]*dZ
    L_rempli -= dZ


#Calcul de deltap_canal
deltap_canal = 12*eta*Q*((I1/(2*w1) + I2/w2)/h1**3 + (I3/w3 + I4/w4)/h3**3)


#Détermination de la pression Q
compt = 0
sol = h(Q)[0]
LM, LH, LQ = [],[],[]
plt.figure()
GRAD = []
qa, qb = Q1, Q2
qm = (qa + qb)/2
plt.figure()
while abs(sol)>eps and compt <100 :
    compt+=1
    hqm = h(qm)
    ha,hb,hm = h(qa)[0], h(qb)[0], hqm[0]
    if ha*hb > 0:
        print("pbm")
    if ha*hm > 0:
        qa = qm
    else :
        qb = qm
    qm = (qa+qb)/2
    sol = hm
    LM += [hqm[1]]
    LH+=[hqm[0]]
    LQ += [qm]
    if compt%1 == 0 :
        DP = [fr(i,qm) for i in DZ]
        plt.plot(DZ,DP, label = "p (MPa) pour Q = {}".format(qm))
    print(compt)
plt.plot([L-L_rempli]*2,[0,20],label ="Longueur remplie à obtenir")
plt.ylim(0,16)
plt.grid()
plt.legend()
plt.xlabel("Z (mm)")
plt.title("Evolution pression le long des vis pour différents débits")
plt.figure()
plt.scatter(LQ, LM)
plt.scatter(LQ, LH)
plt.show()
print("Débit obtenu :", qm)

