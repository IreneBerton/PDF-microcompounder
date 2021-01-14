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
# D0e = 14
# D0i = 11
# D1e = 6
# D1i = 4
L = 110
Lmax = 20
B = 22.8
CI0 = D0e/2+D0i/2
CI1 = D1e/2+D1i/2
# CI0 = 12.5
# CI1 = 5
r = 3.5

#Paramètres de discrétisation
pdZ = 0.5
peps = 0.1

#Paramètres du fluide/de la machine
peta = 1.00E-3
pQ = 0.15E3
pQ1 = 0.1E3
pQ2 = 0.65E3
pm = 3
prho = 0.9E-3
pN = 200
e = 10.32


###

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


def debit(m = pm, N = pN, eta = peta,Q1 = pQ1,rho = prho, dZ = pdZ, eps = peps, trace = True):
    pN = N
    pQ1 = Q1
    prho = rho
    pm = m
    peta = eta
    pdZ, peps = dZ, eps

    def V(Z):
        #Vitesse relative du fourreau Vf
        return (pi*pN/60)*De(Z)*cos(theta(Z))

    def dz(Z):
        #Incrément dz obtenu pour un incrément dZ en Z
        return pdZ/sin(theta(Z))

    def Sl(Z):
        #Surface libre
        C = CI0 + (CI1-CI0)*Z/L
        print(Z,C,De(Z), C/De(Z))
        phi = acos(C/De(Z))
        alpha = pi/2 - 2*phi
        Sf = (pi - phi) * (De(Z)**2)/2 + (C*De(Z)*sin(phi))/2
        Sv = 2 *(phi*C**2 - C*De(Z)*sin(phi)/2) + alpha * ((De(Z)**2)/4 +(C-De(Z)/2)**2)
        return Sf, Sv, Sf - 2*Sv, phi*180/pi, alpha*180/pi

    def dp(Z, q=pQ):
        #correspond à (dp(Z)/dz)*dz(Z) (delta de pression en Z)
        dpdz = 12*peta*(V(Z)*W(Z)*H(Z)/2 - q/3)/(W(Z)*H(Z)**3)
        return dpdz*dz(Z)

    def deltap_canal(q=pQ):
        return 12*peta*q*((I1/(2*w1) + I2/w2)/h1**3 + (I3/w3 + I4/w4)/h3**3)

    def fr(Z, q=pQ):
        #Fonction à annuler : deltapcanal - deltapchenal(Z), quand on est nul on est à Z correspondant à la surface libre
        if Z < L:
            return fr(Z+pdZ, q) - dp(Z, q)
        else :
            return deltap_canal(q) - dp(L, q)

    def h(q):
        #Fonction à minimiser : (L-L_rempli) - Z, quand on est nul le débit q donne le bon remplissage.
        P = deltap_canal(q) - dp(L, q)
        c=0
        a, b = Lmax, L
        m = (a+b)/2
        while abs(a-b)>peps/40 and c<2000 :
            fa,fb,fm = fr(a,q), fr(b, q), fr(m, q)
            if fa*fb > 0:
                print(q,a,b,fa,fb)
                print("pbm dans L_rempli")
            if fa*fm > 0:
                a = m
            else :
                b = m
            m = (a+b)/2
        return L-m-L_rempli, m, fm


    Vcanal = ((w1*(2*I1) + w2*I2)*h1 + 5*h3*(w2 + (w3-w2)/2) + 5*(h1-h3)*w2/2 + (w3 - w2)*5*(h1-h3)/2 + h3*(w3*I3 + (w4+(w3-w4)/2)*5 + w4*I4))
    mmini = Vcanal*prho
    if m < mmini :
        print("masse trop faible ! Il faut m >", mmini)
    DZ = np.arange(0,L,dZ)
    #Détermination de L_rempli
    V_rempli = pm/prho - Vcanal

    V_rempli = pm/prho - 1700
    L_rempli, Vcalc, Z = L, 0, 0

    while Z < L :
        Vcalc += Sl(Z)[2]*dZ
        Z+=dZ
    print(Vcalc)
    while Vcalc > V_rempli :
        Vcalc -= Sl(L-L_rempli)[2]*dZ
        L_rempli -= dZ
    print(V_rempli, L_rempli)

    Q2 = pQ1
    frl = fr(L, Q2)
    while fr(Lmax,Q2)*fr(L,Q2)<0 :
        Q2+=1
    Q2-=1
    print(Q2)

    #Détermination de la pression Q
    compt = 0

    LM, LH, LQ = [],[],[]
    plt.figure()
    GRAD = []
    qa, qb = Q1, Q2
    qm = (qa + qb)/2
    hqm = h(qm)
    sol = hqm[0]
    if trace :
        plt.figure()
    while abs(sol)>eps and compt <100 :
        compt+=1
        print("qm",qm)
        ha,hb,hm = h(qa)[0], h(qb)[0], hqm[0]
        print(qa, qb, ha, hb)
        if ha*hb > 0:
            print("pbm dans l'intervalle de Q")
            return False
        if ha*hm > 0:
            qa = qm
        else :
            qb = qm
        qm = (qa+qb)/2
        hqm = h(qm)
        sol = hqm[0]
        LQ += [qm]
        DP = [fr(i,qm) for i in DZ]
        if trace :
            plt.plot(DZ,DP, label = "p (MPa) pour Q = {}mm^3/s".format(qm))
        print(compt)
    if compt == 100 :
        print("Pas de convergence...")
    if trace :
        plt.plot([L-L_rempli]*2,[0,20],label ="Longueur remplie à obtenir")
        plt.ylim(0,25)
        plt.grid()
        plt.legend()
        plt.xlabel("Z (mm)")
        plt.title("Evolution pression le long des vis, m ={}g, N={}tour/min, eta={}MPa/s, rho={}g/mm^3".format(m,N,eta,rho))
        plt.show()
    return qm, V_rempli, L_rempli, DP

##Affichage débit/pression...

a= debit()
DZ = np.arange(0,L,pdZ)
#q1 = debit(m = 1.4, trace = False)
q2 = debit(m = 2, trace = False)
q3 = debit(m = 3, trace = False)
q4 = debit(m = 3.4, trace = False)
#DP1 = q1[-1]
DP2 = q2[-1]
DP3 = q3[-1]
DP4 = q4[-1]
plt.figure()
#plt.plot(DZ, DP1, label = "Q = {}mm^3/s".format(q1[0]))
plt.plot(DZ, DP2, color = 'b',label = "Q = {}mm^3/s, m = 2 g".format(q2[0]))
plt.plot(DZ, DP3, color = 'r',label = "Q = {}mm^3/s, m = 3 g".format(q3[0]))
plt.plot(DZ, DP4, color='g',label = "Q = {}mm^3/s, m = 3.4 g".format(q4[0]))
#plt.plot([L-q1[2]]*2,[0,20],label ="L_remplie pour Q = {}mm^3/s".format(q1[0]))
plt.plot([L-q2[2]]*2,[-2,2],color = 'b',label ="L_remplie pour m = {}g".format(2))
plt.plot([L-q3[2]]*2,[-2,2],color = 'r', label ="L_remplie pour m = {}g".format(3))
plt.plot([L-q4[2]]*2,[-2,2],color = 'g',label ="L_remplie pour m = {}g".format(3.4))
plt.grid()
plt.xlabel("Z (mm)")
plt.ylabel("pression (MPa)")
plt.legend()
plt.show()

## Affichage géométrie

DDe = [De(i) for i in DZ]
DDi = [Di(i) for i in DZ]
Dtheta = [theta(i)*180/pi for i in DZ]
DW = [W(i) for i in DZ]
DH = [H(i) for i in DZ]
DV = [V(i) for i in DZ]
Ddz = [dz(i) for i in DZ]
DSl = [Sl(i)[2] for i in DZ]
DSf = [Sl(i)[0] for i in DZ]
DSv = [Sl(i)[1] for i in DZ]
Dalpha = [Sl(i)[4] for i in DZ]
Dpsi = [Sl(i)[3] for i in DZ]

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
plt.plot(DZ,Dtheta, label ="Theta (deg)")
plt.plot(DZ,Dalpha, label ="alpha (deg)")
plt.plot(DZ,Dpsi, label ="psi (deg)")
plt.legend()
plt.xlabel("Z(mm)")
plt.title("Evolution angle de filet, alpha et psi")
plt.grid()
plt.subplot(2,3,4)
plt.plot(DZ,DV, label ="Vf (mm/s)")
plt.legend()
plt.xlabel("Z(mm)")
plt.title("Evolution vitesse relative du fourreau")
plt.grid()
plt.subplot(2,3,5)
plt.plot(DZ,DSl, label ="Sl (mm²)")
plt.plot(DZ,DSf, label ="Sf (mm²)")
plt.plot(DZ,DSv, label ="Sv (mm²)")
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
