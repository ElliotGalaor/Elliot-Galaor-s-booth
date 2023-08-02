# -*- coding: utf-8 -*-
"""
Created on Fri Dec 24 16:29:46 2021

@author: ellio
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.animation as animation
import os
import sqlite3
from mpl_toolkits import mplot3d

clear_table=1

conn = sqlite3.connect('données_mesure.db')
cursor = conn.cursor()

if clear_table==1:
    cursor.execute("""DROP TABLE valeurs;""")

cursor.execute("""
CREATE TABLE IF NOT EXISTS valeurs(
     condition FlOAT,
     temps FLOAT,
     horizontal FLOAT,
     vertical FLOAT,
     pressure FLOAT)
""")


#on créer une nouvelle classe pour les conditions aux limites afin de pouvoir résoudre les équations aux dérivées partielles

class Boundary:
    def __init__(self,boundary_type,boundary_value):
        self.DefineBoundary(boundary_type,boundary_value)
        
    def DefineBoundary(self,boundary_type,boundary_value):
        self.type=boundary_type
        self.value=boundary_value
        
class Space:
    def __init__(self):
        pass
    
    def créer_grille(self,lignes,colonnes):
        self.lignes=lignes      #grille de domaine
        self.colonnes=colonnes
        self.u=np.zeros((self.lignes+2,self.colonnes+2))    #matrice des vitesses
        self.v=np.zeros((self.lignes+2,self.colonnes+2))
        self.u_star=np.zeros((self.lignes+2,self.colonnes+2))
        self.v_star=np.zeros((self.lignes+2,self.colonnes+2))
        self.u_next=np.zeros((self.lignes+2,self.colonnes+2))
        self.v_next=np.zeros((self.lignes+2,self.colonnes+2))
        self.u_c=np.zeros((self.lignes,self.colonnes))
        self.v_c=np.zeros((self.lignes,self.colonnes))
        self.p=np.zeros((self.lignes+2,self.colonnes+2))    #matrice de pression
        self.p_c=np.zeros((self.lignes,self.colonnes))
        self.choix_terme_source()     #choix du terme de source par défaut    
        
    def choix_deltas(self,largeur,longueur):
        self.dx=longueur/(self.colonnes-1)
        self.dy=largeur/(self.lignes-1)
    def choix_u_initiale(self,U):
        self.u=U*self.u
        
    def choix_v_initiale(self,V):
        self.v=V*self.v
        
    def choix_pression_initiale(self,P):
        self.p=P*self.p
    def choix_terme_source(self,S_x=0,S_y=0):
        self.S_x=S_x
        self.S_y=S_y
        
class Fluid:
    def __init__(self,rho,mu):
        self.choix_propriétés_du_fluide(rho,mu)
    
    def choix_propriétés_du_fluide(self,rho,mu):
        self.rho=rho
        self.mu=mu
 

# !!!! changer les conditions initiales selon le modèle (à effacer quand c'est fait)

#choix des conditions initiales de la vélocités selon x aux frontières 
        
def choix_limites_u(space,left,right,top,bottom):
    if(left.type=="D"):     #il y a 2 types de conditions aux limites, les limites de Dirichlet et de Neumann
        space.u[:,0]=left.value     #les limites de Neumann sont des valeurs de la dérivée temporelle de la variable (dépendante de la frontière) à la frontière
    elif(left.type=="N"):   #les limites de Dirichlet sont des valeurs de la variable (dépendante de la frontière) à la frontière
        space.u[:,0]=-left.value*space.dx+space.u[:,1]
    if(right.type=="D"):
        space.u[:,-1]=right.value
    elif(right.type=="N"):
        space.u[:,-1]=right.value*space.dx+space.u[:,-2]
    if(top.type=="D"):
        space.u[-1,:]=2*top.value-space.u[-2,:]
    elif(top.type=="N"):
        space.u[-1,:]=-top.value*space.dy+space.u[-2,:]
    if(bottom.type=="D"):
        space.u[0,:]=2*bottom.value-space.u[1,:]
    elif(bottom.type=="N"):
        space.u[0,:]=bottom.value*space.dy+space.u[1,:]

#choix des conditions initiales de la vélocités selon y aux frontières
        
def choix_limites_v(space,left,right,top,bottom):
    if(left.type=="D"):
        space.v[:,0]=2*left.value-space.v[:,1]
    elif(left.type=="N"):
        space.v[:,0]=-left.value*space.dx+space.v[:,1]
    if(right.type=="D"):
        space.v[:,-1]=2*right.value-space.v[:,-2]
    elif(right.type=="N"):
        space.v[:,-1]=right.value*space.dx+space.v[:,-2]
    if(top.type=="D"):
        space.v[-1,:]=top.value
    elif(top.type=="N"):
        space.v[-1,:]=-top.value*space.dy+space.v[-2,:]
    if(bottom.type=="D"):
        space.v[0,:]=bottom.value
    elif(bottom.type=="N"):
        space.v[0,:]=bottom.value*space.dy+space.v[1,:]
        
#choix des conditions initiales de la pression
        
def choix_limites_pression(space,left,right,top,bottom):
    if(left.type=="D"):
        space.p[:,0]=left.value
    elif(left.type=="N"):
        space.p[:,0]=-left.value*space.dx+space.p[:,1]
    if(right.type=="D"):
        space.p[1,-1]=right.value
    elif(right.type=="N"):
        space.p[:,-1]=right.value*space.dx+space.p[:,-2]
    if(top.type=="D"):
        space.p[-1,:]=top.value
    elif(top.type=="N"):
        space.p[-1,:]=-top.value*space.dy+space.p[-2,:]
    if(bottom.type=="D"):
        space.p[0,:]=bottom.value
    elif(bottom.type=="N"):
        space.p[0,:]=bottom.value*space.dy+space.p[1,:]
        
def choix_pas_temporel(CFL,space,fluid):
    with np.errstate(divide='ignore'):
        dt=CFL/np.sum([np.amax(space.u)/space.dx,np.amax(space.v)/space.dy])
    if np.isinf(dt):    #évite l'erreur si dt est infini du à une vitesse initiale nulle
        dt=CFL*(space.dx+space.dy)
    space.dt=dt
    
    
def obtenir_vitesses_etoile(space,fluid):
    rows=int(space.lignes)  #on enregistre les paramètres des objets comme des variables locales
    cols=int(space.colonnes)
    u=space.u.astype(float,copy=False)
    v=space.v.astype(float,copy=False)
    dx=float(space.dx)
    dy=float(space.dy)
    dt=float(space.dt)
    S_x=float(space.S_x)
    S_y=float(space.S_y)
    rho=float(fluid.rho)
    mu=float(fluid.mu)
    u_star=u.copy()
    v_star=v.copy()
    u1_y=(u[2:,1:cols+1]-u[0:rows,1:cols+1])/(2*dy) #calcul des dérivées de u et v avec la méthode des différences finies
    u1_x=(u[1:rows+1,2:]-u[1:rows+1,0:cols])/(2*dx)
    u2_y=(u[2:,1:cols+1]-2*u[1:rows+1,1:cols+1]+u[0:rows,1:cols+1])/(dy**2)
    u2_x=(u[1:rows+1,2:]-2*u[1:rows+1,1:cols+1]+u[1:rows+1,0:cols])/(dx**2)
    v_face=(v[1:rows+1,1:cols+1]+v[1:rows+1,0:cols]+v[2:,1:cols+1]+v[2:,0:cols])/4
    u_star[1:rows+1,1:cols+1]=u[1:rows+1,1:cols+1]-dt*(u[1:rows+1,1:cols+1]*u1_x+v_face*u1_y)+(dt*(mu/rho)*(u2_x+u2_y))+(dt*S_x)   
    v1_y=(v[2:,1:cols+1]-v[0:rows,1:cols+1])/(2*dy)
    v1_x=(v[1:rows+1,2:]-v[1:rows+1,0:cols])/(2*dx)
    v2_y=(v[2:,1:cols+1]-2*v[1:rows+1,1:cols+1]+v[0:rows,1:cols+1])/(dy**2)
    v2_x=(v[1:rows+1,2:]-2*v[1:rows+1,1:cols+1]+v[1:rows+1,0:cols])/(dx**2)
    u_face=(u[1:rows+1,1:cols+1]+u[1:rows+1,2:]+u[0:rows,1:cols+1]+u[0:rows,2:])/4
    v_star[1:rows+1,1:cols+1]=v[1:rows+1,1:cols+1]-dt*(u_face*v1_x+v[1:rows+1,1:cols+1]*v1_y)+(dt*(mu/rho)*(v2_x+v2_y))+(dt*S_y)
    space.u_star=u_star.copy()  #on enregistre les vitesses etoile sur l'objet space
    space.v_star=v_star.copy()    
    
    
#on résoud itérativement l'équation de Poisson de la pression avec les vitesses etoiles pour calculer la pression à t+dt
def résoudre_Poisson(space,fluid,left,right,top,bottom):
    rows=int(space.lignes)
    cols=int(space.colonnes)
    u_star=space.u_star.astype(float,copy=False)
    v_star=space.v_star.astype(float,copy=False)
    p=space.p.astype(float,copy=False)
    dx=float(space.dx)
    dy=float(space.dy)
    dt=float(space.dt)
    rho=float(fluid.rho)
    factor=1/(2/dx**2+2/dy**2)
    error=1     #on définit l'erreur et la tolérence pour la convergence
    tol=1e-3
    ustar1_x=(u_star[1:rows+1,2:]-u_star[1:rows+1,0:cols])/(2*dx)   #on évalue la dérivée des vitesses etoiles
    vstar1_y=(v_star[2:,1:cols+1]-v_star[0:rows,1:cols+1])/(2*dy)
    i=0
    while(error>tol):
        i+=1
        p_old=p.astype(float,copy=True) #on enregistre la pression actuelle sous le nom de p_old
        p2_xy=(p_old[2:,1:cols+1]+p_old[0:rows,1:cols+1])/dy**2+(p_old[1:rows+1,2:]+p_old[1:rows+1,0:cols])/dx**2   #on évalue la dérivée seconde de la pression avec p_old
        p[1:rows+1,1:cols+1]=(p2_xy)*factor-(rho*factor/dt)*(ustar1_x+vstar1_y) #calcul de la nouvelle pression
        error=np.amax(abs(p-p_old)) #on trouve l'erreur maximale entre l'ancienne et la nouvelle pression
        choix_limites_pression(space,left,right,top,bottom)   #on applique les conditions limites de la pression
        if(i>500):  #condition pour sortir du while si on obtient rien après 500 itération
            tol=10*tol
            
    
#on calcule la vitesse a t+dt en utilisant la pression a t+dt et la vitesse etoile
def résolution_équation_quantité_mouvement(space,fluid):
    rows=int(space.lignes)
    cols=int(space.colonnes)
    u_star=space.u_star.astype(float,copy=False)
    v_star=space.v_star.astype(float,copy=False)
    p=space.p.astype(float,copy=False)
    dx=float(space.dx)
    dy=float(space.dy)
    dt=float(space.dt)
    rho=float(fluid.rho)
    u=space.u.astype(float,copy=False)
    v=space.v.astype(float,copy=False)
    p1_x=(p[1:rows+1,2:]-p[1:rows+1,0:cols])/(2*dx) #on évalue la dérivée par rapport à x de la pression
    u[1:rows+1,1:cols+1]=u_star[1:rows+1,1:cols+1]-(dt/rho)*p1_x    #on calcule u avec un pas de dt
    p1_y=(p[2:,1:cols+1]-p[0:rows,1:cols+1])/(2*dy) #on évalue la dérivée par arpport à y de la pression
    v[1:rows+1,1:cols+1]=v_star[1:rows+1,1:cols+1]-(dt/rho)*p1_y    #on calcule v avec un pas de dt

def centre_PUV(space):
    space.p_c=space.p[1:-1,1:-1]
    space.u_c=space.u[1:-1,1:-1]
    space.v_c=space.v[1:-1,1:-1]
    

def créer_répertoire_des_résultats(wipe=False):
    cwdir=os.getcwd()   #obtient le chemin du répertoire des résultats
    dir_path=os.path.join(cwdir,"Result")
    if not os.path.isdir(dir_path):     #si le répertoire n'existe pas, le créer
        os.makedirs(dir_path,exist_ok=True)
    else:
        if wipe:    #si wipe est True i.e. on souhaite effacer le répertoire, l'efface 
            os.chdir(dir_path)
            filelist=os.listdir()
            for file in filelist:
                os.remove(file)
    os.chdir(cwdir)
            
    
def écrit_les_données(space,iteration,interval):
    if(iteration%interval==0):
        dir_path=os.path.join(os.getcwd(),"Result")
        filename="PUV{0}.txt".format(iteration)
        path=os.path.join(dir_path,filename)
        with open(path,"w") as f:
            for i in range(space.lignes):
                for j in range(space.colonnes):
                    f.write("{}\t{}\t{}\n".format(space.p_c[i,j],space.u_c[i,j],space.v_c[i,j]))

#Valeurs spatiales et temporelles
longueur=12     #longueur du domaine choisi selon x
largeur=4       #largeur du domaine choisi selon y
colonnes=64    #nombre de point selon x dans la grille 
lignes=64     #nombre de point selon y dans la grille 
cavity=Space()  #créer un objet de classe Space
cavity.créer_grille(lignes,colonnes)
cavity.choix_deltas(largeur,longueur)


#Propriétés du fluide
rho=1               #densité du fluide
mu=1.8*10**(-2)     #viscosité du fluide
water=Fluid(rho,mu) #créer un objet de classe Fluid

#Propriétés des limites
gauche=Boundary("D",1)
droite=Boundary("N",0)
haut=Boundary("D",0)
bas=Boundary("D",0)
pgauche=Boundary("N",0)
pressureatm=Boundary("D",0)

#Compteurs initiaux
t=0
i=0

#Paramètres de la mesure
u_mes=[]
v_mes=[]
p_mes=[]
t_mes=[]

#Paramètres de simulation
time=100        # "temps" du décompte de la simulation
CFL_number=0.8  #réduire ceci si la solution diverge,mais pas trop. Valeur de base : 0.8
file_flag=1     #1 pour enregistrer les résultats, 0 pour ne pas les garder
interval=100    #enregistre les valeurs par 100 nombre d'itération 

#Lancer la simulation
print("######## Beginning FlowPy Simulation ########")
print("#############################################")
print("# Simulation time: {0:.2f}".format(time))
print("# Mesh: {0} x {1}".format(colonnes,lignes))
print("# Re/u: {0:.2f}\tRe/v:{1:.2f}".format(rho*longueur/mu,rho*largeur/mu))
print("# Save outputs to text file: {0}".format(bool(file_flag)))
créer_répertoire_des_résultats(wipe=True)

while(t<time):
    sys.stdout.write("\rSimulation time left: {0:.2f}".format(time-t))  #affiche le temps restant
    sys.stdout.flush()
    choix_pas_temporel(CFL_number,cavity,water)
    timestep=cavity.dt
    choix_limites_u(cavity,gauche,droite,haut,bas)
    choix_limites_v(cavity,gauche,droite,haut,bas)
    choix_limites_pression(cavity,pgauche,pgauche,pressureatm,pgauche)
    obtenir_vitesses_etoile(cavity,water)
    résoudre_Poisson(cavity,water,pgauche,pgauche,pressureatm,pgauche)
    résolution_équation_quantité_mouvement(cavity,water)
    centre_PUV(cavity)  #enregistre les variables et écrit les données
    if(file_flag==1):
        écrit_les_données(cavity,i,interval)
    u_mes.append(cavity.u_c[10,32])
    v_mes.append(cavity.v_c[10,32])
    p_mes.append(cavity.p_c[10,32])
    t+=timestep     #augmente le temps d'un pas
    i+=1            #incrémente le compteur


cwdir=os.getcwd()   #obtiens le chemin du repertoire des résultats
dir_path=os.path.join(cwdir,"Result")
os.chdir(dir_path)
filenames=[]    #enregistre le nom des données obtenues
iterations=[]
for root,dirs,files in os.walk(dir_path):
    for datafile in files:
        if "PUV" in datafile:
            filenames.append(datafile)
            no_ext_file=datafile.replace(".txt","").strip()
            iter_no=int(no_ext_file.split("V")[-1])
            iterations.append(iter_no)
initial_iter=np.amin(iterations)    #discerne l'itération finale des intervalles
final_iter=np.amax(iterations)
inter=(final_iter - initial_iter)/(len(iterations)-1)
number_of_frames=len(iterations)
sorted_iterations=np.sort(iterations)


def read_datafile(iteration):
    filename="PUV{0}.txt".format(iteration) #choisi le nom des données et le lien selon l'itération
    filepath=os.path.join(dir_path,filename)
    arr=np.loadtxt(filepath,delimiter="\t") #charge les données texte comme des array numpy
    rows,cols=arr.shape
    p_p=np.zeros((lignes,colonnes)) #défini des arrays vides pour la pression et les vitesses
    u_p=np.zeros((lignes,colonnes))
    v_p=np.zeros((lignes,colonnes))
    p_arr=arr[:,0]  #organise les array importées en variables
    u_arr=arr[:,1]
    v_arr=arr[:,2]
    p_p=p_arr.reshape((lignes,colonnes))    #transforme les données 1d en 2d
    u_p=u_arr.reshape((lignes,colonnes))
    v_p=v_arr.reshape((lignes,colonnes))
    return p_p,u_p,v_p


x=np.linspace(0,longueur,colonnes)  #créer une grille pour les valeurs de X et Y pour la figure
y=np.linspace(0,largeur,lignes)
[X,Y]=np.meshgrid(x,y)


index_cut_x=int(colonnes/10)    #determine l'indexation pour le tracé de flux (seulemnt 10 points)
index_cut_y=int(lignes/10)


fig=plt.figure(figsize=(16,8))
ax=plt.axes(xlim=(0,longueur),ylim=(0,largeur))
p_p,u_p,v_p=read_datafile(0)    #créer le contour initial et le tracé de flux et la barre de couleur
ax.set_xlim([0,longueur])
ax.set_ylim([0,largeur])
ax.set_xlabel("$x$",fontsize=12)
ax.set_ylabel("$y$",fontsize=12)
ax.set_title("Frame No: 0")
cont=ax.contourf(X,Y,p_p)
stream=ax.streamplot(X[::index_cut_y,::index_cut_x],Y[::index_cut_y,::index_cut_x],u_p[::index_cut_y,::index_cut_x],v_p[::index_cut_y,::index_cut_x],color="k")
fig.colorbar(cont,label='Pression')
fig.tight_layout()


def animate(i):
    sys.stdout.write("\rFrames remaining: {0:03d}".format(len(sorted_iterations)-i))    #affiche le nombre d'images restantes à être afficher durant l'animation
    sys.stdout.flush()
    iteration=sorted_iterations[i]  #obtient des itération de sorted_itertion de manière séquentiel 
    p_p,u_p,v_p=read_datafile(iteration)    #utilise la fonction read_datafile pour obtenir la pression et les vitesses
    ax.clear()  #efface les plots précédents pour afficher le nouveau
    ax.set_xlim([0,longueur])
    ax.set_ylim([0,largeur])
    ax.set_xlabel("$x$",fontsize=12)
    ax.set_ylabel("$y$",fontsize=12)
    ax.set_title("Frame No: {0}".format(i))
    cont=ax.contourf(X,Y,p_p)
    stream=ax.streamplot(X[::index_cut_y,::index_cut_x],Y[::index_cut_y,::index_cut_x],u_p[::index_cut_y,::index_cut_x],v_p[::index_cut_y,::index_cut_x],color="k")
    return cont,stream


print("######## Making FlowPy Animation ########")
print("#########################################")
#anim=animation.FuncAnimation(fig,animate,frames=number_of_frames,interval=50,blit=False)
t_mes=[i*float(cavity.dt) for i in range(0,len(u_mes)-1)]
cursor = conn.cursor()
for k in range(0,len(u_mes)-1):
    cursor.execute("""INSERT INTO valeurs(condition ,temps ,horizontal ,vertical ,pressure) VALUES (?,?,?,?,?)""",(mu,t_mes[k],u_mes[k],v_mes[k],p_mes[k]))

fig=cursor.execute("""SELECT condition,temps,horizontal,vertical,pressure from valeurs""")
img0=[]
img1=[]
img2=[]
img3=[]
img4=[]
for row in fig:
    img0.append(row[0])
    img1.append(row[1])
    img2.append(row[2])
    img3.append(row[3])
    img4.append(row[4])


plt.close('all')
img1=[img1[i]*20 for i in range(len(img1)-1)]
vitesse=[(img2[i]**2+img3[i]**2)**(1/2) for i in range(len(img1))]
pression=[img4[i]/8+1 for i in range(len(img1))]

matplotlib.rc('xtick', labelsize=20) 
matplotlib.rc('ytick', labelsize=20)
""
N=10000
Q=2.098
w0=1.192
a=((w0/Q)**2+4*(w0)**4)**(1/2)
tau=14
A=1.05
n=100
t=[i*n/N for i in range(N)]
f=lambda x:A*np.exp(-2*x/tau)*(np.sin(2*x/a))+1
y=[f(i) for i in t]
plt.plot(t,y,color='r')
plt.plot(img1,pression,color='blue')
plt.legend(['Modèle correspondant','Vitesse simulée'],fontsize=30)
plt.title("Vitesse en (2,2)",fontsize=40)
plt.xlabel("$Temps$",fontsize=30)
plt.ylabel(r"$\frac{V}{V_{0}} $",fontsize=30,rotation=0,labelpad=20)
plt.xlim(0,100)
plt.ylim(0.6,1.8)
plt.grid()
plt.show()
""
"""
plt.plot(img1,vitesse,color='blue')
plt.title("Pression en (2,2)",fontsize=40)
plt.xlabel("$Temps$",fontsize=30)
plt.ylabel(r"$\frac{P}{P_{0}} $",fontsize=30,rotation=0,labelpad=20)
plt.xlim(0,100)
plt.ylim(0.035,0.2)
plt.grid()
plt.show()
"""
conn.commit()

















