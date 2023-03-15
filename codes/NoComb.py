from turtle import color
import numpy as np
import matplotlib.pyplot as plt
pi = np.pi

R=287
Ru=8.314
LHV=42.588e6
P_atm=101325
rho_atm=1.2923
gamma=1.4
cp=1.005e3
cv=cp/gamma
T_atm=300
S=80*10**-3 #Stroke
Rad=S/2 #radius of crank
L=S*1.75 #length of connecting rod
r=18.5 #compression ratio
Bore=76.5*10**-3 # bore of cylinder
LIFT_MAX=4.5*10**-3 # Maximum Lift
Vd=(pi*(Bore**2)*Rad)/2
Vcl=Vd/(r-1)
#S=(Vd*4)/((Bore**2)*pi)

N=5000 #engine rpm

W= (2*pi*N)/60
CM=(N*2*Rad)/30 # Mean piston Speed
IVO=pi/12 # IVO Before TDC
IVC=pi/9 # IVC After BDC
EVO=pi/9 # EVO before BDC
EVC=IVO # EVC after TDC
D=30*10**-3 # Inlet Valve Diameter
D2=2*D # Inlet manifold diameter
A1=pi*0.25*D**2
A2=pi*0.25*D2**2
FIS=pi/36 #fuel injection starts before TDC
FIC=pi/7.2 #fuel injector closes after TDC
MWair=0.02987
dtheta=0.0002
a=0.137
b=0.0387e-3
stoichiometric_af_ratio=14.5
#--------------------------------------------------------------------------#
#--------------Modelling of Intake-----------------------------------------#
#--------------------------------------------------------------------------#
theta_i = np.arange(-IVO,pi + IVC + dtheta, dtheta)
N_i = len(theta_i)
P_i = np.zeros(N_i)
P_i=np.zeros(N_i)
dP_i=np.zeros(N_i)
T_i=np.zeros(N_i)
V_i=Vd*((r/(r-1)) - ((1-np.cos(theta_i - pi))/2) + (L/S) - 0.5*np.sqrt((2*L/S)**2 - (np.sin(theta_i - pi))**2))
dV_i=np.gradient(V_i,dtheta)
m_i=np.zeros(N_i)
avg_cyl_gas_velocity=2.28*CM
dm_i=np.zeros(N_i)
Twall=400
T_im=300
P_m=20e3
P_im=101325
A_i= (pi*Bore**2)/2 + (pi*Bore*S/2)*(1 - np.cos(theta_i) + (L/Rad) - np.sqrt((L/Rad)**2 - (np.sin(theta_i))**2))
P_i[0]=P_m
T_i[0]=T_im
m_i[0]=(P_i[0]*V_i[0])/(T_i[0]*(Ru/MWair))
for i in range(0,N_i-1):
    cv = 0.2*T_i[i] + 700
    cp = R+cv
    gamma=cp/cv
    CD=(107.78*((LIFT_MAX/D)**4))-(77.204*((LIFT_MAX/D)**3))+(14.1*((LIFT_MAX/D)**2))-(1.01*(LIFT_MAX/D))+0.6687
    CURTAREA=(pi/4)*D**2
    Cht=3.26*(10**-3)*(P_i[i]**(0.8))*(T_i[i]**(-0.55))*(Bore**(-0.2))*(avg_cyl_gas_velocity**(0.8))
    if P_i[i]<=P_im:
        if (P_i[i]/P_im)<=(2/(gamma+1))**(gamma/(gamma-1)):
            dm_i[i]=(P_im/(np.sqrt(R*T_im)))*((np.sqrt(gamma))*((2*gamma/(gamma+1))**((gamma+1)/2*(gamma-1))))*CD*CURTAREA/W
            m_i[i+1]=(dtheta*dm_i[i])+m_i[i]
        else:
            dm_i[i]=(P_im/(np.sqrt(R*T_im)))*(((P_i[i]/P_im)**(1/gamma))*np.sqrt((2*gamma/(gamma-1))*(1 - (P_i[i]/P_im)**((gamma-1)/gamma))))*CD*CURTAREA/W
            m_i[i+1]=(dtheta*dm_i[i])+m_i[i]
        n=m_i[i]/MWair
        dP_i[i] = ((W*dm_i[i]*cp*T_i[i])-A_i[i]*Cht*(T_i[i]-Twall) - W*dV_i[i]*(P_i[i]*(1 + cv/(R)) + 2*a*b*n*n*n*cv/(R*((V_i[i])**3)) - a*n*n*cv/(R*((V_i[i])**2))))/(W*(V_i[i]*cv/(R) - cv*b*n/R))
        P_i[i+1]=(dtheta*dP_i[i])+P_i[i]
        T_i[i+1]=(P_i[i+1]*V_i[i+1])/(m_i[i+1]*(Ru/MWair))
    else:
        if (P_im/P_i[i])<=(2/(gamma+1))**(gamma/(gamma-1)):
            dm_i[i]=-(P_i[i]/np.sqrt(R*T_i[i]))*((np.sqrt(gamma))*((2*gamma/(gamma+1))**((gamma+1)/2*(gamma-1))))*CD*CURTAREA/W
            m_i[i+1]=(dtheta*dm_i[i])+m_i[i]
        else:
            dm_i[i]=(P_i[i]/np.sqrt(R*T_i[i]))*(((P_im/P_i[i])**(1/gamma))*np.sqrt((2*gamma/(gamma-1))*(1 - (P_im/P_i[i])**((gamma-1)/gamma))))*CD*CURTAREA/W
            m_i[i+1]=(dtheta*dm_i[i])+m_i[i]
        n=m_i[i]/MWair
        dP_i[i] = (W*dm_i[i]*cp*T_i[i] - W*dV_i[i]*(P_i[i]*(1 + cv/(R)) + 2*a*b*n*n*n*cv/(R*((V_i[i])**3)) - a*n*n*cv/(R*((V_i[i])**2))))/(W*(V_i[i]*cv/(R) - cv*b*n/R))
        P_i[i+1]=(dtheta*dP_i[i])+P_i[i]
        T_i[i+1]=(P_i[i+1]*V_i[i+1])/(m_i[i+1]*(Ru/MWair))
mass_of_air_inducted = m_i[N_i-1]
#--------------------------------------------------------------------------#
#----------------------Modelling of Compression----------------------------#
#--------------------------------------------------------------------------#
theta_c = np.arange(pi + IVC,2*pi - FIS + dtheta, dtheta)
N_c = len(theta_c)
P_c = np.zeros(N_c)
P_c=np.zeros(N_c)
dP_c=np.zeros(N_c)
T_c=np.zeros(N_c)
V_c=Vd*((r/(r-1)) - ((1-np.cos(theta_c - pi))/2) + (L/S) - 0.5*np.sqrt((2*L/S)**2 - (np.sin(theta_c - pi))**2))
dV_c=np.gradient(V_c,dtheta)
A_c= (pi*Bore**2)/2 + (pi*Bore*S/2)*(1 - np.cos(theta_c) + (L/Rad) - np.sqrt((L/Rad)**2 - (np.sin(theta_c))**2))
Twall=400
avg_cyl_gas_velocity=2.28*CM
n = m_i[N_i-1]/MWair
T_c[0]=T_i[N_i-1]
P_c[0]=P_i[N_i-1]
for i in range(0,N_c-1):
    cv = 0.2*T_c[i] + 700
    Cht=3.26*(10**-3)*(P_c[i]**(0.8))*(T_c[i]**(-0.55))*(Bore**(-0.2))*(avg_cyl_gas_velocity**(0.8))
    dP_c[i] = -(A_c[i]*Cht*(T_c[i]-Twall) + W*dV_c[i]*(P_c[i]*(1 + cv/(R)) + 2*a*b*n*n*n*cv/(R*((V_c[i])**3)) - a*n*n*cv/(R*((V_c[i])**2))))/(W*(V_c[i]*cv/(R) - cv*b*n/R))
    P_c[i+1]=(dtheta*dP_c[i])+P_c[i]
    T_c[i+1]=(P_c[i+1]*V_c[i+1])/(Ru*n)
#--------------------------------------------------------------------------#
#---------------------Modelling of Combustion------------------------------#
#--------------------------------------------------------------------------#
theta_b = np.arange(2*pi - FIS,2*pi + FIC +0.15 + dtheta,dtheta)
N_b = len(theta_b)
P_b = np.zeros(N_b)
P_b=np.zeros(N_b)
dP_b=np.zeros(N_b)
T_b=np.zeros(N_b)
V_b=Vd*((r/(r-1)) - ((1-np.cos(theta_b - pi))/2) + (L/S) - 0.5*np.sqrt((2*L/S)**2 - (np.sin(theta_b - pi))**2))
dV_b=np.gradient(V_b,dtheta)
A_c= (pi*Bore**2)/2 + (pi*Bore*S/2)*(1 - np.cos(theta_b) + (L/Rad) - np.sqrt((L/Rad)**2 - (np.sin(theta_b))**2))
Twall=700
avg_cyl_gas_velocity=2.28*CM
n = m_i[N_i-1]/MWair
T_b[0]=T_c[N_c-1]
P_b[0]=P_c[N_c-1]
for i in range(0,N_b-1):
    cv = 0.2*T_b[i] + 900
    Cht=3.26*(10**-3)*(P_b[i]**(0.8))*(T_b[i]**(-0.55))*(Bore**(-0.2))*(avg_cyl_gas_velocity**(0.8))
    dP_b[i] = -(A_c[i]*Cht*(T_b[i]-Twall) + W*dV_b[i]*(P_b[i]*(1 + cv/(R)) + 2*a*b*n*n*n*cv/(R*((V_b[i])**3)) - a*n*n*cv/(R*((V_b[i])**2))))/(W*(V_b[i]*cv/(R) - cv*b*n/R))
    P_b[i+1]=(dtheta*dP_b[i])+P_b[i]
    T_b[i+1]=(P_b[i+1]*V_b[i+1])/(Ru*n)
#--------------------------------------------------------------------------#
#---------------------Modelling of Expansion-------------------------------#
#--------------------------------------------------------------------------#
theta_e=np.arange(2*pi + FIC+0.15,dtheta+ 3*pi - EVO,dtheta)
N_e = len(theta_e)
P_e=np.zeros(N_e)
dP_e=np.zeros(N_e)
T_e=np.zeros(N_e)
V_e=Vd*((r/(r-1)) - ((1-np.cos(theta_e - pi))/2) + (L/S) - 0.5*np.sqrt((2*L/S)**2 - (np.sin(theta_e - pi))**2))
dV_e=np.gradient(V_e,dtheta)
A_e= (pi*Bore**2)/2 + (pi*Bore*S/2)*(1 - np.cos(theta_e) + (L/Rad) - np.sqrt((L/Rad)**2 - (np.sin(theta_e))**2))
Twall=1100
avg_cyl_gas_velocity=2.28*CM
n = m_i[N_i-1]/MWair
T_e[0]=T_b[N_b-1]
P_e[0]=P_b[N_b-1]
for i in range(0,N_e-1):
    cv = 0.2*T_e[i] + 700
    Cht=3.26*(10**-3)*(P_e[i]**(0.8))*(T_e[i]**(-0.55))*(Bore**(-0.2))*(avg_cyl_gas_velocity**(0.8))
    dP_e[i] = -(A_e[i]*Cht*(T_e[i]-Twall) + W*dV_e[i]*(P_e[i]*(1 + cv/(R)) + 2*a*b*n*n*n*cv/(R*((V_e[i])**3)) - a*n*n*cv/(R*((V_e[i])**2))))/(W*(V_e[i]*cv/(R) - cv*b*n/R))
    P_e[i+1]=(dtheta*dP_e[i])+P_e[i]
    T_e[i+1]=(P_e[i+1]*V_e[i+1])/(Ru*n)
#--------------------------------------------------------------------------#
#---------------------Modelling of Exhuast---------------------------------#
#--------------------------------------------------------------------------#
theta_ex=np.arange(3*pi - EVO,dtheta+ 4*pi + EVC,dtheta)
N_ex = len(theta_ex)
P_ex=np.zeros(N_ex)
dP_ex=np.zeros(N_ex)
T_ex=np.zeros(N_ex)
V_ex=Vd*((r/(r-1)) - ((1-np.cos(theta_ex - pi))/2) + (L/S) - 0.5*np.sqrt((2*L/S)**2 - (np.sin(theta_ex - pi))**2))
dV_ex=np.gradient(V_ex,dtheta)
m_ex=np.zeros(N_ex)
avg_cyl_gas_velocity=6.18*CM
dm_ex=np.zeros(N_ex)
Twall=700
T_exm=1000
P_exm=101325
A_ex= (pi*Bore**2)/2 + (pi*Bore*S/2)*(1 - np.cos(theta_ex) + (L/Rad) - np.sqrt((L/Rad)**2 - (np.sin(theta_ex))**2))
P_ex[0]=P_e[N_e-1]
T_ex[0]=T_e[N_e-1]
m_ex[0]=m_i[N_i-1]

for i in range(0,N_e-1):
    cv = 0.2*T_ex[i] + 700
    cp = R+cv
    gamma=cp/cv
    CD=(107.78*((LIFT_MAX/D)**4))-(77.204*((LIFT_MAX/D)**3))+(14.1*((LIFT_MAX/D)**2))-(1.01*(LIFT_MAX/D))+0.6687
    CURTAREA=(pi/4)*D**2
    Cht=3.26*(10**-3)*(P_ex[i]**(0.8))*(T_ex[i]**(-0.55))*(Bore**(-0.2))*(avg_cyl_gas_velocity**(0.8))
    if P_exm<=P_ex[i]:
        if (P_exm/P_ex[i])<=(2/(gamma+1))**(gamma/(gamma-1)):
            dm_ex[i]=-(P_ex[i]/(np.sqrt(R*T_ex[i])))*((np.sqrt(gamma))*((2*gamma/(gamma+1))**((gamma+1)/2*(gamma-1))))*CD*CURTAREA/W
            m_ex[i+1]=(dtheta*dm_ex[i])+m_ex[i]
        else:
            dm_ex[i]=(P_ex[i]/(np.sqrt(R*T_ex[i])))*(((P_exm/P_ex[i])**(1/gamma))*np.sqrt((2*gamma/(gamma-1))*(1 - (P_exm/P_ex[i])**((gamma-1)/gamma))))*CD*CURTAREA/W
            m_ex[i+1]=(dtheta*dm_ex[i])+m_ex[i]
        n=m_ex[i]/MWair
        dP_ex[i] = ((W*dm_ex[i]*cp*T_ex[i])-A_ex[i]*Cht*(T_ex[i]-Twall)-W*dV_ex[i]*(P_ex[i]*(1 + cv/(R)) + 2*a*b*n*n*n*cv/(R*((V_ex[i])**3)) - a*n*n*cv/(R*((V_ex[i])**2))))/(W*(V_ex[i]*cv/(R) - cv*b*n/R))
        P_ex[i+1]=(dtheta*dP_ex[i])+P_ex[i]
        T_ex[i+1]=(P_ex[i+1]*V_ex[i+1])/(m_ex[i+1]*(Ru/MWair))
    else:
        if (P_ex[i]/P_exm)<=(2/(gamma+1))**(gamma/(gamma-1)):
            dm_ex[i]=(P_exm/np.sqrt(R*T_exm))*((np.sqrt(gamma))*((2*gamma/(gamma+1))**((gamma+1)/2*(gamma-1))))*CD*CURTAREA/W
            m_ex[i+1]=(dtheta*dm_ex[i])+m_ex[i]
        else:
            dm_ex[i]=(P_exm/np.sqrt(R*T_exm))*(((P_ex[i]/P_exm)**(1/gamma))*np.sqrt((2*gamma/(gamma-1))*(1 - (P_ex[i]/P_exm)**((gamma-1)/gamma))))*CD*CURTAREA/W
            m_ex[i+1]=(dtheta*dm_ex[i])+m_ex[i]
        n=m_ex[i]/MWair
        dP_ex[i] = (W*dm_ex[i]*cp*T_ex[i]-W*dV_ex[i]*(P_ex[i]*(1 + cv/(R)) + 2*a*b*n*n*n*cv/(R*((V_ex[i])**3)) - a*n*n*cv/(R*((V_ex[i])**2))))/(W*(V_ex[i]*cv/(R) - cv*b*n/R))
        P_ex[i+1]=(dtheta*dP_ex[i])+P_ex[i]
        T_ex[i+1]=(P_ex[i+1]*V_ex[i+1])/(m_ex[i+1]*(Ru/MWair))


theta = np.concatenate((theta_i,theta_c,theta_b,theta_e,theta_ex))
theta = theta*180/pi
theta = theta-360
P = np.concatenate((P_i,P_c,P_b,P_e,P_ex))
T = np.concatenate((T_i,T_c,T_b,T_e,T_ex))
V = np.concatenate((V_i,V_c,V_b,V_e,V_ex))
plt.plot(V,P,color='r',label='noComb')

#############################################################################3

R=287
Ru=8.314
LHV=42.588e6
P_atm=101325
rho_atm=1.2923
gamma=1.4
cp=1.005e3
cv=cp/gamma
T_atm=300
S=80*10**-3 #Stroke
Rad=S/2 #radius of crank
L=S*1.75 #length of connecting rod
r=18.5 #compression ratio
Bore=76.5*10**-3 # bore of cylinder
LIFT_MAX=4.5*10**-3 # Maximum Lift
Vd=(pi*(Bore**2)*Rad)/2
Vcl=Vd/(r-1)
#S=(Vd*4)/((Bore**2)*pi)

N=5000 #engine rpm

W= (2*pi*N)/60
CM=(N*2*Rad)/30 # Mean piston Speed
IVO=pi/12 # IVO Before TDC
IVC=pi/9 # IVC After BDC
EVO=pi/9 # EVO before BDC
EVC=IVO # EVC after TDC
D=30*10**-3 # Inlet Valve Diameter
D2=2*D # Inlet manifold diameter
A1=pi*0.25*D**2
A2=pi*0.25*D2**2
FIS=pi/36 #fuel injection starts before TDC
FIC=pi/7.2 #fuel injector closes after TDC
MWair=0.02987
dtheta=0.0002
a=0.137
b=0.0387e-3
stoichiometric_af_ratio=14.5
#--------------------------------------------------------------------------#
#--------------Modelling of Intake-----------------------------------------#
#--------------------------------------------------------------------------#
theta_i = np.arange(-IVO,pi + IVC + dtheta, dtheta)
N_i = len(theta_i)
P_i = np.zeros(N_i)
P_i=np.zeros(N_i)
dP_i=np.zeros(N_i)
T_i=np.zeros(N_i)
V_i=Vd*((r/(r-1)) - ((1-np.cos(theta_i - pi))/2) + (L/S) - 0.5*np.sqrt((2*L/S)**2 - (np.sin(theta_i - pi))**2))
dV_i=np.gradient(V_i,dtheta)
m_i=np.zeros(N_i)
avg_cyl_gas_velocity=2.28*CM
dm_i=np.zeros(N_i)
Twall=400
T_im=300
P_m=20e3
P_im=101325
A_i= (pi*Bore**2)/2 + (pi*Bore*S/2)*(1 - np.cos(theta_i) + (L/Rad) - np.sqrt((L/Rad)**2 - (np.sin(theta_i))**2))
P_i[0]=P_m
T_i[0]=T_im
m_i[0]=(P_i[0]*V_i[0])/(T_i[0]*(Ru/MWair))
for i in range(0,N_i-1):
    cv = 0.2*T_i[i] + 700
    cp = R+cv
    gamma=cp/cv
    CD=(107.78*((LIFT_MAX/D)**4))-(77.204*((LIFT_MAX/D)**3))+(14.1*((LIFT_MAX/D)**2))-(1.01*(LIFT_MAX/D))+0.6687
    CURTAREA=(pi/4)*D**2
    Cht=3.26*(10**-3)*(P_i[i]**(0.8))*(T_i[i]**(-0.55))*(Bore**(-0.2))*(avg_cyl_gas_velocity**(0.8))
    if P_i[i]<=P_im:
        if (P_i[i]/P_im)<=(2/(gamma+1))**(gamma/(gamma-1)):
            dm_i[i]=(P_im/(np.sqrt(R*T_im)))*((np.sqrt(gamma))*((2*gamma/(gamma+1))**((gamma+1)/2*(gamma-1))))*CD*CURTAREA/W
            m_i[i+1]=(dtheta*dm_i[i])+m_i[i]
        else:
            dm_i[i]=(P_im/(np.sqrt(R*T_im)))*(((P_i[i]/P_im)**(1/gamma))*np.sqrt((2*gamma/(gamma-1))*(1 - (P_i[i]/P_im)**((gamma-1)/gamma))))*CD*CURTAREA/W
            m_i[i+1]=(dtheta*dm_i[i])+m_i[i]
        n=m_i[i]/MWair
        dP_i[i] = ((W*dm_i[i]*cp*T_i[i])-A_i[i]*Cht*(T_i[i]-Twall) - W*dV_i[i]*(P_i[i]*(1 + cv/(R)) + 2*a*b*n*n*n*cv/(R*((V_i[i])**3)) - a*n*n*cv/(R*((V_i[i])**2))))/(W*(V_i[i]*cv/(R) - cv*b*n/R))
        P_i[i+1]=(dtheta*dP_i[i])+P_i[i]
        T_i[i+1]=(P_i[i+1]*V_i[i+1])/(m_i[i+1]*(Ru/MWair))
    else:
        if (P_im/P_i[i])<=(2/(gamma+1))**(gamma/(gamma-1)):
            dm_i[i]=-(P_i[i]/np.sqrt(R*T_i[i]))*((np.sqrt(gamma))*((2*gamma/(gamma+1))**((gamma+1)/2*(gamma-1))))*CD*CURTAREA/W
            m_i[i+1]=(dtheta*dm_i[i])+m_i[i]
        else:
            dm_i[i]=(P_i[i]/np.sqrt(R*T_i[i]))*(((P_im/P_i[i])**(1/gamma))*np.sqrt((2*gamma/(gamma-1))*(1 - (P_im/P_i[i])**((gamma-1)/gamma))))*CD*CURTAREA/W
            m_i[i+1]=(dtheta*dm_i[i])+m_i[i]
        n=m_i[i]/MWair
        dP_i[i] = (W*dm_i[i]*cp*T_i[i] - W*dV_i[i]*(P_i[i]*(1 + cv/(R)) + 2*a*b*n*n*n*cv/(R*((V_i[i])**3)) - a*n*n*cv/(R*((V_i[i])**2))))/(W*(V_i[i]*cv/(R) - cv*b*n/R))
        P_i[i+1]=(dtheta*dP_i[i])+P_i[i]
        T_i[i+1]=(P_i[i+1]*V_i[i+1])/(m_i[i+1]*(Ru/MWair))
mass_of_air_inducted = m_i[N_i-1]
#--------------------------------------------------------------------------#
#----------------------Modelling of Compression----------------------------#
#--------------------------------------------------------------------------#
theta_c = np.arange(pi + IVC,2*pi - FIS + dtheta, dtheta)
N_c = len(theta_c)
P_c = np.zeros(N_c)
P_c=np.zeros(N_c)
dP_c=np.zeros(N_c)
T_c=np.zeros(N_c)
V_c=Vd*((r/(r-1)) - ((1-np.cos(theta_c - pi))/2) + (L/S) - 0.5*np.sqrt((2*L/S)**2 - (np.sin(theta_c - pi))**2))
dV_c=np.gradient(V_c,dtheta)
A_c= (pi*Bore**2)/2 + (pi*Bore*S/2)*(1 - np.cos(theta_c) + (L/Rad) - np.sqrt((L/Rad)**2 - (np.sin(theta_c))**2))
Twall=400
avg_cyl_gas_velocity=2.28*CM
n = m_i[N_i-1]/MWair
T_c[0]=T_i[N_i-1]
P_c[0]=P_i[N_i-1]
for i in range(0,N_c-1):
    cv = 0.2*T_c[i] + 700
    Cht=3.26*(10**-3)*(P_c[i]**(0.8))*(T_c[i]**(-0.55))*(Bore**(-0.2))*(avg_cyl_gas_velocity**(0.8))
    dP_c[i] = -(A_c[i]*Cht*(T_c[i]-Twall) + W*dV_c[i]*(P_c[i]*(1 + cv/(R)) + 2*a*b*n*n*n*cv/(R*((V_c[i])**3)) - a*n*n*cv/(R*((V_c[i])**2))))/(W*(V_c[i]*cv/(R) - cv*b*n/R))
    P_c[i+1]=(dtheta*dP_c[i])+P_c[i]
    T_c[i+1]=(P_c[i+1]*V_c[i+1])/(Ru*n)
#--------------------------------------------------------------------------#
#---------------------Modelling of Combustion------------------------------#
#--------------------------------------------------------------------------#
theta_b=np.arange(2*pi - FIS,2*pi + FIC +0.15 + dtheta,dtheta)
N_b = len(theta_b)
P_b=np.zeros(N_b)
dP_b=np.zeros(N_b)
T_b=np.zeros(N_b)
dqc=np.zeros(N_b)
V_b=Vd*((r/(r-1)) - ((1-np.cos(theta_b - pi))/2) + (L/S) - 0.5*np.sqrt((2*L/S)**2 - (np.sin(theta_b - pi))**2))
dV_b=np.gradient(V_b,dtheta)
A_b= (pi*Bore**2)/2 + (pi*Bore*S/2)*(1 - np.cos(theta_b) + (L/Rad) - np.sqrt((L/Rad)**2 - (np.sin(theta_b))**2))
Twall=1200
avg_cyl_gas_velocity=2.28*CM
n = m_i[N_i-1]/MWair
lam=1.4 #equivalence air-fuel ratio
mfuel=m_i[N_i-1]/(stoichiometric_af_ratio*lam)
qav = LHV*mfuel
ab=6.908
comb_dur=FIC+FIS+0.15
mb=2.36
T_b[0]=T_c[N_c-1]
P_b[0]=P_c[N_c-1]
for i in range(0,N_b-1):
    cv = 0.2*T_b[i] + 700
    Cht=3.26*(10**-3)*(P_b[i]**(0.8))*(T_b[i]**(-0.55))*(Bore**(-0.2))*(avg_cyl_gas_velocity**(0.8))
    dqc[i] = ab*W*(mb+1)*((qav*pi/(180*comb_dur))**mb)*(((theta_b[i] - (2*pi-FIS))/comb_dur)**mb)*np.exp(-ab*(((theta_b[i] - (2*pi-FIS))/comb_dur)**(mb+1)))
    dP_b[i] = (dqc[i] - A_b[i]*Cht*(T_b[i]-Twall) - W*dV_b[i]*(P_b[i]*(1 + cv/(R)) + 2*a*b*n*n*n*cv/(R*((V_b[i])**3)) - a*n*n*cv/(R*((V_b[i])**2))))/(W*(V_b[i]*cv/(R) - cv*b*n/R))
    P_b[i+1]=(dtheta*dP_b[i])+P_b[i]
    T_b[i+1]=(P_b[i+1]*V_b[i+1])/(Ru*n)
#--------------------------------------------------------------------------#
#---------------------Modelling of Expansion-------------------------------#
#--------------------------------------------------------------------------#
theta_e=np.arange(2*pi + FIC+0.15,dtheta+ 3*pi - EVO,dtheta)
N_e = len(theta_e)
P_e=np.zeros(N_e)
dP_e=np.zeros(N_e)
T_e=np.zeros(N_e)
V_e=Vd*((r/(r-1)) - ((1-np.cos(theta_e - pi))/2) + (L/S) - 0.5*np.sqrt((2*L/S)**2 - (np.sin(theta_e - pi))**2))
dV_e=np.gradient(V_e,dtheta)
A_e= (pi*Bore**2)/2 + (pi*Bore*S/2)*(1 - np.cos(theta_e) + (L/Rad) - np.sqrt((L/Rad)**2 - (np.sin(theta_e))**2))
Twall=1100
avg_cyl_gas_velocity=2.28*CM
n = m_i[N_i-1]/MWair
T_e[0]=T_b[N_b-1]
P_e[0]=P_b[N_b-1]
for i in range(0,N_e-1):
    cv = 0.2*T_e[i] + 700
    Cht=3.26*(10**-3)*(P_e[i]**(0.8))*(T_e[i]**(-0.55))*(Bore**(-0.2))*(avg_cyl_gas_velocity**(0.8))
    dP_e[i] = -(A_e[i]*Cht*(T_e[i]-Twall) + W*dV_e[i]*(P_e[i]*(1 + cv/(R)) + 2*a*b*n*n*n*cv/(R*((V_e[i])**3)) - a*n*n*cv/(R*((V_e[i])**2))))/(W*(V_e[i]*cv/(R) - cv*b*n/R))
    P_e[i+1]=(dtheta*dP_e[i])+P_e[i]
    T_e[i+1]=(P_e[i+1]*V_e[i+1])/(Ru*n)
#--------------------------------------------------------------------------#
#---------------------Modelling of Exhuast---------------------------------#
#--------------------------------------------------------------------------#
theta_ex=np.arange(3*pi - EVO,dtheta+ 4*pi + EVC,dtheta)
N_ex = len(theta_ex)
P_ex=np.zeros(N_ex)
dP_ex=np.zeros(N_ex)
T_ex=np.zeros(N_ex)
V_ex=Vd*((r/(r-1)) - ((1-np.cos(theta_ex - pi))/2) + (L/S) - 0.5*np.sqrt((2*L/S)**2 - (np.sin(theta_ex - pi))**2))
dV_ex=np.gradient(V_ex,dtheta)
m_ex=np.zeros(N_ex)
avg_cyl_gas_velocity=6.18*CM
dm_ex=np.zeros(N_ex)
Twall=700
T_exm=1000
P_exm=101325
A_ex= (pi*Bore**2)/2 + (pi*Bore*S/2)*(1 - np.cos(theta_ex) + (L/Rad) - np.sqrt((L/Rad)**2 - (np.sin(theta_ex))**2))
P_ex[0]=P_e[N_e-1]
T_ex[0]=T_e[N_e-1]
m_ex[0]=m_i[N_i-1]

for i in range(0,N_e-1):
    cv = 0.2*T_ex[i] + 700
    cp = R+cv
    gamma=cp/cv
    CD=(107.78*((LIFT_MAX/D)**4))-(77.204*((LIFT_MAX/D)**3))+(14.1*((LIFT_MAX/D)**2))-(1.01*(LIFT_MAX/D))+0.6687
    CURTAREA=(pi/4)*D**2
    Cht=3.26*(10**-3)*(P_ex[i]**(0.8))*(T_ex[i]**(-0.55))*(Bore**(-0.2))*(avg_cyl_gas_velocity**(0.8))
    if P_exm<=P_ex[i]:
        if (P_exm/P_ex[i])<=(2/(gamma+1))**(gamma/(gamma-1)):
            dm_ex[i]=-(P_ex[i]/(np.sqrt(R*T_ex[i])))*((np.sqrt(gamma))*((2*gamma/(gamma+1))**((gamma+1)/2*(gamma-1))))*CD*CURTAREA/W
            m_ex[i+1]=(dtheta*dm_ex[i])+m_ex[i]
        else:
            dm_ex[i]=(P_ex[i]/(np.sqrt(R*T_ex[i])))*(((P_exm/P_ex[i])**(1/gamma))*np.sqrt((2*gamma/(gamma-1))*(1 - (P_exm/P_ex[i])**((gamma-1)/gamma))))*CD*CURTAREA/W
            m_ex[i+1]=(dtheta*dm_ex[i])+m_ex[i]
        n=m_ex[i]/MWair
        dP_ex[i] = ((W*dm_ex[i]*cp*T_ex[i])-A_ex[i]*Cht*(T_ex[i]-Twall)-W*dV_ex[i]*(P_ex[i]*(1 + cv/(R)) + 2*a*b*n*n*n*cv/(R*((V_ex[i])**3)) - a*n*n*cv/(R*((V_ex[i])**2))))/(W*(V_ex[i]*cv/(R) - cv*b*n/R))
        P_ex[i+1]=(dtheta*dP_ex[i])+P_ex[i]
        T_ex[i+1]=(P_ex[i+1]*V_ex[i+1])/(m_ex[i+1]*(Ru/MWair))
    else:
        if (P_ex[i]/P_exm)<=(2/(gamma+1))**(gamma/(gamma-1)):
            dm_ex[i]=(P_exm/np.sqrt(R*T_exm))*((np.sqrt(gamma))*((2*gamma/(gamma+1))**((gamma+1)/2*(gamma-1))))*CD*CURTAREA/W
            m_ex[i+1]=(dtheta*dm_ex[i])+m_ex[i]
        else:
            dm_ex[i]=(P_exm/np.sqrt(R*T_exm))*(((P_ex[i]/P_exm)**(1/gamma))*np.sqrt((2*gamma/(gamma-1))*(1 - (P_ex[i]/P_exm)**((gamma-1)/gamma))))*CD*CURTAREA/W
            m_ex[i+1]=(dtheta*dm_ex[i])+m_ex[i]
        n=m_ex[i]/MWair
        dP_ex[i] = (W*dm_ex[i]*cp*T_ex[i]-W*dV_ex[i]*(P_ex[i]*(1 + cv/(R)) + 2*a*b*n*n*n*cv/(R*((V_ex[i])**3)) - a*n*n*cv/(R*((V_ex[i])**2))))/(W*(V_ex[i]*cv/(R) - cv*b*n/R))
        P_ex[i+1]=(dtheta*dP_ex[i])+P_ex[i]
        T_ex[i+1]=(P_ex[i+1]*V_ex[i+1])/(m_ex[i+1]*(Ru/MWair))


theta = np.concatenate((theta_i,theta_c,theta_b,theta_e,theta_ex))
theta = theta*180/pi
theta = theta-360
P = np.concatenate((P_i,P_c,P_b,P_e,P_ex))
T = np.concatenate((T_i,T_c,T_b,T_e,T_ex))
V = np.concatenate((V_i,V_c,V_b,V_e,V_ex))
# plt.plot(V,P,color='b',label='withComb')
plt.xlabel('Vol ($m^3$)')
plt.ylabel('Pressure (Pa)')
plt.legend()
plt.grid()
plt.savefig('graphs/NoComb/P_V_noComb')