function[] = SI_Engine_Simulation()
clear all
clc
R=287;
g=9.81;
Ru=8314;
LHV=44.42*10^6; %in J/kg
P_atm=1E-5;
rho_atm=1.2;
gamma=1.3;
T_atm=300;
xr=0.063; %input('enter the residual gas fraction');
display('Enter dimensions in m')
Rad=34.62*10^-3; % radius of crank
L=136.5*10^-3; % length of connecting rod
r=8.5; % compression ratio
Bore=95.25*10^-3; % bore of cylinder
LIFT_MAX=4.5*10^-3; % Maximum Lift
Vd=(pi*(Bore^2)*Rad)/2;
VOL_CLEARANCE=Vd/(r-1);
display('Enter any speed between 2000 to 4500 RPM')
N=2000;
W=(2*pi*N)/60;
CM=(N*2*Rad)/30; % Mean Piston Speed
IVO=30; % IVO Before TDC
IVC=60; % IVC After BDC
EVO=0;
EVC=IVO;
D=41.2*10^-3; % Inlet Valve Diameter
D2=2*D; % Inlet manifold diameter
A1=pi*0.25*D^2;
A2=pi*0.25*D2^2;
T_m=((-0.043624*(N/1000))+1.2953)*T_atm; % Charge Heating
P_m=90000 % manifold pressure
SA=335 % Spark Timing angle
%--------------------------Global Chemical Reaction---------------%
fi=0.91; % Equivalence ratio
x=8;y=18; %x=no of carbon atoms in a mole of fuel & y=no of hydrogen atoms in a mole of fuel
Nrfuel=1;
Nro2=1*((x+(y/4))/fi);
Nrn2=3.76*((x+(y/4))/fi);
Npco2=x;
Nph2o=y/2;
Npo2=((x+(y/4))/fi)-x-(y/4);
Npn2=3.76*((x+(y/4))/fi);
Nr=Nrfuel+ Nro2+ Nrn2; % Number of moles of reactant
Np=Npco2+ Nph2o+Npo2+ Npn2;
Xrfuel=Nrfuel/Nr; % mole fraction of fuel
Xro2= Nro2/Nr;
Xrn2=Nrn2/Nr;
Xpco2=Npco2/Np;
Xph2o=Nph2o/Np;
Xpo2=Npo2/Np;
Xpn2=Npn2/Np;
MWc=12; MWh=1; MWo2=32; MWn2=28;MWrfuel=(x*MWc+y*MWh);
MWr=sqrt(Xrfuel*power((x*MWc+y*MWh),2)+Xro2*power(MWo2,2)+Xrn2*power(MWn2,2)); % Molecular weight of reactant in kg/kmol
MWp=sqrt(Xpco2*power((MWc+MWo2),2)+Xph2o*power(((2*MWh)+(MWo2/2)),2)+Xpo2*power(MWo2,2)+Xpn2*power(MWn2,2));
cprfuel=0.85; cpro2=0.918; cprn2=1.040; cppco2=0.85;cpph2o=4.18;cppo2=0.918;cppn2=1.042;
cpr=((Xrfuel*(x*MWc+y*MWh)*cprfuel+Xro2*MWo2*cpro2+Xrn2*MWn2*cprn2))*10^3 ; % specific heat of reactant in J/kmol-K
cvr=cpr-Ru;
cpp=((Xpco2*(MWc+MWo2)*cppco2+Xph2o*((2*MWh)+(MWo2/2))*cpph2o+Xpo2*MWo2*cppo2+Xpn2*MWn2*cppn2))*10^3; % specific heat of product in J/kmol-K
cvp=cpp-Ru;
%--------------------------------------------------------------------------%
%--------------Modelling of Intake Stroke----------------------------------%
%--------------------------------------------------------------------------%
gamma1=cpr/cvr;
U_m=((pi/4)*Bore^2*CM)/(A2); % Approximate velocity of air in inlet manifold
U_v=((pi/4)*Bore^2*CM)/(A1); % Approximate velocity of air at inlet Valve
U_Sound=sqrt(gamma1*(Ru/MWr)*T_m);
Mat=U_v/U_Sound;
Mam=U_m/U_Sound;
Pexhaust=101325;%input('Enter the exhaust manifold pressure=');
rho_m=rho_atm/((1+((gamma1-1)*0.5*Mam^2))^(1/(gamma1-1)));
rhot=rho_atm/((1+((gamma1-1)*0.5*Mat^2))^(1/(gamma1-1)));
SP_m=P_m+(0.5*rho_m*U_m^2); % Stagnation pressure
ST_m=T_m+((0.5*U_m^2)/(cpr/MWr)); % Stagnation Temperature
dtheta=0.0008;
if IVO>=0
theta1=(0-IVO):dtheta:0;
theta2=dtheta:dtheta:(180+IVC);
N1=floor((IVO/dtheta)+1);
N2=floor(((180+IVC-dtheta)/dtheta))+(N1+1);
else
theta1=0:dtheta:abs(IVO);
theta2=(abs(IVO)+dtheta):dtheta:(180+IVC);
N1=floor((abs(IVO)/dtheta)+1);
N2=floor((((180+IVC)-((abs(IVO)+dtheta)))/dtheta))+(N1+1);
end
theta(1,1:N1)=theta1;
theta(1,N1+1:N2)=theta2;
V_v=zeros(N2,1);
dV=zeros(N2,1);
for i=1:N1 % crank angle from IVO to TDC
V_v(i,1)=((Vd/(r-1))+((Vd/2)*(1+(L/Rad)-cosd(theta1(1,i))-((L/Rad)^2-(sind(theta1(1,i)))^2)^0.5))); % Incylinder volume (Cubic meter ) calculation
end
for i=(N1+1):N2 % crank angle from TDC to IVC
V_v(i,1)=((Vd/(r-1))+((Vd/2)*(1+(L/Rad)-cosd(theta2(1,(i-N1)))-((L/Rad)^2-(sind(theta2(1,i-N1)))^2)^0.5)));
end
dV=gradient(V_v,dtheta);
P_v=zeros(N2,1);
T_v=zeros(N2,1);
m=zeros(N2,1);
rho_v=zeros(N2,1);
dm=zeros(N2,1);
dP_v=zeros(N2,1);
P_v(1,1)=P_m; % Initial Value
T_v(1,1)=T_m;
m(1,1)=(P_v(1,1)*V_v(1,1))/((Ru/MWr)*T_v(1,1));
rho_v(1,1)=m(1,1)/V_v(1,1);
Uavg_v=0;                   
CF_Avg=0;
for i=1:N2-1
si=(pi*(IVO-IVC+540+(2*theta(1,i))))/(IVO+IVC+180);
LIFT=LIFT_MAX*(1+cos(si))/2;
CD=(107.78*((LIFT/D)^4))-(77.204*((LIFT/D)^3))+(14.1*((LIFT/D)^2))-(1.01*(LIFT/D))+0.6687;
CURTAREA=pi*D*LIFT;
if LIFT>=D/4
CURTAREA=(pi/4)*D^2;
end
PORTAREA=(pi/4)*D^2;
CF=CURTAREA/PORTAREA;
CF_Avg=(CF*dtheta)+CF_Avg;
if (P_v(i,1)<=SP_m)
if (P_v(i,1)/SP_m)>=(2/(gamma1+1))^(gamma1/(gamma1-1)) % subsonic condition
dm(i,1)=((((CD*CURTAREA*SP_m)/(sqrt(R*ST_m)))*((P_v(i,1)/SP_m)^(1/gamma1))*(((2*gamma1)/(gamma1-1))*(1-((P_v(i,1)/SP_m)^((gamma1-1)/gamma1))))^0.5)/W);
m(i+1,1)=(dtheta*dm(i,1))+m(i,1);
U_v=(dm(i,1)*W)/(rhot*A1);
U_m=(dm(i,1)*W)/(rho_m*A2);
Mat=U_v/U_Sound;
Mam=U_m/U_Sound;
rhot=rho_atm/((1+((gamma1-1)*0.5*Mat^2))^(1/(gamma1-1)));
rho_m=rho_atm/((1+((gamma1-1)*0.5*Mam^2))^(1/(gamma1-1))); 
SP_m=P_m+(0.5*rho_m*U_m^2);
end
if (P_v(i,1)/SP_m)<(2/(gamma1+1))^(gamma1/(gamma1-1)) %sonic condition
dm(i,1)=(((CD*CURTAREA*SP_m)/(sqrt(R*ST_m)))*((gamma1)^0.5)*(2/(gamma1+1))^((gamma1+1)/(2*(gamma1-1))))/W;
m(i+1,1)=(dtheta*dm(i,1))+m(i,1);
U_v=(dm(i,1)*W)/(rhot*A1);
U_m=(dm(i,1)*W)/(rho_m*A2);
Mat=U_v/U_Sound;
Mam=U_m/U_Sound;
rhot=rho_atm/((1+((gamma1-1)*0.5*Mat^2))^(1/(gamma1-1)));
rho_m=rho_atm/((1+((gamma1-1)*0.5*Mam^2))^(1/(gamma1-1)));
SP_m=P_m+(0.5*rho_m*U_m^2);
end
dP_v(i,1)=(((gamma1-1)*(dm(i,1)*(cpr/MWr)*T_v(i,1)))-(gamma1*P_v(i,1)*dV(i,1)))/V_v(i,1);
P_v(i+1,1)=(dtheta*dP_v(i,1))+P_v(i,1);
T_v(i+1,1)=(P_v(i+1,1)*V_v(i+1,1))/(m(i+1,1)*(Ru/MWr));
rho_v(i+1,1)= (P_v(i+1,1)*MWr)/(Ru*T_v(i+1,1));
else
if (SP_m/P_v(i,1))>=(2/(gamma1+1))^(gamma1/(gamma1-1)) % subsonic condition
dm(i,1)=((((CD*CURTAREA*P_v(i,1))/(sqrt(R*T_v(i,1))))*((SP_m/P_v(i,1))^(1/gamma1))*(((2*gamma1)/(gamma1-1))*(1-((SP_m/P_v(i,1))^((gamma1-1)/gamma1))))^0.5)/W);
m(i+1,1)=(dtheta*-dm(i,1))+m(i,1);
U_v=(-dm(i,1)*W)/(rhot*A1);
U_m=(-dm(i,1)*W)/(rho_m*A2);
Mat=U_v/U_Sound;
Mam=U_m/U_Sound;
rhot=rho_atm/((1+((gamma1-1)*0.5*Mat^2))^(1/(gamma1-1)));
rho_m=rho_atm/((1+((gamma1-1)*0.5*Mam^2))^(1/(gamma1-1)));
SP_m=P_m+(0.5*rho_m*U_m^2);
end
if (SP_m/P_v(i,1))<(2/(gamma1+1))^(gamma1/(gamma1-1)) %sonic condition
dm(i,1)=(((CD*CURTAREA*P_v(i,1))/(sqrt(R*T_v(i,1))))*((gamma1)^0.5)*(2/(gamma1+1))^((gamma1+1)/(2*(gamma1-1))))/W;
m(i+1,1)=(dtheta*-dm(i,1))+m(i,1);
U_v=(-dm(i,1)*W)/(rhot*A1);
U_m=(-dm(i,1)*W)/(rho_m*A2);
Mat=U_v/U_Sound;
Mam=U_m/U_Sound;
rhot=rho_atm/((1+((gamma1-1)*0.5*Mat^2))^(1/(gamma1-1)));
rho_m=rho_atm/((1+((gamma1-1)*0.5*Mam^2))^(1/(gamma1-1)));
SP_m=P_m+(0.5*rho_m*U_m^2);
end
dP_v(i,1)=(((gamma1-1)*(-dm(i,1)*(cpr/MWr)*T_v(i,1)))-(gamma1*P_v(i,1)*dV(i,1)))/V_v(i,1);
P_v(i+1,1)=(dtheta*dP_v(i,1))+P_v(i,1);
T_v(i+1,1)=(P_v(i+1,1)*V_v(i+1,1))/(m(i+1,1)*(Ru/MWr));
rho_v(i+1,1)= (P_v(i+1,1)*MWr)/(Ru*T_v(i+1,1));
end
Uavg_v=Uavg_v+U_v;
end
Uavg_v=Uavg_v/(N2-1); % Average velocity at the inlet valve for turbulent velocity (Ut) calculations.
CF_Avg=(1/(IVO+180+IVC))*CF_Avg; % average flow coefficient.
MI=(CM/U_Sound)*(Bore/D)^2*(1/CF_Avg); % Mach Index
mivc=m(N2,1);
mf=(mivc)/((14.7/fi)+1);
ma=mivc-mf;
VOL_EFF=ma/(rho_atm*Vd);
P_v(N2)
T_v(N2)
m(N2)
%-------------------Modelling of Compression Stroke-----------------------%
%-------------------------------------------------------------------------%
Tw=400;
dtheta=0.005; % dtheta value should not increase above 0.010
thetac=180+IVC:dtheta:SA;
Nc=floor((SA-(180+IVC))/dtheta)+1;
V_comp=zeros(Nc,1);
Aw=zeros(Nc,1);
Ru=8314;
for i=1:Nc % crank angle from IVC to Spark Timing
V_comp(i,1)=((Vd/(r-1))+((Vd/2)*(1+(L/Rad)-cosd(thetac(1,i))-((L/Rad)^2-(sind(thetac(1,i)))^2)^0.5))); % Incylinder volume (Cubic meter ) calculation
Aw(i,1)=(pi*0.5*Bore^2)+((pi*Bore*Rad)*((L/Rad)+1-cosd(thetac(1,i))+((L/Rad)^2-(sind(thetac(1,i)))^2)^0.5));
end
dVc=gradient(V_comp,dtheta); % gradient of cylinder volume
P_comp=zeros(Nc,1);
T_comp=zeros(Nc,1);
rho=zeros(Nc,1);
A_comp=zeros(Nc,1);
Q_comp=zeros(Nc,1);
Sr=zeros(Nc,1);
rho(1,1)=mivc/V_comp(1,1);
T_comp(1,1)=T_v(N2,1);
P_comp(1,1)=rho(1,1)*(Ru/MWr)*T_comp(1,1);
for i=1:Nc-1 % Compression Loop
if T_comp(i,1)<1000
c1o2=0.03212936E2;c2o2=0.112748E-2;c3o2=-0.057561E-5;c4o2=0.131387E-8;c5o2=-0.0876855E-11;c6o2=-0.1005249E4;c7o2=0.060347E2;
c1n2=0.0329867E2;c2n2=0.140824E-2;c3n2=-0.0396322E-4;c4n2=0.0564151E-7;c5n2=-0.0244485E-10;c6n2=-0.1020899E4;c7n2=0.0395037E2;
else
c1o2=0.036975E2;c2o2=0.0613519E-2;c3o2=-0.1258842E-6;c4o2=0.0177528E-9;c5o2=-0.1136435E-14;c6o2=-0.1233930E4;c7o2=0.0318916E2;
c1n2=0.0292664E2;c2n2=0.1487976E-2;c3n2=-0.0568476E-5;c4n2=0.10097038E-9;c5n2=-0.0675335E-13;c6n2=-0.0922797E4;c7n2=0.0598052E2;
end
if T_comp(i,1)<1396
c1f=-4.20868893E+00;c2f=1.11440581E-01;c3f=-7.91346582E-05;c4f=2.92406242E-08;c5f=-4.43743191E-12;c6f=-2.99446875E+04;c7f=4.49521701E+01;
else
c1f=2.71373590E+01;c2f=3.79004890E-02;c3f=-1.29437358E-05;c4f=2.00760372E-09;c5f=-1.16400580E-13;c6f=-4.07958177E+04;c7f=-1.23277495E+02;
end
cprfuel=Ru*(c1f+(c2f*T_comp(i,1))+(c3f*(T_comp(i,1))^2)+(c4f*(T_comp(i,1))^3)+(c5f*(T_comp(i,1))^4));
cpro2=Ru*(c1o2+(c2o2*T_comp(i,1))+(c3o2*(T_comp(i,1))^2)+(c4o2*(T_comp(i,1))^3)+(c5o2*(T_comp(i,1))^4));
cprn2=Ru*(c1n2+(c2n2*T_comp(i,1))+(c3n2*(T_comp(i,1))^2)+(c4n2*(T_comp(i,1))^3)+(c5n2*(T_comp(i,1))^4));
cpr=((Xrfuel*cprfuel+Xro2*cpro2+Xrn2*cprn2)); % specific heat of reactant in J/kmol-K
cvr=cpr-Ru;
Hrfuel=Ru*T_comp(i,1)*(c1f+((c2f/2)*T_comp(i,1))+((c3f/3)*(T_comp(i,1))^2)+((c4f/4)*(T_comp(i,1))^3)+((c5f/5)*(T_comp(i,1))^4)+(c6f/T_comp(i,1)));
Hro2=Ru*T_comp(i,1)*(c1o2+((c2o2/2)*T_comp(i,1))+((c3o2/3)*(T_comp(i,1))^2)+((c4o2/4)*(T_comp(i,1))^3)+((c5o2/5)*(T_comp(i,1))^4)+(c6o2/T_comp(i,1)));
Hrn2=Ru*T_comp(i,1)*(c1n2+((c2n2/2)*T_comp(i,1))+((c3n2/3)*(T_comp(i,1))^2)+((c4n2/4)*(T_comp(i,1))^3)+((c5n2/5)*(T_comp(i,1))^4)+(c6n2/T_comp(i,1)));
Hr=((Xrfuel*Hrfuel+Xro2*Hro2+Xrn2*Hrn2));
Srfuel=(Ru*((c1f*log(T_comp(i,1)))+(c2f*T_comp(i,1))+((c3f/2)*(T_comp(i,1))^2)+((c4f/3)*(T_comp(i,1))^3)+((c5f/4)*(T_comp(i,1))^4)+c7f))-(Ru*log((Xrfuel*P_comp(i,1))/P_atm));
Sro2=(Ru*((c1o2*log(T_comp(i,1)))+(c2o2*T_comp(i,1))+((c3o2/2)*(T_comp(i,1))^2)+((c4o2/3)*(T_comp(i,1))^3)+((c5o2/4)*(T_comp(i,1))^4)+c7o2))-(Ru*log((Xro2*P_comp(i,1))/P_atm));
Srn2=(Ru*((c1n2*log(T_comp(i,1)))+(c2n2*T_comp(i,1))+((c3n2/2)*(T_comp(i,1))^2)+((c4n2/3)*(T_comp(i,1))^3)+((c5n2/4)*(T_comp(i,1))^4)+c7n2))-(Ru*log((Xrn2*P_comp(i,1))/P_atm));
Sr(i,1)=((Xrfuel*Srfuel+Xro2*Sro2+Xrn2*Srn2));
drhoP=(MWr/(Ru*T_comp(i,1)))*(cvr/(cvr+Ru));
drhoT=(rho(i,1)*cvr)/(Ru*T_comp(i,1));
A=((1/rho(i,1))*(drhoT/drhoP))+(cpr/MWr);
B=1/drhoP;
Ku=0.0243; %thermal conductivity of air
dQ_comp=(Aw(i,1)*(((0.45*Ku*((80000)^0.75)*(T_comp(i,1)-Tw))/Bore)+((4.3*10^-9)*((T_comp(i,1))^4-Tw^4))))/W; % Heat loss in J/radian
f_T1=((B/A)*((-dVc(i,1)/V_comp(i,1))-(dQ_comp/(B*mivc))));
T_comp(i+1,1)=(f_T1*dtheta)+T_comp(i,1); % Temperature during compression
rho(i+1,1)=mivc/V_comp(i+1,1);
P_comp(i+1,1)=(rho(i+1,1)*(Ru/MWr)*T_comp(i+1,1)); % Pressure during compression
end
Srfuel=(Ru*((c1f*log(T_comp(Nc,1)))+(c2f*T_comp(Nc,1))+((c3f/2)*(T_comp(Nc,1))^2)+((c4f/3)*(T_comp(Nc,1))^3)+((c5f/4)*(T_comp(Nc,1))^4)+c7f))-(Ru*log((Xrfuel*P_comp(Nc,1))/P_atm));
Sro2=(Ru*((c1o2*log(T_comp(Nc,1)))+(c2o2*T_comp(Nc,1))+((c3o2/2)*(T_comp(Nc,1))^2)+((c4o2/3)*(T_comp(Nc,1))^3)+((c5o2/4)*(T_comp(Nc,1))^4)+c7o2))-(Ru*log((Xro2*P_comp(Nc,1))/P_atm));
Srn2=(Ru*((c1n2*log(T_comp(Nc,1)))+(c2n2*T_comp(Nc,1))+((c3n2/2)*(T_comp(Nc,1))^2)+((c4n2/3)*(T_comp(Nc,1))^3)+((c5n2/4)*(T_comp(Nc,1))^4)+c7n2))-(Ru*log((Xrn2*P_comp(Nc,1))/P_atm));
Sr(Nc,1)=((Xrfuel*Srfuel+Xro2*Sro2+Xrn2*Srn2));
%--------------------------------------------------------------------%
%----------------Modeling of combustion------------------------------%
%--------------------------------------------------------------------%
dtheta=0.005; %Step size for Combustion
comb_dura=((-1.6189*(N/1000)^2)+(19.886*(N/1000))+39.951);
thetaEOB=(SA)+comb_dura;
thetacb=(SA):dtheta:thetaEOB;
Ncb=floor((thetaEOB-(SA))/dtheta)+1;
V_cb=zeros(Ncb,1);
Aw_cb=zeros(Ncb,1);
double(V_cb);
double(Aw_cb);
for i=1:Ncb
V_cb(i,1)=((Vd/(r-1))+((Vd/2)*(1+(L/Rad)-cosd(thetacb(1,i))-((L/Rad)^2-(sind(thetacb(1,i)))^2)^0.5))); % Incylinder volume (Cubic meter ) calculation
Aw_cb(i,1)=(pi*0.5*Bore^2)+((pi*Bore*Rad)*((L/Rad)+1-cosd(thetacb(1,i))+((L/Rad)^2-(sind(thetacb(1,i)))^2)^0.5));
end
dVcb=gradient(V_cb,dtheta); % gradient of cylinder volume
f=zeros(Ncb,1);
Vu_cb=zeros(Ncb,1);
Vb_cb=zeros(Ncb,1);
double(f);
double(Vu_cb);
double(Vb_cb);
a=5; m1=3;j1=0;j2=0; % for fixed com_dura m=2.1, for variable comb_dura m=2.2, w/o HT m=2.0
for i=1:Ncb
f(i,1)=1-(exp(-a*((thetacb(1,i)-(SA))/comb_dura)^m1));
Vu_cb(i,1)=(1-f(i,1))*V_cb(i,1);
Vb_cb(i,1)=f(i,1)*V_cb(i,1);
if f(i,1)<=0.1
j1=j1+1;
end
if f(i,1)>0.1&&f(i,1)<=0.95
j2=j2+1;
end
end
thetad=(j1-1)*dtheta;
thetab=(j2-1)*dtheta;
df=gradient(f,dtheta);
P_cb=zeros(Ncb,1);rhou=zeros(Ncb,1);
Pu_cb=zeros(Ncb,1);rhob=zeros(Ncb,1);
Tu_cb=zeros(Ncb,1);mu=zeros(Ncb,1);
Pb_cb=zeros(Ncb,1);mb=zeros(Ncb,1);
Tb_cb=zeros(Ncb,1);me=zeros(Ncb,1);
dmb=zeros(Ncb,1);f_Pu1=zeros(Ncb,1);
dmu=zeros(Ncb,1);f_Tu1=zeros(Ncb,1);
f_Pb1=zeros(Ncb,1);Af=zeros(Ncb,1);
f_Tb1=zeros(Ncb,1);Aw2=zeros(Ncb,1);
A_cb=zeros(Ncb,1);Au_cb=zeros(Ncb,1);
Ab_cb=zeros(Ncb,1);Vf=zeros(Ncb,1);
Sr=zeros(Ncb,1);Sp=zeros(Ncb,1);
double(P_cb); double(Pu_cb);double(Pb_cb);double(Tu_cb);
double(Tb_cb);double(rhou);double(rhob); double(mu);
double(mb); double(me);double(dmb); double(dmu);
double(f_Pu1); double(f_Pb1);double(f_Tu1); double(f_Tb1);
rhou(1,1)=rho(Nc,1);
rhob(1,1)=1E-10;
Pu_cb(1,1)=P_comp(Nc,1);
Pb_cb(1,1)=1.1E-10;
Tu_cb(1,1)=T_comp(Nc,1);
Tb_cb(1,1)=T_comp(Nc,1);
P_cb(1,1)=Pu_cb(1,1);
Au_cb(1,1)=A_comp(Nc,1);
Ab_cb(1,1)=1E-10;
A_cb(1,1)=A_comp(Nc,1);
mb(2,1)=f(2,1)*mivc;
mu(2,1)=mivc-mb(2,1);
me(2,1)=mb(2,1);
Tu_cb(2,1)=T_comp(Nc,1)+(1E-5);
Pr=0;
Tad = 2583;
Tb_cb(2,1)=Tad;
Pu_cb(2,1)=((mu(2,1)*(Ru/MWr)*Tu_cb(2,1))/Vu_cb(2,1));
Pb_cb(2,1)=1.1E-5;
P_cb(2,1)=Pu_cb(2,1)+Pb_cb(2,1)+Pr;
Au_cb(2,1)=(1-f(2,1))*A_comp(Nc,1);
Ab_cb(2,1)=f(2,1)*A_comp(Nc,1);
A_cb(2,1)=Au_cb(2,1)+Ab_cb(2,1);
for i=2:Ncb-1 %Starting from SA+dtheta
dVu_cb=((1-f(i,1))*dVcb(i,1))+(V_cb(i,1)*-df(i,1));
dVb_cb=(df(i,1)*V_cb(i,1))+(f(i,1)*dVcb(i,1));
double(dVu_cb); double(dVb_cb);
if mb(i,1)<=1E-9&&Vb_cb(i,1)<=1E-10
rhob(i,1)=1E-8;
else if mb(i,1)<=1E-9
rhob(i,1)=1E-8;
else
rhob(i,1)=(mb(i,1))/Vb_cb(i,1);
end
end
if mu(i,1)<=1E-9&&Vu_cb(i,1)<=1E-10
rhou(i,1)=1E-8;
else if mu(i,1)<=1E-9
rhou(i,1)=1E-8;
else
rhou(i,1)=mu(i,1)/Vu_cb(i,1);
end
end
Vf(i,1)=Vb_cb(i,1)+((me(i,1)-mb(i,1))/rhou(i,1));
if Vf<=0.0
Vf=1E-10;
end
hz=(V_cb(i,1)*4)/(pi*Bore^2); % Chamber height
rc=51.5*10^-3; % spark plug location from the edge of the cylinder.
rf=((thetacb(1,i)-SA-thetad*(1-exp((-(thetacb(1,i)-SA)/thetad))))/thetab)*(Bore/2); % Flame front radius
if rf<=hz
x=rf; % Flame depth
else
x=hz;
end
if rf<=rc
Af(i,1)=2*pi*rf*x;
Aw2(i,1)=0;
end
if rf>rc
y=0:x/100:x;
S=zeros(101,1);
Aw2c=zeros(101,1);
for j=1:101
fy=sqrt(rf^2-(y(1,j))^2);
if fy<=rc
p=2*pi*fy;
q=0;
S(j,1)=((rf*p)/fy)*(x/100);
Aw2c(j,1)=q*(x/100);
end
if fy>=(Bore-rc)
p=0;
q=pi*Bore;
S(j,1)=((rf*p)/fy)*(x/100);
Aw2c(j,1)=q*(x/100);
end
if fy>rc&&fy<(Bore-rc)
p=2*(pi-Bore)*fy;
beta1=acos(1+(((rc/Bore)^2-(fy/Bore)^2)/(0.5-(rc/Bore))));
q=Bore*beta1;
S(j,1)=((rf*p)/fy)*(x/100);
Aw2c(j,1)=q*(x/100);
end
end
Af(i,1)=sum(S);
Aw2(i,1)=sum(Aw2c);
end
Awb_cb=Aw2(i,1);
Awu_cb=Aw_cb(i,1)-Aw2(i,1);
double(Awu_cb);double(Awb_cb);
alpha=2.18-(0.8*(fi-1));
beta=-0.16+(0.22*(fi-1));
Sl0=0.305+(-0.549*(fi-1.21)^2);
if i<=j1
Sl=0.13;
else
Sl=Sl0*((Tu_cb(i,1)/298)^alpha)*((P_cb(i,1)/101300)^beta);
end
if i<=j1
Ut=1E-5;
Lt=1E-5;
else
Ut=0.08*Uavg_v*sqrt((rhou(i,1)*N2)/sum(rho_v));
Lt=0.8*LIFT_MAX*((sum(rho_v)/(N2*rhou(i,1)))^(3/4));
end
dmb(i,1)=((rhou(i,1)*Af(i,1)*Sl)+(((me(i,1)-mb(i,1))*Sl)/Lt))/W;
dme=(rhou(i,1)*Af(i,1)*(Ut+Sl))/W;
mb(i+1,1)=(dmb(i,1)*dtheta)+mb(i,1);
me(i+1,1)=(dme*dtheta)+me(i,1);
dmu(i,1)=-dmb(i,1);
mu(i+1,1)=mu(i,1)+(dmu(i,1)*dtheta);
if mb(i+1,1)<=1E-9&&Vb_cb(i+1,1)<=1E-10
rhob(i+1,1)=1E-8;
else if mb(i+1,1)<=1E-9
rhob(i+1,1)=1E-8;
else
rhob(i+1,1)=(mb(i+1,1))/Vb_cb(i+1,1);
end
end
if mu(i+1,1)<=1E-9&&Vu_cb(i+1,1)<=1E-10
rhou(i+1,1)=1E-8;
else if mu(i+1,1)<=1E-9
rhou(i+1,1)=1E-8;
else
rhou(i+1,1)=mu(i+1,1)/Vu_cb(i+1,1);
end
end
% Unburned Mixture calculation
if Tu_cb(i,1)<1000
c1o2=0.03212936E2;c2o2=0.112748E-2;c3o2=-0.057561E-5;c4o2=0.131387E-8;c5o2=-0.0876855E-11;c6o2=-0.1005249E4;c7o2=0.060347E2;
c1n2=0.0329867E2;c2n2=0.140824E-2;c3n2=-0.0396322E-4;c4n2=0.0564151E-7;c5n2=-0.0244485E-10;c6n2=-0.1020899E4;c7n2=0.0395037E2;
else
c1o2=0.036975E2;c2o2=0.0613519E-2;c3o2=-0.1258842E-6;c4o2=0.0177528E-9;c5o2=-0.1136435E-14;c6o2=-0.1233930E4;c7o2=0.0318916E2;
c1n2=0.0292664E2;c2n2=0.1487976E-2;c3n2=-0.0568476E-5;c4n2=0.10097038E-9;c5n2=-0.0675335E-13;c6n2=-0.0922797E4;c7n2=0.0598052E2;
end
if Tu_cb(i,1)<1396
c1f=-4.20868893E+00;c2f=1.11440581E-01;c3f=-7.91346582E-05;c4f=2.92406242E-08;c5f=-4.43743191E-12;c6f=-2.99446875E+04;c7f=4.49521701E+01;
else
c1f=2.71373590E+01;c2f=3.79004890E-02;c3f=-1.29437358E-05;c4f=2.00760372E-09;c5f=-1.16400580E-13;c6f=-4.07958177E+04;c7f=-1.23277495E+02;
end
cprfuel=Ru*(c1f+(c2f*Tu_cb(i,1))+(c3f*(Tu_cb(i,1))^2)+(c4f*(Tu_cb(i,1))^3)+(c5f*(Tu_cb(i,1))^4));
cpro2=Ru*(c1o2+(c2o2*Tu_cb(i,1))+(c3o2*(Tu_cb(i,1))^2)+(c4o2*(Tu_cb(i,1))^3)+(c5o2*(Tu_cb(i,1))^4));
cprn2=Ru*(c1n2+(c2n2*Tu_cb(i,1))+(c3n2*(Tu_cb(i,1))^2)+(c4n2*(Tu_cb(i,1))^3)+(c5n2*(Tu_cb(i,1))^4));
cpr=((Xrfuel*cprfuel+Xro2*cpro2+Xrn2*cprn2)); % specific heat of reactants in J/kmol-K
cvr=cpr-Ru;
Hrfuel=Ru*Tu_cb(i,1)*(c1f+((c2f/2)*Tu_cb(i,1))+((c3f/3)*(Tu_cb(i,1))^2)+((c4f/4)*(Tu_cb(i,1))^3)+((c5f/5)*(Tu_cb(i,1))^4)+(c6f/Tu_cb(i,1)));
Hro2=Ru*Tu_cb(i,1)*(c1o2+((c2o2/2)*Tu_cb(i,1))+((c3o2/3)*(Tu_cb(i,1))^2)+((c4o2/4)*(Tu_cb(i,1))^3)+((c5o2/5)*(Tu_cb(i,1))^4)+(c6o2/Tu_cb(i,1)));
Hrn2=Ru*Tu_cb(i,1)*(c1n2+((c2n2/2)*Tu_cb(i,1))+((c3n2/3)*(Tu_cb(i,1))^2)+((c4n2/4)*(Tu_cb(i,1))^3)+((c5n2/5)*(Tu_cb(i,1))^4)+(c6n2/Tu_cb(i,1)));
Hr=((Xrfuel*Hrfuel+Xro2*Hro2+Xrn2*Hrn2)); % specific enthalpy of reactants in J/kmol
Srfuel=(Ru*((c1f*log(Tu_cb(i,1)))+(c2f*Tu_cb(i,1))+((c3f/2)*(Tu_cb(i,1))^2)+((c4f/3)*(Tu_cb(i,1))^3)+((c5f/4)*(Tu_cb(i,1))^4)+c7f))-(Ru*log((Xrfuel*Pu_cb(i,1))/P_atm));
Sro2=(Ru*((c1o2*log(Tu_cb(i,1)))+(c2o2*Tu_cb(i,1))+((c3o2/2)*(Tu_cb(i,1))^2)+((c4o2/3)*(Tu_cb(i,1))^3)+((c5o2/4)*(Tu_cb(i,1))^4)+c7o2))-(Ru*log((Xro2*Pu_cb(i,1))/P_atm));
Srn2=(Ru*((c1n2*log(Tu_cb(i,1)))+(c2n2*Tu_cb(i,1))+((c3n2/2)*(Tu_cb(i,1))^2)+((c4n2/3)*(Tu_cb(i,1))^3)+((c5n2/4)*(Tu_cb(i,1))^4)+c7n2))-(Ru*log((Xrn2*Pu_cb(i,1))/P_atm));
Sr(i,1)=((Xrfuel*Srfuel+Xro2*Sro2+Xrn2*Srn2)); % specific entropy of reactants in J/kmol-k
if i==1
drhoPu=(rho(Nc,1)-rho((Nc-1),1))/(P_comp(Nc,1)-P_comp((Nc-1),1));
drhoTu=(rho(Nc,1)-rho((Nc-1),1))/(T_comp(Nc,1)-T_comp((Nc-1),1));
else
if abs(rhou(i,1)-rhou(i-1,1))<=1E-20
drhoPu=1E-10;
drhoTu=1E-10;
else
drhoPu=(rhou(i,1)-rhou((i-1),1))/(Pu_cb(i,1)-Pu_cb((i-1),1));
drhoTu=(rhou(i,1)-rhou((i-1),1))/(Tu_cb(i,1)-Tu_cb((i-1),1));
end
end
A=((1/rhou(i,1))*(drhoTu/drhoPu))+(cpr/MWr);
B=1/drhoPu;
if A>=1E20
A=1E20;
end
if A<=-1E20
A=-1E20;
end
if B>=1E20
B=1E20;
end
if B<=-1E20
B=-1E20;
end
hu=Hr/MWr; %Specific enthalpy of unburned gases in J/kg
Ku=0.0243; %thermal conductivity of air
dQu=(Awu_cb*(((0.45*Ku*((80000)^0.75)*(Tu_cb(i,1)-Tw))/Bore)+((4.3*10^-9)*((Tu_cb(i,1))^4-Tw^4))))/W; % Heat loss in J/radian
f_Tu1(i,1)=((B/A)*((-dVu_cb/Vu_cb(i,1))+((-dQu)/(B*mu(i,1)))+(dmu(i,1)/mu(i,1))));
Tu_cb(i+1,1)=(f_Tu1(i,1)*dtheta)+Tu_cb(i,1); % Temperature of unburned gases
Pu_cb(i+1,1)=rhou(i+1,1)*(Ru/MWr)*Tu_cb(i+1,1); % Pressure of unburned gases
% Burned Mixture Calculation
if Tb_cb(i,1)<1000
c1o2=0.03212936E2;c2o2=0.112748E-2;c3o2=-0.057561E-5;c4o2=0.131387E-8;c5o2=-0.0876855E-11;c6o2=-0.1005249E4;c7o2=0.060347E2;
c1n2=0.0329867E2;c2n2=0.140824E-2;c3n2=-0.0396322E-4;c4n2=0.0564151E-7;c5n2=-0.0244485E-10;c6n2=-0.1020899E4;c7n2=0.0395037E2;
c1co2=0.0227572E2;c2co2=0.0992207E-1;c3co2=-0.1040911E-4;c4co2=0.06866686E-7;c5co2=-0.0211728E-10;c6co2=-0.0483731E6;c7co2=0.1018848E2;
c1h2o=0.0338684E2;c2h2o=0.0347498E-1;c3h2o=-0.0635469E-4;c4h2o=0.0696858E-7;c5h2o=-0.0250658E-10;c6h2o=-0.0302081E6;c7h2o=0.0259023E2;
else
c1o2=0.036975E2;c2o2=0.0613519E-2;c3o2=-0.1258842E-6;c4o2=0.0177528E-9;c5o2=-0.1136435E-14;c6o2=-0.1233930E4;c7o2=0.0318916E2;
c1n2=0.0292664E2;c2n2=0.1487976E-2;c3n2=-0.0568476E-5;c4n2=0.10097038E-9;c5n2=-0.0675335E-13;c6n2=-0.0922797E4;c7n2=0.0598052E2;
c1co2=0.0445362E2;c2co2=0.0314016E-1;c3co2=-0.1278410E-5;c4co2=0.0239399E-8;c5co2=-0.1669033E-13;c6co2=-0.0489669E6;c7co2=-0.0955395E1;
c1h2o=0.02672145E2;c2h2o=0.0305629E-1;c3h2o=-0.0873026E-5;c4h2o=0.1200996E-9;c5h2o=-0.06391618E-13;c6h2o=-0.02989921E6;c7h2o=0.06862817E2;
end
cppco2=Ru*(c1co2+(c2co2*Tb_cb(i,1))+(c3co2*(Tb_cb(i,1))^2)+(c4co2*(Tb_cb(i,1))^3)+(c5co2*(Tb_cb(i,1))^4));
cpph2o=Ru*(c1h2o+(c2h2o*Tb_cb(i,1))+(c3h2o*(Tb_cb(i,1))^2)+(c4h2o*(Tb_cb(i,1))^3)+(c5h2o*(Tb_cb(i,1))^4));
cppo2=Ru*(c1o2+(c2o2*Tb_cb(i,1))+(c3o2*(Tb_cb(i,1))^2)+(c4o2*(Tb_cb(i,1))^3)+(c5o2*(Tb_cb(i,1))^4));
cppn2=Ru*(c1n2+(c2n2*Tb_cb(i,1))+(c3n2*(Tb_cb(i,1))^2)+(c4n2*(Tb_cb(i,1))^3)+(c5n2*(Tb_cb(i,1))^4));
cpp=((Xpco2*cppco2+Xph2o*cpph2o+Xpo2*cppo2+Xpn2*cppn2)); % specific heat of product in J/kmol-K
cvp=cpp-Ru;
Hpo2=Ru*Tb_cb(i,1)*(c1o2+((c2o2/2)*Tb_cb(i,1))+((c3o2/3)*(Tb_cb(i,1))^2)+((c4o2/4)*(Tb_cb(i,1))^3)+((c5o2/5)*(Tb_cb(i,1))^4)+(c6o2/Tb_cb(i,1)));
Hpn2=Ru*Tb_cb(i,1)*(c1n2+((c2n2/2)*Tb_cb(i,1))+((c3n2/3)*(Tb_cb(i,1))^2)+((c4n2/4)*(Tb_cb(i,1))^3)+((c5n2/5)*(Tb_cb(i,1))^4)+(c6n2/Tb_cb(i,1)));
Hpco2=Ru*Tb_cb(i,1)*(c1co2+((c2co2/2)*Tb_cb(i,1))+((c3co2/3)*(Tb_cb(i,1))^2)+((c4co2/4)*(Tb_cb(i,1))^3)+((c5co2/5)*(Tb_cb(i,1))^4)+(c6co2/Tb_cb(i,1)));
Hph2o=Ru*Tb_cb(i,1)*(c1h2o+((c2h2o/2)*Tb_cb(i,1))+((c3h2o/3)*(Tb_cb(i,1))^2)+((c4h2o/4)*(Tb_cb(i,1))^3)+((c5h2o/5)*(Tb_cb(i,1))^4)+(c6h2o/Tb_cb(i,1)));
Hp=((Xpco2*Hpco2+Xph2o*Hph2o+Xpo2*Hpo2+Xpn2*Hpn2)); % specific enthalpy of products in J/kmol
Spo2=(Ru*((c1o2*log(Tb_cb(i,1)))+(c2o2*Tb_cb(i,1))+((c3o2/2)*(Tb_cb(i,1))^2)+((c4o2/3)*(Tb_cb(i,1))^3)+((c5o2/4)*(Tb_cb(i,1))^4)+c7o2))-(Ru*log((Xpo2*Pb_cb(i,1))/P_atm));
Spn2=(Ru*((c1n2*log(Tb_cb(i,1)))+(c2n2*Tb_cb(i,1))+((c3n2/2)*(Tb_cb(i,1))^2)+((c4n2/3)*(Tb_cb(i,1))^3)+((c5n2/4)*(Tb_cb(i,1))^4)+c7n2))-(Ru*log((Xpn2*Pb_cb(i,1))/P_atm));
Spco2=(Ru*((c1co2*log(Tb_cb(i,1)))+(c2co2*Tb_cb(i,1))+((c3co2/2)*(Tb_cb(i,1))^2)+((c4co2/3)*(Tb_cb(i,1))^3)+((c5co2/4)*(Tb_cb(i,1))^4)+c7co2))-(Ru*log((Xpco2*Pb_cb(i,1))/P_atm));
Sph2o=(Ru*((c1h2o*log(Tb_cb(i,1)))+(c2h2o*Tb_cb(i,1))+((c3h2o/2)*(Tb_cb(i,1))^2)+((c4h2o/3)*(Tb_cb(i,1))^3)+((c5h2o/4)*(Tb_cb(i,1))^4)+c7h2o))-(Ru*log((Xph2o*Pb_cb(i,1))/P_atm));
Sp(i,1)=((Xpco2*Spco2+Xph2o*Sph2o+Xpo2*Spo2+Xpn2*Spn2)); % specific entropy of products in J/kmol-K
if (Tb_cb(i,1)-Tb_cb((i-1),1))<=1E-5
DTb=1E-5;
else if(Tb_cb(i,1)-Tb_cb((i-1),1))<=-1E-5
DTb=-1E-5;
else
DTb=(Tb_cb(i,1)-Tb_cb((i-1),1));
end
end
if i==1
drhoPb=1E-10;
drhoTb=1E-10;
else
if abs(rhob(i,1)-rhob(i-1,1))<=1E-20
drhoPb=1E-10;
drhoTb=1E-10;
else
drhoPb=(rhob(i,1)-rhob((i-1),1))/(Pb_cb(i,1)-Pb_cb((i-1),1));
drhoTb=(rhob(i,1)-rhob((i-1),1))/DTb;
end
end
A=((1/rhob(i,1))*(drhoTb/drhoPb))+(cpp/MWp); %modification is done
B=1/drhoPb;
if A>=1E20
A=1E20;
end
if A<=-1E20
A=-1E20;
end
if B>=1E20
B=1E20;
end
if B<=-1E20
B=-1E20;
end
hb=Hp/MWp; %Specific enthalpy of burned gases in J/kg
Ku=0.0243; %thermal conductivity of air
dQb=(Awb_cb*(((0.45*Ku*((80000)^0.75)*(Tb_cb(i,1)-Tw))/Bore)+((4.3*10^-9)*((Tb_cb(i,1))^4-Tw^4))))/W; % Heat loss in J/radian
f_Tb1(i,1)=((B/A)*((-dVb_cb/Vb_cb(i,1))+((-dQb+(dmb(i,1)*(hu-hb)))/(B*mb(i,1)))+(dmb(i,1)/mb(i,1))));
Tb_cb(i+1,1)=(f_Tb1(i,1)*dtheta)+Tb_cb(i,1); % Burned gas Temperature during combustion.
Pb_cb(i+1,1)=rhob(i+1,1)*(Ru/MWp)*Tb_cb(i+1,1); % Burned gas Pressure during combustion.
P_cb(i+1,1)=Pu_cb(i+1,1)+Pb_cb(i+1,1); % Total cylinder pressure.
xb=mb(i+1,1)/mivc;
rhoT=(xb*rhob(i+1,1))+((1-xb)*rhou(i+1,1));
MWT=((Np*MWp)+(Nr*MWr))/(Nr+Np);
T_cb=(P_cb(i+1,1)*(MWT))/(rhoT*Ru); % total cylinder temperature.
if i>(j1+j2)
thetarEOB=(i-1)*dtheta;
Nrcb=i;
break;
end
end
%---------------------------------------%
% Start of Expansion %
%---------------------------------------%
thetaexp=(SA+thetarEOB):dtheta:(540-EVO);
Nexp=floor(((540-EVO-SA-thetarEOB)/dtheta)+1);
V_exp=zeros(Nexp,1);
Aw=zeros(Nexp,1);
Tw=1300;
for i=1:Nexp
V_exp(i,1)=((Vd/(r-1))+((Vd/2)*(1+(L/Rad)-cosd(thetaexp(1,i))-((L/Rad)^2-(sind(thetaexp(1,i)))^2)^0.5))); % Incylinder volume (Cubic meter ) calculation
Aw(i,1)=(pi*0.5*Bore^2)+((pi*Bore*Rad)*((L/Rad)+1-cosd(thetaexp(1,i))+((L/Rad)^2-(sind(thetaexp(1,i)))^2)^0.5));
end
dVexp=gradient(V_exp,dtheta); % gradient of cylinder volume
P_exp=zeros(Nexp,1);
T_exp=zeros(Nexp,1);
rho_exp=zeros(Nexp,1);
Sp=zeros(Nexp,1);
A_exp=zeros(Nexp,1);
P_exp(1,1)=P_cb(Nrcb,1);
T_exp(1,1)=T_cb;
rho_exp(1,1)=rhoT;
A_exp(1,1)=A_cb(Nrcb,1);
cvpavg=0;
for i=1:Nexp-1 % expansion Loop
if T_exp(i,1)<1000
c1o2=0.03212936E2;c2o2=0.112748E-2;c3o2=-0.057561E-5;c4o2=0.131387E-8;c5o2=-0.0876855E-11;c6o2=-0.1005249E4;c7o2=0.060347E2;
c1n2=0.0329867E2;c2n2=0.140824E-2;c3n2=-0.0396322E-4;c4n2=0.0564151E-7;c5n2=-0.0244485E-10;c6n2=-0.1020899E4;c7n2=0.0395037E2;
c1co2=0.0227572E2;c2co2=0.0992207E-1;c3co2=-0.1040911E-4;c4co2=0.06866686E-7;c5co2=-0.0211728E-10;c6co2=-0.0483731E6;c7co2=0.1018848E2;
c1h2o=0.0338684E2;c2h2o=0.0347498E-1;c3h2o=-0.0635469E-4;c4h2o=0.0696858E-7;c5h2o=-0.0250658E-10;c6h2o=-0.0302081E6;c7h2o=0.0259023E2;
else
c1o2=0.036975E2;c2o2=0.0613519E-2;c3o2=-0.1258842E-6;c4o2=0.0177528E-9;c5o2=-0.1136435E-14;c6o2=-0.1233930E4;c7o2=0.0318916E2;
c1n2=0.0292664E2;c2n2=0.1487976E-2;c3n2=-0.0568476E-5;c4n2=0.10097038E-9;c5n2=-0.0675335E-13;c6n2=-0.0922797E4;c7n2=0.0598052E2;
c1co2=0.0445362E2;c2co2=0.0314016E-1;c3co2=-0.1278410E-5;c4co2=0.0239399E-8;c5co2=-0.1669033E-13;c6co2=-0.0489669E6;c7co2=-0.0955395E1;
c1h2o=0.02672145E2;c2h2o=0.0305629E-1;c3h2o=-0.0873026E-5;c4h2o=0.1200996E-9;c5h2o=-0.06391618E-13;c6h2o=-0.02989921E6;c7h2o=0.06862817E2;
end
cppco2=Ru*(c1co2+(c2co2*T_exp(i,1))+(c3co2*(T_exp(i,1))^2)+(c4co2*(T_exp(i,1))^3)+(c5co2*(T_exp(i,1))^4));
cpph2o=Ru*(c1h2o+(c2h2o*T_exp(i,1))+(c3h2o*(T_exp(i,1))^2)+(c4h2o*(T_exp(i,1))^3)+(c5h2o*(T_exp(i,1))^4));
cppo2=Ru*(c1o2+(c2o2*T_exp(i,1))+(c3o2*(T_exp(i,1))^2)+(c4o2*(T_exp(i,1))^3)+(c5o2*(T_exp(i,1))^4));
cppn2=Ru*(c1n2+(c2n2*T_exp(i,1))+(c3n2*(T_exp(i,1))^2)+(c4n2*(T_exp(i,1))^3)+(c5n2*(T_exp(i,1))^4));
cpp=((Xpco2*cppco2+Xph2o*cpph2o+Xpo2*cppo2+Xpn2*cppn2)); % specific heat of product in J/kmol-K
cvp=cpp-Ru;
Hpo2=Ru*T_exp(i,1)*(c1o2+((c2o2/2)*T_exp(i,1))+((c3o2/3)*(T_exp(i,1))^2)+((c4o2/4)*(T_exp(i,1))^3)+((c5o2/5)*(T_exp(i,1))^4)+(c6o2/T_exp(i,1)));
Hpn2=Ru*T_exp(i,1)*(c1n2+((c2n2/2)*T_exp(i,1))+((c3n2/3)*(T_exp(i,1))^2)+((c4n2/4)*(T_exp(i,1))^3)+((c5n2/5)*(T_exp(i,1))^4)+(c6n2/T_exp(i,1)));
Hpco2=Ru*T_exp(i,1)*(c1co2+((c2co2/2)*T_exp(i,1))+((c3co2/3)*(T_exp(i,1))^2)+((c4co2/4)*(T_exp(i,1))^3)+((c5co2/5)*(T_exp(i,1))^4)+(c6co2/T_exp(i,1)));
Hph2o=Ru*T_exp(i,1)*(c1h2o+((c2h2o/2)*T_exp(i,1))+((c3h2o/3)*(T_exp(i,1))^2)+((c4h2o/4)*(T_exp(i,1))^3)+((c5h2o/5)*(T_exp(i,1))^4)+(c6h2o/T_exp(i,1)));
Hp=((Xpco2*Hpco2+Xph2o*Hph2o+Xpo2*Hpo2+Xpn2*Hpn2)); % specific enthalpy of products in J/kmol
Spo2=(Ru*((c1o2*log(T_exp(i,1)))+(c2o2*T_exp(i,1))+((c3o2/2)*(T_exp(i,1))^2)+((c4o2/3)*(T_exp(i,1))^3)+((c5o2/4)*(T_exp(i,1))^4)+c7o2))-(Ru*log((Xpo2*P_exp(i,1))/P_atm));
Spn2=(Ru*((c1n2*log(T_exp(i,1)))+(c2n2*T_exp(i,1))+((c3n2/2)*(T_exp(i,1))^2)+((c4n2/3)*(T_exp(i,1))^3)+((c5n2/4)*(T_exp(i,1))^4)+c7n2))-(Ru*log((Xpn2*P_exp(i,1))/P_atm));
Spco2=(Ru*((c1co2*log(T_exp(i,1)))+(c2co2*T_exp(i,1))+((c3co2/2)*(T_exp(i,1))^2)+((c4co2/3)*(T_exp(i,1))^3)+((c5co2/4)*(T_exp(i,1))^4)+c7co2))-(Ru*log((Xpco2*P_exp(i,1))/P_atm));
Sph2o=(Ru*((c1h2o*log(T_exp(i,1)))+(c2h2o*T_exp(i,1))+((c3h2o/2)*(T_exp(i,1))^2)+((c4h2o/3)*(T_exp(i,1))^3)+((c5h2o/4)*(T_exp(i,1))^4)+c7h2o))-(Ru*log((Xph2o*P_exp(i,1))/P_atm));
Sp(i,1)=((Xpco2*Spco2+Xph2o*Sph2o+Xpo2*Spo2+Xpn2*Spn2)); % specific entropy of products in J/kmol
drhoP=(MWp/(Ru*T_exp(i,1)))*(cvp/(cvp+Ru));
drhoT=(rho_exp(i,1)*cvp)/(Ru*T_exp(i,1));
A=((1/rho_exp(i,1))*(drhoT/drhoP))+(cpp/MWp);
B=1/drhoP;
Ka=0.0243; %thermal conductivity of air
dQ_exp=(Aw(i,1)*(((0.45*Ka*((80000)^0.75)*(T_exp(i,1)-Tw))/Bore)+((4.3*10^-9)*((T_exp(i,1))^4-Tw^4))))/W; % Heat loss in J/degree
f_T1=((B/A)*((-dVexp(i,1)/V_exp(i,1))-(dQ_exp/(B*mivc))));
T_exp(i+1,1)=(f_T1*dtheta)+T_exp(i,1); % Temperature during expansion
P_exp(i+1,1)=(mivc*(Ru/MWp)*T_exp(i+1,1))/V_exp(i+1,1);
rho_exp(i+1,1)=(P_exp(i+1,1)*MWp)/(Ru*T_exp(i+1,1));
end
Spo2=(Ru*((c1o2*log(T_exp(Nexp,1)))+(c2o2*T_exp(Nexp,1))+((c3o2/2)*(T_exp(Nexp,1))^2)+((c4o2/3)*(T_exp(Nexp,1))^3)+((c5o2/4)*(T_exp(Nexp,1))^4)+c7o2))-(Ru*log((Xpo2*P_exp(Nexp,1))/P_atm));
Spn2=(Ru*((c1n2*log(T_exp(Nexp,1)))+(c2n2*T_exp(Nexp,1))+((c3n2/2)*(T_exp(Nexp,1))^2)+((c4n2/3)*(T_exp(Nexp,1))^3)+((c5n2/4)*(T_exp(Nexp,1))^4)+c7n2))-(Ru*log((Xpn2*P_exp(Nexp,1))/P_atm));
Spco2=(Ru*((c1co2*log(T_exp(Nexp,1)))+(c2co2*T_exp(Nexp,1))+((c3co2/2)*(T_exp(Nexp,1))^2)+((c4co2/3)*(T_exp(Nexp,1))^3)+((c5co2/4)*(T_exp(Nexp,1))^4)+c7co2))-(Ru*log((Xpco2*P_exp(Nexp,1))/P_atm));
Sph2o=(Ru*((c1h2o*log(T_exp(Nexp,1)))+(c2h2o*T_exp(Nexp,1))+((c3h2o/2)*(T_exp(Nexp,1))^2)+((c4h2o/3)*(T_exp(Nexp,1))^3)+((c5h2o/4)*(T_exp(Nexp,1))^4)+c7h2o))-(Ru*log((Xph2o*P_exp(Nexp,1))/P_atm));
Sp(Nexp,1)=((Xpco2*Spco2+Xph2o*Sph2o+Xpo2*Spo2+Xpn2*Spn2));
%--------------------------------------------------------%
theta=180+IVC:dtheta:540-EVO;
NT=floor(((540-EVO)-(180+IVC))/dtheta)+1;
V=zeros(NT,1);
for i=1:NT
V(i,1)=((Vd/(r-1))+((Vd/2)*(1+(L/Rad)-cosd(theta(1,i))-((L/Rad)^2-(sind(theta(1,i)))^2)^0.5))); % Incylinder volume (Cubic meter ) calculation
end
DV=gradient(V,dtheta);
P=zeros(NT,1);
P(1:Nc,1)=P_comp(:,1);
P(Nc+1:(Nc+Nrcb-1),1)=P_cb(2:Nrcb,1);
P((Nc+Nrcb):(Nrcb+Nc+Nexp-2),1)=P_exp(2:Nexp,1);
figure;plot(theta',P);grid on;
%---Modelling of Exhaust Process------------------%
Pr=(1.05/1.25)*P_atm;
Tr=T_exp(Nexp,1)/((P_exp(Nexp,1)/Pr)^(1/3));
%------------------Power and Torque Calculations------------%
IWD=0;
for i=1:NT
IWD=IWD+(P(i,1)*DV(i,1)*dtheta);
end
IMEP=(IWD/Vd)*10^-5
FMEP=(0.05*(N/1000)^2)+(0.15*(N/1000))+0.97
BMEP=IMEP-FMEP % in Bar
IP=IMEP*Vd*(N/120)*10^2 % in KW
ITRQ=(IP/W)*10^3 % in J
BP=BMEP*Vd*(N/120)*10^2
BTRQ=(BP/W)*10^3
efficicency=(BP*120*10^5)/(mf*LHV*N)
end