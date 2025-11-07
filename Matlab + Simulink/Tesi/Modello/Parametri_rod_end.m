dt = 1e-4; %passo di integrazione (s)
g=9.81; %Acceleazione di gravità (m/s^2)
m_in=0.0215; %Massa anello interno (kg)
m_out=0.1; %Massa anello esterno (kg)  La massa è comprensiva di rod...
Ixx_in=0.0000008; %momento di inerzia anello interno [kgm^2]
Iyy_in=0.0000008;
Izz_in=0.00000011;
Ixx_out=0.0000321; %Momento di inerzia anello esterno [kgm^2]
Iyy_out=0.0000035;
Izz_out=0.0000337;
R_in=0.019; %Raggio anello interno (m)
R_out=0.0195; %Raggio anello esterno (m)
E_in=2.1e11; %Modulo di Young (N/m^2)
E_out=2.1e11; %Modulo di Young (N/m^2)
v_in=0.3; %Coefficiente di Poisson
v_out=0.3; %Coefficiente di Poisson
sigma_k=(1-v_in^2)/(pi*E_in);
sigma_j=(1-v_out^2)/(pi*E_out);
k=4/(3*pi*(sigma_j+sigma_k))*((R_in*R_out)/(R_in+R_out))^0.5; %Rigidezza contatto
%k=5e10;
e=0.9; %Coefficiente di restituzione

v_rel=0:0.0001:5; %Velocità tangenziale relativa (m/s)
%v_st= 0.001;  %Velocità di soglia attrito statico (m/s)
v_st=0.03;
%v_dyn=0.005; %Velocità di soglia attrito dinamico (m/s)
v_dyn=0.04;
mu_st=0.11;
mu_v= 0.003;
mu_dyn= 0.06;
for i=1:length(v_rel)
if v_rel(i)>=0 && v_rel(i)<=v_st
    mu1(i)=Step(v_rel(i),-v_st,-(mu_st-mu_v*v_st),v_st,mu_st-mu_v*v_st)+mu_v*v_rel(i);
end
end

for i=1:length(v_rel)
if v_rel(i)>v_st
    mu2(i-length(mu1))=Step(v_rel(i),v_st,mu_st-mu_v*v_st,v_dyn,mu_dyn-mu_v*v_dyn)+mu_v*v_rel(i);
end
end

mu=[mu1,mu2];
plot(v_rel,mu)
% Kw=5.05e-4; %(mm^2/Nm) Dacapire meglio
% H=1; %(GPa) durezza brinell va bene come unità di misura?
% Var_aus=0; %Variabile ausialiaria per il calcolo dell'usura
% N=1/dt+1; %stop time, per avere il vettore Fc in 2D e non in 3D

% model = createpde;
% gm = importGeometry(model,"Anello interno.stl");
% center = mean(gm.Vertices);
% h = translate(gm,-center); %trasla le geometria g di una distanza s.
% %Così dovrei aver centarto la figura nel centro, c'è un modo per capirlo
% %dalla figura senza guardare punto per punto?
% pdegplot(model)
% mesh_default = generateMesh(model)
% figure
% pdemesh(mesh_default)
% %Bisogna trovare il modo di mettere lo zero della mesh al centro
% %dell'anello

Initial_condition=[0;0;0.000];   %Disallineamento iniziale dell'anello esterno
I=2.3;    %Momento di inerzia della manovella (da definire)(1e8)
m=61;  %Massa della manovella (da definire)(100)
c=55; %Smorzamento rotazionale della manovella 
L=0.37; %[m] Lunghezza a riposo della biella
t=0.38; %[m] Lunghezza del telaio
s=0.11; %[m] Lunghezza della manovella
eta=deg2rad(15); %Angolo di inclinazione del telaio rispetto all'orizzontale
Gamma0=acos((t^2+s^2-L^2)/(2*s*t));  %Angolo tra telaio e manovella all'istante iniziale
Alpha_p0=eta-acos((L^2+t^2-s^2)/(2*L*t)); %Angolo tra biella e l'orizzontale all'istante iniziale
Alpha0=acos((L^2+t^2-s^2)/(2*L*t));
Beta0=acos((L-t*cos(Alpha0))/s);
Beta01=acos((s^2+L^2-t^2)/(2*s*L));

Damping=1e3; %Smorzamento delle molle tridimensionali
Damping_rot=100; 
Stifness=1e8; %Rigidezza delle molle tridimensionali
Stiffness_rot=1e3; %Rigidezza torsionale molla

% Dati=importdata('Risultati_base.mat');
