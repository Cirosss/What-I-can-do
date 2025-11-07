%Progetto Metodi e modelli numerici

%Sezione a 90 gradi circa al centro del soggiorno (vedi foto), qui
%l'interasse è di 300mm. Per cui il tubo si trova al centro di una porzione
%di pavimento con 15mm a dx, sx, sopra e sotto
clear all
close all
clc

h = 0.3;            %interasse di posa dei tubi
d_ext = 0.016;      %diametro esterno dei tubi
d_int = 0.013;      %diametro interno dei tubi

Amax=@(x,y)1e-5;

%Regions per il disegno dei tratti di tubo verticali
R1=regions.rect([h/2,h/2],[d_int,h]);
R2=regions.rect([-h/2,h/2],[d_int,h]);


%Regions che servono per disegnare il semicerchio dei tubi
Ca = regions.circle([0,0],[h/2+d_int/2]);
Cb = regions.circle([0,0],[h/2-d_int/2]);
R3 = regions.rect([0,h/2],[h,h]);

A=Ca-Cb-R3;
R=A+R1+R2;

S=regions.rect([0,h/4-0.05/2],[2*h,h/2+h+0.05]);       %Con questa definizione di S riesco a valutare il tubo come a 5cm dal muro
S =S-R;

%figure
%axis equal
%R.draw()
%S.draw('e')
%A.draw


%% Parametri
K=0.000832;  %m/s è la trasmissività da trovare

lambda_acqua=0.606; %W/(mK)
rho_acqua=1000;     %kg/m^3
cp_acqua=4186;      %J/(Kg*K)

lambda_pavimento=1.3;       %W/(mK) stimato da nostro excel
rho_pavimento=2000;          %kg/m^3
cp_pavimento=1000;            %da google

alpha_aria= 7.7;             %Preso da fisica tecnica
alpha_suolo= 1;          %da fisica tecnica
Taria= 273.15+20;           %nel soggiorno [K]
Tsuolo= 273.15+8;           %da UNI 11300  [K]

Pin=454.15;  %Pa
Pout=0;   %Pa
g=9.81;             %Accelerazione di gravità



%% problema in pressione
%Inserisco le BC di neumann su tutti i bordi perché viene più facile 
%mettere dirichlet poi solo sui due bordi necessari essendo che i raccordi
%hanno tanti bordi

for k=1:length(R.Borders)
    R.Borders(k).Bc(:)=boundaries.neumann(0);
end

R.Borders(1).Bc(37)=boundaries.dirichlet(Pout);
R.Borders(1).Bc(18)=boundaries.dirichlet(Pin);

R.mu=K;   %Termine diffusivo dell'equazione

figure
axis equal
R.draw('bc')

%costruzione della mesh
Me=mesh2D(R,Amax);

%figure
%axis equal
%Me.draw

%Posso risolvere l'equazione prenendo un esempio con solo dirichlet non
%omogeneo

f=@(x,y)0*x;                                %creo forzante nulla da passare alla funzione
A=FemAssembler.BuildStiffness(Me);
b=FemAssembler.BuildForce(Me,f);
bdir=FemAssembler.BuildDirichlet(Me);
b=b+bdir;

u=Me.copyToAllNodes(A\b); %Soluzione trovata
%figure
%Me.draw(u)

%Ora devo calcolare il gradiente della soluzione perchè sto cercando 
%la velocità
[vx,vy]=Me.gradient(u);
vx=-K*vx;   %Calcolo la velocità del fluido lungo x e y
vy=-K*vy;
figure
Me.quiver(vx,vy)
axis equal


%% equazione del calore
Tiniz=312.15;  %temperatura [K] alla quale ariva il flusso nella fase di carica

dt=3600; %[s]
tfin=60*60*24; %[s]  (24h)
time=(0:dt:tfin)';
T0=290;  %Temperatura a cui è il sistema all'istante iniziale

Tin=@(t)315.15;

R.mu=lambda_acqua;
S.mu=lambda_pavimento;

R.rho=rho_acqua*cp_acqua;
S.rho=rho_pavimento*cp_pavimento;

R.sigma = alpha_aria+alpha_suolo;
S.sigma = alpha_aria+alpha_suolo;

R.beta=@(x,y)[Me.interpolate(vx,[x,y]),Me.interpolate(vy,[x,y])];

%Metto le condizioni per avere scambio termico tra tubo e il pavimento
%circostante

R.Borders(1).Bc(:)=boundaries.none;

for k=1:length(S.Borders)
    S.Borders(k).Bc(:)=boundaries.none;
end

TempDomain=[R,S];
TempDomain(2).Borders(2).Bc(1)=boundaries.neumann(0);    
TempDomain(1).Borders(1).Bc(18)=boundaries.dirichlet(Tiniz);
TempDomain(1).Borders(1).Bc(37)=boundaries.neumann(0);
TempDomain(2).Borders(2).Bc(3)=boundaries.neumann(0);
TempDomain(2).Borders(2).Bc(22)=boundaries.neumann(0);
TempDomain(2).Borders(1).Bc(18)=boundaries.neumann(0);
TempDomain(2).Borders(2).Bc(2)=boundaries.periodic(@(x,y)[-x,y]);
TempDomain(2).Borders(2).Bc(23)=boundaries.periodic(@(x,y)[-x,y]);


figure
TempDomain.draw('bc')
%TempDomain.draw('bc')
%TempDomain.draw('e')

Me2=mesh2D(TempDomain,Amax);

%figure
%Me2.draw;



%% CALORE STAZIONARIO


g=@(x,y)alpha_aria*Taria+alpha_suolo*Tsuolo;


At=FemAssembler.BuildStiffness(Me2);
b2=FemAssembler.BuildDirichlet(Me2);
b1=FemAssembler.BuildForce(Me2,g);
bt=b1+b2;

w=Me2.copyToAllNodes(At\bt); %Soluzione trovata
min=min(w);
mean=mean(w);
figure
Me2.draw(w,'h')
shading interp
zlim([290 315.2])



%% Calore problema evolutivo

%figure
%TempDomain2.draw('bc')
%TempDomain2.draw('e')

M=FemAssembler.BuildMass(Me2);  %Costruisco solo la matrice di massa che manca nella parte evolutiva

% Uso Eulero implicito per la risoluzione
ndof=max(Me2.Nodes.Dof);
T=T0*ones(ndof,1);
ut=tfin/dt;  %Sono i passi temporali
minT=length(size(ut));
meanT=length(size(ut));
E=M+dt*At;   %Definisco la matrice di "Eulero"
figure
for i=1:ut

    T=E\(M*T+dt*bt);
    %minT(i)=max(T);
    %meanT(i)=mean(T);
    Me2.draw(Me2.copyToAllNodes(T));
    shading interp
    zlim([290 315.2])

    title(['t=',num2str(i*dt/3600),' h'])
    pause(0.5)
    hold off
end


