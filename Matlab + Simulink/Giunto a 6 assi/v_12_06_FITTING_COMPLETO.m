%MODELLO MASSA-MOLLA-SMORZATORE

close all
clear all
clc
%Dobbiamo mettere fourire di tutte le prove sullo stesso grafico così il
%picco rimane fino alla frequenza fi 5 Hz che è nyquist e dopo si perde
%l'andamento

%% ACQUISIZIONE DATI

% caricamento dati
Folder = 'Acquisizioni';
FileName = 'Acquisizione_20s_1000Hz_mk6_lunghezzaBraccio2,46.txt';

dati = importdata([Folder,filesep,FileName]);

% frequenza di acquisizione
f_acq = 1000; % Hz

% vettore tempi
tt = (0:1/f_acq*ones(size(dati(:,1))):1/f_acq*(length(dati(:,1))-1))';

% segnali acquisiti (V)
acq = dati(:,1:6);

% forze (N) 
Fx = dati(:,7);
Fy = dati(:,8);
Fz = dati(:,9);

% coppie (Nm)
Mx = dati(:,10);
My = dati(:,11);
Mz = dati(:,12);

%% ACQUISIZIONE DATI CON  MATRICE
%Alternativa all'acquisizione dati precedente, qui il file txt ha solo i
%valori di tensione e quindi deve passare attraverso la matrice di
%conversione ed è quello che fa la function che verrà citata

%caricamento direttorio
% Folder = 'Acquisizioni';
% FileName = 'b1_1000hz_dots.txt';
% 
% [Fx,Fy,Fz,Mx,My,Mz,tt]=Acq_con_matrice(Folder,FileName);

%% Parametri da controllare

%Dati misurati sperimentalmente
%m=0.62;
m=2.0417; % massa [kg]
%k=15.39^2*m;
k=484; %rigidezza molla [N/m]
g=9.81; %Accelerazione di gravità [m/s^2] 

f_acq=1000; %[Hz] Almeno 10 volte la frequenza naturale del sistema in modo da prender ein maniera molto precisa l'andamento

w_n=sqrt(k/m); %frequenza naturale 15.39 (rad/s)
T=2*pi/w_n; %Periodo di oscillazione
%OSS: non conoscendo zita utilizziamo w_n poichè per smorzamenti molto
%bassi il suo valore è molto simile a w_d

%Per come è scritto il codice si preferisce avere una forza positiva
if mean(Fy)<0
    Fy=-Fy;
end

%% PRE-ANALISI
%Vengono stimati alcuni parametri di prima ipotesi var_ hp da utilizzare come var_in per fminsearch

%A) Decremento logaritmico

%Filtro
Fy_f=movmean(Fy,max([round(f_acq/50),1]));

%Massimi
[pks,loc]=findpeaks(Fy_f,f_acq,"MinPeakDistance",0.9*T);

figure (3)
plot(tt, Fy, 'r',tt, Fy_f, 'g', loc, pks, 'bo')
legend("Dati sperimentali","Modello con media mobile","Picchi")

n=min([60, size(pks)]); %numero di picchi da considerare (prendo il minimo tra 15 e size pks)
delta=1/n*log(pks(1)/pks(n+1));
zeta_dec_log=delta/(sqrt(4*pi^2+delta^2));
beta_dec_log=2*zeta_dec_log*w_n*m;

%Trasformata di Fourier
Y = fft(Fy);
L=length(Fy);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = f_acq*(0:(L/2))/L;
figure (10)
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

%B) Ampiezza
%Minimi
[pks_min,loc_min]=findpeaks(-Fy_f,f_acq,"MinPeakDistance",0.9*T);
pks_min=-pks_min;
figure (3)
hold on
plot(loc_min, pks_min, 'ko')
legend("Minimi")

x0_hp=(pks(1)-pks_min(1))/k/2; %Doppio meno perchè pks_min è un valore negativo
m_hp=(pks(1)+pks_min(1))/g/2;  

%C) Fase
n=10; %numero di picchi da considerare
Ts_hp=(loc(n+1)-loc(1))/n; %Periodo (media primi 10 picchi) 

ws_hp=2*pi/Ts_hp;
ws_dec_log=w_n*sqrt(1-zeta_dec_log^2);
[max_idx, idx]=max(P1(10:end));
ws_ttf=f(9+idx)*2*pi;

phi_hp=loc(1)*ws_hp-pi/2;
if phi_hp<0
    phi_hp=2*pi+phi_hp; %+ perchè phi è negativo
end

%% FITTING PARAMETRI

%FORZA DEL VINCOLO

%NOTA: var=[m, beta, k, x0, phi]
%Equazione modello con variabile vettoriale
%Primo con mg seconda senza
F_sim=@(var, tt) var(1)*g+var(3)*var(4)*sin(sqrt(var(3)/var(1))*sqrt(1-(var(2)/(2*sqrt(var(3)/var(1))*var(1)))^2).*tt-var(5)).*exp(-(var(2)/(2*sqrt(var(3)/var(1))*var(1))*sqrt(var(3)/var(1)).*tt));
%F_sim=@(var, tt) var(3)*var(4)*sin(sqrt(var(3)/var(1))*sqrt(1-(var(2)/(2*sqrt(var(3)/var(1))*var(1)))^2).*tt-var(5)).*exp(-(var(2)/(2*sqrt(var(3)/var(1))*var(1))*sqrt(var(3)/var(1)).*tt));


%Valori primo tentativo di test
m_in=m*1.01;
beta_in=beta_dec_log;
k_in=k*0.98;
x0_in=x0_hp;
phi_in=phi_hp;

var_in=[m_in,beta_in,k_in,x0_in,phi_in]; %vettore variabili primo tentativo

%Funzione errore - fminsearch
err= @(var) norm(F_sim(var,tt)-Fy);
option= optimset('Display','iter','TolFun',1e-6,'TolX',1e-6,'MaxFunEvals',5000);
var_opt= fminsearch(err,var_in,option); %parametri di fitting
zeta=var_opt(2)/2/(((var_opt(3)/var_opt(1))^(1/2)));
w_n_opt=(var_opt(3)/var_opt(1))^(1/2);
disp("    Confronto var in - opt")
disp(["m    ", "beta", "k   ", "x0  ", "phi ","zeta", "w_n"])
disp([var_in zeta_dec_log w_n;var_opt zeta w_n_opt])
%Comparazione_in_exp_opt=[var_in;var_exp;var_opt]

figure (2)
plot(tt, F_sim(var_opt,tt))
hold on
plot(tt, Fy)
legend('Fitting','Dati sperimentali')

%% FOR VAR_OPT

for j=1:10
    

    %Valori primo tentativo di test
    m_in=m*(100+randn)/100;
    beta_in=beta_dec_log;
    k_in=k*(100-randn)/100;
    x0_in=x0_hp;
    phi_in=phi_hp;
    
    var_in_for(j,:)=[m_in,beta_in,k_in,x0_in,phi_in]; %vettore variabili primo tentativo
    
    %Funzione errore - fminsearch
    err= @(var) norm(F_sim(var,tt)-Fy);
    option= optimset('TolFun',1e-6,'TolX',1e-6);
    var_opt_for(j,:)= fminsearch(err,var_in,option); %parametri di fitting

end

zeta_for=mean(var_opt_for(:,2)/2./(((var_opt_for(:,3)./var_opt_for(:,1)).^(1/2))));
w_n_for=(((var_opt_for(:,3)./var_opt_for(:,1)).^(1/2)));
disp("    Confronto var in - opt CICLO FOR")
disp(["m    ", "beta","k   ", "x0  ", "phi ","zeta", "w_n"])
disp([mean(var_in_for) zeta_dec_log w_n;mean(var_opt_for) zeta_for mean(w_n_for)])

%% MOMENTO

h=(Mx./Fy); %- perchè hanno segni opposti
h_mean=mean(h);