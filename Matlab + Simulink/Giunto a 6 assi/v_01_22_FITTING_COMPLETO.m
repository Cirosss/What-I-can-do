%FITTING GIUNTO A 6 ASSI
%Il codice, partendo dai dati acquisiti nel programma labview, permette di
%identificare i parametri che meglio fittano il modello teorico
%massa-molla-smorzatore

close all
clear all
clc

%% ACQUISIZIONE DATI

%Il file testo generato prevede le prime 6 colonne dedicate ai segnali
%acquisiti in V, dalla 7 alla 9 alle forze in N a seguito della
%trasformazione della matrice (presente all'interno del codice labview),
%dalla 10 alla 12 i momenti riferiti al sistema di riferimento scelto per
%l'intero banco (vedere schema relazione)

% caricamento dati
Folder = 'Acquisizioni';
FileName = 'Acquisizione_20s_1000Hz.txt';

dati = importdata([Folder,filesep,FileName]);

% frequenza di acquisizione
%NOTA: impostare la corretta frequenza di acquisizione
f_acq = 1000; % Hz

% vettore tempi (essendo nota la frequenza di acquisizione è possibile
% generare il vettore dei tempi)
tt = (0:1/f_acq*ones(size(dati(:,1))):1/f_acq*(length(dati(:,1))-1))';

% segnali acquisiti (V)
acq = dati(:,1:6);

% forze (N) 
Fx = dati(:,7);
Fy = dati(:,8); %è la forza utilizzata per il fitting
Fz = dati(:,9);

% coppie (Nm)
Mx = dati(:,10);
My = dati(:,11);
Mz = dati(:,12);

%Per come è scritto il codice si preferisce avere una forza positiva
if mean(Fy)<0
    Fy=-Fy;
end

%% PRE-ANALISI
%Vengono stimati alcuni parametri di prima ipotesi var_ hp da utilizzare come var_in per fminsearch

%0) Dati misurati sperimentalmente
%questi valori possono essere anche solo una stima dei valori reali (servono per ipotizzare un periodo di oscillazione utile nelle successive pre-analisi) 
m=2.0417; % massa [kg]
k=484; %rigidezza molla [N/m]
g=9.81; %Accelerazione di gravità [m/s^2]

w_n=sqrt(k/m); %frequenza naturale 15.39 (rad/s)/ 2.4 Hz
T=2*pi/w_n; %Periodo di oscillazione (stima di massima)
%OSS: non conoscendo ancora zita utilizziamo w_n poichè per smorzamenti molto
%bassi il suo valore è molto simile a w_d


%A) Decremento logaritmico

%Filtro (media mobile)
Fy_f=movmean(Fy,max([round(f_acq/50),1]));

%Identificazione dei massimi
%OSS: si utilizza la funzione findpeaks che necessita del pacchetto matlab
%dedicato. La funzione richiede il segnale, la frequenza di acquisizione
%(permette di restituire loc in tempi invece che in "indici") e viene
%specificata una distanza minima attraverso il periodo T ricavato nel
%punto 0
[pks,loc]=findpeaks(Fy_f,f_acq,"MinPeakDistance",0.9*T);

figure (3)
plot(tt, Fy, 'r',tt, Fy_f, 'g', loc, pks, 'bo')
xlabel('Tempo [s]')
ylabel('F_y [N]')
legend("Dati sperimentali","Modello con media mobile","Picchi")

n=min([30, length(pks)]); %numero di picchi da considerare nel decremento logaritmico
delta=1/n*log(pks(1)/pks(n+1));
zeta_dec_log=delta/(sqrt(4*pi^2+delta^2));
beta_dec_log=2*zeta_dec_log*w_n*m;

%Trasformata di Fourier
%OSS: viene effettuata un'analisi della trasformata di fourier usata poi
%anche per confermare la frequenza di oscillazione

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
%Identificazione dei minimi (stesse condsiderazioni della ricerca dei
%massimi)
[pks_min,loc_min]=findpeaks(-Fy_f,f_acq,"MinPeakDistance",0.9*T);
pks_min=-pks_min;

figure (3)
hold on
plot(loc_min, pks_min, 'ko')
legend("Minimi")

x0_hp=(pks(1)-pks_min(1))/k/2; 
m_hp=(pks(1)+pks_min(1))/g/2; 
%è possibile anche ricavare un valore di primo tentativo per la massa (può essere eventualmente sostituita per la massa di prima ipotesi)

%C) Fase

%Identificazione frequenza di oscillazione (vengono utilizzati diversi
%metodi alternativi)
n=10; %numero di picchi da considerare
Ts_hp=(loc(n+1)-loc(1))/n; %Periodo (media primi 10 picchi) 

%OSS: sono stati implementati 3 metodi alternativi
ws_hp=2*pi/Ts_hp;
ws_dec_log=w_n*sqrt(1-zeta_dec_log^2);
[max_idx, idx]=max(P1(10:end));
ws_ttf=f(9+idx)*2*pi; %OSS: massimo della trasformata di Fourier

phi_hp=loc(1)*ws_hp-pi/2;
if phi_hp<0
    phi_hp=2*pi+phi_hp; %+ perchè phi è negativo
end
%per come è stato impostato il modello è corretto avere una fase positiva

%% FITTING PARAMETRI

%FORZA DEL VINCOLO

%Equazione modello con variabile vettoriale (è richiesto dalla funzione
%fminsearch)
%NOTA: var=[m, beta, k, x0, phi]

F_sim=@(var, tt) var(1)*g+var(3)*var(4)*sin(sqrt(var(3)/var(1))*sqrt(1-(var(2)/(2*sqrt(var(3)/var(1))*var(1)))^2).*tt-var(5)).*exp(-(var(2)/(2*sqrt(var(3)/var(1))*var(1))*sqrt(var(3)/var(1)).*tt));

%Valori primo tentativo di test
%OSS:vengono riportati i valori di primo tentativo (m e k vengono
%leggermente variate dagli eventuali valori attesi)
m_in=m*1.01;
beta_in=beta_dec_log;
k_in=k*0.98;
x0_in=x0_hp;
phi_in=phi_hp;

var_in=[m_in,beta_in,k_in,x0_in,phi_in]; %vettore variabili primo tentativo

%Funzione errore - fminsearch
err= @(var) norm(F_sim(var,tt)-Fy);
option= optimset('Display','iter','TolFun',1e-6,'TolX',1e-6,'MaxFunEvals',5000,'MaxIter',5000);
var_opt= fminsearch(err,var_in,option); %parametri di fitting (optimisation)

zeta_opt=var_opt(2)/2/(((var_opt(3)/var_opt(1))^(1/2)));
w_n_opt=(var_opt(3)/var_opt(1))^(1/2);

disp("    Confronto var in - opt")
disp(["m    ", "beta", "k   ", "x0  ", "phi ","zeta", "w_n"])
disp([var_in zeta_dec_log w_n;var_opt zeta_opt w_n_opt])

figure (2)
plot(tt, F_sim(var_opt,tt))
hold on
plot(tt, Fy)
legend('Fitting','Dati sperimentali')
xlabel('Tempo [s]')
ylabel('F_y [N]')

%% FOR VAR_OPT
% Questa section permette di effettuare j prove di fitting a fronte di
% piccole variazioni dei parametri di primo tentaivo (il procedimento utilizzato è quella della sezione precedente).
% (Nei dati sperimentali 22-23 si ottengono sempre gli stessi valori, ma in generale potrebbe non essere garantita sempre la stessa convergenza) 

for j=1:10
    disp(['Iter ',num2str(j),'/',num2str(10)])

    %Valori primo tentativo di test
    m_in=m*(100+randn)/100;
    beta_in=beta_dec_log;
    k_in=k*(100-randn)/100;
    x0_in=x0_hp;
    phi_in=phi_hp;
    
    var_in_for(j,:)=[m_in,beta_in,k_in,x0_in,phi_in]; %vettore variabili primo tentativo
    
    %Funzione errore - fminsearch
    err= @(var) norm(F_sim(var,tt)-Fy);
    option= optimset('TolFun',1e-6,'TolX',1e-6,'MaxIter',5000,'MaxFunEvals',5000);
    var_opt_for(j,:)= fminsearch(err,var_in,option); %parametri di fitting

end

zeta_for=mean(var_opt_for(:,2)/2./(((var_opt_for(:,3)./var_opt_for(:,1)).^(1/2))));
w_n_for=(((var_opt_for(:,3)./var_opt_for(:,1)).^(1/2)));
disp("    Confronto var in - opt CICLO FOR")
disp(["m    ", "beta","k   ", "x0  ", "phi ","zeta", "w_n"])
disp([mean(var_in_for) zeta_dec_log w_n;mean(var_opt_for) zeta_for mean(w_n_for)])

%% MOMENTO
% il braccio viene identificato come la media del rapporto tra il momento e
% la forza durente tutta la prova

h=(Mx./Fy); 
h_mean=mean(h);