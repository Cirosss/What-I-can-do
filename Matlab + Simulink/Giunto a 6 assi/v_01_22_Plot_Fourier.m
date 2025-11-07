%% FILE PER AVERE FOURIER UNO SOPRA L'ALTRO
%README
%Inserire i nomi dei file da analizzare nella matrice_FileName,
%specificando la grandezza del nome dei file. Modficare matrice_f_acq con
%le diverse frequenze di acquisizione
%Scegliere il tipo di grafico modificando x,ylim

close all
clear all
clc


matrice_Filename(1,1:26)='Acquisizione_20s_2,4Hz.txt';
matrice_Filename(2,1:24)='Acquisizione_20s_5Hz.txt';
matrice_Filename(3,1:25)='Acquisizione_20s_10Hz.txt';
matrice_Filename(4,1:26)='Acquisizione_20s_100Hz.txt';
matrice_Filename(5,1:27)='Acquisizione_20s_1000Hz.txt';

%Definizione alternativa
%matrice_Filename=['Acquisizione_20s_2,4Hz.txt';'Acquisizione_20s_5Hz.txt';'Acquisizione_20s_10Hz.txt';'Acquisizione_20s_100Hz.txt';'Acquisizione_20s_1000Hz.txt'];

%NOTA: Modficare matrice_f_acq con le diverse frequenze di acquisizione
matrice_f_acq=[2.4; 5; 10; 100; 1000];

for i=1:length(matrice_f_acq)

    Folder = 'Acquisizioni';
    FileName = matrice_Filename(i,:);
    
    dati = importdata([Folder,filesep,FileName]);
    
    f_acq = matrice_f_acq(i); % Hz
    
    % vettore tempi
    tt = (0:1/f_acq*ones(size(dati(:,1))):1/f_acq*(length(dati(:,1))-1))';
    
    % segnali acquisiti (V)
    acq = dati(:,1:6);
    
    %OSS: le altre forze e momenti possiamo toglierli?
    % forze (N) 
    Fx = dati(:,7);
    Fy = dati(:,8);
    Fz = dati(:,9);
    
    % coppie (Nm)
    Mx = dati(:,10);
    My = dati(:,11);
    Mz = dati(:,12);
    
    %Trasformata di Fourier
    Y = fft(Fy);
    L=length(Fy);
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = f_acq*(0:(L/2))/L;
    
    figure (1)
    
    %colori 
    option=["r";"g";"b";"m";"k"];

    plot(f,P1,option(i))
    %loglog(f,P1,option(i))

    title('Single-Sided Amplitude Spectrum of X(t)')
    
    xlim([2 2.8])
    ylim([0 2.5])
    
    %loglog
    %xlim([0, 1000])
    %ylim([0.0001 10])
    
    xlabel('f (Hz)')
    ylabel('|P1(f)|')
    
    Leg(i)="f acq = " + num2str(f_acq);
    
    hold on
end

Leg(i+1)= "w_n = 2,4 Hz";
xline(2.4,'LineWidth',1.5,'Color',[0.9290 0.6940 0.1250])

legend(Leg) 
