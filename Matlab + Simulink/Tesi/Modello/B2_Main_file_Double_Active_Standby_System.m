%%% Author:  Andrea Dellacasa
%%%          Andrea Mornacchi
%%%          Polytechnic of Turin
%%% Subject: Aileron Model main file
%%%          Prognostic of Flight Controls
%%% Year: 2013
%%%--------------------------------%%%


% close all
% clear all
% clc

run('slblocks.m');

%% Hydraulic Data
 
Pr=600000;                      %[Pa] Reutn pressure
Ps=35000000;                    %[Pa] Supply pressure 35000000
T=40;                           %[°C] Oil temperature

run('D4_Fluid')                 %computation and storage fluid properties

%% Environmental Data
Patm=101325;                    %[Pa] Atmospheric pressure
load E5_rho_atm.mat;               %Load density profiler

%Wind
vmean_w=40;                     %[m/s]Mean speed wind
kg_w=0.80;                      %[-] Probability of gust [0-1]
load E3_Gust.mat;                  %Load gust profiler 

%Turbolence
TurbProb=3;                     %Probability of exceedance:
                                % 1) 2*10^-2
                                % 2) 10^-1
                                % 3) 10^-2 Light
                                % 4) 10^-3 Moderate
                                % 5) 10^-4
                                % 6) 10^-5 Severe
                                % 7) 10^-6
                                    
                 


%% Noise Switch 1-->ON, 0-->OFF
                       
Fae_on=1;                       %Aerodynamic force
Turb_on=1;                      %Turbulence
Del_on=1;                       %Electric noise

%% Fault switch 1-->ON, 0-->OFF

%% Data-air gap variation
air_gap_fault_time=10000;                   %Time at which the air gap start varying [hours]
first_cycle = 0;
flight_hours= 8;
k_airgap = 0.5/100*0.01;                    %Rate of change of the air gap height 
b_airgap = 0.3;                             %Rate of increase with usage
air_gap0_1 = 0;                             %value at begin of simulation
%airgap_size

%% Data-demagnetization of the torque motor
Temperature_threshold= 0;
m_thermal=0;                               % coefficient establishing the amount of demagnetization with the exposure time
R_th1 =0;                                  
tau_w1=0;
demag0_1 = 0;                              %jet pipe demagnetization at begin of each simulation
%demag_size

%% Data-jet pipe distorsion
enable_jetpipe_distorsion=0;
Gust_threshold=10000;
Gust_stiffness_jetpipe=1e8;
jetpipe_initialpos = 0;                    %jet pipe distortion at begin of each simulation
%jetpipe_size

%% Data-null shift
null_shift_sign = -1+((-1+2*rand(1,1))>=0)*2;
null_shift = -1+2*rand(1,1);

%% Data-crack growth block

% Spring size
b_fbspring = 2E-3;  % spring width [m]
h_fbspring = 1E-3;  % spring thickness [m]
a_svfb0_1 = 1E-3;   % initial size of the defect on side 1 [mm]
a_svfb0_2 = 0;      % initial size of the defect on side 2 [mm]

% Material properties
Rp02_svspring = 1014;   % Rp0.2 of the spring material [MPa] (assumed A666)

% Crack growth parameters
C_svfb = 1;          % stress intensity factor
A_svfb = 79.12.*1E-9.*sqrt(pi);% growth factor [between 1E-9 and 5E-9]
n_svfb = 3;                % growth power

% utilities
sv_crackON1 = 0; % switch for degradation on side 1 mettere a zero
sv_crackON2 = 0; % switch for degradation on side 2



%% Data-nozzle area variation
hyd_fluid_repl= 20;                        %frequency at which the fluid is replaced [hours]
a_jetpipe_area=0;                          %variation of the contamination level with time
k_c0_jetpipe_area=0*1e-1;                    %coefficient which depends on the general contamination level of the hydraulic fluid
jetpipe_passage_area=1;                    %jet pipe nozzle area [m2]
perc_passage_area0_1 = 0;                %percentage of area occlusion at begin of each simulation
%passage_area_size
%% UTAS_simulation parameters
rate_limiter_activation =0.000000;  % 1: rate limiter on ; 0: rate limiter off (ATP configuration)
inadvertant=0.000000;%              % Activation ou non de la simulation inadvertant driving (== +/-1 
deflection_rate=10000.000000;% NaN
%% Load physical data 
 
run ('C1_Actuator_1_Data')
run ('C2_Actuator_2_Data')
load('E2_Cmd_prog')

%% Simulation setting

dt_sim=1e-5;                    %[s] Step integration     1e-4;          
end_sim=Cmd(1,end);             %[s] End simulation time
