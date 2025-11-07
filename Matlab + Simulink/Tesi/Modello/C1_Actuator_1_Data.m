 %%% Author:  Andrea Dellacasa
%%%          Andrea Mornacchi
%%%          Polytechnic of Turin
%%% Subject: Aileron Model and Data
%%%          Prognostic of Flight Controls
%%% Year: 2013
%%%--------------------------------%%%

 
%%% Author:  Andrea Dellacasa
%%%          Andrea Mornacchi
%%%          Polytechnic of Turin
%%% Subject: Aileron Model and Data
%%%          Prognostic of Flight Controls
%%% Year: 2013
%%%--------------------------------%%%
 
 
%%Servovalve Data
% Moog series 30
 
    imax_sv=8e-3;                       %Input max current [A]
    Resistance = 1; % to be set
    Voltage_com = Resistance; % to be set
% Torque motor
    la_sv=41e-3;                        %Distance between left and right pole [m]
 
    lp_sv=28e-3;                        %Length of permanent magnet [m]
    Ap_sv=504e-6;                       %Cross_sectional area of permanent magnet [m2]
    Vp_sv=5.67e3;                       %Magnetomotive force of permanent magnet [A]
    Br_sv=0.28;                         %Residual flux density of permanent magnet [T]
    B0_sv=0.460;                        %Normalizing flux density [T]
    mup_sv=1.1*4*pi/1e7;                %Permeability of magnet [N/A2]
    mua_sv=4*pi/1e7;                    %Permeability of air [N/A2]
 
    Ag_sv=18e-6;                        %Cross_sectional area of on air-gap [m2]
    l0_sv=0.5e-3;                         %Normal value of the four air_gaps [m]
    H_sv=0;                             %Air_gap height error [m]
    h_sv=H_sv/l0_sv;                      
    W_sv=0;                             %Air_gap right and left imbalance [m]
    w_sv=W_sv/l0_sv;
    G_sv=0;                             %Air_gap imbalance incline of armature [m]
    g_sv=G_sv/l0_sv;
    gap_ms=0;                           %gap between sping and spool
    
    n1_sv=1140;                         %Number of coil turns
    n2_sv=1140;                         %Number of coil turns
 
    k_sv=mup_sv*Vp_sv/(2*B0_sv*l0_sv)-1;                        
    T0_sv=la_sv*Ag_sv*Vp_sv^2*mua_sv/(4*(k_sv+1)^2*l0_sv^2);
 
%Hydraulic amplifer
    Gq_sv=0.069;                         %Hydraulic gain [(m3/s)/m] 0.087 0.07 0.069
    Gp_sv=(Ps-Pr)/3.05e-5;              %Pressure gain [Pa/m]
    dfl_sv=5.0654e-5;                   %Equivalent nozzle diameter [m] 
    Afl_sv=dfl_sv^2*pi/4;               %Equivalent nozzle area [m2] 
    Pocl_A=0;
    Pocl_B=0;
    Pott=0;
%     F_act=1;
 
%Flapper
    Ja_sv=1.0575e-5;                    %Rotational mass of armature/flapper [Nm/(m/s2)] 1.9575e-5
    bf_sv=0.45;                         %Damping on armature/flapper [Nm/(m/s)] 0.08 0.4 0.45
    Kf_sv=301.6;                        %Stifness of armature/flapper [Nm/m] 511.6 291.6;
    xf_sv_max=3.05e-5;                  %Max flapper displacement [m]  3.05e-5
 
%Feedback spring
    Kw_sv=14.5;                           %Feedback wire stiffness [Nm/m] 74 24 13 15
    lf = 0.025;                         %Force leverage of feedback wire [m]
%Spool
    Asp_sv=1.35e-5;                     %Spool and area [m2] 1.35e-5
    msp_sv=0.005;                       %Spool mass [Kg] 0.005
    bsp_sv=0.25;                        %Spool damping [N/(m/s)] 0.25 0.0468
    xspmax_sv=0.001250;                 %Max spool displacement [m]
    wsp_sv=0.009334348785340;           %Width of servovalve port [m]
    hr_sv=6e-6;                         %Spool radial gap [m]
 
    Cdmax_sv=0.7392;                    %Discharge coefficient max
    Cdmin_sv=0.6428;                    %Discharge coefficient min
    Remax_sv=5000;                      %Reynold max
    F_att_sv=4;                         %Spool couombian friction [N] cambiata era 1
    rj_sv=1e-6;                         %Rdius of curvature [m]
    load E4_R_sp_C                              %coef. curve c'
%UTAS model additional parameters
    bk1=0.001500;                       % abscisse de début de changement de pente des fentes de servovalve (en courant)
    bk2=0.002500;                       % abscisse de fin de changement de pente des fentes de servovalve (en courant)
    x0=bk1*xspmax_sv/imax_sv;           % abscisse de début de changement de pente des fentes de servovalve (en course tiroir)
    x1=bk2*xspmax_sv/imax_sv;           % abscisse de fin de changement de pente des fentes de servovalve (en course tiroir)
    l1=0.500000;                        % largeur initiale en % de la largeur nominale équivalente Lar
    l2=(4* imax_sv - bk2 - bk1)/(4*imax_sv -2*bk2 -2*bk1);     % largeur finale en % de la largeur nominale équivalente Lar
%Lap
    OV_S1=0*2.344e-05;                    %Spool lap Supply-1 [m]  0.01*1e-3
    OV_R1=0*2.344e-05;                    %Spool lap Return-1 [m] 
    OV_S2=0*2.344e-05;                    %Spool lap Supply-2 [m]  0.01*1e-3
    OV_R2=0*2.344e-05;                    %Spool lap Return-2 [m] 

%% Shut off bypass valve (most of the data have been assessed)
    xmax_sb= 5.5e-3;                                    %valve stroke [m]
    A_ss= pi*(4e-3)^2/4 ;                               %Spool end area [m2] diametro prima 4 mm
    A_sp= pi*(2.5e-3)^2/4 ;                             %pressure area (piston) [m2]diametro prima 2.5 mm
    V0_ss= (A_ss-A_sp)*(xmax_sb+0.1e-3); %[m^3] Volume chambers shutoff bypass valve
    Vo_solpool= A_ss*(0.1e-3);                     %[m^3] Volume chambers shutoff valve spool dynamic
    K_solpool= 13300;                                   % Spring rate spool shut off valve [N/m]

 %Holes
    d_hs= 3e-3;
    N_hs= 2;
    d_hr= 1.5e-3;
    N_hr= 1;
    load E1_Circ_segment                                        %Ellipsoid area of passage seen during opening of valve with circular holes

 %Lap
    OL_sb_S1=0;                                         %Spool lap Supply-1 [m] 
    OL_sb_S2=0;                                         %Spool lap Supply-2 [m]
    OL_sb_R1=0;                                         %Spool lap Return-1 [m]
    OL_sb_R2=0;                                         %Spool lap Return-2 [m]
    hr_ss= 0.0045e-3;                                   %Spool radial gap [m]
    Lport_ss= d_hs;                                %Width of servovalve port [m]
    Lport_s2= d_hr;                                %Width of servovalve port [m]
    rj_ss= 1e-6;                                        %Rdius of curvature [m]

% Shut off valve dynamic
    
    d_hrs= 3e-3;
    Lport_s3= d_hrs;                               %Width of servovalve port [m]
%    V_junc = V0_ss;
    A_ss_2= pi*(9.5e-3)^2/4 ;
    A_sp2= pi*(6e-3)^2/4 ;
    V_junc=(A_ss_2-A_sp2)*(xmax_sb+0.1e-3);
    M_shut_spool= 75e-3;                 %spool mass [kg]
    F_shutSpool= 50;                     % Preload [N]
    b_soff_spool=30;                      % Viscous friction coefficient [N/(m/s)]
% solenoid
    V_al= 28;                            % Alimentazione [V]
    Res= 28;                             % resistenza [Ohm]
    Ind= 0.7;                            % induttanza [H]
    F_h0= 18.3;                          % Preload [N] 20.8
    k_h= 17700;                         % Spring rate [N/m]
    xsol_max= 0.3e-3;                    % Solenoid stroke [m]
    d_sol= 1e-3;                         % Solenoid Spool diameter [m]
    A_sol= pi*d_sol^2/4;                 % Active area [m]
    M_sol= 1e-3;                         % Spool mass [kg]
    Fsb= 2.2;                            % Static friction [N] 0.8... 3.3 per static friction con zero increase factor
    b_sol= 30;                           % Viscous friction coefficient [N/(m/s)]
    K_magn= 31.5;                        % magnetic coefficient [N/A2]

 %Isolation valve
    A_isol= pi*(5e-3^2)/4;               %Spool end area [m2]
    F_isol= 90;                          %Preload [N] 30
    Vo_isol= A_isol*(0.1e-3);            %[m^3] Volume chambers isolation valve spool dynamic
    Fsb_isol= 10;                        % Static friction [N]
    b_isol= 3;                           % Viscous friction coefficient [N/(m/s)]
    xmax_isol= 1.5e-3;                   % Max spool displacement [m] 3e-3
    k_isol= 440;                         % Spring rate [N/m] 640
    M_isol= 1e-3;                        %spool mass [kg]
    
%Accumulator
    A_acc= pi*(10e-2^2)/4;               %area [m2]
    F_acc= 5;                            %Preload [N]
    Vo_acc= A_acc*(0.01);                %[m^3] Volume chambers isolation valve spool dynamic
    Fsb_acc= 5;                          % Static friction [N] 2.2
    b_acc= 30;                           % Viscous friction coefficient [N/(m/s)]
    xmax_acc= 10e-2;                     % Max spool displacement [m]
    P_acc= 1.5e6;
    k_acc= 1000;                         % Spring rate [N/m]
    M_acc= 5e-3;                         %spool mass [kg]
%% Data actuator
%Through rod

%Geometrical data
    d_la=0.047570;                      %Rod diameter [m]
    D_la=0.09765;                       %Cylinder diameter [m]
    A1_la=5.710000e-003;                %Active area 1 [m2]
    A2_la=5.710000e-003;                %Active area 2 [m2]
    L_la1=0.059900;                     %Half cylinder stroke retraction [m]
    L_la2=0.084800;                     %Half cylinder stroke extension [m]
%     L_la1=L_la2;
    Vm_la1=A1_la*L_la1*0.17;            %Chamber1 volume  [m3]
    Vm_la2=A2_la*L_la2*0.17;
    Mrod_la=11.313490;                  %Mass rod [kg]
    Mcyl_la=30.425550;                  %Mass cylinder [kg]
    
%Equilibrium condiction
    x0_la=0;                            %Initial position piston [m]
%     V01_la=Vm_la1+A1_la*(L_la1+x0_la);    %Chamber1 volume0 [m3]
%     V02_la=Vm_la2+A2_la*(L_la2-x0_la);    %Chamber2 volume0 [m3]
    V01_la=(Vm_la1+A1_la*(L_la1+x0_la));    %Chamber1 volume0 [m3]
    V02_la=(Vm_la2+A2_la*(L_la2-x0_la));    %Chamber2 volume0 [m3]
    Pmean=0.5*(Ps-Pr);
    DP_eq=Pmean*(A1_la-A2_la)/(A1_la+A2_la);
    P01_la=Pmean-DP_eq;
    P02_la=Pmean+DP_eq;
    Pbupass0 = P_acc;
    
    Fstall_la=180e3;                    %Maximum force [N]

 %Friction
    run D2_Data_seal

    Kstdi_la=1.5;                       %Ratio of static to dynamic
    Gamma_la=1.667e4;                   %Coefficient of viscous friction [Ns/m]

 %Leakage coefficient
    Kle11=0;                            %External chamber 1 [m3/(sPa)]
    Kle12=0;                            %External chamber1 [m3/(sPa^1/2)]

    Kle21=0;                            %External chamber2 [m3/(sPa)]
    Kle22=0;                            %External chamber2 [m3/(sPa^1/2)]

    Kli1=0;                             %Internal[m3/(sPa)]
    Kli2=0;                             %Internal[m3/(sPa^1/2)]
    
%% Data Load
    J=177.000000;                       %Inertia surface[kg*m^2]
    Bl=0.155000;                        %Max leverage arm [m]
    Meq_la=J/(Bl^2);                    %Load Equivalent Mass [kg]
    mvar_e;                             %Elevator kinematics
    Kt_la=18720000;                     %Cross Stifness [N/m] 3.36 107N/m  era
    Ct_la=0.02*sqrt(Kt_la*Meq_la);      %Cross Damping  [Ns/m]
    Cext1_la=300;                       %External Damping 1 [Ns/m]
    Cext2_la=700;                       %External damping 2 [Ns/m]
   

    Ksa_la=12480000;                    %Stiffness cylinder mount [N/m] 5e7
    Csa_la=0.02*sqrt(Ksa_la*Meq_la);    %Damping cylinder mount [Ns/m] 8e3

    
%%backlash
    b_max= 0*0.0001;                    %total backlash [m]
    b_mplus=b_max/2;
    b_mminus= - b_mplus;
    phi_rel0= 0;                        %unstretched spring lenght [m]
    bEps= 1e-10;                        %minimum backlash [m]

%% Feedback 
%%UTAS AEU    

y0 =0.000000;%  initial conditions[m]

%% AEU
%AEU-signal processing
AEU_sampling_frequency =1/0.01;%1/0.01; 1/0.003%  [Hz] 
to_e=0.01;%1/0.02; 0.003;%  [Hz] 
ST_dsp = 1/AEU_sampling_frequency    ;%   [s]
N_bit =12.000000;%  [NA] Stable bits 
%AEU-bus delays
delay_REU2_ACE2 =2.000000;%  [N° of samples]
delay_REU1_ACE1 =2.000000;%  [N° of samples]
AEU_input_commutation_time =8.000000;%  [N° of samples]
AEU_output_commutation_time =5.000000;%  [N° of samples]
%AEU_Rid LVDT demodulation filter
w_fr = 2*pi*160 ;%  [rad/s]
ksi_fr =0.600000;%  [NA]
NUM_fr = [1] ;%  
DEN_fr = [(1/w_fr)^2 2*ksi_fr*(1/w_fr) 1] ;%  
NF_fr = tf(NUM_fr,DEN_fr);% NaN
NF_D_fr = c2d(NF_fr,1/AEU_sampling_frequency) ;%  Discretize Notch filter
[num_d_fr,den_d_fr] = tfdata(NF_D_fr,'v') ;%  Digital Notch filter
% cu=filt(num_d_fr,den_d_fr)
%AEU_DP-LVDT demodulation filter
w_Dp = 2*pi*125 ;%  [Hz]
NUM_Dp = [1] ;%  
DEN_Dp = [(1/w_Dp) 1] ;%  
NF_Dp = tf(NUM_Dp,DEN_Dp);% NaN
NF_D_Dp = c2d(NF_Dp,1/AEU_sampling_frequency) ;%  Discretize Notch filter
[num_d_Dp,den_d_Dp] = tfdata(NF_D_Dp,'v') ;%  Digital Notch filter
%AEU - signal range
FSrod_range =0.084800;%  [m]
FSehsv_range =0.001250;%  [m]
FSDP_range =44200000.000000;%  [Pa]
%AEU - signal quantification
QI = FSrod_range/2^N_bit    ;%  [m]
QI_p = FSDP_range/2^N_bit  ;%  [Pa]
QI_cmd_ACE = FSrod_range/2^N_bit  ;%  [m]
%AEU- controller
% position_gain =0.530000;%  [A/m]
AEU_current_saturation =0.008000;%  [A]
i0=0.000000;%  % Courant de commande de départ
Offset=0.000000;%    % Offset sur la consigne de position sinusoidale
Ka=0.500000;%     % gain du calculateur S/C (boucle de position) 0.530000
[Curr,Alpha]=calcul_alpha2(imax_sv,xspmax_sv,x0,x1,l1,l2);% NaN
rate_pilot='inf';%  % limitation de la consigne pilote
w0=57.800000;% % Paramètres linéaires filtre trou calé sur 9.2Hz (Rm=3.12e7N/m)
w1= 1.2*w0;% % Paramètres linéaires filtre trou calé sur 9.2Hz (Rm=3.12e7N/m)
z1=0.038800;% % Paramètres linéaires filtre trou calé sur 9.2Hz (Rm=3.12e7N/m)
z2=0.300000;% % Paramètres linéaires filtre trou calé sur 9.2Hz (Rm=3.12e7N/m)
type_c2d='matched';%  %numérique
[numzHA denzHA]=tfdata(c2d(tf([1/w0^2 2*z1/w0 1],[1/w1^2 2*z2/w1 1]),ST_dsp,type_c2d),'v');%  %numérique
Tinteg=10.000000;%   % Constante de temps de l'intégrateur de position
kep3=0.000000;%  % prise en compte intégrateur de compensation (=1) ou non (=0)
to_e=0.01;%0.010000; 0.003;%  % Echantillonnage des calculateurs primaires
% to_em=0.040000;%     % Echantillonnage de la consigne pilote
f1=0.005000;%     % retard émission
f11=0.004750;%     % retard reception
cu = tf(numzHA,denzHA,0.1,'variable','z^-1');
%% MODEL SIMULATION 
IS = 0.01/AEU_sampling_frequency  ;%  [s]
NGM =45.000000;% NaN
NPM =35.000000;% NaN

amp_i=0.000000;% % amplitude du courant de perturbation
freq=0.000000;% % fréquence du courant de perturbation
i0=0.000000;% NaN
Offset=0.000000;% NaNclo
y0=Offset;% NaN
type_signal=2.000000;% NaN
pas=0.000100;% NaN
pas_sauve=pas;% NaN



%%POLITO controller
%LVDT Demodulation filter
    w_fb=310;                            %Filter natural freqency [rad/s]
    z_fb=0.5;                            %Filter damping

%A/D
    fs_AD=350;                           %Sample Frequency [Hz]
    ST_AD= ceil(10000/fs_AD)/10000; %Remote Electronic Unit Sample Rate [s]
    Nbit_AD= 14;                         %DAC bits
    QI_AD= 2*(L_la1+L_la2)/2^Nbit_AD;   %Quantization interval [V]
%% Control

    Kp=8;                               %Proportional gain [A/m]
    Ki=0.1;                             %Integrator gain [A/(m/s)]
    eps_c=0.02e-3;                      %Integrator dead band [m]
    Li_c=4e-3;                          %Integrator saturation [A]
    L_c=imax_sv;                   %Command saturation [A]
   
    
    
    isvdata=[-0.00125000000000	-0.00001116666255
-0.00109375000000	-0.00000954044956
-0.00093750000000	-0.00000791423656
-0.00078125000000	-0.00000650485197
-0.00062500000000	-0.00000509546738
-0.00046875000000	-0.00000314401179
-0.00039062500000	-0.00000208155263
-0.00031250000000	-0.00000143106743
-0.00023437500000	-0.00000095404496
-0.00015625000000	-0.00000056375384
-0.00007812500000	-0.00000030355976
-0.00003125000000	-0.00000008239479
-0.00001562500000	-0.00000002601941
0.00000156250000	0.00000000000000
0.00001562500000	0.00000002168284
0.00003125000000	0.00000008673136
0.00007812500000	0.00000030355976
0.00015625000000	0.00000060711952
0.00023437500000	0.00000095404496
0.00031250000000	0.00000143106743
0.00039062500000	0.00000212491831
0.00046875000000	0.00000303559759
0.00062500000000	0.00000498705318
0.00078125000000	0.00000683009457
0.00093750000000	0.00000823947916
0.00109375000000	0.00001019093475
0.00125000000000	0.00001203397615];