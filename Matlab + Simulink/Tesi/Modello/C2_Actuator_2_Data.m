%%% Author:  Andrea Dellacasa
%%%          Andrea Mornacchi
%%%          Polytechnic of Turin
%%% Subject: Aileron Model and Data
%%%          Prognostic of Flight Controls
%%% Year: 2013
%%%--------------------------------%%%

 
%%Servovalve Data
% Moog series 30

    imax_sv_act2=8e-3;                 %Input max current [A]

% Torque motor
    la_sv_act2=41e-3;                   %Distance between left and right pole [m]

    lp_sv_act2=28e-3;                   %Length of permanent magnet [m]
    Ap_sv_act2=504e-6;                  %Cross_sectional area of permanent magnet [m2]
    Vp_sv_act2=5.67e3;                       %Magnetomotive force of permanent magnet [A]
    Br_sv_act2=0.28;                         %Residual flux density of permanent magnet [T]
    B0_sv_act2=0.460;                        %Normalizing flux density [T]
    mup_sv_act2=1.1*4*pi/1e7;                %Permeability of magnet [N/A2]
    mua_sv_act2=4*pi/1e7;                    %Permeability of air [N/A2]

    Ag_sv_act2=18e-6;                        %Cross_sectional area of on air-gap [m2]
    l0_sv_act2=0.5e-3;                       %Normal value of the four air_gaps [m]
    H_sv_act2=0;                             %Air_gap height error [m]
    h_sv_act2=H_sv_act2/l0_sv_act2;                      
    W_sv_act2=0;                             %Air_gap right and left imbalance [m]
    w_sv_act2=W_sv_act2/l0_sv_act2;
    G_sv_act2=0;                             %Air_gap imbalance incline of armature [m]
    g_sv_act2=G_sv_act2/l0_sv_act2;
    gap_ms_act2=0;                           %gap between sping and spool
    
    n1_sv_act2=1140;                         %Number of coil turns
    n2_sv_act2=1140;                         %Number of coil turns

    k_sv_act2=mup_sv_act2*Vp_sv_act2/(2*B0_sv_act2*l0_sv_act2)-1;                        
    T0_sv_act2=la_sv_act2*Ag_sv_act2*Vp_sv_act2^2*mua_sv_act2/(4*(k_sv_act2+1)^2*l0_sv_act2^2);

%Hydraulic amplifer
    Gq_sv_act2=0.087;                        %Hydraulic gain [(m3/s)/m]
    Gp_sv_act2=(Ps-Pr)/3.05e-5;              %Pressure gain [Pa/m]
    dfl_sv_act2=5.0654e-5;                   %Equivalent nozzle diameter [m] 
    Afl_sv_act2=dfl_sv_act2^2*pi/4;          %Equivalent nozzle area [m2] 
    Pocl_A_act2=0;
    Pocl_B_act2=0;
    Pott_act2=0;
    F_act_act2=1;

%Flapper
    Ja_sv_act2=1.9575e-5;                    %Rotational mass of armature/flapper [Nm/(m/s2)]
    bf_sv_act2=0.08;                         %Damping on armature/flapper [Nm/(m/s)]
    Kf_sv_act2=511.6;                        %Stifness of armature/flapper [Nm/m]
    xf_sv_max_act2=3.05e-5;                  %Max flapper displacement [m]

%Feedback spring
    Kw_sv_act2=74;                           %Feedback wire stiffness [Nm/m]
    lf = 0.025;                         %Force leverage of feedback wire [m]
%Spool
    Asp_sv_act2=1.35e-5;                     %Spool and area [m2]
    msp_sv_act2=0.005;                       %Spool mass [Kg]
    bsp_sv_act2=0.25;                        %Spool damping [N/(m/s)]
    xspmax_sv_act2=0.001250;                 %Max spool displacement [m]
    wsp_sv_act2=0.009334348785340;           %Width of servovalve port [m]
    hr_sv_act2=6e-6;                         %Spool radial gap [m]

    Cdmax_sv_act2=0.7392;                    %Discharge coefficient max
    Cdmin_sv_act2=0.6428;                    %Discharge coefficient min
    Remax_sv_act2=5000;                      %Reynold max
    F_att_sv_act2=1;                         %Spool couombian friction [N] cambiata era 1
    rj_sv_act2=1e-6;                         %Rdius of curvature [m]
    load E4_R_sp_C                              %coef. curve c'
    
%Lap
    OV_S1_act2=0*2.344e-05;                    %Spool lap Supply-1 [m]  0.01*1e-3
    OV_R1_act2=0*2.344e-05;                    %Spool lap Return-1 [m] 
    OV_S2_act2=0*2.344e-05;                    %Spool lap Supply-2 [m]  0.01*1e-3
    OV_R2_act2=0*2.344e-05;                    %Spool lap Return-2 [m] 

%% Shut off bypass valve (most of the data have been assessed)
    xmax_sb_act2= 5.5e-3;                                    %valve stroke [m]
    A_ss_act2= pi*(4e-3)^2/4 ;                               %Spool end area [m2] diametro prima 4 mm
    A_sp_act2= pi*(2.5e-3)^2/4 ;                             %pressure area (piston) [m2]diametro prima 2.5 mm
    V0_ss_act2= (A_ss_act2-A_sp_act2)*(xmax_sb_act2+0.1e-3); %[m^3] Volume chambers shutoff bypass valve
    Vo_solpool_act2= A_ss_act2*(0.1e-3);                     %[m^3] Volume chambers shutoff valve spool dynamic
    K_solpool_act2= 13300;                                   % Spring rate spool shut off valve [N/m]

 %Holes
    d_hs_act2= 3e-3;
    N_hs_act2= 2;
    d_hr_act2= 1.5e-3;
    N_hr_act2= 1;
    load E1_Circ_segment                                        %Ellipsoid area of passage seen during opening of valve with circular holes

 %Lap
    OL_sb_S1_act2=0;                                         %Spool lap Supply-1 [m] 
    OL_sb_S2_act2=0;                                         %Spool lap Supply-2 [m]
    OL_sb_R1_act2=0;                                         %Spool lap Return-1 [m]
    OL_sb_R2_act2=0;                                         %Spool lap Return-2 [m]
    hr_ss_act2= 0.0045e-3;                                   %Spool radial gap [m]
    Lport_ss_act2= d_hs_act2;                                %Width of servovalve port [m]
    Lport_s2_act2= d_hr_act2;                                %Width of servovalve port [m]
    rj_ss_act2= 1e-6;                                        %Rdius of curvature [m]

% Shut off valve dynamic
    
    d_hrs_act2= 3e-3;
    Lport_s3_act2= d_hrs_act2;                               %Width of servovalve port [m]
%    V_junc _act2= V0_ss;
    A_ss_2_act2= pi*(9.5e-3)^2/4 ;
    A_sp2_act2= pi*(6e-3)^2/4 ;
    V_junc_act2=(A_ss_2_act2-A_sp2_act2)*(xmax_sb_act2+0.1e-3);
    M_shut_spool_act2= 75e-3;                 %spool mass [kg]
    F_shutSpool_act2= 50;                     % Preload [N]
    b_soff_spool_act2=30;                      % Viscous friction coefficient [N/(m/s)]
% solenoid
    V_al_act2= 28;                            % Alimentazione [V]
    Res_act2= 28;                             % resistenza [Ohm]
    Ind_act2= 0.7;                            % induttanza [H]
    F_h0_act2= 18.3;                          % Preload [N] 20.8
    k_h_act2= 17700;                         % Spring rate [N/m]
    xsol_max_act2= 0.3e-3;                    % Solenoid stroke [m]
    d_sol_act2= 1e-3;                         % Solenoid Spool diameter [m]
    A_sol_act2= pi*d_sol_act2^2/4;                 % Active area [m]
    M_sol_act2= 1e-3;                         % Spool mass [kg]
    Fsb_act2= 2.2;                            % Static friction [N] 0.8... 3.3 per static friction con zero increase factor
    b_sol_act2= 30;                           % Viscous friction coefficient [N/(m/s)]
    K_magn_act2= 31.5;                        % magnetic coefficient [N/A2]

 %Isolation valve
    A_isol_act2= pi*(5e-3^2)/4;               %Spool end area [m2]
    F_isol_act2= 30;                          %Preload [N]
    Vo_isol_act2= A_isol_act2*(0.1e-3);            %[m^3] Volume chambers isolation valve spool dynamic
    Fsb_isol_act2= 10;                        % Static friction [N]
    b_isol_act2= 3;                           % Viscous friction coefficient [N/(m/s)]
    xmax_isol_act2= 1.5e-3;                   % Max spool displacement [m] 3e-3
    k_isol_act2= 440;                         % Spring rate [N/m] 640
    M_isol_act2= 1e-3;                        %spool mass [kg]
    
%Accumulator
    A_acc_act2= pi*(5e-3^2)/4;                %area [m2]
    F_acc_act2= 5;                            %Preload [N]
    Vo_acc_act2= A_acc_act2*(0.1e-3);              %[m^3] Volume chambers isolation valve spool dynamic
    Fsb_acc_act2= 5;                          % Static friction [N] 2.2
    b_acc_act2= 30;                           % Viscous friction coefficient [N/(m/s)]
    xmax_acc_act2= 10e-3;                     % Max spool displacement [m]
    P_acc_act2= 2.5e6;                        % was 1.5e6
    k_acc_act2= 1000;                         % Spring rate [N/m]
    M_acc_act2= 5e-3;                         %spool mass [kg]
%% Data actuator
%Through rod

%Geometrical data
    d_la_act2=0.047570;                       %Rod diameter [m]
    D_la_act2=0.109950;                       %Cylinder diameter [m]
    A1_la_act2=5.710000e-003;                 %Active area 1 [m2]
    A2_la_act2=5.710000e-003;                 %Active area 2 [m2]
    L_la1_act2=0.059900;                      %Half cylinder stroke retraction [m]
    L_la2_act2=0.084800;                      %Half cylinder stroke extension [m]
    Vm_la1_act2=A1_la_act2*L_la1_act2*0.17;
    Vm_la2_act2=A2_la_act2*L_la2_act2*0.17;
    Mrod_la_act2=11.313490;                   %Mass rod [kg]
    Mcyl_la_act2=30.425550;                   %Mass cylinder [kg]
    
%Equilibrium condiction
    x0_la_act2=0;                            %Initial position piston [m]
    V01_la_act2=Vm_la1_act2+A1_la_act2*(L_la1_act2+x0_la_act2);    %Chamber1 volume0 [m3]
    V02_la_act2=Vm_la2_act2+A2_la_act2*(L_la2_act2-x0_la_act2);    %Chamber2 volume0 [m3]
      
    Pmean_act2=0.5*(Ps-Pr);
    DP_eq_act2=Pmean_act2*(A1_la_act2-A2_la_act2)/(A1_la_act2+A2_la_act2);
    P0_damped=2.5e6;
    
    
    Ksa_la_act2=5e7;                         %Stiffness cylinder mount [N/m]
    Csa_la_act2=8e3;                         %Damping cylinder mount [Ns/m]

    Fstall_la_act2=180e3;                    %Maximum force [N]

 %Friction
    run D2_Data_seal

    Kstdi_la_act2=1.5;                       %Ratio of staic to dinamic
    Gamma_la_act2=1.667e4;                   %Coefficient of viscous friction [Ns/m]

 %Leakage coefficient
    Kle11_act2=0;                            %External chamber 1 [m3/(sPa)]
    Kle12_act2=0;                            %External chamber1 [m3/(sPa^1/2)]

    Kle21_act2=0;                            %External chamber2 [m3/(sPa)]
    Kle22_act2=0;                            %External chamber2 [m3/(sPa^1/2)]

    Kli1_act2=0;                             %Internal[m3/(sPa)]
    Kli2_act2=0;                             %Internal[m3/(sPa^1/2)]
    
%% Data Load

    Kt_la_act2=Kt_la;                         %Cross Stifness [N/m]
    Ct_la_act2=Ct_la;                         %Cross Damping  [Ns/m]
    Cext1_la_act2=300;                       %External Damping 1 [Ns/m]
    Cext2_la_act2=700;                       %External damping 2 [Ns/m]
   
%%backlash
    b_max_act2= 0*0.0001;                          %total backlash [m]
    b_mplus_act2=b_max_act2/2;
    b_mminus_act2= - b_mplus_act2;
    phi_rel0_act2= 0;                       %unstretched spring lenght [m]
    bEps_act2= 1e-10;                       %minimum backlash [m]
  


%% Feedback 
    
%LVDT Demodulation filter
    w_fb_act2=310;                            %Filter natural freqency [rad/s]
    z_fb_act2=0.5;                            %Filter damping

%A/D
    fs_AD_act2=350;                           %Sample Frequency [Hz]
    ST_AD_act2= ceil(10000/fs_AD_act2)/10000;      %Remote Electronic Unit Sample Rate [s]
    Nbit_AD_act2= 14;                         %DAC bits
    QI_AD_act2= 2*(L_la1_act2+L_la2_act2)/2^Nbit_AD_act2;             %Quantization interval [V]
%% Control

    Kp_act2=8;                               %Proportional gain [A/m]
    Ki_act2=0.1;                             %Integrator gain [A/(m/s)]
    eps_c_act2=0.02e-3;                      %Integrator dead band [m]
    Li_c_act2=4e-3;                          %Integrator saturation [A]
    L_c_act2=imax_sv_act2;                        %Command saturation [A]
   