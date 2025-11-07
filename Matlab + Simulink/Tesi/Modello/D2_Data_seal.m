%%% Author:  Andrea Mornacchi/ Andrea Dellacasa
%%%          Polytechnic of Turin
%%% Subject: Seal Data
%%%          Prognostic of Flight Controls
%%% Year: 2012
%%%--------------------------------%%%
%%% Version 2-25/10/2012
%%%--------------------------------%%%

%% Data Seal Actuator1

De_se=[53.80 53.80 53.80 97.65 53.80 53.80 53.80]';        %Seal external diameter [mm] 
Di_se=[47.57 47.57 47.57 59.88 47.57 47.57 47.57]';        %Seal internal diameter [mm]
Dm_se=(Di_se+De_se)/2;                                     %Seal meam diameter [mm]
Dd_se=[47.57 47.57 47.57 97.65 47.57 47.57 47.57]';        %Sliding seal diameter [mm]
W_se=[3.53 3.53 3.53 3.53 3.53 3.53 3.53]';                %Seal cross-section [mm]
Tw_se=[12.18 12.18 12.18 12.18 12.18 12.18 12.18]';        %Squeeze seal [%]
Hs_se=[70 70 70 70 70 70 70]';                             %Shore hardness [sh]

%% Estimate FH
%FH equation
%FH=0.175*pi*Dd*(-0.884+0.0206*Hs-0.0001Hs^2)*Tw

m=length(De_se);                                                                           %Matix dimension
FH=0;
for ii=1:m;
    FH=FH+0.175*pi*Dd_se(ii)*(-0.884+0.0206*Hs_se(ii)-0.0001*(Hs_se(ii))^2)*Tw_se(ii);     %Estimate FH [N]
end


%% Estimate FP constant part
%FP equation
%FP=7.82*10^(-2)*pi*Dm*W

m=3;                                                                                                                                    %Matix dimension
KFPi=zeros(m,1);                                                                          %Preallocation KFPi
for ii=1:m; 
    KFPi(ii,1)=7.82e-2*pi*Dm_se(ii+2)*W_se(ii+2);                                          %Estimate KFPi for each seal [N]
end

clear ('De_se','Di_se','Dm_se','De_se','Dd_se','W_se','Tw_se','Hs_se','ii','m')
