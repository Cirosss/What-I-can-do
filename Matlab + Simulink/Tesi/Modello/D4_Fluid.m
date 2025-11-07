%%% Author:  Andrea Mornacchi
%%%          Polytechnic of Turin
%%% Year: 2012
%%%
%%%Fluid: Phosphate ester type IV/V:
%%%       Skydrol500
%%%--------------------------------%%%
%%% Version 2-07/12 At the bottom release notes
%%%--------------------------------%%%


%% Expression of fluid properties (function of temperature)

%Density [kg/m3]:
    ro=(-0.85*T +1015.4)*(1+20/8000);  %Density ro(T) [kg/m3]
                      
%Kinematic viscosity [cS]: Log[Log(ni-cv)]=av-bvT
    ni=1e-6*10^(0.3829*log10(T +51)^4-1.8481*log10(T +51)^3+2.2079*log10(T +51)^2-1.1052*log10(T +51)+2.9544); %  [m^2/s] oil viscosity 
    
%Absolute viscosity [Pa*s]:
    mu=ro*ni;                            %Absolute viscosity mu(T) [Pa*s]
      
%Equivalent bulk modulus [Pa]: 
    beta=800000000;              
    
%%% Release note
%%% 
%%% Version 1 08/Apr/2012 Firt Relase
%%% Version 2 07/Dec/2012 fr function of oil temperature and pressure
