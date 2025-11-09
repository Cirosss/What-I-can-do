close all
clear all
clc

t=0.38; %[m] lunghezza del telaio
s=0.11; %[m] Lunghezza della manovella
l=0.37; %[m] Lunghezza a riposo della biella
delta_l=-0.06:0.0001:0.06;

L=l+delta_l;
alpha=acos((L.^2+t^2-s^2)./(2*L.*t));
beta=acos((L-t.*cos(alpha))./s);
gamma_0=pi-alpha(601)-beta(601);
gamma=pi-alpha-beta;
theta=gamma-gamma_0;

figure
plot(delta_l,theta*180/pi,'LineWidth',1)
title('Gear Ratio')
xlabel("\DeltaL [m]")
ylabel("\theta [°]")
grid on


I=2.3;    %Momento di inerzia della manovella (da definire)(1e8)
m=61;  %Massa della manovella (da definire)(100)
meq=zeros(1,length(theta));
for i=1:length(theta)
    meq(i)=(m+I/(s^2))*(delta_l(i)/theta(i));
end

mean(meq);