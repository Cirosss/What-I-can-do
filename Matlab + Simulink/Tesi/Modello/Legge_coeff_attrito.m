%%Calcolo del coefficiente di attrito a partire dalla velocità relativa
clear all 
clc

v_rel=0:0.0001:5; %Velocità tangenziale relativa (m/s)
v_st= 0.001;  %Velocità di soglia attrito statico (m/s)
v_dyn=0.005; %Velocità di soglia attrito dinamico (m/s)
mu_st=0.11;
mu_v= 0.003;
mu_dyn= 0.06;
for i=1:length(v_rel)
if v_rel(i)>=0 & v_rel(i)<=v_st
    mu1(i)=Step(v_rel(i),-v_st,-(mu_st-mu_v*v_st),v_st,mu_st-mu_v*v_st)+mu_v*v_rel(i);
end
end

for i=1:length(v_rel)
if v_rel(i)>v_st
    mu2(i-length(mu1))=Step(v_rel(i),v_st,mu_st-mu_v*v_st,v_dyn,mu_dyn-mu_v*v_dyn)+mu_v*v_rel(i);
end
end

mu=[mu1,mu2];
plot(v_rel,mu)