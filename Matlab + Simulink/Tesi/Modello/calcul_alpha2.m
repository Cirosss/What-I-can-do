function [Curr,Alpha]=calcul_alpha(imax,xmax,x0,x1,l1,l2)

Curr=[0 0.5 1 2 3 4 6 8]*1e-3;
X=Curr.*xmax/imax;
for i=1:length(Curr)    
    
% Loi Exacte
     if X(i)<(x0/2)        
         Alpha(i)=1/l1;  
     elseif X(i)<((l2*(x1-x0)+l1*(x1+x0))/2)
         Alpha(i)=(2*(l2*x0-l1*x1)+sqrt(abs(4*(l1*x1-l2*x0)^2-4*(l2-l1)*((l2-l1)*x0^2-2*(x1-x0)*X(i)))))/(2*(l2-l1)*X(i)); 
     else 
         Alpha(i)=(2*X(i)+(l2-l1)*(x1+x0))/(2*l2*X(i));  
     end

% Loi approchée


%     if X(i)<(x0/2)        
%         Alpha(i)=1/l1;  
%     else
%         Alpha(i)=(X(i)+x0*(l2-l1))/(l2*X(i));
%     end


end
        
Alpha=max(Alpha,1);