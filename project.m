clc;
clear;
close all;
%1 
Ebsilon=[];
Itration=[];
[phi_new11 ,ebsilon1 ,itteration1]=WSOR(1,201,301);
delta_x=5/200;
delta_y=3/300;
x=0:delta_x:5;
y=0:delta_y:3;
contour(x,y,phi_new11);
surf(x,y,phi_new11);
Ebsilon(1)=ebsilon1;
Itration(1)=itteration1;
%2 
[phi_new22,ebsilon2 ,itteration2]=WSOR(1.2,201,301);
Ebsilon(2)=ebsilon2;
Itration(2)=itteration2;

%3
[phi_new33 ,ebsilon3 ,itteration3]=WSOR(1.4,201,301);
Ebsilon(3)=ebsilon3;
Itration(3)=itteration3;



%4
[phi_new44 ,ebsilon4 ,itteration4]=WSOR(1.6,201,301);
Ebsilon(4)=ebsilon4;
Itration(4)=itteration4;



%5
[phi_new55 ,ebsilon5 ,itteration5]=WSOR(1.8,201,301);
 Ebsilon(5)=ebsilon5;
 Itration(5)=itteration5; 
 
% plot of normale-iteration
 for i=1:5
     plot(1:Itration(i),Ebsilon(i))
     hold on;
 end
 % tabel of drfferent method
 result = table;
 result.name=['w=1','w=1.2','w=1.4','w=1.6','w=1.8'];
 result.valu=[Itration(1),Itration(1) ,Itration(3) ,Itration(4) ,Itration(5)];
 
 % stability
 
 
%[phi_exact] = exact_solotion(201,301);
%[phi_new1 ,ebsilon1 ,itteration1 ] = WLSOR_rowbyrow(1,201,301) ;
[phi_new2 ,ebsilon2 ,itteration2 ] = WLSOR_columnbycolumn(1,201,301) ;



function [phi_new ,ebsilon ,itteration] = WSOR(w,n_x,n_y) 

delta_x = 5/(n_x-1);
delta_y = 3/(n_y-1);
roh_y = (delta_x^2/(2*(delta_x^2 + delta_y^2)));
roh_x = (delta_y^2/(2*(delta_x^2 + delta_y^2)));
roh = ((delta_x^2 * delta_y^2)/(2*(delta_x^2 + delta_y^2)));
phi = zeros(n_y,n_x);
% here i define the boundry condition %
   phi(n_y,:)= 1;        

   for i = 1:n_x
    x = (i-1)*delta_x;
    phi(1,i) = x^2 - 10*x+2 ;
   end
   
% after defining boundry condition I want to solve thr equation 
phi_new = phi;
itteration = 1;
ebsilon =[1];
while ebsilon(itteration) > 10^(-12)
    itteration=itteration +1;
    for j = n_y-1:-1:2 
        for i = 2:1:n_x-1
            X = (i-1)*delta_x;
            Y = 3-(j-1)*delta_y;
            if i==1
                  phi_new(j,i) = phi(j,i)+ (w/2)*(roh_y*(phi(j+1,i)+phi_new(j-1,i))+ roh_x*(phi(j,i+1)-2*(-10*Y/3)*delta_x+phi(j,i+1))- roh*(((-34*(pi^2))/225)*cos((X*pi)/5)*sin((Y*pi)/3)+2*Y/3)-phi(j,i));  
            elseif i==n_x
                  phi_new(j,i) = phi(j,i) +(w/2)*(roh_y*(phi(j+1,i)+phi_new(j-1,i))+ roh_x*(phi_new(j,i-1)+phi_new(j,i-1))- roh*(((-34*(pi^2))/225)*cos((X*pi)/5)*sin((Y*pi)/3)+2*Y/3)-phi(j,i));
            else
                 phi_new(j,i) = phi(j,i) +(w/2)*(roh_y*(phi(j+1,i)+phi_new(j-1,i))+ roh_x*(phi_new(j,i-1)+phi(j,i+1))- roh*(((-34*(pi^2))/225)*cos((X*pi)/5)*sin((Y*pi)/3)+2*Y/3)-phi(j,i));
                
            end   

        end
    end
    delta =phi_new - phi ;
    ebsilon(itteration) = sqrt(sum(sum(delta.*delta))/(n_x*n_y));
    phi=phi_new;
end   

end

function [phi_exact] = exact_solotion(n_x,n_y)
delta_x = 5/(n_x-1);
delta_y = 3/(n_y-1);
phi_exact=zeros(n_y,n_x);
for j = n_y:-1:1
    for i = 1:1:n_x
        x =(i-1)*delta_x;
        Y = 3-(j-1)*delta_y;
        phi_exact(j,i)=cos(x*pi/5)*sin(Y*pi/3)+ (Y/3)*(x^2-10*x+1)+1;      
    end
end
end 


function [phi_new ,ebsilon ,itteration ] = WLSOR_rowbyrow(w,n_x,n_y) 
   % discritizing both dimension
delta_x = 5/(n_x-1);
delta_y = 3/(n_y-1);
   % calculation the coefficient of each equation
roh_y = (delta_x^2/(2*(delta_x^2 + delta_y^2)));
roh_x = (delta_y^2/(2*(delta_x^2 + delta_y^2)));
roh=((delta_x^2 * delta_y^2)/(2*(delta_x^2 + delta_y^2)));
   % meshing the domain 
phi = zeros(n_y,n_x);
% here i define the boundry condition %
   % down boundry condition
   phi(n_y,:)= 1;        
   % up boundry condition
   for i = 1:n_x
    x = (i-1)*delta_x;
    phi(1,i) = x^2 - 10*x+2 ;
   end
   
   % after defining boundry condition I want to solve thr equation 
phi_new = phi;
itteration = 1;
ebsilon = [1] ;
   % define two matrix which are requier for tomas algurithem
ans_mat = [];
coefficiant_mat=-eye(n_x);
vector=repmat(w*roh_x ,n_x-1,1);
mat1=diag(vector,1);
mat2 =diag(vector,-1);
coefficiant_mat = coefficiant_mat + mat1 +mat2;
coefficiant_mat(1,2)=2*w*roh_x;
coefficiant_mat(n_x,n_x-1)=2*w*roh_x;
% solving core
while ebsilon(itteration) > 10^(-12)
        itteration=itteration +1;
    for j = n_y-1:-1:2 
        for i = 1:1:n_x
            X = (i-1)*delta_x;
            Y = 3-(j-1)*delta_y;
            if i==1
                ans_mat(i) = (w-1)*phi(j,i) + w*(roh*((-34*pi^2/225)*cos(X*pi/5)*sin(Y*pi/3)+2*Y/3)-roh_y*(phi_new(j-1,i)+phi(j+1,i))+2*roh_x*(-10*Y/3)*delta_x);
            elseif i ==n_x
                ans_mat(i) = (w-1)*phi(j,i) + w*(roh*((-34*pi^2/225)*cos(X*pi/5)*sin(Y*pi/3)+2*Y/3)-roh_y*(phi_new(j-1,i)+phi(j+1,i))-2*roh_x*(-10*Y/3)*delta_x);
            else   
                ans_mat(i) = (w-1)*phi(j,i) + w*(roh*((-34*pi^2/225)*cos(X*pi/5)*sin(Y*pi/3)+2*Y/3)-roh_y*(phi_new(j-1,i)+phi(j+1,i)));
            end
        end
       ans=tomas(n_x+2,ans_mat,coefficiant_mat);
       for i=1:n_x
           phi_new(j,i)=ans(i);
       end
    end
        delta =phi_new - phi ;
       ebsilon(itteration) = sqrt(sum(sum(delta.*delta))/(n_x*n_y));
       phi=phi_new ; 
end
end

function [phi_new ,ebsilon ,itteration ] = WLSOR_columnbycolumn(w,n_x,n_y) 
   % discritizing both dimension
delta_x = 5/(n_x-1);
delta_y = 3/(n_y-1);
   % calculation the coefficient of each equation
roh_y = (delta_x^2/(2*(delta_x^2 + delta_y^2)));
roh_x = (delta_y^2/(2*(delta_x^2 + delta_y^2)));
roh=((delta_x^2 * delta_y^2)/(2*(delta_x^2 + delta_y^2)));
   % meshing the domain 
phi = zeros(n_y,n_x);
% here i define the boundry condition %
   % down boundry condition
   phi(n_y,:)= 1;        
   % up boundry condition
   for i = 1:n_x
    x = (i-1)*delta_x;
    phi(1,i) = x^2 - 10*x+2 ;
   end
   % right boundry condition
   for i = 1:n_y
       y = 3-(i-1)*delta_y;
       phi(i,1) = -10*y/3 ;
   end 
   % left boundry condition 
   phi(:,n_x)=0; 
   % after defining boundry condition I want to solve thr equation 
phi_new = phi;
itteration = 1;
ebsilon = [1] ;
   % define two matrix which are requier for tomas algurithem
ans_mat = [];
coefficiant_mat=-eye(n_y-2);
vector=repmat(w*roh_y ,n_y-3,1);
mat1=diag(vector,1);
mat2 =diag(vector,-1);
coefficiant_mat = coefficiant_mat + mat1 +mat2;
% solving core
while ebsilon(itteration) > 10^(-12)
        itteration=itteration +1;
    for i = 1:n_x
        for j = n_y-1:-1:2 
            X = (i-1)*delta_x;
            Y = 3-(j-1)*delta_y;
            if i==1
                ans_mat(j-1) = (w-1)*phi(j,i) + w*(roh*((-34*pi^2/225)*cos(X*pi/5)*sin(Y*pi/3)+2*Y/3)-roh_x*(phi_new(j,i+1)+2*(10*Y/3)*delta_x+phi(j,i+1)));
            elseif i==n_x
                ans_mat(j-1) = (w-1)*phi(j,i) + w*(roh*((-34*pi^2/225)*cos(X*pi/5)*sin(Y*pi/3)+2*Y/3)-roh_x*(phi_new(j,i-1)+phi(j,i-1)));
            else
                ans_mat(j-1) = (w-1)*phi(j,i) + w*(roh*((-34*pi^2/225)*cos(X*pi/5)*sin(Y*pi/3)+2*Y/3)-roh_x*(phi_new(j,i-1)+phi(j,i+1)));
            end
        end
       ans=tomas(n_y,ans_mat,coefficiant_mat);
       for j = n_y-1:-1:2 
           phi_new(j,i)=ans(j-1);
       end
    end
        delta =phi_new - phi ;
       ebsilon(itteration) = sqrt(sum(sum(delta.*delta))/(n_x*n_y));
       phi=phi_new;
end
end

function x=tomas(n_x,ans_mat,coefficiant_mat)
n=n_x-2;
A=coefficiant_mat;
d=ans_mat;
w=zeros(n,1);
q=zeros(n,1);
g=zeros(n,1);
w(1,1)=A(1,1);
g(1,1)=d(1,1)/w(1,1);
    for i=2:n
        q(i-1,1)=A(i-1,i)/w(i-1,1);
        w(i,1)=A(i,i)-A(i,i-1)*q(i-1,1);
        g(i,1)=(d(i,1)-A(i,i-1)*g(i-1,1))/w(i,1);
    end
x(n,1)=g(n,1);
    for i=n-1:-1:1
      x(i,1)=g(i,1)-q(i,1)*x(i+1,1);
    end
end