clear
clc
close all

B=3.2;
H=0.9;
W=0.85;
g=9.81; %m/s^2
x=[2*W -W/2 0; 2*W W/2 0;2*W 0 H;0 0 H;0 -B H;0 B H; W 0 H]; %m
E_barra=75000*10^6; %Pa
Dens_barra=3350; %kg/m^3
D1=0.018; %m
d1=0.0075; %m
A_barra=pi*((D1/2)^2-(d1/2)^2); %m^2
In_barra=pi/4*((D1/2)^4-(d1/2)^4); %m^4
E_cable=147000*10^6; %Pa
Dens_cable=950; %kg/m^3
D2=0.003; %m
A_cable=pi*(D2/2)^2; %m^2
M=150; %kg
syms T;
D=T;


Tn=[1 5; 1 2; 2 3; 1 3; 3 5; 5 7; 4 5; 1 4; 1 7; 2 4; 2 7; 3 7; 4 7; 4 6; 6 7; 3 6; 2 6];
alpha=23*10^-6;
delta_T=0;
%F0=A*E*alpha*delta_T;
[numRows,numCols] = size(Tn);
Td=rand(numRows,numCols*3);
[numRowsx,numColsx] = size(x);
for i=1:numRows
    for j=1:numCols
    Td(i,3*j-2)=Tn(i,j)*3-2;   
    Td(i,3*j-1)=Tn(i,j)*3-1;
  
    Td(i,3*j)=Tn(i,j)*3;    
        
    end
    
end
K_G=zeros(numRowsx*numColsx);

F=zeros(max(max(Td)),1);
l_e=zeros(numRows,1);
m=zeros(numRows,1);
Momento_y=0;
for e=1:numRows
    x_1e=x(Tn(e,1),1);
    y_1e=x(Tn(e,1),2);
    x_2e=x(Tn(e,2),1);
    y_2e=x(Tn(e,2),2);
    z_1e=x(Tn(e,1),3);
    z_2e=x(Tn(e,2),3);
    l_e(e,:)=sqrt((x_2e-x_1e)^2+(y_2e-y_1e)^2+(z_2e-z_1e)^2);
    Ret=1/l_e(e)*[x_2e-x_1e y_2e-y_1e z_2e-z_1e 0 0 0; 0 0 0 x_2e-x_1e y_2e-y_1e z_2e-z_1e];
    if e==1 || e==8 || e==9 || e==10 || e==11 || e==17  
        K_e_prima=((A_cable*E_cable)/l_e(e))*[1 -1;-1 1];
        m(e,:)=Dens_cable*A_cable*l_e(e);
    else
        K_e_prima=((A_barra*E_barra)/l_e(e))*[1 -1;-1 1];
        m(e,:)=Dens_barra*A_barra*l_e(e);
    end
    K_e=Ret.'*K_e_prima*Ret;
    [numRows_ke,numCols_ke] = size(K_e);
    for i=1:numRows_ke
        for j=1:numCols_ke
            posRowsTd=e;
            posRowsK=Td(e,i);
            posColsK=Td(e,j);
            K_G(posRowsK,posColsK)=K_G(posRowsK,posColsK)+K_e(i,j);
        end   
    end
    Fm=[0; 0; -m(e)*g/2; 0; 0; -m(e)*g/2];
    for i=1:6
       F(Td(e,i))=F(Td(e,i))+Fm(i); 
    end
    
    
end
L=-1*(sum(F)-M*g);
F_ext=[T/2;0;-M*g/2;T/2;0;-M*g/2;-D/5;0;L/5;-D/5;0;L/5;-D/5;0;L/5;-D/5;0;L/5;-D/5;0;L/5];
Dist_nodo4=[-H; -W/2; -2*W; -H; W/2; -2*W; 0; 0; -2*W; 0; 0; 0; 0; -B; 0; 0; B; 0; 0; 0; -W]; 
F_total=F_ext+F;
for i=1:21
    
    Momento_y=Momento_y+F_total(i)*Dist_nodo4(i);
    Momento_vec(i,:)=F_total(i)*Dist_nodo4(i);
end

eqn=Momento_y==0;
T=double(solve(eqn));
Momento_vec=double(subs(Momento_vec));
F_total=double(subs(F_total));
F_ext=double(subs(F_ext));
D=T;
F_massa=m*g;

V_R=[10;11;12;13;18;21];%dof restringidos
U_R=[0;0;0;0;0;0];
V_L=[1;2;3;4;5;6;7;8;9;14;15;16;17;19;20];

[numRows_R,~] = size(V_R);
[numRows_L,~] = size(V_L);

K_RR=zeros(numRows_R, numRows_R);
K_LL=zeros(numRows_L, numRows_L);
K_RL=zeros(numRows_R, numRows_L);
K_LR=zeros(numRows_L, numRows_R);

for i=1:numRows_R
    for j=1:numRows_R
        K_RR(i,j)=K_G(V_R(i),V_R(j));
    end
end

for i=1:numRows_L
    for j=1:numRows_L
        K_LL(i,j)=K_G(V_L(i),V_L(j));
    end
end

for i=1:numRows_R
    for j=1:numRows_L
        K_RL(i,j)=K_G(V_R(i),V_L(j));
    end
end

for i=1:numRows_L
    for j=1:numRows_R
        K_LR(i,j)=K_G(V_L(i),V_R(j));
    end
end

F_L=zeros(numRows_L,1);
for i=1:numRows_L
   F_L(i,1)=F_total(V_L(i)); 
end


F_R=zeros(numRows_R,1);
for i=1:numRows_R
   F_R(i,1)=F_total(V_R(i)); 
end


%K_LL_inv=inv(K_LL);
U_L=K_LL\(F_L-(K_LR*U_R));

R_R=K_RR*U_R+K_RL*U_L-F_R;

U=zeros(numRows,1);
for i=1:numRows_L
    U(V_L(i,1),1)=U_L(i,1);
end
for i=1:numRows_R
    U(V_R(i,1),1)=U_R(i,1);
end

Sigma_cr=zeros(numRows,1);
u_e=zeros(4,1);
Sigma_e=zeros(numRows,1);
Epsilon_e=zeros(numRows,1);

for e=1:numRows
    x_1e=x(Tn(e,1),1);
    y_1e=x(Tn(e,1),2);
    x_2e=x(Tn(e,2),1);
    y_2e=x(Tn(e,2),2);
    z_1e=x(Tn(e,1),3);
    z_2e=x(Tn(e,2),3);
    l_e(e,:)=sqrt((x_2e-x_1e)^2+(y_2e-y_1e)^2+(z_2e-z_1e)^2);
    Ret=1/l_e(e)*[x_2e-x_1e y_2e-y_1e z_2e-z_1e 0 0 0; 0 0 0 x_2e-x_1e y_2e-y_1e z_2e-z_1e];

    if e==1 || e==8 || e==9 || e==10 || e==11 || e==17  
        Sigma_cr(e,1)=0;
    else
        Sigma_cr(e,1)=(pi^2*E_barra*In_barra)/(l_e(e).^2*A_barra);
    end
    
    
    
    for i=1:2*3
        I=Td(e,i);
        u_e(i,1)=U(I);
        
    end
    u_e_prima=Ret*u_e;
    Epsilon_e(e)=(1/l_e(e))*([-1 1]*u_e_prima);
    
    
    if e==1 || e==8 || e==9 || e==10 || e==11 || e==17 
       
    Sigma_e(e,1)=E_cable*Epsilon_e(e);
    else
    Sigma_e(e,1)=E_barra*Epsilon_e(e);
    end
end

for i=1:numRows
    if Sigma_e(i)<0 && abs(Sigma_e(i))>Sigma_cr(i)
        disp(i)
    end
end


plotBarStress3D(x,Tn,U,Sigma_e,50)














