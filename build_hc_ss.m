function [ A, rhs ] = build_hc_ss( Told, Tmod, hc, powdens, therm, geom )

% nb of non zero elements in the tridiagonal heat conduction matrix
dim_tq = geom.Nfuel + geom.Nclad;
nonzeros = 3*dim_tq-2;
% heat conduction matrix
A = spalloc(dim_tq,dim_tq,nonzeros);
% rhs of the heat conduction system
rhs = zeros(dim_tq,1);

%
Rfuel = therm.r_fuel;
Rclad = therm.r_clad;
Nfuel = geom.Nfuel;
Nclad = geom.Nclad;
d = Nfuel+Nclad;

hg = therm.h_gap;
% hc = hDB( hydro.h_in, hydro.pressure, hydro.speed, hydro.Dhy);

%
% create the matrix A(Told) for the fixed point method
%
% aux fuel data
auxC = zeros(Nfuel,1);
Q    = powdens*(Rfuel/(Nfuel-1))^2;
for i=1:Nfuel-1
    t = 0.5*(Told(i)+Told(i+1));
    auxC(i)=kf(t)*(2*i-1)/2;
    %auxC(i)=kkf*(2*i-1)/2;
end
% 1st fuel ring
A(1,1) = auxC(1);
A(1,2) = -A(1,1);
rhs(1) = Q/8;
% fuel, i=2..Nfuel-1
for i=2:Nfuel-1
    A(i,i)   = auxC(i)+auxC(i-1);
    A(i,i-1) = -auxC(i-1);
    A(i,i+1) = -auxC(i);
    rhs(i)   = Q*(i-1);
end
% fuel, last mesh : Nfuel
A(Nfuel,Nfuel-1) = -auxC(Nfuel-1);
A(Nfuel,Nfuel  ) = auxC(Nfuel-1)+Rfuel*hg;
A(Nfuel,Nfuel+1) = -Rfuel*hg;
rhs(Nfuel) = Q*(Nfuel-1-1/4)/2;

% aux clad data
clear auxC;
auxC = zeros(Nclad,1);
for j=1:Nclad-1
    t = 0.5*(Told(Nfuel+j)+Told(Nfuel+j+1));
    auxC(j) = kc(t)*(Rfuel/(Rclad-Rfuel)*(Nclad-1)+(2*j-1)/2);
    %auxC(j)=kkc*(Rfuel/(Rclad-Rfuel)*(Nclad-1)+(2*j-1)/2);
end
% j=1
A(Nfuel+1,Nfuel+1) = auxC(1)+Rfuel*hg;
A(Nfuel+1,Nfuel  ) = -Rfuel*hg;
A(Nfuel+1,Nfuel+2) = -auxC(1);
% j=2..Nclad-1
for j=2:Nclad-1
    A(Nfuel+j,Nfuel+j  ) = auxC(j)+auxC(j-1);
    A(Nfuel+j,Nfuel+j-1) = -auxC(j-1);
    A(Nfuel+j,Nfuel+j+1) = -auxC(j);
end
% j=Nclad
A(Nfuel+Nclad,Nfuel+Nclad  ) = auxC(Nclad-1)+Rclad*hc;
A(Nfuel+Nclad,Nfuel+Nclad-1) = -auxC(Nclad-1);

rhs(d) = Tmod*Rclad*hc;


