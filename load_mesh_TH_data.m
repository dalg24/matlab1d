function [MSH, therm, hydro] = load_mesh_TH_data( i_mat, ref_nz);


% nbr of axial regions
[nreg] = length(i_mat);
if(nreg<=0)
    error(sprintf('wrong length of i_mat input: nreg =%i',nreg))
end
% 1D problem
MSH.dim = 1;

% fuel assembly size, cm
fa_size = 21.504;
% core height, cm
core_height  = 400   ;


% computational mesh size (in cm's because neutronics XS and inverse cm's)
dz = core_height/(ref_nz*nreg);

% mesh to medium array
z2med=zeros(nreg*ref_nz,1);
for ireg=1:nreg
    i1=(ireg-1)*ref_nz+1;
    i2=ireg*ref_nz;
    z2med(i1:i2,1)=i_mat(ireg);
end

% nb of spatial axial cells
nz = nreg *ref_nz;
% create axial positions (in cm)
z  = linspace(0,core_height,nz+1);

% create mesh structure
MSH.z  = z;
MSH.nz = nz;
MSH.dz = dz;
MSH.z2med  = z2med;

[therm,hydro]=load_therm_hydr( fa_size, core_height );

therm.Watt0 = therm.powass; % Watts of thermal power at ss

%%%%%%%%%%%%%%%%%%%%%%%%%%% HEAT CONDUCTION GEOMETRY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           ------------------------
% # of fuel nodes , of clad nodes   
Nfuel = 4;
Nclad = 2;
MSH.radial.Nfuel = Nfuel ;   
MSH.radial.Nclad = Nclad ; 
MSH.radial.Ntot  = Nclad + Nfuel ; 

% create fuel and clad mesh
for i=1:Nfuel,
    ri(i) = (i-1) * therm.r_fuel / ( Nfuel-1 );
end

for j=1:Nclad,
    rj(j) = therm.r_fuel + (j-1) *( therm.r_clad - therm.r_fuel ) / ( Nclad - 1 );
end
r=ri; 
r(Nfuel+1:Nfuel+Nclad)=rj;

MSH.radial.r=r;
