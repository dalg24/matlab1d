function [MSH, XS, therm, hydro ] = load_input(ref_nz);

% geometry input
% we enter the XS ID for the material, this creates the 1D material mesh
% the ID numbers correspond to the Casmo XS of the MSLB benchmark form the
% OECD
i_mat=[24 289 290 291 292 293 294 295 295 295 295 296 296 296 296 298 298 298 298 -299 -300 -301 -302 -303 26 ]';
% i_mat=289*ones(20,1);

[MSH, therm, hydro] = load_mesh_TH_data( i_mat, ref_nz );

% load XS: we can either read the entire file containing all of the X
% for all materials (set read_entire_file=true) or, if you have already
% read it, you can simply read the smaller file that contains only the one
% in your current geometry (read_entire_file=false).
% Note: every time you change the geometry, you should reset
% read_entire_file to true

read_entire_file=false;
XS =load_XS_all(read_entire_file,MSH.z2med);
% perform a check
if( (read_entire_file==false) && (any(i_mat-XS.z2medORI)))
    disp('you have changed the geometry, the entire XS file will be read');
    disp('and the XS needed for your new geomtrical layout will be saved');
    read_entire_file=true;
    XS =load_XS_all(read_entire_file,MSH.z2med);
end
% make the rodded XS less absorbing
blending=true;
if(blending)
    weight = 0.5;
    XS = blend_XS( XS, MSH, weight );
end

imat2(1:MSH.nz,1)= XS.imat;
MSH.z2mat = imat2;


