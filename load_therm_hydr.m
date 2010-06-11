function [thermiq,hydro]=load_therm_hydr(fa_size, height );

%----------------------------%
% heat conduction data
%----------------------------%
% units: meter
% fuel pellet radius, clad radius
thermiq.r_fuel = 0.410E-2; 
thermiq.r_clad = 0.475E-2; 

thermiq.rho_fuel = 10.8E3;
thermiq.rho_clad =  7.6E3;

% water-hole radius, assembly pitch
thermiq.r_wh      = 0.6025E-2;
thermiq.ass_pitch = fa_size/100; % convert to meter

% compute power density W/m^3
Ptot   = 2750E6 ; % typical total power (in W) of a PWR
Nass   = 157;     % that contains 157 fuel assemblies
Npins  = 264;     % each FA contains 264 fuel pins
Nwh    = 25;      % and 25 water holes
thermiq.height = height/100; % convert to meter
thermiq.powass  = Ptot/Nass;  % average nominal power in a single FA
thermiq.pinsurf = pi*(thermiq.r_fuel)^2; % pin area
% average nominal power density in a single fule element
thermiq.powdens = thermiq.powass/(Npins*height*thermiq.pinsurf);  % units: W/m3

thermiq.pindens = Npins/fa_size; % pin per cm of radial direction


% gap conductance
thermiq.h_gap = 1E4; % W/m^2

% T effective in fuel = rw(1) T fuel centerline + rw(2) T fuel surface
thermiq.rw(1:2) = [4, 5]/9;


%----------------------------%
% hydraulical parameters
%----------------------------%
% coolant inlet temperature, pressure, inlet speed
hydro.T_in     = 292;     % in Celsius
hydro.pressure = 155E+05; % Pa
hydro.speed    = 4.9;     % m/s

% compute the inlet enthalpy
[hydro.h_in]   = hlsat(hydro.pressure , hydro.T_in);
% compute the inlet density
[hydro.rho_in] = rholiq(hydro.h_in, hydro.pressure);
% compute the mass flux, in kg/(m^2.s)
hydro.mass_flux = hydro.rho_in * hydro.speed; 
% compute the hydraulic flow area: it is the area of the sqaure FA 
% minus the area of the water holes minus the area of the fuel elements
% Shy=a*a-264*pi*Rclad^2-25*pi*Rtg^2;
hydro.S_hy      = thermiq.ass_pitch^2 ...
    - pi* (Nwh*thermiq.r_wh^2 + Npins*thermiq.r_clad^2);
% compute the mass flow rate, in kg/s
hydro.mass_flow_rate = hydro.mass_flux * hydro.S_hy ;

% compute the wetted parameter for an average pin
hydro.pwet =2*pi*( Npins*thermiq.r_clad + Nwh*thermiq.r_wh);
% compute the hydraulic diameter
hydro.Dhy  = 4* hydro.S_hy / hydro.pwet;

% % wall-fluid exchange coefficient
% thermiq.h_conv = hDB( hydro.h_in, hydro.pressure, hydro.speed, hydro.Dhy);

