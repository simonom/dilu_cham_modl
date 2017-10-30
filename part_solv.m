% Copyright Notice: This code is in Copyright.  Any use leading to 
% publication or 
% financial gain is prohibited without the permission of the authors Simon 
% O'Meara : simon.omeara@manchester.ac.uk.  First published 2017.

% This file is part of dilu_cham_modl

% dilu_cham_modl 
% is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% dilu_cham_modl
% is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

% You should have received a copy of the GNU General Public License
% along with dilu_cham_modl
% (see the LICENSE file).  If not, see 
% <http://www.gnu.org/licenses/>.
% -------------------------------------------------------


% Eq. numbers correspond to Jacobson 2005 (Fundamentals of 
% Atmospheric Modelling)
function part_solv

% temperature (K)
Temp = 298.15;
% Volume of gas phase (m3)
VgT = 1.0;
% time span to solve ode over (s)
tspan = [0.0 20.0];
% initial gas-phase concentration of VBS bin (ug/m3)
Cgqt0 = 20.0;
% radius of particle (m)
rad = 1.0e-7;
% volume of one particle (m3)
Vi0 = (4.0/3.0)*pi*rad^3.0;
% number concentration of particles in size bin (#/m3)
Npb = 10.0;
% VBS molecular weights (g/mol) of condensed phase
% from CRC online handbook
MW = [200.0];
rho = [0.8]; % VBS densities (g/m3)
% initial mass fraction of VBS bins in particle phase
MF0 = [0.5];
% initial volume of VBS in one particle (m3)
rho_bar = sum(rho*MF0); % average initial density (g/m3)
Mit = (rho_bar*Vi0); % mass of particle (g)
Vqi0 = (MF0*Mit)/rho; % volume of each component (m3)

% accommodation coefficient (fraction)
alpha=1.0;
% air pressure (hPa)
Pa = 1013.25;
% molecular weight dry air (g/mol)
MWda = 28.966;
% gas constant for dry air (m3.hPa/g.K)
Rda = 8.31451e-2/MWda;
% density of dry air (g/m3)
rho_da = Pa/(Rda*Temp);
% gas-phase diffusion coefficient of VBS bins (m2/s) (16.17) corrected
% for non-continuum effects
Dq = Dg_calc(MW, Temp, rho_da, MWda, rad, alpha);
% effective saturation vapour pressure of components (Pa) corrected for
% curvature and solute effects
P_qs = P_qs_calc(P_qs_ref, Temp, MF0, MW, rad, rho_bar);


% call on ode solver to find new volume in particle
% (ug/m3)
[~, Vqit] = ode45(@(tspan, Vqi0) part_ode(tspan, ...
                Vqi0, rad, Npb, Dq, Vi0, Vqi0, Cgqt0, rho, VgT), tspan, Vqi0);
            
% [~, Cgqt] = ode45(@(tspan, Cgqt0) part_ode(tspan, ...
%                 Cgqt0, rad, Npb, Dgi, Vi0, Vqi0, rho, VgT), tspan, Cgqt0);

% convert to gas-phase concentration (g/m3)
Cgqt = Cgqt0-(((Vqit-Vqi0)*rho*Npb)/VgT);
            
disp(Cgqt)

    % function that defines the ode
    function V1 = part_ode(~, Vqit, rad, Npb, Dq, Vit, Vqith, Cgqt, rho, VgT)
    %function C1 = part_ode(~, Cgqt, rad, Npb, Dgi, Vit, Vqith, rho, VgT)
        % ----------------------------------------------------
        % inputs:
        % t - time to solve over (s)
        % C_g_i_t - initial concentrations in VBS bins in gas  
        % phase (ug/m3)
        % rad - particle radii (m)
        % Npb - number per bin (#/m3)
        % Vit - initial volume of one particle (m3)
        % Vqi0 - initial volume of components in one 
        % particle (m3)
        % Vqit - initial guess for new component volume 
        % (m3) in one particle
        % rho - component density (g/m3)
        % VgT - total volume of gas phase
        % ----------------------------------------------------


        % mass transfer rate (/s)
        %kimt = part_tran_rate(rad, Npb, Dgi, Vit, Vqith, Vqit, rho_da);
        % effective saturation concentration (ug/m3)
        Cstar_imt = 10.0;
        kimt=1.0;
        % code to solve the partitioning equation using an ode 
        % solver.  Eq. 16.51.  Gives change in volume of each 
        % component in particle phase and therefore loss/gain of volume
        % from/to gas phase for condensation/evaporation
        V1 = kimt*((Cgqt-(((Vqit)*rho*Npb)/VgT))-Cstar_imt);
        %C1 = -kimt*(Cgqt-Cstar_imt);
        
    end
    
    % function to calculate gas-phase diffusion coefficient 
    % (16.17)
    function Dq = Dg_calc(MW, Temp, rho_da, MWda, rad, alpha)
        
        % -------------------------------------------------
        % inputs:
        % MW - molecular weight of VBS (g/mol)
        % Temp - Temperature (K)
        % rho_da - denisty of dry air (g/m3)
        % MWda - molecular weight of dry air (g/m3)
        % rad - particle radius (m)
        % alpha - accommodation coefficient (fraction)
        
        % notes - p528 (Jacobson 2005)
        % -------------------------------------------------
        
        % gas phase diffusion coefficient (m2/s)
        Dq = ((5.0/(16.0*6.0221367e23*1.0e-9^2.0*rho_da))* ...
            ((((8.31451e3*Temp*MWda)/(2.0*pi))*(MW+MWda)/MW)^(0.5)));
        % mean free path of particles (m) (16.22)
        lam = ((MWda/(pi*6.0221367e23*1.0e-9^2.0*rho_da))* ...
            (((MWda)/(MWda+MW))^(0.5)));
        % Knudsen number (dimensionless) (16.20)
        Kn = lam/rad;
        % correction factor (16.19)
        omega = ((1.0+((1.33+0.71*Kn^-1.0)/(1.0+Kn^-1.0)+...
            (4.0*(1.0-alpha))/(3.0*alpha))*Kn)^-1.0);
        % correction factor for venting effects 
        % (only necessary for relatively large particles) (16.24)
        Fq = 1.0;
        % Corrected gas-phase diffusion coefficient (m2/s) (16.18)
        Dq = Dq*omega*Fq;
        
    end

    % estimate effective saturation vapour pressure
    function P_qs = P_qs_calc(P_qs_ref, Temp, MF0, MW, rad, rho_bar)
    
        % -------------------------------------------------
        % inputs:
        % P_qs_ref - saturation vapour pressure (Pa)
        % Temp - temperature (K)
        % MF0 - mass fraction of VBS in particle phase
        % MW - VBS molecular weights (g/mol)
        % rad - particle radius (m)
        % rho_bar - average particle density (g/m3)
        % -------------------------------------------------
        
        % estimate surface tension (14.19) (see p533 if 
        % you want to include the effects of organics and 
        % inorganic ions)
        if Temp-273.15>=0
            sigma = 76.1-0.155*Temp;
        end
        % average particle molecular weight (g/mol)
        mp_bar = sum(MF0*MW);
        % curvature (Kelvin) term (16.33)
        curv = 1.0+(2.0*sigma*mp_bar)/(rad*8.31451e3*Temp*rho_bar);
        % mole fraction in condensed phase
        x_q = C_q/C_T;
        % correct saturation vapour pressure to account for curvature and
        % solute effect (eq. 3 Riipinen et al. 2010 for solute effect, i.e.
        % Raoult's law
        % combined with 16.33 of Jacobson 2005 for for curvature effect)
        P_qs = (P_qs_ref*(x_q))*(1+curv);
        
        
        
    end

    % function that estimates the vapour-particle mass 
    % transfer rate
    function k_i_m_t = part_tran_rate(rad, Npb, Dgi, Vit, Vqith, Vqit)

        % -------------------------------------------------
        % inputs:
        % rad - particle radii (m)
        % Npb - number per bin (#/m3)
        % Dg - gas-phase diffusion coefficients per bin 
        % (m2/s)
        % Vit - new volume of components in one particle (m3)
        % Vqith - old volume of one particle (m3)
        % Vqit - old volume of components in one particle (m3)
        % -------------------------------------------------

        % new volume of components in one particle (m3)
        
        % square brackets in 16.51 
        k_i_m_t_part1 = (48.0*pi^2.0*(Vqit+Vit-Vqith))^(1.0/3.0);
        % mass transfer rate (/s)
        k_i_m_t = Dgi*k_i_m_t_part1;
        k_i_m_t = 4.0*pi*rad*Npb*k_i_m_t;
        
    end
    % -----------------------------------------------------
    

end
