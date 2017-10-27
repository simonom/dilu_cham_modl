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
Cgqt = 20.0;
% radius of particle (m)
rad = 1.0e-7;
% volume of one particle (m3)
Vi0 = (4.0/3.0)*pi*rad^3.0;
% number concentration of particles in size bin (#/m3)
Npb = 10.0;
% VBS molecular weights (g/mol)
% from CRC online handbook
MW = [200.0];
rho = [0.8]; % VBS densities (g/m3)
% initial mass fraction of VBS
MF0 = [0.5];
% initial volume of VBS in one particle (m3)
rho_bar = np.sum(rho*MF0); % average initial density (g/m3)
Mit = (rho_bar*Vit); % mass of particle (g)
Vqi0 = (MF0*Mit)/rho; % volume of each component (m3)
% gas-phase diffusion coefficient of VBS bins (m2/s)
Dgi = Dg_calc(MW, rho, Temp);

% call on ode solver to find new gas_phase concentration
% (ug/m3)
[~, Cgqt] = ode45(@(tspan, Vqi0) part_ode(tspan, ...
                Vqi0, rad, Npb, Dgi, Vi0, Vqi0, rho), tspan, Vqi0);

disp(Cgqt)

    % function that defines the ode
    function Cgqt = part_ode(~, Vqit, rad, Npb, Dgi, Vit, Vqith, rho, VgT)

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
        k_i_m_t = part_tran_rate(rad, Npb, Dgi, Vit, Vqith, Vqit);
        % effective saturation concentration (ug/m3)
        Cstar_i_m_t = 10.0;

        % code to solve the partitioning equation using an ode 
        % solver.  Eq. 16.51.  Gives change in volume of each 
        % component in particle phase and therefore loss/gain of volume
        % from/to gas phase for condensation/evaporation
        Vqit = -k_i_m_t(1)*(C_g_i_t(1)-Cstar_i_m_t(1));
        % convert to gas-phase concentration change (g/m3)
        Cgqt = (Vqit*rho*Npb)/VgT;
        
    end
    
    % function to calculate gas-phase diffusion coefficient 
    % (p44 of thesis2)
    function Dgi = Dg_calc(MW, rho, Temp)
        
        % -------------------------------------------------
        % inputs:
        % MW - molecular weight of VBS (g/mol)
        % rho - density of VBS (g/m3)
        % Temp - Temperature (K)
        
        % notes - equation reference is p44 of thesis2
        % -------------------------------------------------
        
        % (p44 of thesis2):
        % molecular volume (m3/molecule)
        Mvol = (MW/rho)/6.022e23; 
        % molecular diameter (m)
        sigma = (((3*Mvol)/(4*pi))^(1.0/3.0))*2;
        % molecular mass (kg/molecule) (convert from g/mol)
        Mmolec = (MW/6.022e23)*1e-3;
        % mean molecular speed in gas phase (m/s)
        molv = ((3.0*Temp*1.381e-23)./Mmolec).^(0.5);
        % air pressure (kg/m.s2)
        Pair = 101325.0;
        % ideal gas constant (kg/m2.s2.mol.K)
        R = 8.31;
        % number density (molecules/m3)
        n_rho = ((Pair)/(R*Temp))*6.022e23;  
        % mean free path of molecules (m)
        lam = 0.707/(pi*n_rho*(sigma^2));
        % gas phase diffusion coefficient (m2/s) 
        % (p44 of thesis2)
        Dgi = (molv*lam)/3.0;
        
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
