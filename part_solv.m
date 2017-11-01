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
% This version assumes just one VBS that comprises the
% entire particle (mole fraction=1, mass fraction=1)

% Eq. and table numbers correspond to Jacobson 2005 
% (Fundamentals of Atmospheric Modelling) unless otherwise 
% stated
function part_solv


% ---------------------------------------------------------
% Environmental Variables:

% temperature (K)
Temp = 298.15;
% air pressure (hPa)
Pa = 1013.25;
% time span to solve ode over (s)
tspan = [0.0 20.0];
% molecular weight dry air (g/mol) (Table 16.1)
MWda = 28.966;
% gas constant for dry air (m3.hPa/g.K) (2.24)
Rda = 8.31451e-2/MWda;
% density of dry air (g/m3 (air)) (2.23)
rho_da = Pa/(Rda*Temp);
% total volume of chamber (m3)
VT = 1.0;


% ---------------------------------------------------------
% Particle Properties:

% radius of particle (m)
rad = 1.0e-7;
% volume of one particle (m3)
Vi = (4.0/3.0)*pi*rad^3.0;
% number concentration of particles per size bin (#/m3)
Npb = 1.0e6;


% ---------------------------------------------------------
% VBS Properties:

% VBS molecular weights (g/mol (particle))
MW = 200.0;
rho = 8.0e6; % densities (g/m3)
% saturation vapour pressure (Pa (in SI base units: kg/m.s2)
P_q_ref = 1.0e-3;
% saturation concentration (ug/m3 (air)) (eq. 1 O'Meara 
% 2014), R has units kg.m2/s2.mol.K (p. 533 Jacobson 2005)
C_qg_ref = (1.0e6*MW*P_q_ref)/(8.31451*Temp);
% gas-phase concentration (ug/m3 (air))
C_gq0 = 0.0;
C_gq00 = C_gq0;
% mass fraction of particle phase comprised per VBS
MF = 1.0;
% accommodation coefficient (fraction)
alpha = 1.0;
% gas-phase diffusion coefficient of VBS bins (m2/s) (16.17) 
% corrected for non-continuum effects
Dq = Dg_calc(MW, Temp, rho_da, MWda, rad, alpha);


% ---------------------------------------------------------
% VBS and Particle Properties:

% average initial density of particle (g/m3 (particle))
rho_bar = sum(rho*MF); 
Mi0 = (rho_bar*Vi)*1.0e6; % mass of one particle (ug)
% total mass of particles (ug)
Mi0T = Mi0*Npb*VT;
% concentration of VBS in particle phase (ug/m3 (particle))
C_qp = (MF*Mi0)/Vi;
% surface concentration of components (ug/m3 (air)) corrected for
% curvature (16.33) and solute (16.36) effects
[crv1, crv2, C_sqt] = C_sqt_calc(C_qg_ref, Temp, MF, MW, rad, rho_bar, C_qp);


% ---------------------------------------------------------
% call on ode solver to find new concentration in gas
% (ug/m3 (gas)) (eq.3 Zaveri et al. 2008)
%opts = odeset('RelTol',1e-4,'AbsTol',1e-4); % solver tolerances

[~, C_gqt1] = ode45(@(tspan, C_gq0) part_ode(tspan, ...
                C_gq0, rad, Npb, Dq, rho, ...
            crv1, Mi0, VT, C_qg_ref, crv2, C_gq00, C_sqt), tspan, C_gq0);
            
% mass gained in gas phase (ug)
delta_mg = ((C_gqt1(end))-C_gq0)*VT;

% check in case everything evaporates in this time step
if delta_mg>Mi0T
    C_gqt1(end) = Mi0*Npb;
    delta_mg = ((C_gqt1(end))-C_gq0)*VT;
end

% mass fraction remaining in particle phase (ug)
mfr = (Mi0T-delta_mg)/Mi0T;
% display result
disp(mfr)

    % function that defines the ode (eq. 3 Zaveri et al. 2008)
    function C1 = part_ode(~, C_gqt, rad, Npb, Dq, rho, ...
            crv1, m_p, VT, C_qg_ref, crv2, C_gq0, C_sqt)
        
        % ----------------------------------------------------
        % inputs:
        % t - time to solve over (s)
        % C_g_i_t - initial concentrations in VBS bins in gas  
        % phase (ug/m3)
        % rad - particle radii (m)
        % Npb - number per bin (#/m3)
        % Dq - diffusion coefficient (m2/s)
        % rho - component density (g/m3)
        % C_sqt - corrected surface concentration (ug/m3 (air))
        % m_p - one particle mass (g)
        % VT - total gas-phase volume (m3)
        % C_qg_ref - reference saturation concentration (ug/m3 (air))
        % Temp - temperature (K)
        % MF - particle-phase mass fractions from VBS
        % MW - VBS molecular weight (g/mol)
        % rho_bar - average particle density (g/m3)
        % C_qp
        % ----------------------------------------------------


        % mass transfer rate (m3/s) (eq. 5 Zaveri et al. 2008)
%         k_qt = part_tran_rate(rad, Dq, Npb);
        
        
        % partitioning equation.  Gives change in concentration 
        % in gas phase (eq. 3 Zaveri et al. 2008)
        % this first equation for C1 accounts for partitioning effect on
        % particle size and therefore curvature effect and mass transfer
        % rate
        % the first line is the mass transfer rate (/s)
        % second line is the gas-phase concentration
        % third line is the saturation concentration corrected for solute
        % effect
        % fourth line is correction of saturation concentration due to
        % curvature effect
        C1 = (-4.0*pi*(((3.0/(4.0*pi))*((m_p-(((C_gqt-C_gq0)*VT)/Npb))*(1.0/(rho*1.0e6))))^(1.0/3.0))*Npb*Dq)*...
            (C_gqt-...
            C_qg_ref*1.0*...
            (exp(crv1/(crv2*(((3.0/(4.0*pi))*((m_p-(((C_gqt-C_gq0)*VT)/Npb))*(1.0/(rho*1.0e6))))^(1.0/3.0))))));
        % equation used in MOSAIC (partitioning effect on radius assumed
        % negligible and therefore mass transfer rate and curvature effect
        % assumed constant)
%         C1 = -k_qt*(C_gqt-C_sqt); 
        
        % can show radius change with time
        %disp((((3.0/(4.0*pi))*((m_p-(((C_gqt-C_gq0)*VT)/Npb))*(1.0/(rho*1.0e6))))^(1.0/3.0)))
        
        
    end
    
    % function to calculate gas-phase diffusion coefficient (16.17)
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
        % (only necessary for relatively large particles) 
        % (16.24)
        Fq = 1.0;
        % corrected gas-phase diffusion coefficient (m2/s) 
        % (16.18)
        Dq = Dq*omega*Fq;
        
    end

    % estimate effective saturation vapour pressure
    function [crv1, crv2, C_sqt] = C_sqt_calc(C_q_ref, Temp, MF, MW, ...
            rad, rho_bar, C_qp)
    
        % -------------------------------------------------
        % inputs:
        % P_qs_ref - saturation vapour pressure (Pa)
        % Temp - temperature (K)
        % MF0 - mass fraction of VBS in particle phase
        % MW - VBS molecular weights (g/mol)
        % rad - particle radius (m)
        % rho_bar - average particle density (g/m3)
        % C_qp - concentration of VBS in particle phase 
        % (ug/m3 (particle))
        % -------------------------------------------------
        
        % estimate surface tension (14.19) (see p533 if 
        % you want to include the effects of organics and 
        % inorganic ions)
        if Temp-273.15>=0
            sigma = 76.1-0.155*(Temp-273.15);
        end
        % average particle molecular weight (g/mol)
        mp_bar = sum(MF*MW);
        % curvature (Kelvin) term (16.33)
        curv = exp((2.0*sigma*mp_bar)/(rad*8.31451e3*Temp*rho_bar));
        % mole fraction in condensed phase (16.36)
        x_q = C_qp/sum(C_qp);
        % correct saturation vapour pressure to account for curvature and
        % solute effect (16.33, 16.36 and combined effect in 16.39)
        C_sqt = (C_q_ref*x_q*curv);
        crv1 = 2.0*sigma*mp_bar;
        crv2 = 8.31451e3*Temp*rho_bar;
        
        
        
    end

    % function that estimates the vapour-particle mass 
    % transfer rate (eq. 5 Zaveri et al. 2008)
    function k_qt = part_tran_rate(rad, Dq, Npb)

        % -------------------------------------------------
        % inputs:
        % rad - particle radii (m)
        % Npb - number per bin (#/m3)
        % Dg - gas-phase diffusion coefficients per bin 
        % (m2/s)
        % -------------------------------------------------

        % mass transfer rate (m3/s) (eq. 5 Zaveri et al. 2008)
        k_qt = (4.0*pi*rad*Npb*Dq);
        
    end
    % -----------------------------------------------------
    

end
