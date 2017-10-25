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
function part_solv_swe

% time span to solve ode over (s)
tspan = [0.0 20.0];
% initial gas-phase concentration of VBS bin (ug/m3)
C_g_i_t = 20.0;

% call on ode solver to find new gas_phase concentration
% (ug/m3)
[~, C_g_i_t] = ode45(@part_ode, tspan, C_g_i_t);
disp(C_g_i_t)

    % function that defines the ode
    function C1 = part_ode(~, C_g_i_t)

        % ----------------------------------------------------
        % inputs:
        % t - time to solve over (s)
        % C_g_i_t - initial concentrations in VBS bins in gas  
        % phase (ug/m3)
        % ----------------------------------------------------

        % number of size bins (i)
        sbn = length(C_g_i_t);
        % mass transfer rate (/s)
        k_i_m_t = zeros(sbn);
        k_i_m_t(:) = 1.0;
        % effective saturation concentration (ug/m3)
        Cstar_i_m_t = 10.0;

        % code to solve the partitioning equation using an ode 
        % solver.  Eq. 16.51
        C1 = -k_i_m_t(1)*(C_g_i_t(1)-Cstar_i_m_t(1));
    end
end
