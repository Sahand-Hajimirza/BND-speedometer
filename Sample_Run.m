function Sample_Run

% This clears command window 
clc
clear

%============================================================
% Input variables, from supplementary table 1.

P_H2O = 100;                % Saturation pressure (MPa) [50-300]
T  = 860;                   % Temperature (C)  [700-1100]
crys = 0;                   % Crystal content [0 1]
theta = 155;                % Heterogeneous nucleation contact angle (degrees) [0 180]
MDR = 10e7;                  % Mass discharge rate (Kg/s) 

% Vary conduit radius until BND matches observations
R_conduit = 30;             % Conduit radius (m)

Input = table(MDR,P_H2O,T,MDR,crys,theta,R_conduit);
%============================================================


[output] = BND_sim(Input); 

disp('%===================================')
disp('Final bubble number density (m^-3):')
disp(output.BND(end))
disp('%===================================')
disp('Time-averaged decompression rate (MPa/s):')
dpdt = (output.pm(1) - output.pm(end))/output.t(end)/1e6;
disp(dpdt)
disp('%===================================')

end

