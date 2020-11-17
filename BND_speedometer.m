
function [output] = BND_sim(Input)

% By Sahand Hajimirza (sahand@rice.edu)
% This code calculates bubble nucleation and growth during Plinian silicic
% eruptions.
% See (Hajimirza et. al. 2020, EarthArxiv) for details on equations that are solved.

Par.Pi = Input.P_H2O*1e6;
Par.Pf = 1e5;
Par.T = Input.T + 273.15;
% Nucleation factor
Par.alpha = (1+cos(Input.theta*pi/180)).^2 .* (2-cos(Input.theta*pi/180)) / 4;

Par.MDR = Input.MDR;
Par.R_conduit = Input.R_conduit;



load('magma_properties.mat');
Par.solubility = solubility;
Par.viscosity = viscosity;
Par.DH2O = DH2O;
Par.fugacity = fugacity;


Par.C0 = Par.solubility(Par.Pi,Par.T);

%==========================================================================
% Constants
%==========================================================================
Par.Mw = 18.0528e-3;                    %Molar Mass of water [Kg/mol]
Par.KB = 1.38e-23;                      %Boltzman constant [J/K]
Par.AV = 6.022e23;                      %Avogadro number [1/mol]
Par.Rgas = 8.314;                       %Gas constant [J/(K.mol)]
Par.rho_melt = 2324;                    %Melt density [Kg/m^3](Ochs & Lange, Science 1999)
Par.ND = 1;                             %Minimum BND in the system
Par.g = 9.81;                           %Gravity accelaration 
%==========================================================================
Par.mu = Par.viscosity(Par.T,Par.C0);
Par.crys = Input.crys/.6;
Par.mu = Par.mu * (1 - Par.crys)^(-2.5);

Par.dpdt = Par.MDR/(pi*Par.R_conduit^2) *...
           (Par.g + Par.MDR/(pi*Par.R_conduit^2)*Par.rho_melt^(-2)...
           *(8*Par.mu/Par.R_conduit^2));
       
Par.tdec = (Par.Pi-1e5)/Par.dpdt;   

%==========================================================================


%==========================================================================
% Decompression Starts, No Bubble
%==========================================================================
option1=odeset('AbsTol',1e-5,'RelTol',1e-5,'InitialStep',1e-20,...
    'Events',@EventsFcn);
M0 = 0;     % Initial BND
[t,y,te,~,~] = ode15s(@before_nucleation_onset,[0 Par.tdec],M0,option1,Par);

[~,par1]=before_nucleation_onset(t,y,Par);
clear t y


if isempty(te)
    error('*** No bubbles formed. Try increasing saturation pressure ***')
end
 
M0 = par1.BND(end);   
r = 2 * par1.rc(end);
M1 = M0 * r;
M2 = M0 * r^2;
M3 = M0 * r^3;
pg = par1.pm(end) + 2*par1.sigma(end)/par1.rc(end);
mg = (pg*4/3*pi*r^3*Par.Mw) / (Par.Rgas*Par.T)*M0;
Cm = Par.C0 - M0/Par.rho_melt * mg;
pm = par1.pm(end);



%==========================================================================
% Bubble nucleation and growth
%==========================================================================

           
options=odeset('AbsTol',[1e-25,1e-3,1e-12,1e-20,1e-28,1e-2],...
               'RelTol',1e-5,...
               'Events',@EventsFcn3);
           


y0 = [mg,M0,M1,M2,M3,pm]';

[t,y] = ode15s(@after_nucleation_onset,[te Par.tdec],y0,options,Par);



%==========================================================================
% Output
%==========================================================================

[~,par2] = after_nucleation_onset(t,y',Par);


output.t = [par1.t; par2.t];
output.pm = [par1.pm; par2.pm];
output.cm = [par1.cm; par2.cm];
output.phi_g = [par1.phi_g; par2.phi_g];
output.J = [par1.J; par2.J];   
output.BND = [par1.BND; par2.BND];  
output.dpdt = [par1.dpdt; par2.dpdt];    
 

if abs(output.pm(end) - 1e5) < 1
    warning('*** Fragmentation condition is not reached ***')
end
    


end



function [dydt,out]=before_nucleation_onset(t,y,Par)

    M0 = y;


    pm1 = Par.Pi - Par.dpdt * t(t<Par.tdec);
    pm2 = Par.Pf * ones(size(t,1)-size(pm1,1),1);
    pm = [pm1;pm2];


    Cm = Par.C0 * ones(size(t));
    D = Par.DH2O(Cm,pm,Par.T*ones(size(pm)));    % Diffusivity
    psat = Par.Pi * ones(size(t));    % Saturation pressure
    pb = findPb(pm,psat,Par.T,Par);       % Nucleus pressure

    sigma = SurfaceTension(psat,pb,pm,Par); % Surface tension
    %=============================================
    % Nucleation rate
    rc = 2 * sigma ./ (pb-pm);
    n0 = Cm * Par.rho_melt * Par.AV / Par.Mw;
    a0 = n0.^(-1/3);
    Vw = Molecular_VH2O(pm,Par.T);
    A = (2*Vw .* n0.^2 .* D) ./ (a0) .* sqrt(sigma./(Par.KB*Par.T));
    
    
    Wcl = (16*pi * sigma.^3) ./ (3 * (pb-pm).^2);
    J = A .* exp(-Wcl * Par.alpha /(Par.KB*Par.T));
    %=============================================
    dydt = J;
    
    % Output variables
    out.t = t;
    out.pm = pm;
    out.cm = Cm;
    out.phi_g = zeros(size(t));
    out.BND = M0;
    out.dpdt = Par.dpdt * ones(size(t));
    out.J = J;
    out.u = Par.MDR ./ (pi*Par.R_conduit^2*Par.rho_melt);
    out.sigma = sigma;
    out.rc = rc;


end


function [dydt,out]=after_nucleation_onset(t,y,Par)



t = t;
mg = [y(1,:)]';

M0 = [y(2,:)]';
M1 = [y(3,:)]';
M2 = [y(4,:)]';
M3 = [y(5,:)]';
pm = [y(6,:)]';


% Conservation of mass
Cm = Par.C0 - mg/Par.rho_melt;


D = Par.DH2O(Cm,pm,Par.T*ones(size(pm)));    %Diffusivity of volatile in melt [m^2/s]

psat = find_psat(Cm,Par);
pb = findPb(pm,psat,Par.T,Par);


sigma = SurfaceTension(psat,pb,pm,Par); % Surface tension
%=============================================
% Nucleation rate
rc = 2 * sigma ./ (pb-pm);
n0 = Cm * Par.rho_melt * Par.AV / Par.Mw;
a0 = n0.^(-1/3);
Vw = Molecular_VH2O(pm,Par.T);
A = (2*Vw .* n0.^2 .* D) ./ (a0) .* sqrt(sigma/(Par.KB*Par.T));


Wcl = (16*pi*sigma.^3) ./ (3*(pb-pm).^2);
J = A .* exp(-Wcl*Par.alpha/(Par.KB*Par.T));

%=============================================
% Bubble growth

r = M1 ./ M0;  % Mean radius

% Gas pressure
V = 4/3 * pi * r.^3;
m = mg./M0;
rho_g = m./V;
pg = EoS_H2O(rho_g,Par.T);

% Growth rate
cin = Par.solubility(pg,Par.T*ones(size(pg)));
D_cin = Par.DH2O(cin,pg,Par.T*ones(size(pm)));
mu = Par.viscosity(Par.T*ones(size(cin)), cin);
drdt = (r./(4*mu)) .* (pg-pm-2*sigma./r);

mc = pb .* 4/3 * pi .* (rc).^3 / (Par.Rgas/Par.Mw * Par.T);

% Diffusion rate
dC_diff = 4*pi*r .* D_cin *Par.rho_melt .* (Cm-cin) .* M0; % dmg/dt
dC_nuc = mc .* J;

%=============================================
% Magma ascent
%Viscosity
mu = Par.viscosity(Par.T*ones(size(Cm)), Cm);
mu = mu * (1 - Par.crys)^(-2.5);

% Porisity
porosity = M3 ./ (M3 + 3/(4*pi));
dporosity = (3/(4*pi)) ./ (M3 + 3/(4*pi)).^2 ...
           .* (3 * drdt .* M2 + J .* (rc).^3);

% Magma density
rho_magma = porosity .* rho_g + (1-porosity) .* Par.rho_melt;

drhog = 3/(4*pi)*1 ./ (M0 .* r.^3) .* dC_diff...
           - rho_g ./ M0 .* J...
           - 3 * rho_g./r .* drdt;
       
drho_magma = dporosity .* (rho_g-Par.rho_melt) + drhog .* porosity;

% Decompression rate
dPdt = Par.MDR/(pi*Par.R_conduit^2) *...
           (Par.g + Par.MDR/(pi*Par.R_conduit^2)*rho_magma.^(-2)...
           .*(8 * mu/Par.R_conduit^2 - drho_magma));
       



dydt=[(dC_diff + dC_nuc)
    J
    drdt .* M0 + J .* rc
    2 * drdt .* M1 + J .* (rc).^2
    3 * drdt .* M2 + J .* (rc).^3
    -dPdt];


% Output variables
    out.t = t;
    out.pm = pm;
    out.cm = Cm;
    out.phi_g = porosity;
    out.BND = M0;
    out.dpdt = dPdt;
    out.u = Par.MDR ./ (pi*Par.R_conduit^2*rho_magma);
    out.J = J;
    
    
    
% Fragmentation criteria
    dudz = -1./rho_magma .* drho_magma;
    out.frag = dudz - .01 * 3e9./mu;        % Based on (Papale, Nature 1999) criterion
        


end




function psat = find_psat(cm,Par)


psat = zeros(size(cm));


for i = 1:length(cm)
 
  fun =  @(x) abs(Par.solubility(x,Par.T*ones(size(x))) - cm(i)); 
  
 psat(i) = fminbnd(fun,Par.Pf,Par.Pi,optimset('TolX',1e-2));   
    
    
end

end
% --------------------------------------------------------------


function pb = findPb(pm,psat,T,Par)

KB = 1.38e-23;                      %Boltzman constant [J/K]

fw_psat = Par.fugacity(psat,T*ones(size(psat)));

pb = zeros(size(pm));

for i = 1:length(pm)
    Vw = Molecular_VH2O(pm(i),T);
const = exp(Vw/KB/T*(pm(i)-psat(i)))...
    .* fw_psat(i);

p1 = pm(i);
p2 = psat(i);


fun =  @(x) abs(Par.fugacity(x,T) - const); 
pb(i) = fminbnd(fun,p1,p2,optimset('TolX',1e-1)); 


end

end



function Pg = EoS_H2O(rho,T)

% Modified Redlich and Kwong EoS for water vapor from Holloway 1977


R = 83.12;          % Gas constant cm^3.bar/(deg mole)
M = 18.01528e-3;    % Molar mass of water kg/mol

TC = T - 273.15;    % Degree Cel

ao = 35e6;
b = 14.6;
a = 166.8e6 - 193080*TC + 186.4*TC.^2 - 0.071288*TC.^3;



V = (1./rho)*M*1e6;      % molar volume (cm^3/mol)

Pg = R*T./(V-b) - a./(sqrt(T).*V.*(V+b));

Pg = Pg * 1e5;          % bar to Pa

end


function Vw = Molecular_VH2O(P,T)

AV = 6.022e23;      % Avogadro number
T = T - 273.15;     % K to C
P =  P/1e5;         % Pa to bar

% Molar volume of water (Ochs & Lange, Science 1999)
V = 22.9 + 9.5e-3*(T-1000) - 3.2e-4*(P-1);  % molar volume of water (cm^3)

Vw = V * 1e-6 / AV;

end


function ST = SurfaceTension(psat,pb,pm,Par)
% Hajimirza et al. J. Geophysical Res. Solid Earth

delta = .3200e-3;   % micro meter
alpha = .51;


Psat = psat/1e6;    % Pa to MPa
Pb = pb/1e6;
Pm = pm/1e6;    % Pa to MPa
T = Par.T - 273.15;         % K to C

ST_B = 1.21e-1 * exp(-2.24e-2*Psat)...
     + 1.47e-1 * exp(-1.90e-3*Psat)...
     + 7.5e-5 * (T-1000);
 
 ST_inf = ST_B * (1-alpha);

ST = ST_inf + delta * (Pb-Pm);

end



function [position,isterminal,direction] = EventsFcn(t,y,Par)

position = log10(y(1))-log10(Par.ND); % The value that we want to be zero
isterminal = 1;  % Halt integration 
direction = 1;   % Event function increasing
end



function [position,isterminal,direction] = EventsFcn3(t,y,Par)

mg = y(1,:);

M0 = y(2,:);
M1 = y(3,:);
M2 = y(4,:);
M3 = y(5,:);
pm = y(6,:);


% Conservation of mass
Cm = Par.C0 - mg/Par.rho_melt;
D = Par.DH2O(Cm,pm,Par.T*ones(size(pm)));    %Diffusivity of volatile in melt [m^2/s]

psat = find_psat(Cm,Par);
pb = findPb(pm,psat,Par.T,Par);


sigma = SurfaceTension(psat,pb,pm,Par); % Surface tension
%=============================================
% Nucleation rate
rc = 2 * sigma ./ (pb-pm);
n0 = Cm * Par.rho_melt * Par.AV / Par.Mw;
a0 = n0.^(-1/3);
Vw = Molecular_VH2O(pm,Par.T);
A = (2*Vw .* n0.^2 .* D) ./ (a0) .* sqrt(sigma/(Par.KB*Par.T));


Wcl = (16*pi*sigma.^3) ./ (3*(pb-pm).^2);
J = A .* exp(-Wcl*Par.alpha/(Par.KB*Par.T));

%=============================================
% Bubble growth

r = M1 ./ M0;  % Mean radius

% Gas pressure
V = 4/3 * pi * r.^3;
m = mg./M0;
rho_g = m./V;
pg = EoS_H2O(rho_g,Par.T);

% Growth rate
cin = Par.solubility(pg,Par.T*ones(size(pg)));
D_cin = Par.DH2O(cin,pg,Par.T*ones(size(pm)));
mu = Par.viscosity(Par.T*ones(size(cin)), cin);
drdt = (r./(4*mu)) .* (pg-pm-2*sigma./r);



% Diffusion rate
dC_diff = 4*pi*r .* D_cin *Par.rho_melt .* (Cm-cin) .* M0; % dmg/dt


%=============================================
% Magma ascent
%Viscosity
mu = Par.viscosity(Par.T*ones(size(Cm)), Cm);
mu = mu * (1 - Par.crys)^(-2.5);

% Porisity
porosity = M3 ./ (M3 + 3/(4*pi));
dporosity = (3/(4*pi)) ./ (M3 + 3/(4*pi)).^2 ...
           .* (3 * drdt .* M2 + J .* (rc).^3);

% Magma density
rho_magma = porosity .* rho_g + (1-porosity) .* Par.rho_melt;

drhog = 3/(4*pi)*1 ./ (M0 .* r.^3) .* dC_diff...
           - rho_g ./ M0 .* J...
           - 3 * rho_g./r .* drdt;
       
drho_magma = dporosity .* (rho_g-Par.rho_melt) + drhog .* porosity;
   
% Fragmentation criteria
    dudz = -1./rho_magma .* drho_magma;
    frag = dudz - .01 * 3e9./mu;        % Based on (Papale, Nature 1999) criterion
        

% The value that we want to be zero
% Pressure-atmospheric pressure and Papale 1999's fragmentation criterion
position = [y(6)-1e5 frag]; 
isterminal = [1 1];  % Halt integration 
direction = [-1 1];   % The zero is approached from up to down
end










