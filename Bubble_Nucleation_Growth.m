function [BND,varargout] = Bubble_Nucleation_Growth_tData(In)

% This code is written by Sahand Hajimirza. For detail decription of model
% see "Hajimirza et. al. 2019, JGR-Solid Earth".

% For scientific and non-commercial use only.

% Please cite "Hajimirza et. al. 2019, JGR-Solid Earth" upon using this code.

% This code estimates bubble nucleation and growth in rhyolitic melt during
% decompression in nucleation experiments. The volatile phase is water.

% Input: A table that includes experimental conditions. The required
% experimental conditions are initial and final pressures in [MPa], a
% constant decompression rate in [MPa/s], temperature, only 825, 850, 875 [C],
% and decompression and post decompression times in [s]. 

% Output: Final bubble number density, BND, in number/m^3 melt. If two
% outputs is called, the evolution of model parameters in time will be the
% second output. 


global fw_interpolant


%==========================================================================
% Experimental conditions
%==========================================================================
Par.Pinitial = In.Pi*1e6;           % Initial pressure [Pa] 
Par.pf = In.Pf*1e6;                 % Finale pressure [Pa]
Par.dpdt = -In.dPdt * 1e6;          % Decompression rate [Pa/s] 
Par.tdec = In.td;                   % Decompression time [s]
Par.tpost = In.tpost;               % Post decompression time [s]
Par.T = In.T + 273.15;              % Temperature [K]

Par.Xc = 0;                         % Mole fraction of CO2. Will be developed in future versions. 
Par.C0 = solubility(Par.Pinitial, Par.Xc, Par.T);           % Initial water concentration
Par.Cf = solubility(Par.pf, Par.Xc, Par.T);                 % Final water concentration

%==========================================================================
% Sample Properties & Constants
%==========================================================================
Par.Mw = 18.0528e-3;                % Molar Mass of water [Kg/mol]
Par.KB = 1.38e-23;                  % Boltzman constant [J/K]
Par.AV = 6.022e23;                  % Avogadro number [1/mol]
Par.Rgas = 8.314;                   % Gas constant [J/(K.mol)]
Par.rho_melt = 2324;                % Melt density [Kg/m^3](Ochs & Lange, Science 1999)
Par.ND = 1;                         % Detectability Threshold [m^-3](Cluzel etal, 2008)
%==========================================================================


% To accelerate the run, fugacity coefficient of water is calculated for
% three temperatures, 825, 850, and 875 C. The model uses interpolation of
% the calculated data.

if strcmp(num2str(Par.T-273.15),'825')
    load('fw_interpolant.mat','fw_interpolant_825C')
    fw_interpolant = fw_interpolant_825C;
    clear fw_interpolant_825C;
elseif strcmp(num2str(Par.T-273.15),'850')
    load('fw_interpolant.mat','fw_interpolant_850C')
    fw_interpolant = fw_interpolant_850C;
    clear fw_interpolant_850C;
elseif strcmp(num2str(Par.T-273.15),'875')
    load('fw_interpolant.mat','fw_interpolant_875C')
    fw_interpolant = fw_interpolant_875C;
    clear fw_interpolant_875C;
else
    error('Temperature is not in the defined Ts')
end




%==========================================================================
% Decompression Starts, No Bubble
%==========================================================================
% The first step in the code is before bubbles form. 

option1=odeset('AbsTol',1e-5,'RelTol',1e-5,'InitialStep',1e-20,...
    'Events',@EventsFcn);
M0 = 0;     % Initial BND
[t,y,te,~,~] = ode15s(@before_nucleation,[0 Par.tdec+Par.tpost],M0,option1,Par);

 
[~,par1]=before_nucleation(t,y,Par);
clear t y

% If no bubble form in the system during experimental duration, code
% returns 0 bubble number density
if isempty(te)
    BND = 0;
    Out = par1;
    Out.Csol = solubility(Par.pm,Par.Xc,Par.T);
    
    return
end
 
% Calculating initil condtions for the next step
M0 = par1.M0(end);   
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
% The second step is after bubbles for in the system. 
           
options=odeset('AbsTol',[1e-25,1e-3,1e-12,1e-20,1e-28,1e-2],...
               'RelTol',1e-5);

          
y0 = [mg,M0,M1,M2,M3,pm]';      % Initial condition for the second step

[t,y] = ode15s(@after_nucleation,[te Par.tdec+Par.tpost],y0,options,Par);


%==========================================================================
% Output
%==========================================================================
BND = y(end,2);

if nargout > 1

    [~,par2] = after_nucleation(t,y',Par);
    parnames = fieldnames(par1);

    for i = 1:length(parnames)
        Out.(char(parnames(i))) = [par1.(char(parnames(i))) ; par2.(char(parnames(i)))];
    end
    Out.Csol = solubility(Out.pm,Par.Xc,Par.T);
    varargout{1} = Out;

end



end



function [dydt,out]=before_nucleation(t,y,Par)
% ODE function for the first step, before bubbles form.

global fw_interpolant

    out.t = t;
    out.M0 = y;

    % Melt pressure
    pm1 = Par.Pinitial + Par.dpdt * t(t<Par.tdec);
    pm2 = Par.pf * ones(size(t,1)-size(pm1,1),1);
    out.pm = [pm1;pm2];


    out.Cm = Par.C0 * ones(size(t));
    out.D = DH2Orhyolite(out.Cm,out.pm,Par.T);    %Diffusivity of volatile in melt [m^2/s]
    out.psat = Par.Pinitial * ones(size(t));      % Saturation pressure
    out.pb = findPb(out.pm,out.psat,Par.T,fw_interpolant);  % Nucleus pressure
    %=============================================

    out.sigma = SurfaceTension(out,Par.T); 
    out.rc = 2 * out.sigma ./ (out.pb - out.pm);
    n0 = out.Cm * Par.rho_melt * Par.AV / Par.Mw;
    a0 = n0.^(-1/3);
    Vw = Molecular_VH2O(out.pm,Par.T);
    A = (2*Vw .* n0.^2 .* out.D) ./ (a0) .* sqrt(out.sigma./(Par.KB*Par.T));
    Wcl = (16*pi * out.sigma.^3) ./ (3 * (out.pb-out.pm).^2);
   
    out.J = A .* exp(-Wcl./(Par.KB*Par.T));     % Nucleation rate from CNT
    
    dydt = out.J;
    
    out.cin = NaN * ones(size(t));
    out.pg = NaN * ones(size(t));
    out.dC_diff = zeros(size(t));
    out.dC_nuc = zeros(size(t));
    out.mg = zeros(size(t));
    out.M1 = zeros(size(t));
    out.M2 = zeros(size(t));
    out.M3 = zeros(size(t));

end


function [dydt,out]=after_nucleation(t,y,Par)
% ODE function for the first step, when bubbles present in the system.
global fw_interpolant

    out.t = t;
    out.mg = [y(1,:)]';

    out.M0 = [y(2,:)]';
    out.M1 = [y(3,:)]';
    out.M2 = [y(4,:)]';
    out.M3 = [y(5,:)]';
    out.pm = [y(6,:)]';
    out.pm(out.pm<1e5) = 1e5;


    out.Cm = Par.C0 - out.mg/Par.rho_melt;
    out.dPdt = -Par.dpdt;

    % Nucleation rate
    out.D = DH2Orhyolite(out.Cm,out.pm,Par.T);    %Diffusivity of volatile in melt [m^2/s]
    out.psat = find_psat(out.Cm,Par.T,Par.Xc,Par.Pinitial);

    out.pb = findPb(out.pm,out.psat,Par.T,fw_interpolant);

    out.sigma = SurfaceTension(out,Par.T); 

    out.rc = 2*out.sigma ./ (out.pb-out.pm);
    n0 = out.Cm * Par.rho_melt * Par.AV / Par.Mw;
    a0 = n0.^(-1/3);
    Vw = Molecular_VH2O(out.pm,Par.T);
    A = (2*Vw .* n0.^2 .* out.D) ./ (a0) .* sqrt(out.sigma/(Par.KB*Par.T));
    out.Wcl = (16*pi*out.sigma.^3) ./ (3*(out.pb-out.pm).^2);

    out.J = A .* exp(-out.Wcl/(Par.KB*Par.T));

    % Bubble growth rate
    r = out.M1 ./ out.M0;

    V = 4/3 * pi * r.^3;
    m = out.mg./out.M0;
    rho = m./V;
    out.pg = EoS_H2O(rho,Par.T);
    out.cin = solubility(out.pg,Par.Xc,Par.T);
    out.D_cin = DH2Orhyolite(out.cin,out.pg,Par.T);
    mu = viscosity(Par.T,out.cin);
    drdt = (r./(4*mu)) .* (out.pg-out.pm-2*out.sigma./r);

    % Diffusion rate
    mc = out.pb .* 4/3 * pi .* (out.rc).^3 / (Par.Rgas/Par.Mw * Par.T);
    out.dC_diff = 4*pi*r .* out.D_cin *Par.rho_melt .* (out.Cm-out.cin) .* out.M0; % dmg/dt
    out.dC_nuc = mc .* out.J;

    dydt=[(out.dC_diff + out.dC_nuc)
           out.J
           drdt .* out.M0 + out.J .* out.rc
           2 * drdt .* out.M1 + out.J .* (out.rc).^2
           3 * drdt .* out.M2 + out.J .* (out.rc).^3
           -out.dPdt];

end



function [Cw,Cc] = solubility(P,Xc,T)
% Solubility of water in rhyolite. (Liu et al 2005)
    Xw = 1 - Xc;
    Pw = Xw .* P/1e6;
    Pc = Xc .* P/1e6;
  
    a1 = 354.94;
    a2 = 9.623;
    a3 = -1.5223;
    a4 = 0.0012439;
    a5 = -1.084e-4;
    a6 = -1.362e-5;
    %
    b1 = 5668;
    b2 = 0.4133;
    b3 = 2.041e-3;
    b4 = -55.99;
    
    
    Cw = (a1*Pw.^(1/2) + a2*Pw + a3*Pw.^(3/2))./T + a4*Pw.^(3/2) +...
      Pc.*(a5*Pw.^(1/2) + a6*Pw); %H2O content in wt.%
    Cc = b1*Pc/T + Pc.*(b2*Pw.^(1/2) + b3*Pw.^(3/2)) + b4*Pc.*Pw/T;
  
    Cw = Cw/100; % fraction
    Cc = Cc/1e6; % fraction
    
end


function psat = find_psat(cm,T,Xc,Pinitial)
% Calculates saturation pressre from solubility curve.

psat = zeros(size(cm));

p2 = Pinitial;
p1 = 1e5;

for i = 1:length(cm)
 
  fun =  @(x) abs(solubility(x,Xc,T) - cm(i)); 
  
 psat(i) = fminbnd(fun,p1,p2,optimset('TolX',1e-2));   
    
    
end


end


function eta = viscosity(T,H2O)
% Viscosity of rhyolite (Hui & Zhang 2007)

SiO2 = 76.53e-2;
TiO2 = .06e-2;
Al2O3 = 13.01e-2;
FeO = .79e-2;
MnO = .08e-2;
MgO = .02e-2;
CaO = .74e-2;
Na2O = 3.87e-2;
K2O = 4.91e-2;

mO = 15.9994;
mSiO2 = 28.0855 + 2*mO;
mTiO2 = 47.88 + 2*mO;
mAl2O3 = 2*26.98154 + 3*mO;
mFeO = 55.847 + mO;
mMnO = 54.9380 + mO;
mMgO = 24.305 + mO;
mCaO = 40.08 + mO;
mNa2O = 2*22.98977 + mO;
mK2O = 2*39.0983 + mO;
mH2O = 2*1.00794 + mO;

nSiO2 = SiO2 / mSiO2;
    nTiO2 = TiO2 / mTiO2;
    nAl2O3 = Al2O3 / mAl2O3;
    nFeO = FeO / mFeO;
    nMnO = MnO / mMnO;
    nMgO = MgO / mMgO;
    nCaO = CaO / mCaO;
    nNa2O = Na2O / mNa2O;
    nK2O = K2O / mK2O;
%    nP2O5 = P2O5 / mP2O5;
    nH2O = H2O / mH2O;
    nFeMnO = nFeO + nMnO;
    
    nNaK = 2*nNa2O + 2*nK2O;
    nAl = 2*nAl2O3;
    if nNaK <= nAl
        nNaKAlO2 = nNaK;
        nAl2O3ex = (nAl-nNaK)/2;
        nNaK2Oex = 0;
    else
        nNaKAlO2 = nAl;
        nNaK2Oex = (nAl-nNaK)/2;
        nAl2O3ex = 0;
    end
    
    nmol = nSiO2 + nTiO2 + nFeMnO + ...
        nMgO + nCaO + nNaKAlO2 + nAl2O3ex + nNaK2Oex + nH2O;
    %
    XSiO2 = nSiO2./nmol;
    XTiO2 = nTiO2./nmol;
    XFeMnO = nFeMnO./nmol;
    XMgO = nMgO./nmol;
    XCaO = nCaO./nmol;
 %   XP2O5 = nP2O5./nmol;
 XP2O5 = 0;
    XNaKAlO2 = nNaKAlO2./nmol;
    XAl2O3ex = nAl2O3ex./nmol;
    XNaK2Oex = nNaK2Oex./nmol;
    XH2O = nH2O./nmol;
    %
    Z = XH2O.^(1./(1 + 185.797./T));
    
    eta = 10.^( ...
        ( ...
        -  6.83*XSiO2 ...
        - 170.79*XTiO2 ...
        - 14.71*XAl2O3ex ...
        - 18.01*XMgO ...
        - 19.76*XCaO ...
        - 34.31*XNaK2Oex ...
        - 140.38*Z ...
        + 159.26*XH2O ...
        - 8.43*XNaKAlO2 ...
        ) ...
        + ...
        ( ...
        + 18.14*XSiO2 ...
        + 248.93*XTiO2 ...
        + 32.61*XAl2O3ex ...
        + 25.96*XMgO ...
        + 22.64*XCaO ...
        - 68.29*XNaK2Oex ...
        + 38.84*Z ...
        - 48.55*XH2O ...
        + 16.12*XNaKAlO2 ...
        )* 1000./T ...
        + ...
        exp( ...
        ( ...
        + 21.73*XAl2O3ex ...
        - 61.98*XFeMnO ...
        - 105.53*XMgO ...
        - 69.92*XCaO ...
        - 85.67*XNaK2Oex ...
        + 332.01*Z ...
        - 432.22*XH2O ...
        - 3.16*XNaKAlO2 ...
        ) ...
        + ...
        ( ...
        + 2.16*XSiO2 ...
        - 143.05*XTiO2 ...
        - 22.10*XAl2O3ex ...
        + 38.56*XFeMnO ...
        + 110.83*XMgO ...
        + 67.12*XCaO ...
        + 58.01*XNaK2Oex ...
        + 384.77*XP2O5...
        - 404.97*Z ...
        + 513.75*XH2O ...
        )*1000./T ...
        ) ...
        );
end


function D = DH2Orhyolite(concentration,pm,T)
% Diffusivity of water in rhyolite, (Zhang & Behrens 2007) or (Ni & Zhang 2008)

diffusivity = 'Zhang_Behrens_2000';
%
if strcmp(diffusivity,'Zhang_Behrens_2000')
    %H2O bulk diffusivity (D in m^2/s) in rhyolitic melt after Zhang and Behrens (2000)
    %
    %T is in Kelvin
    P = pm/1e6; %from Pa to MPa
    C = 100*concentration; % from fractional to percent
    X = (C/18.015) ./ (C/18.015 + (100-C)/32.49); %mole fraction of total H2O on a single oxygen basis
    m = -20.79 - 5030./T - 1.4*P./T;

    D = 1e-12 .* X .* exp(m) .* ( ...
        1 +exp( ...
        56 + m + X.*(-34.1 + 44620./T + 57.3*P./T) ...
        - sqrt(X).*(0.091 + 4.77e6./T.^2) ...
        ) ...
        );
elseif strcmp(diffusivity,'Ni_Zhang_2008')
    %H2O bulk diffusivity (D in m^2/s) in rhyolitic melt after
    %Ni and Zhang, Chemical Geology, 250, 68-78 (2008)
    %
    %T is in Kelvin
    P = pm/1e9; %from Pa to GPa
    C = 100*concentration; % from fractional to percent
    X = (C/18.015) ./ (C/18.015 + (100-C)/32.49); %mole fraction of total H2O on a single oxygen basis

    D = 1e-12 * X.*exp( ...
        13.47 - 49.996*X + 7.0827*sqrt(X) + 1.8875*P ...
        - (9532.3 - 91933*X + 13403*sqrt(X) + 3625.6*P)./T ...
        );


end

end


% --------------------------------------------------------------




function pb = findPb(pm,psat,T,fw_interpolant)
% Finds the pressure inside bubbles. 

KB = 1.38e-23;                      %Boltzman constant [J/K]

pb = zeros(size(pm));

for i = 1:length(pm)
    Vw = Molecular_VH2O(pm(i),T);
const = exp(Vw/KB/T*(pm(i)-psat(i)))...
    .* fw_interpolant(psat(i));

p1 = pm(i);
p2 = psat(i);

fun =  @(x) abs(fw_interpolant(x) - const); 
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
% Molar volume of water (Ochs & Lange, Science 1999)
AV = 6.022e23;      % Avogadro number
T = T - 273.15;     % K to C
P =  P/1e5;         % Pa to bar


V = 22.9 + 9.5e-3*(T-1000) - 3.2e-4*(P-1);  % molar volume of water (cm^3)
Vw = V * 1e-6 / AV;

end





% --------------------------------------------------------------


function [position,isterminal,direction] = EventsFcn(t,y,Par)
% Event function to detect detectibility threshold
position = log10(y(1))-log10(Par.ND); % The value that we want to be zero
isterminal = 1;  % Halt integration 
direction = 1;   % Event function increasing
end





function ST = SurfaceTension(Par,T)
% Surface tension (Hajimirza et al 2019)

delta = .3200e-3;       % micro meter
alpha = .51;

Psat = Par.psat/1e6;    % Pa to MPa
Pb = Par.pb/1e6;
Pm = Par.pm/1e6;        % Pa to MPa
T = T - 273.15;         % K to C

ST_B = 1.21e-1 * exp(-2.24e-2*Psat)...
     + 1.47e-1 * exp(-1.90e-3*Psat)...
     + 7.5e-5 * (T-1000);
 
 ST_inf = ST_B * (1-alpha);

ST = ST_inf + delta * (Pb-Pm);

end





