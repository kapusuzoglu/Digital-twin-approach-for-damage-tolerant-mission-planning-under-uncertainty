%==========================================================================    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%   Copyright: Berkcan Kapusuzoglu %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This function solves the 2D crack growth optimization problem using
%   PARIS LAW for biaxial loading, % da/dN = C*(K_eff)^m - Tanaka (Pandey
%   and Patel paper 1984).
%
%
%   RBDO (RELIABILITY BASED DESIGN OPTIMIZATION)
%   Minimize P[ a(Nf) > a_crit ]
%   subject to F <= 500
%
%   RDO (ROBUST DESIGN OPTIMIZATION)
%   Minimize c*mean(a) + (1-c)*std(a)
%   subject to N >= 5000
%     
%   PARAMETERS:
%    
%       N: number of cyles
%       2a: the crack length
%       a_crit: critical crack size
%       dK: Stress intensity factor = Kmax - Kmin = F*sqrt(pi*a)
%       [MPa.m^0.5] (for uniaxial loading)
%
%       7075-T6 Aluminum alloy
%       ======================
%       Constants in Paris eq.
%       C = 4.9e-10; m = 2.69;
%       Poissons ratio, nu = 0.333
%       Yield strength, sigma_ys = 503 [MPa]
%       Shear modulus, mu = 26.9 [GPa]
%       Shear strength = 331 [MPa]
%       Fracture toughness, K_Ic = 25 [MPa]
%==========================================================================    
clc; clear; close all;

rng default % For reproducibility

% Initial crack size (half of the size of the crack since crack
% size is actually 2a in [m])
% mu_a0 = 3; cov_a0 = 0.1;

% a0 = icdf('Normal', lhsdesign(1,1), mu_a0, mu_a0*cov_a0);

% % 1st mission
% a0 = 4.826;
% a0 = icdf('Normal', lhsdesign(1,1), mu_a0, 1) % 1 mm std  -> 3.510720918736531
% 2nd mission -> obtained using InitCrackSize_Estim... file
% a0 = 6.666292293350;
% a0 = icdf('Normal', lhsdesign(1,1), mu_a0, 1) % 1 mm std  -> 2.936631361669532
% 3rd mission
% a0 = 8.590680307865;
% a0 = icdf('Normal', lhsdesign(1,1), mu_a0, 1) % 1 mm std
% 4th mission
a0=14.646725536198; % mm
% a0 = icdf('Normal', lhsdesign(1,1), mu_a0, 1) % 1 mm std

%------------------------------------------------------------------
% For material AL 7075-6:

% Paris law material parameter
m = 3.42;

% Number of samples for Monte Carlo Simulation
n_samples = 30;

time1=tic;

% Run in parallel
poolobj = parpool;

[Fsol,history] = BlockLoad(a0,m,n_samples);

delete(poolobj)

toc(time1)

save 30SamplesC_LocalOpt_Miss4
fid = fopen('Results_30SamplesC_LocalOpt_Miss4.txt', 'wt+');
format short
cl = size(history.X, 2);
for i = 1:size(history.X, 1)
  fprintf(fid, ' %.3f \t',history.X(i,1:cl));
  fprintf(fid, ' %.3f \t',history.Fval(i,1));
  fprintf(fid, ' \n');
end
fclose(fid);




