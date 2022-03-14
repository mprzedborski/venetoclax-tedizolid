%==========================================================================
% This is the script for running the optimal drug timing simulations for
% MOLM-13 R2 AML cells
%==========================================================================
% Note: Here, timescale is minutes (or units of t_scale), and concentration 
% scale is nM.
%==========================================================================
function rSq = Run_drug_timing_simulations()
rng default
rSq = 0.0;
close all

%--------------------------------------------------------------------------
% Initial drug levels:
%--------------------------------------------------------------------------
ABT_0 = 400.0;  %nM
Ted_0 = 5000.0; %nM
t_scale = 12.0*60.0;       % t-scale is the number of minutes in 12 hours.

%--------------------------------------------------------------------------
% Bcl-2 family protein parameter values:
%--------------------------------------------------------------------------
myPars = [13.669755696602916, 485.697261747894288, 1044.141691325814236,...
          1153.185096086854855, 7.192553700040385, 82.687472599671310,...
          452.154079358190415, 0.161810685573990, 3808.779122436302714,...
          129.902673437969696, 2.320021253551643, 11.452565945560234,...
          84.164288956134570, 1.729787988695099, 3.151231128658153,...
          0.602022914394669, 3.520635994758134, 1.440200545088872,...
          20.192519980890044, 14.705005047381947, 25.637206631177463,...
          1.950200147859286, 0.735094861939094, 0.846708465943897,...
          1.726707485864304, 0.530205713378790, 1.638885619999066,...
          1.858928933236780, 3.438405173378365, 8.080797703919686,...
          10.536611750584026, 3.706640794850193, 2.721147803974464,...
          3.922800957032131, 2.377366569413438, 4.620362679300571,...
          6.817649060433075, 2.060800458146430, 10.993215713647745,...
          0.521150510527861, 2.002345966157941, 2.727210293394596,...
          2.271195899462174, 6.699538751856998, 3.484730980617899,...
          3.721176559333263, 1.025447266055947, 1.702891500541281,...
          0.140738806893586, 1.019986323110426, 0.014757110602751];
      
%--------------------------------------------------------------------------
% Switches:
proceed = 0;
proceed_2 = 0;
proceed_3 = 0;

%==========================================================================
% (1) Read in the experimental data
%==========================================================================
% Protein data file has the following row order: MCL-1; BCL-2; BCL-XL; BIM; 
% BAX; BAK; aCASP-3; 
num_protocols = 4;         % number of treatment protocols + untreated case.
fname_protein = 'Protein.dat';         % Protein data file.

data_format = '%f';                             % data file formats (or e.g. %16f)
fid_protein = fopen(fname_protein, 'r');

m = 1;                     % loop counter variable
num_proteins = 8;          % number of proteins in the Western Blot data file.
num_protocols_protein = 6; % x2 because avg + std.

protein_expt = zeros(num_proteins,num_protocols_protein);
while true                % infinite loop until end of all files is reached
   data_protein = textscan(fid_protein, data_format, num_protocols_protein); % (num_protocols-1) data columns 
   if isempty(data_protein{1}(:)) 
       break
   end
      
   if ~isempty(data_protein{1}(:)) 
      for i=1:num_protocols_protein
        protein_expt(m,i) = data_protein{1}(i);            % Protein file goes as: ABT-199, Tedizolid, Combination.
      end
   end                 
   m = m + 1;
end

fclose(fid_protein);

%==========================================================================
% 1. Drug and transcription factor parameter values
%==========================================================================
TFPars = [1.159914, 0.027393, 0.972924, 0.106661, 11.425474,...
          2.146729, 3.284086, 4.849403, 2.669784, 3.7034320,...
          1.682787, 1.568424, 6.943936, 4.620715, 0.0001039,...
          2.423171, 208.89961, 0.863131, 2.285103, 0.499996,...
          0.418970, 382.49187, 1.288237, 20.22993, 5.736299,...
          0.920794, 1.5560653, 10.74388, 8.771861, 0.634600,...
          2.107549, 113.76057];     
      
%--------------------------------------------------------------------------
%(2) Other parameter values:
%--------------------------------------------------------------------------
dp_myc_max_tilde=log(2)/15*t_scale; dp_myc_min_tilde=log(2)/180*t_scale;
dp_chop_max_tilde=log(2)/(30)*t_scale; dp_chop_min_tilde=log(2)/(4*60)*t_scale;
%--------------------------------------------------------------------------
dA_tilde = TFPars(1); dE_ABT = TFPars(2); %Already scaled by t_Scale
dT_tilde = TFPars(3); dE_Ted = TFPars(4);  %Already scaled by t_Scale
kA_uptake = dE_ABT*(dE_ABT/dA_tilde)^(dA_tilde/(dE_ABT-dA_tilde)); 
kT_uptake = dE_Ted*(dE_Ted/dT_tilde)^(dT_tilde/(dE_Ted-dT_tilde)); 
%--------------------------------------------------------------------------
kmax_myc_tilde=TFPars(5);
Kf1_myc_A_tilde=TFPars(6); nf1_myc_A=TFPars(7);
Kf1_myc_T_tilde=TFPars(8); nf1_myc_T=TFPars(9);
Kg1_myc_A_tilde=TFPars(10); ng1_myc_A=TFPars(11);
Kg1_myc_T_tilde=TFPars(12); ng1_myc_T=TFPars(13);
dp_myc_0_tilde=TFPars(14); kb_myc_tilde=dp_myc_0_tilde;
%--------------------------------------------------------------------------
kmax_chop_tilde=TFPars(22); 
Kf1_chop_A_tilde=TFPars(23); nf1_chop_A=TFPars(24);
Kg1_chop_A_tilde=TFPars(25); ng1_chop_A=TFPars(26);
Kf1_chop_T_tilde=TFPars(27); nf1_chop_T=TFPars(28);
Kg1_chop_T_tilde=TFPars(29); ng1_chop_T=TFPars(30);
dp_chop_0_tilde=TFPars(31); kb_chop_tilde=dp_chop_0_tilde;
%--------------------------------------------------------------------------
%Weights for combination treatment: (value from 0 to 1)
w1_myc=TFPars(15); w2_myc=TFPars(16); w3_myc=1; w4_myc=1; 
w5_myc=1; w6_myc=1; w7_myc=1; w8_myc=1;
w1_chop=1; w2_chop=1; w3_chop=1; w4_chop=1; 
w5_chop=1; w6_chop=1; w7_chop=1; w8_chop=1;
%--------------------------------------------------------------------------
q_myc=TFPars(17); %cooperativity for cMyc (for production) 
q_chop=TFPars(32); %cooperativity for Chop (for production)
%--------------------------------------------------------------------------
lags = TFPars(18:21); 

%--------------------------------------------------------------------------
% (2) Transcription factor parameter array:
%--------------------------------------------------------------------------
TF_params = [kA_uptake, dA_tilde, dE_ABT,...
             kT_uptake, dT_tilde, dE_Ted,...
             kb_myc_tilde, kmax_myc_tilde,...
             Kf1_myc_A_tilde, nf1_myc_A,...
             Kf1_myc_T_tilde, nf1_myc_T,...
             Kg1_myc_A_tilde, ng1_myc_A,...
             Kg1_myc_T_tilde, ng1_myc_T,...
             dp_myc_max_tilde, dp_myc_min_tilde, dp_myc_0_tilde,...
             kb_chop_tilde, kmax_chop_tilde,... 
             Kf1_chop_A_tilde, nf1_chop_A,...
             Kf1_chop_T_tilde, nf1_chop_T,...   
             Kg1_chop_A_tilde, ng1_chop_A,...
             Kg1_chop_T_tilde, ng1_chop_T,...
             dp_chop_max_tilde, dp_chop_min_tilde, dp_chop_0_tilde,...
             w1_myc, w2_myc, w3_myc, w4_myc,...
             w5_myc, w6_myc, w7_myc, w8_myc,...
             w1_chop, w2_chop, w3_chop, w4_chop,...
             w5_chop, w6_chop, w7_chop, w8_chop,...
             q_myc, q_chop];

%==========================================================================
% (3) Bcl-2 family protein parameter values:
%==========================================================================
% (a) Caspase cleavage parameters:
%--------------------------------------------------------------------------
kc_tilde=myPars(1); kcc_tilde=myPars(2); 
%--------------------------------------------------------------------------
% (b) Target steady-state untreated protein levels: (in nM)
%--------------------------------------------------------------------------
mcl1_0=myPars(3); bcl2_0=myPars(4); 
bim_0=myPars(5); bax_0=myPars(6); bak_0=myPars(7);
casp3_0=myPars(8);
%--------------------------------------------------------------------------
% (c) Background protein production rates:
%--------------------------------------------------------------------------
kb_mcl1_tilde=myPars(9)/mcl1_0; kb_bcl2_tilde=myPars(10)/bcl2_0; 
kb_bim_tilde=myPars(11)/bim_0; kb_bax_tilde=myPars(12)/bax_0; kb_bak_tilde=myPars(13)/bak_0;
%--------------------------------------------------------------------------
% (d) Protein production half saturation constants:
%--------------------------------------------------------------------------
K_myc_mcl1=myPars(14); K_myc_bcl2=myPars(15); 
K_myc_bim=myPars(16); K_myc_bax=myPars(17); K_myc_bak=myPars(18);
K_chop_bcl2=myPars(19);
K_chop_bim=myPars(20); K_chop_bax=myPars(21); K_chop_bak=myPars(22);
%--------------------------------------------------------------------------
% (e) Protein decay rates: 
%--------------------------------------------------------------------------
%https://cancerres.aacrjournals.org/content/canres/suppl/2013/11/19/73.2.519.DC1/methfiglegtab1-9.pdf
% Free proteins:
dp_mcl1_tilde=log(2)/45*t_scale; dp_bcl2_tilde=log(2)/300*t_scale; 
dp_bim_tilde=log(2)/240*t_scale; 
dp_bax_tilde=log(2)/(24*60)*t_scale; dp_bak_tilde=log(2)/(24*60)*t_scale;
%google:
dp_casp3_tilde=log(2)/(8*60)*t_scale; kb_casp3_tilde=dp_casp3_tilde;
% Protein complexes: 
dp_bcl2_bim_tilde=log(2)/75*t_scale; dp_bclxl_bim_tilde=log(2)/75*t_scale; 
dp_mcl1_bim_tilde=log(2)/150*t_scale;
dp_bax_bim_tilde=log(2)/(24*60)*t_scale; dp_bak_bim_tilde=log(2)/(24*60)*t_scale;
dp_bcl2_bax_tilde=log(2)/(24*60)*t_scale; dp_bcl2_bak_tilde=log(2)/(24*60)*t_scale;
dp_mcl1_bax_tilde=log(2)/(24*60)*t_scale; dp_mcl1_bak_tilde=log(2)/(24*60)*t_scale;
%--------------------------------------------------------------------------
% (f) Binding and dissociation constants for protein complexs
%--------------------------------------------------------------------------
%https://cancerres.aacrjournals.org/content/canres/suppl/2013/11/19/73.2.519.DC1/methfiglegtab1-9.pdf
kf_mcl1_bim_tilde = 1.3e-3*60*t_scale; kr_mcl1_bim_tilde = 2.6e-4*60*t_scale;
kf_mcl1_bax_tilde = 2.6e-8*60*t_scale; kr_mcl1_bax_tilde = 2.6e-4*60*t_scale;
kf_mcl1_bak_tilde = 3.25e-5*60*t_scale; kr_mcl1_bak_tilde = 2.6e-4*60*t_scale;
kf_bcl2_bim_tilde = 3.0e-5*60*t_scale; kr_bcl2_bim_tilde = 1.4e-4*60*t_scale;
kf_bcl2_bax_tilde = 9.33e-6*60*t_scale; kr_bcl2_bax_tilde = 1.4e-4*60*t_scale;
kf_bcl2_bak_tilde = 1.4e-8*60*t_scale; kr_bcl2_bak_tilde = 1.4e-4*60*t_scale;
kf_bax_bim_tilde = 2.57e-6*60*t_scale; kr_bax_bim_tilde = 2.57e-4*60*t_scale; 
kappa_bax_bim_tilde = 1.16e-1*60*t_scale;
kf_bak_bim_tilde = 2.57e-6*60*t_scale; kr_bak_bim_tilde = 2.57e-4*60*t_scale; 
kappa_bak_bim_tilde = 1.16e-1*60*t_scale;

%--------------------------------------------------------------------------
% (g) Drug-related parameters:
%--------------------------------------------------------------------------
%From dissociation constant of Bcl-XL with ABT-737 (was the lowest)
%https://cancerres.aacrjournals.org/content/canres/suppl/2013/11/19/73.2.519.DC1/methfiglegtab1-9.pdf
%pg. 18
kf_abt_bcl2_tilde=3.85e-4*60*t_scale; kr_abt_bcl2_tilde=1.93e-4*60*t_scale; 
kted_tilde=0.0; kted_2_tilde=0.0;
K_ABT_bax_tilde=myPars(23);
%--------------------------------------------------------------------------
% (h) Cell viability parameters:
%--------------------------------------------------------------------------
lambda_bcl2_tilde=myPars(24); lambda_myc_tilde=myPars(25); 
lambda_0_tilde=myPars(26); delta_casp3_tilde=myPars(27); 
Bcl2_thresh=myPars(28); q_Ted_tilde=myPars(46); 
n_q_Ted=myPars(47);
Myc_thresh=myPars(39); Casp3_thresh=myPars(40);

%Test the parameter values for biological relevancy:
if (lambda_bcl2_tilde > log(2)/(8*60)*t_scale) 
   rSq = rSq + exp(lambda_bcl2_tilde - log(2)/(8*60)*t_scale);
elseif (lambda_bcl2_tilde < log(2)/(28*60)*t_scale)
   rSq = rSq + exp(log(2)/(28*60)*t_scale - lambda_bcl2_tilde);
end
%-----------
if (lambda_0_tilde > log(2)/(8*60)*t_scale) 
   rSq = rSq + exp(lambda_0_tilde - log(2)/(8*60)*t_scale);
elseif (lambda_0_tilde < log(2)/(28*60)*t_scale)
   rSq = rSq + exp(log(2)/(28*60)*t_scale - lambda_0_tilde);
end
%-----------
if (lambda_myc_tilde > log(2)/(4*60)*t_scale) 
   rSq = rSq + exp(lambda_myc_tilde - log(2)/(4*60)*t_scale);
elseif (lambda_myc_tilde < log(2)/(2*24*60)*t_scale)
   rSq = rSq + exp(log(2)/(2*24*60)*t_scale - lambda_myc_tilde);
end
%-----------
if (delta_casp3_tilde > log(2)/(4*60)*t_scale) 
   rSq = rSq + exp(delta_casp3_tilde - log(2)/(4*60)*t_scale);
elseif (delta_casp3_tilde < log(2)/(2*24*60)*t_scale)
   rSq = rSq + exp(log(2)/(2*24*60)*t_scale - delta_casp3_tilde);
end

%--------------------------------------------------------------------------
n1=myPars(29); n2=myPars(30); n3=myPars(31); n4=myPars(32); n5=myPars(33);
n6=myPars(34); n7=myPars(35); n8=myPars(36); n9=myPars(37); n10=myPars(38);
%--------------------------------------------------------------------------
K_chop_mcl1=myPars(41); n11=myPars(42); n12=myPars(43);
nv_1=myPars(44); nv_2=myPars(45);

n13=myPars(48); K_myc_bim_2=myPars(49);
n14=myPars(50); K_myc_bcl2_2=myPars(51);

%--------------------------------------------------------------------------
% (i) Calculated protein production rates:
%--------------------------------------------------------------------------
k_mcl1_tilde = (dp_mcl1_tilde - kb_mcl1_tilde)*(1.0+K_myc_mcl1^n1)*(1.0+1.0/K_chop_mcl1^n11);        
if k_mcl1_tilde < 0
   disp('Error: Negative Mcl-1 production rate')
   rSq = rSq + 1.0e8 + exp(-k_mcl1_tilde);
   return
end   
%--------------------------------------------------
k_bcl2_tilde = (dp_bcl2_tilde - kb_bcl2_tilde)*(1.0+1.0/K_myc_bcl2^n2)*(1.0+1.0/K_chop_bcl2^n3)*(1.0 + K_myc_bcl2_2^n14);
if k_bcl2_tilde < 0
   disp('Error: Negative Bcl-2 production rate')
   rSq = rSq + 1.0e8 + exp(-k_bcl2_tilde);
   return
end   
%--------------------------------------------------         
k_bim_tilde = (dp_bim_tilde - kb_bim_tilde)/( (1.0 - 1.0/(K_myc_bim^n5 + 1.0))*1.0/(K_chop_bim^n4 + 1.0)/(1.0 + (1.0/K_myc_bim_2)^n13)...
    + 1.0/(K_myc_bim^n5 + 1.0) );    
if k_bim_tilde < 0
   disp('Error: Negative Bim production rate')
   rSq = rSq + 1.0e8 + exp(-k_bim_tilde);
   return
end   
%--------------------------------------------------                                               
k_bax_tilde = (dp_bax_tilde - kb_bax_tilde)*(1.0+1.0/K_chop_bax^n6 + 1.0/K_myc_bax^n7 + 1.0/(K_chop_bax^n6*K_myc_bax^n7))...
                                           /(1.0/K_chop_bax^n6 + 1.0/K_myc_bax^n7 + 1.0/(K_chop_bax^n6*K_myc_bax^n7));
if k_bax_tilde < 0
%   disp('Error: Negative Bax production rate')
   rSq = rSq + 1.0e8 + exp(-k_bax_tilde);
   return
end   
%--------------------------------------------------        
k_bak_tilde = (dp_bak_tilde - kb_bak_tilde)*(1.0+1.0/K_chop_bak^n8+1.0/K_myc_bak^n9+1.0/(K_chop_bak^n8*K_myc_bak^n9))...
                                           /(1.0/K_chop_bak^n8+1.0/K_myc_bak^n9+1.0/(K_chop_bak^n8*K_myc_bak^n9));
if k_bak_tilde < 0
%   disp('Error: Negative Bak production rate')
   rSq = rSq + 1.0e8 + exp(-k_bak_tilde);
   return
end           

%==========================================================================
% (1) Equilibrium simulations:
%==========================================================================
%--------------------------------------------------------------------------
% (ii) Equilibrium initial conditions protein array:
%--------------------------------------------------------------------------
% MCL-1/mcl_0, BCL-2/bcl2_0,  
% BIM/bim_0, BAX/bax_0, BAK/bak_0, BAX*/bax_0, BAK*/bak_0, 
% CASP-3/casp3_0, CASP-3*/casp3_0, MCL-1-/mcl1_0, 
% MCL-1:BIM/mcl1_0, BCL-2:BIM/bcl2_0, 
% BIM:BAX/bax_0, BIM:BAK/bak_0, 
% MCL-1:BAX*/mcl1_0, MCL-1:BAK*/mcl1_0, BCL-2:BAX*/bcl2_0, BCL-2:BAK*/bcl2_0,  
% BCL-2-/bcl2_0
IC_protein_0 = [1.0, 1.0,...
                1.0, 1.0, 1.0, 0, 0,...
                1.0, 0, 0,... 
                0, 0,... 
                0, 0,... 
                0, 0, 0, 0,... 
                0];
IC_protein_0 = IC_protein_0.';

%--------------------------------------------------------------------------
% (iii) Equilibrium parameter array:
%--------------------------------------------------------------------------
DE_params_equil = [kc_tilde, kcc_tilde,... 
                   k_mcl1_tilde, k_bcl2_tilde,... 
                   k_bim_tilde, k_bax_tilde, k_bak_tilde,...
                   kb_mcl1_tilde, kb_bcl2_tilde,...  
                   kb_bim_tilde, kb_bax_tilde, kb_bak_tilde,...
                   kb_casp3_tilde,... 
                   K_myc_mcl1, K_myc_bcl2,... 
                   K_myc_bim, K_myc_bax, K_myc_bak,...
                   K_chop_bcl2,...
                   K_chop_bim, K_chop_bax, K_chop_bak,...
                   dp_mcl1_tilde, dp_bcl2_tilde,... 
                   dp_bim_tilde, dp_bax_tilde, dp_bak_tilde,...
                   dp_bcl2_bim_tilde, dp_mcl1_bim_tilde,...
                   dp_bax_bim_tilde, dp_bak_bim_tilde,...
                   dp_bcl2_bax_tilde, dp_bcl2_bak_tilde,...
                   dp_mcl1_bax_tilde, dp_mcl1_bak_tilde,...
                   kf_mcl1_bim_tilde, kr_mcl1_bim_tilde,...
                   kf_mcl1_bax_tilde, kr_mcl1_bax_tilde,...
                   kf_mcl1_bak_tilde, kr_mcl1_bak_tilde,...
                   kf_bcl2_bim_tilde, kr_bcl2_bim_tilde,...
                   kf_bcl2_bax_tilde, kr_bcl2_bax_tilde,...
                   kf_bcl2_bak_tilde, kr_bcl2_bak_tilde,...
                   kf_bax_bim_tilde, kr_bax_bim_tilde,... 
                   kappa_bax_bim_tilde,...
                   kf_bak_bim_tilde, kr_bak_bim_tilde,... 
                   kappa_bak_bim_tilde,...
                   n1, n2, n3, n4, n5, n6, n7, n8, n9,...
                   K_chop_mcl1, n11,...
                   n13, K_myc_bim_2,...
                   n14, K_myc_bcl2_2];               
               
%--------------------------------------------------------------------------
% (iv) Equilibrium transcription factor protein array:
%--------------------------------------------------------------------------
TF_0 = zeros(2,1);
TF_0(1) = 1.0;    %protein/protein_0
TF_0(2) = 1.0; 

P_0 = [mcl1_0; bcl2_0;...  
       bim_0; bax_0; bak_0; casp3_0];                             

%--------------------------------------------------------------------------
% (v) Solve for equilibrium untreated protein levels:
%--------------------------------------------------------------------------
% An event function is included here to ensure the system reaches an 
% equilibrium state. Equilibrium is defined to be a steady state in which 
% all protein levels are constant, with activated caspase-3 levels very low, 
% e.g. < 10% of the overall concentration of caspase-3. This way, the 
% equilibrium state in untreated cells reflects a growing cellular 
% population (and the cells won't spontaneously undergo apoptosis). 
tmin = 0.0;
tmax_equil = 1000000/t_scale;
num_equil_points = 2000000;
tspan_equil=linspace(tmin,tmax_equil,num_equil_points);

XX = [mcl1_0; bcl2_0; bim_0; bax_0; bak_0;...
      bax_0; bak_0; casp3_0; casp3_0;... 
      mcl1_0; mcl1_0; bcl2_0;... 
      bax_0; bak_0;...       
      mcl1_0; mcl1_0;... 
      bcl2_0; bcl2_0;... 
      bcl2_0]; 

% Pass the parameters to the event function first:
ParamEventHdl = @(t_equil, x_equil)EquilibriumEventFcn_no_feedback(t_equil,x_equil,DE_params_equil,TF_0,P_0,t_scale,XX); 
opts_equil = odeset('RelTol',1e-10,'AbsTol',1e-12,'Events',ParamEventHdl); % Add the events function here. 
%--------------------------------------------------------------------------
% Run the equilibration simulation:
%--------------------------------------------------------------------------
% Using a stiff solver to make sure the results are equivalent (stiff 
% solver is faster for our equations):
[t_equil,x_equil] = ode15s(@(t_equil,x_equil)...
       integratingfunction_equilibrium_no_feedback(t_equil,x_equil,DE_params_equil,TF_0,P_0),tspan_equil,IC_protein_0,opts_equil);
%--------------------------------------------------------------------------
% Confirm that the event function was called: 
%--------------------------------------------------------------------------
% Ensure that an appropriate steady state solution was found. This will not 
% be the case (i.e. an event will not be found) max(t_equil) = tmax_equil. 
% Thus a control statement is included which prints a warning message and 
% if an appropriate steady state was not found for the given parameter 
% values, so that new parameter values can be tried. 
max_ind = max(size(t_equil));
Calculated_protein_SS = x_equil(max_ind,:);

if max(t_equil)==tmax_equil                        
%   msg = 'Error: no steady state found in alloted time for initial protein levels from simulations. Max(t_equil) = ';   
%   disp(msg)
   proceed = 0;
   rSq = 1.0e9;                     
else  % Steady state was found in alloted time.
   acasp_frac_equil = Calculated_protein_SS(9)/(Calculated_protein_SS(8)+Calculated_protein_SS(9));     % Casp-3*/total Casp-3.
   if acasp_frac_equil >= 0.15
%      msg = 'Warning: Parameter values lead to apoptotic steady-state. acasp_frac = ';
%      disp(msg)
      proceed = 0;
      rSq = 1.0e9 + exp(acasp_frac_equil);
   else
      proceed = 1;
   end
end

% If a non-apoptotic steady state was found, simulate the untreated system
% for 5 days and calculate the corresponding live cell numbers and cell
% viabilities:
%--------------------------------------------------------------------------
if proceed == 1
%--------------------------------------------------------------------------
   Bcl2_0 = Calculated_protein_SS(2);
   IC_N_0 = [50000, 813];   %Initial viability set to average of the media viabilities = 98.4
   IC_N_0 = IC_N_0.';

   DE_params_N = [Bcl2_0, TF_0(1), acasp_frac_equil,...
                  lambda_bcl2_tilde, lambda_myc_tilde,... 
                  lambda_0_tilde, Bcl2_thresh, Myc_thresh,... 
                  delta_casp3_tilde, Casp3_thresh, n12,...
                  nv_1, nv_2];
              
   % Run the simulation:                             
   tmin = 0.0;
   tmax = 5.0*24*60/t_scale;
   num_points = 25000;
   tspan=linspace(tmin,tmax,num_points);

   opts = odeset('RelTol',1e-10,'AbsTol',1e-12);  

   [t_N,x_N] = ode15s(@(t_N,x_N)...
       integratingfunction_untreated_N(t_N,x_N,DE_params_N),tspan,IC_N_0,opts);

   if((max(size(t_N))) < num_points)
     rSq = 1.0e9;
   else  
   
   % Compare cell counts and cell viability with experimental data: 
   days = 24.0*60.0/t_scale*[1, 2, 3, 4, 5];  
   day_ind = zeros(1,5);
   
%--------------------------------------------------------------------------
   % Experimental ratio of N on days 1-5 to N on day 0:
   Expt_N_0 = [2.853333333, 5.786666667, 39.0, 90.0, 180.0];
   Std_N_0 = [0.435583899, 1.892758129, 5.196152423, 0.0, 0.0];   
%-------------------------------------------------------------------------- 
   % Experimental viability data on days 1-5
   Expt_V_0 = [97.76666667, 98.33333333, 98.7, 98.63333333, 98.8];
   Std_V_0 = [0.56862407, 0.4163332, 0.3, 0.450924975, 0.264575131];

   N_Sim_0 = zeros(1,5);
   V_Sim_0 = zeros(1,5);
   N_sq_0 = 0;
   V_sq_0 = 0;

   for i=1:5
       [dist, ind] = min(abs(t_N-days(i)));
       day_ind(i) = ind;
       N_Sim_0(i) = x_N(day_ind(i),1)/x_N(1,1);
       L_N = Expt_N_0(i) - 1.5*Std_N_0(i); 
       U_N = Expt_N_0(i) + 1.5*Std_N_0(i);
      
       if (N_Sim_0(i)>=L_N)
           if (N_Sim_0(i)<=U_N)
               N_sq = 0.0;
           else 
               N_sq = abs(N_Sim_0(i)-U_N)/U_N;     
           end    
       else
           N_sq = abs(N_Sim_0(i)-L_N)/L_N;
       end                              
       N_sq_0 = N_sq_0 + N_sq;
%---------------------------------       
       V_Sim_0(i) = 100.0*x_N(day_ind(i),1)/(x_N(day_ind(i),1) + x_N(day_ind(i),2));
       L_V = Expt_V_0(i) - 1.5*Std_V_0(i);
       U_V = Expt_V_0(i) + 1.5*Std_V_0(i);

       if (V_Sim_0(i)>=L_V)
           if (V_Sim_0(i)<=U_V)
               V_sq = 0.0;
           else 
               V_sq = abs(V_Sim_0(i)-U_V)/U_V;     
           end    
       else
           V_sq = abs(V_Sim_0(i)-L_V)/L_V;
       end                              
       V_sq_0 = V_sq_0 + V_sq;
   end
   
   proceed_2 = 1;
   end   
end   

% Now proceed to ABT-199 and Tedizolid treatment simulations:
%--------------------------------------------------------------------------
if proceed_2 == 1
%--------------------------------------------------------------------------
% Determine total protein levels in untreated cells. The total protein 
% array has the following order: MCL-1/mcl1_0; BCL-2/bcl2_0; BIM/bim_0; 
% BAX/bax_0; BAK/bak_0; aCASP-3/casp3_0; c-Myc/cMyc_0; Chop/Chop_0.
   Total_protein_equil = zeros(num_proteins,1);  

   %Mcl-1: 
   Total_protein_equil(1) = Calculated_protein_SS(1) + Calculated_protein_SS(10)...
                          + Calculated_protein_SS(11)...
                          + Calculated_protein_SS(15) + Calculated_protein_SS(16);   
   %Bcl-2:
   Total_protein_equil(2) = Calculated_protein_SS(2) + Calculated_protein_SS(19)...
                          + Calculated_protein_SS(12)...
                          + Calculated_protein_SS(17) + Calculated_protein_SS(18);                      
   %Bim:
   Total_protein_equil(3) = Calculated_protein_SS(3)...
                          + Calculated_protein_SS(11)*mcl1_0/bim_0 + Calculated_protein_SS(12)*bcl2_0/bim_0...
                          + Calculated_protein_SS(13)*bax_0/bim_0 + Calculated_protein_SS(14)*bak_0/bim_0;
   %Bax:
   Total_protein_equil(4) = Calculated_protein_SS(4) + Calculated_protein_SS(6)...
                          + Calculated_protein_SS(13)...
                          + Calculated_protein_SS(15)*mcl1_0/bax_0 + Calculated_protein_SS(17)*bcl2_0/bax_0;
   %Bak:
   Total_protein_equil(5) = Calculated_protein_SS(5) + Calculated_protein_SS(7)...
                          + Calculated_protein_SS(14)...
                          + Calculated_protein_SS(16)*mcl1_0/bak_0 + Calculated_protein_SS(18)*bcl2_0/bak_0;
   %aCasp-3:                
   Total_protein_equil(6) = Calculated_protein_SS(9);  %Calculated_protein_SS(8) +               
   
   %c-Myc:
   Total_protein_equil(7) = 1.0;
   
   %Chop:
   Total_protein_equil(8) = 1.0;

%==========================================================================
% (2) Protein levels in treated cells:
%==========================================================================
% Treatment parameter array:
Drug_params = [kf_abt_bcl2_tilde, kr_abt_bcl2_tilde,... 
               kted_tilde, kted_2_tilde,...
               K_ABT_bax_tilde];

Via_params = [lambda_bcl2_tilde, lambda_myc_tilde,... 
              lambda_0_tilde, delta_casp3_tilde,...
              Bcl2_thresh,...
              q_Ted_tilde];
                   
DE_params = [DE_params_equil(1:53), TF_params, Drug_params, Via_params,...
             n1, n2, n3, n4, n5, n6, n7, n8, n9, n10,...
             n_q_Ted, Myc_thresh, Casp3_thresh,...
             K_chop_mcl1, n11, n12, nv_1, nv_2,...
             n13, K_myc_bim_2, n14, K_myc_bcl2_2];               
     
%--------------------------------------------------------------------------
% Time span for simulations:
   tmin = 0.0;
   num_points = 25000;
   tmax = 5.0*24*60/t_scale;
   tspan_treatment=linspace(tmin,tmax,num_points);
%--------------------------------------------------------------------------
% (a) ABT-199 400 nM treatment:
%--------------------------------------------------------------------------
% (i) Initial conditions protein array for calculation
%--------------------------------------------------------------------------
% MCL-1/mcl1_0, BCL-2/bcl2_0, 
% BIM/bim_0, BAX/bax_0, BAK/bak_0, BAX*/bax_0, BAK*/bak_0, 
% CASP-3/casp3_0, CASP-3*/casp3_0, MCL-1-/mcl1_0, 
% MCL-1:BIM/mcl1_1, BCL-2:BIM/bcl2_0,  
% BIM:BAX/bax_0, BIM:BAK/bak_0, 
% MCL-1:BAX*/mcl1_0, MCL-1:BAK*/mcl1_0, BCL-2:BAX*/bcl_2, BCL-2:BAK*/bcl2_0, 
% BCL-2-/bcl2_0, 
% c-Myc/cMyc_0, Chop/Chop_0
% extracellular ABT concentration/ABT_0, free intracellular ABT concentration/ABT_0, 
% ABT:Bcl-2 complex/ABT_0,
% live cell number, dead cell number=x(26);
   IC_protein_ABT = [Calculated_protein_SS,...
                     1.0, 1.0,...
                     1.0, 0.0, 0.0,...
                     50000, 813]; % Initial viability set to 98.4 = average of untreated
   opts_ABT = odeset('RelTol',1e-8,'AbsTol',1e-10,'Stats','off');%,'Events',ParamEventHdl); % Add the events function here. 
   P_0_ABT = [mcl1_0; bcl2_0;...  
              bim_0; bax_0; bak_0; casp3_0;...
              ABT_0];   

%--------------------------------------------------------------------------
% (b) Tedizolid 5000 nM treatment:
%--------------------------------------------------------------------------
% (i) Initial conditions protein array for calculation
%--------------------------------------------------------------------------
% MCL-1/mcl1_0, BCL-2/bcl2_0, 
% BIM/bim_0, BAX/bax_0, BAK/bak_0, BAX*/bax_0, BAK*/bak_0, 
% CASP-3/casp3_0, CASP-3*/casp3_0, MCL-1-/mcl1_0, 
% MCL-1:BIM/mcl1_0, BCL-2:BIM/bcl2_0,  
% BIM:BAX/bax_0, BIM:BAK/bak_0, 
% MCL-1:BAX*/mcl1_0, MCL-1:BAK*/mcl1_0, BCL-2:BAX*/bcl2_0, BCL-2:BAK*/bcl2_0, 
% BCL-2-/bcl2_0, 
% c-Myc/cMyc_0, Chop/Chop_0
% external Tedizolid concentration/Ted_0, 
% intracellular Tedizolid concentration/Ted_0, 
% live cell number, dead cell number
   IC_protein_Ted = [Calculated_protein_SS,...
                     1.0, 1.0,...
                     1.0, 0.0,...
                     50000, 813];  %Initial viability set to average of untreated
   opts_Ted = odeset('RelTol',1e-8,'AbsTol',1e-10,'Stats','off');%,'Events',ParamEventHdl); % Add the events function here. 
   P_0_Ted = [mcl1_0; bcl2_0;...  
              bim_0; bax_0; bak_0; casp3_0;...
              Ted_0];         

%--------------------------------------------------------------------------
% (c) Combination ABT-199 and Tedizolid treatment:
%--------------------------------------------------------------------------
% (i) Initial conditions protein array for calculation
%--------------------------------------------------------------------------
% MCL-1/mcl1_0, BCL-2/bcl2_0, 
% BIM/bim_0, BAX/bax_0, BAK/bak_0, BAX*/bax_0, BAK*/bak_0, 
% CASP-3/casp3_0, CASP-3*/casp3_0, MCL-1-/mcl1_0, 
% MCL-1:BIM/mcl1_0, BCL-2:BIM/bcl2_0,  
% BIM:BAX/bax_0, BIM:BAK/bak_0, 
% MCL-1:BAX*/bax_0, MCL-1:BAK*/bak_0, BCL-2:BAX*/bax_0, BCL-2:BAK*/bak_0, 
% BCL-2-/bcl2_0, 
% c-Myc/cMyc_0, Chop/Chop_0
% extracellular ABT concentration/ABT_0, free intracellular ABT concentration/ABT_0, 
% ABT:Bcl-2 complex/ABT_0,
% extracellular Tedizolid concentration/Ted_0, 
% intracellular Tedizolid concentration/Ted_0, 
% live cell number, dead cell number=x(26);
   IC_protein_AT = [Calculated_protein_SS,...
                    1.0, 1.0,...
                    1.0, 0.0, 0.0,...
                    1.0, 0.0,...
                    50000, 813];
   opts_AT = odeset('RelTol',1e-8,'AbsTol',1e-10,'Stats','off');%,'Events',ParamEventHdl); % Add the events function here. 
   P_0_AT = [mcl1_0; bcl2_0;...  
             bim_0; bax_0; bak_0; casp3_0;...
             ABT_0; Ted_0];   

%--------------------------------------------------------------------------
% Run simulations          
%--------------------------------------------------------------------------
   %ABT-199 monotherapy
   sol_ABT = dde15s_new(@integratingfunction_ABT199_no_feedback,lags,IC_protein_ABT,tspan_treatment,opts_ABT,DE_params,P_0_ABT);   
   sol_ABT = deval(sol_ABT, tspan_treatment);       
   t_ABT = tspan_treatment;
   x_ABT = sol_ABT.';

   %Tedizolid monotherapy
   sol_Ted = dde15s_new(@integratingfunction_Tedizolid_no_feedback,lags,IC_protein_Ted,tspan_treatment,opts_Ted,DE_params,P_0_Ted);
   sol_Ted = deval(sol_Ted, tspan_treatment);
   t_Ted = tspan_treatment;
   x_Ted = sol_Ted.';

   %Combination therapy
   sol_AT = dde15s_new(@integratingfunction_combination_ABT_Ted_no_feedback,lags,IC_protein_AT,tspan_treatment,opts_AT,DE_params,P_0_AT);
   sol_AT = deval(sol_AT, tspan_treatment);
   t_AT = tspan_treatment;
   x_AT = sol_AT.';

   Mat_size = [max(size(x_ABT)), max(size(x_Ted)), max(size(x_AT))];

%--------------------------------------------------------------------------   
%--------------------------------------------------------------------------   
% Compare simulations with experimental data:
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
if ( Mat_size(1) < num_points) || ( Mat_size(2) < num_points) || ( Mat_size(3) < num_points)
   rSq = 1.0e9;
%   disp('Problem!')
else       
%   disp('No problem!')
%--------------------------------------------------------------------------
% Protein levels on day 3:
%--------------------------------------------------------------------------
  ABT_protein_3d = x_ABT(day_ind(3),:);
  Ted_protein_3d = x_Ted(day_ind(3),:);    
    
%--------------------------------------------------------------------------
%(1) ABT-199 cell numbers and cell viability on days 1, 2, 3, 4, 5.
%--------------------------------------------------------------------------
   % Experimental ratio of N on days 1-5 to N on day 0:
   Expt_N_A = [2.26, 5.793333333, 30.56, 55.0, 180.0];
   Std_N_A = [0.207846097, 2.715903778, 5.182740588, 8.660254038, 0.0];      
   %-------------------------------------------------------------------------- 
   % Experimental viability data on days 1-5
   Expt_V_A = [98.2, 98.4, 97.73333333, 95.56666667, 95.1];   
   Std_V_A = [0.360555128, 0.360555128, 0.37859389, 0.808290377, 1.228820573];
   
   N_Sim_A = zeros(1,5);
   V_Sim_A = zeros(1,5);
   N_sq_A = 0;
   V_sq_A = 0;
   for i=1:5
       [dist, ind] = min(abs(t_ABT-days(i)));
       day_ind(i) = ind;
       N_Sim_A(i) = x_ABT(day_ind(i),25)/x_ABT(1,25);
       L_N = Expt_N_A(i) - 1.0*Std_N_A(i); 
       U_N = Expt_N_A(i) + 1.0*Std_N_A(i);
      
       if (N_Sim_A(i)>=L_N)
           if (N_Sim_A(i)<=U_N)
               N_sq = 0.0;
           else 
               N_sq = abs(N_Sim_A(i)-U_N)/U_N;     
           end    
       else
           N_sq = abs(N_Sim_A(i)-L_N)/L_N;
       end                              
       N_sq_A = N_sq_A + N_sq;
%---------------------------------
       V_Sim_A(i) = 100.0*x_ABT(day_ind(i),25)/(x_ABT(day_ind(i),25) + x_ABT(day_ind(i),26));
       L_V = Expt_V_A(i) - 1.25*Std_V_A(i);
       U_V = Expt_V_A(i) + 1.25*Std_V_A(i);

       if (V_Sim_A(i)>=L_V)
           if (V_Sim_A(i)<=U_V)
               V_sq = 0.0;
           else 
               V_sq = abs(V_Sim_A(i)-U_V)/U_V;     
           end    
       else
           V_sq = abs(V_Sim_A(i)-L_V)/L_V;
       end                              
       V_sq_A = V_sq_A + V_sq;
   end
   
%--------------------------------------------------------------------------
%(2) Tedizolid cell numbers and cell viability on days 1, 2, 3, 4, 5.
%--------------------------------------------------------------------------
% Ratio of N on days 1-5 to N on day 0:
   Expt_N_T = [2.113333333, 4.706666667, 22.08, 25.24, 75.0];
   Std_N_T = [0.620429958, 1.18816385, 3.862434465, 4.994557037, 25.98076211];   
%----------------------------------
   Expt_V_T = [94.73333333, 94.5, 94.3, 93.5, 94.03333333]; 
   Std_V_T = [0.901849951, 0.458257569, 0.458257569, 0.264575131, 0.838649708];
   
   N_Sim_T = zeros(1,5);
   V_Sim_T = zeros(1,5);
   N_sq_T = 0;
   V_sq_T = 0;
   for i=1:5
       [dist, ind] = min(abs(t_Ted-days(i)));
       day_ind(i) = ind;
       N_Sim_T(i) = x_Ted(day_ind(i),24)/x_Ted(1,24);
       L_N = Expt_N_T(i) - 1.0*Std_N_T(i); 
       U_N = Expt_N_T(i) + 1.0*Std_N_T(i);

       if (N_Sim_T(i)>=L_N)
           if (N_Sim_T(i)<=U_N)
               N_sq = 0.0;
           else 
               N_sq = abs(N_Sim_T(i)-U_N)/U_N;     
           end    
       else
           N_sq = abs(N_Sim_T(i)-L_N)/L_N;
       end
       if i~=3                              
          N_sq_T = N_sq_T + N_sq;
       end   
%---------------------------------
       V_Sim_T(i) = 100.0*x_Ted(day_ind(i),24)/(x_Ted(day_ind(i),24) + x_Ted(day_ind(i),25));
       L_V = Expt_V_T(i) - 1.25*Std_V_T(i);
       U_V = Expt_V_T(i) + 1.25*Std_V_T(i);

       if (V_Sim_T(i)>=L_V)
           if (V_Sim_T(i)<=U_V)
               V_sq = 0.0;
           else 
               V_sq = abs(V_Sim_T(i)-U_V)/U_V;     
           end    
       else
           V_sq = abs(V_Sim_T(i)-L_V)/L_V;
       end                              
       V_sq_T = V_sq_T + V_sq;
   end
   
%==========================================================================
% (3) Calculate rSq protein value for both treatment protocols:
%==========================================================================
%--------------------------------------------------------------------------
% (i) Calculate the total steady-state protein levels: 
%--------------------------------------------------------------------------
   Proteins_treatment = zeros(3,21);
   Proteins_treatment(1,:) = ABT_protein_3d(1:21);
   Proteins_treatment(2,:) = Ted_protein_3d(1:21);
   
% The total protein array has the following order: 
% MCL-1/mcl1_0; BCL-2/bcl2_0; BIM/bim_0; BAX/bax_0; BAK/bak_0; 
% aCASP-3/casp3_0; c-Myc/cMyc_0; Chop/Chop_0.
   Total_protein = zeros(num_proteins, num_protocols);             % Total protein levels on day 3 after treatment. 
   Total_protein(:,1) = Total_protein_equil;                       % Untreated equilibrium levels calculated previously.
   Protein_ratio = zeros(num_proteins, num_protocols-1); 
   
   Protein_comparison = zeros(num_proteins, num_protocols-1);      % Protein ratios in the same order as Western Blot experimental data

%--------------------------------------------------------------------------
% Calculate total protein levels and rSq value for each treatment 
% condition:
%--------------------------------------------------------------------------
  for i=2:3
       Protein_ratio_rSq = 0.0;
       
       if i==2
          Bcl2_complex = ABT_protein_3d(24);     % ABT:Bcl-2 complex level.
       else
          Bcl2_complex = 0.0;
       end   
      
       %Mcl-1: 
       Total_protein(1,i) = Proteins_treatment(i-1,1)...
                          + Proteins_treatment(i-1,10)...
                          + Proteins_treatment(i-1,11)...
                          + Proteins_treatment(i-1,15)...
                          + Proteins_treatment(i-1,16); 
       %Bcl-2:
       Total_protein(2,i) = Proteins_treatment(i-1,2)...
                          + Proteins_treatment(i-1,19)...
                          + Proteins_treatment(i-1,12)...
                          + Proteins_treatment(i-1,17)...
                          + Proteins_treatment(i-1,18)... 
                          + Bcl2_complex*ABT_0/bcl2_0;     % ABT:Bcl-2 complex level.
       %Bim:
       Total_protein(3,i) = Proteins_treatment(i-1,3)...
                          + Proteins_treatment(i-1,11)*mcl1_0/bim_0...
                          + Proteins_treatment(i-1,12)*bcl2_0/bim_0...
                          + Proteins_treatment(i-1,13)*bax_0/bim_0...
                          + Proteins_treatment(i-1,14)*bak_0/bim_0; 
       %Bax: 
       Total_protein(4,i) = Proteins_treatment(i-1,4)...
                          + Proteins_treatment(i-1,6)...
                          + Proteins_treatment(i-1,13)...
                          + Proteins_treatment(i-1,15)*mcl1_0/bax_0...
                          + Proteins_treatment(i-1,17)*bcl2_0/bax_0; 
       %Bak:
       Total_protein(5,i) = Proteins_treatment(i-1,5)...
                          + Proteins_treatment(i-1,7)...
                          + Proteins_treatment(i-1,14)...
                          + Proteins_treatment(i-1,16)*mcl1_0/bak_0...
                          + Proteins_treatment(i-1,18)*bcl2_0/bak_0; 
       %aCasp-3:                
       Total_protein(6,i) = Proteins_treatment(i-1,9);     % Is actually active casp-3 level. 

       %c-Myc:
       Total_protein(7,i) = Proteins_treatment(i-1,20);
       
       %Chop:
       Total_protein(8,i) = Proteins_treatment(i-1,21);
       
       %disp(Total_protein)
       
       Protein_ratio(:,i-1) = Total_protein(:,i)./Total_protein(:,1);
  

% Cleaved Caspase-3 levels should be roughly the same as the untreated 
% levels for ABT-199 and Tedizolid treatment. So we check to make sure that 
% the ratio of cleaved Casp-3 during treatment is between 0.5 and 2.5. If 
% it falls outside this range, then we will impose a large penalty on the 
% rSq value.
      for k=1:num_proteins
          if k==6  % cleaved Casp-3
             acasp3_ratio = Protein_ratio(6,i-1); 
             if i==2 
                 if (acasp3_ratio>=1.0)
                     if (acasp3_ratio<=2.5)
                         diff_sq = 0.0;
                     else 
                         diff_sq = abs(acasp3_ratio-2.5)/2.5;     
                     end    
                 else
                     diff_sq = abs(1.0-acasp3_ratio)/1.0;
                 end                 
             else    
                 if (acasp3_ratio>=0.5)
                     if (acasp3_ratio<=2.5)
                         diff_sq = 0.0;
                     else 
                         diff_sq = abs(acasp3_ratio-2.5)/2.5;     
                     end    
                 else
                     diff_sq = abs(0.5-acasp3_ratio)/0.5;
                 end                 
             end
             Protein_comparison(k,i-1) = acasp3_ratio;
          else
             L_P = protein_expt(k,2*(i-1)-1) - protein_expt(k,2*(i-1)); 
             U_P = protein_expt(k,2*(i-1)-1) + protein_expt(k,2*(i-1));             
             if (Protein_ratio(k,i-1)>=L_P)
                 if (Protein_ratio(k,i-1)<=U_P)
                     diff_sq = 0.0;
                 else 
                     diff_sq = abs(Protein_ratio(k,i-1)-U_P)/U_P;        
                 end    
             else
                 diff_sq = abs(Protein_ratio(k,i-1)-L_P)/L_P;        %Can change this to actual rsq value if you want (e.g. (Protein_sim(i,j)-L_P)^2 )
             end                              
             Protein_comparison(k,i-1) = Protein_ratio(k,i-1);                
          end
          Protein_ratio_rSq = Protein_ratio_rSq + diff_sq;
      end      
      rSq = rSq + Protein_ratio_rSq;   
  end
  
end

%==========================================================================
% Now do combination treatment:
%==========================================================================
%--------------------------------------------------------------------------
% Protein levels on day 3:
%--------------------------------------------------------------------------
  [dist, d3_ind] = min(abs(t_AT-3*24*60/t_scale));
   AT_protein_3d = x_AT(d3_ind,:);
%--------------------------------------------------------------------------
%(1) Cell numbers and cell viability on days 1, 2, 3, 4, 5.
%--------------------------------------------------------------------------
   % Compare cell counts and cell viability with experimental data: 
   
%--------------------------------------------------------------------------
   % Ratio of N on days 1-5 to N on day 0:
   Expt_N_AT = [2.366666667, 2.953333333, 9.48, 7.56, 1.88];
   Std_N_AT = [0.574224114, 1.372054421, 1.383907511, 0.634980315, 0.069282032];   
%-------------------------------------------------------------------------- 
   % Experimental viability data on days 1-5
   Expt_V_AT = [92.26666667, 80.73333333, 64.03333333, 35.5, 15.4];
   Std_V_AT = [1.171893055, 1.942506971, 1.601041328, 1.5, 0.871779789];
   
   N_Sim_AT = zeros(1,5);
   V_Sim_AT = zeros(1,5);
   N_sq_AT = 0;
   V_sq_AT = 0;
   for i=1:5
       [dist, ind] = min(abs(t_AT-days(i)));
       day_ind(i) = ind;
       N_Sim_AT(i) = x_AT(day_ind(i),27)/x_AT(1,27);
       N_AT = x_AT(end,27)/50000;
       L_N = Expt_N_AT(i) - 1.0*Std_N_AT(i); 
       U_N = Expt_N_AT(i) + 1.0*Std_N_AT(i);
      
       if (N_Sim_AT(i)>=L_N)
           if (N_Sim_AT(i)<=U_N)
               N_sq = 0.0;
           else 
               N_sq = abs(N_Sim_AT(i)-U_N)/U_N;     
           end    
       else
           N_sq = abs(N_Sim_AT(i)-L_N)/L_N;
       end                              
       N_sq_AT = N_sq_AT + N_sq;
%---------------------------------       
       V_Sim_AT(i) = 100.0*x_AT(day_ind(i),27)/(x_AT(day_ind(i),27) + x_AT(day_ind(i),28));
       L_V = Expt_V_AT(i) - 1.5*Std_V_AT(i);
       U_V = Expt_V_AT(i) + 1.5*Std_V_AT(i);
 
      if (V_Sim_AT(i)>=L_V)
          if (V_Sim_AT(i)<=U_V)
               V_sq = 0.0;
           else 
               V_sq = abs(V_Sim_AT(i)-U_V)/U_V;     
           end    
       else
           V_sq = abs(V_Sim_AT(i)-L_V)/L_V;
       end                              
       V_sq_AT = V_sq_AT + V_sq;
   end

   Sim_Via_AT = 100.0*x_AT(:,27)./(x_AT(:,27) + x_AT(:,28));

   Via_0 = 98.4;
   Via_Diff_AT = Via_0 - Sim_Via_AT;
   R_AT = sum(Via_Diff_AT(1:end))/max(size(Via_Diff_AT));  %We redefined this in the revised manuscript (see below)
   
%==========================================================================
% (3) Calculate rSq protein value for the combination treatment protocols:
%==========================================================================
%--------------------------------------------------------------------------
% (i) Calculate the total steady-state protein levels: 
%--------------------------------------------------------------------------
   Proteins_treatment(3,:) = AT_protein_3d(1:21);
   
% The total protein array has the following order: 
% MCL-1/mcl1_0; BCL-2/bcl12_0; BIM/bim_0; BAX/bax_0; BAK/bak_0; 
% aCASP-3/casp3_0; c-Myc/cMyc_0; Chop/Chop_0.
%--------------------------------------------------------------------------
% Calculate total protein levels and rSq value for each treatment 
% condition:
%--------------------------------------------------------------------------
   Protein_ratio_rSq = 0.0;
       
   Bcl2_complex = AT_protein_3d(24);     % ABT:Bcl-2 complex level.

   %Mcl-1: 
   Total_protein(1,4) = Proteins_treatment(3,1)...
                      + Proteins_treatment(3,10)...
                      + Proteins_treatment(3,11)...
                      + Proteins_treatment(3,15)...
                      + Proteins_treatment(3,16); 
   %Bcl-2:
   Total_protein(2,4) = Proteins_treatment(3,2)...
                      + Proteins_treatment(3,19)...
                      + Proteins_treatment(3,12)...
                      + Proteins_treatment(3,17)...
                      + Proteins_treatment(3,18)... 
                      + Bcl2_complex*ABT_0/bcl2_0;     % ABT:Bcl-2 complex level.
   %Bim:
   Total_protein(3,4) = Proteins_treatment(3,3)...
                      + Proteins_treatment(3,11)*mcl1_0/bim_0...
                      + Proteins_treatment(3,12)*bcl2_0/bim_0...
                      + Proteins_treatment(3,13)*bax_0/bim_0...
                      + Proteins_treatment(3,14)*bak_0/bim_0; 
   %Bax: 
   Total_protein(4,4) = Proteins_treatment(3,4)...
                      + Proteins_treatment(3,6)...
                      + Proteins_treatment(3,13)...
                      + Proteins_treatment(3,15)*mcl1_0/bax_0...
                      + Proteins_treatment(3,17)*bcl2_0/bax_0; 
   %Bak:
   Total_protein(5,4) = Proteins_treatment(3,5)...
                      + Proteins_treatment(3,7)...
                      + Proteins_treatment(3,14)...
                      + Proteins_treatment(3,16)*mcl1_0/bak_0...
                      + Proteins_treatment(3,18)*bcl2_0/bak_0; 
   %aCasp-3:                
   Total_protein(6,4) = Proteins_treatment(3,9);     % Is actually active casp-3 level. 

   %c-Myc:
   Total_protein(7,4) = Proteins_treatment(3,20);
       
   %Chop:
   Total_protein(8,4) = Proteins_treatment(3,21);
       
       
   Protein_ratio(:,3) = Total_protein(:,4)./Total_protein(:,1);
  
% Protein data file has the following row order: MCL-1; BCL-2; BIM; BAX; 
% BAK; aCASP-3; c-Myc; Chop
   for k=1:num_proteins
       L_P = protein_expt(k,5) - protein_expt(k,6); 
       U_P = protein_expt(k,5) + protein_expt(k,6);             
       if (Protein_ratio(k,3)>=L_P)
          if (Protein_ratio(k,3)<=U_P)
              diff_sq = 0.0;
          else 
              diff_sq = abs(Protein_ratio(k,3)-U_P)/U_P;        
          end    
       else
          diff_sq = abs(Protein_ratio(k,3)-L_P)/L_P;        %Can change this to actual rsq value if you want (e.g. (Protein_sim(i,j)-L_P)^2 )
       end                              
       Protein_comparison(k,3) = Protein_ratio(k,3);                     
       Protein_ratio_rSq = Protein_ratio_rSq + diff_sq;
   end
   rSq = rSq + Protein_ratio_rSq;

%--------------------------------------------------------------------------
% Calculate total rSq value:
  protein_sq = rSq; %disp(protein_sq)
  via_sq = V_sq_0 + V_sq_A + V_sq_T + V_sq_AT; %disp(via_sq)  
  num_sq = N_sq_0 + N_sq_A + N_sq_T + N_sq_AT; %disp(num_sq)
  rSq = protein_sq + 0.2*via_sq + 0.05*num_sq;

  %disp(rSq)
  %disp(protein_sq)
proceed_3 = 1;
end

% Proceed to alternative treatment schedules:
if proceed_3 == 1
      tspan_treatment=linspace(tmin,2*tmax,2*num_points);
      sol_AT = dde15s_new(@integratingfunction_combination_ABT_Ted_no_feedback,lags,IC_protein_AT,tspan_treatment,opts_AT,DE_params,P_0_AT);

   %(1) Give ABT-199 first followed by Tedizolid/Combination treatment some 
   % hours later, then calculate the time averaged response for each
   % treatment.
    
   R_Treatment_Schedules = zeros(2,96);   %96

   Delay = zeros(1,96);
   N_live_dead_scatter = zeros(4,96);
   %Nlive max and Ndead max for ABT first, followed by Ted first 

   for i = 1:max(size(Delay))
       Delay(i) = i;       
   end
   Delay = 60.0*Delay; %Min
   Delay = Delay/t_scale; %t_scale units
   
   N_Live_plot = zeros(2,max(size(Delay)));
   NFinal = zeros(2,max(size(Delay))); 
 
   for i = 1:max(size(Delay))   
       tmin = 0.0;
       num_points = 25000;
       
       tmax_delay = Delay(i);
       tmax_total = 5*24*60/t_scale + tmax_delay;
       num_points_delay = ceil(num_points/(5*24*60/t_scale)*tmax_delay);
       num_points = num_points + num_points_delay;

       tspan_treatment=linspace(tmin,tmax_delay,num_points_delay);
 
       IC_delay = [Calculated_protein_SS,...
                   1.0, 1.0,...
                   1.0, 0.0, 0.0,...
                   50000, 813]; % Initial viability set to 98.4 = average of untreated
   
       P_0 = [mcl1_0; bcl2_0;...  
              bim_0; bax_0; bak_0; casp3_0;...
              ABT_0];   
   
       x_delay_1 = dde15s_new(@integratingfunction_ABT199_no_feedback,lags,IC_delay,tspan_treatment,opts_AT,DE_params,P_0);
       x_delay_1 = deval(x_delay_1, tspan_treatment);        
       x_delay_1 = x_delay_1.';       
 
       num_points_finish = num_points-num_points_delay; %round((tmax_total-tmax_delay)/time_step,0); 
       tspan_treatment=linspace(tmax_delay,tmax_total,num_points_finish);
       np_total = num_points_delay + num_points_finish;
       
       %To be consistent with experiments, we perform drug washout,
       %then replenish both drugs (internal concentrations should
       %remain, but external should be reset, i.e IC_finish looks like: 
       IC_finish = [x_delay_1(end,1:21),...
                    1.0, x_delay_1(end,23:24),...
                    1.0, 0.0,...
                    x_delay_1(end,25:26)]; 

       P_0 = [mcl1_0; bcl2_0;...  
              bim_0; bax_0; bak_0; casp3_0;...
              ABT_0; Ted_0];   
                 
        x_delay_2 = dde15s_new(@integratingfunction_combination_ABT_Ted_no_feedback,lags,IC_finish,tspan_treatment,opts_AT,DE_params,P_0);
        x_delay_2 = deval(x_delay_2, tspan_treatment);        
        x_delay_2 = x_delay_2.';
        
        N_live = [x_delay_1(:,25); x_delay_2(:,27)];
        N_dead = [x_delay_1(:,26); x_delay_2(:,28)];
        N_Live_plot(1,i) = N_live(end);
        N_live_dead_scatter(1,i) = max(N_live);
        N_live_dead_scatter(2,i) = max(N_dead);     
 
        NFinal(1,i) = N_live(end); 
        tspan_treatment=linspace(tmin,(tmax_total-tmax_delay),num_points_finish);

        sol_AT2 = deval(sol_AT, tspan_treatment);
        x_AT = sol_AT2.';
        N_live_0 = x_AT(:,27);
        Casp_temporal_0 = x_AT(:,9);
        Max_N_live_0_ABT = max(N_live_0);

        R_Treatment_Schedules(1,i) = (max(N_live_0)*N_live_0(end))/(max(N_live)*N_live(end));

   end

   %(2) Give Tedizolid first followed by combination treatment a number of 
   % hours later. Then calculate the time-averaged response and statistics 
   % of other quantities, such as endpoint cell numbers, for each delay to 
   % determine the optimal treatment strategy.
   max_Bak = zeros(1,96);
   max_Bax = zeros(1,96);
   max_Casp = zeros(1,96);
   for i = 1:max(size(Delay))   
       tmin = 0.0;
       num_points = 25000;
       
       tmax_delay = Delay(i);    %This is in units of t_scale now.
       tmax_total = 5*24*60/t_scale + tmax_delay;
       num_points_delay = ceil(num_points/(5*24*60/t_scale)*tmax_delay);
       num_points = num_points + num_points_delay;

       tspan_treatment=linspace(tmin,tmax_delay,num_points_delay);

       IC_delay = [Calculated_protein_SS,...
                   1.0, 1.0,...
                   1.0, 0.0,...
                   50000, 813];  %Initial viability set to average of untreated

       P_0 = [mcl1_0; bcl2_0;...  
              bim_0; bax_0; bak_0; casp3_0;...
              Ted_0];      

       x_delay_1 = dde15s_new(@integratingfunction_Tedizolid_no_feedback,lags,IC_delay,tspan_treatment,opts_AT,DE_params,P_0);
       x_delay_1 = deval(x_delay_1, tspan_treatment);              
       x_delay_1 = x_delay_1.';
        
       num_points_finish = num_points-num_points_delay;
       tspan_treatment=linspace(tmax_delay,tmax_total,num_points_finish);
       np_total = num_points_delay + num_points_finish;
                    
       %To be consistent with experiments, we perform drug washout, then 
       %replenish both drugs (i.e IC_finish looks like):
       IC_finish = [x_delay_1(end,1:21),...
                    1.0, 0.0, 0.0,...
                    1.0,...
                    x_delay_1(end,23:25)];    
      
       P_0 = [mcl1_0; bcl2_0;...  
              bim_0; bax_0; bak_0; casp3_0;...
              ABT_0; Ted_0];   

        x_delay_2 = dde15s_new(@integratingfunction_combination_ABT_Ted_no_feedback,lags,IC_finish,tspan_treatment,opts_AT,DE_params,P_0);
        x_delay_2 = deval(x_delay_2, tspan_treatment);        
        x_delay_2 = x_delay_2.';
        
        N_live = [x_delay_1(:,24); x_delay_2(:,27)];
        N_dead = [x_delay_1(:,25); x_delay_2(:,28)];
        N_Live_plot(2,i) = N_live(end);
        N_live_dead_scatter(3,i) = max(N_live);
        N_live_dead_scatter(4,i) = max(N_dead);

        Sim_Via = 100.0*N_live./(N_live + N_dead);   
        NFinal(2,i) = N_live(end); 

        tspan_treatment=linspace(tmin,(tmax_total-tmax_delay),num_points_finish);
        sol_AT2 = deval(sol_AT, tspan_treatment);
        x_AT = sol_AT2.';
        N_live_0 = x_AT(:,27);
        Casp_temporal_0 = x_AT(:,9);
        Max_N_live_0_Ted = max(N_live_0);
        
        R_Treatment_Schedules(2,i) = (max(N_live_0)*N_live_0(end))/(max(N_live)*N_live(end));

        Bak_temporal = [x_delay_1(:,7); x_delay_2(:,7)]; 
        Bax_temporal = [x_delay_1(:,6); x_delay_2(:,6)];
        Casp_temporal = [x_delay_1(:,9); x_delay_2(:,9)];
        max_Bak(i) = max(Bak_temporal/Bak_temporal(1));
        max_Bax(i) = max(Bax_temporal/Bax_temporal(1));
        max_Casp(i) = max(Casp_temporal/Casp_temporal(1));
        
        if i == 36
           Casp_temporal_36 = Casp_temporal(num_points_delay+1:np_total);
        end
   end  
   figure; 
   hold on;
   Delay = [0, Delay];
   R_AT = 1;
   R_Treatment_Schedules = [R_AT, R_Treatment_Schedules(1,:);...
                            R_AT, R_Treatment_Schedules(2,:)];                        
   plot(Delay/60*t_scale, R_Treatment_Schedules(1,:), 'ko', 'MarkerSize',10,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','k');
   plot(Delay/60*t_scale, R_Treatment_Schedules(2,:), 'bo', 'MarkerSize',10,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor','b');
   box on;
   set(gca,'LineWidth',2,'TickLength',[0.015 0.015]);
   set(gca,'XMinorTick','on','YMinorTick','on')
   set(gca,'fontsize',26);
   xlabel('Delay (hours)');
   ylabel('Response');
   legend('ABT-199 before Tedizolid', 'Tedizolid before ABT-199');
   title('Response with different treatment schedules of single doses of 5,000 nM Tedizolid and 400 nM ABT-199');


   %Plot the end-point and maximal absolute cell numbers:
   figure; 
   hold on;
   N_Live_plot = [N_AT*50000, N_Live_plot(1,:); N_AT*50000, N_Live_plot(2,:)];
   N_live_dead_scatter_plot = [Max_N_live_0_ABT, N_live_dead_scatter(1,:);...
                               Max_N_live_0_Ted, N_live_dead_scatter(3,:)];
     
   plot(Delay/60*t_scale, N_Live_plot(1,:), 'ko', 'MarkerSize',10,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','k');
   plot(Delay/60*t_scale, N_Live_plot(2,:), 'bo', 'MarkerSize',10,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor','b');
   plot(Delay/60*t_scale, N_live_dead_scatter_plot(1,:), 'ro', 'MarkerSize',10,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor','r');
   plot(Delay/60*t_scale, N_live_dead_scatter_plot(2,:), 'mo', 'MarkerSize',10,...
    'MarkerEdgeColor','m',...
    'MarkerFaceColor','m');
   box on;
   set(gca,'LineWidth',2,'TickLength',[0.015 0.015]);
   set(gca,'XMinorTick','on','YMinorTick','on')
   set(gca,'fontsize',26);
   xlabel('Delay (hours)');
   ylabel('Live cell number');
   legend('Endpoint, ABT-199 pre-treat', 'Endpoint, Tedizolid pre-treat','Peak, ABT-199 pre-treat', 'Peak, Tedizolid pre-treat');
   title('Cell numbers with different treatment schedules of single doses of 5000 nM Tedizolid and 400 nM ABT-199');

%--------------------------------------------------------------------------
%Plot active Caspase-3 levels:
   figure
   t_0 = linspace(0,120,max(size(Casp_temporal_0)));
   t_36 = linspace(0,120,max(size(Casp_temporal_36)));
   plot(t_0, Casp_temporal_0/Casp_temporal_0(1), 'k-')
   hold on
   plot(t_36, Casp_temporal_36/Casp_temporal_0(1), 'b-')  %Do not change these indices
   box on;
   set(gca,'LineWidth',2,'TickLength',[0.015 0.015]);
   set(gca,'XMinorTick','on','YMinorTick','on')
   set(gca,'fontsize',26);
   xlabel('Time (hours)');
   ylabel('Casp-3* level');
   legend('No pre-treat', 'Tedizolid pre-treat');
   title('Caspase 3 levels with different treatment schedules of 5,000 nM Tedizolid and 400 nM ABT-199 doses');
end
end