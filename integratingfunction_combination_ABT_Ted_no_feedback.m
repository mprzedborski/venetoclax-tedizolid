function dxdt = integratingfunction_combination_ABT_Ted_no_feedback(t,x,Z,DE_params,P)
%--------------------------------------------------------------------------
% Variable definitions:
%--------------------------------------------------------------------------
%mcl/mcl_0=x(1); bcl2/bcl2_0=x(2); bim/bim_0=x(3); bax/bax_0=x(4);
%bak/bak_0=x(5); abax/bax_0=x(6); abak/bak_0=x(7); casp3/casp3_0=x(8); 
%acasp3/casp3_0=x(9); 
%mcli/mcl1_0=x(10); mclbim/mcl1_0=x(11); bcl2bim/bcl2_0=x(12); 
%bimbax/bax_0=x(13); bimbak/bak_0=x(14);        
%mclabax/mcl1_0=x(15); mclabak/mcl1_0=x(16); 
%bcl2abax/bcl2_0=x(17); bcl2abak/bcl2_0=x(18); 
%bcl2i/bcl2_0=x(19); 
%cmyc/cMyc_0=x(20); chop/Chop_0=x(21)
%--------------------------------------------------------------------------
%external ABT concentration/ABT_0=x(22); free internal ABT concentration/ABT_0=x(23); 
%ABT-Bcl-2 complex/ABT_0=x(24);
%--------------------------------------------------------------------------
%external Tedizolid concentration/Ted_0=x(25); 
%intracellular Tedizolid concentration/Ted_0=x(26); 
%--------------------------------------------------------------------------
%live cell number=x(27); dead cell number=x(28);
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Parameters
%--------------------------------------------------------------------------
% Read in both drug parameters for both monotherapies.
%--------------------------------------------------------------------------
params=DE_params;
%--------------------------------------------------------------------------
% Native protein parameters:
%--------------------------------------------------------------------------
kc_tilde=params(1); kcc_tilde=params(2); 
%--------------------------------------------------------------------------
k_mcl1_tilde=params(3); k_bcl2_tilde=params(4); 
k_bim_tilde=params(5); k_bax_tilde=params(6); k_bak_tilde=params(7);
%--------------------------------------------------------------------------
kb_mcl1_tilde=params(8); kb_bcl2_tilde=params(9); 
kb_bim_tilde=params(10); kb_bax_tilde=params(11); kb_bak_tilde=params(12);
%--------------------------------------------------------------------------
kb_casp3_tilde=params(13); dp_casp3_tilde=params(13);
%--------------------------------------------------------------------------
K_myc_mcl1=params(14); K_myc_bcl2=params(15); 
K_myc_bim=params(16); K_myc_bax=params(17); K_myc_bak=params(18);
K_chop_bcl2=params(19);
K_chop_bim=params(20); K_chop_bax=params(21); K_chop_bak=params(22);
%--------------------------------------------------------------------------
dp_mcl1_tilde=params(23); dp_bcl2_tilde=params(24); 
dp_bim_tilde=params(25); dp_bax_tilde=params(26); dp_bak_tilde=params(27);
%--------------------------------------------------------------------------
dp_bcl2_bim_tilde=params(28); dp_mcl1_bim_tilde=params(29);
dp_bax_bim_tilde=params(30); dp_bak_bim_tilde=params(31);
dp_bcl2_bax_tilde=params(32); dp_bcl2_bak_tilde=params(33);
dp_mcl1_bax_tilde=params(34); dp_mcl1_bak_tilde=params(35);

kf_mcl1_bim_tilde=params(36); kr_mcl1_bim_tilde=params(37);
kf_mcl1_bax_tilde=params(38); kr_mcl1_bax_tilde=params(39);
kf_mcl1_bak_tilde=params(40); kr_mcl1_bak_tilde=params(41);

kf_bcl2_bim_tilde=params(42); kr_bcl2_bim_tilde=params(43);
kf_bcl2_bax_tilde=params(44); kr_bcl2_bax_tilde=params(45);
kf_bcl2_bak_tilde=params(46); kr_bcl2_bak_tilde=params(47);

kf_bax_bim_tilde=params(48); kr_bax_bim_tilde=params(49); 
kappa_bax_bim_tilde=params(50);
kf_bak_bim_tilde=params(51); kr_bak_bim_tilde=params(52); 
kappa_bak_bim_tilde=params(53);

%--------------------------------------------------------------------------
% Transcription factor parameters:
%--------------------------------------------------------------------------
kA_uptake=params(54); dA_tilde=params(55); dE_ABT=params(56);
kT_uptake=params(57); dT_tilde=params(58); dE_Ted=params(59);
kb_myc_tilde=params(60); kmax_myc_tilde=params(61);
Kf1_myc_A_tilde=params(62); nf1_myc_A=params(63);
Kf1_myc_T_tilde=params(64); nf1_myc_T=params(65);
Kg1_myc_A_tilde=params(66); ng1_myc_A=params(67);
Kg1_myc_T_tilde=params(68); ng1_myc_T=params(69);
dp_myc_max_tilde=params(70); dp_myc_min_tilde=params(71); dp_myc_0_tilde=params(72);
kb_chop_tilde=params(73); kmax_chop_tilde=params(74); 
Kf1_chop_A_tilde=params(75); nf1_chop_A=params(76); 
Kf1_chop_T_tilde=params(77); nf1_chop_T=params(78); 
Kg1_chop_A_tilde=params(79); ng1_chop_A=params(80);
Kg1_chop_T_tilde=params(81); ng1_chop_T=params(82);
dp_chop_max_tilde=params(83); dp_chop_min_tilde=params(84); dp_chop_0_tilde=params(85);
w1_myc=params(86); w2_myc=params(87); w3_myc=params(88); w4_myc=params(89);
w5_myc=params(90); w6_myc=params(91); w7_myc=params(92); w8_myc=params(93);
w1_chop=params(94); w2_chop=params(95); w3_chop=params(96); w4_chop=params(97);
w5_chop=params(98); w6_chop=params(99); w7_chop=params(100); w8_chop=params(101);
q_myc=params(102); q_chop=params(103);

%--------------------------------------------------------------------------
% Other drug-related parameters:
%--------------------------------------------------------------------------
kf_abt_bcl2_tilde=params(104); kr_abt_bcl2_tilde=params(105); 
kted_tilde=params(106); kted_2_tilde=params(107);
K_ABT_bax_tilde=params(108);
%--------------------------------------------------------------------------
% Cell viability parameters
%--------------------------------------------------------------------------
lambda_bcl2_tilde=params(109); lambda_myc_tilde=params(110); 
lambda_0_tilde=params(111); delta_casp3_tilde=params(112); 
Bcl2_thresh=params(113); q_Ted_tilde=params(114);
%--------------------------------------------------------------------------
n1=params(115); n2=params(116); n3=params(117); n4=params(118); n5=params(119);
n6=params(120); n7=params(121); n8=params(122); n9=params(123); n10=params(124);
%--------------------------------------------------------------------------
n_q_Ted=params(125);
Myc_thresh=params(126); Casp3_thresh=params(127);
K_chop_mcl1=params(128); n11=params(129); n12=params(130);
nv_1=params(131); nv_2=params(132);
%--------------------------------------------------------------------------
n13=params(133); K_myc_bim_2=params(134);
n14=params(135); K_myc_bcl2_2=params(136);

%--------------------------------------------------------------------------
mcl1_0=P(1); bcl2_0=P(2); bim_0=P(3); bax_0=P(4); bak_0=P(5); casp3_0=P(6);
ABT_0=P(7); Ted_0=P(8);

%--------------------------------------------------------------------------
xlag1 = Z(:,1); xlag2 = Z(:,2); %xlag5 = Z(:,5);
%Chop:
xlag3 = Z(:,3); xlag4 = Z(:,4);

%--------------------------------------------------------------------------
dxdt = zeros(max(size(x)),1);
%--------------------------------------------------------------------------
% Equations defining output:  
%--------------------------------------------------------------------------
% We assume that production of c-Myc and Chop depend on total internal
% ABT-199 concentration (not just free ABT)
% Total internal ABT concentration:
ABT = x(23) + x(24);

%(1) Mcl-1:
dxdt(1) = k_mcl1_tilde*x(20)^n1/(K_myc_mcl1^n1+x(20)^n1)/(1.0+(x(21)/K_chop_mcl1)^n11)...
    + kb_mcl1_tilde... 
    - kf_mcl1_bim_tilde*bim_0*x(1)*x(3) + kr_mcl1_bim_tilde*x(11)... 
    - (kf_mcl1_bax_tilde*bax_0*x(1)*x(6) + kf_mcl1_bak_tilde*bak_0*x(1)*x(7))...
    + kr_mcl1_bax_tilde*x(15) + kr_mcl1_bak_tilde*x(16)...
    - kcc_tilde*casp3_0*x(1)*x(9)...
    - dp_mcl1_tilde*x(1)...
    - kted_2_tilde*Ted_0*x(26)*x(1);  %drug interaction

%(2) Bcl-2:
dxdt(2) = k_bcl2_tilde*x(20)^n14/(K_myc_bcl2_2^n14+x(20)^n14)*(1+x(20)^n2/K_myc_bcl2^n2)^(-1)*(1+x(21)^n3/K_chop_bcl2^n3)^(-1)...
    + kb_bcl2_tilde...
    - kf_bcl2_bim_tilde*bim_0*x(2)*x(3) + kr_bcl2_bim_tilde*x(12)...
    - (kf_bcl2_bax_tilde*bax_0*x(2)*x(6) + kf_bcl2_bak_tilde*bak_0*x(2)*x(7))...
    + kr_bcl2_bax_tilde*x(17) + kr_bcl2_bak_tilde*x(18)...
    - kcc_tilde*casp3_0*x(2)*x(9)...
    - dp_bcl2_tilde*x(2)...
    - kf_abt_bcl2_tilde*x(2)*ABT_0*x(23) + kr_abt_bcl2_tilde*ABT_0*x(24)/bcl2_0;   % Direct effects of drugs.   

%(3) Bim:
dxdt(3) = k_bim_tilde*( (1.0 - x(20)^n5/(K_myc_bim^n5 + x(20)^n5))*x(21)^n4/(K_chop_bim^n4 + x(21)^n4)/(1.0 + (x(20)/K_myc_bim_2)^n13)...
    + x(20)^n5/(K_myc_bim^n5 + x(20)^n5) )...
    + kb_bim_tilde...
    - dp_bim_tilde*x(3)... 
    - kf_mcl1_bim_tilde*mcl1_0*x(1)*x(3) + kr_mcl1_bim_tilde*x(11)*mcl1_0/bim_0...
    - kf_bcl2_bim_tilde*bcl2_0*x(2)*x(3) + kr_bcl2_bim_tilde*x(12)*bcl2_0/bim_0...
    - kf_bax_bim_tilde*bax_0*x(3)*x(4) + kr_bax_bim_tilde*x(13)*bax_0/bim_0...
    + kappa_bax_bim_tilde*x(13)*bax_0/bim_0...
    - kf_bak_bim_tilde*bak_0*x(3)*x(5) + kr_bak_bim_tilde*x(14)*bak_0/bim_0...
    + kappa_bak_bim_tilde*x(14)*bak_0/bim_0...
    + kted_tilde*Ted_0*x(26)*x(11)*mcl1_0/bim_0;       % Tedizolid inactivates Mcl-1 which causes unbinding of Mcl-1:Bim complex.

%(4) Bax:
dxdt(4) = k_bax_tilde*(x(21)^n6/K_chop_bax^n6 + x(20)^n7/K_myc_bax^n7 + x(21)^n6*x(20)^n7/(K_chop_bax^n6*K_myc_bax^n7))...
                     /(1.0+x(21)^n6/K_chop_bax^n6 + x(20)^n7/K_myc_bax^n7 + x(21)^n6*x(20)^n7/(K_chop_bax^n6*K_myc_bax^n7))...
                     /(1.0+ABT^n10/K_ABT_bax_tilde^n10)...
    + kb_bax_tilde...
    - dp_bax_tilde*x(4)... 
    - kf_bax_bim_tilde*bim_0*x(3)*x(4) + kr_bax_bim_tilde*x(13);

%(5) Bak:
dxdt(5) = k_bak_tilde*(x(21)^n8/K_chop_bak^n8 + x(20)^n9/K_myc_bak^n9 + x(21)^n8*x(20)^n9/(K_chop_bak^n8*K_myc_bak^n9))...
                     /(1.0+x(21)^n8/K_chop_bak^n8 + x(20)^n9/K_myc_bak^n9 + x(21)^n8*x(20)^n9/(K_chop_bak^n8*K_myc_bak^n9))...
    + kb_bak_tilde...
    - dp_bak_tilde*x(5)... 
    - kf_bak_bim_tilde*bim_0*x(3)*x(5) + kr_bak_bim_tilde*x(14);

%(6) Bax*:
dxdt(6) =  -dp_bax_tilde*x(6)... 
    + kappa_bax_bim_tilde*x(13)...
    - kf_bcl2_bax_tilde*bcl2_0*x(2)*x(6) + kr_bcl2_bax_tilde*x(17)*bcl2_0/bax_0...
    - kf_mcl1_bax_tilde*mcl1_0*x(1)*x(6) + kr_mcl1_bax_tilde*x(15)*mcl1_0/bax_0...
    + kted_tilde*Ted_0*x(26)*x(15)*mcl1_0/bax_0;        % Tedizolid inactivates Mcl-1 which causes unbinding of Mcl-1:aBax complex.

%(7) Bak*:
dxdt(7) =  -dp_bak_tilde*x(7)... 
    + kappa_bak_bim_tilde*x(14)...
    - kf_bcl2_bak_tilde*bcl2_0*x(2)*x(7) + kr_bcl2_bak_tilde*x(18)*bcl2_0/bak_0...
    - kf_mcl1_bak_tilde*mcl1_0*x(1)*x(7) + kr_mcl1_bak_tilde*x(16)*mcl1_0/bak_0...
    + kted_tilde*Ted_0*x(26)*x(16)*mcl1_0/bak_0;        % Tedizolid inactivates Mcl-1 which causes unbinding of Mcl-1:aBak complex.

%(8) Casp-3:
dxdt(8) =  kb_casp3_tilde...
    - dp_casp3_tilde*x(8)... 
    - kc_tilde*(bax_0*x(6)*x(8) + bak_0*x(7)*x(8));

%(9) Casp-3*:
dxdt(9) = -dp_casp3_tilde*x(9)... 
    + kc_tilde*(bax_0*x(6)*x(8) + bak_0*x(7)*x(8));

%(10) Mcl-1-:
dxdt(10) = -dp_mcl1_tilde*x(10) + kcc_tilde*casp3_0*x(1)*x(9)...
         + kted_2_tilde*Ted_0*x(26)*x(1)...
         + kted_tilde*Ted_0*x(26)*x(11)...
         + kted_tilde*Ted_0*x(26)*x(15)...
         + kted_tilde*Ted_0*x(26)*x(16);    

%(11) Mcl-1:Bim
dxdt(11) = kf_mcl1_bim_tilde*bim_0*x(1)*x(3) - kr_mcl1_bim_tilde*x(11)...
         - dp_mcl1_bim_tilde*x(11)...
         - kted_tilde*Ted_0*x(26)*x(11);
     
%(12) Bcl-2:Bim
dxdt(12) = kf_bcl2_bim_tilde*bim_0*x(2)*x(3) - kf_bcl2_bim_tilde*x(12)...
         - dp_bcl2_bim_tilde*x(12);

%(13) Bim:Bax
dxdt(13) = kf_bax_bim_tilde*bim_0*x(3)*x(4) - kr_bax_bim_tilde*x(13)...
         - kappa_bax_bim_tilde*x(13) - dp_bax_bim_tilde*x(13);

%(14) Bim:Bak
dxdt(14) = kf_bak_bim_tilde*bim_0*x(3)*x(5) - kf_bak_bim_tilde*x(14)...
         - kappa_bak_bim_tilde*x(14) - dp_bak_bim_tilde*x(14);   

%(15) Mcl-1:aBax
dxdt(15) = kf_mcl1_bax_tilde*bax_0*x(1)*x(6) - kr_mcl1_bax_tilde*x(15)...
         - dp_mcl1_bax_tilde*x(15)...
         - kted_tilde*Ted_0*x(26)*x(15);
         
%(16) Mcl-1:aBak
dxdt(16) = kf_mcl1_bak_tilde*bak_0*x(1)*x(7) - kr_mcl1_bak_tilde*x(16)...
         - dp_mcl1_bak_tilde*x(16)...
         - kted_tilde*Ted_0*x(26)*x(16);
          
%(17) Bcl-2:aBax
dxdt(17) = kf_bcl2_bax_tilde*bax_0*x(2)*x(6) - kr_bcl2_bax_tilde*x(17)...
         - dp_bcl2_bax_tilde*x(17);
     
%(18) Bcl-2:aBak
dxdt(18) = kf_bcl2_bak_tilde*bak_0*x(2)*x(7) - kr_bcl2_bak_tilde*x(18)...
         - dp_bcl2_bak_tilde*x(18);

%(19) Bcl-2i
dxdt(19) = kcc_tilde*casp3_0*x(2)*x(9) - dp_bcl2_tilde*x(19);

%--------------------------------------------------------------------------
% We assume that production of c-Myc and Chop depend on total internal
% ABT-199 concentration (not just free ABT)
% Total internal ABT concentration:
ABT = x(23) + x(24);
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% c-Myc Functions
%--------------------------------------------------------------------------
% f functions are for transcription (1- increased transcription; 2- 
% decreased transcription). 
% g functions are for stability effects (1- increased stability; 
% 2- decreased stability). 

% Combination therapy:
rf1_myc_A = ((xlag1(23)+xlag1(24))/Kf1_myc_A_tilde)^nf1_myc_A;
rf1_myc_T = (xlag2(26)/Kf1_myc_T_tilde)^nf1_myc_T;
c_myc_f1 = w1_myc*rf1_myc_A + w2_myc*rf1_myc_T + w1_myc*rf1_myc_A*w2_myc*rf1_myc_T;

rg1_myc_A = (ABT/Kg1_myc_A_tilde)^ng1_myc_A;
rg1_myc_T = (x(26)/Kg1_myc_T_tilde)^ng1_myc_T;
c_myc_g1 = w5_myc*rg1_myc_A + w6_myc*rg1_myc_T + q_myc*w5_myc*rg1_myc_A*w6_myc*rg1_myc_T;


%--------------------------------------------------------------------------
% Chop Functions
%--------------------------------------------------------------------------
% f functions are for stability (1- increased stability; 2- decreased 
% stability). 

% g functions are for stability effects (1- increased stability; 
% 2- decreased stability). 

% Combination therapy:
rf1_chop_A = (ABT/Kf1_chop_A_tilde)^nf1_chop_A;
rf1_chop_T = (xlag4(26)/Kf1_chop_T_tilde)^nf1_chop_T;
c_chop_f1 = w1_chop*rf1_chop_A + w2_chop*rf1_chop_T...
          + q_chop*w1_chop*rf1_chop_A*w2_chop*rf1_chop_T;
rg1_chop_A = (ABT/Kg1_chop_A_tilde)^ng1_chop_A;
rg1_chop_T = (xlag3(26)/Kg1_chop_T_tilde)^ng1_chop_T;
c_chop_g1 = w5_chop*rg1_chop_A + w6_chop*rg1_chop_T + w5_chop*rg1_chop_A*w6_chop*rg1_chop_T;  


%(20) cMyc levels: (Add weights to the combination treatment)
dxdt(20) = kb_myc_tilde*(1.0 + kmax_myc_tilde*c_myc_f1/(1.0 + c_myc_f1))... 
        - (dp_myc_0_tilde...
        + (dp_myc_max_tilde - dp_myc_0_tilde)*c_myc_g1/(1.0 + c_myc_g1))*x(20);
           
%(21) Chop levels:
dxdt(21) = kb_chop_tilde*(1.0...
         + kmax_chop_tilde*c_chop_f1/(1.0 + c_chop_f1))... 
         - (dp_chop_0_tilde...
         + (dp_chop_max_tilde - dp_chop_0_tilde)*c_chop_g1/(1.0 + c_chop_g1))*x(21);
 
%--------------------------------------------------------------------------
%(22) External ABT-199 concentration:
dxdt(22) = -dE_ABT*x(22);

%(23) Free internal ABT-199 concentration:
dxdt(23) = kA_uptake*x(22)...
         - kf_abt_bcl2_tilde*bcl2_0*x(2)*x(23) + kr_abt_bcl2_tilde*x(24)...
         - dA_tilde*x(23);   % Direct effects of drugs.   

%(24) ABT-199:bcl-2 complex concentration:
dxdt(24) = kf_abt_bcl2_tilde*bcl2_0*x(2)*x(23) - kr_abt_bcl2_tilde*x(24)...
         - dA_tilde*x(24);  
%--------------------------------------------------------------------------
%(25) External Tedizolid concentration:
dxdt(25) = -dE_Ted*x(25);

%(26) Intracellular Tedizolid concentration:
dxdt(26) = kT_uptake*x(25) - dT_tilde*x(26);
%--------------------------------------------------------------------------
% Active Caspase-3 fraction:
aCasp3_frac = x(9)/(x(8)+x(9));

%(27) Number of live cells:
dxdt(27) = ((lambda_0_tilde + lambda_bcl2_tilde*x(2)^nv_1/(Bcl2_thresh^nv_1+x(2)^nv_1)...
         + lambda_myc_tilde*x(20)^nv_2/(Myc_thresh^nv_2+x(20)^nv_2))/(1.0+(x(26)/q_Ted_tilde)^n_q_Ted)...
         - delta_casp3_tilde*aCasp3_frac^n12/(Casp3_thresh^n12+aCasp3_frac^n12))*x(27);      

%(28) Dead cell count:
dxdt(28) = (delta_casp3_tilde*aCasp3_frac^n12/(Casp3_thresh^n12+aCasp3_frac^n12))*x(27);

end
