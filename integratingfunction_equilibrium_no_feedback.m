%===============================================================
% This is the system of first order ODEs for the ABT-199 and   
% Tedizolid combination treatment.
%===============================================================

function dxdt = integratingfunction_equilibrium_no_feedback(t,x,DE_params_equil,TF,P)

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

%--------------------------------------------------------------------------
% Parameters
%--------------------------------------------------------------------------
params=DE_params_equil;
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
n1=params(54); n2=params(55); n3=params(56); n4=params(57); n5=params(58);
n6=params(59); n7=params(60); n8=params(61); n9=params(62);
%--------------------------------------------------------------------------
K_chop_mcl1=params(63); n11=params(64);
%--------------------------------------------------------------------------
n13=params(65); K_myc_bim_2=params(66);
n14=params(67); K_myc_bcl2_2=params(68);
%--------------------------------------------------------------------------
myc_tilde = TF(1);
chop_tilde = TF(2);

mcl1_0=P(1); bcl2_0=P(2); bim_0=P(3); bax_0=P(4); bak_0=P(5); casp3_0=P(6);

%--------------------------------------------------------------------------
% Equations defining output:  
%--------------------------------------------------------------------------
dxdt = zeros(max(size(x)),1);

%(1) Mcl-1:
dxdt(1) = k_mcl1_tilde*myc_tilde^n1/(K_myc_mcl1^n1+myc_tilde^n1)/(1.0+(chop_tilde/K_chop_mcl1)^n11)...
    + kb_mcl1_tilde... 
    - dp_mcl1_tilde*x(1)...
    - kf_mcl1_bim_tilde*bim_0*x(1)*x(3)...
    + kr_mcl1_bim_tilde*x(11)... 
    - (kf_mcl1_bax_tilde*bax_0*x(1)*x(6) + kf_mcl1_bak_tilde*bak_0*x(1)*x(7))...
    + kr_mcl1_bax_tilde*x(15) + kr_mcl1_bak_tilde*x(16)...
    - kcc_tilde*casp3_0*x(1)*x(9);
    
%(2) Bcl-2:
dxdt(2) = k_bcl2_tilde*myc_tilde^n14/(K_myc_bcl2_2^n14+myc_tilde^n14)*(1+myc_tilde^n2/K_myc_bcl2^n2)^(-1)*(1+chop_tilde^n3/K_chop_bcl2^n3)^(-1)...
    + kb_bcl2_tilde...
    - kf_bcl2_bim_tilde*bim_0*x(2)*x(3) + kr_bcl2_bim_tilde*x(12)...
    - (kf_bcl2_bax_tilde*bax_0*x(2)*x(6) + kf_bcl2_bak_tilde*bak_0*x(2)*x(7))...
    + kr_bcl2_bax_tilde*x(17) + kr_bcl2_bak_tilde*x(18)...
    - kcc_tilde*casp3_0*x(2)*x(9)...
    - dp_bcl2_tilde*x(2);

%(3) Bim:
dxdt(3) = k_bim_tilde*( (1.0 - myc_tilde^n5/(K_myc_bim^n5 + myc_tilde^n5))*chop_tilde^n4/(K_chop_bim^n4 + chop_tilde^n4)/(1.0 + (myc_tilde/K_myc_bim_2)^n13)...
    + myc_tilde^n5/(K_myc_bim^n5 + myc_tilde^n5) )...
    + kb_bim_tilde...
    - dp_bim_tilde*x(3)... 
    - kf_mcl1_bim_tilde*mcl1_0*x(1)*x(3) + kr_mcl1_bim_tilde*x(11)*mcl1_0/bim_0...
    - kf_bcl2_bim_tilde*bcl2_0*x(2)*x(3) + kr_bcl2_bim_tilde*x(12)*bcl2_0/bim_0...
    - kf_bax_bim_tilde*bax_0*x(3)*x(4) + kr_bax_bim_tilde*x(13)*bax_0/bim_0...
    + kappa_bax_bim_tilde*x(13)*bax_0/bim_0...
    - kf_bak_bim_tilde*bak_0*x(3)*x(5) + kr_bak_bim_tilde*x(14)*bak_0/bim_0...
    + kappa_bak_bim_tilde*x(14)*bak_0/bim_0;

%(4) Bax:
dxdt(4) = k_bax_tilde*(chop_tilde^n6/K_chop_bax^n6 + myc_tilde^n7/K_myc_bax^n7 + chop_tilde^n6*myc_tilde^n7/(K_chop_bax^n6*K_myc_bax^n7))...
                     /(1.0+chop_tilde^n6/K_chop_bax^n6 + myc_tilde^n7/K_myc_bax^n7 + chop_tilde^n6*myc_tilde^n7/(K_chop_bax^n6*K_myc_bax^n7))...
    + kb_bax_tilde...
    - dp_bax_tilde*x(4)... 
    - kf_bax_bim_tilde*bim_0*x(3)*x(4) + kr_bax_bim_tilde*x(13);

%(5) Bak:
dxdt(5) = k_bak_tilde*(chop_tilde^n8/K_chop_bak^n8 + myc_tilde^n9/K_myc_bak^n9 + chop_tilde^n8*myc_tilde^n9/(K_chop_bak^n8*K_myc_bak^n9))...
                     /(1.0+chop_tilde^n8/K_chop_bak^n8 + myc_tilde^n9/K_myc_bak^n9 + chop_tilde^n8*myc_tilde^n9/(K_chop_bak^n8*K_myc_bak^n9))...
    + kb_bak_tilde...
    - dp_bak_tilde*x(5)... 
    - kf_bak_bim_tilde*bim_0*x(3)*x(5) + kr_bak_bim_tilde*x(14);

%(6) Bax*:
dxdt(6) =  -dp_bax_tilde*x(6)... 
    + kappa_bax_bim_tilde*x(13)...
    - kf_bcl2_bax_tilde*bcl2_0*x(2)*x(6) + kr_bcl2_bax_tilde*x(17)*bcl2_0/bax_0...
    - kf_mcl1_bax_tilde*mcl1_0*x(1)*x(6) + kr_mcl1_bax_tilde*x(15)*mcl1_0/bax_0;

%(7) Bak*:
dxdt(7) =  -dp_bak_tilde*x(7)... 
    + kappa_bak_bim_tilde*x(14)...
    - kf_bcl2_bak_tilde*bcl2_0*x(2)*x(7) + kr_bcl2_bak_tilde*x(18)*bcl2_0/bak_0...
    - kf_mcl1_bak_tilde*mcl1_0*x(1)*x(7) + kr_mcl1_bak_tilde*x(16)*mcl1_0/bak_0;

%(8) Casp-3:
dxdt(8) =  kb_casp3_tilde...
    - dp_casp3_tilde*x(8)... 
    - kc_tilde*bax_0*x(6)*x(8)...
    - kc_tilde*bak_0*x(7)*x(8);

%(9) Casp-3*:
dxdt(9) = -dp_casp3_tilde*x(9)... 
    + kc_tilde*bax_0*x(6)*x(8)...
    + kc_tilde*bak_0*x(7)*x(8);

%(10) Mcl-1i:
dxdt(10) = -dp_mcl1_tilde*x(10) + kcc_tilde*casp3_0*x(1)*x(9);

%(11) Mcl-1:Bim normalized to mcl1_0
dxdt(11) = kf_mcl1_bim_tilde*bim_0*x(1)*x(3) - kr_mcl1_bim_tilde*x(11)...
         - dp_mcl1_bim_tilde*x(11);

%(12) Bcl-2:Bim normalized to bcl2_0
dxdt(12) = kf_bcl2_bim_tilde*bim_0*x(2)*x(3) - kr_bcl2_bim_tilde*x(12)...
         - dp_bcl2_bim_tilde*x(12);

%(13) Bim:Bax normalized to bax_0
dxdt(13) = kf_bax_bim_tilde*bim_0*x(3)*x(4) - kr_bax_bim_tilde*x(13)...
         - kappa_bax_bim_tilde*x(13) - dp_bax_bim_tilde*x(13);

%(14) Bim:Bak normalized to bak_0
dxdt(14) = kf_bak_bim_tilde*bim_0*x(3)*x(5) - kr_bak_bim_tilde*x(14)...
         - kappa_bak_bim_tilde*x(14) - dp_bak_bim_tilde*x(14);

%(15) Mcl-1:aBax normalized to mcl1_0
dxdt(15) = kf_mcl1_bax_tilde*bax_0*x(1)*x(6) - kr_mcl1_bax_tilde*x(15)...
         - dp_mcl1_bax_tilde*x(15);

%(16) Mcl-1:aBak normalized to mcl1_0
dxdt(16) = kf_mcl1_bak_tilde*bak_0*x(1)*x(7) - kr_mcl1_bak_tilde*x(16)...
         - dp_mcl1_bak_tilde*x(16);

%(17) Bcl-2:aBax normalized to bcl2_0
dxdt(17) = kf_bcl2_bax_tilde*bax_0*x(2)*x(6) - kr_bcl2_bax_tilde*x(17)...
         - dp_bcl2_bax_tilde*x(17);

%(18) Bcl-2:aBak normazlized to bcl2_0
dxdt(18) = kf_bcl2_bak_tilde*bak_0*x(2)*x(7) - kr_bcl2_bak_tilde*x(18)...
         - dp_bcl2_bak_tilde*x(18);

%(19) Bcl-2i normalized to bcl2_0
dxdt(19) = kcc_tilde*casp3_0*x(2)*x(9) - dp_bcl2_tilde*x(19);

end
