%===============================================================
% This is the system of first order ODEs for the ABT-199 and   
% Tedizolid combination treatment.
%===============================================================

function dxdt = integratingfunction_untreated_N(t,x,DE_params_N)

%--------------------------------------------------------------------------
% Variable definitions:
%--------------------------------------------------------------------------
%N=x(1); number of live cells 
%D=x(2); number of dead cells
%--------------------------------------------------------------------------
% Parameters
%--------------------------------------------------------------------------
params=DE_params_N;

Bcl2_0=params(1);
Myc_0=params(2);
aCasp3_frac_0=params(3);
lambda_bcl2=params(4);
lambda_myc=params(5);
lambda_0=params(6);
Bcl2_thresh=params(7);
Myc_thresh=params(8);
delta_casp3=params(9);
Casp3_thresh=params(10);
n12=params(11);
nv_1=params(12);
nv_2=params(13);

%--------------------------------------------------------------------------
% Equations defining output:  
%--------------------------------------------------------------------------
dxdt = zeros(max(size(x)),1);

%The terms represent cellular respiration + glycolysis - apoptosis,
%where 
%cellular respiration depends on Bcl-2 concentration (with threshold), 
%glycolysis depends on c-Myc levels (with threshold),
%apoptosis depends on active caspase 3 (with threshold).
%we also include a constant background production term 

%(1) Live cell count:
dxdt(1) = (lambda_0 + lambda_bcl2*Bcl2_0^nv_1/(Bcl2_thresh^nv_1+Bcl2_0^nv_1)...
        + lambda_myc*Myc_0^nv_2/(Myc_thresh^nv_2+Myc_0^nv_2)...
        - delta_casp3*aCasp3_frac_0^n12/(Casp3_thresh^n12+aCasp3_frac_0^n12))*x(1);

%(2) Dead cell count:
dxdt(2) = (delta_casp3*aCasp3_frac_0^n12/(Casp3_thresh^n12+aCasp3_frac_0^n12))*x(1);

end