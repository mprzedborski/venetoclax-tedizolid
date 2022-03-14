function [position,isterminal,direction] = EquilibriumEventFcn_no_feedback(t,x,DE_params_equil,TF,P,t_scale,XX)
% See https://www.mathworks.com/help/matlab/math/ode-event-location.html
% for more details.
% This function is used to find an equilibrium steady state for the
% untreated apoptosis-ISR pathway for Molm-13 leukemia cells.
isterminal = 1;
direction = [];
threshold = 0.005;

%--------------------------------------------------------------------------
dxdt = integratingfunction_equilibrium_no_feedback(t,x,DE_params_equil,TF,P);
%--------------------------------------------------------------------------
derivatives_vec=XX.*dxdt/t_scale;
norm_vec=norm(derivatives_vec);

% Logical(cond) evaluates to 1 if cond is true and to 0 if cond is false.
% So we want to make the condition the opposite of what we actually want,
% e.g. check if norm_vec is above the threshold. Then when it is below,
% position will drop to zero and the integration will stop.
if norm_vec > threshold
   position = 1; 
else
   position = 0;
end

end

