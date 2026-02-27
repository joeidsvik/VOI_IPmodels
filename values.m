function [v_cont, v_stop] = values(T_m, C_inj, C_leak, tax, r, t, T)
%This function computes the prospect values for continuing and stopping 
%injection at a particular time t
%   T_m = times of leakage from the topmost layer for 
%         different samples (vector)
%   C_inj = cost of injection per volume
%   C_leak = cost of leakage per volume
%   tax = tax per volume not injected
%   r: volume rate of CO2 injection
%   t = time at which the prospect values are to be computed
%   T = total planned time period of injection
%   v_cont = prospect values at time t for continuing injection for
%            different samples (vector)
%   v_stop = prospect values at time t for stopping injection for
%            different samples (vector)
    v_cont = - C_inj*(T-t)*r - C_leak*max(T-T_m,0)*r;
    v_stop = - tax*(T-t)*r;
end

