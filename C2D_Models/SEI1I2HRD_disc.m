function [xk1] = SEI1I2HRD_disc(T, xk, uk, param) 
%SEDISC discretization for the continuous SEI1I2HRD model
% * This model is a discrete Euler approximated model of the 
% * nonlinear, continuous SEI1I2HRD model
% * INPUTS
% *     T - the sampling period
% *     xk - contains states:
% *         1. S 
% *         2. E
% *         3. I1
% *         4. I2
% *         5. H
% *         6. R
% *         7. D
% *     uk - parameter coefficients as input
% *         1. u_beta
% *         2. u_gamma2
% *         3. p1
% *         4. p2
% *     param - parameter matrix, containing the following parameters:
% *         1. beta
% *         2. c
% *         3. gamma1
% *         4. gamma2
% *         5. gamma3
% *         6. gamma4
% * OUTPUTS
% *     yk -  states 
% *     Xk1 - next states, also as output, because we assume that every state
% *           may be observed
% *
    N = sum(xk);
    Sk1 = min([max([(1-T*uk(1)*param(1)*(xk(3)+xk(4)+param(2)*xk(5)))*xk(1), 0]), N]);
    Ek1 = min([max([(1-T*param(3))*xk(2) + T*uk(1)*param(1)*xk(1)*(xk(3)+xk(4)+param(2)*xk(5)), 0]), N]);
    I1k1 = min([max([(1-T*uk(2)*param(4))*xk(3) + T*param(3)*xk(2),0]), N]);
    I2k1 = min([max([T*uk(2)*param(4)*(1-uk(3))*xk(3) + (1-T*param(6))*xk(4), 0]),N]);
    Hk1 = min([max([T*uk(2)*param(4)*uk(3)*xk(3) + (1-T*param(5))*xk(5), 0]), N]);
    Rk1 = min([max([xk(6) + T*param(6)*xk(4) + T*param(5)*(1-uk(4))*xk(5), 0]), N]);
    Dk1 = min([max([xk(7) + T*param(5)*uk(4)*xk(5),0]),N]);
    
    xk1 = [Sk1, Ek1, I1k1, I2k1, Hk1, Rk1, Dk1]';
end

