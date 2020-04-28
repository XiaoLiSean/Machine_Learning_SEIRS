function [x_dot]=SEI1I2HRD_cont(t,x,beta,c,gamma1,gamma2,gamma3,gamma4,p1,p2)
% SEI1I2HRD_cont: nonlinear, complete, continuous model for the SEI1I2HRD model         
% * INPUTS: 
% *     t - time, just for ode solvers
% *     x - state vector. x(1) to x(7) are S, E, I1, I2, H, R and D,
% *         respectively.
% *     beta, c, ... p2 - parameters in the SEI1I2HRD model
% * OUTPUTS:
% *     xdot - deriative of the state vector
    [S,E,I1,I2,H,R,D]=deal(x(1),x(2),x(3),x(4),x(5),x(6),x(7));
    Sdot=-beta*S*(I1+I2+c*H);
    Edot=beta*S*(I1+I2+c*H)-gamma1*E;
    I1dot=gamma1*E-gamma2*I1;
    I2dot=gamma2*(1-p1)*I1-gamma4*I2;
    Hdot=gamma2*p1*I1-gamma3*H;
    Rdot=gamma4*I2+gamma3*(1-p2)*H;
    Ddot=gamma3*p2*H;
    x_dot=[Sdot,Edot,I1dot,I2dot,Hdot,Rdot,Ddot]';
end