global N r beta c;
N       = 1.1e6; % Total population in our imaginary city
r       = 0.5; % Number of contacts per person per day
beta    = 0.5; % Infection rate if contact patients
c       = 0.1; % Proportion of severe patients who could not be hospitalized

X0      = [N; 4; 1e4; 0; 0; 0; 0]; % Initial State
tSpan   = 0:0.1:200; % Simulation time span
[T,X]   = ode45(@infectious_ODE, tSpan, X0);
save('Infectious_data', 'T', 'X');