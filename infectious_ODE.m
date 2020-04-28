function X_dot = infectious_ODE(T, X)
    % ---------------------------------------------------------------------
    % Infectious ODE model: Credit to "Epidemiological parameter review 
    % and comparative dynamics of influenza, respiratory syncytial virus,
    % rhinovirus, human coronavirus, and adenovirus"
    % ---------------------------------------------------------------------
    % Model States:
    % ---------------------------------------------------------------------
    % S     Number of susceptible individuals
    % E     Number of exposed (not infectious) individuals
    % I1    Number of initially infectious individuals
    % I2    Number of infected, non-hospitalized individuals
    % H     Number of hospitalized individuals
    % R     Number of recovered individuals
    % D     Number of dead individuals
    % ---------------------------------------------------------------------
    
    global N r beta c;
    gamma1  = 0.1961; % per capita rate of progress from exposed to infectious state
    gamma2  = 0.1176; % per capita rate of progress through initial infectious state
    gamma3  = 0.0286; % per capita rate of progress through hospitalized state
    gamma4  = 0.1818; % per capita rate of progress through non-hospitalized infectious state
    p1      = 0.138; % Proportion of severe patients
    p2      = 0.5; % Death rate of severe patients 
    
    [S,E,I1,I2,H,R,D] = deal(X(1),X(2),X(3),X(4),X(5),X(6),X(7));
    
    Sdot    = -r*beta/N*S*(I1+I2+c*H);
    Edot    = r*beta/N*S*(I1+I2+c*H) - gamma1*E;
    I1dot   = gamma1*E - gamma2*I1;
    I2dot   = gamma2*(1-p1)*I1 - gamma4*I2;
    Hdot    = gamma2*p1*I1 - gamma3*H;
    Rdot    = gamma4*I2 + gamma3*(1-p2)*H;
    Ddot    = gamma3*p2*H;
    X_dot   = [Sdot,Edot,I1dot,I2dot,Hdot,Rdot,Ddot]';
end

