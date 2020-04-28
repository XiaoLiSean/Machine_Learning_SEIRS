clear
clc
close all
dT = 1; 
T_pass = 200;
N=1100; % total people in the area
r=0.1; % average people each one contact everyday
beta_=0.5; % the probability to be infected if contacted
beta=r*beta_/N; % the probability to be infected 
% beta=(1-(1-beta_)^r)/N;
c=0.2; % the proportion of severe cannot be hospitalized
gamma1=1/5.1; % the inverse of incubation period
gamma2=1/8.5; % the inverse of the period from mild to severe
gamma3=1/35; % the inverse of the treatment time in hospital
gamma4=1/5.5; % the inverse of (tretement time - the period from mild to severe)
u_beta = 1;
u_gamma2 = 1;
p1=0.14; % proportion of severe people amoung all the patients
p2=0.5; % proportion of hospitalized patients that die
n=1; % initial infected people is about 10000
[S0,E0,I10,I20,H0,R0,D0]=deal(N-n,0,n,0,0,0,0); % the disease start from n people in an N-people city
x0=[S0,E0,I10,I20,H0,R0,D0]';
u0=[u_beta, u_gamma2, p1, p2]';
parameter_matrix = [beta, c, gamma1, gamma2, gamma3, gamma4]';

% * generate continuous results
[T_cont,X_cont] = ode45(@(t,x) SEI1I2HRD_cont(t,x,beta,c,gamma1,gamma2,gamma3,gamma4,p1,p2), [0 T_pass], x0);
Ntot = T_pass/dT; % * total steps

X_disc = zeros(length(x0), Ntot+1);
X_disc(:, 1) = x0;
T_disc = 0:dT:T_pass;
% * generate discrete results
for idx = 1:Ntot
    X_disc(:, idx+1) = SEI1I2HRD_disc(dT, X_disc(:, idx), u0, parameter_matrix);
end
%%
dS=zeros(T_pass,1);
for i=1:T_pass-1
    dS(i)=(X_disc(1,i+1)-X_disc(1,i))/dT;
end
for i=1:T_pass-1
    if dS(i+1)>dS(i) && abs(dS(i))<0.005
       T_final=i; 
       break
    end 
end


%% * plots
figure
plot(T_cont,X_cont,'LineWidth',1);
hold on;
plot(T_disc, X_disc, '--', 'LineWidth', 2)
title({['\beta = ' num2str(beta_,2) ', r = ' num2str(r,2) ', c = ' num2str(c,2) ', \gamma_{1} = ' ...
    num2str(gamma1,2) ', \gamma_{2} = ' num2str(gamma2,2) ', \gamma_{3} = ' ...
    num2str(gamma3,2) ', \gamma_{4} = ' num2str(gamma4,2) ','], ['p_{1} = ' num2str(p1,2)...
    ', p_{2} = ' num2str(p2,2)]})
legend({'S, cont','E, cont','I_{1}, cont','I_{2}, cont','H, cont','R, cont','D, cont',...
        'S, disc','E, disc','I_{1}, disc','I_{2}, disc','H, disc','R, disc','D, disc'},...
        'NumColumns', 2);
ylim([0,1200])


