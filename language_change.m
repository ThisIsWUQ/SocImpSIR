% this MATLAB code is adapted from Alexandre Celestino's Epidemic-Models/SIR
% https://github.com/AlexandreCelestino/Epidemic-Models/tree/main/SIR

%% define the set of SIR equations with beta adjusted with Social Impact Theory
function dydt = SocImp_SIR(t, Y, beta0, gamma, status0, intimacy0, k)
    
    Susceptible = Y(1);
    Infected = Y(2);
    Recovered = Y(3);

    % define social status: S(t) <- str0 + f(t) where f(t), a change in the influence over time 't' is optional
    status = status0 * exp(-0.05 * t);

    % define social intimacy: I(t) <- imd0 + g(t) where g(t) is optional
    intimacy = intimacy0;

    % formulate beta informed by Social Impact Theory
    beta = beta0 * status * intimacy * (Infected^k);

    % SIR equations
    dSdt = (-1) * beta * Susceptible * Infected;
    dIdt = (beta * Susceptible * Infected) - (gamma * Infected);
    dRdt = gamma * Infected;

    dydt = [dSdt; dIdt; dRdt];

end

%% determine intial conditions

N = 1000; % population size

Infected = 25; % given that there are 25 people who learn a slang word
Susceptible = N - Infected; % the rest of the people probably adopt the new word
Recovered = 0; % given that, at the initial state, no one is giving up using the term

SIR = [Susceptible Infected Recovered];
%% determin parameters

beta0 = 0.2; % transmission rate
gamma = 0.05; % recovery rate
status0 = 0.1; % baseline influence of social status
intimacy0 = 0.025; % baseline influence of social intimacy
k = 0.25; % scaling exponential

%% determine timeframe

T = 100; % evaluation time
Tspam = [0:0.1:T]; % time interval

%% solve the equations

[T,Y] = ode45(@(t,Y) SocImp_SIR(t,Y,beta0,gamma,status0,intimacy0,k),Tspam,SIR);

S=Y(:,1); % Solution S
I=Y(:,2); % Solution I
R=Y(:,3); % Solution R

%% plot
plot(T,S,'g', 'LineWidth', 2);
hold on;
grid on;
plot(T,I,'r--', 'LineWidth', 2);
plot(T,R,'b-.', 'LineWidth', 2);
hold off;
title(['SIR model with parameters: \beta= ',num2str(beta0),', \gamma= ',num2str(gamma), ', N=',num2str(N), ', Status0=',num2str(status0)])
xlabel('Time')
ylabel('Number of Individuals')
legend('S','I','R','Location','best')
