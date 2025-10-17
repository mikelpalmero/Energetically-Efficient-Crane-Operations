%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Grandfather-Clock optimization (2D)                       %
%------------------------------------------------------------------------------%
%       Mikel Palmero                      ||      7/2/2020    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear



%% Parameters

m = 5;
M = 1000;
w = 0.2;
L = 20;
D = 12;
l0 = 5.5;
lf = 0.5;
g = 9.8;
xmax = 8;
mu = 0.02;
n = 0;

Omega0 = sqrt(g/l0);
Gamma = (lf/l0)^(1/4);
tf = pi/w;
x = linspace(0, xmax, 1e5)';
dx = x(2)-x(1);
tt = linspace(0, tf, 1e5)';


%% Minimization to get optimal shortcut for the functional

%Seed for the minimization
a3_seed = 0;
a4_seed = 0;
% a6_seed = -32.910731083795184;
% a7_seed = 7.791864009772690;
a5_seed = 0;
% 
% %We call to the minimizing function
% [z, min]  = fminsearch(@(z) minimizeBC(z, Gamma, tf, Omega0, m, g, l0, lf), [a6_seed, a7_seed]);
[z, min]  = fminsearch(@(z) minimizeBC(z, Gamma, tf, Omega0, m, g, l0, lf), [a3_seed, a4_seed, a5_seed]);

a6 = z(1);
a7 = z(2);
a8 = z(3);

%Tolerance options for the differential equation
options = bvpset('RelTol',10^(-6));

%solver the differential equation for this specific a6, a7
tspan = [0, tf];
init_cond = [l0 0];

% [T, X] = ode45(@(t, y) l_diff_eq(t, y, a6, a7, Gamma, tf, g, Omega0), tspan, init_cond, options);
[T, X] = ode45(@(t, y) l_diff_eq(t, y, a6, a7, a8, Gamma, tf, g, Omega0), tt, init_cond, options);
    
%Get l and l' from the solution for the given free parameters
l = X(:,1);
l_dt = X(:,2);
l_ddt = gradient(l_dt, tt);



%% Solve the Equations of motion

% Solve the first equation to get theta
%init_cond = [pi/12, 0];
init_cond = [pi/9, 0];

% [T_th, X_th] = ode45(@(t, y) theta_EM(t, y, l, T, dl_t, dx, g), tspan, init_cond, options);
[T_th, X_th] = ode45(@(t, y) theta_EM(t, y, l, tt, l_dt, dx, g), tt, init_cond, options);

%Get theta and theta' from the solution for the given free parameters
theta = X_th(:,1);
theta_dt = X_th(:,2);


%Get F_{k} from the second second equation of motion


% Fr = -mu*M*dx_ho;


Fk = Fk_EM(l, l_dt, l_ddt, theta, theta_dt, m, M, g, mu);


% figure
% subplot(3,1,1)
% plot(tt, theta)
% title('Angle')
% ylabel('$\theta (t)$', 'Interpreter', 'latex')
% xlabel('$t$', 'Interpreter', 'latex')
% 
% % figure(2)
% subplot(3,1,2)
% plot(tt, theta_dt)
% title('Derivative of the angle')
% ylabel('$\dot{\theta} (t)$', 'Interpreter', 'latex')
% xlabel('$t$', 'Interpreter', 'latex')
% 
% % figure(3)
% subplot(3,1,3)
% plot(tt, Fk)
% title('Kick Force')
% ylabel('$F_k (t)$', 'Interpreter', 'latex')
% xlabel('$t$', 'Interpreter', 'latex')

%% Energy Consumption

% Power of the kick force

P = Fk.*l_dt;
P_plus = max(P, 0);
P_minus = -max(-P, 0);


E_cons = trapz(tt, P_plus) + n*trapz(tt, P_minus)
E_exc = m*g*abs(l(end)*cos(theta(end))-lf) + 1/2*m*(l_dt(end)^2 + l(end)^2*theta_dt(end)^2);

E_p = m*g*(l(1)-l(end));
(E_cons-E_p)/E_p*100

%% Plot some figures

% figure
% 
% % subplot(3,2,1)
% % plot(x_ho, h)
% % title('design of the rail')
% % ylabel('$h(x)$', 'Interpreter', 'latex')
% % xlabel('$x$', 'Interpreter', 'latex')
% 
% subplot(3,2,3)
% plot(T, l)
% title('Distance between the pulley and the hanging mass')
% ylabel('$l (t)$', 'Interpreter', 'latex')
% xlabel('$t$', 'Interpreter', 'latex')
% 
% subplot(3,2,5)
% plot(T, l_dt)
% title('Speed of the hanging mass')
% ylabel('$\dot{l} (t)$', 'Interpreter', 'latex')
% xlabel('$t$', 'Interpreter', 'latex')
% 
% 
% subplot(3,2,2)
% plot(tt, theta)
% title('Angle')
% ylabel('$\theta (t)$', 'Interpreter', 'latex')
% xlabel('$t$', 'Interpreter', 'latex')
% 
% % figure(2)
% subplot(3,2,4)
% plot(tt, theta_dt)
% title('Derivative of the angle')
% ylabel('$\dot{\theta} (t)$', 'Interpreter', 'latex')
% xlabel('$t$', 'Interpreter', 'latex')
% 
% % figure(3)
% subplot(3,2,6)
plot(tt, Fk)
title('Kick Force')
ylabel('$F_k (t)$', 'Interpreter', 'latex')
xlabel('$t$', 'Interpreter', 'latex')
% 
% sgtitle('Benchmark; $\cos$ ansatz (3 param); $\mu = 0.02$', 'Interpreter', 'latex')


%% Function definitions
%Function for the minimization

function min = minimizeBC(z, Gamma, tf, Omega0, m, g, l0, lf)
% [min, l, l_dt, T] = minimizeBC(z, Gamma, tf, Omega0, m, g, l0, lf)

    %initialize the free parameters
    a3 = z(1);
    a4 = z(2);
    a5 = z(3);
    % a8 = 0;


    %Tolerance options for the differential equation
    options = bvpset('RelTol',10^(-6));
    
    %solver the differential equation for this specific a6, a7
    tspan = [0, tf];
    init_cond = [l0 0];
    [T, X] = ode45(@(t, y) l_diff_eq(t, y, a3, a4, a5, Gamma, tf, g, Omega0), tspan, init_cond, options);
    % [T, X] = ode45(@(t, y) l_diff_eq(t, y, a6, a7, Gamma, tf, g, Omega0), tspan, init_cond);
    % [T, X] = ode45(@(t, y) l_diff_eq(t, y, a6, a7, Gamma, tf, g, Omega0), tspan, init_cond, options);
    
    
    %Get l and l' from the solution for the given free parameters
    l = X(:,1);
    l_dt = X(:,2);
    
    %Function to be minimized. Energy difference with respect to the target
    %configuration
    min = m*g*abs(l(end)-lf) + 1/2*m*l_dt(end)^2;
    
    % sprintf('a3 = %0.4f, a4 = %0.4f, a5 = %0.4f, min = %0.2s ', a3, a4, a5, min)

end


%Function for the Differential Equation that designs the shortcut for l

function z = l_diff_eq(t, y, a3, a4, a5, Gamma, tf, g, Omega0)


    %3 free parameters
    rho = (1 + Gamma)/2 +... 
          1/16*(9 + 32*a3 + 80*a4 + 144*a5 - 9*Gamma)*cos(pi*t/tf) +... 
          1/16*(-1 - 48*a3 - 96*a4 - 160*a5 + Gamma)*cos(3*pi*t/tf) +... 
          a3*cos(5*pi*t/tf) + a4*cos(7*pi*t/tf) + a5*cos(9*pi*t/tf);
    rho2 = -pi^2/(16*tf^2)*(9 + 32*a3 + 80*a4 + 144*a5 - 9*Gamma)*cos(pi*t/tf) +... 
          (-9*pi^2)/(16*tf^2)*(-1 - 48*a3 - 96*a4 - 160*a5 + Gamma)*cos(3*pi*t/tf) +... 
          (-25*pi^2)/tf^2*a3*cos(5*pi*t/tf) + (-49*pi^2)/tf^2*a4*cos(7*pi*t/tf) + (-81*pi^2)/tf^2*a5*cos(9*pi*t/tf);


    z = zeros(2,1);
    z(1) = y(2);
    z(2) = g-(Omega0^2/rho^4-rho2/rho)*y(1);
    
end


%Function for the Differential Equation that solves the first equation of motion

function z = theta_EM(t, y, l, T, dl_t, dx, g)

    l = interp1(T, l, t);
    dl = interp1(T, dl_t, t);


    z = zeros(2,1);
    z(1) = y(2);
    % z(2) = g/l*sin(y(1))-2*dl/l.^2*dx*y(2);
    % z(2) = -g*(l-l0)/l.^2*sin(y(1))-2*dl./l*dx*y(2);
    z(2) = -g./l*sin(y(1))-2*dl./l*dx*y(2);
    
end


%Function for the second equation of motion

function Fk = Fk_EM(l, dl_t, ddl_t, theta, dtheta, m, M, g, mu)

    Fk = (m+M).*ddl_t - m*l.*dtheta - m*g.*cos(theta) + mu*M.*dl_t;
    
end


%3 point first derivative

function df = first_deriv(f, dx)

    Df0 = (-3*f(1)+4*f(2)-f(3))/(2*dx);
    Dff = (f(end-2)-4*f(end-1)+3*f(end))/(2*dx);
    Df1 = (-f(1:end-2)+f(3:end))/(2*dx);
    df = [Df0; Df1; Dff];

end


%5 point second derivative

function ddf = second_deriv(f, dx)

    ddf1 = (35*f(1)-104*f(2)+114*f(3)-56*f(4)+11*f(5))/(12*dx^2);
    ddf2 = (11*f(1)-20*f(2)+6*f(3)+4*f(4)-f(5))/(12*dx^2);
    ddfint = (-f(1:end-4)+16*f(2:end-3)-30*f(3:end-2)+16*f(4:end-1)-f(5:end))/(12*dx^2);
    ddfend1 = (-f(end-4)+4*f(end-3)+6*f(end-2)-20*f(end-1)+11*f(end))/(12*dx^2);
    ddfend2 = (11*f(end-4)-56*f(end-3)+114*f(end-2)-104*f(end-1)+35*f(end))/(12*dx^2);

    ddf = [ddf1; ddf2; ddfint; ddfend1; ddfend2];  
    
end
