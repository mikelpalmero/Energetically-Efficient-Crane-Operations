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
Omega0 = sqrt(g/l0);
Gamma = (lf/l0)^(1/4);
tf = pi/w;
x = linspace(0, xmax, 1e5)';
dx = x(2)-x(1);
tt = linspace(0, tf, 1e5)';


%% Minimization to get optimal shortcut for the functional

%Seed for the minimization
% a6_seed = 0;
% a7_seed = 0;
a6_seed = -32.910731083795184;
a7_seed = 7.791864009772690;
a8_seed = 0;
% 
% %We call to the minimizing function
% [z, min]  = fminsearch(@(z) minimizeBC(z, Gamma, tf, Omega0, m, g, l0, lf), [a6_seed, a7_seed]);
[z, min]  = fminsearch(@(z) minimizeBC(z, Gamma, tf, Omega0, m, g, l0, lf), [a6_seed, a7_seed, a8_seed]);

a6 = z(1);
a7 = z(2);
% a8 = 0;
a8 = z(3);

E_exc = min;

%Tolerance options for the differential equation
options = bvpset('RelTol', 1e-6);

%solver the differential equation for this specific a6, a7
tspan = [0, tf];
init_cond = [l0 0];

% [T, X] = ode45(@(t, y) l_diff_eq(t, y, a6, a7, Gamma, tf, g, Omega0), tspan, init_cond, options);
[T, X] = ode45(@(t, y) l_diff_eq(t, y, a6, a7, a8, Gamma, tf, g, Omega0), tt, init_cond, options);
    
%Get l and l' from the solution for the given free parameters
l = X(:,1);
l_dt = X(:,2);

%We solve the dynamics of the mass in the harmonic oscillator to get x(t)

% x_ho = xmax/2*(1+cos(w*T));
% dx_ho = -xmax*w/2*sin(w*T);
% ddx_ho = -xmax*w^2/2*cos(w*T);
x_ho = xmax/2*(1+cos(w*tt));
dx_ho = -xmax*w/2*sin(w*tt);
ddx_ho = -xmax*w^2/2*cos(w*tt);

%Function of the rail (h(t))
h = sqrt((L-l).^2-(D-x_ho).^2);

% figure
% 
% % figure(1)
% subplot(3,1,1)
% plot(x_ho, h)
% title('design of the rail')
% ylabel('$h(x)$', 'Interpreter', 'latex')
% xlabel('$x$', 'Interpreter', 'latex')
% 
% % figure(2)
% subplot(3,1,2)
% plot(T, l)
% title('Distance between the pulley and the hanging mass')
% ylabel('$l (t)$', 'Interpreter', 'latex')
% xlabel('$t$', 'Interpreter', 'latex')
% 
% % figure(3)
% subplot(3,1,3)
% plot(T, l_dt)
% title('Speed of the hanging mass')
% ylabel('$\dot{l} (t)$', 'Interpreter', 'latex')
% xlabel('$t$', 'Interpreter', 'latex')


%% Interpolation of the l function and numerical derivatives

%l as a function of x
l_x = interp1(x_ho, l, x);

%Numerical derivatives
dl_x = first_deriv(l_x, dx);
dl_x_b = gradient(l_x, dx);
ddl_x = second_deriv(l_x, dx);
ddl_x_b = first_deriv(dl_x,dx);
ddl_x_c = gradient(dl_x_b, dx);

%Interpolate back as a function of T
% dl_t = interp1(x, dl_x, x_ho);
% ddl_t = interp1(x, ddl_x, x_ho);

dl_t = interp1(x, dl_x_b, x_ho);
ddl_t = interp1(x, ddl_x_c, x_ho);


%% Solve the Equations of motion

% Solve the first equation to get theta
init_cond = [pi/12, 0];
% init_cond = [pi/12, 0];

% [T_th, X_th] = ode45(@(t, y) theta_EM(t, y, l, T, dl_t, dx, g), tspan, init_cond, options);
[T_th, X_th] = ode45(@(t, y) theta_EM(t, y, l, l0, T, dl_t, dx, g), tt, init_cond, options);

%Get theta and theta' from the solution for the given free parameters
theta = X_th(:,1);
theta_dt = X_th(:,2);


%Get F_{k} from the second second equation of motion

Fr = -0.02*M*dx_ho;
theta_int = interp1(T_th, theta, T);
theta_int_dt = interp1(T_th, theta_dt, T);
% theta_int = theta;
% theta_int_dt = theta_dt;

Fk = Fk_EM(l, l0, dl_t, ddl_t, theta_int, theta_int_dt, m, M, g, Fr, w, xmax, x_ho, dx_ho, ddx_ho);


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

P = Fk.*dx_ho;
P_plus = max(P, 0);
P_minus = -max(-P, 0);

n = 0;

E_cons = trapz(tt, P_plus) + n*trapz(tt, P_minus);



%% Plot some figures


figure

subplot(3,2,1)
plot(x_ho, h)
title('design of the rail')
ylabel('$h(x)$', 'Interpreter', 'latex')
xlabel('$x$', 'Interpreter', 'latex')

subplot(3,2,3)
plot(T, l)
title('Distance between the pulley and the hanging mass')
ylabel('$l (t)$', 'Interpreter', 'latex')
xlabel('$t$', 'Interpreter', 'latex')

subplot(3,2,5)
plot(T, l_dt)
title('Speed of the hanging mass')
ylabel('$\dot{l} (t)$', 'Interpreter', 'latex')
xlabel('$t$', 'Interpreter', 'latex')


subplot(3,2,2)
plot(tt, theta_int)
title('Angle')
ylabel('$\theta (t)$', 'Interpreter', 'latex')
xlabel('$t$', 'Interpreter', 'latex')

% figure(2)
subplot(3,2,4)
plot(tt, theta_int_dt)
title('Derivative of the angle')
ylabel('$\dot{\theta} (t)$', 'Interpreter', 'latex')
xlabel('$t$', 'Interpreter', 'latex')

% figure(3)
subplot(3,2,6)
plot(tt, Fk)
title('Kick Force')
ylabel('$F_k (t)$', 'Interpreter', 'latex')
xlabel('$t$', 'Interpreter', 'latex')

sgtitle('$E_{exc}$ Optimization; Poly ansatz; $\omega = 0.2$', 'Interpreter', 'latex')


%% Function definitions
%Function for the minimization

function min = minimizeBC(z, Gamma, tf, Omega0, m, g, l0, lf)
% [min, l, l_dt, T] = minimizeBC(z, Gamma, tf, Omega0, m, g, l0, lf)

    %initialize the free parameters
    a6 = z(1);
    a7 = z(2);
    a8 = z(3);
%     a8 = 0;


    %Tolerance options for the differential equation
    options = bvpset('RelTol',10^(-6));
    
    %solver the differential equation for this specific a6, a7
    tspan = [0, tf];
    init_cond = [l0 0];
    [T, X] = ode45(@(t, y) l_diff_eq(t, y, a6, a7, a8, Gamma, tf, g, Omega0), tspan, init_cond, options);
%     [T, X] = ode45(@(t, y) l_diff_eq(t, y, a6, a7, Gamma, tf, g, Omega0), tspan, init_cond);
    % [T, X] = ode45(@(t, y) l_diff_eq(t, y, a6, a7, Gamma, tf, g, Omega0), tspan, init_cond, options);
    
    
    %Get l and l' from the solution for the given free parameters
    l = X(:,1);
    l_dt = X(:,2);
    
    %Function to be minimized. Energy difference with respect to the target
    %configuration
    min = m*g*abs(l(end)-lf) + 1/2*m*l_dt(end)^2;
    
%     sprintf('a6 = %0.2f, a7 = %0.2f, a8 = %0.2f, min = %0.2s ', a6, a7, a8, min)

end


%Function for the Differential Equation that designs the shortcut for l

function z = l_diff_eq(t, y, a6, a7, a8, Gamma, tf, g, Omega0)


    %2 free parameters
    % rho = 1 + (-10 - a6 - 3*a7 + 10*Gamma)*t^3/tf^3 +...
    %      (15 + 3*a6 + 8*a7 - 15*Gamma)*t^4/tf^4 -...
    %      3*(2 + a6 + 2*a7 - 2*Gamma)*t^5/tf^5 + a6*t^6/tf^6 + a7*t^7/tf^7;
    % rho2 = 6*(-10 - a6 - 3*a7 + 10*Gamma)*t/tf^3 *...
    %        12*(15 + 3*a6 + 8*a7 - 15*Gamma)*t^2/tf^4 -...
    %        60*(2 + a6 + 2*a7 - 2*Gamma)*t^3/tf^5 + 30*a6*t^4/tf^6 + 42*a7*t^5/tf^7;

    %3 free parameters
    rho = 1 + (-10 - a6 - 3*a7 - 6*a8 + 10*Gamma)*t^3/tf^3 +...
         (15 + 3*a6 + 8*a7 + 15*a8 - 15*Gamma)*t^4/tf^4 +...
         (-6 - 3*a6 - 6*a7 - 10*a8 + 6*Gamma)*t^5/tf^5 +...
         a6*t^6/tf^6 + a7*t^7/tf^7 + a8*t^8/tf^8;
    rho2 = 6*(-10 - a6 - 3*a7 - 6*a8 + 10*Gamma)*t/tf^3 +...
           12*(15 + 3*a6 + 8*a7 + 15*a8 - 15*Gamma)*t^2/tf^4 +...
           20*(-6 - 3*a6 - 6*a7 - 10*a8 + 6*Gamma)*t^3/tf^5 +...
           30*a6*t^4/tf^6 + 42*a7*t^5/tf^7 + 56*a8*t^6/tf^8;


    z = zeros(2,1);
    z(1) = y(2);
    z(2) = g-(Omega0^2/rho^4-rho2/rho)*y(1);
    
end


%Function for the Differential Equation that solves the first equation of motion

function z = theta_EM(t, y, l, l0, T, dl_t, dx, g)

    l = interp1(T, l, t);
    dl = interp1(T, dl_t, t);


    z = zeros(2,1);
    z(1) = y(2);
    % z(2) = g/l*sin(y(1))-2*dl/l.^2*dx*y(2);
    % z(2) = -g*(l-l0)/l.^2*sin(y(1))-2*dl./l*dx*y(2);
    z(2) = -g./l*sin(y(1))-2*dl./l*dx*y(2);
    
end


%Function for the second equation of motion

function Fk = Fk_EM(l, l0, dl_t, ddl_t, theta, dtheta, m, M, g, Fr, w, xmax, x, dx, ddx)


    % Fk = m*dl_t.*ddl_t.*dx.^2 + m*dl_t.^2.*ddx - m*l.*dl_t.*dtheta.^2 +...
    %      m*g*dl_t.*cos(theta) + M*ddx + M*w^2*(x-xmax/2) - Fr;
    % Fk = m*dl_t.*ddl_t.*dx.^2 + (m*dl_t.^2 + M).*ddx - m*l.*dl_t.*dtheta.^2 -...
    %      m*g*(dl_t-l0).*cos(theta) + M*w^2*(x-xmax/2) - Fr;
    Fk = m*dl_t.*ddl_t.*dx.^2 + (m*dl_t.^2 + M).*ddx - m*l.*dl_t.*dtheta.^2 -...
         m*g*(dl_t).*cos(theta) + M*w^2*(x-xmax/2) - Fr;
    
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
