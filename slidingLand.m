clc;clear;close all


%% PLANET PARAMS %%
GM = 813730.7;
a = 20000;
b = 7000;
c = 6500;
rRate = 3.3118/10000;

%% LANDER IC %%
r0 = [25 -1 -1]'*1000;
v0 = [0 0 1]';
m_wet = 500;
t0 = 0;

%% LANDER TC %%
rd = [a 0 0]';
vd = [0 0 0]';
tf = 7778.3;
n = 0.5;
tf2 = n*tf;

%% SMC IC %%
lam = 1.5; %{1.5,2,4,8}
s10 = r0 - rd;
s1_dot0 = v0 - vd;
s20 = s1_dot0 + lam*s10/tf;
phi = abs(s20)/tf2;

%% DE solve %%
N = 10000;
tSpan = linspace(t0,tf,N);
init = [r0' v0' m_wet']';

[t,x] = ode45(@odefun, tSpan, init);

%% PLOTS %%

[X,Y,Z] = ellipsoid(0,0,0,a,b,c);

% TRAJECTORY %
figure(1)
plot3(x(:,1),x(:,2),x(:,3), 'LineWidth', 1.5);
hold on
plot3(x(1,1),x(1,2),x(1,3),'ks')
hold on
plot3(x(N,1),x(N,2),x(N,3),'ro')
hold on
surf(X,Y,Z);
grid on

% POSITION %
figure(2)
subplot(3,1,1)
plot(t,x(:,1))
grid on

subplot(3,1,2)
plot(t,x(:,2))
grid on

subplot(3,1,3)
plot(t,x(:,3))
grid on

% VELOCITY %
figure(3)
subplot(3,1,1)
plot(t,x(:,4))
grid on

subplot(3,1,2)
plot(t,x(:,5))
grid on

subplot(3,1,3)
plot(t,x(:,6))
grid on

% MASS %
figure(4)
plot(t,x(:,7))
grid on


%% DE FUNC %%
function dx = odefun(t,x)
    dx = zeros(7,1);
    
    %%% PLANET PARAMS %%%
    GM = 813730.7;
    a = 20000;
    b = 7000;
    c = 6500;
    rRate = 3.3118/10000;
    
    %%% LANDER IC %%%
    r0 = [25 -1 -1]'*1000;
    v0 = [0 0 1]';
    
    %%% LANDER TC %%%
    rd = [a 0 0]';  %%% desired position
    vd = [0 0 0]';  %%% desired velocity. generally, [0 0 0] is too restrictive,
                    ... and vz ~~ 1-2m/s on earth
    tf = 7778.3;  %%% fixed terminal time. "CAN LOOK INTO FREE TERMINAL TIME PROBLEMS"
    n = 0.5;  %%% rate at which 2nd sliding surface converges, must n \in [0,1]
    tf2 = n*tf;
    
    %%% GRAVITY AND THRUST PARAMS %%%
    g = grav(GM, a, b, c, x(1), x(2), x(3));  %%% gravity
    
    g0 = grav(GM, a, b, c, r0(1), r0(2), r0(3));  %%% gravity at desired position
    Isp = 300;  %%% specific impulse

    %%% SLIDING MODE CONTROL %%%
    v = [x(4) x(5) x(6)]';  %%% current velocity
    r = [x(1) x(2) x(3)]';  %%% current position
    
    lam = 1.5;  %%% must be more than 1
    
    %%% SMC IC %%%
    s10 = (r0 - rd);
    s1_dot0 = v0 - vd;
    s1d_dot0 = -lam*s10/tf;
    s20 = s1_dot0 - s1d_dot0;
    
    s1 = r - rd;  %%% first sliding surface to bring altitude to desired altitude(rd).
    s1_dot = v - vd;  %%% A
    s1d_dot = -lam*s1/(tf - t);  %%% B
    s2 = s1_dot - s1d_dot;  %%% second sliding surface for velocity profile(A) to track desired velocity profile(B)
    
    phi = abs(s20)/tf2;  %%% must be positive always
    
    %%% THRUST GENERATION %%%
    T(1) = -x(7)*(2*rRate*x(5) + x(1)*rRate^2 + g(1) + (lam*s1_dot(1)*(tf - t) + lam*s1(1))/((tf - t)^2) + phi(1)*sign(s2(1)));
    T(2) = -x(7)*(-2*rRate*x(4) + x(2)*rRate^2 + g(2) + (lam*s1_dot(2)*(tf - t) + lam*s1(2))/((tf - t)^2) + phi(2)*sign(s2(2)));
    T(3) = -x(7)*(g(3) + (lam*s1_dot(3)*(tf - t) + lam*s1(3))/((tf - t)^2) + phi(3)*sign(s2(3)));

    %%% ODE EQNS %%%
    dx(1) = x(4);%x
    dx(2) = x(5);%y
    dx(3) = x(6);%z
    dx(4) = 2*rRate*x(5) + x(1)*rRate^2 + g(1) + T(1)/x(7);%vx
    dx(5) = -2*rRate*x(4) + x(2)*rRate^2 + g(2) + T(2)/x(7);%vy
    dx(6) = g(3) + T(3)/x(7);%vz
    dx(7) = -norm(T)/(Isp*norm(g0));%m

end

