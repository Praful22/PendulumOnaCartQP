% A generic template for computing Cart-Pole Dynamics and implementing QP
% control.

syms l M m g real  % Parameters/Parametric variables of the dynamic system
syms x theta dx dtheta u real % State and input variables of the system

[A,B,ddq] = computeRobotDynamicsSymbolic();

% Substituting the system parameters with physical measurement values
% Mass of cart(M) = 1.0 [kg]
% Mass of pole(m) = 0.2 [kg]
% length of pole(l) = 0.3 [m]
% acceleration due to gravity(g) = 9.81 [m/s^2]

A = subs(A,[M,m,l,g],[1.0,0.2,0.3,9.81]);

B = subs(B,[M,m,l,g],[1.0,0.2,0.3,9.81]);

ddq = subs(ddq,[M,m,l,g],[1.0,0.2,0.3,9.81]);

A_fun = matlabFunction(A);

A_Mat = A_fun(0,0,0); % Desired States dtheta, dx, x, theta

%disp(A_afterComputation);

B_fun = matlabFunction(B);

B_vec = B_fun(0);

ddq_fun = matlabFunction(ddq);

% LQR control
% Compute Gain
Qa = diag([10,10,10,10]); 

Ra = 1;

global params

params.Ka = lqr(A_Mat,B_vec,Qa,Ra);

x0=[0;pi/6;0;0]; % x0 is the intial state of the system

tspan= 0:0.01:2; % simulation time

[t,x]=ode45(@(t,x) sys_dynamics(t,x,ddq_fun),tspan,x0);

% recreate control inputs
 for i=1:length(t)
     u(:,i)=controller(t(i),x(i,:)');
 end

 % Plot x(t) and theta(t)
figure;
subplot(3, 1, 1);
plot(t, x(:,1));
xlabel('Time (s)');
ylabel('x (m)');
title('Cart Position over Time with QP controller.');
grid on;

subplot(3, 1, 2);
plot(t, x(:,2));
xlabel('Time (s)');
ylabel('\theta (rad)');
title('Pendulum Angle over Time with QP controller.');
grid on;

subplot(3,1,3);
plot(t,u);
title('control input');
title('Control input with QP controller.');
xlabel('Time (s)');
ylabel('u(N)')
grid on;

% Animate the motion

%animate_cart_pendulum(t, x(:, 1), x(:, 2), 0.3);

function dx=sys_dynamics(t,x,ddq_fun)

u=controller(t,x);

dq = x(3:4);

ddq = ddq_fun(x(4),x(2),u);

dx= [dq;ddq];
end


function u=controller(t,x)
% QP CONTROL
global params
%Control input
    umax = 10;

    umin = -10;

    H = 1;

    f = params.Ka * x; % f = - ud = - (-Ka*x) = Ka*x
    
    Aineq = [1;-1];

    bineq = [umax;-umin];
    options = optimoptions('quadprog','Display','off');

    u = quadprog(H,f,Aineq,bineq,[],[],[],[],[],options);
    
end


function [A,B,ddq] = computeRobotDynamicsSymbolic()

% Symbolically note down the variables

syms l M m g real % Parameters/Parametric variables of the dynamic system
syms x theta dx dtheta u real % State and input variables of the system

% State vectors; State of the dynamic system: X = [q;dq]
% =[x;theta;dx;dtheta]

q = [x; theta];
dq = [dx; dtheta];
u_vec = [u;0];

% x_cm_pole = x - l*sin(theta); y_cm_pole = l*cos(theta); 

vxpole = dx - l * cos(theta) * dtheta;

vypole = -l*sin(theta)*dtheta;

% Kinetic Energy
T = (0.5 * M * dx^2) + (0.5 * m * (vxpole^2+vypole^2)); % K.E of cart + K.E. of pole

% Potential Energy(P.E)
% P.E of cart = 0, so U = 0 + P.E of pole 
U = m*g*l*cos(theta);

% Lagrangian
L_q_dq = simplify(T - U);

f_q_dq = simplify(jacobian(L_q_dq,dq));

D_q = simplify(jacobian(f_q_dq,dq)); % Inertia Matrix % Say Internal components affecting states.

N_q_dq = simplify(jacobian(f_q_dq,q)*dq-jacobian(L_q_dq,q)');% Say External components affecting states.

ddq = simplify(D_q\(u_vec-N_q_dq));

f_x_u = [dq;ddq];

X = [q;dq];

A = simplify(jacobian(f_x_u,X));

B = simplify(jacobian(f_x_u,u));

end