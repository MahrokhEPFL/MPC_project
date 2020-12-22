function [ctrl, traj] = ctrl_NMPC(quad)

import casadi.*

opti = casadi.Opti(); % Optimization problem

N = 50; % MPC horizon [SET THIS VARIABLE]

% ---- decision variables ----
X = opti.variable(12,N+1); % state trajectory variables
U = opti.variable(4,N); % control variables (four propeller speeds)

X0 = opti.parameter(12,1); % initial state
REF = opti.parameter(4,1); % reference position [x,y,z,yaw]

%%%%%%%%%%%%%%%%%%%%%%%%
%%%% YOUR CODE HERE %%%%

% Continuous-time dynamics
f=@(x,u) quad.f(x,u);

% Input constraints 
UMIN=0;
UMAX=1.5;

% Alpha constraint 
AlphaMAX=deg2rad(2);

% Beta constraint
BetaMAX=deg2rad(2);

% Discrete-time model
Ts=1/5; % sample period
F=@(x,u) RK4_quad(x,u,Ts,f);

% Split the state into its part 
[omega,theta,vel,pos]=quad.parse_state(X);

% Objective
opti.minimize(100*sum((theta(3,:)-REF(4)).^2)+...
              100*sum((pos(1,:)-REF(1)).^2)+...
              100*sum((pos(2,:)-REF(2)).^2)+...
              100*sum((pos(3,:)-REF(3)).^2)+...
              75*omega(1,:)*omega(1,:)'+75*omega(2,:)*omega(2,:)'+75*omega(3,:)*omega(3,:)'+...
              100*theta(1,:)*theta(1,:)'+100*theta(2,:)*theta(2,:)'+...
              75*vel(1,:)*vel(1,:)'+75*vel(2,:)*vel(2,:)'+75*vel(3,:)*vel(3,:)'+...
              sum(U(1,:).^2)+sum(U(1,:).^2)+sum(U(2,:).^2)+sum(U(3,:).^2)+sum(U(4,:).^2));
 
% Dynamic constraints
for i=1:N
    opti.subject_to(X(:,i+1)==F(X(:,i),U(:,i)));
end

% Input constraints
opti.subject_to(UMIN<=U<=UMAX);

% States constraints 
%opti.subject_to(-AlphaMAX<=theta(1,:)<=AlphaMAX);
%opti.subject_to(-BetaMAX<=theta(2,:)<=BetaMAX);

% Initial conditions
opti.subject_to(X(:,1)==X0);
%%%% YOUR CODE HERE %%%%
%%%%%%%%%%%%%%%%%%%%%%%%

ctrl = @(x,ref) eval_ctrl(x, ref, opti, X0, REF, X, U);
end

function u = eval_ctrl(x, ref, opti, X0, REF, X, U)
% ---- Set the initial state and reference ----
opti.set_value(X0, x);
opti.set_value(REF, ref);

% ---- Setup solver NLP ----
ops = struct('ipopt', struct('print_level',0,'tol',1e-3), 'print_time', false);
opti.solver('ipopt', ops);

% ---- Solve the optimization problem ----
sol = opti.solve();
assert(sol.stats.success == 1, 'Error computing optimal input');

u = opti.value(U(:,1));

% Use the current solution to speed up the next optimization
opti.set_initial(sol.value_variables());
opti.set_initial(opti.lam_g, sol.value(opti.lam_g));
end

