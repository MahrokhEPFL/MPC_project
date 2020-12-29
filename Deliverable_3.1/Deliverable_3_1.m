%Deliverable 3.1
close all
clear all
addpath('C:\gurobi910\win64\matlab')


set(0, 'DefaultLineLineWidth', 1)

Ts=1/5;
quad=Quad(Ts);
[xs,us]=quad.trim();
sys=quad.linearize(xs,us);
[sys_x,sys_y,sys_z,sys_yaw]=quad.decompose(sys,xs,us);

tol=1e-3; %tolerance for convergence

%% Yaw

% Design MPC controller
mpc_yaw=MPC_Control_yaw(sys_yaw,Ts);

sol_yaw.x(:,1)=[0;deg2rad(45)]; % x0

% Simulate the closed-loop system

i=1;
while norm(sol_yaw.x(:,end))>tol
    
% Get the optimal control input
sol_yaw.u(i)=mpc_yaw.get_u(sol_yaw.x(:,i));

% Apply the optimal control input to the system
sol_yaw.x(:,i+1)=mpc_yaw.A*sol_yaw.x(:,i)+mpc_yaw.B*sol_yaw.u(i);
i=i+1;
end

settling_step=size(sol_yaw.x,2)-1;
t=Ts*[0:settling_step]; %discrete times for which the system is simulated

figure
plot(t,sol_yaw.x(1,:),'--.k','markersize', 18)
xlabel('$t$ [s]','Interpreter','latex')
ylabel('$\dot{\gamma}$','Interpreter','latex')
set(gcf,'color','white')
set(gca,'fontsize',14)

figure
plot(t,sol_yaw.x(2,:),'--.b','markersize', 18)
xlabel('$t$ [s]','Interpreter','latex')
ylabel('$\gamma$','Interpreter','latex')
set(gcf,'color','white')
set(gca,'fontsize',14)

figure
stairs(t(1:end-1),sol_yaw.u,'k')
yline(0.2,'--r','$M_{\gamma} max$','LineWidth',2,'Interpreter','latex')
yline(-0.2,'--r','$M_{\gamma} min$','LineWidth',2,'Interpreter','latex')
xlabel('$t$ [s]','Interpreter','latex')
ylabel('$M_\gamma$','Interpreter','latex')
ylim([-0.25 0.25])
set(gcf,'color','white')
set(gca,'fontsize',14)

% Plot invariant set
Q=10*eye(2);
R=1;
sysMPT=LTISystem('A',mpc_yaw.A,'B',mpc_yaw.B);
sysMPT.u.max=0.2;
sysMPT.u.min=-0.2;
sysMPT.x.penalty=QuadFunction(Q);
sysMPT.u.penalty=QuadFunction(R);
Xf=sysMPT.LQRSet;
P=sysMPT.LQRPenalty.weight;
K=sysMPT.LQRGain;

figure
plot(Xf,'alpha',0.5,'color','r')
xlabel('$\dot{\gamma}$','Interpreter','latex')
ylabel('$\gamma$','Interpreter','latex')
set(gcf,'color','white')
set(gca,'fontsize',14)
%% X

%Design MPC controller
mpc_x=MPC_Control_x(sys_x,Ts);

%Get control inputs
sys_x_x0=[0 0 0 2]'; %x0

sol_x.states(:,1)=sys_x_x0;

% Simulate the closed-loop system
k=1;
while norm(sol_x.states(:,end))>tol
    
sol_x.u(:,k)=mpc_x.get_u(sol_x.states(:,k));

sol_x.states(:,k+1)=mpc_x.A*sol_x.states(:,k)+mpc_x.B*sol_x.u(:,k);
k=k+1;
end

green = [0.13, 0.54, 0.13];
orange = [0.8500, 0.3250, 0.0980];

beta_lim=deg2rad(2); % -2°<=beta<=2°

settling_step=size(sol_x.states,2)-1;
t=Ts*[0:settling_step]; %discrete times for which the system is simulated

figure
plot(t,sol_x.states(1,:),'--.k','markersize', 18)
xlabel('$t$ [s]','Interpreter','latex')
ylabel('$\dot{\beta}$','Interpreter','latex')
set(gcf,'color','white')
set(gca,'fontsize',14)

figure
plot(t,sol_x.states(2,:),'--.b','markersize', 18)
xlabel('$t$ [s]','Interpreter','latex')
ylabel('$\beta$','Interpreter','latex')
yline(beta_lim,'--r','$\beta_{max}$','LineWidth',2,'Interpreter','latex')
yline(-beta_lim,'--r','$\beta_{min}$','LineWidth',2,'Interpreter','latex')
set(gcf,'color','white')
set(gca,'fontsize',14)

figure
plot(t,sol_x.states(3,:),'--.','Color',green,'markersize', 18)
xlabel('$t$ [s]','Interpreter','latex')
ylabel('$\dot{x}$','Interpreter','latex')
set(gcf,'color','white')
set(gca,'fontsize',14)

figure
plot(t,sol_x.states(4,:),'--.','Color',orange,'markersize', 18)
xlabel('$t$ [s]','Interpreter','latex')
ylabel('$x$','Interpreter','latex')
set(gcf,'color','white')
set(gca,'fontsize',14)

figure
stairs(t(1:end-1),sol_x.u,'k')
yline(0.3,'--r','$M_{\beta} max$','LineWidth',2,'Interpreter','latex')
yline(-0.3,'--r','$M_{\beta} min$','LineWidth',2,'Interpreter','latex')
xlabel('$t$ [s]','Interpreter','latex')
ylabel('$M_\beta$','Interpreter','latex')
ylim([-0.35 0.35])
set(gcf,'color','white')
set(gca,'fontsize',14)


% Plot invariant set
Q=diag([5 15 5 15]);
R=1;
      
Sys=LTISystem('A',mpc_x.A,'B',mpc_x.B);
Sys.x.min=[-Inf;-beta_lim;-Inf;-Inf];  %-2°<=x2<=2°
Sys.x.max=[Inf;beta_lim;Inf;Inf];    
Sys.u.min=-0.3; %-0.3<=M_beta<=0.3
Sys.u.max=0.3;
Sys.x.penalty=QuadFunction(Q);
Sys.u.penalty=QuadFunction(R);

K=Sys.LQRGain; % Terminal controler
P=Sys.LQRPenalty.weight; %Terminal cost
Xf=Sys.LQRSet; %Terminal invariant set

figure
Xf.projection(1:2).plot('alpha',0.5,'color','r');
xlabel('$\dot{\beta}$','Interpreter','latex')
ylabel('$\beta$','Interpreter','latex')
set(gcf,'color','white')
set(gca,'fontsize',14)


figure
Xf.projection(2:3).plot('alpha',0.5,'color','r');
xlabel('$\beta$','Interpreter','latex')
ylabel('$\dot{x}$','Interpreter','latex')
set(gcf,'color','white')
set(gca,'fontsize',14)


figure
Xf.projection(3:4).plot('alpha',0.5,'color','r');
xlabel('$\dot{x}$','Interpreter','latex')
ylabel('$x$','Interpreter','latex')
set(gcf,'color','white')
set(gca,'fontsize',14)


%% Y

%Design MPC controller
mpc_y=MPC_Control_y(sys_y,Ts);

%Get control inputs
sys_y_x0=[0 0 0 2]'; %initial state

sol_y.states(:,1)=sys_y_x0;

% Simulate the closed-loop system
k=1;
while norm(sol_y.states(:,end))>tol
    
sol_y.u(:,k)=mpc_y.get_u(sol_y.states(:,k));

sol_y.states(:,k+1)=mpc_y.A*sol_y.states(:,k)+mpc_y.B*sol_y.u(:,k);
k=k+1;
end

green = [0.13, 0.54, 0.13];
orange = [0.8500, 0.3250, 0.0980];

alpha_lim=deg2rad(2); % -2°<=alpha<=2°

settling_step=size(sol_y.states,2)-1;
t=Ts*[0:settling_step]; %discrete times for which the system is simulated

figure
plot(t,sol_y.states(1,:),'--.k','markersize', 18)
xlabel('$t$ [s]','Interpreter','latex')
ylabel('$\dot{\alpha}$','Interpreter','latex')
set(gcf,'color','white')
set(gca,'fontsize',14)

figure
plot(t,sol_y.states(2,:),'--.b','markersize', 18)
xlabel('$t$ [s]','Interpreter','latex')
ylabel('$\alpha$','Interpreter','latex')
yline(alpha_lim,'--r','$\alpha_{max}$','LineWidth',2,'Interpreter','latex')
yline(-alpha_lim,'--r','$\alpha_{min}$','LineWidth',2,'Interpreter','latex')
set(gcf,'color','white')
set(gca,'fontsize',14)

figure
plot(t,sol_y.states(3,:),'--.','Color',green,'markersize', 18)
xlabel('$t$ [s]','Interpreter','latex')
ylabel('$\dot{y}$','Interpreter','latex')
set(gcf,'color','white')
set(gca,'fontsize',14)

figure
plot(t,sol_y.states(4,:),'--.','Color',orange,'markersize', 18)
xlabel('$t$ [s]','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
set(gcf,'color','white')
set(gca,'fontsize',14)

figure
stairs(t(1:end-1),sol_y.u,'k')
yline(0.3,'--r','$M_{\alpha} max$','LineWidth',2,'Interpreter','latex')
yline(-0.3,'--r','$M_{\alpha} min$','LineWidth',2,'Interpreter','latex')
xlabel('$t$ [s]','Interpreter','latex')
ylabel('$M_\beta$','Interpreter','latex')
ylim([-0.35 0.35])
set(gcf,'color','white')
set(gca,'fontsize',14)


% Plot invariant set
Q=diag([5 15 5 15]);
R=1;
      
Sys=LTISystem('A',mpc_y.A,'B',mpc_y.B);
Sys.x.min=[-Inf;-alpha_lim;-Inf;-Inf];  %-2°<=x2<=2°
Sys.x.max=[Inf;alpha_lim;Inf;Inf];    
Sys.u.min=-0.3; %-0.3<=M_alpha<=0.3
Sys.u.max=0.3;
Sys.x.penalty=QuadFunction(Q);
Sys.u.penalty=QuadFunction(R);

K=Sys.LQRGain; % Terminal controler
P=Sys.LQRPenalty.weight; %Terminal cost
Yf=Sys.LQRSet; %Terminal invariant set

figure
Yf.projection(1:2).plot('alpha',0.5,'color','b');
xlabel('$\dot{\alpha}$','Interpreter','latex')
ylabel('$\alpha$','Interpreter','latex')
hold on
set(gcf,'color','white')
set(gca,'fontsize',14)


figure
Yf.projection(2:3).plot('alpha',0.5,'color','b');
xlabel('$\alpha$','Interpreter','latex')
ylabel('$\dot{y}$','Interpreter','latex')
hold on
set(gcf,'color','white')
set(gca,'fontsize',14)


figure
Yf.projection(3:4).plot('alpha',0.5,'color','b');
xlabel('$\dot{y}$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
hold on
set(gcf,'color','white')
set(gca,'fontsize',14)


%% Z

% Design MPC controller
mpc_z=MPC_Control_z(sys_z,Ts);

sol_z.x(:,1)=[0;2]; % x0

% Simulate the closed-loop system
i=1;
while norm(sol_z.x(:,end))>tol    
sol_z.u(i)=mpc_z.get_u(sol_z.x(:,i));

sol_z.x(:,i+1)=mpc_z.A*sol_z.x(:,i)+mpc_z.B*sol_z.u(i);
i=i+1;
end

settling_step=size(sol_z.x,2)-1;
t=Ts*[0:settling_step]; %discrete times for which the system is simulated

figure
plot(t,sol_z.x(1,:),'--.k','markersize', 18)
xlabel('$t$ [s]','Interpreter','latex')
ylabel('$\dot{z}$','Interpreter','latex')
set(gcf,'color','white')
set(gca,'fontsize',14)

figure
plot(t,sol_z.x(2,:),'--.b','markersize', 18)
xlabel('$t$ [s]','Interpreter','latex')
ylabel('$z$','Interpreter','latex')
set(gcf,'color','white')
set(gca,'fontsize',14)

figure
stairs(t(1:end-1),sol_z.u,'k')
yline(0.3,'--r','$F_{max}$','LineWidth',2,'Interpreter','latex')
yline(-0.2,'--r','$F_{min}$','LineWidth',2,'Interpreter','latex')
xlabel('$t$ [s]','Interpreter','latex')
ylabel('$F$','Interpreter','latex')
ylim([-0.25 0.35])
set(gcf,'color','white')
set(gca,'fontsize',14)

% Plot invariant set
Q=diag([10 20]);
R=1;
sysMPT=LTISystem('A',mpc_z.A,'B',mpc_z.B);
sysMPT.u.max=0.3;
sysMPT.u.min=-0.2;
sysMPT.x.penalty=QuadFunction(Q);
sysMPT.u.penalty=QuadFunction(R);
Xf=sysMPT.LQRSet;
P=sysMPT.LQRPenalty.weight;

figure
plot(Xf,'alpha',0.5,'color','b')
xlabel('$\dot{z}$','Interpreter','latex')
ylabel('$z$','Interpreter','latex')
set(gcf,'color','white')
set(gca,'fontsize',14)
