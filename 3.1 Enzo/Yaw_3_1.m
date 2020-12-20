%% 3.1 Yaw
Ts=1/5;
quad=Quad(Ts);
[xs,us]=quad.trim();
sys=quad.linearize(xs,us);
[sys_x,sys_y,sys_z,sys_yaw]=quad.decompose(sys,xs,us);

% Design MPC controller
mpc_yaw=MPC_Control_yaw(sys_yaw,Ts);

sol_yaw.x(:,1)=[0;deg2rad(45)]; % x0

% Simulate the closed-loop system
for i=1:50 %i=40 corresponds to 8s (settling time)
% Get the optimal control input
sol_yaw.u(i)=mpc_yaw.get_u(sol_yaw.x(:,i));

% Apply the optimal control input to the system
sol_yaw.x(:,i+1)=mpc_yaw.A*sol_yaw.x(:,i)+mpc_yaw.B*sol_yaw.u(i);
end

figure
plot(Ts*[0:50],sol_yaw.x(2,:),'-ok')
hold on
plot(Ts*[0:50],sol_yaw.x(1,:),'-ob')
hold off
xlabel('$t$ [s]','Interpreter','latex')
legend('$\gamma$','$\dot{\gamma}$','Interpreter','latex')

figure
stairs(Ts*[0:49],sol_yaw.u,'k')
hold on
plot(Ts*[0:49],-0.2*ones(1,50),'-.r')
plot(Ts*[0:49],0.2*ones(1,50),'-.r')
hold off
xlabel('$t$ [s]','Interpreter','latex')
ylabel('$M_\gamma$','Interpreter','latex')
ylim([-0.25 0.25])