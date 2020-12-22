%Deliverable 3.2
close all
clear 

set(0, 'DefaultLineLineWidth', 1)

Ts=1/5;
quad=Quad(Ts);
[xs,us]=quad.trim();
sys=quad.linearize(xs,us);
[sys_x,sys_y,sys_z,sys_yaw]=quad.decompose(sys,xs,us);

%% Design MPC controllers
fprintf('\n\nTuning controllers...\n')
fprintf('X direction:\n')
mpc_x = MPC_Control_x(sys_x, Ts);
fprintf('\nY direction:\n')
mpc_y = MPC_Control_y(sys_y, Ts);
fprintf('\nZ direction:\n')
mpc_z = MPC_Control_z(sys_z, Ts);
fprintf('\nYaw:\n')
mpc_yaw = MPC_Control_yaw(sys_yaw, Ts);

%% Solve problem
ref = [-2, -2, -2, deg2rad(45)];

% Simulating the closed-loop system
subsystems = ["x", "y", "z", "yaw"];    % subsystems names
tol = 1e-3;                             % convergence criterion
settling_time = zeros(4, 1);            % settling time of each subsystem

% simulate for each subsystem
for sys_num = 1:4
    fprintf(strcat('\n\nSimulating subsystem ', {' '}, subsystems(sys_num)));
    
    % select subsystem
    switch sys_num
        case 1
            subsys = sys_x;
        case 2
            subsys = sys_y;
        case 3
            subsys = sys_z;
        case 4
            subsys = sys_yaw;
        otherwise
            fprintf('\n[Error] subsystem number out of range.');
    end
    
    % set initial condition
    x0 = zeros(size(subsys.B));
    clear sol
    sol.x(:,1) = x0;
    
    % select controller
    switch sys_num
        case 1
            mpc_cont = mpc_x;
        case 2
            mpc_cont = mpc_y;
        case 3
            mpc_cont = mpc_z;
        case 4
            mpc_cont = mpc_yaw;
        otherwise
            fprintf('\n[Error] subsystem number out of range.')
    end
    
    % simulate at each time step
    i = 1;
    maxiter = 100;
    try
        % Simulate until convergence
        xs = [zeros(length(subsys.A)-1,1);ref(sys_num)];
        while norm(sol.x(:,end) - xs)>tol
            % Solve MPC problem for current state
            sol.u(:,i) = mpc_cont.get_u(sol.x(:,i), ref(sys_num));   
            % Apply the optimal input to the system
            sol.x(:,i+1) = mpc_cont.A*sol.x(:,i) + mpc_cont.B*sol.u(:,i);
            % Go to next time step
            i = i + 1;
            % Maximum iterations limit
            if(i>maxiter)
                fprintf('\nDid not converge\n');
                break
            end
        end
    catch
        error('---> Initial state is outside the feasible set <---\n');
    end
    
    % find settling time
    settled = prod(abs(sol.x - xs)<=0.05,1);% settled[i] = 1 iff all states are settled at index i
    settled = find(settled);                % indexes when all states are settled
    settling_time(sys_num)=Ts*(settled(1)-1);% settling time in seconds
    msg = strcat("\nSettling time for this subsystem = ", {' '}, num2str(settling_time(sys_num)), "s.");
    fprintf(msg);
    
    % save results
    switch sys_num
        case 1
            solution.x = sol;
        case 2
            solution.y = sol;
        case 3
            solution.z = sol;
        case 4
            solution.yaw = sol;
        otherwise
            fprintf('[Error] subsystem number out of range.');
    end
    clear sol
end

fprintf(strcat("\n\nMaximum settling time of all subsystems = ", {' '}, num2str(max(settling_time)), "s.\n"));

%% Plot
% define input names
input_names = ["$M_\alpha$", "$M_\beta$", "$F$", "$M_\gamma$"];
sys_x.StateName(1) = {'vel_{pitch}'};
sys_yaw.StateName(1) = {'vel_{yaw}'};
% units?

% produce two figures for each subsystem
for sys_num = 1:4
    % retrive solution
    switch sys_num
        case 1
            sol = solution.x;
            subsys = sys_x;
        case 2
            sol = solution.y;
            subsys = sys_y;
        case 3
            sol = solution.z;
            subsys = sys_z;
        case 4
            sol = solution.yaw;
            subsys = sys_yaw;
        otherwise
            fprintf('[Error] subsystem number out of range.');
    end
    xs = [zeros(length(subsys.A)-1,1);ref(sys_num)];
    
    % get input bounds
    [umin, umax] = get_input_bounds(sys_num);
    
    % solution properties
    num_states = size(sol.x, 1);
    num_sample = size(sol.x, 2);
    
    % constants
    times = 0: Ts: Ts * (num_sample-1);
    o = ones(1, num_sample);
    
    % produce figures for each state
    fig_title = strcat("States Evolution for Subsystem", {' '}, subsystems(sys_num));
    figure('Name', fig_title,'NumberTitle','off');
    hold on; grid on;
    for fig_num = 1 : num_states
        subplot(num_states, 1, fig_num);
        hold on; grid on;
        plot(times, sol.x(fig_num, :),'-ok','markersize',5,'linewidth',2);
        plot(times, xs(fig_num)-0.05*o,'c','linewidth',1)
        plot(times, xs(fig_num)+0.05*o,'c','linewidth',1)
        xline(settling_time(sys_num), '-.r','linewidth',2) 
        legend('state', '$-0.05$ error bound', '$+0.05$ error bound', 'settling time' ,'Interpreter','latex') 
        ylabel(strcat("$", subsys.StateName(fig_num), "$") ,'Interpreter','latex');   
        xlabel('$t$ [s]','Interpreter','latex')
        hold off
    end
    sgtitle(fig_title);
    
    % figure for input
    fig_title = strcat("Applied Inputs to Subsystem", {' '}, subsystems(sys_num));
    figure('Name', fig_title,'NumberTitle','off');
    hold on; grid on;
    stairs(times(1:end-1), sol.u,'k','markersize',20,'linewidth',2);
    plot(times(1:end-1), umin*o(1:end-1),'c','linewidth',2);
    plot(times(1:end-1), umax*o(1:end-1),'c','linewidth',2);
    xline(settling_time(sys_num), '-.r','linewidth',2); 
    xlabel('$t$ [s]','Interpreter','latex')
    ylim(1.1*[umin umax]);
    ylabel(input_names(fig_num), 'Interpreter','latex');
    legend('input', 'lower bound', 'upper bound', 'settling time', 'Interpreter','latex');
    title(fig_title, 'Interpreter','latex');
end

%% Assistive functions
% find u_min and u_max
function [umin, umax] = get_input_bounds(sys_num)
    switch sys_num
        case 1
            umin = -0.3; umax = 0.3;
        case 2
            umin = -0.3; umax = 0.3;
        case 3
            umin = -0.2; umax = 0.3;
        case 4
            umin = -0.2; umax = 0.2;
        otherwise
            fprintf('[Error] subsystem number out of range.');
    end
end
