function [x_next]=RK4_quad(x,u,h,f)
%% Compute discrete-time dynamics with Runge-Kutta 4 integrator
k1=f(x,u);
k2=f(x+h/2*k1,u);
k3=f(x+h/2*k2,u);
k4=f(x+h*k3,u);
x_next=x+h/6*(k1+2*k2+2*k3+k4);
end