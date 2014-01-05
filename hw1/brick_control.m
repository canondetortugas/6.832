function [tout,xout]=brick_control()
global A B K

IC=[2; 1]; %initial conditions

% Set up the system: xdot = Ax + Bu
A = [0, 1; 0, 0];
B = [0, 0; 0, 1];

% LQR Parameters
Q = 0.25*eye(2); % Yields a policy similar to bang-bang with low
                 % torque limit
%Q = 1000*[1, 0;0, 0]; % This yields a policy similar to bang-bang
                       % with high torque limit
R = 10*eye(2);                          % Control cost

% Solve algebraic ricatti equation for our system and cost
% matrices, yielding steady-state / infinite-horizon LQR cost-to-go
% matrix (S) and the control gain (K)
[K, S] = lqr(A, B, Q, R);

% simulate system
dt=.01;
T=20;
[tout,xout]=ode45(@control_dyn,[0:.01:T],IC);

% Calculate time taken 
time = 0;
for idx = 1:length(tout)
    x = xout(idx, :);
    if abs(x(1)) < 0.05 && abs(x(2)) < 0.05
        time = tout(idx);
        break;
    end
end

disp(sprintf('Time to reach (0, 0): [ %f ].', time));

% Phase plot
plot_phase(xout);

% Animate brick with the trajectory we solved for
plotdt=.1;
for i=1:plotdt/dt:T/dt
    draw((i-1)*dt,xout(i,:));
    pause(.001)
end

end

%dynamics function
function xdot=control_dyn(t,x)
global A B

u=lqr_control(x);
%u=min_time_control(x); %this controller may take a bit longer to simulate. that's fine

xdot=A*x+B*u;
end

%implement the LQR controller in this function
function u=lqr_control(x)
global K
u = -K*x;
end

% implement the minimum time controller in this function
% Derived with Pontryagin
% u is constrained to lie in [-1, 1]
function u=min_time_control(x)
    
    if x(2) > min_time_control_surface(x(1))
        c = -1;
    else
        c = 1;
    end
    u = [0;c];
end

function qdot = min_time_control_surface(q)
    qdot = -sign(q).*sqrt(2*sign(q).*q);
end

% ==============================================================
% This is the draw function.
% ==============================================================
function draw(t,x)

persistent hFig blockx blocky;

if (isempty(hFig))
  hFig = figure(25);
  set(hFig,'DoubleBuffer','on');
  blockx = [-1, -1, 1, 1, -1];
  blocky = [0, 0.5, 0.5, 0, 0];
end

figure(hFig);
clf;

% draw the mass
brickcolor=[.75 .6 .5];
fill(blockx+repmat(x(1),1,5),blocky,brickcolor);
hold on

faintline=[.6 .8 .65]*1.1;
plot(min(blockx)+[0 0],[-5 5],'k:','Color',faintline);
plot(max(blockx)+[0 0],[-5 5],'k:','Color',faintline);

% draw the ground
line([-5, 5], [0, 0],'Color',[.3 .5 1],'LineWidth',1);
axis([-5 5 -1 2]);
%grid on
axis equal;
title(['t = ', num2str(t)]);

drawnow;
end

function plot_phase(x)
    figure; hold on;
    
    s = -5:.1:5;
    sdot = min_time_control_surface(s);
    plot(s, sdot, 'k');
        
    scatter(x(:,1), x(:,2), 1, jet(length(x)));
    title('Phase Plot');
    xlabel('q'); ylabel('qdot');
    legend('Minimum Time Switching Surface', 'Trajectory');
end