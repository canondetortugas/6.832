function [J,dJdalpha] = cartpolefun(alpha_in)
% dynamics dt
dt = 0.01; T = 1.5;

% pendulum parameters
global mp mc g l xdes;
% mc = 10; mp = 1; l = 0.5; g = 9.8;
mc = 1; mp = 1; l = 0.5; g = 1;
xdes = [0 pi 0 0]'; % the desired final state

N = floor(T/dt)+1;
xtape = zeros(4,N);
utape = zeros(1,N);
alpha = zeros(N,1);
if nargin>0
    alpha = alpha_in;
end

% Simulate forward
IC = [0 0 0 0]';
x = IC; % arbitrary (but fixed) initial condition
for i=1:N
    xtape(:,i) = x;
    u = alpha(i);
    utape(i) = u;
    x = x + dynamics(x,u).*dt;
end

% Plot the trajectory
figure(24)
subplot(2,1,1);
plot(xtape(1,:),xtape(3,:)); xlabel('x'); ylabel('xdot');
subplot(2,1,2);
plot(xtape(2,:),xtape(4,:)); xlabel('theta'); ylabel('thetadot');
drawnow;

dJdalpha = compute_gradients(xtape,utape,dt);
J = sum(cost(xtape,utape,dt)) + finalCost(xtape(:,N));
% sum(cost(xtape,utape,dt))
% finalCost(xtape(:,N))

end % of pendfun


% =========================================================
% This function returns the gradients by
% integrating the adjoint equations
% =========================================================
function dJdalpha = compute_gradients(x,u,dt)
global xdes; %desired x location
N = size(x,2);

[Q,R,Qend] = get_QR;
Q = Q.*dt; R = R.*dt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INTEGRATE ADJOINT EQUATIONS TO PERFORM BPTT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[dfdx,dfdu] = gradients(x(:,N),u(N)); 
dgdu = 2*u(N)*R;                        % Note: Added the 2 term
                                        % compared to original
                                        % hw. Makes sense given
                                        % that g is quadratic in u.
dudalpha = zeros(1,N); dudalpha(N) = 1; %gradient of u
                                        %w.r.t. parameters for open
                                        %loop policy at N. We will
                                        %decrement the index of the
                                        %non-zero column as we
                                        %integrate back in time.

F_alpha = dfdu*dudalpha;
G_alpha = dgdu*dudalpha;

y = [-2*Qend(1,1)*x(1,end);  2*Qend(2,2)*(mod(x(2,end)+pi, 2*pi)-pi); ...
     -2*Qend(3,3)*x(3,end); -2*Qend(4,4)*x(4,end)]; %terminal
                                                 %condition for y
                                                 %(note: this form
                                                 %will be correct
                                                 %if Q is not diagonal
dJdalpha = (G_alpha'-F_alpha'*y).*dt; % dJdalpha for first time step
for n=N-1:-1:1 %integrate adjoint equations backwards in time
    dgdx = zeros(1,4); %gradient of cost with respect to current state 
    dgdu = 2*u(n)*R; %gradient of cost with respect to current action
    [dfdx,dfdu] = gradients(x(:,n),u(n)); %gradient of f w.r.t. current position, action
    F_x = dfdx;
    G_x = dgdx;
    dudalpha = zeros(1,N); dudalpha(n) = 1; %gradient of u w.r.t. parameters for open loop policy at current time
    F_alpha = dfdu*dudalpha;
    G_alpha = dgdu*dudalpha;
    y = y - (G_x - y'*F_x)'.*dt; %solve for y
    dJdalpha = dJdalpha + (G_alpha' - F_alpha'*y).*dt; %add this step's contribution to dJdalpha
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

% ============================================================
% This function evaulates the gradients at a particular x,u
% ============================================================
% NOTE: These expressions were computed with mc = mp =1, l = 0.5
% They won't work if these parameters are changed
function [dfdx,dfdu] = gradients(x,u)
% parameters
global mc mp l g;

% Sin and cosine
s = sin(x(2)); c = cos(x(2));
thetadot = x(4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Give the gradients of the dynamics with respect to 
% a particular state and action
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
df3dthetadot = thetadot*s/(1 + s^2);

df4dthetadot = -2*thetadot*c*s/(1+s^2);

df3dtheta = (c^2 - c*s^2)*((thetadot^2)/2 + g*c)/((1+s^2)^2) - ...
    g*s*(s/(1+s^2));

df4dtheta = 2*( -2*s*c*(-u*c - c*s*(thetadot^2)/2 - 2*g*s)/((1 + s^2)^2) ...
                + (u*s - (c^2 - s^2)*(thetadot^2)/2 - 2*g*c)/(1 + s^2));

dfdx = [zeros(2,2), eye(2,2); 0, df3dtheta, 0, df3dthetadot; 0, ...
        df4dtheta, 0, df4dthetadot ];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gradients of f wrt u
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
df1du = 0;                              % xdot
df2du = 0;                               % thetadot
df3du = 1/(1 + s^2);                    % xdotdot
df4du = -2*c/(1 + s^2);                 % thetadotdot

dfdu = [df1du, df2du, df3du, df4du]';

end

% =============================================================
% This function defines the continuous dynamics of the pendulum
% =============================================================
function xdot = dynamics(x,u)
    global mc mp l g;
    
    s = sin(x(2)); c = cos(x(2));

    %    H = [mc+mp, mp*l*c; mp*l*c, mp*l^2];
    %    C = [0 -mp*x(4)*l*s; 0 0];
    %    G = [0; mp*g*l*s];
    %    B = [1; 0];
    %    xdot = [x(3:4); inv(H)*[B*u - C*x(3:4) - G]];

    xddot = [u + mp*s*(l*x(4)^2 + g*c)]/[mc+mp*s^2];
    tddot = [-u*c - mp*l*x(4)^2*c*s - (mc+mp)*g*s]/[l*(mc+mp*s^2)];
    xdot = [x(3:4); xddot; tddot];
end



% =============================================================
% This function defines the instantaneous cost (i.e. g(x,u))
% =============================================================
function C = cost(X,u,dt)
global xdes;
[Q,R] = get_QR;
Q = Q*dt;
R = R*dt;

Xerr = X - repmat(xdes,1,size(X,2));
Xerr(2,:) = mod(Xerr(2,:)+pi,2*pi)-pi;

% implement a quadratic cost
C = diag( u'*R*u);
end

% =============================================================
% Implements a final cost
% =============================================================

function C = finalCost(X)
global xdes;
[~,~,Qend] = get_QR;

Xerr = X - repmat(xdes,1,size(X,2));
Xerr(2,:) = mod(Xerr(2,:)+pi,2*pi)-pi;  % Wrap pendulum angle error

% implement the final cost cost
C = Xerr(:,end)'*Qend*Xerr(:,end);
end

% ============================================================
% Returns the cost matrices
% ============================================================
function [Q,R, Qend] = get_QR

% penalty matrices
Q = zeros(4,4);
Qend = 500*eye(4,4);
R = 0.01*eye(1,1);

end