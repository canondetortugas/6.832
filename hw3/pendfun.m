function [J,dJdalpha] = pendfun(alpha_in)
% dynamics dt
dt = 0.01; T = 2.5;

% pendulum parameters
global m g l I b xdes;
m=1; g = 9.8; l = 1; I = m*l*l; b = 0.1;
xdes = [pi 0]'; % the desired final state

N = floor(T/dt)+1;
xtape = zeros(2,N);
utape = zeros(1,N);
alpha = zeros(N,1);
if nargin>0
    alpha = alpha_in;
end

% Simulate forward
IC = [0 0]';
x = IC; % arbitrary (but fixed) initial condition
for i=1:N
    xtape(:,i) = x;
    u = alpha(i);
    utape(i) = u;
    x = x + dynamics(x,u).*dt;
end

figure(24); hold off;
plot(xtape(1,:),xtape(2,:)); drawnow;

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

y = [ 2*Qend(1,1)*(mod(x(1,end)+pi, 2*pi)-pi); -2*Qend(2,2)*x(2,end)]; %terminal
                                                 %condition for y
                                                 %(note: this form
                                                 %will be correct
                                                 %if Q is not diagonal
dJdalpha = (G_alpha'-F_alpha'*y).*dt; % dJdalpha for first time step
for n=N-1:-1:1 %integrate adjoint equations backwards in time
    dgdx = zeros(1,2); %gradient of cost with respect to current state 
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
function [dfdx,dfdu] = gradients(x,u)
% pend parameters
global m g l I b;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Give the gradients of the dynamics with respect to 
% a particular state and action
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dfdx = [0, 1; -(g/l)*cos(x(1)), -b/(m*l^2)];
dfdu = [0; 1/(m*l^2)];

end

% =============================================================
% This function defines the continuous dynamics of the pendulum
% =============================================================
function xdot = dynamics(x,u)
% pendulum parameters
global m g l I b;

xdot = [x(2,:); (u-m*g*l*sin(x(1,:))-b*x(2,:))./I];
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
Xerr(1,:) = mod(Xerr(1,:)+pi,2*pi)-pi;

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
Xerr(1,:) = mod(Xerr(1,:)+pi,2*pi)-pi;

% implement the final cost cost
C = Xerr(:,end)'*Qend*Xerr(:,end);
end

% ============================================================
% Returns the cost matrices
% ============================================================
function [Q,R, Qend] = get_QR

% penalty matrices
Q = zeros(2,2);
Qend = 500*eye(2,2);
R = 0.01*eye(1,1);

end