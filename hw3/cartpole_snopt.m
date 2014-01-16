function w = cartpole_snopt()
% pendulum parameters
    global mc mp g l xdes N;
    % mc = 10; mp = 1; l = 0.5; g = 9.8;

    % NOTE: This function is broken if parameters other than
    % gravity are changed
    mc = 1; mp = 1; l = 0.5; g = 1;
    T = 1.5;
    dt=.01;
    N = floor(T/dt)+1;
    xdes = [0 pi 0 0]'; % the desired final state [x, theta, x, thetadot

    snsummary off;
    snseti     ('Major Iteration limit', 100);

    %don't worry if SNOPT says "Failed to find optimal solution" when it
    %terminates.  So long as the error has dropped below this tolerance,
    %the solution is "optimal enough" for our purposes.
    snsetr     ('Major optimality tolerance',1e-4);

    snprint('cartpole_snopt.log');

    [alpha,alphalow,alphaupp,alphamul,alphastate,Flow,Fupp,Fmul,Fstate,ObjAdd,ObjRow,    ...
     A,iAfun,jAvar,iGfun,jGvar] = penddata;

    [alpha,F,wmul,Fmul,inform]= snsolve( alpha, alphalow, alphaupp, ...
                                         alphamul, alphastate,    ...
                                         Flow, Fupp, Fmul, Fstate,       ...
                                         ObjAdd, ObjRow, A, iAfun, jAvar,...
                                         iGfun, jGvar, 'cartpolefun');

    snset('Defaults');
    

    % Play back the policy. We will use this trajectory to create
    % an LQR stabilizer
    x = [0 0 0 0]';
    xtape = x;
    for i=1:N
        control = alpha(i);
        x = x + dynamics(x, control).*dt;
        xtape = [xtape, x];
    end

    %  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Build LQR stabilizer
    %  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Q = 50*eye(4,4).*dt;
    R = 0.01*eye(1,1).*dt;

    [A, B] = linearize_cartpole_around_sequence(xtape, alpha');
    [K, S] = lqr_stabilizer(A, B, Q, R, dt);

    % keyboard;

    % Simulate the policy with noise
    control_stdev = 16;
    
    x = [0 0 0 0]';
    z = x;
    xtraj = x;
    ztraj = z;
    ctraj = [];
    for i=1:N
        draw(z, (i-1)*dt);
        zerr = z - xtape(:,i);
        control = alpha(i) - K{i}*zerr;
        ctraj = [ctraj, control];
        e = randn(1)*control_stdev;                         % Gaussian RV
        % e = (rand(1)-0.5)*control_stdev;                         % Uniform RV
        x = x + dynamics(x, alpha(i) + e).*dt;
        z = z + dynamics(z, control + e).*dt;
        xtraj = [xtraj, x];
        ztraj = [ztraj, z];
    end

    [J,dJdalpha] = cartpolefun(alpha);
    fprintf('\nCost of the found solution: %3.2f\n',J)
    fprintf('Final stabilized state: [ %f, %f ].\n\n', z(1), z(2));

    figure; hold on;
    plot([1:length(alpha)]*dt, alpha, 'b');
    plot([1:length(alpha)]*dt, ctraj, 'g');
    xlabel('t'); ylabel('\pi_{\alpha{}}(t)');
    title('Policies'); legend('Open Loop', 'Stabilized');
    
    figure(24);
    subplot(2,1,1);
    hold on;
    plot(xtraj(1,:),xtraj(3,:), 'r'); 
    plot(ztraj(1,:),ztraj(3,:), 'g'); 
    xlabel('x'); ylabel('xdot');
    title('Trajectories');
    legend('Target', 'Unstabilized', 'Stabilized');
    subplot(2,1,2);
    hold on;
    plot(xtraj(2,:),xtraj(4,:), 'r'); 
    plot(ztraj(2,:),ztraj(4,:), 'g'); 
    xlabel('theta'); ...
        ylabel('thetadot');
    legend('Target', 'Unstabilized', 'Stabilized');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [alpha,alphalow,alphaupp,alphamul,alphastate,Flow,Fupp,Fmul,Fstate,ObjAdd,ObjRow, ...
	  A,iAfun,jAvar,iGfun,jGvar] = penddata()
    global N;

    ObjRow = 1;
    ObjAdd = 0;

    alpha      = ones(N,1);
    alphalow   = -inf*ones(N,1);
    alphaupp   = inf*ones(N,1);
    alphamul   = zeros(N,1);
    alphastate = zeros(N,1);

    Flow   = -inf;
    Fupp   = inf;
    Fmul   = 0;
    Fstate = 0;

    A     = [];
    iAfun = [];
    jAvar = [];

    iGfun = ones(N,1); jGvar = [1:N]';
end

% =============================================================
% This function defines the continuous dynamics of the cart pole
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


% ===============================================================
% This is the draw function
%================================================================
function draw(x,t)
    global mc mp l g;

    persistent hFig base a1 raarm wb lwheel rwheel;
    if (isempty(hFig))
        hFig = figure(25);
        set(hFig,'DoubleBuffer', 'on');

        a1 = l+0.25;
        av = pi*[0:.05:1];
        theta = pi*[0:0.05:2];
        wb = .3; hb=.15;
        aw = .01;
        wheelr = 0.05;
        lwheel = [-wb/2 + wheelr*cos(theta); -hb-wheelr + wheelr*sin(theta)]';
        base = [wb*[1 -1 -1 1]; hb*[1 1 -1 -1]]';
        arm = [aw*cos(av-pi/2) -a1+aw*cos(av+pi/2)
               aw*sin(av-pi/2) aw*sin(av+pi/2)]';
        raarm = [(arm(:,1).^2+arm(:,2).^2).^.5, atan2(arm(:,2),arm(:,1))];
    end

    figure(hFig); cla; hold on; view(0,90);
    patch(x(1)+base(:,1), base(:,2),0*base(:,1),'b','FaceColor',[.3 .6 .4])
    patch(x(1)+lwheel(:,1), lwheel(:,2), 0*lwheel(:,1),'k');
    patch(x(1)+wb+lwheel(:,1), lwheel(:,2), 0*lwheel(:,1),'k');
    patch(x(1)+raarm(:,1).*sin(raarm(:,2)+x(2)-pi),-raarm(:,1).*cos(raarm(:,2)+x(2)-pi), 1+0*raarm(:,1),'r','FaceColor',[.9 .1 0])
    plot3(x(1)+l*sin(x(2)), -l*cos(x(2)),1, 'ko',...
          'MarkerSize',10,'MarkerFaceColor','b')
    plot3(x(1),0,1.5,'k.')
    title(['t = ', num2str(t,'%.2f') ' sec']);
    set(gca,'XTick',[],'YTick',[])

    axis image; axis([-2.5 2.5 -1.5*l 1.5*l]);
    drawnow;
end

% NOTE: These expressions were computed with mc = mp =1, l = 0.5
% They won't work if these parameters are changed
function [A, B] = linearize_cartpole_around_sequence(xtraj, utraj)
% parameters
    global mc mp l g;


    % Discard the final x point (where we end up after the final
    % control input is applied)
    for idx = 1:(length(xtraj)-1)
        
        x = xtraj(:,idx);
        u = utraj(:,idx);
        
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

        A{idx} = dfdx;


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Gradients of f wrt u
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        df1du = 0;                              % xdot
        df2du = 0;                               % thetadot
        df3du = 1/(1 + s^2);                    % xdotdot
        df4du = -2*c/(1 + s^2);                 % thetadotdot

        dfdu = [df1du, df2du, df3du, df4du]';

        B{idx} = dfdu;
    end

end