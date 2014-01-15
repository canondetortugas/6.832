function w = cartpole_snopt()
% pendulum parameters
    global mc mp g l xdes N;
    % mc = 10; mp = 1; l = 0.5; g = 9.8;
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

    % playback the learned policy
    x = [0 0 0 0]';
    for i=1:N
        draw(x, (i-1)*dt);
        x = x + dynamics(x,alpha(i)).*dt;
    end

    [J,dJdalpha] = cartpolefun(alpha);
    fprintf('\nCost of the found solution: %3.2f\n\n',J)
    fprintf('Final state: [ %f, %f ].\n', x(1), x(2));

    figure; hold on;
    plot([1:length(alpha)]*dt, alpha);
    xlabel('t'); ylabel('\pi_{\alpha{}}(t)');
    title('Optimal Policy');

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