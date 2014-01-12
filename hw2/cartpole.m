function cartpole

% Parameters:
% mc = 10; mp = 1; l = 0.5; g = 9.8;
    mc = 1; mp = 1; l = 1; g = 1;
    T = 30;
    plant_dt = 1e-3;
    display_dt = 0.1;

    % Linearized dynamics
    dGdq = [0 0; 0 -mp*g*l];
    H = [mc + mp, -mp*l; -mp*l, mp*l^2];
    A = [zeros(2,2), eye(2,2); -inv(H)*dGdq, zeros(2,2)];
    B = [ zeros(2,4); zeros(2,2), inv(H)*[1,0;0,0] ];
    % dGdq = [0 0; 0 mp*g*l];
    % H = [mc + mp, mp*l; mp*l, mp*l^2];
    % A = [zeros(2,2), eye(2,2); -inv(H)*dGdq, zeros(2,2)]
    % B = [ zeros(2,4); zeros(2,2), inv(H)*[1,0;0,0] ]


    % LQR
    Q = diag([100, 100, 100, 100])*plant_dt;
    R = 10*eye(4,4)*plant_dt; % A little silly since only one dimension gets
                              % passed into the system, but the LQR solution
                              % takes care of this

    % TODO: Figure out how discretization in these functions works
    sys = ss(A,B,eye(4,4), zeros(4,4));
    sysd = c2d(sys, plant_dt);
    [K,S] = dlqr(sysd.a,sysd.b,Q,R);
    
    % [K,S] = lqr(A,B,Q,R)


    % Initial Conditions:  
    target_state = [0, pi, 0, 0]'; % [ x, theta, xdot, thetadot]
    target_energy = mp*g*l; % Potential energy at pi with no kinetic energy
                            % target_state = [0,0,0,0]';

    randn('state',sum(100*clock))
    %   x = [x,\theta,\dot{x},\dot\theta]^T 
    %   (set this to [0;0;0;0] to test true swing-up)
    % x = target_state +  0.2*randn(4,1);  
    x = 0.5*randn(4,1); 
    % x = [0, 0.01*randn(1,1), 0, 0.01*randn(1,1)]'; 
    % x = zeros(4,1);


    % Euler Integration Loop:
    last_display_t = -inf;
    for t=0:plant_dt:T
        u = control(x,t);

        if (t>last_display_t + display_dt)
            draw(x,t);
            last_display_t = t;

            % Energy calc (good way to verify eqs. of motion)
            T = 0.5*(mc+mp)*x(3)^2 + mp*x(3)*x(4)*l*cos(x(2)) + 0.5*mp*l^2*x(4)^2;
            U = -mp*g*l*cos(x(2));
            E = T+U
        end
        
        xdot = dynamics(x,u);
        x = x + plant_dt*xdot;
    end
    draw(x,t);


    function u = control(x,t)
    % u = lqr_control(x,t);
    % u = energy_control(x,t);
    % u = -x(1);
    % u = 0;p
    % TODO: Fix this so that it properly triggers when angle is
    % near pi
        r1 = abs(mod(x(2), 2*pi) - target_state(2)) < 1;
        r2 = abs(x(4) - target_state(4)) < 1;

        [r1 r2]

        if r1 && r2
            u = lqr_control(x,t);
            disp('LQR On');
        else
            u = energy_control(x,t);
        end
    end

    function u = lqr_control(x,t)
        z = x - target_state;
        u = -K*z;
        u = u(3); % This should be the only non-zero component
                  % since our R matrix has no coupling 
    end

    function u = energy_control(x,t)
        theta = x(2);
        thetadot = x(4);

        % Gain
        k = 1;
        
        % Energy shaping control
        % TODO: Fix this to depend on coefficients
        E = 0.5*thetadot^2 - cos(theta);
        Eerr = E - target_energy;
        xddotd = k*thetadot*cos(theta)*Eerr;

        % xdottd = 0;
        
        % Collocated partial feedback linearization - forces
        % x double dot equal to our desired input force
        % TODO: Fix this so that it doesn't only work when all
        % coefficients are 1
        u = (mp+mc - mp*(cos(theta))^2)*xddotd ...
            - mp*l*(thetadot^2)*sin(theta) ...
            - mp*g*sin(theta)*cos(theta);
    end
    
    function xdot = dynamics(x,u)
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


    function draw(x,t)
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

        axis image; axis([-10 10 -1.5*l 1.5*l]);
        drawnow;
    end


end


