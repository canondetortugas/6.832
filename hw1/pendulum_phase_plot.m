function pendulum_phase_plot(u, b, x0, type)
    m = 1;                              % mass of ball
    l =1;                               % length
    g = 9.8;                            % gravity
    
    iterations = 1e3;

    if strcmpi(type, 'full')
        [t, x] = ode45(@full_dynamics, [0, iterations], x0);
    elseif strcmpi(type, 'pfl_damping')
        [t, x] = ode45(@simple_pfl, [0, iterations], x0);
    elseif strcmpi(type, 'pfl_2g')
        [t, x] = ode45(@pfl_2g, [0, iterations], x0);
    elseif strcmpi(type, 'pfl_ig')
        [t, x] = ode45(@pfl_ig, [0, iterations], x0);
    end

    figure;
    scatter( x(:,1), x(:,2), 5, jet(length(x)));
    xlabel('theta');ylabel('theta dot');
    title(type);

    % Invert gravity
    function xdot = pfl_ig(t,x)
        xdot = dynamics(u + 2*m*g*l*sin(x(1)), x);
    end
    
    % Double gravity
    function xdot = pfl_2g(t, x)
        xdot = pfl_dynamics(u- m*g*l*sin(x(1)),x);
    end
    
    function xdot = simple_pfl(t, x)
        xdot = pfl_dynamics(u,x);
    end
    
    % Use feedback linearization to eliminate damping term
    function xdot = pfl_dynamics(u, x)
        xdot = dynamics(u + b*x(2), x);
    end

    function xdot = full_dynamics(t, x)
        xdot = dynamics(u, x);
    end

    % Dynamics for a motorized pendulum with damping, mass at tip,
    % and negligible pole mass
    function xdot = dynamics(u,x)
        xdot = [x(2); (u - m*g*l*sin(x(1)) - b*x(2))/(m*l*l)];
    end

end
