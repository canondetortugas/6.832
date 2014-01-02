pendulum_phase_plot(0, 0.5, [pi/4, 0]', 'full'); % Standard system dynamics
pendulum_phase_plot(0, 0.5, [pi/4, 0]', 'pfl_damping'); % Use PFL
                                                        % to remove
                                                        % damping
pendulum_phase_plot(0, 0.5, [pi/4, 0]', 'pfl_2g'); % No damping,
                                                   % double gravity
pendulum_phase_plot(0, 0.5, [pi/4, 0]', 'pfl_ig'); % Inverted
                                                   % gravity with damping