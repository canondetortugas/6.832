% Simulates a simple pendulum and calculates the basin of
% attraction for a given fixed point.
function calc_basin(u, b, fp)
% pendulum parameters
m=1; g = 9.8;
l = 1; I = m*l*l;
dt = 1e-2;
t0 = 0; tf = .1;
nbins = 50;                              % 

th_bins = linspace(-2*pi,2*pi,nbins);
thdot_bins = linspace(-2*g/l,2*g/l,nbins);
[th,thdot] = ndgrid(th_bins,thdot_bins);
basin = zeros(size(th));
max_iterations = 1e2;

for i=1:size(th,1)
    for j=1:size(th,2)
       x0 =  [th(i,j) thdot(i,j)]';
       
       [t, x] = ode45(@dynamics, [0, max_iterations], x0);
       
       xf = x(end,:)';
       if( norm(xf-fp)< 1 )
           basin(i,j) = true;
       end
       
    end
    if (mod(i,5)==0)
        disp([num2str((i/size(th,1))*1e2),' % done']);
    end
end

figure;
imagesc([th(1,1) th(end,1)],[thdot(1,1) thdot(1,end)],basin');
axis xy; colormap('winter');colorbar;
xlabel('theta'); ylabel('theta dot');
plot_title = sprintf('Basin of attraction for [%s], u=%f, b=%f', ...
                     num2str(fp'), u, b);
title(plot_title);

    function xdot = dynamics(t, x)
        xdot = [x(2); (u-m*g*l*sin(x(1))-b*x(2))./I];
    end
end

