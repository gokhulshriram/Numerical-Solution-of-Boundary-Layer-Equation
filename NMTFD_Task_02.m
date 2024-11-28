%%
close all
clear
clc
%% Defining dx and dy
dx = 0.0001;
dy = 0.001;
x = 0:dx:1;
y = 0:dy:0.1;
no_x = size(x,2);
no_y = size(y,2);
%%
Re = 10000;
v_extras = (dy/(2*dx));
u_extra_1 = (dx/(Re*(dy^2)));
u_extra_2 = (dx/(2*dy));
%% Setting Boundary conditions
u = zeros(no_y,no_x);
v = zeros(no_y,no_x);

u(1,:) = 0;     %boundary condition; not needed
v(1,:) = 0;     %boundary condition; not needed
u(end,:) = 1;   %boundary condition
u(:,1) = 1;
u(1,1) = 0;
[X,Y] = meshgrid(x,y);
%%
for j = 1:(no_x - 1) % traverse theough rows
    for i = 2:(no_y - 1) % traverse through columns
        % Update for u(i,j+1)
        u(i, j + 1) = u(i, j) + u_extra_1 * (u(i + 1, j) - 2 * u(i, j) + u(i - 1, j)) / u(i, j) - v(i, j) * u_extra_2 * (u(i + 1, j) - u(i - 1, j)) / (u(i, j));
        % Update for v(i,j+1)
        v(i, j + 1) = v(i - 1, j + 1) - v_extras * (u(i, j + 1) - u(i, j) + u(i - 1, j + 1) - u(i - 1, j));
    end
    % Update v for i = no_y
    v(no_y, j + 1) = v(no_y - 1, j + 1) - 0.5 * dy / dx * (u(no_y, j + 1) - u(no_y, j) + u(no_y - 1, j + 1) - u(no_y - 1, j));
end
%%

contourf(x,y,v,'LineColor', 'none');
colormap('parula');
colorbar;
xlabel('x')
ylabel('y')
title('v velocity contour')