dx = .1;
dy = .1;

x = [0:dx:10];
y = [0:dy:10];
nx = length(x);
ny = length(y);

[X,Y] = meshgrid(x,y)

Z = sin(X/1.4)+cos(Y/1.4)+sqrt(exp(X/4));
figure(1);
contourf(X,Y,Z);
colorbar

figure(2);
surf(X,Y,Z);
colorbar

%crit points
min1 = [6, 4.5]; %minima
min2 = [0, 4.5]; %minima
saddle1 = [6, 0];
saddle2 = [6, 10];

%4
Xmx = .5*(X(1:nx-1,:)+X(2:nx,:)); % xmid = .5*(x(2:nx)+x(1:nx-1))
Ymx = Y(:,2:nx); % taking last 2:nx elements of Y in x-direction

Xmy = .5*(X(2:ny, :)); % taking last 2:nx elements of Y in x-direction
Ymy = Y(:,1:ny-1)+Y(:,2:ny);

dX = diff(X, 1, 2);
dY = diff(Y, 1, 1);

dZx = diff(Z,1,2);
dZy = diff(Z,1,1);

partZx = dZx./dX
partZy = dZy./dY

figure(3);
contourf(Xmx,Ymx,partZx);
colorbar
title('Partial Derivitave wrt x');
figure(4);
surf(Xmx,Ymx,partZx);
colorbar
title('Partial Derivitave wrt x');
figure(5);
contourf(Xmy,Ymy,partZy);
colorbar
title('Partial Derivitave wrt y');
figure(6);
surf(Xmy,Ymy,partZy);
colorbar
title('Partial Derivitave wrt y');

partZxmids = .5.*(partZx(2:end,:)+partZx(1:end-1,:))
partZymids = .5.*(partZy(:,2:end)+partZy(:,1:end-1))

steepness = sqrt((partZxmids.^2)+(partZymids.^2));

figure(6);
contourf(Xmx,Ymy,steepness);
colorbar

figure(7);
surf(Xmx,Ymy,steepness);
colorbar










