% functions
fi = @(x)exp(x);
fii = @(x)x.^2+2
fiii = @(x)sin(x)
fiv = @(x)log(x+1)/(x+1)
% evaluate numerically
dx = .01;

%fi
integ_plot(dx, fi, 5, 1)
subplot(3,1,3)
xb = [0:dx:5];
fia = @(x)exp(x)
plot(xb(2:end),fia(xb(2:end)))

integ_plot(dx, fii, 5, 2)

integ_plot(dx, fiii, 5, 3)

integ_plot(dx, fiv, 5, 4)


function answer = integ_plot(dx, func, b, figi)
    
    xb = [0:dx:b];
    
    nxb = length(xb);
    
    xmidpoints = @(x, nx).5 * (x(1:nx-1)+x(2:nx));
    
    xbmids = xmidpoints(xb, nxb);
    
    ybmids = func(xbmids);
    
    Fix = cumsum(ybmids(1:nxb-1).*dx);

    figure(figi);
    subplot(3, 1, 1);
    plot(xb(2:end), Fix);
    title('F(x) vs x')
    
    subplot(3, 1, 2);
    plot(xb(2:end), ybmids);
    title('f(x) vs x')
    
    end

    

function answer = integrate(dx, func, a, b)
    
    xa = [0:dx:a];
    xb = [0:dx:b];
    nxa = length(xa);
    nxb = length(xb);
    xmidpoints = @(x, nx).5 * (x(1:nx-1)+x(2:nx));
    xamids = xmidpoints(xa, nxa);
    xbmids = xmidpoints(xb, nxb);
    
    yamids = func(xamids);
    ybmids = func(xbmids);
    
    y_int_0_b = sum(ybmids(1:nxb-1).*dx);
    y_int_0_a = sum(yamids(1:nxa-1).*dx);
    y_int_a_b = y_int_0_b - y_int_0_a;
    
    answer = y_int_a_b
end


    