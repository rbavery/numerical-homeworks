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