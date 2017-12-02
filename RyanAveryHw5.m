load('socal_elevation_v2.mat');

[LON, LAT] = meshgrid(lon, lat);

% 1.
figure(1)
surf(LON, LAT, hgt)
title('3D Elevation Plot')
xlabel('Longitude')
ylabel('Latitude')
zlabel('Elevation')

dZlon = diff(hgt,1,2);
dZlat = diff(hgt,1,1);

dLONlon = diff(LON, 1, 2);
dLATlat = diff(LAT, 1, 1);

pZlon = dZlon./dLONlon;
pZlat = dZlat./dLATlat;

LONmid = .5*(LON(:,2:end) + LON(:,1:end-1));
LATtruncbyLON = LAT(:, 2:end);

LATmid = .5*(LAT(2:end,:) + LAT(1:end-1,:));
LONtruncbyLAT = LON(2:end, :);

figure(2)
subplot(2, 1, 1)
contour(LONmid, LATtruncbyLON, pZlon)
title('Partial LON Derivitive of Height Plot')
xlabel('Longitude')
ylabel('Latitude')
zlabel('Steepness by LON')
subplot(2, 1, 2)
contour(LONtruncbyLAT, LATmid, pZlat)
title('Partial LAT Derivitive of Height Plot')
xlabel('Longitude')
ylabel('Latitude')
zlabel('Steepness by LAT')

% generally steeper in west east direction (santa ynez mountains)
figure(3)
lonmid = .5.*(lon(2:end)+lon(1:end-1));
latmid = .5.*(lat(2:end)+lat(1:end-1));
[LONmid, LATmid] = meshgrid(lonmid,latmid);
pZlatmids = .5.*(pZlat(:,1:end-1)+pZlat(:,2:end));
pZlonmids = .5.*(pZlon(1:end-1,:)+pZlon(2:end,:));
aplusbsquared = pZlatmids.^2+pZlonmids.^2;
steepness = sqrt(aplusbsquared);
contour(LONmid, LATmid, steepness)
title('Steepness Plot')
xlabel('Longitude')
ylabel('Latitude')
zlabel('Steepness')

% 4 
p2Zlon = diff(pZlon,2,2);
p2Zlat = diff(pZlat,2,1);

p2Zlon = (p2Zlon./dLONlon(:, 2:end-1));
p2Zlat = (p2Zlat./dLATlat(2:end-1,:));

p2Zlatmids = .5.*(p2Zlat(:,1:end-1)+p2Zlat(:,2:end));
p2Zlonmids = .5.*(p2Zlon(1:end-1,:)+p2Zlon(2:end,:));

LONtrunc2 = LON(2:end,2:end);
LATtruncbyLON2 = LATtruncbyLON(2:end,:);

LAT2 = LAT(4:end, 2:end)
LONtruncbyLAT2 = LONtruncbyLAT(3:end,2:end)

figure(4)

subplot(2, 1, 1)
contourf(LONtrunc2, LATtruncbyLON2, p2Zlonmids)
title('Roughness LON Direction Plot')
xlabel('Longitude')
ylabel('Latitude')
zlabel('Roughness')
subplot(2, 1, 2)
contourf(LONtruncbyLAT2, LAT2, p2Zlatmids)
title('Roughness LAT Direction Plot')
xlabel('Longitude')
ylabel('Latitude')
zlabel('Roughness')


















