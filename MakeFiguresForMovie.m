%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% MakeFiguresForMovie
%
% Reads output from EddyFieldWithNIO.m and make cool figures for a movie
%
% Jeffrey J. Early
% jeffrey@jeffreyearly.com


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Open the model output, create a new WVT
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = "dmd-eddy-tide.nc";
[wvt, ncfile] = WVTransform.waveVortexTransformFromFile(filename,iTime=2);
t_full = ncfile.readVariables('t');

%%

% iTime=20;
% wvt.initFromNetCDFFile(ncfile,iTime=iTime)

xIndices = 1:wvt.Nx;
yIndices = floor(wvt.Ny/2);
zIndices = 1:wvt.Nz;
horzAxis = wvt.y/1e3;
vertAxis = wvt.z;

eta = wvt.eta_w;
rv = wvt.diffX(wvt.v_g) - wvt.diffY(wvt.u_g);
zeta_limits = [-0.3 0.3];

figure
tiledlayout('flow')
colormap('gray')

nexttile
p3 = pcolor(wvt.x/1e3, wvt.y/1e3, squeeze(wvt.zeta_z(:,:,end)).'/wvt.f); shading interp, axis equal
xlim([min(wvt.x) max(wvt.x)]/1e3), ylim([min(wvt.y) max(wvt.y)]/1e3)
cb2 = colorbar('eastoutside');
cb2.Label.String = 'f_0';
clim(zeta_limits)
set(gca, 'XTickLabel', [])


nexttile
p2 = pcolor(horzAxis,vertAxis,squeeze(wvt.zeta_z(xIndices,yIndices,zIndices)).'/wvt.f); shading interp, hold on
% title('rv')
cb2 = colorbar('eastoutside');
cb2.Label.String = 'f_0';
clim(zeta_limits)



return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Relative vorticity figures
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xIndices = floor(wvt.Nx/2);
yIndices = 1:wvt.Ny;
zIndices = 1:wvt.Nz;
horzAxis = wvt.y/1e3;
vertAxis = wvt.z;

fig1 = figure('Units', 'points', 'Position', [50 50 860 400]);
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');

t = ncfile.readVariables('t');
t = t-t(1);
zeta_range = 0.15;
for iTime = 1:1 %length(t)
    clf
    wvt.initFromNetCDFFile(ncfile,iTime=iTime)

    eta = wvt.eta_w;
    rv = wvt.diffX(wvt.v_g) - wvt.diffY(wvt.u_g);

    p3 = pcolor(horzAxis,vertAxis,squeeze(100*eta(xIndices,yIndices,zIndices)).'); shading interp, hold on
    colormap(crameri('cork'))
    contour(horzAxis,vertAxis,squeeze(rv(xIndices,yIndices,zIndices)).'/wvt.f,[-0.01 -0.01],'k',linewidth=1.5)
    contour(horzAxis,vertAxis,squeeze(rv(xIndices,yIndices,zIndices)).'/wvt.f,[0.01 0.01],'k',linewidth=1.5)
    clim([-1 1]*600)
    title('\eta (isopycnal deviation)')
    cb3 = colorbar('eastoutside');
    cb3.Label.String = 'cm';

    % print(sprintf('movie-figures/t-%03d.png',iTime),'-r300','-dpng')
end