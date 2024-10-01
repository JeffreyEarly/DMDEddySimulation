%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% MakeMovieFigures
%
% Jeffrey J. Early
% jeffrey@jeffreyearly.com

filename = "dmd-eddy-tide.nc";
[wvt, ncfile] = WVTransform.waveVortexTransformFromFile(filename,iTime=1);

outputVar = WVVariableAnnotation('ssh_w',{'x','y'},'m', 'sea-surface height, wave component');
wvt.addOperation(WVOperation('ssh_w', outputVar,@(wvt) wvt.p_w(:,:,end)/(wvt.rho0*wvt.g)));

outputVar = WVVariableAnnotation('ssh_pv',{'x','y'},'m', 'sea-surface height, pv component');
wvt.addOperation(WVOperation('ssh_pv', outputVar,@(wvt) (wvt.p_mda(:,:,end)+wvt.p_g(:,:,end))/(wvt.rho0*wvt.g)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Relative vorticity figures
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig1 = figure('Units', 'points', 'Position', [50 50 860 600]);
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');

t = ncfile.readVariables('t');
t = t-t(1);
zeta_range = 0.15;
for iTime = 1:length(t)
    clf
    wvt.initFromNetCDFFile(ncfile,iTime=iTime)
    depth = 0;
    depthIndex = find(wvt.z <= 0,1,'last');
    rv_g = wvt.diffX(wvt.v_g) - wvt.diffY(wvt.u_g);
    rv_igw = wvt.diffX(wvt.v_w) - wvt.diffY(wvt.u_w);

    tl = tiledlayout(2,3,'TileSpacing','Compact');
    title(tl,sprintf('DMD simulation, day %d, %d:%02d',floor(t(iTime)/86400),floor(mod(t(iTime)/3600,24)),floor(mod(t(iTime)/60,60))));

    Ejk_g = wvt.transformToRadialWavenumber(wvt.A0_TE_factor.* (wvt.A0.*conj(wvt.A0)));
    Ejk_w = wvt.transformToRadialWavenumber(wvt.Apm_TE_factor .* (abs(wvt.Ap).^2 + abs(wvt.Am).^2));

    %transformFromDFTGridToWVGrid
    ssh_g_bar = wvt.horizontalModes.transformFromDFTGridToWVGrid(wvt.horizontalModes.transformFromSpatialDomain(wvt.ssh_pv));
    ssh_w_bar = wvt.horizontalModes.transformFromDFTGridToWVGrid(wvt.horizontalModes.transformFromSpatialDomain(wvt.ssh_w));
    Sssh_g = wvt.transformToRadialWavenumber(abs(ssh_g_bar).^2);
    Sssh_w = wvt.transformToRadialWavenumber(abs(ssh_w_bar).^2);

    kAxis = wvt.kRadial; % radians/m
    lambda = kAxis*1e3/(2*pi); % cycles/km
    dk = kAxis(2) - kAxis(1);
    SWOT_error = 2 + 0.00125*(lambda).^(-2); % 'cm^2/(cycle/km)'
    Sconv = 1e4*2*pi/1e3; % convert from power in m^2/(radians/m) to cm^2/(cycles/km)

    sp1 = nexttile;
    pcolor(wvt.x/1000,wvt.y/1000,rv_g(:,:,depthIndex).'/wvt.f), shading interp
    axis equal, xlim([min(wvt.x) max(wvt.x)]/1000), ylim([min(wvt.y) max(wvt.y)]/1000)
    % colormap("gray")
    colormap(sp1,cmocean('balance'))
    ylabel('km'), title('surface vorticity (geostrophic)')
    cb1 = colorbar('southoutside');
    % cb1.Label.String = 'vorticity (f)';
    set(gca,'XTickLabel',[]);
    clim(zeta_range*[-1 1])

    sp2 = nexttile;
    pcolor(wvt.x/1000,wvt.y/1000,rv_igw(:,:,depthIndex).'/wvt.f), shading interp
    axis equal, xlim([min(wvt.x) max(wvt.x)]/1000), ylim([min(wvt.y) max(wvt.y)]/1000)
    % colormap("gray")
    colormap(sp2,cmocean('balance'))
    % ylabel('km')
    title('surface vorticity (igw)')
    cb2 = colorbar('southoutside');
    % cb1.Label.String = 'vorticity (f)';
    set(gca,'XTickLabel',[]);
    set(gca,'YTickLabel',[]);
    clim(zeta_range*[-1 1])


    sp3 = nexttile;
    plot(kAxis,2*pi*1e-3*squeeze(sum(Ejk_g,1))/dk,LineWidth=2), hold on
    plot(kAxis,2*pi*1e-3*squeeze(sum(Ejk_w,1))/dk,LineWidth=2)
    xlog, ylog
    ylabel('m^3/s^2/(cycles/km)')
    set(gca,'XTickLabel',[]);
    % xlabel('1/m')
    set(gca,'YAxisLocation','right')
    ylim([1e0 2e5])
    xlim([kAxis(2) max(kAxis)/3])
    title('depth-integrated total energy spectrum')
    legend('geostrophic', 'igw')

    sp4 = nexttile;
    pcolor(wvt.y/1000,wvt.z/1000,squeeze(rv_g(1,:,:)).'/wvt.f), shading interp
    colormap(sp4,cmocean('balance'))
    clim(zeta_range*[-1 1])
    xlabel('km'), ylabel('km'), title('interior vorticity (geostrophic)')

    sp5 = nexttile;
    pcolor(wvt.y/1000,wvt.z/1000,squeeze(rv_igw(1,:,:)).'/wvt.f), shading interp
    colormap(sp5,cmocean('balance'))
    clim(zeta_range*[-1 1])
    xlabel('km'),title('interior vorticity (igw)')
    set(gca,'YTickLabel',[]);

    sp6 = nexttile;
    plot(kAxis,Sconv*squeeze(Sssh_g)/dk,LineWidth=2), hold on
    plot(kAxis,Sconv*squeeze(Sssh_w)/dk,LineWidth=2)
    plot(kAxis,SWOT_error,"Color","k",LineWidth=2)
    xlog, ylog
    ylabel('cm^2/(cycle/km)')
    set(gca,'XTickLabel',[]);
    set(gca,'YAxisLocation','right')
    xlim([kAxis(2) max(kAxis)/3])
    ylim([1e-1 1e4])
    title('ssh spectrum')
    legend('geostrophic', 'igw','swot error')
    xticks(2*pi./(1e3*[1500 150 15]))
    labels = cell(3,1); labels{1} = '1500 km';  labels{2} = '150 km';  labels{3} = '15 km';
    xticklabels(labels)

    % sp3 = nexttile;
    % plot(wvt.kRadial,sum(Ekj_g,2)/(wvt.kRadial(2)-wvt.kRadial(1))), hold on
    % plot(wvt.kRadial,sum(Ekj_w,2)/(wvt.kRadial(2)-wvt.kRadial(1)))
    % xlog, ylog
    % ylabel('m^3/s^2')
    % xlabel('1/m')

    

    print(sprintf('movie-figures/t-%03d.png',iTime),'-r300','-dpng')
end
% 
