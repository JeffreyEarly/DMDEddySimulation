%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% MakeMovieFigures
%
% Jeffrey J. Early
% jeffrey@jeffreyearly.com

filename = "dmd-eddy-tide.nc";
[wvt, ncfile] = WVTransform.waveVortexTransformFromFile(filename,iTime=1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Extract mooring time series
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = ncfile.readVariables('t');
t = t-t(1);
cv = zeros(length(t),wvt.Nz);
cv_g = zeros(length(t),wvt.Nz);
cv_w = zeros(length(t),wvt.Nz);
for iTime = 1:length(t)
    clf
    wvt.initFromNetCDFFile(ncfile,iTime=iTime)
    
    cv(iTime,:) = wvt.u(1,1,:) + sqrt(-1)*wvt.v(1,1,:);
    cv_g(iTime,:) = wvt.u_g(1,1,:) + sqrt(-1)*wvt.v_g(1,1,:);
    cv_w(iTime,:) = wvt.u_w(1,1,:) + sqrt(-1)*wvt.v_w(1,1,:);
end

max(abs(wvt.u(:) - wvt.u_g(:) - wvt.u_w(:)))
max(abs(wvt.v(:) - wvt.v_g(:) - wvt.v_w(:)))
max(abs(cv(:) - cv_g(:) - cv_w(:)))

%%
figure
tiledlayout(2,1)
nexttile, plot(t,[real(cv_w(:,end)),imag(cv_w(:,end))])
nexttile, plot(t,[real(cv_g(:,end)),imag(cv_g(:,end))])

%%
figure
plot(t/86400,[real(cv(:,end)),imag(cv(:,end))]), hold on
plot(t/86400,[real(cv_w(:,end)+cv_g(:,end)),imag(cv_w(:,end)+cv_g(:,end))])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Compute the spectrum
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dt = t(2)-t(1);
S = zeros(length(t),wvt.Nz);
S_g = zeros(length(t),wvt.Nz);
S_w = zeros(length(t),wvt.Nz);
[psi,lambda]=sleptap(length(t),3);
% psi = [];
for iDepth = 1:wvt.Nz
    [~, Spp, Snn, Spn] = mspec(dt,cv(:,iDepth),psi,'nodemean');
    S(:,iDepth) = (1/(2*pi))*[flipud(vmean(Snn,2)); vmean(Spp(2:end,:),2)];
    [~, Spp, Snn, Spn] = mspec(dt,cv_g(:,iDepth),psi,'nodemean');
    S_g(:,iDepth) = (1/(2*pi))*[flipud(vmean(Snn,2)); vmean(Spp(2:end,:),2)];
    [omega_p, Spp, Snn, Spn] = mspec(dt,cv_w(:,iDepth),[]);
    S_w(:,iDepth) = (1/(2*pi))*[flipud(vmean(Snn,2)); vmean(Spp(2:end,:),2)];
end
omega = [ -flipud(omega_p(2:end)); omega_p];

S_g( S_g<1e-8 ) = nan;
S_w( S_w<1e-8 ) = nan;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Make a figure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


depths = [0 -200 -1000];
depthIndex = zeros(size(depths));
for iDepth=1:length(depths)
    depthIndex(iDepth) = find(wvt.z<depths(iDepth),1,'last');
end

fig1 = figure('Units', 'points', 'Position', [50 50 600 600]);
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');


tl = tiledlayout(length(depths),1);
for iDepth=1:length(depths)
    nexttile
    semilogy(omega*86400/(2*pi),S(:,depthIndex(iDepth)), 'LineWidth', 2,'Color',0*[1 1 1]), hold on
    xlim([-12 12]), ylim([5e-4 3e3])
    plot(-86400/wvt.inertialPeriod*[1 1],ylim,'LineWidth',1,'Color',0*[1 1 1])
    plot(86400/wvt.inertialPeriod*[1 1],ylim,'LineWidth',1,'Color',0*[1 1 1])
    if iDepth < length(depths)
        xtick([])        
    end
    title(sprintf('%d meters',abs(round(depths(iDepth)))))
end
title(tl,'horizontal velocity power spectrum')
ylabel(tl,'m^2/s^2')
xlabel(tl,'frequency (cycles per day)')

% print('MooringSpectaNotSeparated.eps','-depsc2')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Make a figure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig1 = figure('Units', 'points', 'Position', [50 50 600 600]);
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');

tl = tiledlayout(length(depths),1);
for iDepth=1:length(depths)
    nexttile
    semilogy(omega*86400/(2*pi),S(:,depthIndex(iDepth)), 'LineWidth', 2,'Color',0*[1 1 1]), hold on
    ax = gca; ax.ColorOrderIndex = 1;
    semilogy(omega*86400/(2*pi),S_g(:,depthIndex(iDepth)), 'LineWidth', 2)
    semilogy(omega*86400/(2*pi),S_w(:,depthIndex(iDepth)), 'LineWidth', 2)
    xlim([-5 5]), ylim([5e-4 3e3])
    plot(-86400/wvt.inertialPeriod*[1 1],ylim,'LineWidth',1,'Color',0*[1 1 1])
    plot(86400/wvt.inertialPeriod*[1 1],ylim,'LineWidth',1,'Color',0*[1 1 1])
    if iDepth < length(depths)
        xtick([])        
    end
    title(sprintf('%d meters',abs(round(depths(iDepth)))))
end
title(tl,'horizontal velocity power spectrum')
ylabel(tl,'m^2/s^2')
xlabel(tl,'frequency (cycles per day)')

% print('MooringSpectaSeparated.eps','-depsc2')
