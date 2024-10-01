%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CyprusEddyRedux
%
% Repeats the Cyprus eddy problem (Lelong et. al, (2020)) with slightly
% different stratification and initial conditions
%
% Jeffrey J. Early
% jeffrey@jeffreyearly.com


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Specify the problem dimensions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nxy = 256;
Nz = 30;
lat = 27;

N0 = 3*2*pi/3600; % buoyancy frequency at the surface, radians/seconds
L_gm = 1300; % thermocline exponential scale, meters
Lz = 4000;
N2 = @(z) N0*N0*exp(2*z/L_gm);

im = InternalModesWKBSpectral(N2=N2,zIn=[-Lz 0],latitude=lat);
[~,~,~,k_sd] = im.ModesAtFrequency(2*pi/(12.420602*3600));

k_sd1 = k_sd(1);
L = 5*(2*pi/k_sd1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Spin up a geostrophic field
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('spun-up-qg-field-128.nc','file')
    wvt_qg = WVTransformHydrostatic([L, L, Lz],[128, 128, Nz], N2=N2,latitude=lat);

    % First barolinic mode plus a barotropic mode so that the bottom velocity
    % is zero.
    u0 = 0.10;
    wvt_qg.setGeostrophicModes(k=0,l=5,j=1,phi=0,u=u0);
    wvt_qg.setGeostrophicModes(k=0,l=5,j=0,phi=0,u=max(max(wvt_qg.u(:,:,1))));
    A0Force = wvt_qg.A0;
    wvt_qg.removeAll();

    % Now add some noise
    wvt_qg.addRandomFlow('geostrophic',uvMax=0.005);

    nlFlux = WVNonlinearFluxQGForced(wvt_qg,uv_damp=0.20,r=1/(200*86400));
    nlFlux.setGeostrophicForcingCoefficients(A0Force,tau0=50*86400);
    model = WVModel(wvt_qg,nonlinearFlux=nlFlux);
    model.integrateToTime(2000*86400);

    wvt_qg.writeToFile('spun-up-qg-field-128.nc');
else
    wvt_qg = WVTransform.waveVortexTransformFromFile('spun-up-qg-field-128.nc');
end

%%
if ~exist('spun-up-qg-field-256.nc','file')
    wvt_qg = wvt_qg.waveVortexTransformWithResolution([256 256 Nz]);
    model = WVModel(wvt_qg);
    model.setupIntegrator(relTolerance=1e-3);
    model.integrateToTime(wvt_qg.t + 50*86400);

    wvt_qg.writeToFile('spun-up-qg-field-256.nc');
end

Ekj = wvt_qg.transformToRadialWavenumber( wvt_qg.A0_TE_factor .* (abs(wvt_qg.A0).^2));
figure, plot(wvt_qg.kRadial,sum(Ekj,1)), xlog, ylog, xlabel('k'), ylabel('energy'), title(sprintf('geostrophic energy spectrum, day %d',round(wvt_qg.t/86400)))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Now build the wave forcing
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wvt = WVTransformHydrostatic([L, L, Lz],[Nxy, Nxy, Nz], N2=N2,latitude=lat);

outputVar = WVVariableAnnotation('ssh_w',{'x','y'},'m', 'sea-surface height, wave component');
wvt.addOperation(WVOperation('ssh_w', outputVar,@(wvt) wvt.p_w(:,:,end)/(wvt.rho0*wvt.g)));

outputVar = WVVariableAnnotation('ssh_pv',{'x','y'},'m', 'sea-surface height, pv component');
wvt.addOperation(WVOperation('ssh_pv', outputVar,@(wvt) (wvt.p_mda(:,:,end)+wvt.p_g(:,:,end))/(wvt.rho0*wvt.g)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Select a wave band around the semi-diurnal frequency +/ dPeriod
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dPeriod = 300;
omega_min =2*pi/(12.420602*3600+dPeriod);
omega_max =2*pi/(12.420602*3600-dPeriod);
Omega = wvt.Omega;
omega_sd = Omega > omega_min & Omega < omega_max & wvt.J == 1;
fprintf('Found %d modes within +/- %d seconds of the semi-diurnal period after removing the aliased modes.\n',sum(omega_sd(:)),dPeriod);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Restrict the Garrett-Munk spectrum to semi-diurnal frequencies
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wvt.removeAll;

% wvt.initWithGMSpectrum();
% wvt.Ap(~omega_sd) = 0;
% wvt.Am(~omega_sd) = 0;
% 
% ApTide = 6*wvt.Ap;
% AmTide = 0*wvt.Am;

theta0 = pi*rand(1);
j = 1;
k0 = k_sd1;
L = 75e3;
dTheta = 0.6;
sshAmplitude = 2.2e-2;

% randomize amplitudes, and restrict to the wavenumbers (with a
% gaussian) and vertical mode of interest.
[Ap_rand,~,~] = wvt.waveComponent.randomAmplitudes();
Ap = zeros(wvt.spectralMatrixSize);
Ap(j+1,:) = Ap_rand(j+1,:);
Theta = atan2(wvt.L,wvt.K);

Ap = exp(-((wvt.Kh-k0)*L).^2 - ((Theta-theta0)/dTheta).^2) .* Ap;

% Now rescale to match the requested amplitude.
p = wvt.transformToSpatialDomainWithF(Apm=wvt.NAp.*Ap);
ssh = p(:,:,end);
Ap = Ap*(sshAmplitude/max(abs(ssh(:))))/2;

forcingIndices = find( abs(Ap) > 1e-5*max(abs(Ap(:))) );
ApTide = zeros(wvt.spectralMatrixSize);
AmTide = zeros(wvt.spectralMatrixSize);
ApTide(forcingIndices) = Ap(forcingIndices);
AmTide(forcingIndices) = -Ap(forcingIndices);

wvt.Ap = ApTide;
wvt.Am = AmTide;

fprintf('Tidal amplitudes are umax: %.1f cm/s and ssh_max: %.1f cm\n',wvt.uvMax*100,100*max(abs(wvt.ssh(:))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Select a near-inertial wavenumber, and force near the surface
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wvt.removeAll;

wvt.initWithGMSpectrum();

if 0
    wvt.Ap(wvt.Kh > 3*wvt.dk) = 0;
    wvt.Am(wvt.Kh > 3*wvt.dk) = 0;
    wvt.Ap(wvt.Kh == 0) = 0;
    wvt.Am(wvt.Kh == 0) = 0;

    Ld = 500;
    taper = exp((wvt.Z/Ld));
    wvt.initWithUVEta(wvt.u .* taper, wvt.v .* taper, wvt.eta .* taper );

    wvt.Ap(wvt.Kh > 3*wvt.dk) = 0;
    wvt.Am(wvt.Kh > 3*wvt.dk) = 0;
    wvt.Ap(wvt.Kh == 0) = 0;
    wvt.Am(wvt.Kh == 0) = 0;
else
    wvt.Ap(wvt.Kh > 0) = 0;
    wvt.Am(wvt.Kh > 0) = 0;
    Ld = 500;
    taper = exp((wvt.Z/Ld));
    wvt.initWithUVEta(wvt.u .* taper, wvt.v .* taper, wvt.eta .* taper );
end

ApIO = wvt.Ap;
AmIO = wvt.Am;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Set up a model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wvt.removeAll;
wvt.initWithGMSpectrum();
Ap_forcing = ApTide + ApIO;
Am_forcing = AmTide + AmIO;
wvt.Ap(abs(Ap_forcing) > 0) = Ap_forcing(abs(Ap_forcing) > 0);
wvt.Am(abs(Am_forcing) > 0) = Am_forcing(abs(Am_forcing) > 0);

% Get rid of the inertial stuff...
% wvt.Ap(wvt.Kh == 0) = 0;
% wvt.Am(wvt.Kh == 0) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Add an eddy
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 1
    wvt.A0 = wvt_qg.A0;
else
    % Consider a "shallow eddy" where the density anomaly sits close to the
    % surface. This example was constructed in Early, Hernández-Dueñas, Smith,
    % and Lelong (2024), https://arxiv.org/abs/2403.20269

    Le = 80e3;
    He = wvt.Lz/5;
    U = 0.25; % m/s
    x0 = (max(wvt.x)-min(wvt.x))/2;
    y0 = (max(wvt.y)-min(wvt.y))/2;

    H = @(z) exp(-(z/He/sqrt(2)).^2 );
    F = @(x,y) exp(-((x-x0)/Le).^2 -((y-y0)/Le).^2);
    psi = @(x,y,z) U*(Le/sqrt(2))*exp(1/2)*H(z).*(F(x,y) - (pi*Le*Le/(wvt.Lx*wvt.Ly)));
    wvt.setGeostrophicStreamfunction(psi);
end

wvt.summarizeEnergyContent;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Build flux
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nlFlux = WVNonlinearFluxForced(wvt,uv_damp=wvt.uvMax,r=1/(200*86400));
nlFlux.setWaveForcingCoefficients(ApTide + ApIO,AmTide + AmIO,tauP=10*86400,tauM=10*86400);
% nlFlux.setWaveForcingCoefficients(ApTide + ApIO,AmTide + AmIO);

[energyFrequency,omegaVector] = wvt.convertFromWavenumberToFrequency;
figure, plot(omegaVector/wvt.f,squeeze(sum(energyFrequency,1))), xlog, ylog
figure, pcolor(wvt.ssh.'), shading interp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Run the model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model = WVModel(wvt,nonlinearFlux=nlFlux);

% See how things look
Ekj = wvt.transformToRadialWavenumber( wvt.Apm_TE_factor .* (abs(wvt.Ap).^2 + abs(wvt.Am).^2));
figure(name="Forced inertial-tidal energy and velocity variance")
tl1 = tiledlayout(2,2);
nexttile(tl1), plot(wvt.kRadial,sum(Ekj,1)), xlog, ylog, xlabel('k'), ylabel('energy'), title(sprintf('wave energy spectrum, day %d',round(wvt.t/86400)))
nexttile(tl1), plot(100*sqrt(squeeze(mean(mean(wvt.u.^2 + wvt.v.^2,1),2))),wvt.z), ylabel('z (m)'), xlabel('velocity variance (cm/s)'), title('velocity variance, day 0')

tRecord = wvt.t;
totalEnergy = wvt.waveEnergy;
totalEnstrophy = wvt.totalEnstrophy/wvt.f/wvt.f;
sshVariance = sqrt(mean(wvt.ssh(:).^2));
uvVariance = sqrt(mean(wvt.u(:).^2 + wvt.v(:).^2));
uvMax = wvt.uvMax;
zetaMax = max(abs(wvt.zeta_z(:)))/wvt.f;

figure(name="Forced inertial-tidal surface vorticity during spinup")
tl = tiledlayout("flow",TileSpacing="compact");
title(tl,'surface vorticity')

figure(name="Forced inertial-tidal instability: spinup metrics")
tl2 = tiledlayout(6,1,TileSpacing="tight");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Run the model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

spinupTime = 8*86400;
nIncrements = 9;
dIncrement = spinupTime/(nIncrements-1);
t0 = wvt.t;
for i=0:8
    model.integrateToTime(t0 + i*dIncrement);

    tRecord(end+1) = wvt.t;
    totalEnergy(end+1) = wvt.waveEnergy;
    totalEnstrophy(end+1) = wvt.totalEnstrophy/wvt.f/wvt.f;
    sshVariance(end+1) = sqrt(mean(wvt.ssh(:).^2));
    uvVariance(end+1) = sqrt(mean(wvt.u(:).^2 + wvt.v(:).^2));
    uvMax(end+1) = wvt.uvMax;
    zetaMax(end+1) = max(abs(wvt.zeta_z(:)))/wvt.f;

    nexttile(tl2,1)
    plot(tRecord/86400,totalEnergy)
    ylabel('energy (m^2/s)'), xtick([])
    nexttile(tl2,2)
    plot(tRecord/86400,totalEnstrophy)
    ylabel('enstrophy (f^2)'), xtick([])
    nexttile(tl2,3)
    plot(tRecord/86400,sshVariance)
    ylabel('rms ssh'), xtick([])
    nexttile(tl2,4)
    plot(tRecord/86400,uvVariance)
    ylabel('uv-rms'), xtick([])
    nexttile(tl2,5)
    plot(tRecord/86400,uvMax)
    ylabel('uv-max'), xtick([])
    nexttile(tl2,6)
    plot(tRecord/86400,zetaMax)
    ylabel('zeta_z-max')
    xlabel('time (days)')
       
    nexttile(tl)
    pcolor(wvt.x/1e3, wvt.y/1e3, wvt.zeta_z(:,:,end).'), shading interp
    colormap("gray")
    title(sprintf('%d days',round(wvt.t/86400)))
    xtick([]), ytick([])
    pause(0.1);
end

Ekj = wvt.transformToRadialWavenumber( wvt.Apm_TE_factor .* (abs(wvt.Ap).^2 + abs(wvt.Am).^2));
nexttile(tl1,3), plot(wvt.kRadial,sum(Ekj,1)), xlog, ylog, xlabel('k'), ylabel('energy'), title(sprintf('wave energy spectrum, day %d',round(wvt.t/86400)))
nexttile(tl1,4), plot(100*sqrt(squeeze(mean(mean(wvt.u.^2 + wvt.v.^2,1),2))),wvt.z), ylabel('z (m)'), xlabel('velocity variance (cm/s)'), title(sprintf('velocity variance, day %d',round(wvt.t/86400)))

%%

model.createNetCDFFileForModelOutput("dmd-eddy-tide.nc",shouldOverwriteExisting=0,outputInterval=1800);
model.addNetCDFOutputVariables('ssh_w','ssh','ssh_pv');
model.integrateToTime(t0+80*86400)