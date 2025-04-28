
%Below command lets you see variables stored in file
%ncdisp(file); 


file = "C:\Users\benny\Downloads\2 step pulse 2w 3w 25TW.exo";
figures = 1;
setlimits = 0;
Measured = 1;


%Import data from file. 
Time  = ncread(file,'time_whole');
Radius  = ncread(file,'zone_boundaries');
Density  = ncread(file,'mass_density');
ElectronDensity  = ncread(file,'elec_density');
IonTemp  = ncread(file,'ion_temperature');
ElecTemp  = ncread(file,'elec_temperature');
Mass=ncread(file,'zone_mass');
Pressure=ncread(file, 'ion_pressure')+ncread(file, 'elec_pressure'); %in J/ cm^3
Velocity=ncread(file, 'fluid_velocity')./100000; %Converted to km/s from cm/s
Volume = Mass./Density;
try
Flux = ncread(file, 'Freq. resolved net flux at region interfaces [J per (cm2.sec.eV)]');
Groups = ncread(file, 'photon_energy_group_bounds');

LaserEnergy=ncread(file,'LaserEnDeliveredTimeInt');
LaserPower=[0; diff(LaserEnergy)./diff(Time)];
SimulationPower = LaserPower;




DepositedLaserPowerZones = ncread(file,'LaserPwrSrc'); %Laser energy deposition rates [J per (g.sec)].
DepositedLaserPower = sum(DepositedLaserPowerZones.*Mass);
DepositedEnergy1 = trapz(Time, DepositedLaserPower);
DepositedEnergy = ncread(file,'LaserEnTimeIntg'); %Time-integrated laser energy deposition [J per (cm**X.zone)] .
DepositedEnergy = ncread(file,'LaserEnTimeIntg').*Volume; %Time-integrated laser energy deposition [J per (cm**X.zone)] .
DepositedEnergy = sum(DepositedEnergy(:, end));


%Laser statistics calculated. Hyades has a 1cm^2 area, so multiply by area
%in cm^2 to get actual power. This also means that laser power is actually
%Laser Intensity
AdjustedEnergy = (LaserEnergy) * pi()*(0.06/2)^2 %Multiply by area to get real energy
LaserIntensity = LaserPower; %This is intensity
RealLaserPower = LaserIntensity * pi()*(0.06/2)^2 %Multiply by area to get real power
DepositedLaserPower = (DepositedLaserPower) * pi() * (0.02)^2;
MaxLaserIntensity = max(LaserIntensity);
catch
end
IonTempKelvin = 1.16E4 .*IonTemp;

% FoamDensity = 0.253;
FoamDensity = Density(end,1);

%Define shock front as last zone in increasing radius in the sim where
%there is a substantial difference between density now and original density
%- since density doesn't really change until shock reaches it. Doesn't work
%so well for gold/ start of al where there is preheat.
Shock = (Density - repmat(Density(:,1), 1, length(Time))) > 0.1;
%Shock = (Density ./ repmat(Density(:,1), 1, length(Time))) > 1.5;

%To avoid it picking up changing density in aluminised plastic files, I
%have also specified that these layers shouldn't be included.
FoamMaxRadius = max(Radius(:,1));
BulkRegions = Radius(1:end-1,1)<FoamMaxRadius;
Shock = Shock.*BulkRegions;
ShockIndex = arrayfun(@(x)find(Shock(:,x),1,'last'),1:size(Shock,2), 'UniformOutput',false);
tf = cellfun('isempty',ShockIndex); % true for empty cells
ShockIndex(tf) = {1};
ShockIndex= [ShockIndex{:}];
ShockIndexReshaped=sub2ind(size(Radius),ShockIndex,1:length(Time));
ShockRadius=Radius(ShockIndexReshaped);

InterfaceIndex = find(abs(Density(2:end,1) - Density(1:end-1))>1, 1, 'first');
InterfaceRadius=Radius(InterfaceIndex, :);

% ShockVelocity = smooth(diff(ShockRadius)./diff(Time.'.*10^9.'), 30);
% InterfaceVelocity = diff(InterfaceRadius)./diff(Time.'.*10^9.');
% ParticleVelocity = Velocity(InterfaceIndex, :)./10000;

ShockVelocity = smooth(diff(ShockRadius)./diff(Time.'), 10); %cm/s
InterfaceVelocity = diff(InterfaceRadius)./diff(Time.'); %cm/s
ParticleVelocity = Velocity(InterfaceIndex, :); %km/s - conversion performed when imported

figure;
DensityPlot = surf(Time.*10^9,Radius(1:end-1,:).*10000, Density);
%xlim([(MaxDensityTime-.5e-9) (MaxDensityTime+.5e-9)]);
%xlim([xlowerlim xupperlim]);
%ylim([0 0.01]);
%zlim([0 MaxDensity]);
colorbar;
colormap(jet);
title('Density Plot CH+Al+SiO2+Aerogel')
xlabel('Time (ns)');
ylabel('Radius (\mu m)')
zlabel('Mass Density (g/cm^3)');
shading interp
view(2)
set(gca,'colorscale','log')
caxis([0.01 20]);
ylim([-100 1200])
xlim([0 10])
% hold on
% plot3(Time.*10^9, ShockRadius.*10000, repmat(max(max(Density)), 1, length(Time)))
% plot3(Time.*10^9, InterfaceRadius.*10000, repmat(max(max(Density)), 1, length(Time)))
% hold off

figure;
DensityPlot = surf(Time.*10^9,Radius(1:end-1,:).*10000, Density);
%xlim([(MaxDensityTime-.5e-9) (MaxDensityTime+.5e-9)]);
%xlim([xlowerlim xupperlim]);
%ylim([0 0.01]);
%zlim([0 MaxDensity]);
colorbar;
colormap(jet);
title('Density Plot CH+Al+SiO2+Aerogel')
xlabel('Time (ns)');
ylabel('Radius (\mu m)')
zlabel('Mass Density (g/cm^3)');
shading interp
view(2)
set(gca,'colorscale','log')
caxis([0.01 20]);
ylim([-100 1200])
xlim([0 10])
hold on
plot3(Time.*10^9, ShockRadius.*10000, repmat(max(max(Density)), 1, length(Time)))
plot3(Time.*10^9, InterfaceRadius.*10000, repmat(max(max(Density)), 1, length(Time)))
hold off

figure
plot(Time.*10^9, ShockRadius.*10000)
hold on
plot(Time.*10^9, InterfaceRadius.*10000)
hold off
xlabel('Time (ns)');
ylabel('Radius (\mu m)')
legend('Shock radius', 'Interface radius')

% figure
% plot(Time(1:end-1).*10^9, ShockVelocity.*1000)
% hold on
% plot(Time(1:end-1).*10^9, InterfaceVelocity.*1000)
% hold off
% xlabel('Time (ns)');
% ylabel('Velocity (m/s)')
% ylim([0 0.5])
% legend('Shock velocity', 'Interface velocity')

figure
plot(Time(1:end-1).*10^9, smooth(ShockVelocity./(10^5)))
hold on
plot(Time(1:end-1).*10^9, InterfaceVelocity./(10^5))
hold off
xlabel('Time (ns)');
ylabel('Velocity (km/s)')
legend('Shock velocity', 'Interface velocity')


figure
plot(Time.*10^9, ParticleVelocity)
hold on
plot(Time(1:end-1).*10^9, InterfaceVelocity./(10^5))
hold off
xlabel('Time (ns)');
ylabel('Velocity (km/s)')
% ylim([0 5e-4])
legend('Particle velocity', 'Interface velocity')

figure;
DensityPlot = surf(Time.*10^9,Radius(1:end-1,:).*10000, IonTempKelvin);
colorbar;
colormap(jet);
title('Ion Temp Plot');
xlabel('Time (ns)');
ylabel('Radius (cm)');
zlabel('Mass Density (g/cm^3)');
shading interp
set(gca,'colorscale','log')
view(2)
ylim([0 1295])


%Plot zoning
figure
MassDifference = 100*diff(Mass(:,1))./Mass(2:end,1);
plot(MassDifference)
ylabel('Percentage Mass Difference')
yyaxis right;
plot(Mass(:,1), ':');
ylabel('Zone Mass')
title('Zoning Plot')
xlabel('Zone')

figure
title('Laser Power and input energy')
yyaxis left
plot(Time.*10^9, LaserIntensity./10^12)
ylabel('Laser Power (TW)')
yyaxis right
plot(Time.*10^9, AdjustedEnergy)
xlabel('Time (ns)')
ylabel('Laser Energy (J)')






