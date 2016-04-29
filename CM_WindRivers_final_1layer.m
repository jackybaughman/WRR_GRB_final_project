%% Wind River mountains and green basin evolution
% incorporates rock path over time with temperature
% uses hillslope diffusion and geothermal gradient
% vertical thrust fault
% Computer Modeling, 2016, Final project, JSB

clear all
figure(3)
clf
%% Initialize 

% Constants and variables
Rho_r = 2700; % rock density kg/m3
Rho_s = 1600; % soil density kg/m3
kappa=0.002; %m^2/yr 
k=kappa*Rho_s; % efficiency 
wgdotnaught = 100 * (10^-5); % m/year initial weathering rate granite basement
wmdotnaught = 100 *(10^-5); % m/yr initial weathering rate meso cover
Hstar = 2; % m
Qm = .07; % heat flux
kT=2; % T conductivity
Ts = 15; % surface T

% Thrust faulting - 3 time slices
thrust = 4 * (10^-5); % total thrust1
uplift = .75*thrust; % m/yr uplift rate;
downlift = -.25*thrust; % downlift rate;

thrust2 = 33 * (10^-5);% total thrust 2
uplift2 = .75*thrust2; % m/yr uplift rate;
downlift2 = -.25*thrust2; % downlift rate;

thrust3 = 1 * (10^-5);% total thrust 3
uplift3 = .75*thrust3; % m/yr uplift rate;
downlift3 = -.25*thrust3; % downlift rate;

xstar = 1*300; %km - deflection parameter

% Arrays

% z array
dz = 1; % z step m
zmax = 4000; % zmax m
zT = 0:dz:zmax;% z array
T = ones(size(zT))*Ts + ((Qm*zT)/kT); %initial T based on geotherm

% time array setup
dt = 100; % time step (years)
tmax = 35000 * 1000; % max time (ka)
t = 0:dt:tmax; % time array
imax = length(t)-200000; %thrust1
jmax = length(t)-100000; %thrust2
kmax = length(t); %thrust3

% distance array setup
dx = 10; % distance step (10s m)
xmax = 2000; % max dist (m) % min dist (10s m)
x = 0:dx:xmax; % dist array

% Initial condtions

% Initial mobile regolith
Hnaught = 0; % H at time zero 
H = Hnaught * ones(size(x)); %initial H

% Initial topography 
zg = 100*rand(size(x))-500; % randomized initial granite topo
zm = zg + 750; % initial meso cover
z = zg+(zm-zg)+H; % Initial hill slope topo

% Sea level constants
meansl = 0; % mean sea level (m)
A = 80; % sea level amplitude (m)
P = 47 * 1000; % sea level change period (ka)
sealevel = zeros(size(x)); % initial array size x

% initial rock package I care about
rock = -3000*ones(size(x));
rock(120)= -2850;
rockinitial = rock(120);
surfaceinitial = z(120);

% plot animation
n=150; %number of plots
tplot = tmax/n;


%% RUN

% In the time loop

% uplift1
for i = 1:imax;
% sealevel - sinusoidal change with time
sealevel = 0*x+meansl+A*sin(2*pi*(t(i)/P)); % absolute change of sea level sinusoidal    

% rock uplift
rock(120) = rock(120)+(uplift*dt); % absolute rock uplift

% surface uplift  
%granite uplift
zg(100:140)=zg(100:140)+(uplift*dt); %thrust uplift
zg(1:100) = zg(1:100)+downlift*dt*exp((x(1:100)-1000)/xstar);%thrust down + deflection

%meso uplift
zm(100:140)=zm(100:140)+(uplift*dt); 
zm(1:100) = zm(1:100)+downlift*dt*exp((x(1:100)-1000)/xstar);

% weathering rate
wgdot = wgdotnaught*exp(-(H/Hstar)); % weathering flux exponential function granite
wmdot = wmdotnaught*exp(-(H/Hstar)); % weathering flux exponential function meso

% calculate slope at each point
dzdx = diff(z)/dx; % slope

% calculate Q and diff Q
Q = dzdx.*-k.*(1-exp(-(H(2:end)./Hstar))); % flux inputs
Q = [Q(1) Q Q(end)]; % filling in the Q array
dQdx = diff(Q)/dx; % rate of change of flux

% Determine mobile regolith over time
dHdtg = ((Rho_r/Rho_s)*wgdot)-((1/Rho_s)*dQdx); % rate of change of mobile regolith granite
dHdtm = ((Rho_r/Rho_s)*wmdot)-((1/Rho_s)*dQdx); % rate of change of mobile regolith meso

% Decide which layer to turn to mobile regolith - only the surface exposed layer
if zm > zg % when meso > granite
    H = H + (dHdtm*dt);  % erode meso
else
    H = H + (dHdtg*dt); % if granite>meso, erode granite
end

H = max(H,0); % cannot have negative mobile regolith

% Calculate bedrock over time
% decide which layer to weather
if zm > zg % meso>granite
    zm = zm-(wmdot*dt); % change of bedrock due to weathering
else % granite>meso
    zg = zg-(wgdot*dt); % change of bedrock due to weathering
    zm = zm-(wmdot*dt); % change of bedrock due to weathering

end

zm = max(zm,zg); % can't have negative meso

% recalculate hillslope topo
z = zg+(zm-zg)+H; % hillslope topography

% Find where lake would form
water = zeros(size(x)); % water size of x array
waterfill = find(sealevel>z); % find where sealevel > surface
water(waterfill)=sealevel(waterfill); % water = sealevel where > surface

% rock package temperature at changing times 
rockz=z(120)-rock(120); % distance of rock package from surfaace
rockz = round(rockz); % round that distance to the nearest integer
rockT = T(abs(rockz)); % calculate T at that distance from the surface


%% PLOT
if (rem(t(i),tplot)==0) % decrease plotting frequency - speed up animation
figure(3)
subplot('position',[.1 .5 .8 .45]);
plot(x/10,zm,'Color',[0.1 0.5 0.1],'linewidth',3)
hold on
plot(x/10,z,'r','linewidth',2) % plot hill slope over time
hold on
plot(x/10,zg,'k','linewidth',3) % plot bedrock over time
hold on
plot(x/10,sealevel, '-.b', 'linewidth',2)
hold on
plot((x(120))/10,rock(120),'o','markerfacecolor','r','markersize',12)

hold off

   xlabel('Distance (km)','fontname','arial','fontsize',21) % x label
   ylabel('Elevation (m)','fontname','arial','fontsize',21) % y label
   set(gca,'fontsize',18,'fontname','arial') % axes number labels
   title(['Mountain and basin evolution at ',num2str(75-(t(i)/1000000)),' Ma ']) % title - accumulates model time
   axis([0 1350/10 -3000 3800]) % hold axes constant
   legend('Mesozoic rock','Regolith','Granitic rock','Sealevel','Rock package','location','northwest')
   
        %Subplot of Rock T
        subplot('position',[.1 .07 .8 .1]);
        plot(T,zT,'k','linewidth',2)
        hold on
        plot(rockT,rockz, 'o','markerfacecolor','r','markersize',10)
        hold off
        title('Rock package temperature with time','fontname','arial','fontsize',18);
        xlabel('T (C)','fontname','arial','fontsize',16);
        ylabel('Depth (m)','fontname','arial','fontsize',16);
        set(gca,'fontsize',14,'fontname','arial') % axes number labels
        set(gca,'YDIR','reverse')
        axis([0 200 0 4000])
        pause(0.02)
        
        
        %Subplot of flux
        subplot('position',[.1 .27 .8 .1]);
        plot(x/10,Q(1:end-1),'k','linewidth',2)
        title('Flux with time','fontname','arial','fontsize',18);
        xlabel('distance along profile,m','fontname','arial','fontsize',16);
        ylabel('flux,m^2','fontname','arial','fontsize',16);
        set(gca,'fontsize',14,'fontname','arial') % axes number labels
        axis([0 1350/10 0 10])
        hold all;
        pause(0.02)
        
        hold off
   
end

end


% uplift2
for i = imax:jmax;
sealevel = 0*x+meansl+A*sin(2*pi*(t(i)/P)); % absolute change of sea level sinusoidal    

rock(120) = rock(120)+(uplift2*dt);

%uplift    
zg(100:140)=zg(100:140)+(uplift2*dt);
zg(1:100) = zg(1:100)+downlift2*dt*exp((x(1:100)-1000)/xstar);
%zg(140:end) = zg(140:end)+downlift*dt*exp((-x(140:end)+1400)/xstar);

zm(100:140)=zm(100:140)+(uplift2*dt);
zm(1:100) = zm(1:100)+downlift2*dt*exp((x(1:100)-1000)/xstar);
%zm(140:end) = zm(140:end)+downlift*dt*exp((-x(140:end)+1400)/xstar);
% weathering rate

wgdot = wgdotnaught*exp(-(H/Hstar)); % weathering flux exponential function granite
wmdot = wmdotnaught*exp(-(H/Hstar)); % weathering flux exponential function meso

% calculate slope at each point
dzdx = diff(z)/dx; % slope

% calculate Q and diff Q
%Q = dzdx*-k;
%Q = dzdx.*-k.*H(2:end);
Q = dzdx.*-k.*(1-exp(-(H(2:end)./Hstar))); % flux inputs
Q = [Q(1) Q Q(end)]; % filling in the Q array
dQdx = diff(Q)/dx; % rate of change of flux

% Determine mobile regolith over time
dHdtg = ((Rho_r/Rho_s)*wgdot)-((1/Rho_s)*dQdx); % rate of change of mobile regolith granite
dHdtm = ((Rho_r/Rho_s)*wmdot)-((1/Rho_s)*dQdx); % rate of change of mobile regolith meso

if zm > zg
    H = H + (dHdtm*dt); 
else
    H = H + (dHdtg*dt);
end

H = max(H,0); % cannot have negative mobile regolith

% Calculate bedrock over time

if zm > zg
    zm = zm-(wmdot*dt); % change of bedrock due to weathering
else
    zg = zg-(wgdot*dt); % change of bedrock due to weathering
    zm = zm-(wmdot*dt); % change of bedrock due to weathering

end

zm = max(zm,zg);
% recalculate hillslope topo
z = zg+(zm-zg)+H; % hillslope topography

water = zeros(size(x));
waterfill = find(sealevel>z);
water(waterfill)=sealevel(waterfill);


rockz=(z(120))-(rock(120));
rockz = round(rockz);
rockT = T(abs(rockz));
%% PLOT
if (rem(t(i),tplot)==0) % decrease plotting frequency - speed up animation
figure(3)
subplot('position',[.1 .5 .8 .45]);
plot(x/10,zm,'Color',[0.1 0.5 0.1],'linewidth',3)
hold on
plot(x/10,z,'r','linewidth',2) % plot hill slope over time
hold on
plot(x/10,zg,'k','linewidth',3) % plot bedrock over time
hold on
plot((x(120))/10,rock(120),'o','markerfacecolor','r','markersize',12)
hold on
plot((x(waterfill))/10,water(waterfill), '-.b', 'linewidth',2)
hold off

   xlabel('Distance (km)','fontname','arial','fontsize',21) % x label
   ylabel('Elevation (m)','fontname','arial','fontsize',21) % y label
   set(gca,'fontsize',18,'fontname','arial') % axes number labels
   title(['Mountain and basin evolution at ',num2str(75-(t(i)/1000000)),' Ma ']) % title - accumulates model time
   axis([0 1350/10 -3000 3800]) % hold axes constant
   legend('Mesozoic rock','Regolith','Granitic rock','Rock package','location','northwest')
   
        %Subplot of Rock T
        subplot('position',[.1 .07 .8 .1]);
        plot(T,zT,'k','linewidth',2)
        hold on
        plot(rockT,rockz, 'o','markerfacecolor','r','markersize',10)
        hold off
        title('Rock package temperature with time','fontname','arial','fontsize',18);
        xlabel('T (C)','fontname','arial','fontsize',16);
        ylabel('Depth (m)','fontname','arial','fontsize',16);
        set(gca,'fontsize',14,'fontname','arial') % axes number labels
        set(gca,'YDIR','reverse')
        axis([0 200 0 4000])
        pause(0.02)
        
        
        %Subplot of flux
        subplot('position',[.1 .27 .8 .1]);
        plot(x/10,Q(1:end-1),'k','linewidth',2)
        title('Flux with time','fontname','arial','fontsize',18);
        xlabel('distance along profile,m','fontname','arial','fontsize',16);
        ylabel('flux,m^2','fontname','arial','fontsize',16);
        set(gca,'fontsize',14,'fontname','arial') % axes number labels
        axis([0 1350/10 0 10])
        hold all;
        pause(0.02)
        
        hold off
end


end


%uplift3
for i = jmax:kmax;
sealevel = 0*x+meansl+A*sin(2*pi*(t(i)/P)); % absolute change of sea level sinusoidal    

rock(120) = rock(120)+(uplift3*dt);

%uplift    
zg(100:140)=zg(100:140)+(uplift3*dt);
zg(1:100) = zg(1:100)+downlift3*dt*exp((x(1:100)-1000)/xstar);
%zg(140:end) = zg(140:end)+downlift*dt*exp((-x(140:end)+1400)/xstar);

zm(100:140)=zm(100:140)+(uplift3*dt);
zm(1:100) = zm(1:100)+downlift3*dt*exp((x(1:100)-1000)/xstar);
%zm(140:end) = zm(140:end)+downlift*dt*exp((-x(140:end)+1400)/xstar);
% weathering rate

wgdot = wgdotnaught*exp(-(H/Hstar)); % weathering flux exponential function granite
wmdot = wmdotnaught*exp(-(H/Hstar)); % weathering flux exponential function meso

% calculate slope at each point
dzdx = diff(z)/dx; % slope

% calculate Q and diff Q
%Q = dzdx*-k;
%Q = dzdx.*-k.*H(2:end);
Q = dzdx.*-k.*(1-exp(-(H(2:end)./Hstar))); % flux inputs
Q = [Q(1) Q Q(end)]; % filling in the Q array
dQdx = diff(Q)/dx; % rate of change of flux

% Determine mobile regolith over time
dHdtg = ((Rho_r/Rho_s)*wgdot)-((1/Rho_s)*dQdx); % rate of change of mobile regolith granite
dHdtm = ((Rho_r/Rho_s)*wmdot)-((1/Rho_s)*dQdx); % rate of change of mobile regolith meso

if zm > zg
    H = H + (dHdtm*dt); 
else
    H = H + (dHdtg*dt);
end

H = max(H,0); % cannot have negative mobile regolith

% Calculate bedrock over time

if zm > zg
   zm = zm-(wmdot*dt); % change of bedrock due to weathering
else
   zg = zg-(wgdot*dt); % change of bedrock due to weathering
   zm = zm-(wmdot*dt); % change of bedrock due to weathering

end

zm = max(zm,zg);
% recalculate hillslope topo
z = zg+(zm-zg)+H; % hillslope topography

water = zeros(size(x));
waterfill = find(sealevel>z);
water(waterfill)=sealevel(waterfill);

rockz=(z(120))-(rock(120));
rockz = round(rockz);
rockT = T(abs(rockz));    


%% PLOT
if (rem(t(i),tplot)==0) % decrease plotting frequency - speed up animation
figure(3)
subplot('position',[.1 .5 .8 .45]);
plot(x/10,zm,'Color',[0.1 0.5 0.1],'linewidth',3)
hold on
plot(x/10,z,'r','linewidth',2) % plot hill slope over time
hold on
plot(x/10,zg,'k','linewidth',3) % plot bedrock over time
hold on
plot((x(120))/10,rock(120),'o','markerfacecolor','r','markersize',12)
hold on
plot((x(waterfill))/10,water(waterfill), '-.b', 'linewidth',2)


hold off

   xlabel('Distance (km)','fontname','arial','fontsize',21) % x label
   ylabel('Elevation (m)','fontname','arial','fontsize',21) % y label
   set(gca,'fontsize',18,'fontname','arial') % axes number labels
   title(['Mountain and basin evolution at ',num2str(75-(t(i)/1000000)),' Ma ']) % title - accumulates model time
   axis([0 1350/10 -3000 3800]) % hold axes constant
   legend('Mesozoic rock','Regolith','Granitic rock','Rock package','location','northwest')
   
        %Subplot of Rock T
        subplot('position',[.1 .07 .8 .1]);
        plot(T,zT,'k','linewidth',2)
        hold on
        plot(rockT,rockz, 'o','markerfacecolor','r','markersize',10)
        hold off
        title('Rock package temperature with time','fontname','arial','fontsize',18);
        xlabel('T (C)','fontname','arial','fontsize',16);
        ylabel('Depth (m)','fontname','arial','fontsize',16);
        set(gca,'fontsize',14,'fontname','arial') % axes number labels
        set(gca,'YDIR','reverse')
        axis([0 200 0 4000])
        pause(0.02)
        
        
        %Subplot of flux
        subplot('position',[.1 .27 .8 .1]);
        plot(x/10,Q(1:end-1),'k','linewidth',2)
        title('Flux with time','fontname','arial','fontsize',18);
        xlabel('distance along profile,m','fontname','arial','fontsize',16);
        ylabel('flux,m^2','fontname','arial','fontsize',16);
        set(gca,'fontsize',14,'fontname','arial') % axes number labels
        axis([0 1350/10 0 10])
        hold all;
        pause(0.02)
        
        hold off
end
end
rockfinal = rock(120);
rockuplift = rockfinal-rockinitial

surfacefinal = z(120);
surfaceuplift = surfacefinal-surfaceinitial

exhumation = rockuplift-surfaceuplift