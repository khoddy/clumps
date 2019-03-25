% Script to calculate combined errors from replicates and calibration using
% MC approach. Currently based on Kelson 2017 calibration.

clear all; close all

% set initial parameters
m = 0.3:0.01:0.7; %range of D47 for plot
n = 1000; % number of MC iterations per D47


%D47err = 0.0332; %Standard dev! This is measurement error from zeros or standards
D47err = 0.017; %Standard dev! This is measurement error from zeros or standards
Ameas = 0.0417;
Aerr = 0.0013; %Standard dev?
Bmeas = 0.139;
Berr = 0.014; %Standard dev?
T=zeros(n,1);

%% Loops to MC measurement and calibration error
for j = 1:length(m)

    for i = 1:n
    D47 = random('Normal',m(j),D47err);
    A = random('Normal',Ameas,Aerr);
    B = random('Normal',Bmeas,Berr);
    T(i)=sqrt((A*1e6)/(D47-B));
    end

    T = T - 273.15; %Convert to C
    
    % Calculate best estimate of temps (no MC) and stats on MC results
    Tmeas(j) = sqrt((Ameas*1e6)/(m(j)-Bmeas)) - 273.15;
    Tstd(j) = std(T);
    
end

TstdA = Tstd;

a=plot(Tmeas,m,'k');
hold on
b=plot(Tmeas - Tstd,m,'b');
plot(Tmeas + Tstd,m,'b')
c=plot(Tmeas - 2*Tstd,m,'r');
plot(Tmeas + 2*Tstd,m,'r')
xlabel('Temperature (C)')
ylabel('D47 (permil)')
legend([a b c], 'Calculated temperature','1SD error','2SD error')
title('D47 and T relationship with measurement and calibration errors')

figure
plot(Tmeas,2*Tstd,'b');hold on
plot(Tmeas,4*Tstd,'r')
xlabel('Calculated temperature (C)')
ylabel('Error (C)')
legend('1SD error','2SD error','location','northwest')
title('Measurement and calibration errors with temperature')
axis([0 250 0 300])


%% Loops to MC just measurement errors

for j = 1:length(m)
    
    for i = 1:n
    D47 = random('Normal',m(j),D47err);
    T(i)=sqrt((Ameas*1e6)/(D47-Bmeas));
    end

    T = T - 273.15;
    Tmeas(j) = sqrt((Ameas*1e6)/(m(j)-Bmeas)) - 273.15;
    Tstd(j) = std(T);
    
end

TstdB = Tstd;

figure
a=plot(Tmeas,m,'k');
hold on
b=plot(Tmeas - Tstd,m,'b');
plot(Tmeas + Tstd,m,'b')
c=plot(Tmeas - 2*Tstd,m,'r');
plot(Tmeas + 2*Tstd,m,'r')
xlabel('Temperature (C)')
ylabel('D47 (permil)')
legend([a b c], 'Calculated temperature','1SD error','2SD error')
title('D47 and T relationship with measurement errors')

figure
plot(Tmeas,2*Tstd,'b');hold on
plot(Tmeas,4*Tstd,'r')
xlabel('Calculated temperature (C)')
ylabel('Error (C)')
legend('1SD error','2SD error','Location','northwest')
title('Measurement errors with temperature')
axis([0 250 0 300])

%% Loops to MC just Calibration error

for j = 1:length(m)

    for i = 1:n
    A = random('Normal',Ameas,Aerr);
    B = random('Normal',Bmeas,Berr);
    T(i)=sqrt((A*1e6)/(m(j)-B));
    end

    T = T - 273.15;
    Tmeas(j) = sqrt((Ameas*1e6)/(m(j)-Bmeas)) - 273.15;
    Tstd(j) = std(T);
        
end
 
TstdC = Tstd;

figure
a=plot(Tmeas,m,'k');
hold on
b=plot(Tmeas - Tstd,m,'b');
plot(Tmeas + Tstd,m,'b')
c=plot(Tmeas - 2*Tstd,m,'r');
plot(Tmeas + 2*Tstd,m,'r')
xlabel('Temperature (C)')
ylabel('D47 (permil)')
legend([a b c], 'Calculated temperature','1SD error','2SD error')
title('D47 and T relationship with calibration errors')

figure
plot(Tmeas,2*Tstd,'b');hold on
plot(Tmeas,4*Tstd,'r')
xlabel('Calculated temperature (C)')
ylabel('Error (C)')
legend('1SD error','2SD error','Location','northwest')
title('Calibration errors with temperature')
axis([0 250 0 300])
 
%% Comparison Figure
p=polyfit(Tmeas,TstdA,2);
errFitY=polyval(p,Tmeas);


figure
a = plot(Tmeas,TstdA);
hold on
b = plot(Tmeas,TstdB);
c = plot(Tmeas,TstdC);
d = plot(Tmeas,errFitY,'LineWidth',2);
xlabel('Calculated temperature (C)')
ylabel('Error (C)')
legend('Combined','Analytical','Calibration','Best Fit','Location','northwest')
title('Comparison of error magnitudes (1 SD)')
axis([0 250 0 80])

    
    
