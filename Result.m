clear all
close all
clc

%% PARAMETERS %%
T_fr=800; %K
T_ad= 2000; %K
T_act=1e4; %K
tau_ch=10^-6; %s
tol=0.001;

%% Q7: curve comparison %%
T=linspace (T_fr,T_ad,10000); % values of T_br for the horizontal axis
tau_s=[0.0008,0.001035,0.002,0.00615,0.02]; % values of tau_s to test

chem=(1/tau_ch)*exp(-T_act./T).*(T_ad-T);
plot(T,chem), hold on
for k=1:length(tau_s)
    conv=(T-T_fr)/tau_s(k);
    plot(T,conv)
    ylim([0,10e5])
end


%% Q8: S-curve %%
tau_s=linspace(0.0001,0.01,2000); % values of tau_s to test

chem=(1/tau_ch)*exp(-T_act./T).*(T_ad-T); % chemical term

tau_s_pos=[]; % initialize list of results
T_pos=[]; % initialize list of results

for k=1:length(tau_s)
    conv=(T-T_fr)/tau_s(k); % convective term
    dist=abs(chem-conv); % distance between the two curves
    bool=(dist<tol*chem); % selection of temperature indices where the two curves are close enough (crossings)
    T_pos=[T_pos,T(bool)]; % add the list of steady state temperatures for this tau_s to the previous ones
    tau_s_pos=[tau_s_pos,tau_s(k)*ones(size(T(bool)))]; % add the corresponding tau_s to make the horizontal axis
end

figure;
hold on
plot([tau_s(1),tau_s(end)],[T_fr,T_fr]) % minimum temperature line
plot([tau_s(1),tau_s(end)],[T_ad,T_ad]) % maximum temperature line
scatter(tau_s_pos,T_pos) % possible steady state temperatures
legend('T_{fr}','T_{ad}','Steady state temperatures','Location','east')
ylim([T_fr-100,T_ad+100])
xlabel('Residence time (s)')
ylabel('PSR temperature (K)')


%% Q10: Extinction limits %%
phi=[0.6,0.7,0.8,0.9,1,1.05,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2];
T_ad=[1670,1850,2000,2150,2220,2225,2205,2150,2050,1970,1900,1830,1760,1690,1610,1570];

% resolution of 2nd order equation for T
delta=(T_act*(T_fr+T_ad)).^2+4*(T_act*T_fr.*T_ad.*(T_fr-T_ad-T_act));
T_ign=(-(T_fr*T_act+T_ad.*T_act)+sqrt(delta))./(2*(T_fr-T_ad-T_act)); % colder solution: ignition
T_ext=(-(T_fr*T_act+T_ad.*T_act)-sqrt(delta))./(2*(T_fr-T_ad-T_act)); % hotter solution: extinction

tau_rat_ign=exp(-T_act./T_ign).*(T_ad-T_ign)./(T_ign-T_fr); % ratio of tau_ch/tau_s at ignition
tau_rat_ext=exp(-T_act./T_ext).*(T_ad-T_ext)./(T_ext-T_fr); % ratio of tau_ch/tau_s at extinction


figure;
plot(phi,tau_rat_ign,phi,tau_rat_ext)
legend('Ignition Limit','Extinction Limit')
xlabel('Equivalence ratio (-)')
ylabel('$\frac{\tau_{ch}}{\tau_s}$ (-)','Interpreter','latex','FontSize',18)
