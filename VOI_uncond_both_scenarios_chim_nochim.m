close all; clear all; clc;

rng('default');
r_inj = 1e7;
T_total = 40; % Total planned time period of injection
drho = 400;
%eta = log(1e5)*ones(1,9);
eta = [log(1e5)*ones(1,8) log(5e5)]; %for strong cap rock
tau = 0.5*ones(1,9);
n = 10000;
m = length(eta); %number of stratigraphic layers
f = 0.5; % Fraction of injected volume migrating to layer 9 via chimney
%radii = compute_radii(n);
load radii1.mat; load radii2.mat;
[T1,C_th1,h1,a1,b1,c1] = event_times_Sleipner(r_inj, drho, eta, tau, n/2, radii1);
[T2,C_th2,h2,a2,b2,c2] = event_times_Sleipner_chim19(r_inj, drho, eta, tau, n/2, f, radii2);
T=[T1;T2]; C_th=[C_th1;C_th2]; h=[h1;h2]; a=[a1;a2]; b=[b1;b2]; c=[c1;c2];

sn_high = 4; %standard deviation of high noise data
sn_low = 2; %standard deviation of low noise data
thres_lq = 10; %minimum height that can be detected by low quality data
thres_hq = 1; %minimum height that can be detected by high quality data

survey_time = 1:1:40;

v_cont = nan(n,numel(survey_time));
v_stop = nan(1,numel(survey_time));

v_cont_fit_perfect_hq = nan(n,numel(survey_time));
v_cont_fit_perfect_lq = nan(n,numel(survey_time));
v_cont_fit_hq_low = nan(n,numel(survey_time));
v_cont_fit_hq_high = nan(n,numel(survey_time));
v_cont_fit_lq_high = nan(n,numel(survey_time));

VOI_perfect = nan(1,numel(survey_time));
VOI_perfect_hq = nan(1,numel(survey_time));
VOI_perfect_lq = nan(1,numel(survey_time));
VOI_hq_low = nan(1,numel(survey_time));
VOI_hq_high = nan(1,numel(survey_time));
VOI_lq_high = nan(1,numel(survey_time));

inst_h = cell(1,numel(survey_time)); %instantaneous heights
perfect_data_hq = cell(1,numel(survey_time));
perfect_data_lq = cell(1,numel(survey_time));
hq_data_low_noise = cell(1,numel(survey_time));
hq_data_high_noise = cell(1,numel(survey_time));
lq_data_high_noise = cell(1,numel(survey_time));
coeffs_perfect_hq = cell(1,numel(survey_time));
coeffs_perfect_lq = cell(1,numel(survey_time));
coeffs_hq_low = cell(1,numel(survey_time));
coeffs_hq_high = cell(1,numel(survey_time));
coeffs_lq_high = cell(1,numel(survey_time));

high_noise = sn_high*randn(n,m);
low_noise = sn_low*randn(n,m);

for j = 1:numel(survey_time)
    inst_h{j} = h.*(survey_time(j) > T(:,2:end));
    for i=1:n/2
        if any(T(i,:)>survey_time(j))
            ind = find(T(i,:)>survey_time(j), 1);
            V_t = r_inj * (survey_time(j) - T(i,ind-1));
            inst_h{j}(i,ind-1) = compute_h(V_t,a(i,ind-1),b(i,ind-1),c(i,ind-1),0,c(i,ind-1));
        end
    end
    for i=(n/2+1):n
        if any(T(i,:)>survey_time(j))
            ind = find(T(i,:)>survey_time(j), 1);
            if ind == 10
                V_t_9 = r_inj*(survey_time(j)-T(i,ind-1)) + f*r_inj*T(i,ind-1);
                inst_h{j}(i,9) = compute_h(V_t_9,a(i,9),b(i,9),c(i,9),0,c(i,9));
            else
                V_t = (1-f)*r_inj*(survey_time(j)-T(i,ind-1));
                V_t_9 = f*r_inj*survey_time(j);
                inst_h{j}(i,ind-1) = compute_h(V_t,a(i,ind-1),b(i,ind-1),c(i,ind-1),0,c(i,ind-1));
                inst_h{j}(i,9) = compute_h(V_t_9,a(i,9),b(i,9),c(i,9),0,c(i,9));
            end
        end
    end
    
    perfect_data_hq{j} = max(0, 0.*(inst_h{j} < thres_hq) + inst_h{j} .* (inst_h{j} >= thres_hq));
    perfect_data_lq{j} = max(0, 0.*(inst_h{j} < thres_lq) + inst_h{j} .* (inst_h{j} >= thres_lq));
    hq_data_low_noise{j} = max(0, low_noise.*(inst_h{j} < thres_hq) + (inst_h{j} + low_noise) .* (inst_h{j} >= thres_hq));
    hq_data_high_noise{j} = max(0, high_noise.*(inst_h{j} < thres_hq) + (inst_h{j} + high_noise) .* (inst_h{j} >= thres_hq));
    lq_data_high_noise{j} = max(0, high_noise.*(inst_h{j} < thres_lq) + (inst_h{j} + high_noise) .* (inst_h{j} >= thres_lq));
    
    [v_cont(:,j), v_stop(j)] = values(T(:,end), 2, 5, 3, r_inj, survey_time(j), T_total);
    
    %Linear regression
    coeffs_perfect_hq{j} = regress(v_cont(:,j),[ones(n,1),perfect_data_hq{j}]);
    coeffs_perfect_lq{j} = regress(v_cont(:,j),[ones(n,1),perfect_data_lq{j}]);
    coeffs_hq_low{j} = regress(v_cont(:,j),[ones(n,1),hq_data_low_noise{j}]);
    coeffs_hq_high{j} = regress(v_cont(:,j),[ones(n,1),hq_data_high_noise{j}]);
    coeffs_lq_high{j} = regress(v_cont(:,j),[ones(n,1),lq_data_high_noise{j}]);
    v_cont_fit_perfect_hq(:,j) = [ones(n,1),perfect_data_hq{j}]*coeffs_perfect_hq{j};
    v_cont_fit_perfect_lq(:,j) = [ones(n,1),perfect_data_lq{j}]*coeffs_perfect_lq{j};
    v_cont_fit_hq_low(:,j) = [ones(n,1),hq_data_low_noise{j}]*coeffs_hq_low{j};
    v_cont_fit_hq_high(:,j) = [ones(n,1),hq_data_high_noise{j}]*coeffs_hq_high{j};
    v_cont_fit_lq_high(:,j) = [ones(n,1),lq_data_high_noise{j}]*coeffs_lq_high{j};
    
    PV = max(mean(v_cont(:,j)),v_stop(j));
    PoV_perfect_hq = mean(max(v_cont_fit_perfect_hq(:,j),v_stop(j)));
    PoV_perfect_lq = mean(max(v_cont_fit_perfect_lq(:,j),v_stop(j)));
    PoV_hq_low = mean(max(v_cont_fit_hq_low(:,j),v_stop(j)));
    PoV_hq_high = mean(max(v_cont_fit_hq_high(:,j),v_stop(j)));
    PoV_lq_high = mean(max(v_cont_fit_lq_high(:,j),v_stop(j)));
    VOI_perfect_hq(j) = PoV_perfect_hq - PV;
    VOI_perfect_lq(j) = PoV_perfect_lq - PV;
    VOI_hq_low(j) = PoV_hq_low - PV;
    VOI_hq_high(j) = PoV_hq_high - PV;
    VOI_lq_high(j) = PoV_lq_high - PV;
    
    fprintf('Completed computation for survey time = %f years\n',survey_time(j));
end

tt = 1:numel(survey_time);
tt2 = tt;
figure(1); plot(survey_time(tt),VOI_perfect_hq(tt),'r'); hold on;
plot(survey_time(tt),VOI_perfect_lq(tt),'m');
plot(survey_time(tt),VOI_hq_low(tt),'b');
plot(survey_time(tt),VOI_hq_high(tt),'g');
plot(survey_time(tt),VOI_lq_high(tt),'c');
legend('Perfect data HQ','Perfect data LQ','HQ Low noise','HQ High noise','LQ High noise');
xlabel('Time of survey (years)'); ylabel('VOI');
set(findall(figure(1),'-property','FontSize'),'FontSize',16)
set(findall(figure(1),'-property','FontName'),'FontName','Cambria')
set(findall(figure(1),'-property','LineWidth'),'LineWidth',2)

figure(2); plot(survey_time(tt2),VOI_hq_low(tt2),'b'); hold on;
plot(survey_time(tt2),VOI_hq_high(tt2),'g');
plot(survey_time(tt2),VOI_lq_high(tt2),'r');
legend('Low detection threshold, Low noise variance','Low detection threshold, High noise variance','High detection threshold, High noise variance','Location','Northwest');
xlabel('Time of survey (years)'); ylabel('VOI');
set(findall(figure(2),'-property','FontSize'),'FontSize',16)
set(findall(figure(2),'-property','FontName'),'FontName','Cambria')
set(findall(figure(2),'-property','LineWidth'),'LineWidth',2)

[f_T1,xi_T1] = ksdensity(T(:,2));
[f_T2,xi_T2] = ksdensity(T(:,3));
[f_T3,xi_T3] = ksdensity(T(:,4));
[f_T4,xi_T4] = ksdensity(T(:,5));
[f_T5,xi_T5] = ksdensity(T(:,6));
[f_T6,xi_T6] = ksdensity(T(:,7));
[f_T7,xi_T7] = ksdensity(T(:,8));
[f_T8,xi_T8] = ksdensity(T(:,9));
[f_T9,xi_T9] = ksdensity(T(:,10));
figure(3);
plot(xi_T1,f_T1); hold on; plot(xi_T2,f_T2); plot(xi_T3,f_T3); ...
    plot(xi_T4,f_T4); plot(xi_T5,f_T5); plot(xi_T6,f_T6); ...
    plot(xi_T7,f_T7); plot(xi_T8,f_T8); plot(xi_T9,f_T9);
xlim([0,200]);
legend('Migration from layer 1 to 2', 'Migration from layer 2 to 3', 'Migration from layer 3 to 4', ...
    'Migration from layer 4 to 5', 'Migration from layer 5 to 6', 'Migration from layer 6 to 7', ...
    'Migration from layer 7 to 8', 'Migration from layer 8 to 9', 'Leakage from layer 9');
xlabel('Time (years)'); ylabel('Probability density');
set(findall(figure(3),'-property','FontSize'),'FontSize',16)
set(findall(figure(3),'-property','FontName'),'FontName','Cambria')
set(findall(figure(3),'-property','LineWidth'),'LineWidth',2)
