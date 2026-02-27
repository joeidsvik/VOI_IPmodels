close all; clear all; clc;
rng('default');
r_inj = 1e7;
T_plan = 40; % total planned time period of injection 
drho = 400;
eta = [log(7.6e5),log(7.4e5),log(7.2e5),log(7.9e5)];
tau = [0.1,0.08,0.07,0.12];
mlna = [log(1000),log(950),log(900),log(970)];
slna = [0.05,0.05,0.05,0.05];
mlnb = [log(1000),log(950),log(900),log(970)];
slnb = [0.05,0.05,0.05,0.05];
mlnc = [log(1000),log(950),log(900),log(970)];
slnc = [0.05,0.05,0.05,0.05];
n = 10000;
m = length(eta); %number of stratigraphic layers
[T,C_th,h,a,b,c] = event_times(r_inj, drho, eta, tau, mlna, slna, mlnb, slnb, mlnc, slnc, n);
sn_high = 50; %standard deviation of high noise data
sn_low = 25; %standard deviation of low noise data
thres = 5; %minimum height that can be detected by data

survey_time = 1:60;
v_cont = nan(n,numel(survey_time));
v_stop = nan(1,numel(survey_time));
v_cont_fit_perfect = nan(n,numel(survey_time));
v_cont_fit_high = nan(n,numel(survey_time));
v_cont_fit_low = nan(n,numel(survey_time));
VOI_perfect = nan(1,numel(survey_time));
VOI_imperfect_high = nan(1,numel(survey_time));
VOI_imperfect_low = nan(1,numel(survey_time));
inst_h = cell(1,numel(survey_time)); %instantaneous heights
perfect_data = cell(1,numel(survey_time));
imperfect_data_high_noise = cell(1,numel(survey_time));
imperfect_data_low_noise = cell(1,numel(survey_time));
coeffs_perfect = cell(1,numel(survey_time));
coeffs_high = cell(1,numel(survey_time));
coeffs_low = cell(1,numel(survey_time));

high_noise = sn_high*randn(n,m);
low_noise = sn_low*randn(n,m);

for j = 1:numel(survey_time)
    inst_h{j} = h.*(survey_time(j) > T(:,2:end));
    for i=1:n
        if any(T(i,:)>survey_time(j))
            ind = find(T(i,:)>survey_time(j), 1);
            V_t = r_inj * (survey_time(j) - T(i,ind-1));
            inst_h{j}(i,ind-1) = compute_h(V_t,a(i,ind-1),b(i,ind-1),c(i,ind-1),0,c(i,ind-1));
        end
    end
    
    perfect_data{j} = max(0, 0.*(inst_h{j} < thres) + inst_h{j} .* (inst_h{j} >= thres));
    imperfect_data_high_noise{j} = max(0, high_noise.*(inst_h{j} < thres) + (inst_h{j} + high_noise) .* (inst_h{j} >= thres));
    imperfect_data_low_noise{j} = max(0, low_noise.*(inst_h{j} < thres) + (inst_h{j} + low_noise) .* (inst_h{j} >= thres));

    [v_cont(:,j), v_stop(j)] = values(T(:,end), 2, 5, 3, r_inj, survey_time(j), T_plan);
    %{
    plot(v_cont)
    hold on
    plot(v_stop*ones(size(v_cont)))
    %}
    
    %Linear regression
    coeffs_perfect{j} = regress(v_cont(:,j),[ones(n,1),perfect_data{j}]);
    coeffs_high{j} = regress(v_cont(:,j),[ones(n,1),imperfect_data_high_noise{j}]);
    coeffs_low{j} = regress(v_cont(:,j),[ones(n,1),imperfect_data_low_noise{j}]);
    v_cont_fit_perfect(:,j) = [ones(n,1),perfect_data{j}]*coeffs_perfect{j};
    v_cont_fit_high(:,j) = [ones(n,1),imperfect_data_high_noise{j}]*coeffs_high{j};
    v_cont_fit_low(:,j) = [ones(n,1),imperfect_data_low_noise{j}]*coeffs_low{j};
    
    PV = max(mean(v_cont(:,j)),v_stop(j));
    %PV_perfect = max(mean(v_cont_fit_perfect(:,j)),v_stop(j));
    %PV_imperfect_high = max(mean(v_cont_fit_high(:,j)),v_stop(j));
    %PV_imperfect_low = max(mean(v_cont_fit_low(:,j)),v_stop(j));
    %PoV_perfect = mean(max(v_cont(:,j),v_stop(j)));
    PoV_perfect = mean(max(v_cont_fit_perfect(:,j),v_stop(j)));
    PoV_imperfect_high = mean(max(v_cont_fit_high(:,j),v_stop(j)));
    PoV_imperfect_low = mean(max(v_cont_fit_low(:,j),v_stop(j)));
    VOI_perfect(j) = PoV_perfect - PV;
    VOI_imperfect_high(j) = PoV_imperfect_high - PV;
    VOI_imperfect_low(j) = PoV_imperfect_low - PV;
    
    fprintf('Completed computation for survey time = %f years\n',survey_time(j));
end

figure(1); plot(survey_time,VOI_perfect,'r'); hold on;
plot(survey_time,VOI_imperfect_low,'g');
plot(survey_time,VOI_imperfect_high,'b');
yl = ylim; ylim([0,yl(2)]);
legend('Perfect data','Imperfect data (low noise)','Imperfect data (high noise)','Location','best');
xlabel('Time of survey (years)'); ylabel('VOI ($)');
%title('VOI for 4 layer case');
set(findall(figure(1),'-property','FontSize'),'FontSize',24)
set(findall(figure(1),'-property','FontName'),'FontName','Cambria')
set(findall(figure(1),'-property','LineWidth'),'LineWidth',2)

[f_T1,xi_T1] = ksdensity(T(:,2));
[f_T2,xi_T2] = ksdensity(T(:,3));
[f_T3,xi_T3] = ksdensity(T(:,4));
[f_T4,xi_T4] = ksdensity(T(:,5));
figure(2);
plot(xi_T1,f_T1,'g'); hold on; plot(xi_T2,f_T2,'b'); plot(xi_T3,f_T3,'r'); plot(xi_T4,f_T4,'m');
legend('Migration from layer 1 to 2', 'Migration from layer 2 to 3', 'Migration from layer 3 to 4', 'Leakage from layer 4');
xlabel('Time (years)'); ylabel('Probability density');
set(findall(figure(2),'-property','FontSize'),'FontSize',24)
set(findall(figure(2),'-property','FontName'),'FontName','Cambria')
set(findall(figure(2),'-property','LineWidth'),'LineWidth',2)

h1 = cell2mat(cellfun(@(x) x(:,1), inst_h, 'UniformOutput', false));
h2 = cell2mat(cellfun(@(x) x(:,2), inst_h, 'UniformOutput', false));
h3 = cell2mat(cellfun(@(x) x(:,3), inst_h, 'UniformOutput', false));
h4 = cell2mat(cellfun(@(x) x(:,4), inst_h, 'UniformOutput', false));
ti = 1:2:numel(survey_time);

figure(3);
boxplot(h1(:,ti));
xticklabels(survey_time(ti))
ax = gca;
labels = string(ax.XAxis.TickLabels);
labels(2:2:end) = ' ';
ax.XAxis.TickLabels = labels;
xlabel('Time (years)'); ylabel('Height (m)');
%title('Column height in first layer with time');
set(findall(figure(3),'-property','FontSize'),'FontSize',40)
set(findall(figure(3),'-property','FontName'),'FontName','Cambria')
set(findall(figure(3),'-property','LineWidth'),'LineWidth',2)

figure(4);
boxplot(h2(:,ti));
xticklabels(survey_time(ti))
ax = gca;
labels = string(ax.XAxis.TickLabels);
labels(2:2:end) = ' ';
ax.XAxis.TickLabels = labels;
xlabel('Time (years)'); ylabel('Height (m)');
%title('Column height in second layer with time');
set(findall(figure(4),'-property','FontSize'),'FontSize',40)
set(findall(figure(4),'-property','FontName'),'FontName','Cambria')
set(findall(figure(4),'-property','LineWidth'),'LineWidth',2)

figure(5);
boxplot(h3(:,ti));
xticklabels(survey_time(ti))
ax = gca;
labels = string(ax.XAxis.TickLabels);
labels(2:2:end) = ' ';
ax.XAxis.TickLabels = labels;
xlabel('Time (years)'); ylabel('Height (m)');
%title('Column height in third layer with time');
set(findall(figure(5),'-property','FontSize'),'FontSize',40)
set(findall(figure(5),'-property','FontName'),'FontName','Cambria')
set(findall(figure(5),'-property','LineWidth'),'LineWidth',2)

figure(6);
boxplot(h4(:,ti));
xticklabels(survey_time(ti))
ax = gca;
labels = string(ax.XAxis.TickLabels);
labels(2:2:end) = ' ';
ax.XAxis.TickLabels = labels;
xlabel('Time (years)'); ylabel('Height (m)');
%title('Column height in fourth layer with time');
set(findall(figure(6),'-property','FontSize'),'FontSize',40)
set(findall(figure(6),'-property','FontName'),'FontName','Cambria')
set(findall(figure(6),'-property','LineWidth'),'LineWidth',2)


