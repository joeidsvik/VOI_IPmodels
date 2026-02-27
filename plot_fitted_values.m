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

survey_time = 35;
v_cont = nan(n,numel(survey_time));
v_stop = nan(1,numel(survey_time));
v_cont_fit = nan(n,numel(survey_time));
inst_h = cell(1,numel(survey_time)); %instantaneous heights
imperfect_data = cell(1,numel(survey_time));
coeffs = cell(1,numel(survey_time));

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
    
    imperfect_data{j} = max(0, low_noise.*(inst_h{j} < thres) + (inst_h{j} + low_noise) .* (inst_h{j} >= thres));

    [v_cont(:,j), v_stop(j)] = values(T(:,end), 2, 5, 3, r_inj, survey_time(j), T_plan);
    %{
    plot(v_cont)
    hold on
    plot(v_stop*ones(size(v_cont)))
    %}
    
    %Linear regression
    coeffs{j} = regress(v_cont(:,j),[ones(n,1),imperfect_data{j}]);
    v_cont_fit(:,j) = [ones(n,1),imperfect_data{j}]*coeffs{j};
    
end

figure(1); plot(imperfect_data{1}(:,1),v_cont_fit,'.'); hold on;
plot(imperfect_data{1}(:,1),v_stop*ones(size(v_cont_fit)),'r');
xlabel('Data on heights in layer 1 (m)'); ylabel('Value ($)');
set(findall(figure(1),'-property','FontSize'),'FontSize',40)
legend('Fitted values for continuing','Value for stopping','Location','Best','FontSize',30);
set(findall(figure(1),'-property','FontName'),'FontName','Cambria')
set(findall(figure(1),'-property','LineWidth'),'LineWidth',2)

figure(2); plot(imperfect_data{1}(:,2),v_cont_fit,'.'); hold on;
plot(imperfect_data{1}(:,2),v_stop*ones(size(v_cont_fit)),'r');
xlabel('Data on heights in layer 2 (m)'); ylabel('Value ($)');
set(findall(figure(2),'-property','FontSize'),'FontSize',40)
legend('Fitted values for continuing','Value for stopping','Location','Best','FontSize',30);
set(findall(figure(2),'-property','FontName'),'FontName','Cambria')
set(findall(figure(2),'-property','LineWidth'),'LineWidth',2)

figure(3); plot(imperfect_data{1}(:,3),v_cont_fit,'.'); hold on;
plot(imperfect_data{1}(:,3),v_stop*ones(size(v_cont_fit)),'r');
xlabel('Data on heights in layer 3 (m)'); ylabel('Value ($)');
set(findall(figure(3),'-property','FontSize'),'FontSize',40)
legend('Fitted values for continuing','Value for stopping','Location','Best','FontSize',30);
set(findall(figure(3),'-property','FontName'),'FontName','Cambria')
set(findall(figure(3),'-property','LineWidth'),'LineWidth',2)

figure(4); plot(imperfect_data{1}(:,4),v_cont_fit,'.'); hold on;
plot(imperfect_data{1}(:,4),v_stop*ones(size(v_cont_fit)),'r');
xlabel('Data on heights in layer 4 (m)'); ylabel('Value ($)');
set(findall(figure(4),'-property','FontSize'),'FontSize',40)
legend('Fitted values for continuing','Value for stopping','Location','Northeast','FontSize',30);
set(findall(figure(4),'-property','FontName'),'FontName','Cambria')
set(findall(figure(4),'-property','LineWidth'),'LineWidth',2)

