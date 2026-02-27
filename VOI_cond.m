close all; clear all; clc;

% Set parameter values
sn = 2; %standard deviation of noise (high = 4, low = 2)
thres = 1; %minimum height that can be detected by data (high = 10, low = 1)
true_data_case = 2; %1 = max, 2 = p percentile (medium data), 3 = min
monitor_time = 5;
p = 75; % percentile of medium data

rng('default');
r_inj = 1e7;
T_total = 40; % Total planned time period of injection
drho = 400;
eta = log(1e5)*ones(1,9);
tau = 0.5*ones(1,9);
n = 5000;
load radii1.mat;
[T,C_th,h,a,b,c] = event_times_Sleipner(r_inj, drho, eta, tau, n, radii1);

g = 9.81; %acceleration due to gravity
m = length(eta); %number of stratigraphic layers

survey_time = 1:1:20;

v_cont = nan(n-1,numel(survey_time));
v_stop = nan(1,numel(survey_time));

v_cont_fit= nan(n-1,numel(survey_time));

VOI = nan(1,numel(survey_time));

inst_h = cell(1,numel(survey_time)); %instantaneous heights
data = cell(1,numel(survey_time));
true_data = cell(1,numel(survey_time));
dataset = cell(1,numel(survey_time)); % dataset = {data, true_data}
coeffs = cell(1,numel(survey_time));

noise = sn*randn(n,m);

inst_h_mon_time = h.*(monitor_time > T(:,2:end));
for i = 1:n
    if any(T(i,:) > monitor_time)
        ind = find(T(i,:) > monitor_time, 1);
        V_t = r_inj * (monitor_time - T(i,ind-1));
        inst_h_mon_time(i,ind-1) = compute_h(V_t,a(i,ind-1),b(i,ind-1),c(i,ind-1),0,c(i,ind-1));
    end
end

if true_data_case == 1
    [true_val, true_ind] = max(sum(inst_h_mon_time(:,2:end),2));
elseif true_data_case == 2
    inst_h_mon_time_sort = sort(sum(inst_h_mon_time(:,2:end),2));
    p_prctile_ind = round((p/100)*length(inst_h_mon_time_sort));
    true_val = inst_h_mon_time_sort(p_prctile_ind);
    true_ind = find(sum(inst_h_mon_time(:,2:end),2) == true_val, 1);
elseif true_data_case == 3
    [true_val, true_ind] = min(sum(inst_h_mon_time(:,2:end),2));
else
    error('Invalid true data case');
end

data_ind = setdiff(1:size(inst_h_mon_time,1), true_ind);

for j = 1:numel(survey_time)
    inst_h{j} = h.*(survey_time(j) > T(:,2:end));
    for i=1:n
        if any(T(i,:)>survey_time(j))
            ind = find(T(i,:)>survey_time(j), 1);
            V_t = r_inj * (survey_time(j) - T(i,ind-1));
            inst_h{j}(i,ind-1) = compute_h(V_t,a(i,ind-1),b(i,ind-1),c(i,ind-1),0,c(i,ind-1));
        end
    end
    
    dataset{j} = max(0, noise.*(inst_h{j} < thres) + (inst_h{j} + noise) .* (inst_h{j} >= thres));
    
    true_data{j} = dataset{j}(true_ind,:);
    data{j} = dataset{j}(data_ind,:);
    
    [v_cont(:,j), v_stop(j)] = values(T(data_ind,end), 2, 5, 3, r_inj, survey_time(j), T_total);
    
    %Linear regression
    coeffs{j} = regress(v_cont(:,j),[ones(n-1,1),data{j}]);
    v_cont_fit(:,j) = [ones(n-1,1),data{j}]*coeffs{j};
    
    PV = max(mean(v_cont(:,j)),v_stop(j));
    PoV= mean(max(v_cont_fit(:,j),v_stop(j)));
    VOI(j) = PoV - PV;
    
    if survey_time(j) == monitor_time
         % Update model parameters
        x = [a(data_ind,:)';b(data_ind,:)';c(data_ind,:)';C_th(data_ind,:)'];
        d = data{j}';
        td = true_data{j}';
        x_post = update_x(x,d,td);
        a(data_ind,:) = x_post(1:m,:)';
        b(data_ind,:) = x_post(m+1:2*m,:)';
        c(data_ind,:) = x_post(2*m+1:3*m,:)';
        C_th(data_ind,:) = x_post(3*m+1:4*m,:)';

        for i=1:m
            h(data_ind,i) = C_th(data_ind,i)/(drho*g);
            T(data_ind,i+1) = ((pi.*a(data_ind,i).*b(data_ind,i))./(r_inj.*3.*c(data_ind,i).^2))...
                .*(h(data_ind,i).^2).*(3.*c(data_ind,i)-h(data_ind,i)) + T(data_ind,i);
        end
    end
    
    fprintf('Completed computation for survey time = %f years\n',survey_time(j));
end

tt = 1:numel(survey_time);
tt2 = tt;

figure(1); plot(survey_time(tt2),VOI(tt2));
xlabel('Time of survey (years)'); ylabel('VOI');
title(['True data percentile = ', num2str(p)]);
set(findall(figure(1),'-property','FontSize'),'FontSize',16)
set(findall(figure(1),'-property','FontName'),'FontName','Cambria')
set(findall(figure(1),'-property','LineWidth'),'LineWidth',2)
