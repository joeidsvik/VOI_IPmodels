function [T,C_th,h,a,b,c] = event_times(r, drho, eta, tau, mlna, slna, mlnb, slnb, mlnc, slnc, n)
    % r: volume rate of CO2 injection
    % drho: density difference between brine and CO2
    % eta: vector of means of ln C_th for different layers
    % tau: vector of standard deviations of ln C_th for different layers
    % mlna: vector of means of ln a for different layers
    % slna: vector of standard deviations of ln a for different layers
    % mlnb: vector of means of ln b for different layers
    % slnb: vector of standard deviations of ln b for different layers
    % mlnc: vector of means of ln c for different layers
    % slnc: vector of standard deviations of ln c for different layers
    % n: number of samples to simulate
    % T: vector of event times
    % C_th: vector of capillary pressure thresholds
    g = 9.81; %acceleration due to gravity
    m = length(eta); % number of layers
    T = zeros(n,m+1);
    C_th = zeros(n,m);
    a = zeros(n,m);
    b = zeros(n,m);
    c = zeros(n,m);
    h = zeros(n,m);
    for i = 1:m
        ln_C_th = eta(i)+tau(i).*randn(n,1);
        ln_a = mlna(i)+slna(i).*randn(n,1);
        ln_b = mlnb(i)+slnb(i).*randn(n,1);
        ln_c = mlnc(i)+slnc(i).*randn(n,1);
        C_th(:,i) = exp(ln_C_th);
        a(:,i) = exp(ln_a);
        b(:,i) = exp(ln_b);
        c(:,i) = exp(ln_c);
        h(:,i) = C_th(:,i)/(drho*g);
        T(:,i+1) = ((pi.*a(:,i).*b(:,i))./(r.*3.*c(:,i).^2)).*(h(:,i).^2).*(3.*c(:,i)-h(:,i)) + T(:,i);
    end
end