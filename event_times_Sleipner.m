function [T,C_th,h,a,b,c] = event_times_Sleipner (r, drho, eta, tau, n, radii)
    % r: volume rate of CO2 injection
    % drho: density difference between brine and CO2
    % eta: vector of means of ln C_th for different layers
    % tau: vector of standard deviations of ln C_th for different layers
    % n: number of samples to simulate
    % T: vector of event times
    % C_th: vector of capillary pressure thresholds
    % h: maximum heights at which migration happens
    % a, b, c: radii of fitted ellipsoids
    % radii: cell array containing samples of radii for each layer
    g = 9.81; %acceleration due to gravity
    m = length(eta); % number of layers
    T = zeros(n,m+1);
    C_th = zeros(n,m);
    h = zeros(n,m);
    a = zeros(n,m);
    b = zeros(n,m);
    c = zeros(n,m);
    
    for i = 1:m
        a(:,i) = radii{i}(:,1);
        b(:,i) = radii{i}(:,2);
        c(:,i) = radii{i}(:,3);
        ln_C_th = eta(i)+tau(i).*randn(n,1);
        C_th(:,i) = exp(ln_C_th);
        h(:,i) = C_th(:,i)/(drho*g);
        T(:,i+1) = ((pi.*a(:,i).*b(:,i))./(r.*3.*c(:,i).^2)).*(h(:,i).^2).*(3.*c(:,i)-h(:,i)) + T(:,i);
    end
end