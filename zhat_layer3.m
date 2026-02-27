function F = zhat_layer3(a,data)
    x = data(:,1);
    y = data(:,2);
    F = 139.3115 - a(3) + a(3)*sqrt(1-((x-4775).^2)./a(1)^2 - ((y-1525).^2)./a(2)^2);
end

