function F = zhat_layer4(a,data)
    x = data(:,1);
    y = data(:,2);
    F = 163.5615 - a(3) + a(3)*sqrt(1-((x-825).^2)./a(1)^2 - ((y-1425).^2)./a(2)^2);
end

