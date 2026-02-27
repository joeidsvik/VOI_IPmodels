function F = zhat_layer6(a,data)
    x = data(:,1);
    y = data(:,2);
    F = 197.1512 - a(3) + a(3)*sqrt(1-((x-1525).^2)./a(1)^2 - ((y-2375).^2)./a(2)^2);
end

