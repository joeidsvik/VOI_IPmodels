function F = zhat_layer9(a,data)
    x = data(:,1);
    y = data(:,2);
    F = 270.0783 - a(3) + a(3)*sqrt(1-((x-775).^2)./a(1)^2 - ((y-1475).^2)./a(2)^2);
end

