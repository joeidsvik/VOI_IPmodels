function h = compute_h(V,a,b,c,min_h,max_h)
%This function computes the height of an ellipsoidal cap 
%   h: height of the ellipsoidal cap
%   V: volume of the ellipsoidal cap
%   a, b, c: axes of the ellipsoid
%   min_h, max_h: minimum and maximum possible values of h
    c1 = (pi*a*b)/(3*c^2);
    c2 = (-pi*a*b)/c;
    c3 = 0;
    c4 = V;
    h = roots([c1,c2,c3,c4]);
    h = h(h>=min_h & h<=max_h);
end

