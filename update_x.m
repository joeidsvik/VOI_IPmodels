function x_posterior = update_x(x,d,dtrue)

%   x: matrix containing prior samples of model parameters, having realizations in columns and variables in rows

%   d: matrix containing synthetic data, having realizations in columns and variables in rows

%   x_posterior: matrix containing posterior samples of model parameters, having realizations in columns and variables in rows

    xc = x-mean(x,2);

    dc = d-mean(d,2);

    B = size(x,2);

    S_d = (1/B)*(dc*dc');

    S_xd = (1/B)*(xc*dc');

    x_posterior = x+(S_xd/S_d)*(dtrue-d);

end