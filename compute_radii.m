function radii = compute_radii(n)
    %n: number of samples of radii to compute for each stratigraphic layer.
    %radii: 1xm cell array, where m is the number of layers.
    % Each cell of radii is an nx3 matrix, the columns of which contain
    % samples of radii a, b and c for that layer.
    %% Load data
    zcorn = load('zcorn_expanded');
    nx = 64; ny = 118; nz = 263;

    L9_top = zcorn((8*nx*ny*10 + 1) : (8*nx*ny*10 + 4*64*118));
    L8_top = zcorn((8*nx*ny*47 + 1) : (8*nx*ny*47 + 4*64*118));
    L7_top = zcorn((8*nx*ny*66 + 1) : (8*nx*ny*66 + 4*64*118));
    L6_top = zcorn((8*nx*ny*85 + 1) : (8*nx*ny*85 + 4*64*118));
    L5_top = zcorn((8*nx*ny*105 + 1) : (8*nx*ny*105 + 4*64*118));
    L4_top = zcorn((8*nx*ny*123 + 1) : (8*nx*ny*123 + 4*64*118));
    L3_top = zcorn((8*nx*ny*149 + 1) : (8*nx*ny*149 + 4*64*118));
    L2_top = zcorn((8*nx*ny*167 + 1) : (8*nx*ny*167 + 4*64*118));
    L1_top = zcorn((8*nx*ny*205 + 1) : (8*nx*ny*205 + 4*64*118));

    L9_top = reshape(L9_top, 2*nx, 2*ny);
    L8_top = reshape(L8_top, 2*nx, 2*ny);
    L7_top = reshape(L7_top, 2*nx, 2*ny);
    L6_top = reshape(L6_top, 2*nx, 2*ny);
    L5_top = reshape(L5_top, 2*nx, 2*ny);
    L4_top = reshape(L4_top, 2*nx, 2*ny);
    L3_top = reshape(L3_top, 2*nx, 2*ny);
    L2_top = reshape(L2_top, 2*nx, 2*ny);
    L1_top = reshape(L1_top, 2*nx, 2*ny);

    fun = @(block_struct) mean(block_struct.data(:));
    L9_top_means = blockproc(L9_top, [2 2], fun);
    L8_top_means = blockproc(L8_top, [2 2], fun);
    L7_top_means = blockproc(L7_top, [2 2], fun);
    L6_top_means = blockproc(L6_top, [2 2], fun);
    L5_top_means = blockproc(L5_top, [2 2], fun);
    L4_top_means = blockproc(L4_top, [2 2], fun);
    L3_top_means = blockproc(L3_top, [2 2], fun);
    L2_top_means = blockproc(L2_top, [2 2], fun);
    L1_top_means = blockproc(L1_top, [2 2], fun);

    dx = 50; dy = 50;
    x1 = dx/2; y1 = dy/2;
    xx = x1:dx:(x1+(nx-1)*dx);
    yy = y1:dy:(y1+(ny-1)*dy);
    [X,Y] = meshgrid(yy,xx);
    Z = {L1_top_means, L2_top_means, L3_top_means, L4_top_means, L5_top_means, L6_top_means, L7_top_means, L8_top_means, L9_top_means};
    for Zi = 1:numel(Z)
        Z{Zi} = -Z{Zi};
    end
    global_min_z = min(Z{1}(:));
    for Zi = 1:numel(Z)
        Z{Zi} = Z{Zi}-global_min_z;
    end
    radii = cell(1,numel(Z));
    opts = optimset('Display','off');

    %% Fit ellipsoids
    % Layer 1
    for count = 1:n
    j = 1;
    
    ind_x_min = randi([17,21],1);
    ind_x_max = randi([36,39],1);
    ind_y_min = randi([76,81],1);
    ind_y_max = randi([114,117],1);
    ind_x = ind_x_min : ind_x_max;
    ind_y = ind_y_min : ind_y_max;
    
    x = X(ind_x, ind_y);
    y = Y(ind_x, ind_y);
    z = Z{j}(ind_x, ind_y);

    radius = [1500, 2500, 1000];

    x=x(:);y=y(:);z=z(:);
    data=[x,y];
    [radii{j}(count,:), resnorm] = lsqcurvefit(@zhat_layer1, radius, data, z, [], [], opts);
    
    end
    
    %%
    % Layer 2
    for count = 1:n
    j = 2;
    
    ind_x_min = randi([18,21],1);
    ind_x_max = randi([33,36],1);
    ind_y_min = randi([8,11],1);
    ind_y_max = randi([22,25],1);
    ind_x = ind_x_min : ind_x_max;
    ind_y = ind_y_min : ind_y_max;
    x = X(ind_x, ind_y);
    y = Y(ind_x, ind_y);
    z = Z{j}(ind_x, ind_y);

    radius = [1500, 2500, 1000];

    x=x(:);y=y(:);z=z(:);
    data=[x,y];
    [radii{j}(count,:), resnorm] = lsqcurvefit(@zhat_layer2, radius, data, z, [], [], opts);
    
    end
    %%
    % Layer 3
    for count = 1:n
    j = 3;
    
    ind_x_min = randi([17,21],1);
    ind_x_max = randi([36,39],1);
    ind_y_min = randi([76,81],1);
    ind_y_max = randi([114,117],1);
    ind_x = ind_x_min : ind_x_max;
    ind_y = ind_y_min : ind_y_max;
    
    x = X(ind_x, ind_y);
    y = Y(ind_x, ind_y);
    z = Z{j}(ind_x, ind_y);

    radius = [1500, 2500, 1000];

    x=x(:);y=y(:);z=z(:);
    data=[x,y];
    [radii{j}(count,:), resnorm] = lsqcurvefit(@zhat_layer3, radius, data, z, [], [], opts);
    
    end
    %%
    % Layer 4
    for count = 1:n
    j = 4;
    
    ind_x_min = randi([18,21],1);
    ind_x_max = randi([33,36],1);
    ind_y_min = randi([8,11],1);
    ind_y_max = randi([22,25],1);
    ind_x = ind_x_min : ind_x_max;
    ind_y = ind_y_min : ind_y_max;
    
    x = X(ind_x, ind_y);
    y = Y(ind_x, ind_y);
    z = Z{j}(ind_x, ind_y);

    radius = [1500, 2500, 1000];

    x=x(:);y=y(:);z=z(:);
    data=[x,y];
    [radii{j}(count,:), resnorm] = lsqcurvefit(@zhat_layer4, radius, data, z, [], [], opts);
    
    end
    %%
    % Layer 5
    for count = 1:n
    j = 5;
    
    ind_x_min = randi([35,37],1);
    ind_x_max = randi([47,50],1);
    ind_y_min = randi([16,20],1);
    ind_y_max = randi([42,44],1);
    ind_x = ind_x_min : ind_x_max;
    ind_y = ind_y_min : ind_y_max;
    
    x = X(ind_x, ind_y);
    y = Y(ind_x, ind_y);
    z = Z{j}(ind_x, ind_y);
    
    radius = [1500, 2500, 1000];

    x=x(:);y=y(:);z=z(:);
    data=[x,y];
    [radii{j}(count,:), resnorm] = lsqcurvefit(@zhat_layer5, radius, data, z, [], [], opts);
    
    end
    %%
    % Layer 6
    for count = 1:n
    j = 6;
    
    ind_x_min = randi([38,40],1);
    ind_x_max = randi([51,53],1);
    ind_y_min = randi([19,22],1);
    ind_y_max = randi([45,50],1);
    ind_x = ind_x_min : ind_x_max;
    ind_y = ind_y_min : ind_y_max;
    
    x = X(ind_x, ind_y);
    y = Y(ind_x, ind_y);
    z = Z{j}(ind_x, ind_y);
    
    radius = [1500, 2500, 1000];

    x=x(:);y=y(:);z=z(:);
    data=[x,y];
    [radii{j}(count,:), resnorm] = lsqcurvefit(@zhat_layer6, radius, data, z, [], [], opts);
    
    end
    %%
    % Layer 7
    for count = 1:n
    j = 7;
    
    ind_x_min = randi([21,24],1);
    ind_x_max = randi([35,38],1);
    ind_y_min = randi([8,11],1);
    ind_y_max = randi([22,26],1);
    ind_x = ind_x_min : ind_x_max;
    ind_y = ind_y_min : ind_y_max;
    
    x = X(ind_x, ind_y);
    y = Y(ind_x, ind_y);
    z = Z{j}(ind_x, ind_y);
    
    radius = [1500, 2500, 1000];

    x=x(:);y=y(:);z=z(:);
    data=[x,y];
    [radii{j}(count,:), resnorm] = lsqcurvefit(@zhat_layer7, radius, data, z, [], [], opts);
    
    end
    %%
    % Layer 8
    for count = 1:n
    j = 8;
    
    ind_x_min = randi([21,24],1);
    ind_x_max = randi([35,38],1);
    ind_y_min = randi([8,11],1);
    ind_y_max = randi([22,26],1);
    ind_x = ind_x_min : ind_x_max;
    ind_y = ind_y_min : ind_y_max;
    
    x = X(ind_x, ind_y);
    y = Y(ind_x, ind_y);
    z = Z{j}(ind_x, ind_y);
    
    radius = [1500, 2500, 1000];

    x=x(:);y=y(:);z=z(:);
    data=[x,y];
    [radii{j}(count,:), resnorm] = lsqcurvefit(@zhat_layer8, radius, data, z, [], [], opts);
    
    end
    %%
    % Layer 9
    for count = 1:n
    j = 9;
    
    ind_x_min = randi([23,26],1);
    ind_x_max = randi([37,39],1);
    ind_y_min = randi([6,10],1);
    ind_y_max = randi([23,26],1);
    ind_x = ind_x_min : ind_x_max;
    ind_y = ind_y_min : ind_y_max;
    
    x = X(ind_x, ind_y);
    y = Y(ind_x, ind_y);
    z = Z{j}(ind_x, ind_y);
    
    radius = [1500, 2500, 1000];

    x=x(:);y=y(:);z=z(:);
    data=[x,y];
    [radii{j}(count,:), resnorm] = lsqcurvefit(@zhat_layer9, radius, data, z, [], [], opts);
    
    end
end
