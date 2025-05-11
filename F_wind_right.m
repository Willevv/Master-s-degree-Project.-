function func = F_wind_right(x, R)
    fixed_bin_width = 0.01;   % Bin width for the combined histogram
    numpoints = 10000;       % Number of points for curve discretization
   
    s1 = x(1:21);  % Stagger for the first posts
    s2 = x(22:42); % Stagger for the second posts
    sl = x(43:end);
%     sl = 65 * ones(1, 21);
    Counts = cell(length(R), 2);
    
    for i = 1:length(R)
        r = R(i);
        alpha = sl(i)/r;
        delta_x = (r + s1(i)) * cos(0) - (r + s2(i)) * cos(alpha);
        delta_y = (r + s2(i)) * sin(alpha);
        k = -delta_y / delta_x;
        C = k * (r + s1(i));
        theta = linspace(0, alpha, numpoints);
        D = -C ./ (sin(theta) - k * cos(theta));
        cord_length = sqrt(delta_x^2 + delta_y^2); % Length of the Cord connecting the Posts 
        x_e = linspace(0, cord_length, numpoints);
        e = [];
        %--------Eq.10--LAROBOK-KONTAKTLEDNING--------------------
        S = cord_length;
        p = 6.5;
        F = 15000;
        %------------------------------------------------------
        D_w = zeros(1,length(theta));
        x_w = zeros(1,length(theta));
        y_w = zeros(1,length(theta));
        % e_upd = (s2-s1).*x_e/S + s1 + x_e.*(x_e-S)*(1/R + p/F)/2;
        e_upd = -x_e.*(x_e-S)*(p/F)/2;      
        D_w_pos = -C./(sin(theta) - k.*cos(theta))-e_upd - r; %Just add Vindavdrift to this 
        dist = D_w_pos;
        min(dist);
        max_dist = max(dist);
        max_dist=round(max_dist, 2); %rounding to closest 0.01 decimal-----------------------------------------
        min_dist = min(dist);
        min_dist = round(min_dist,2); %rounding to closest 0.01 decimal-----------------------------------------
        % Histogram Calculation
        binEdges = min_dist:fixed_bin_width:max_dist;
        [counts, edges] = histcounts(dist, 'BinEdges', binEdges);
        binCenters = edges(1:end-1) + diff(edges) / 2;
        Counts{i, 1} = counts;
        Counts{i, 2} = binCenters;
    end
    
    % Define combined bin edges that fully span [min(dist), max(dist)]
    global_min = min(cellfun(@(c) min(c), Counts(:, 2))) - fixed_bin_width / 2;
    global_max = max(cellfun(@(c) max(c), Counts(:, 2))) + fixed_bin_width / 2;

    common_bin_edges = global_min:fixed_bin_width:global_max;
    combined_counts = zeros(1, length(common_bin_edges) - 1);
    

    % --------------------------Combine histograms----------------------------
    for i = 1:length(Counts)
        counts = Counts{i, 1};
        bin_centers = Counts{i, 2};
       
        % Assign bin indices with the correct treatment of edges
        [~, bin_indices] = histc(bin_centers, [-inf, common_bin_edges(2:end-1), inf]); % Extend edges to include all points

        for j = 1:length(counts)
            if bin_indices(j) > 0 && bin_indices(j) <= length(combined_counts)
                combined_counts(bin_indices(j)) = combined_counts(bin_indices(j)) + counts(j);
            end
        end
    end
    
    
    %- - - THE- BELOW- SECTION- CREATES- MIRRORED- DISTRIBUTION- - -- - - -
    smoothing_window = 7; 
    smoothed_counts = smoothdata(combined_counts, 'gaussian', smoothing_window); %Not normalized
    mirrored = smoothed_counts(end:-1:1);
    if abs(global_max) > abs(global_min)
        global_min = -abs(global_max);
        common_bin_edges = global_min:fixed_bin_width:global_max;
        ORIGINAL = zeros(1, length(common_bin_edges)-1);
        MIRRORED = zeros(1, length(common_bin_edges)-1);
        ORIGINAL(length(ORIGINAL)-length(smoothed_counts)+1:end) = smoothed_counts;
        MIRRORED(1:length(smoothed_counts)) = mirrored;
        
    elseif abs(global_min) > abs(global_max)
        global_max = abs(global_min);
        common_bin_edges = global_min:fixed_bin_width:global_max;
        ORIGINAL = zeros(1, length(common_bin_edges)-1);
        MIRRORED = zeros(1, length(common_bin_edges)-1);
        ORIGINAL(1:length(smoothed_counts)) = smoothed_counts;
        MIRRORED(length(ORIGINAL)-length(smoothed_counts)+1:end) = mirrored;
    else
        ORIGINAL = zeros(1, length(common_bin_edges)-1);
        MIRRORED = zeros(1, length(common_bin_edges)-1);
        ORIGINAL(1:length(smoothed_counts)) = smoothed_counts;
        MIRRORED(1:length(smoothed_counts)) = mirrored;
    end
    
    FINAL = ORIGINAL + MIRRORED; % NOT NORMALIZED
    FINAL = FINAL/(sum(FINAL)*fixed_bin_width);

    sum(FINAL*fixed_bin_width);
    index_smaller = find(common_bin_edges <- 0.4);
    index_bigger = find(common_bin_edges > 0.4) -1;
    func = sum(FINAL(index_smaller)*fixed_bin_width) + sum(FINAL(index_bigger)*fixed_bin_width);
    
end