function HIST_UPDATED(x, R)
    fixed_bin_width = 0.001;   % Bin width for the combined histogram
    numpoints = 100000;       % Number of points for curve discretization

    s1 = x(1:21);  % Stagger for the first posts
    s2 = x(22:42); % Stagger for the second posts
    sl = x(43:end);
%     sl = 65 * ones(1, 21);
    Counts = cell(length(R), 2);
    
     %
     %andel = ones(1,21);
    andel = [0.00568576716340912,0.00254721525960437,0.000762984433473785,0.00229885808384167,0.00208917171143278,0.00383452099413679,0.00626677754404362,0.00359722767219466,0.00831812141241369,0.0215639779464729,0.0262012078778011,0.0261087983558902,0.0131481785103290,0.00434409049010133,0.00796860900172006,0.00396591742951061,0.00782309298147047,0.0143170483239631,0.00427296571554052,0.00819030755827701,0.0282362193905313];
    andel = andel./sum(andel);
    sum(andel);
    SUM = [];

    for i = 1:length(R)
        r = R(i);
        alpha = sl(i)/r;
        delta_x = (r + s1(i)) * cos(0) - (r + s2(i)) * cos(alpha);
        delta_y = (r + s2(i)) * sin(alpha);
        k = -delta_y / delta_x;
        C = k * (r + s1(i));
        theta = linspace(0, alpha, numpoints);
        D = -C ./ (sin(theta) - k * cos(theta));
        dist = D - r;
        min(dist);
        max_dist = max(dist);
        max_dist=round(max_dist, 2); %rounding to closest 0.01 decimal-----------------------------------------
        min_dist = min(dist);
        min_dist = round(min_dist,2); %rounding to closest 0.01 decimal-----------------------------------------
        % Histogram Calculation
        binEdges = min_dist:fixed_bin_width:max_dist;
        [counts, edges] = histcounts(dist, 'BinEdges', binEdges);
        SUM = [SUM sum(counts)];
        binCenters = edges(1:end-1) + diff(edges) / 2;
        Counts{i, 1} = counts;
        Counts{i, 2} = binCenters;
    end
    minSUM = min(SUM);

    % Define combined bin edges that fully span [min(dist), max(dist)]
    global_min = min(cellfun(@(c) min(c), Counts(:, 2))) - fixed_bin_width / 2;
    global_max = max(cellfun(@(c) max(c), Counts(:, 2))) + fixed_bin_width / 2;

    common_bin_edges = global_min:fixed_bin_width:global_max;
    combined_counts = zeros(1, length(common_bin_edges) - 1);
    
    % --------------------------Combine histograms----------------------------
    for i = 1:length(Counts)
        counts = Counts{i, 1};
        counts = counts.*minSUM/sum(counts);

        bin_centers = Counts{i, 2};
       
        % Assign bin indices with the correct treatment of edges
        [~, bin_indices] = histc(bin_centers, [-inf, common_bin_edges(2:end-1), inf]); % Extend edges to include all points

        for j = 1:length(counts)
            if bin_indices(j) > 0 && bin_indices(j) <= length(combined_counts)
                combined_counts(bin_indices(j)) = combined_counts(bin_indices(j)) + andel(i)*counts(j);
            end
        end
    end
    
 
    % - - - - - - - - - - -Spegla- Båda Hållen - - - - - - - - - - - - - - - - - - - - - -- 
    mirrored_counts = combined_counts(end:-1:1);
    %-- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    
    if abs(global_max) > abs(global_min)
        global_min = -abs(global_max);
        common_bin_edges = global_min:fixed_bin_width:global_max;
        fall = 1;
    elseif abs(global_min) > abs(global_max)
        global_max = abs(global_min);
        common_bin_edges = global_min:fixed_bin_width:global_max;
        fall = 2;
    else
        fall = 0;
    end
    ORIGINAL_DATA_PLACEHOLDER = zeros(1, length(common_bin_edges) - 1);
    MIRRORED_DATA_PLACEHOLDER = zeros(1, length(common_bin_edges) - 1);
    
    if fall == 1
        MIRRORED_DATA_PLACEHOLDER(1:length(mirrored_counts)) = mirrored_counts;
        ORIGINAL_DATA_PLACEHOLDER(length(common_bin_edges)-length(combined_counts):end) = combined_counts;
        
    elseif fall == 2
        ORIGINAL_DATA_PLACEHOLDER(1:length(combined_counts)) = combined_counts;
        MIRRORED_DATA_PLACEHOLDER(length(common_bin_edges)-length(mirrored_counts):end) = mirrored_counts;    
      
    elseif fall == 0
        ORIGINAL_DATA_PLACEHOLDER(1:length(combined_counts)) = combined_counts;
        MIRRORED_DATA_PLACEHOLDER(1:length(combined_counts)) = mirrored_counts;
        
    end
    
    TOTAL = ORIGINAL_DATA_PLACEHOLDER + MIRRORED_DATA_PLACEHOLDER;
    
%      % - -  - - - - - - - - - - -
%     TOTAL(1) = TOTAL(1) + 1000;
%     % - - - - - - - - - - - - -
    
    %---------------SMOTHING--------------------
    TOTAL = TOTAL/(sum(TOTAL)*fixed_bin_width);
    smoothing_window = 100; 
    TOTAL = smoothdata(TOTAL, 'gaussian', smoothing_window);
    TOTAL = TOTAL/(sum(TOTAL)*fixed_bin_width);
    sum(TOTAL*fixed_bin_width)
   
    bar(common_bin_edges(1:end-1), TOTAL, 'histc') %--2
    xlim([min(common_bin_edges) - 0.02, max(common_bin_edges)+0.02]);
    ylim([0,3])
%     title('Resulitng combined distribution for all radii')
%     xlabel('Distance from center(m)')
%     ylabel('Normalized combined counts')  

end