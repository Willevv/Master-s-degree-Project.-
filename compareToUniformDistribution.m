function similarity= compareToUniformDistribution(counts, x)
    % compareToUniformDistribution compares a normalized distribution to a uniform distribution.
    %
    % INPUTS:
    fixed_bin_width = 0.01;
    %   counts - a vector representing the counts or frequencies of the distribution
    %   x - a vector representing the corresponding x-values for the counts
    %
    % OUTPUT:
    %   similarity - a measure of similarity between the distribution and
    %                the uniform distribution over [-0.3, 0.3].
    
    % Normalize the counts to form a probability distribution
    bin_Edges = x - 0.005;
    bin_Edges = [bin_Edges max(x)+0.005];
    bin_centers = x;
    normalized_counts = counts;
    
    % Define the uniform distribution over [-0.3, 0.3]
    uniform_min = -0.3;
    uniform_max = 0.3;
    global_min = min(uniform_min, bin_Edges(1));
    global_max = max(uniform_max, bin_Edges(end));
    % Create a uniform PDF with the same x-resolution

    uniform_pdf = zeros(1, round((global_max-global_min)/fixed_bin_width));
%     in_range = (x >= uniform_min) & (x <= uniform_max);
%     uniform_pdf(in_range) = 1 / Nbins; % Uniform height within range
    %---------------------------------------------------------------------------
    if max(bin_Edges) < uniform_max
        dist_to_cover = uniform_max - max(bin_Edges);
        number_of_bins_to_add = dist_to_cover/fixed_bin_width;
        %adding zeros to x in order to covver full range
        x_adder = zeros(1,round(number_of_bins_to_add));
        normalized_counts = [normalized_counts x_adder];
        bin_centers = [bin_centers [max(bin_centers)+fixed_bin_width:fixed_bin_width:global_max]];
    end
    
    if min(bin_Edges) > uniform_min
        dist_to_cover = abs(min(bin_Edges)-uniform_min);
        number_of_bins_to_add = dist_to_cover/fixed_bin_width;
        %adding zeros to x in order to covver full range
        x_adder = zeros(1,round(number_of_bins_to_add));
        normalized_counts = [x_adder normalized_counts];
        bin_centers = [[global_min+fixed_bin_width:fixed_bin_width:min(bin_centers)] bin_centers];
    end
    
    bin_Edges = bin_centers - 0.005;
    bin_Edges = [bin_Edges max(bin_centers)+0.005];
    index_bigger = find(bin_centers >uniform_min & bin_centers < uniform_max); %------- CHANGE THIS TO COMPARE TO OTHRT U()-DIST----
    uniform_pdf(index_bigger) = 1/(uniform_max-uniform_min); 
   
    similarity = norm(normalized_counts-uniform_pdf); %L2-norm
end
