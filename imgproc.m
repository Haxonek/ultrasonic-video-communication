disp('Starting');

im = imread('img/lower_res.png');
im = uint8(im);
im = rgb2gray(im);
[m, n] = size(im);

% this referes to the number of different pre-determined plots we'll create
% and test. 
num_of_con = 2048;
% this referes to how many times the image will be broken up in each
% concentration, it must match with what is hard coded in the function
squares_height = 6;
squares_width = 4;
squares_in_con = squares_height*squares_width;

% create map of conentrations of pixels
base = 32;
concentrations = ones(num_of_con,squares_in_con);
for mc = 1:num_of_con
    for i = squares_in_con+4:4:base
        index = randi(squares_in_con);
        concentrations(mc,index) = concentrations(mc,index) + 4;
    end
    
    if (sum(concentrations(mc,:)) ~= base)
        fprintf('Issue with concentration %i\n', mc);
    end
    
    % this is so rare now we may as well remove it
    if (rowMatch(concentrations(mc,:), concentrations, mc) ~= 0)
        % we want to generate a new row instead and overwrite this one
        fprintf('Removing matching row at %i from concentrations\n', mc);
        mc = mc - 1;
    end
end

% Create value map
% 1269 is for 15fps, .19% of pixels
disp('Building pre-set tables')
tic
value_maps = zeros(num_of_con,m,n,'uint8');
for mc = 1:num_of_con
    value_maps(mc,:,:) = generate_partitioned_pepper2(m, n, 7812, squeeze(concentrations(mc,:)), squares_height, squares_width, base);
%     value_maps(mc,:,:) = generate_partitioned_pepper(m, n, 7812/2, squeeze(concentrations(mc,:)));
end
toc
%% Create Sobel edge map
% create salt map from image using sobel edge detection
image_salt_map = uint8(edge(im, 'sobel'));
image_salt_map = uint8(image_salt_map);

%% decide on which value map to use

% comparing image salt map with value maps
% higher values = more similarity
disp('Comparing all values to value')
tic
ssim_comparison = 1;
ssim_worst = -1;
% value_map = squeeze(value_maps(1,:,:));

% edge pixels / edge pixels in max square
% salt_ratio = generate_salt_ratio(image_salt_map, squares_height, squares_width);
% min_conc_ratio = min(concentrations(1,:)); % should basically always be 1
% conc_ratios = zeros(num_of_con, squares_in_con);
% 
% for i = 1:num_of_con
%     conc_ratios(i,:) = (concentrations(i,:)-min_conc_ratio)./(max(concentrations(i,:))-min_conc_ratio);
% end
% 
% [match, index] = find_closest_match(conc_ratios, salt_ratio);
% 
% value_map = squeeze(value_maps(index,:,:));
% NEW method; rank them in order, fit ranks with optimization model

% generate salt rank order
salt_rank = generate_salt_rank(image_salt_map, squares_height, squares_width);

% generate concentration rank order
conc_rank = generate_conc_rank(concentrations, squares_height, squares_width);

% finds the best match to the original salt rank
% change to a linear programming model later, but just using immse for now
best_index = 1;
best_immse_value = 1;
for i = 1:num_of_con
    % get immse value
    tmp_immse_value = immse(salt_rank, squeeze(conc_rank(i,:)));
    % if it's better than old value, replace index and immse value
    if tmp_immse_value < best_immse_value
        best_index = i;
    end
end
disp('found best value map')
% final should be set as value_map
value_map = squeeze(value_maps(best_index,:,:));
toc

% for i = 1:num_of_con
%     % create concentration ratio
%     conc_ratio = concentrations(i,:)./max(concentrations(i,:));
%     
%     % Run linear program to compare
%     
%     
% %     tmp = ssim(image_salt_map, squeeze(value_maps(i,:,:)));
% %     tmp = immse(image_salt_map, squeeze(value_maps(i,:,:)));
% 
%     if (tmp < ssim_comparison)
%         ssim_comparison = tmp;
%         value_map = squeeze(value_maps(i,:,:));
%     end
%     
%     if (tmp > ssim_worst)
%         ssim_worst = tmp;
%         worst_value_map = squeeze(value_maps(i,:,:));
%     end
% end
% toc
% disp('SSIM values: '); 
% disp(ssim_comparison);
% disp(ssim_worst);

figure();
imshowpair(image_salt_map, value_map, 'montage');
title('salt image and best value map');
disp('Showing salt image vs best value map')
% 
% figure();
% subplot(1,1,1);
% imshowpair(image_salt_map, worst_value_map, 'montage');
% title('salt image and worst value map');
%% Rebuild image and display it (comparing it to the other images)

% values I would send, but in this case we're just leaving them in the
% vector
compressed_map = generate_compressed_image(value_map, im);

disp('Starting reconstruction of image')
tic
% now we want to rebuild the image from the compressed map, 
[reconstructed_map_1, new_value_map_1] = reconstruct_image(compressed_map, value_map, 7);
disp('first computed')
toc
% repeat with new values
[reconstructed_map_2, new_value_map_2] = reconstruct_image(reconstructed_map_1, new_value_map_1, 9);
disp('second computed')
toc
% repeat with new values
[reconstructed_map_3, new_value_map_3] = reconstruct_image(reconstructed_map_2, new_value_map_2, 11);
disp('third computed')
toc
% repeat with new values
[reconstructed_map_3, new_value_map_3] = reconstruct_image(reconstructed_map_3, new_value_map_3, 17);
disp('fourth computed')
toc
% repeat with new values
[reconstructed_map_3, new_value_map_3] = reconstruct_image(reconstructed_map_3, new_value_map_3, 23);
disp('fifth computed')
toc

figure();
imshowpair(im, reconstructed_map_3,'montage');
title('Original im and final reconstructed')

% figure()
% imshowpair(reconstructed_map_1, reconstructed_map_2,'montage');
% title('reconstructed 1 and 2')





% Generates a matrix of ratios, where the values are pixel ratio / max
% square pixel ratio, the max being 1 and the smallest (likely) being 1 /
% max.
function ratio2 = generate_salt_rank(map, rows, cols)
    ratio = zeros(rows,cols);
    
    [mp,np] = size(map);
    
    m = ones(1,rows).*round(mp/rows);
    for i = 2:rows-1
        m(i) = i*m(1);
    end
    m(rows) = mp;
    
    n = ones(1,cols).*round(np/cols);
    for i = 2:cols-1
        n(i) = i*n(1);
    end
    n(cols) = np;
        
    % iterate through each partitan, assign number of pixels to it equal to
    % the ratio assigned in concentrations
    for i = 1:rows
        for j = 1:cols
            
            if i == 1 && j == 1
                cur = map(1:m(i), 1:n(j));
                ratio(i,j) = sum(sum(cur));
            elseif i == 1
                cur = map(1:m(i), n(j-1)+1:n(j));
                ratio(i,j) = sum(sum(cur));
            elseif j == 1
                cur = map(m(i-1)+1:m(i), 1:n(j));
                ratio(i,j) = sum(sum(cur));
            else % neither == 1
                cur = map(m(i-1)+1:m(i), n(j-1)+1:n(j));
                ratio(i,j) = sum(sum(cur));
            end
        end
    end
    
    % flatten 
    ratio2 = zeros(1, rows*cols);
    for i = 1:rows
        for j = 1:cols
            ratio2(i*j) = ratio(i,j);
        end
    end
    
    % remove straglers
    for i = 1:rows*cols
        % set all values that might be misinterpreted as a rank as the max
        % value
        if (ratio2(i) < rows*cols)
            ratio2(i) = rows*cols+1;
        end
    end
    
    % order by rank
    rank = 1;
    for i = 1:rows*cols
        cur_max = max(ratio2);
        for j = 1:rows*cols
            if (ratio2(j) == cur_max)
                ratio2(j) = rank;
                rank = rank + 1;
                break;
            end
        end
    end
    disp('done')
end

function [conc_rank] = generate_conc_rank(concentrations, squares_height, squares_width)

    [num,len] = size(concentrations);
    
    conc_rank = zeros(num,len);

    for i = 1:1 % num -- just do one for now
        m_max_queue = java.util.LinkedList;
        n_max_queue = java.util.LinkedList;
        max_queue = java.util.LinkedList;
        
        % create max heap
        % add every value
        % for each value above 0, itr through
        
        for cur_m = 1:squares_height
            tmp_total = sum(squeeze(concentrations(i,(((cur_m-1)*squares_width)+1):cur_m*squares_width)));
            has_added_row_sum = false; % has not yet added row, obvi
            if m_max_queue.isEmpty()
                m_max_queue.add([cur_m,tmp_total]);
            else
                tmp_index = [0,0];
                for s = 0:m_max_queue.size()-1
                    % 1 is index, 2 is sum value
                    cur_vals = m_max_queue.get(s);
                    if (cur_vals(2) < tmp_total && ~has_added_row_sum)
                        % assign to temp to be shifted
%                         tmp_index = cur_vals;
                        % set sum to this spot
                        tmp_index = m_max_queue.set(s,[cur_m,tmp_total]);
                        % record that the row sum has already been added
                        has_added_row_sum = true;
                    elseif has_added_row_sum
                        % shifting all values down
                        tmp_index = m_max_queue.set(s,tmp_index);
                    end
                end
                
                if has_added_row_sum
                    % add final value, extending size of ll by one
                    m_max_queue.add(tmp_index);
                else
                    m_max_queue.add([cur_m,tmp_total]);
                end
            end
        end
        
        % do the same for n row 
        for cur_n = 1:squares_width
            tmp_total = sum(squeeze(concentrations(i,(((cur_n-1)*squares_height)+1):cur_n*squares_height)));
            has_added_row_sum = false; % has not yet added row, obvi
            if n_max_queue.isEmpty()
                n_max_queue.add([cur_n,tmp_total]);
            else
                tmp_index = [0,0];
                for s = 0:n_max_queue.size()-1
                    % 1 is index, 2 is sum value
                    cur_vals = n_max_queue.get(s);
                    if (cur_vals(2) < tmp_total && ~has_added_row_sum)
                        % assign to temp to be shifted
%                         tmp_index = cur_vals;
                        % set sum to this spot
                        tmp_index = n_max_queue.set(s,[cur_n,tmp_total]);
                        % record that the row sum has already been added
                        has_added_row_sum = true;
                    elseif has_added_row_sum
                        % shifting all values down
                        tmp_index = n_max_queue.set(s,tmp_index);
                    end
                end
                
                if has_added_row_sum
                    % add final value, extending size of ll by one
                    n_max_queue.add(tmp_index);
                else
                    % has never added row sum, adding now
                    n_max_queue.add([cur_n,tmp_total]);
                end
            end
        end
        
        % now we need to take those two lists and recreate 
        display(m_max_queue.size());
        display(n_max_queue.size());
        disp('sizes (n then m)')
        disp(concentrations(i,1:4));
        disp(concentrations(i,5:8));
        disp(concentrations(i,9:12));
        disp(concentrations(i,13:16));
        disp(concentrations(i,17:20));
        disp(concentrations(i,21:24));
        disp('printed map')
        
        % build the rank index        
        for j = 1:squares_height
            cur_m = m_max_queue.get(j-1);
            
            % add for each n value that's above 1, loop
            for k = 1:squares_width
                cur_n = n_max_queue.get(k-1);

                cur_val = cur_m(2) * cur_n(2);
                
                %%%% add sorted %%%%
                
                has_added_row_sum = false; % has not yet added row, obvi
                if max_queue.isEmpty()
                    max_queue.add([cur_m(1),cur_n(1),cur_val]);
                else
                    tmp_index = [0,0,0];
                    for s = 0:max_queue.size()-1
                        % 1 and 2 are the m and n index, 3 is sum value
                        cur_vals = max_queue.get(s);
                        if (cur_vals(3) < cur_val && ~has_added_row_sum)
                            % assign to temp to be shifted
    %                         tmp_index = cur_vals;
                            % set sum to this spot
                            tmp_index = max_queue.set(s,[cur_m(1),cur_n(1),cur_val]);
                            % record that the row sum has already been added
                            has_added_row_sum = true;
                        elseif has_added_row_sum
                            % shifting all values down
                            tmp_index = max_queue.set(s,tmp_index);
                        end
                    end

                    if has_added_row_sum
                        % add final value, extending size of ll by one
                        max_queue.add(tmp_index);
                    else
                        max_queue.add([cur_m(1),cur_n(1),cur_val]);
                    end
                end
                %%%%  %%%%
            end
        end
        
        disp('size of max_queue')
        disp(max_queue.size())
        
        counter = 1;
        for c = 0:(squares_width*squares_height - 1)
            cur_val = max_queue.get(c);
            cur_index = squares_width*(cur_val(1) - 1) + cur_val(2);
            
            if conc_rank(i,cur_index) ~= 0
                disp('ERROR - value already exists')
            end
            
            conc_rank(i,cur_index) = counter;
            counter = counter + 1;
            
        end
        
        % end
        disp(conc_rank(i,1:4));
        disp(conc_rank(i,5:8));
        disp(conc_rank(i,9:12));
        disp(conc_rank(i,13:16));
        disp(conc_rank(i,17:20));
        disp(conc_rank(i,21:24));
        disp('done with one')
    end
end

% Generates a matrix of ratios, where the values are pixel ratio / max
% square pixel ratio, the max being 1 and the smallest (likely) being 1 /
% max.
function ratio2 = generate_salt_ratio(map, rows, cols)
    ratio = zeros(rows,cols);
    
    [mp,np] = size(map);
    
    m = ones(1,rows).*round(mp/rows);
    for i = 2:rows-1
        m(i) = i*m(1);
    end
    m(rows) = mp;
    
    n = ones(1,cols).*round(np/cols);
    for i = 2:cols-1
        n(i) = i*n(1);
    end
    n(cols) = np;
        
    % iterate through each partitan, assign number of pixels to it equal to
    % the ratio assigned in concentrations
    for i = 1:rows
        for j = 1:cols
            
            if i == 1 && j == 1
                cur = map(1:m(i), 1:n(j));
                ratio(i,j) = sum(sum(cur));
            elseif i == 1
                cur = map(1:m(i), n(j-1)+1:n(j));
                ratio(i,j) = sum(sum(cur));
            elseif j == 1
                cur = map(m(i-1)+1:m(i), 1:n(j));
                ratio(i,j) = sum(sum(cur));
            else % neither == 1
                cur = map(m(i-1)+1:m(i), n(j-1)+1:n(j));
                ratio(i,j) = sum(sum(cur));
            end
        end
    end
    
    ratio = ratio ./ (max(ratio));
    
    % flatten 
    ratio2 = zeros(1, rows*cols);
    for i = 1:rows
        for j = 1:cols
            ratio2(i*j) = ratio(i,j);
        end
    end
end

% This function takes in the compressed image and the value map and spits
% out a reconstructed map by convolving the image and calculating the mean
% square
function [reconstructed_map,new_value_map] = reconstruct_image(compressed_map, value_map, kernel)
    if (mod((kernel - 1),2) ~= 0)
        disp('kernel must be an odd numer'); % also must be >= 3
        return;
    end

    new_value_map = value_map; % will include new values to repeat
    reconstructed_map = compressed_map; % 
    
    [mr, nr] = size(value_map);
    k = (kernel - 1) / 2;
    
    for mc = 1:mr
        for nc = 1:nr
            if (value_map(mc,nc) == 0)
                % missing data in compressed_map
                mc_b = mc - k;
                if (mc_b <= 0)
                    mc_b = 1;
                end
                
                mc_t = mc + k;
                if (mc_t > mr)
                    mc_t = mr;
                end
                
                nc_b = nc - k;
                if (nc_b <= 0)
                    nc_b = 1;
                end
                
                nc_t = nc + k;
                if (nc_t > nr)
                    nc_t = nr;
                end
                
                % calculate reduced maps of just the area we're convolving
                reduced_compression_map = compressed_map(mc_b:mc_t, nc_b:nc_t);
                reduced_value_map = value_map(mc_b:mc_t, nc_b:nc_t);
                
                intensity = mean_square_value(reduced_compression_map, reduced_value_map, mc, nc);
                if (intensity >= 0)
                    % found nearby values
                    new_value_map(mc,nc) = 1;
                    reconstructed_map(mc,nc) = sum(intensity);
                end
            end
        end
    end
end

% This function calculates the intensity of the gray color pixel based off
% values around it TODO - this algorithm is crazy slow, don't recalculate
% every value each time
function intensity = mean_square_value(reduced_compression_map,reduced_value_map, mx, ny)
    intensity = -1;
    is_empty = 0;
    
    [mi,ni] = size(reduced_value_map);
    
    
    d = zeros(mi*ni); % equal in length to all values being valid
    i = zeros(mi*ni);
    counter = 1; % tracks position in d
    dt = 0;

    % iterate through maps, find values, calculate intensity, return
    for mc = 1:mi
        for nc = 1:ni
            if (reduced_value_map(mc,nc) == 1)
                is_empty = is_empty + 1; % no longer empty
                % calculate the distance, add it to dt, add it to d
                d(counter) = sqrt((mx-mc)^2 + (ny-nc)^2);
                i(counter) = reduced_compression_map(mc,nc); % stores gray value
                dt = dt + d(counter);
                counter = counter + 1; % increment counter
            end
        end
    end
    
    if (is_empty == 0)
        % must have usable values
        return;
    end
    
    % sum the intensity
    intensity = sum(i.*d)/dt;
end

% takes in two vectors of equal size and swaps the values from the second
% into the first vector when the first vector has a value of 1
function im_map = generate_compressed_image(pepper_map, im)
    im_map = zeros(size(pepper_map), 'uint8');
    [mi, ni] = size(pepper_map);
    
    if (size(pepper_map) ~= size(im))
        disp('ERROR: maps must be of equal size');
        return;
    end
    
    % look up which direction matlab stores arrays; may be able to make it
    % faster by swapping height and width
    for mc = 1:mi
        for nc = 1:ni
            % if 1, swap value with im
            if (pepper_map(mc,nc) == 1)
                im_map(mc,nc) = im(mc,nc);
            else 
                im_map(mc,nc) = 255;
            end
        end
    end
end

function pepper_map = generate_partitioned_pepper2(mp, np, total, concentrations, rows, cols, base)

    if (size(concentrations) ~= rows*cols)
        disp('concentrations must equal to rows*cols')
    end
    
    pepper_map = zeros(mp,np,'uint8');
    
    m = ones(1,rows,'uint16').*round(mp/rows);
    for i = 2:rows-1
        m(i) = i*m(1);
    end
    m(rows) = mp;
    
    n = ones(1,cols,'uint16').*round(np/cols);
    for i = 2:cols-1
        n(i) = i*n(1);
    end
    n(cols) = np;
    
%     base = 80;
    
    % iterate through each partitan, assign number of pixels to it equal to
    % the ratio assigned in concentrations
    for i = 1:rows
        for j = 1:cols
            if i == 1 && j == 1
                pepper_map(1:m(i), 1:n(j)) = generate_binary_pepper(m(i),n(j),round(total*concentrations(i*j)/base));
            elseif i == 1
                pepper_map(1:m(i), n(j-1)+1:n(j)) = generate_binary_pepper(m(i),n(j)-n(j-1),round(total*concentrations(i*j)/base));
            elseif j == 1
                pepper_map(m(i-1)+1:m(i), 1:n(j)) = generate_binary_pepper(m(i)-m(i-1),n(j),round(total*concentrations(i*j)/base));
            else % neither == 1
                pepper_map(m(i-1)+1:m(i), n(j-1)+1:n(j)) = generate_binary_pepper(m(i)-m(i-1),n(j)-n(j-1),round(total*concentrations(i*j)/base));
            end
        end
    end
    
end

% first value should be vector of the correct size, second value should be
% how many pixels you want to allocate. May change value_map into just size
% results
function pepper_map = generate_binary_pepper(mp, np, total)
    count = 0;
    pepper_map = zeros(mp, np, 'uint8');
    
    while (count < total)
        mc = randi(mp);
        nc = randi(np);
        
        if (pepper_map(mc,nc) > 0)
            % don't use, already has value
            continue;
        else
            % set pepper, increment counter
            pepper_map(mc, nc) = 1;
            count = count + 1;
        end
    end
end


function bool = rowMatch(needle, haystack, index)
    bool = 0;
    
    [mr,nr] = size(haystack);
    
    for mc = 1:(index-1)
        for nc = 1:nr
            if (needle(nc) ~= haystack(mc,nc))
                % different
                break;
            end
            
            % check if we're at the end, if so, fail row
            if (nc == nr)
                disp(mc)
                % match
                bool = mc;
                return;
            end
        end
    end
end

% returns the vector with the closest match, mean square
function [match, index] = find_closest_match(set, key)

    [m,n] = size(set);

    min_val = 100;
    
    for i = 1:m
        tmp = immse(squeeze(set(i,:)), key);
        if tmp < min_val
            min_val = tmp;
            index = i;
        end
    end
    
    match = squeeze(set(index,:));
    
end


