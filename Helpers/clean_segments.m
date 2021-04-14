function peaks = clean_segments(maxk, maxk_loc, mink, mink_loc)
%CLEAN_SEGMENTS Summary of this function goes here
%   Detailed explanation goes here
peaks = []; % [max, maxloc, min, minloc]

for i = 1:length(maxk)

    if (i < length(maxk)) % before last maximum
        min_interval = mink_loc(mink_loc > maxk_loc(i) & mink_loc < maxk_loc(i+1));
        
        if (length(min_interval) == 1) % Exactly one minimum
            min_idx = find(mink_loc==min_interval);
            peaks = [peaks; maxk(i), maxk_loc(i), mink(min_idx), mink_loc(min_idx)];

        elseif (length(min_interval) > 1) % More than one minimum
            min_idx = find(mink==max(mink(ismember(mink_loc,min_interval))));
            peaks = [peaks; maxk(i), maxk_loc(i), mink(min_idx), mink_loc(min_idx)];
        end
    else % last maximum
        min_idx = find(mink_loc>maxk_loc(i), 1);
        if (~isempty(min_idx))
            peaks = [peaks; maxk(i), maxk_loc(i), mink(min_idx), mink_loc(min_idx)];
        end
    end
end

peaks(:,3) = -peaks(:,3);

end

