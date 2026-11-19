function cloud = cloud_screening(P_nr,P_fr,z,merge_ids,settings)
    overlap = merge_ids(1):merge_ids(2);
    scale = median(P_fr(overlap,:)./P_nr(overlap,:),1,'omitnan');
    P = merge_profiles(P_nr.*scale,P_fr,merge_ids(1),merge_ids(2));
    P(:,isnan(scale)) = P_fr(:,isnan(scale));

    frame = settings.frame_length;
    [~,coeff] = sgolay(2,frame);
    C = repmat({coeff},1,numel(z));
    smoothed = movmean(P,frame,1);
    gradient = smooth_derivative(log(smoothed.*z.^2),C,z(2)-z(1));

    in_range = z >= settings.range_km(1)*1e3 & z <= settings.range_km(2)*1e3;
    cloud.mask = false(size(P));
    cloud.layers_km = cell(1,size(P,2));
    cloud.gradient = gradient;

    for i = 1:size(P,2)
        g = gradient(:,i);
        edges = diff([false; in_range & isfinite(g); false]);
        starts = find(edges == 1);
        stops = find(edges == -1)-1;
        layers = zeros(0,2);

        for j = 1:numel(starts)
            ids = (starts(j)+1:stops(j)-1)';
            rising = g(ids) > settings.positive_threshold & g(ids) > g(ids-1) & g(ids) >= g(ids+1);
            falling = g(ids) <= -settings.negative_threshold & g(ids) < g(ids-1) & g(ids) <= g(ids+1);
            positive_peak = [];

            for k = ids(rising | falling)'
                if g(k) > 0
                    positive_peak = k;
                elseif ~isempty(positive_peak)
                    base = find(g(starts(j):positive_peak-1) <= 0,1,'last') + starts(j)-1;
                    top = find(g(k+1:stops(j)) > -0.002,1) + k;
                    if ~isempty(base) && ~isempty(top)
                        cloud.mask(base:top,i) = true;
                        layers(end+1,:) = z([base top])'/1e3;
                    end
                    positive_peak = [];
                end
            end
        end
        cloud.layers_km{i} = layers;
    end
end
