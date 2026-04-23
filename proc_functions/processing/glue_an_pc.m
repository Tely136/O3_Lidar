function [data_glued, scaling_coeffs, corrcoefs, valid] = glue_an_pc(an_data,pc_data,min_toggle,max_toggle)
    toggle_mask = pc_data > min_toggle & pc_data < max_toggle;

    pc_glue_valid = false(size(pc_data));
    pc_glue_valid(toggle_mask) = true;
    
    data_glued = pc_data;
    scaling_coeffs = NaN(size(pc_data,2),2);
    corrcoefs = NaN(size(pc_data,2),1);
    valid = false(size(pc_data));
    
    for profile_idx = 1:size(an_data,2)
        an_temp = an_data(:,profile_idx);
        pc_temp = pc_data(:,profile_idx);
    
        valid_id = find(pc_glue_valid(:,profile_idx));
        [b,r] = scale_analog(an_temp(valid_id), pc_temp(valid_id));
        scaling_coeffs(profile_idx,:) = b;
        corrcoefs(profile_idx) = r;

        an_scaled_temp = an_temp.*b(2) + b(1);
        data_glued(pc_temp > max_toggle,profile_idx) = an_scaled_temp(pc_temp > max_toggle);
        valid(valid_id,profile_idx) = true;
    end
end