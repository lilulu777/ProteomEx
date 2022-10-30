function [measLength_std, RMSerror_std_all_AVG, RMSerror_std_all_STD] = ...
    compute_group_X_Y(in_group_indices, measData)

% Collect X Y for each group
measLength_std = measData(in_group_indices(1)).standardized_XY.measLength;
RMSerror_std_all = zeros(length(in_group_indices),length(measData(in_group_indices(1)).standardized_XY.RMSerror));
for i = 1:length(in_group_indices)
    RMSerror_std_all(i,:) = measData(in_group_indices(i)).standardized_XY.RMSerror;
end

% Create group statistics
RMSerror_std_all_AVG = mean(RMSerror_std_all,1);
RMSerror_std_all_STD = std(RMSerror_std_all,0,1);