function y = calc_ppc_20140923( theta )
%
% usage:
%
% function to calculate the pairwise phase consistency

num_angles = length(theta);
if num_angles > size(theta, 2); theta = theta'; end   % make sure theta is a column vector

angle_matrix = zeros(num_angles-1);

running_sum = 0;
for i_angle1 = 1 : num_angles - 1
%     tempArray1 = circshift(theta, i_angle1 - 1);
%     tempArray1 = tempArray1(1:end-1);
    for i_angle2 = i_angle1 + 1 : num_angles
        
%         tempArray2 = circshift(theta, i_angle2 - 1);
%         tempArray2 = tempArray2(1:end-1);
        
        f = cos(theta(i_angle1)).*cos(theta(i_angle2)) + sin(theta(i_angle1)).*sin(theta(i_angle2));
        running_sum = running_sum + f;
    end

end

y = 2 * running_sum / (num_angles * (num_angles-1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

