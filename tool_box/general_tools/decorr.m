function [A_decorr] = decorr(A, B)
%     % Decorrelates vector A with respect to vector B using Gram-Schmidt process.
%     %
%     % Inputs:
%     % A - The first vector (will be decorrelated with respect to B)
%     % B - The second vector (remains unchanged)
%     %
%     % Outputs:
%     % A_decorr - The decorrelated version of A
%     % B_unchanged - B (returned unchanged for convenience)
%     % Morici Juan Facundo 09/2024
    
    % Check if vectors A and B are same length
    if length(A) ~= length(B)
        error('Vectors A and B must have the same length');
    end
    
    B_norm = B / norm(B);% Normalize vector B
    
    projection_A_on_B = dot(A, B_norm) * B_norm;% Project vector A onto vector B
    
    % Subtraction of the projection from A
    % Decorrelation A from B
    A_decorr = A - projection_A_on_B;

end