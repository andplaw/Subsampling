function xy_sub = MFNUS(xy,fc,K)
% Bengt's Non-Uniform Sample Elimination


arguments
    xy
    fc = 1.5
    K = 10
end

if size(xy,1)<size(xy,2)
    xy = xy';
end

%% algorithm

N = length(xy(:,1)); % Get the number of its dots

xy = sortrows(xy,2); % Sort dots from bottom and up
                       % Create nearest neighbor pointers and distances
[Idx,D] = knnsearch(xy,xy,'K',K);
for k = 1:N % Loop over nodes from bottom and up
    if Idx(k,1)~=0 % Check if node already eliminated
        ind = find(D(k,2:end) < fc*D(k,2));
        ind2 = Idx(k,ind+1);
        ind2(ind2<k) = []; % Mark nodes above present one, and which
        Idx(ind2,1) = 0;   % are within the factor fc of the closest one
    end                % in the original node set
end
xy_sub = xy(Idx(:,1)~=0,:); % Eliminate these marked nodes

end