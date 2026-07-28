function m = repeatedMedian(d)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

x = d(:,1);
y = d(:,2);

N = length(x);

slopeMat = zeros(N, N-1);

for j = 1:size(slopeMat, 1)

        slopeMat(j, :) = (y(1:length(y) ~= j) - y(j)) ./ ...
            (x(1:length(x) ~= j) - x(j));

end

m = median(median(slopeMat, 2, "omitmissing"), 1, 'omitmissing');

end