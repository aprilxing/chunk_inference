% Example matrices
matrixA = rand(161, 15);
matrixB = rand(161, 15);

% Find the maximum value and corresponding row index for each column
[maxValues, rowIndices] = max(gamma');

% Display the results
disp('Maximum values in each column:');
disp(maxValues);

disp('Corresponding row indices in matrixB:');
disp(rowIndices);

% If you want the actual values in matrixB corresponding to the maximum values
maxValuesInB = matrixB(sub2ind(size(matrixB), rowIndices, 1:size(matrixB, 2)));
disp('Actual values in matrixB corresponding to maximum values in matrixA:');
disp(maxValuesInB);
%%

testt = zeros(161, 5);

for i = 1:161
    testt(i, :) = chunks(rowIndices(i), :)';

end 