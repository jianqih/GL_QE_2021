cfunction[corr] = correlationCoefficient(X,Y,n)
    sum_X = 0; sum_Y = 0; sum_XY = 0;
    squareSum_X = 0; squareSum_Y = 0;
    
    for i=1:n
        % // sum of elements of array X.
        sum_X = sum_X + X(i);
    
        %/ sum of elements of array Y.
        sum_Y = sum_Y + Y(i);
    
        %// sum of X[i] * Y[i].
        sum_XY = sum_XY + X(i) * Y(i);
    
       % // sum of square of array elements.
        squareSum_X = squareSum_X + X(i) * X(i);
        squareSum_Y = squareSum_Y + Y(i) * Y(i);
    end
    
    % // use formula for calculating correlation coefficient.
    corr = (n * sum_XY - sum_X * sum_Y) / sqrt((n * squareSum_X - sum_X * sum_X) * (n * squareSum_Y - sum_Y * sum_Y));
end