function ma3ECG_1 = moving_average_filter(N,x)
    num_points = size(x, 2);
    ma3ECG_1 = zeros(size(x));

    for i = 1:num_points
        y_i = 0;
        for k = 0:N-1
            if i-k>0
                y_i = y_i + x(i-k);
            end
        end
        ma3ECG_1(i) = y_i/ N;
    end
end

