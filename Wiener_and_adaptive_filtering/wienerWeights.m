function w0 = wienerWeights(order, Y, N)

    expected_YYT = zeros(order,order);
    expected_Yy = zeros(order,1);
    for k = order:length(Y)-order+1
        temp_y = (Y(k:k+order-1))';
        expected_YYT = expected_YYT + temp_y*temp_y';
        expected_Yy = expected_Yy + temp_y * temp_y(end);
    end
    expected_YYT = expected_YYT / (length(Y)-2*order+2);
    expected_Yy = expected_Yy / (length(Y)-2*order+2);

    expected_NNT = zeros(order,order);
    for k = order:length(N)-order+1
        temp_n = (N(k:k+order-1))';
        expected_NNT = expected_NNT + temp_n*temp_n';
    end
    expected_NNT = expected_NNT/ (length(N)-2*order+2);
    
    w0 = pinv(expected_YYT + expected_NNT)*expected_Yy;
end
