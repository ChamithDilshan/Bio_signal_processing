function MSE = mean_squared_error(noise_free_signal,filtered_signal)
    MSE = sum(power(noise_free_signal-filtered_signal,2));
end
    

