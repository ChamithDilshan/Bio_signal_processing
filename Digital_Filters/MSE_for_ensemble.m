function MSE_k = MSE_for_ensemble(k,epochs,template)
    N = size(epochs,1);
    y_hat_k = mean(epochs(:,(1:k)),2);
    MSE_k = sqrt((sum(power(template-y_hat_k,2)))/N);
end