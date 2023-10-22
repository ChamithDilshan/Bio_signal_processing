function linearModel = createLinearModel(len)

linearModel = []

for i = 1:len
    
    if (1<=i && i<21)
        linearModel = [linearModel 0]
    elseif (21<=i && i<24)
        linearModel = [linearModel 0.1(i-21)]
    end
end