using DynamicalSystems

function logisticmap(dx, x, p, n)
    dx[1] = p[1]*x[1]*(1-x[1])
    return
end

data_size = 3000

ds                      = DiscreteDynamicalSystem(logisticmap, rand(1), 4.0)
Data_Det                = trajectory(ds,data_size; Ttr=10000)[1][:,1][2:end]
