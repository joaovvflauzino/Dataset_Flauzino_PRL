using Distributions, DynamicalSystems, Statistics

function logisticmap(dx, x, p, n)
    dx[1] = p[1]*x[1]*(1-x[1])
    return
end

function NormalWhiteNoise(n, variancia)    
    noise = rand(Normal(0.0,sqrt(variancia)), n)
    return noise
end

function ARgenerator(data_size,x0,a,varAR)

    Transient = 1000
    x=x0
    ARseries = zeros(data_size)

	varNoise = varAR*(1 - a^2)

	Noise = NormalWhiteNoise(data_size+Transient,varNoise)

    for loop_i=1:Transient
        x = a*x + Noise[loop_i]
    end
    for loop_i=1:data_size
        x = a*x + Noise[loop_i]
        ARseries[loop_i] = x
    end

    return ARseries
end

function Normalize(V_Mean_Window_A)
    Maximum_Value,Local_Number=findmax(V_Mean_Window_A)
    Minimum_Value,Local_Number=findmin(V_Mean_Window_A)
    if ((Maximum_Value-Minimum_Value) != 0.0)
        V_Mean_Window_A.=((V_Mean_Window_A.-Minimum_Value)./(Maximum_Value-Minimum_Value))
    end
    return V_Mean_Window_A
end

data_size = 3000

ds                      = DiscreteDynamicalSystem(logisticmap, rand(1), 4.0)
Data_Det                = trajectory(ds,data_size; Ttr=10000)[1][:,1][2:end]

phi = 0.5
variance_noise = 0.1*var(Data_Det)

Data_Noise = ARgenerator(data_size, rand(), phi, variance_noise)
Data_In    = Normalize(Data_Det .+ Data_Noise)