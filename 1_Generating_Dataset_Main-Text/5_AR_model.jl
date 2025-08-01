using Distributions

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

data_size = 3000
phi = 0.5

ARgenerator(data_size,rand(),phi,0.1)