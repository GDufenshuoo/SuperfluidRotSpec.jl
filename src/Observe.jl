function measure(vars, obs, weights, config) 
    # obs: prototype of the observables for each integral
        x, bin = vars #unpack the variables
        obs[1][bin[1]] += weights[1] # circle
        obs[2][bin[1]] += weights[2] # sphere
end