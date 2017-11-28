using FerrenbergSwendsenReweighting
using Base.Test


#Generate a multiple histogram reweighting object:
@test begin
    series = 3
    λs = linspace(10,1,series)
    τints = ones(series)
    xs = [0.1*i .+ 0.01.*randn(10000+i) for i = 1:series]
    rw = ReweightObj(λs, xs)
    w = evaluate(rw, 5.0)
    true
end
