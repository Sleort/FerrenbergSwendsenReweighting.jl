
module AutocorrelationTime

import StatsBase:autocor
export integrated_autocorrelation_time



"""
Find the integrated autocorrelation time of a vector of measurements

Based on the tips given in http://www.stat.unc.edu/faculty/cji/Sokal.pdf
"""
function integrated_autocorrelation_time(samples::Vector)
  c = 6.0 #Some constant determining the tradeoff between large bias (small c) and large variance (large c). Okay value should be around c ∈ [4,10]
  Δtstart = 1
  Δtstop = 30
  ac = autocor(samples, Δtstart:Δtstop)

  τint = 1/2
  M = Δtstart #Window width (to be increased until minimum width is reached)

  for i=1:30
    if c*τint < M
      break
    end
    τint += ac[M - Δtstart + 1]
    M+=1

    if M>Δtstop #More autocorrelation calculation needed
      Δtstart = Δtstop + 1
      Δtstop = 3Δtstop - Δtstart + 1
      ac = autocor(samples, Δtstart:Δtstop)
    end
  end

  return τint
end







end # module
