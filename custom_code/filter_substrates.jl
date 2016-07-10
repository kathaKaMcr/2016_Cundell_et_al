using DataFrames
using JSON
using DataStructures
using XPP

function filter_substrates(infile, outfile, q, rate1, rate2)
  S = readtable(infile)
  dir, fname  = splitdir(outfile)
  fname = split(fname, ".csv")[1]

  #Instantiate log for filtering
  filter_log = OrderedDict()
  filter_log["FILTER: total number of entries"] = size(S)[1]

  S = S[S[:o] .<= 0.005, :]
  filter_log["FILTER: total number of entries with o <= 0.005"] = size(S)[1]

  S = S[S[symbol(rate1)] .> 0, :];
  filter_log["FILTER: total number of entries with $rate1  > 0"] = size(S)[1]

  S = S[S[symbol(rate1)] .> S[symbol(rate2)], :]
  filter_log["FILTER: total number of entries with $rate1 > $rate2"] = size(S)[1]

  S = S[S[symbol(rate1)] .> quantile(S[symbol(rate1)], q), :]
  filter_log["FILTER: total number of entries with $rate1 in $q percentile of dataset"] = size(S)[1]

  sort!(S, cols = symbol(rate1), rev = true)
  writetable(joinpath(dir, "$fname\.csv"), S)
  writetable(joinpath(dir, "$fname\_top.csv"), S[1:20,:])

  s = S[:Sequence]
  p = S[:Position]
  pAA = []
  for i in 1:size(s)[1]
    if length(s[i]) >= p[i]
      pAA = [pAA; s[i][p[i]]]
    else
      pAA = [pAA; "_"]
    end
  end

  S_s = S[pAA .== 'S', :]
  writetable(joinpath(dir, "$fname\_pS.csv"), S_s)
  S_t = S[pAA .== 'T', :]
  writetable(joinpath(dir, "$fname\_pT.csv"), S_t)

  f = open(joinpath(dir,"$fname\.json"), "w")
  write(f, JSON.json(filter_log))
  close(f)
end
