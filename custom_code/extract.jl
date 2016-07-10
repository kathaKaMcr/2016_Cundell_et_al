using DataFrames
using JSON
using DataStructures
using XPP

function describe_data(D)
    De = DataFrame()
    De[:x] = ["N", "Mean", "Min", "Q25", "Q50", "Q75", "Max", "NAs", "NA%", "unique", "error"]
    for t in names(D)[2:end]
        d = D[t]
        dnona = dropna(d)
        if length(dnona) >= 1
            qs = quantile(dnona, [0, 0.25, 0.5, 0.75, 1])
            err = std(dnona)/length(dnona)
        else
            qs = [NaN, NaN, NaN, NaN, NaN]
            err = NaN
        end
        De[t] = [length(d), mean(dnona), qs, sum(isna(d)), round(sum(isna(d))*100/length(d), 2), length(unique(dnona)), err]
    end
    return(De)
end

function get_timecourse_data(data_table, n, t_hdrs, extracted_folder, key = "Ratio_H_L_normalized")
  ix = round(Int, ceil(n/2))
  if mod(n,2) == 1
    id = "S"
    pattern = specs.sampletypes["Supernatant"]
  elseif mod(n,2) == 0
    id = "P"
    pattern = specs.sampletypes["Pellet"]
  end
  columns = PPSIM.queryColumns(pattern, names(data_table))
  columns = PPSIM.queryColumns(Regex(key), columns)
  experiments = PPSIM.queryColumns(pattern, [symbol(e) for e in specs.experiment_ids], level = "inner")
  cm = Symbol[symbol(cm) for cm in specs.crossmixes]
  data = DataFrame()
  data[:x] = [string(e) for e in experiments]
  for s in [t_hdrs; cm]
    data[s] = zeros(length(experiments)) .* NA
  end
  for (i, e) in enumerate(data[:x])
    columns_exp = PPSIM.queryColumns(Regex(e),columns)
    for cond in specs.conditions
      columns_cond = PPSIM.queryColumns(Regex(cond),columns_exp)
      for t in specs.timepoints["hdr"]
        columns_t = PPSIM.queryColumns(Regex(t),columns_cond)
        try
          data[i, symbol("$cond\_" * t[2:end])] = data_table[ix, columns_t[1]]
        end
      end
    end
    for cm in specs.crossmixes
      columns_cm = PPSIM.queryColumns(Regex(cm),columns_exp)
      try
        data[i, symbol(cm)] = data_table[ix, columns_cm[1]]
      end
    end
  end
  summary_stat = describe_data(data)
  for i in 1:size(summary_stat)[1]
    push!(data, convert(Array, summary_stat[i, :]))
  end
  writetable(joinpath(extracted_folder, "$ix\_$id.csv"), data)
  i_data = findfirst(data[:x], "Mean")
  d_data = data[i_data, 2:end]
  i_error = findfirst(data[:x], "error")
  d_error = data[i_error, 2:end]
  return(d_data, d_error)
end

function extract(sourcefile, sinkfile)
  extracted_folder = joinpath(specs.sink, "extracted")
  try
    run(`mkdir $extracted_folder`)
  end
  data = readtable(sourcefile)
  hdrs = [:RecordsIndex, :Gene_names, :Protein, :Position]
  t_hdrs = Symbol[]
  for c in specs.conditions
    for t in specs.timepoints["hdr"]
      t = t[2:end]
      t_hdrs = [t_hdrs; symbol("$c\_$t")]
    end
  end
  cm_hdrs = Symbol[symbol(cm) for cm in specs.crossmixes]
  hdrs = [hdrs; t_hdrs; cm_hdrs]
  error_hdrs = map(x -> symbol("error_" * string(x)), hdrs[5:end])
  hdrs = [hdrs; error_hdrs]
  extracted = DataFrame()
  for h in hdrs[1:3]
    extracted[h] = ["_" for i in 1:(2 * size(data)[1])]
  end
  for h in hdrs[4:end]
    extracted[h] = zeros(2 * size(data)[1]) .* NA
  end
  writetable("test.csv", extracted)
  for i in 1:(2 * size(data)[1])
    println(i)
    ix = round(Int, ceil(i/2))
    extracted[i, :RecordsIndex] = string(round(Int, ceil(i/2))) * "_" * (mod(i,2) == 1 ? "S": "P")
    for h in [:Gene_names, :Position, :Protein]
      extracted[i,h] = data[ix,h]
    end
    extracted[i, [t_hdrs; cm_hdrs]], extracted[i, error_hdrs] = get_timecourse_data(data, i, t_hdrs, extracted_folder)
  end
  writetable(sinkfile, extracted)
end
