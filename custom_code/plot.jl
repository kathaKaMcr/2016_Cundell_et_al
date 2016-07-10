using PyPlot

function plotSet(i, D; folder = Union{}, error = true, scale = true, subplot = false, model = false, vars = ["pS"])
  if folder == Union{}
      pygui(true)
  else
      pygui(false)
      try
          run(`mkdir $folder`)
      end
  end
  if subplot == false
    fig = figure(figsize= (4,4))
    subplot = subplot2grid([1,1], [0,0], 1, 1)
  end
  for cond in specs.conditions
    data = get_data(D, i, cond, scale = scale)
    subplot[:scatter](data[!isnan(data[:,2]),1], data[!isnan(data[:,2]),2], color = specs.colors[cond], label = cond, marker = "x", linewidth = 2)
    if error == true
      subplot[:errorbar](data[!isnan(data[:,3]),1], data[!isnan(data[:,3]),2], data[!isnan(data[:,3]),3], ecolor = specs.colors[cond],elinewidth = 1.5, capthick = 2, linewidth = 0)
    end
  end
  if model != false
    k_ds = D[i,:k_ds]
    k_bg = D[i, :k_bg]
    o = D[i, :o]
    substrate_data = [c => get_data(D, i, c, scale = true) for c in specs.conditions]
    for condition in specs.conditions
      d = substrate_data[condition]
      x = d[!isnan(d[:,2]),1]
      y = d[!isnan(d[:,2]),2]
      model = change_state!(model, condition, ensa_data, substrate_data, depletion_efficiency, b_total)
      model.M.pars["k_bg"] = D[i, :k_bg]
      model.M.pars["k_ds"] = D[i, :k_ds]
      simulate!(model.M, "$condition", collect(0.:45.))
      plotModel(model.M, "$condition"; t = "t", pars = false, vars = vars,  colors = [v => specs.colors[condition] for v in vars], sp = subplot, fig = false)
    end
    subplot[:annotate]("k_bg = $(D[i, :k_bg])", xy = [30,1.32])
    subplot[:annotate]("k_ds = $(D[i, :k_ds])", xy = [30,1.2])
    subplot[:annotate]("o = $(round(D[i, :o],3))", xy = [30,1.08])
  end

  ax = subplot[:axes]
  subplot[:set_ylim]([0,1.5])
  subplot[:set_xlim]([0,47])
  subplot[:set_xlabel]("Time (min)")
  subplot[:set_ylabel]("Phosphorylation (H/L ratio)")
  ax[:spines]["top"][:set_visible](false) # Hide the top edge of the axis
  ax[:spines]["right"][:set_visible](false) # Hide the top edge of the axis
  ax[:xaxis][:set_ticks_position]("bottom") # Set the x-ticks to only the bottom
  ax[:yaxis][:set_ticks_position]("left") # Set the y-ticks to only the left
  ax[:patch][:set_visible](false)

  name = split(D[i, :Gene_names], ";")[1]
  n = D[i, :RecordsIndex]
  p = round(Int, D[i, :Position])
  sequence = D[i, :Sequence_window]
  if sequence == "_"
      fname = name*"_noSequenceFound_"*"$p __$n"
      plottitle = name * " _noSequenceFound"
      subplot[:set_title](plottitle)
  else
      aa = D[i, :Sequence][p]
      sw = sequence
      fname = name*"_p$aa"*"$p __$n"
      plottitle = name * " p$aa" * "$p\n" * sw
      subplot[:set_title](plottitle)
  end
  subplots_adjust(left = 0.16, right = 0.94, hspace = 0.33, wspace = 0.10, top = 0.84, bottom = 0.16)
  if folder != Union{}
      savefig(folder*"/"*fname*".png")
      close()
  end
  return(subplot)
end
