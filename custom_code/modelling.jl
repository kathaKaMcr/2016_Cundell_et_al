using DataFrames
using JSON
using DataStructures
using XPP
using Optim

"""
Change state of the model (parameters and initial conditions) based on experimental data
"""
function change_state!(model::PPSIM.Model, condition, ensa_data, substrate_data, depletion_efficiency, b_total)
  #Dict mapping experimental data to model state for each condition
  model_states = Dict(
    "Control" => Dict(
      "init" => Dict(
        "pEt"  => ensa_data["Control"][1,2],
        "pS"  => substrate_data["Control"][1,2],
        ),
      "pars" => Dict(
        "Bt"  => b_total,
        )),
    "B55" => Dict(
      "init" => Dict(
        "pEt"  => ensa_data["B55"][1,2],
        "pS"  => substrate_data["B55"][1,2],
        ),
      "pars" => Dict(
        "Bt"  => b_total * depletion_efficiency,
        )),
    "GWL" => Dict(
      "init" => Dict(
        "pEt"  => ensa_data["GWL"][1,2],
        "pS"  => substrate_data["GWL"][1,2],
        ),
      "pars" => Dict(
        "Bt"  => b_total,
        ))
      )

  # Loop over state definition and change the corresponding model field
  state_definition = model_states[condition]
  for (p,d) in state_definition
    for (k,v) in d
      model.M.(symbol(p))[k] = v
    end
  end
  #Ensure that pEB is bound
  model.M.init["pEB"] = minimum([model.M.init["pEt"], model.M.pars["Bt"]])
  return(model)
end


function objective_function(y_sim, y_exp, scale_exp)
    diff = (y_sim - y_exp) ./ scale_exp
    lsq = sum(diff .^2) / (2 * length(y_exp))
    return(lsq)
end


function concatenate_data(data_dict)
  y_exp = Float64[]
  t_exp = Float64[]
  scale_exp = Float64[]
  for c in specs.conditions
    d = data_dict[c]
    t = d[!isnan(d[:,2]),1]
    y = d[!isnan(d[:,2]),2]
    err = d[!isnan(d[:,2]),3]
    scale = ones(length(t)) * y[1]
    y_exp = [y_exp; y]
    t_exp = [t_exp; t]
    scale_exp = [scale_exp; scale]
  end
  return(t_exp, y_exp, scale_exp)
end

"""
Return a vector of simulation data for the variable *key* that is aligned in time with experimental timepoints as specified by x
"""
function get_y_data_simulation(M, name, x, key)
  # Get simulation time and simulation data for "key"
  all_t = M.sims[name].D["t"]
  all_y = M.sims[name].D[key]
  # Find the indices in the simulation time vector corresponding to the timepoints specified by x
  indices = [findfirst(all_t, xi) for xi in x]
  # Get the corresponding simulation data
  y = [all_y[i] for i in indices]
  return(y)
end

"""
Simulate ENSA and PRC1 data for a given parameter set and return the value of the objective function for the simulation given the data
"""
function simulate_ensa_prc1(model::PPSIM.Model, k_ass::Float64, k_diss::Float64, k_cat::Float64, b_total::Float64, depletion_efficiency::Float64, k_ds::Float64, k_bg::Float64, ensa_data::Dict, substrate_data::Dict; plot = false)

  #Concatenate experimental data for evalualtion of the objective function
  x_ensa, y_ensa, scale_ensa = concatenate_data(ensa_data)
  x_substrate, y_substrate, scale_substrate = concatenate_data(substrate_data)
  x_exp = [x_ensa; x_substrate]
  y_exp = [y_ensa; y_substrate]
  scale_exp = [scale_ensa; scale_substrate]

  #Instantiate a list holding the ENSA-simulations and the Substrate-simulations
  y_ensa = []
  y_substrate = []

  if plot == true
    figure()
    sp1 = subplot(121)
    sp2 = subplot(122)
  end

  #Loop over conditions
  for condition in specs.conditions
    #Get vectors of t and y for non-missing data
    d_s = substrate_data[condition]
    x_s = d_s[!isnan(d_s[:,2]),1]
    y_s = d_s[!isnan(d_s[:,2]),2]
    d_e = ensa_data[condition]
    x_e = d_e[!isnan(d_e[:,2]),1]
    y_e = d_e[!isnan(d_e[:,2]),2]

    #Change the model state (initial conditions based on experimental data)
    model = change_state!(model, condition, ensa_data, substrate_data, depletion_efficiency, b_total)

    #Change the parameters
    model.M.pars["k_ass"] = k_ass
    model.M.pars["k_diss"] = k_diss
    model.M.pars["k_cat"] = k_cat
    model.M.pars["k_bg"] = k_bg
    model.M.pars["k_ds"] = k_ds

    #Simulate
    # runSimulation!(model.M, "$condition")
    simulate!(model.M, "$condition", collect(0:0.1:46))
    if plot == true
      sp1[:scatter](x_e, y_e, color = specs.colors[condition])
      plotModel(model.M, "$condition", vars = ["pEt"], colors = ["pEt" => specs.colors[condition]], fig = false, sp = sp1)
      sp2[:scatter](x_s, y_s, color = specs.colors[condition])
      plotModel(model.M, "$condition", vars = ["pS"], colors = ["pS" => specs.colors[condition]], fig = false, sp = sp2)
    end

    #Get simulation data for timepoints specified by x
    y_ensa = [y_ensa; get_y_data_simulation(model.M, "$condition", x_e, "pEt")]
    y_substrate = [y_substrate; get_y_data_simulation(model.M, "$condition", x_s, "pS")]
  end

  #Concatenate ENSA- and Substrate data from simulations
  y_sim = [y_ensa; y_substrate]

  #Evaluate objective function
  o = objective_function(y_sim, y_exp, scale_exp)
  return(o)
end
simulate_ensa_prc1(model::PPSIM.Model, ensa_data::Dict, substrate_data::Dict; plot = false) = simulate_ensa_prc1(model, model.M.pars["k_ass"], model.M.pars["k_diss"], model.M.pars["k_cat"],model.M.pars["Bt"], model.M.pars["depletion_efficiency"], model.M.pars["k_ds"], model.M.pars["k_bg"], ensa_data::Dict, substrate_data::Dict; plot = plot)

"""
Fitting routine to find best parameter set of model given ENSA and PRC1 data
"""
function fit_routine_ensa(PRC1, model, ensa_data)
  #Instantiate columns holding the substrate-sepecific parameters
  PRC1[:k_ds] = zeros(size(PRC1)[1])
  PRC1[:k_bg] = zeros(size(PRC1)[1])
  PRC1[:o] = zeros(size(PRC1)[1])

  #Define box and starting values for parameters in fitting routine
  initial_guess = [8, 0.01, 0.25, 0.5, 0.051, 0.0002, 0.0002]
  lower_bound = [0., 0.0001, 0., 0.001, 0.05, 0.0001, 0.]
  upper_bound = [50., 0.1, 1., 1., 1., 1., 1.]

  #Get the data for PRC1
  substrate_data = [c => get_data(PRC1, 1, c, scale = true) for c in specs.conditions]

  #Convert simulation-function into local function only depending of parameter vector and make it differentiable
  f(p) = simulate_ensa_prc1(model, p[1], p[2], p[3], p[4], p[5], p[6], p[7], ensa_data, substrate_data)
  g = DifferentiableFunction(f)

  #Run box-minimization algorithm according to specifications above
  r = fminbox(g, initial_guess, lower_bound,  upper_bound, show_trace = true)

  #Print results
  println(r.minimum)
  println("="^40)

  #Wrap results in dictionary
  pars = Dict(
   "k_ass" => round(r.minimum[1],4),
   "k_diss" => round(r.minimum[2],4),
   "k_cat" => round(r.minimum[3],4),
   "Bt" => round(r.minimum[4],4),
   "depletion_efficiency" => round(r.minimum[5],4),
   "o" => round(r.f_minimum,4)
  )
  PRC1[:k_ds] = round(r.minimum[6],4)
  PRC1[:k_bg] = round(r.minimum[7],4)
  PRC1[:o] = round(r.f_minimum,4)
  return(pars, PRC1)
end

"""
FUNCTION FOR pS(t)
"""
function substrate(p, t, Bt, Si)
  #Approximate the integral of Bt over time via finite sum (rectangle method)
  dt = [0; t[2:end] - t[1:end-1]]
  T = dt .* Bt
  FB = cumsum(T)
  pS = exp(- p[2] * t - p[1] .* FB) .* Si
  return(pS)
end

function get_substrate_simulation_data(pS, t_sim, t_exp)
  indices = [findfirst(t_sim, t) for t in t_exp]
  y = [pS[i] for i in indices]
  return(y)
end

function simulate_substrate(p, substrate_data, b55_sim, t_sim)
  y_sim = Float64[]
  t_exp_concat, y_exp, scale_exp = concatenate_data(substrate_data)
  for condition in specs.conditions
    Si = substrate_data[condition][1,2]
    pS = substrate(p, t_sim[condition], b55_sim[condition], Si)
    d = substrate_data[condition]
    t_exp = d[!isnan(d[:,2]),1]
    y_sim = [y_sim; get_substrate_simulation_data(pS, t_sim[condition], t_exp)]
  end
  o = objective_function(y_sim, y_exp, scale_exp)
  return(o)
end

function fit_routine_substrate(D, model, depletion_efficiency, b_total, ensa_data, key)

  #Simulate ENSA-B55 module to get b55-data and simulation-time to approximate integral of B55 over time
  for condition in specs.conditions
    d = ensa_data[condition]
    x = d[!isnan(d[:,2]),1]
    y = d[!isnan(d[:,2]),2]
    model = change_state!(model, condition, ensa_data, ensa_data, depletion_efficiency, b_total)
    simulate!(model.M, "$condition", collect(0.:0.1:45.))
  end

  b55_sim = [c => convert(Array{Float64,1}, deepcopy(model.M.sims[c].D["B"])) for c in specs.conditions]
  t_sim = [c => convert(Array{Float64,1}, deepcopy(model.M.sims[c].D["t"])) for c in specs.conditions]

  D[:k_ds] = zeros(size(D)[1])
  D[:k_bg] = zeros(size(D)[1])
  D[:o] = zeros(size(D)[1])

  for s in 1:size(D)[1]

    println(s, "\t", D[s,:RecordsIndex], "\t", D[s, :Gene_names], "\t", D[s, :Position])

    substrate_data = [c => get_data(D, s, c, scale = true) for c in specs.conditions]
    f(p) = simulate_substrate(p, substrate_data, b55_sim, t_sim)
    g = DifferentiableFunction(f)
    r = fminbox(g, [0.0001, 0.0001], [1e-5,1e-5],  [0.3, 0.3], ftol = 1e-5)
    k_ds = round(r.minimum[1],4)
    k_bg = round(r.minimum[2],4)
    o = round(r.f_minimum,4)
    println("k_ds = $k_ds \t k_bg = $k_bg \t o = $o")
    println("="^40)
    D[s, :k_ds] = k_ds
    D[s, :k_bg] = k_bg
    D[s, :o] = o
  end
  writetable(joinpath(specs.sink, "fitted.csv"), D)
end


function parametrise(sourcefile, sinkfile)
  # Setup fitting of ENSA and PRC1 pT470 -data
  dir = splitdir(sinkfile)[1]
  S = readtable(sourcefile);
  ENSA = S[S[:Gene_names] .== "ENSA;ARPP19", :]
  ensa_data = [c => get_data(ENSA, 1, c, scale = true) for c in specs.conditions]
  PRC1 = S[S[:Gene_names] .== "PRC1", :]
  PRC1 = PRC1[PRC1[:Position] .== 470, :]
  prc1_data = [c => get_data(PRC1, 1, c, scale = true) for c in specs.conditions]

  pars, PRC1 = fit_routine_ensa(PRC1, model, ensa_data)
  depletion_efficiency = model.M.pars["depletion_efficiency"]
  b_total = pars["Bt"]
  for (p,v) in pars
    model.M.pars[p] = v
    ENSA[symbol(p)] = v
  end
  simulate_ensa_prc1(model, ensa_data, prc1_data, plot = true)

  ensa_sim = [c => convert(Array{Float64,1}, deepcopy(model.M.sims[c].D["pEt"])) for c in specs.conditions]
  t_sim = [c => convert(Array{Float64,1}, deepcopy(model.M.sims[c].D["t"])) for c in specs.conditions]
  for cond in specs.conditions
    for t in collect(0:2.5:45)
      hdr = "Sim_" * cond * replace(string(t), ".", "_")
      ENSA[symbol(hdr)] = zeros(size(ENSA)[1])
      ind = findfirst(t_sim[cond], t)
      ENSA[symbol(hdr)] = ensa_sim[cond][ind]
    end
  end

  writetable(joinpath(dir, "ENSA.csv"), ENSA)

  b55_sim = [c => convert(Array{Float64,1}, deepcopy(model.M.sims[c].D["B"])) for c in specs.conditions]
  t_sim = [c => convert(Array{Float64,1}, deepcopy(model.M.sims[c].D["t"])) for c in specs.conditions]
  fit_routine_substrate(S, model, depletion_efficiency, b_total, ensa_data, "pS")


  global depletion_efficiency = model.M.pars["depletion_efficiency"]
  global b_total = model.M.pars["Bt"]

  k_ass = model.M.pars["k_ass"]
  k_diss = model.M.pars["k_diss"]
  k_cat = model.M.pars["k_cat"]

  for cond in specs.conditions
    for t in collect(0:2.5:45)
      hdr = "Sim_" * cond * replace(string(t), ".", "_")
      S[symbol(hdr)] = zeros(size(S)[1])
    end
  end

  for i in 1:size(S)[1]
    pS = [c => substrate([S[i, :k_ds], S[i,:k_bg]], t_sim[c], b55_sim[c], S[i,symbol(c * "_0")]) for c in specs.conditions]
    for cond in specs.conditions
      for t in collect(0:2.5:45)
        hdr = "Sim_" * cond * replace(string(t), ".", "_")
        ind = findfirst(t_sim[cond], t)
        S[i, symbol(hdr)] = pS[cond][ind]
      end
    end
  end
  writetable(sinkfile, S)
end
