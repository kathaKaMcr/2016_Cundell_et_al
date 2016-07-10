using DataFrames
using JSON
using DataStructures
using XPP

function rescale(S, i, condition)
  if condition == "Control"
    scale_factor = 1
  elseif condition == "B55"
    scale_factor = (S[i, findfirst(names(S), :CM_Con_L_B55_H_)] + 1 / S[i, findfirst(names(S), :CM_Con_H_B55_L_)]) / 2
  elseif condition == "GWL"
    scale_factor = (S[i, findfirst(names(S), :CM_Con_L_GWL_H_)] + 1 / S[i, findfirst(names(S), :CM_Con_H_GWL_L_)]) / 2
  end
  return(scale_factor)
end

function dataRange(S, condition)
  timecourses = Dict(
    "Control" => findfirst(names(S), :Control_0):findfirst(names(S), :Control_45),
    "B55" => findfirst(names(S), :B55_0):findfirst(names(S), :B55_45),
    "GWL" => findfirst(names(S), :GWL_0):findfirst(names(S), :GWL_45),
    )
  return timecourses[condition]
end

function errorRange(S, condition)
  errors = Dict(
    "Control" => findfirst(names(S), :error_Control_0):findfirst(names(S), :error_Control_45),
    "B55" =>  findfirst(names(S), :error_B55_0):findfirst(names(S), :error_B55_45),
    "GWL"=> findfirst(names(S), :error_GWL_0):findfirst(names(S), :error_GWL_45),
    )
  return errors[condition]
end

function get_data(S, i, condition; scale = false)
  i = i
  S = S
  t = convert(Array{Float64,1}, specs.timepoints["values"])
  data_range = dataRange(S, condition)
  y = vec(convert(Array, S[i, data_range]))
  error_range = errorRange(S, condition)
  y_err = vec(convert(Array, S[i, error_range]))
  if scale == true
    scale_factor = rescale(S, i, condition)
    y = y .* scale_factor
    y_err = y_err .* scale_factor
  end
  return([t y y_err])
end
