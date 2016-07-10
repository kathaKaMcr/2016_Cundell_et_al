using DataFrames
using JSON
using DataStructures
using XPP

function trimSequence(seq, pos::Int64, leftflank; rightflank = Union{})
    if rightflank == Union{}
      rightflank = leftflank
    end
    l = maximum([1, pos - leftflank])
    u = minimum([length(seq), pos + rightflank])
    n_headingx = -minimum([0, pos - leftflank - 1])
    n_tailingx = -minimum([0, length(seq) - (pos + rightflank)])
    trimmedseq = ("_" ^ n_headingx) * seq[l:u] * ("_" ^ n_tailingx)
    if length(trimmedseq) != (rightflank + leftflank + 1)
        "Error: $seq"
    end
    return(trimmedseq)
end

function clean(sourcefile, sinkfile, seqfile, ontology, logfile = "cleaning_data_log.json", cutoff_n_data = 4)
  D = readtable(sourcefile)
  cleaning_log = OrderedDict()
  cleaning_log["EXTRACTED: total number of entries"] = size(D)[1]
  cleaning_log["EXTRACTED: total number of unique peptides"] = length(unique([i[1:end-1] for i in D[:RecordsIndex]]))
  cleaning_log["EXTRACTED: entries missing T0 in Control"] = sum(isnan(D[:Control_0]))
  cleaning_log["EXTRACTED: entries missing T0 in B55"] = sum(isnan(D[:B55_0]))
  cleaning_log["EXTRACTED: entries missing T0 in GWL"] = sum(isnan(D[:GWL_0]))
  D = D[D[:Control_0] .>= 0, :]
  D = D[D[:B55_0] .>= 0, :]
  D = D[D[:GWL_0] .>= 0, :]
  cleaning_log["CLEANING: total number of entries left after removing entries missing T0"] = size(D)[1]
  cleaning_log["CLEANING: total number of unique peptides left after removing entries missing T0"] = length(unique([i[1:end-1] for i in D[:RecordsIndex]]))
  D = D[!isna(D[:Gene_names]), :]
  sufficient_data = Int64[]
  for i in 1:size(D)[1]
    n_data = Bool[]
    for c in specs.conditions
      l = findfirst(names(D), symbol(c * "_0"))
      u = findfirst(names(D), symbol(c * "_45"))
      n_data = [n_data; sum(convert(Array, D[i, l:u]) .>= 0) >= cutoff_n_data]
    end
    if all(n_data)
      sufficient_data = [sufficient_data; i]
    end
  end
  D = D[sufficient_data, :]
  cleaning_log["CLEANING: total number of entries left after removing entries with less than $(cutoff_n_data)"] = size(D)[1]
  cleaning_log["CLEANING: total number of unique peptides left after removing entries with less than $(cutoff_n_data)"] = length(unique([i[1:end-1] for i in D[:RecordsIndex]]))
  D = D[D[:CM_Con_L_B55_H_] .>= 0, :]
  D = D[D[:CM_Con_H_B55_L_] .>= 0, :]
  D = D[D[:CM_Con_L_GWL_H_] .>= 0, :]
  D = D[D[:CM_Con_H_GWL_L_] .>= 0, :]

  cleaning_log["CLEANING: total number of entries left after removing entries missing CMs involving Control"] = size(D)[1]
  cleaning_log["CLEANING: total number of unique peptides left after removing entries missing CMs involving Control"] = length(unique([i[1:end-1] for i in D[:RecordsIndex]]))

  data_ok = Bool[]
  for i in 1:size(D)[1]
    d = vec(convert(Array, D[i, findfirst(names(D), :Control_0):findfirst(names(D), :GWL_45)]))
    data_ok = [data_ok; !any(d .> 10)]
  end
  D = D[data_ok,:]
  cleaning_log["CLEANING: total number of entries left after removing entries with timepoints > 10"] = size(D)[1]
  cleaning_log["CLEANING: total number of unique peptides left after removing entries with timepoints > 10"] = length(unique([i[1:end-1] for i in D[:RecordsIndex]]))
  D[:Gene_names] = convert(Array, D[:Gene_names], "_")
  D[:Protein] = convert(Array, D[:Protein], "")
  D[:Sequence] = ["_" for i in 1:size(D)[1]]
  D[:Sequence_window] = ["_" for i in 1:size(D)[1]]

  f = open(sequences, "r")
  seqs = JSON.parse(f)
  close(f)
  background_set = []
  for i in 1:size(D)[1]
    pid = D[i, :Protein]
    if pid in keys(seqs)
      seq = seqs[pid]
      pos = round(Int, D[i, :Position])
      background_set = [background_set; trimSequence(seqs[pid], pos, 25, rightflank = 10)]
    end
  end

  f = open(joinpath(specs.sink, "background_set.txt"), "w")
  write(f, join(background_set, "\n"))
  close(f)

  f = open(joinpath(specs.sink,"background_set_unique.txt"), "w")
  write(f, join(unique(background_set), "\n"))
  close(f)

  for i in 1:size(D)[1]
    pid = D[i, :Protein]
    if pid in keys(seqs)
      D[i, :Sequence] = seqs[pid]
      pos = round(Int, D[i, :Position])
      D[i, :Sequence_window] = trimSequence(seqs[pid], pos, 25, rightflank = 10)
    end
  end

  D = D[D[:Gene_names] .!= "KRT7", :]
  D = D[D[:Gene_names] .!= "KRT8", :]
  D = D[D[:Gene_names] .!= "KRT18", :]
  D = D[D[:Gene_names] .!= "_", :]
  cleaning_log["CLEANING: total number of entries left after removing Keratin and peptides lacking gene_names"] = size(D)[1]
  cleaning_log["CLEANING: total number of unique peptides left after removing Keratin and peptides lacking gene_names"] = length(unique([i[1:end-1] for i in D[:RecordsIndex]]))



  panther = readtable(ontology);
  #annotate fitted-table with panther go-terms
  D[:, :Molecular_function] = ["_" for i in 1:size(D)[1]]
  D[:, :Biological_process] = ["_" for i in 1:size(D)[1]]
  D[:, :Cellular_component] = ["_" for i in 1:size(D)[1]]

  for i in 1:size(D)[1]
    gene_name = split(D[i,:Gene_names], ";")[1]
    pn = panther[panther[:Gene_names] .== gene_name, :]
    if size(pn)[1] > 0
      for c in [:Molecular_function,:Biological_process, :Cellular_component]
        D[i,c] = isna(pn[1,c]) ? "_" : pn[1,c]
      end
    end
  end


  writetable(sinkfile, D)
  f = open(joinpath(specs.sink, logfile), "w")
  write(f, JSON.json(cleaning_log))
  close(f)
  println(cleaning_log)

end
