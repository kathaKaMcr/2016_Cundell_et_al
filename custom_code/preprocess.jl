using PPSIM
using DataFrames

"""
**Preprocessing Routine**
Takes path to sourcefile, and path to sinkfile as arguments.
Queries column names of table and re-structures them so that numbers at the end
of column headers are  attached to the first part of the header preceding the
experiment id, as specified in specs.headerpattern
"""
function preprocess(sourcefile, sinkfile)
  data = readtable(sourcefile, separator = '\t')
  exp_indices = PPSIM.queryColumns(specs.headerpattern, data)
  # println(exp_indices)
  #Find the ones that have a tailing 1-3 preceded by an underscore
  indices = PPSIM.queryColumns(r"_[1-3]\b", exp_indices)
  for index in indices
    index = string(index)
    #get the number
    n = index[end]
    #remove the tail
    tail_removed = split(index, r"_[1-3]\b")[1]
    #Get part preceding specs.headerpattern
    firstpart = split(index, specs.headerpattern)[1]
    experiment = split(tail_removed, Regex(firstpart))[2]
    newname = firstpart * "_$n\_" * experiment
    rename!(data, symbol(index), symbol(newname))
  end
  writetable(sinkfile, data)
end
