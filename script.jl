using PPSIM

specs, streams, model = PPSIM.initialise(pwd());

processed = "/home/luki/Documents/Data/Phosphoproteomics/table_processed.csv"
extracted = joinpath(specs.sink, "extracted.csv")
sequences = joinpath(specs.sink, "sequences.json")
ontology = joinpath(specs.sink, "pantherGeneList.tsv")
cleaned =  joinpath(specs.sink, "cleaned.csv")
fitted = joinpath(specs.sink, "fitted.csv")
substrates = joinpath(specs.sink, "substrates.csv")
independent_substrates = joinpath(specs.sink, "independent_substrates.csv")


close("all")
preprocess(specs.source, processed)
extract(processed, extracted)
clean(extracted, cleaned, sequences, ontology)
S = readtable(cleaned)
parametrise(cleaned, fitted)
filter_substrates(fitted, substrates, 0.85, "k_ds", "k_bg")
filter_substrates(fitted, independent_substrates, 0.85, "k_bg", "k_ds")
