run_dec = function(
    treefile, geofile, multiplierfile, adjacencyfile=NULL, maxrange=2,
  timesfile, distmatrixfile=NULL, resultsfile, section=F, areasfile=NULL
    ){
  library(rexpokit)
  library(cladoRcpp)
  library(BioGeoBEARS)
  library(ape)
  trfn = treefile
  geogfn = geofile
  
  #######################################################
  # Look at your geographic range data:
  tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
  # Maximum range size observed:
  max(rowSums(dfnums_to_numeric(tipranges@df)))
  
  # Set the maximum number of areas any species may occupy; this cannot be larger 
  # than the number of areas you set up, but it can be smaller.
  max_range_size = maxrange
  
  #######################################################
  # Run DEC
  #######################################################
  
  # Intitialize a default model (DEC model)
  BioGeoBEARS_run_object = define_BioGeoBEARS_run()
  
  # Give BioGeoBEARS the location of the phylogeny Newick file
  BioGeoBEARS_run_object$trfn = trfn
  
  # Give BioGeoBEARS the location of the geography text file
  BioGeoBEARS_run_object$geogfn = geogfn
  
  # Input the maximum range size
  BioGeoBEARS_run_object$max_range_size = max_range_size
  
  BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
  BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
  # Also: search script on "include_null_range" for other places to change
  
  # Set up a time-stratified analysis:
  # 1. Here, un-comment ONLY the files you want to use.
  # 2. Also un-comment "BioGeoBEARS_run_object = section_the_tree(...", below.
  # 3. For example files see (a) extdata_dir, 
  #  or (b) http://phylo.wikidot.com/biogeobears#files
  #  and BioGeoBEARS Google Group posts for further hints)
  #
  # Uncomment files you wish to use in time-stratified analyses:
  if(section==T){
  BioGeoBEARS_run_object$timesfn = timesfile
  }
  BioGeoBEARS_run_object$dispersal_multipliers_fn = multiplierfile
  #BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
  #BioGeoBEARS_run_object$areas_adjacency_fn = adjacencyfile
  #BioGeoBEARS_run_object$distsfn = distmatrixfile
  # See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.
  
  # Speed options and multicore processing if desired
  BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
  BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
  BioGeoBEARS_run_object$use_optimx = TRUE    # if FALSE, use optim() instead of optimx();
  # if "GenSA", use Generalized Simulated Annealing, which seems better on high-dimensional
  # problems (5+ parameters), but seems to sometimes fail to optimize on simple problems
  BioGeoBEARS_run_object$num_cores_to_use = 1
  # (use more cores to speed it up; this requires
  # library(parallel) and/or library(snow). The package "parallel" 
  # is now default on Macs in R 3.0+, but apparently still 
  # has to be typed on some Windows machines. Note: apparently 
  # parallel works on Mac command-line R, but not R.app.
  # BioGeoBEARS checks for this and resets to 1
  # core with R.app)
  
  # Sparse matrix exponentiation is an option for huge numbers of ranges/states (600+)
  # I have experimented with sparse matrix exponentiation in EXPOKIT/rexpokit,
  # but the results are imprecise and so I haven't explored it further.
  # In a Bayesian analysis, it might work OK, but the ML point estimates are
  # not identical.
  # Also, I have not implemented all functions to work with force_sparse=TRUE.
  # Volunteers are welcome to work on it!!
  BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
  
  # This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
  # (It also runs some checks on these inputs for certain errors.)
  BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
  
  # Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
  if(section==T){
    BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
  }
  # The stratified tree is described in this table:
  #BioGeoBEARS_run_object$master_table
  
  # Good default settings to get ancestral states
  BioGeoBEARS_run_object$return_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run
  
  # Set up DEC model
  # (nothing to do; defaults)
  
  # Look at the BioGeoBEARS_run_object; it's just a list of settings etc.
  #BioGeoBEARS_run_object
  
  # This contains the model object
  #BioGeoBEARS_run_object$BioGeoBEARS_model_object
  
  # This table contains the parameters of the model 
  #BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table
  
  # Run this to check inputs. Read the error messages if you get them!
  check_BioGeoBEARS_run(BioGeoBEARS_run_object)
  
  # For a slow analysis, run once, then set runslow=FALSE to just 
  # load the saved result.
  runslow = TRUE
  resfn = resultsfile
  if (runslow)
  {
    res = bears_optim_run(BioGeoBEARS_run_object)
    res    
    
    save(res, file=resfn)
    resDEC = res
  } else {
    # Loads to "res"
    load(resfn)
    resDEC = res
  }
  
}

get_table = function(objs){
  tbl = sapply(objs, function(obj){
    obj = obj$optim_result
    this_row = c(obj$p1, obj$p2, obj$value)
  })
  tbl = t(tbl)
  colnames(tbl) = c("d", "e", "lnL")
  return(tbl)
}

plot_models = function(DEC.fit, title){
  library(phytools)
  ## subdivide our plotting area using layout 
  layout(matrix(1:2,1,2),widths=c(0.2,0.8)) ## set plotting parameters 
  par(mar=c(4.1,0.1,3.1,0.1),cex=0.8)
  ## plot legend and tree 
  plot_BioGeoBEARS_results(DEC.fit, analysis_titletxt=title, plotlegend=TRUE,
                           tipcex=0.6,statecex=0.4, plotwhat = "pie")
}
