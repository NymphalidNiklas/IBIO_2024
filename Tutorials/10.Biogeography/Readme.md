# Ancestral area reconstruction

Historical biogeography is a field in which we try to understand the evolutionary history of a clade as it relates to the changes in distribution of the species (and their ancestors) over time. In essence, we estimate the ancestral species ranges on a time-calibrated tree. 

We will carry out such an analysis on our Felidae dataset using an R package [BioGeoBEARS](https://github.com/nmatzke/BioGeoBEARS) and the so-called DEC model (dispersal-extinction-cladogenesis) ([Ree et al. 2005](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.0014-3820.2005.tb00940.x), [Ree & Smith 2008](https://academic.oup.com/sysbio/article/57/1/4/1703014)).

To carry out such an analysis, in addition to a time-calibrated tree, we also need the information on the current distribution of the species in our clade of interest, as well as some information on the possibility of dispersal between the different areas relevant for our dataset.

BioGeoBEARS requires a table with the species distribution data formatted in a specific way:

```
42	5 (O	AT	PA	N	NT)
CM017348_Lynx_canadensis	00010
HM185183_Prionailurus_bengalensis	10100
JF357970_Panthera_tigris_sumatrae	10000
JF357972_Panthera_tigris_corbetti	10000
KJ508413_Panthera_tigris	10100
KJ866876_Panthera_pardus_japonensis	10000
KM236783_Panthera_onca	00001
KP001495_Panthera_leo	01000
KP001507_Panthera_pardus	01100
```

This format loosely resembles the PHYLIP file format that you used with your alignments: in the first line we see the dimensions of the data, meaning the number of taxa that this dataset includes and the areas defined. Finally, within the parantheses, we see the labels for the five geographical areas that we chose. We are going with these large biogeographical areas because many of the cat species are quite widely distributed, and we can leave out Australia because there are no naturally occurring cats there.

They are:
+ O: Oriental
+ AT: Afrotropical
+ PA: Palaearctic
+ N : Nearctic
+ NT: Neotropical


Apart from these two files, we need some other files that define some of the parameters or constraints of the models we are going to use. BioGeoBEARS has several options for creating different models for how ancestral ranges evolved in a group of organisms. For example, one can define which areas are allowed by including an [areas-allowed matrix](http://phylo.wikidot.com/biogeobears#areas_allowed_matrix2). One can also have an [area-adjacency matrix](http://phylo.wikidot.com/biogeobears#areas_adjacency_matrix2). What we will use in our analysis are the so-called [dispersal multipliers](http://phylo.wikidot.com/biogeobears#dispersal_multipliers2) in the form of a _dispersal multiplier matrix_. This is a table that modifies the probability of dispersal between different areas. In our case, for example, we don't expect that the probability of dispersal from the Oriental Region into the Neotropical Region is the same as the probability of dispersal from the Nearctic into the Neotropics. 

Example dispersal multiplier matrix (made pretty for readability):

```
O         AT      PA       N      NT
1       0.01       1    0.01    0.01
0.01       1    0.01    0.01    0.01
1       0.01       1       1    0.01
0.01    0.01       1       1       1
0.01    0.01    0.01       1       1

END

```

Furthermore, we can also create a time-stratified model. We know that the position of continents and the existence of various barriers or land bridges has changed over time. This information can be incorporated into our model by using a _time periods_ file as well as _dispersal multiplier matrices_ that are different for each time period. For our dataset, we will divide our tree into three sections, going from the root of the tree: 20-5 million years ago (mya), 5-2 mya, and 2-0 mya. The file with this information is called _time periods_.

Download the time_periods.txt file from [here](../../Data/Day4/Biogeography/) and save it into your working directory. Open the file in a text editor to see how it is organized. These time periods are made based on the time-calibrated tree that we inferred in Tutorial 6.

We also need to decide on the maximum number of areas that one species (as well as the ancestors) is allowed to occupy in our models. Looking at our tip data, there is one species that is present in all areas and that is the domestic cat _Felis catus_. However, this species became cosmopolitan with the help of humans and it makes this analysis more difficult so we decided to remove that tip (and the associated entry from the species distribution data). Looking at all of the other species, at most they are present in two areas so we want to set the maximum number of areas to 2. That way we limit the total number of area combinations that need to be explored by the model to 16 versus 32.

You can download the pruned tree without the cat tip and all the other files we will need for this tutorial from [here](../../Data/Day4/Biogeography/). 

You can also download the species distribution file from [here](../../Data/Day4/Biogeography/). The file is named Felidae_biogeo_data.txt. Open the file in a text editor to see how it is organized.


--------------------

## Using BioGeoBEARS

### Installation 

First we should start by installing the required R package if you haven't done so already. To do that, open RStudio (or launch R the way you normally use it) and type the following command into your console:

```R
install.packages("rexpokit")
install.packages("cladoRcpp")
install.packages("remotes")
```

Then we will install BioGeoBEARS from GitHub.

```R
remotes::install_github(repo="nmatzke/BioGeoBEARS")
```

If R asks you if you want to update some packages first, you can just hit enter without typing anything and it should skip updating packages.
This installation will try to compile some packages so it is likely that it will take some time or fail altogether. If you are stuck at this step, just ask for help and we will try to solve the problem together.

After a successful installation, you should ***close your RStudio program and reopen it*** before we continue.


### Running DEC models

We should start by running the following line of code by copying and pasting. This will allow you to use a wrapper function to use BioGeoBEARS in a more user friendly way.

```R
source("https://raw.githubusercontent.com/NymphalidNiklas/EB2_2024/main/Data/Day4/Biogeography/run_dec_f.R")
```

This wrapper function has several arguments we need to pass values to. Have a look at the code below to get an understanding of how this works (but don't run this, this is just for clarification).

```R
run_dec = function(treefile, ## Path to your tree file
                   geofile,  ## Path to your species distribution data file
                   multiplierfile, ## Path to file that include multiplier matrices
                   adjacencyfile=NULL, ## Not used
                   maxrange=2, ## Maximum number of areas a species allowed to occupy during modeling
                   timesfile, ## File with time periods defined to use it in time-stratified analyses
                   distmatrixfile=NULL, ## Not used
                   resultsfile, ## filename to save your model fit results, should end with ".RData"
                   section=FALSE, ## Set this to true if you are doing a time-stratified analysis
                   areasfile=NULL ## Not used
    )
```

#### Model testing to answer relevant questions

We can use BioGeoBEARS to answer specific hypotheses about our study group. In the case of cats, we could ask the question of how good the cats are at dispersal? What kind of signal is there in our data? 

We can design two models, one in which dispersal between the continents doesn't depend on the geography (DEC0) and one in which we make dispersal between some continents easier than between others based on the current continet positions (DEC1). Our expectation is that dispersal for cats is much more difficult from e.g. the Nearctic to the Oriental Region than from the Nearctic to the Neotropics and this is reflected in the dispersal multiplier matrix for DEC1.

Download the DEC0 and DEC1 dispersal multiplier matrices from [here](../../Data/Day4/Biogeography/) into your working directory and open them in a text editor to examine them. Do the differences between these models make sense given the current position of continents?

#### Comparison of DEC0 and DEC1 models - does limiting dispersal between some continents improve model fit?

Have a look at the values in the tables and think about whether they make sense given the current position of continents and the large barriers to dispersal (e.g., oceans).

DEC0

|    | O  | AT | PA | N  | NT |
|----|----|----|----|----|----|
| O  | 1  | 1  | 1  | 1  | 1  |
| AT | 1  | 1  | 1  | 1  | 1  |
| PA | 1  | 1  | 1  | 1  | 1  |
| N  | 1  | 1  | 1  | 1  | 1  |
| NT | 1  | 1  | 1  | 1  | 1  |

DEC1

|    | O    | AT   | PA   | N    | NT   |
|----|------|------|------|------|------|
| O  | 1    | 0.1  | 1    | 0.01 | 0.01 |
| AT | 0.1  | 1    | 1    | 0.01 | 0.01 |
| PA | 1    | 1    | 1    | 0.1  | 0.01 |
| N  | 0.01 | 0.01 | 0.1  | 1    | 1    |
| NT | 0.01 | 0.01 | 0.01 | 1    | 1    |


To model DEC0, you can use the wrapper function like the example below (assuming that you have downloaded all the required files).

```R
run_dec(treefile = "FelidaeTimes_pruned_no_F_catus.tre", geofile = "Felidae_biogeo_data.txt",
        multiplierfile = "DEC0_dispersal_multipliers_no_diff.txt", maxrange = 2, timesfile = "time_periods.txt",
        resultsfile = "DEC0_result.Rdata", section=FALSE)
```

After successfully running the code above, you should see a new file saved to your working directory named `DEC0_result.Rdata`. 

Then repeat this for DEC1.

```R
run_dec(treefile = "FelidaeTimes_pruned_no_F_catus.tre", geofile = "Felidae_biogeo_data.txt",
        multiplierfile = "DEC1_dispersal_multipliers_no_time.txt", maxrange = 2, timesfile = "time_periods.txt",
        resultsfile = "DEC1_result.Rdata", section=FALSE)
```

#### DEC2 - a model with time stratification

As you know, the position of continents has changed over time. Felidae is a relatively young group (we inferred the root of the tree being ca. 16 MYA) and the position of continets hasn't changed so much since then, but there are some changes and we'll explore this by adding time stratification into our models and providing different dispersal multiplier matrices for the different time periods.

We will divide the time into two periods: one between 20 and 12 mya and one between 12 and 0 mya. For the most recent time period, we will use the same dispersal multipliers as in DEC1. For the older period (20-12 mya), there are several differences in the dispersal multipliers. Download the DEC2 dispersal multiplier matrix from [here](../../Data/Day4/Biogeography/) and have a look at the dispersal rate multipliers for the older period. Can you spot the differences in values between some regions? Which regions are these and why are there such differences?

Let us run this final DEC model using the wrapper function.

We need to specify `section=TRUE` to state that the tree will be sectioned according to our time periods.

```R
run_dec(treefile = "FelidaeTimes_pruned_no_F_catus.tre", geofile = "Felidae_biogeo_data.txt",
        multiplierfile = "DEC2_dispersal_multipliers_time_strat.txt", maxrange = 2, timesfile = "time_periods.txt",
        resultsfile = "DEC2_result.Rdata", section=TRUE)
```


### Summarizing and visualizing the results

When you run the code block has `source()` above, you should get two additional functions apart from the wrapper `run_dec()` function. These are `plot_models()` and `get_table()`. We will use these functions to compare our models. But first, we need to load the results from the results files. Assuming that you specified `resultsfile="DEC0_result.Rdata"`, `resultsfile="DEC1_result.Rdata"`, etc for DEC0, DEC1:

```R
load("DEC0_result.Rdata")
res_DEC0 = res

load("DEC1_result.Rdata")
res_DEC1 = res

load("DEC2_result.Rdata")
res_DEC2 = res

```

Try this:

```R
results = list(DEC0=res_DEC0, DEC1=res_DEC1, DEC2=res_DEC2)
get_table(results) 
```

you should see a table like below.

```
              d          e        lnL
DEC0 0.03489101 0.01608103 -100.33974
DEC1 0.10078322 0.01672941  -89.98125
DEC2 0.10142983 0.01658964  -89.76118

```

Now, we can statistically compare DEC1 and DEC2 models to our DEC0 model. For this comparison we will use the AIC values. AIC stands for Akaike Information Criterio. You can read more about it [here](https://www.scribbr.com/statistics/akaike-information-criterion/#:~:text=To%20compare%20models%20using%20AIC%2C%20you%20need%20to%20calculate%20the,calculating%20log%2Dlikelihood%20is%20complicated!). Basically, if a model has an AIC value that is more than two points lower, than that model is considered better.


```R
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(res_DEC0)
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(res_DEC1)

stats = AICstats_2models(LnL_1, LnL_2, 2, 2)

stats$AIC1
stats$AIC2
```

```R
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(res_DEC0)
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(res_DEC2)

stats$AIC1
stats$AIC2
```


Can you change the code to compare DEC1 and DEC2? Which model is the best fit to our data?

<!--
```R
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(res_DEC2)
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(res_DEC0)

stats = AICstats_2models(LnL_1, LnL_2, 2, 2)
stats
```
-->


To plot the ancestral area resonctructions, we can use `plot_models()`:

```R
plot_models(res_DEC0, "DEC0")
```


You can also save this to a PDF file by:

```R
pdf("DEC0_plot.pdf", 10, 10)
plot_models(res_DEC0, "DEC0")
dev.off()
```

Plot each of our four models with the function `plot_models()` and save them to their own pdf files.

### Discussion

Hopefully you managed to run the analyses and you got your tree file with the ancestral areas inferred. If not, you can get these files from [here](../../Data/Day4/Biogeography/Results). There are two kinds of plots: i) where the most likely area or area combination is shown at each node and ii) with pie charts showing the probability of the most likely states.

Questions:

1. Which model is the best fit to the data?
2. Are the ancestral area reconstructions the same across the different models?
3. Can you identify some range expansions (i.e., dispersal events)?
4. And can you identify some range extinctions, i.e. where the ancestor was widely distributed but it then went extinct in one of the ranges (range contraction)?
5. Finally, can you identify an inferred vicariance event?
