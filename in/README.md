# Subject Lists

The lists are created according to [Jo's post](https://3.basecamp.com/3758557/buckets/3792852/messages/4106700214). For example, `subjects_wave12_good.txt` contains the subjects with high quality data in wave 1 and 2. `subjects_wave12_prob.txt` is removed.

For `code/timeseries/null_2rpm/run_deconvolve.sh`, we will use the lists `subjects_wave12_all.txt` (30 subjects, 20 good, 6 OK, 4 bad), `subjects_wave23_all.txt` (1 subject, good), and `subjects_wave13_all.txt` (1 subject, bad).

Then in `code/_constants.R` we only look at Stroop and exclude three subjects from further processing:

- DMCC5009144 has no wave 1 Stroop baseline
- 197449 has poor data in Stroop baseline
- 178243 wave 1 Stroop baseline run 2 ended a bit early

Also note that DMCC6418065 and DMCC6671683 have missing Stroop behavioral data.