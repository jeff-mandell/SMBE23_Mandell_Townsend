# See project repository (https://github.com/jeff-mandell/SMBE23_Mandell_Townsend) for instructions
# to install the exact versions of cancereffectsizeR and ces.refset.hg38 used in this analysis.
# Alternatively, future versions of these packages will likely work, too (again, see repo).

library(cancereffectsizeR)
library(data.table)
library(ggplot2)
library(ggrepel)
library(ggstatsplot)
library(irr) # for calculation of intraclass correlation


# Get TCGA SKCM MAF data. Naming it to make clear that we're including metastatic samples,
# but excluding patients that have multiple sequenced samples.  
mela_maf = 'TCGA-SKCM_with_met_exclude_multisample.maf.gz'
if (! file.exists(mela_maf)) {
  tmp_maf = paste0(tempfile(), '.maf')
  get_TCGA_project_MAF(project = 'TCGA-SKCM', filename = tmp_maf, exclude_TCGA_nonprimary = FALSE)
  
  # 2 patients have 2 samples each. For simplicity, we'll these patients.
  maf_to_edit = fread(tmp_maf)
  stopifnot(maf_to_edit[multisample_patient == T, 
                        uniqueN(Tumor_Sample_Barcode) == 4 && uniqueN(Unique_Patient_Identifier) == 2])
  maf_to_edit = maf_to_edit[multisample_patient == FALSE]
  fwrite(maf_to_edit, mela_maf, sep = "\t")
  unlink(tmp_maf)
}

# Load SKCM data into CESAnalysis
mela = preload_maf(maf = mela_maf, refset = 'ces.refset.hg38')
cesa = CESAnalysis('ces.refset.hg38')
cesa = load_maf(cesa = cesa, maf = mela) 


## Implement dNdScv's indel model, just for DBS

# Organize covariates data
gr_genes = ces.refset.hg38$gr_genes
cv = ces.refset.hg38$covariates$SKCM$rotation
pid_to_gene = unique(as.data.table(GenomicRanges::mcols(gr_genes)))
setnames(pid_to_gene, 'names', 'pid')
current_cv_genes = rownames(cv)
present_genes = intersect(pid_to_gene$gene, rownames(cv)) # some may be missing, which gets dealt with later
present_pid_to_gene = pid_to_gene[present_genes, on = "gene"]
cv = cv[present_pid_to_gene$gene, ]
rownames(cv) = present_pid_to_gene$pid

# Load "known_cancergenes" from dNdScv
data(list=sprintf("cancergenes_cgc81"), package="dndscv")
cancer_genes = as.data.table(ces.refset.hg38$gr_genes)[known_cancergenes, unique(names), on = 'gene', nomatch = NULL]

# Annotate DBS with overlapping protein_id
dbs_gr = makeGRangesFromDataFrame(cesa$maf[variant_type == 'dbs', .(seqnames = Chromosome,
                                                                    start = Start_Position, end = Start_Position + 1, variant_id)],
                                  keep.extra.columns = TRUE)
overlaps = findOverlaps(dbs_gr, gr_genes, type = 'within')
dbs_by_pid = cbind(as.data.table(dbs_gr[queryHits(overlaps)]), as.data.table(gr_genes[subjectHits(overlaps)]))
dbs_by_pid = dbs_by_pid[, .(pid_prev = .N, gene_name = gene[1]), by = 'names']
setnames(dbs_by_pid, 'names', 'pid')

# Add in rows for pid that don't have any DBS.
pid_with_none = setdiff(names(ces.refset.hg38$RefCDS), dbs_by_pid$pid)
for_merge = data.table(pid = pid_with_none, pid_prev = 0, gene_name = sapply(ces.refset.hg38$RefCDS[pid_with_none], '[[', 'real_gene_name'))
dbs_by_pid = rbind(dbs_by_pid, for_merge)

# Assemble data for model: Gene name, number of DBS appearing in the CDS (coding sequence; called n_indused to match
# dNdScv's code), length of the CDS, and covariates, cancer gene status. ( Could use max 1 DBS (indel) per sample per
# gene, as dNdScv default does, but not doing that for now.)
for_model = dbs_by_pid[, .(gene_name, n_indused = pid_prev, 
                           cds_length = sapply(ces.refset.hg38$RefCDS[pid], '[[', 'CDS_length'), pid)]
for_model[, is_known_cancer := gene_name %in% known_cancergenes] # known_cancergenes from dNdScv

unif_site_rate = for_model[is_known_cancer == FALSE, sum(n_indused)/sum(cds_length)]
for_model[, exp_unif := cds_length * unif_site_rate]

## Leave out pid with no covariates.
for_model = for_model[pid %in% rownames(cv)]
covs = cv[for_model$pid,]

# Cancer genes are excluded from model fitting. All genes will then get predictions from the fitted model.
# I assume nbrdf is negative binomial regression data frame.
nbrdf_all = cbind(for_model[, .(n_indused, exp_unif)], covs) 
nbrdf = nbrdf_all[! for_model$is_known_cancer]

nb = MASS::glm.nb(n_indused ~ offset(log(exp_unif)) + . , data = nbrdf)

# average mutation rates per DBS site within each pid
for_model[, all_rates := exp(predict(nb,nbrdf_all))]
num_samples = uniqueN(cesa$maf$Unique_Patient_Identifier)

# 9 DBS are possible per CDS site (each site can be first base of a DBS, and there are 9 DBS for any dinucleotide).
for_model[, final_rate := all_rates / (cds_length * 9) / num_samples]

# Get relative trinuc rates across the cohort. Add pseudocounts becasue a couple contexts that don't appear.
dbs = copy(cesa@mutations$dbs)
codon_change = copy(cesa@mutations$dbs_codon_change)

# This was supposed to be exported by cancereffectsizeR; will get fixed in a future update. For now,
# extracting the internal data that lists the 78 COSMIC DBS classes.
cosmic_dbs_classes = cancereffectsizeR:::cosmic_dbs_classes
cosmic_class_count = table(factor(dbs$cosmic_dbs_class, levels = cosmic_dbs_classes))
if(0 %in% cosmic_class_count) {
  cosmic_class_count = cosmic_class_count + 1
}
dbs_prop = cosmic_class_count / sum(cosmic_class_count)

# Likewise, pulling internal tabulation of dinucleotide contexts per gene from refset.
dbs_exposure = cancereffectsizeR:::get_ref_data(cesa, 'cds_dbs_exposure')

# Get context-specific mutation rates for each gene.
gene_site_rates = mapply(
  function(pid, avg_rate) {
    (avg_rate * dbs_prop) / sum(dbs_prop * dbs_exposure[[pid]])
  },
  for_model$pid, for_model$final_rate, SIMPLIFY = FALSE)
names(gene_site_rates) = for_model$pid

# Likelihood function for clonal selection model
dbs_lik = function(rate, num_with, num_without) {
  fn = function(gamma) {
    gamma = unname(gamma) # math faster on unnamed vectors
    sum_log_lik = 0
    if (num_without > 0) {
      sum_log_lik = -1 * gamma * num_without * rate
    }
    if (num_with > 0) {
      sum_log_lik = sum_log_lik + sum(log(1 - exp(-1 * gamma * rep.int(rate, num_with))))
    }
    
    # convert to negative loglikelihood and return
    return(-1 * sum_log_lik)
  }
  # Set default values for gamma (SI), which ces_variant will use to set starting value of optimization
  formals(fn)[["gamma"]] = 1
  bbmle::parnames(fn) = "selection_intensity"
  return(fn)
}

## Collect DBS single-codon changes and run through selection inference
# Restrict to PIDs with covariates, then take one PID per gene to avoid redundancy.
easy_codon_change = codon_change[pid %in% names(gene_site_rates)]
easy_codon_change = easy_codon_change[, .(.N, gene = gene[1]), by = 'pid'][order(-N)][!(duplicated(gene))]
easy_codon_change = codon_change[pid %in% easy_codon_change$pid] # 5,152 distinct records remain

# Get specific rates for DBS codon changes (summing multiple DBS where necessary)
dbs_key = copy(cesa@mutations$aac_dbs_key)[dbs_aac_id %in% easy_codon_change$dbs_aac_id]
dbs_key[dbs, cosmic_dbs_class := cosmic_dbs_class, on = 'dbs_id']
dbs_key[easy_codon_change, pid := pid, on = 'dbs_aac_id']

dbs_key[, rate := mapply(function(x, y) gene_site_rates[[x]][y], dbs_key$pid, dbs_key$cosmic_dbs_class)]
dbs_count = cesa$maf[variant_type == 'dbs', .N, by = 'variant_id']
dbs_key[dbs_count, N := N, on = c(dbs_id = 'variant_id')]
dbs_key[is.na(N), N := 0]
info_by_id = dbs_key[, .(N = sum(N), rate = sum(rate)), by = 'dbs_aac_id']

# Run the selection inference. May take a minute or two.
dbs_res = rbindlist(mapply(
  function(curr_dbs, num_with, only_rate) {
    num_without = num_samples - num_with
    fn = dbs_lik(only_rate, num_with, num_without)
    par_init = formals(fn)[[1]]
    names(par_init) = bbmle::parnames(fn)
    # find optimized selection intensities
    # the selection intensity for any stage that has 0 variants will be on the lower boundary; will muffle the associated warning
    withCallingHandlers(
      {
        fit = bbmle::mle2(fn, method="L-BFGS-B", start = par_init, vecpar = T, lower=1e-3, upper=1e9)
      },
      warning = function(w) {
        if (startsWith(conditionMessage(w), "some parameters are on the boundary")) {
          invokeRestart("muffleWarning")
        }
        if (grepl(x = conditionMessage(w), pattern = "convergence failure")) {
          # a little dangerous to muffle, but so far these warnings are
          # quite rare and have been harmless
          invokeRestart("muffleWarning") 
        }
      }
    )
    
    selection_intensity = bbmle::coef(fit)
    loglikelihood = as.numeric(bbmle::logLik(fit))
    return(list(variant_id = curr_dbs, selection_intensity = selection_intensity, ll = loglikelihood,
                num_with = num_with, rate = only_rate))
    
  }, info_by_id$dbs_aac_id, info_by_id$N, info_by_id$rate, SIMPLIFY = FALSE))

# Annotate output.
dbs_res[codon_change, c("gene", "aachange") := .(gene, aachange), on = c(variant_id = 'dbs_aac_id')]
dbs_res[, sbs_aac := paste0(gene, '_', aachange)]

# See what DBS codon change variants also have SBS AACs in the data set.
can_check = intersect(dbs_res$sbs_aac, cesa$variants$variant_name) # 171 can be checked.

# Run steps for SBS selection inference. Only running the variants that we need to compare to DBS.
cesa = gene_mutation_rates(cesa, covariates = 'SKCM')
signature_exclusions = suggest_cosmic_signature_exclusions('SKCM', T, T)
cesa = trinuc_mutation_rates(cesa, 'COSMIC_v3.2', signature_exclusions = signature_exclusions)
cesa = ces_variant(cesa, variants = cesa$variants[can_check, on = 'variant_name'])


# Merge SBS and DBS inferences and compare
dbs_subset = dbs_res[can_check, on = 'sbs_aac']
merged = merge.data.table(dbs_subset, cesa$selection$selection.1, by.x = 'sbs_aac', by.y = 'variant_name',
                          suffixes = c('.dbs', '.snv'))

clean = merged[, .(variant_name = sbs_aac, dbs_id = variant_id.dbs, N_dbs = num_with, N_snv = included_with_variant, si_dbs = selection_intensity.dbs, 
                   si_snv = selection_intensity.snv)]

# Prepare Figure 2. (Figure 1 was reprinted from another publication.)
# We'll label recurrent DBS and anything getting high enough effect to have space for a label.
clean[, color := 'gray70']
clean[N_dbs > 1, color := 'black']
clean[, label := variant_name]
clean[, gene := gsub('_.*', '', variant_name)]
clean[, is_kc := gene %in% known_cancergenes]
clean[si_snv < 1e3 & si_dbs < 1e4 & N_dbs == 1, label := '']
clean[label != '', label := gsub('_', ' ', label)]

label_fn = function(x) format(x, big.mark = ",", scientific = F)
f2 = ggscatterstats(data = clean,
                    point.args = list(size = 2.5, stroke = 0, color = clean$color), 
                    x = si_snv, y = si_dbs, type = 'parametric', marginal = F, bf.message = F) + 
  scale_x_log10(labels = label_fn, breaks = 10^(1:6)) + 
  scale_y_log10(labels = label_fn, breaks = 10^(2:6), limits = c(99, NA)) +
  geom_text_repel(aes(label = label), size = 3, nudge_x = .3, nudge_y = .04, seed = 919) +
  xlab('SBS-induced amino-acid substitution') + ylab('DBS-induced amino-acid substitution') + 
  ggtitle('Effects of amino-acid changes by nucleotide substitution type') + 
  theme_classic() + 
  theme(plot.title = element_text(size = 16), axis.title = element_text(size = 16))

# Picking and choosing which correlation test outputs to include in subtitle
f2$labels$subtitle = as.call(as.list(f2$labels$subtitle)[c(1, 4:6, 3)])
ggsave(plot = f2, filename = 'figure2.png', width = 8, height = 8, units = 'in')



## For comparison, let's look at cases where distinct SBS induce the same codon change as each other.
# Find AAC that have multiple distinct SNVs in the data
all_aac = cesa$variants[variant_type == 'aac', unique(variant_id)]
multi_snv_aac = cesa@mutations$aac_snv_key[all_aac, .(snv_id, multi = length(snv_id) > 1), 
                                            on = 'aac_id', by = 'aac_id'][multi == TRUE, .(aac_id, snv_id)]
multi_snv_aac[, in_maf := snv_id %in% cesa$maf$variant_id]

# There are no AACs that are inflicted by more than two distinct SBS in the data set.
triple_check = multi_snv_aac[, .(triple = sum(in_maf) > 2), by = 'aac_id']
stopifnot(all(triple_check$triple == F))

to_use = multi_snv_aac[, .(to_use = sum(in_maf) == 2), by = 'aac_id'][to_use == TRUE]
snv_aac_to_use = multi_snv_aac[in_maf == TRUE & aac_id %in% to_use$aac_id]

# Verify that we have gotten exactly 2 SBS per AAC
stopifnot(snv_aac_to_use[, .N, by = 'aac_id'][, unique(N)] == 2)

sbs_aac_variants = select_variants(cesa, variant_ids = snv_aac_to_use$snv_id)
cesa = ces_variant(cesa, variants = sbs_aac_variants, run_name = 'sbs_aac')
sbs_aac_output = cesa$selection$sbs_aac
sbs_aac_output[snv_aac_to_use, aac_id := aac_id, on = c(variant_id = 'snv_id')]

# Verify that every codon (that is, pid/aa_pos) is only used in one pair of variants
sbs_aac_output[cesa$variants, aa_pos := aa_pos, on = c(aac_id = 'variant_id')]
sbs_aac_output[, pid := gsub('_.*', '', aac_id)]
stopifnot(uniqueN(sbs_aac_output[, .(pid, aa_pos)]) == sbs_aac_output[, .N/2])


# For intraclass correlation, doesn't matter how variants are divided into grp1/grp2
grp1 = sbs_aac_output[, .SD[1], by = 'aac_id']
grp2 = sbs_aac_output[, .SD[2], by = 'aac_id']
stopifnot(identical(grp1$aac_id, grp2$aac_id))

# Lots of these "codon changes" are actually synonymous.
# The synonymous cases provide a sort of negative control.
cesa$variants[sbs_aac_output$aac_id, on = 'variant_id'][, table(aa_ref == aa_alt)]
# FALSE  TRUE 
# 224   268

# So syn/nonsyn pairs should be 134 and 112

syn = cesa$variants[sbs_aac_output$aac_id, on = 'variant_id'][aa_ref == aa_alt, unique(variant_id)]
just_syn_grp1 = grp1[syn, on = 'aac_id']
just_syn_grp2 = grp2[syn, on = 'aac_id']
not_syn_grp1 = grp1[! syn, on = 'aac_id']
not_syn_grp2 = grp2[! syn, on = 'aac_id']
stopifnot(just_syn_grp1[, .N] == 134,
          just_syn_grp2[, .N] == 134,
          not_syn_grp1[, .N] == 112,
          not_syn_grp2[, .N] == 112)

nonsyn_sbs = icc(data.table(not_syn_grp1$selection_intensity, not_syn_grp2$selection_intensity), model = 'oneway')
syn_sbs = icc(data.table(just_syn_grp1$selection_intensity, just_syn_grp2$selection_intensity), model = 'oneway')
dbs_sbs = icc(data.table(clean[, .(si_dbs, si_snv)]), model = 'oneway')

icc_table = rbindlist(lapply(list(`DBS/SBS same codon change` = dbs_sbs, 
                                  `SBS/SBS same codon change` = nonsyn_sbs,
                                  `SBS same codon synonymous` = syn_sbs), 
                             function(x) {
                               x[c('subjects', 'value', 'lbound', 'ubound')]
                               }), idcol = 'group')
icc_table = icc_table[, .(group, `number of sites` = subjects, 
                          `intra-class correlation` = round(value, 2),
                          `95% CI` = paste0('(', round(lbound, 2), ', ', round(ubound, 2), ')'))]
fwrite(icc_table, 'dbs_sbs_icc_for_table1.txt', sep = "\t")

# Optionally, save a copy of the CESAnalysis for future work (reload with load_cesa()).
# save_cesa(cesa = cesa, 'for_dbs.rds')
# 
