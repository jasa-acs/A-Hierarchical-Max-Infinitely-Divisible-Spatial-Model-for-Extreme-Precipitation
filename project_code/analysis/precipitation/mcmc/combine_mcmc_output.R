library(stablemix)

#   -----------------------------------------------------------------------
# Description: If MCMC sampler is run in chunks (to deal with e.g. max wall times)
#			   this script will combine the output chunks and store the combined data.
#
# Output: 
#     - Combined datasets for each of the four models are stored in 
#		lnorm_maxid.Rdata
#		lnorm_maxstable.Rdata
#		fixed_maxid.Rdata
#		fixed_maxstable.Rdata
#   -----------------------------------------------------------------------


# Lnorm -------------------------------------------------------------------
# Maxid
out <- lnorm_combine_mcmc(file.path(samples.dir, "lnorm/maxid"), 
                          burnin = 500)
save(out, file = file.path(samples.dir, "lnorm_maxid.Rdata"))
# Maxstable
out <- lnorm_combine_mcmc(file.path(samples.dir, "lnorm/maxstable"), 
                          burnin = 500)
save(out, file = file.path(samples.dir, "lnorm_maxstable.Rdata"))

# Fixed -------------------------------------------------------------------
# Maxid
out <- fixed_combine_mcmc(file.path(samples.dir, "fixed/maxid"), 
                          burnin = 500)
save(out, file = file.path(samples.dir, "fixed_maxid.Rdata"))
# Maxstable
out <- fixed_combine_mcmc(file.path(samples.dir, "fixed/maxstable"),
                          burnin = 500)
save(out, file = file.path(samples.dir, "fixed_maxstable.Rdata"))
