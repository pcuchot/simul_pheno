# -----------------------------------------------------------------------------
# title : 1_simulate_data
# author : 
# date : 
# note : 
# -----------------------------------------------------------------------------


# Packages ----------------------------------------------------------------


# -------------------------------------------------------------------------

# number of sites
n_site <- 100
# number of year (per site?)
n_year <- 5
# number of capture session per site
n_session <- 4
# number of breeding birds per site (catchable adults)
n_breed <- 3
# number of breeding attempts
n_bredd_att <-  1
# mean breeding date per site (vector)
m_ld_site <- runif(n_site,100,200)

# per year_site, create a cumulative vector (one value per day)


