# Todd Sechser and Joshua Alley
# Clean data from YouGov Experiment


# load data
yougov.data <- read_sav("data/UNVA0002_output.sav")
class(yougov.data) <- "data.frame"

# load vignette content 
vignettes.yougov <- read.csv("data/vignette-content-yougov.csv")

# add vignette info
yougov.data <- left_join(yougov.data, vignettes.yougov)

# recode function to switch agreement order
recode.yg <- function(x){
  x <- 6 - x # only works 1-5 Likert
}

yougov.data$Q5[yougov.data$Q5 == 8] <- NA

# clean data 
twdata <- yougov.data %>%
      rename(
        peace.str = Q1s,
        force.worse = Q6s,
        war.unf = Q3s,
        intl.trust = Q4a,
        isolation = Q7,
        us.sup = Q4b,
        us.shame = Q5
      ) %>% 
      mutate(
        # control variables
        republican = ifelse(pid3 == 2, 1, 0),
        democrat = ifelse(pid3 == 1, 1, 0),
        age = 2021 - birthyr,
        female = ifelse(gender == 2, 1, 0),
        
        # dispositional vars
        isolation = recode.yg(isolation),
        isolation.dum = ifelse(isolation > 3, 1, 0),
        mil.assert = (recode.yg(peace.str) +
                     recode.yg(war.unf) - force.worse) +
          3, # scale 0-12
        hawk = ifelse(mil.assert > 6, 1, 0),
        intl.trust.dum = ifelse(intl.trust == 1, 1, 0),
        natl.chauv = us.shame + (4 - us.sup),
        
        # rescale dispositional variables
        mil.assert.rs = mil.assert / max(mil.assert, na.rm = TRUE) * 4,
        isolation.rs = isolation / max(isolation, na.rm = TRUE) * 4,
        natl.chauv.rs = natl.chauv / max(natl.chauv, na.rm = TRUE) * 4,
        
        # casualties dummies
        cas.dum = ifelse(casualties > 0, 1, 0),
        cas.9 = ifelse(casualties == 9, 1, 0),
        cas.50 = ifelse(casualties == 50, 1, 0),
        cas.250 = ifelse(casualties == 250, 1, 0),
        
        # outcome dummy
        use.force = ifelse(Q123 == 1, 1, 0),
        
        # mechanism questions: swap scale 
        reputation = recode.yg(q129s),
        punishment = recode.yg(q130s), # typo issue here
        honor = recode.yg(q131s),
        interests = recode.yg(q132s),
        precedent = recode.yg(q141s),
        
        # military experience
        military.prox = ifelse(
          milstat_1 == 1, 4,
           ifelse(
             milstat_2 == 1, 3,
            ifelse(
              milstat_3 == 1, 2,
              ifelse(
                milstat_4 == 1, 1, 0)
              )))
      ) # end mutate

# edit outcome measures on specific policy
twdata$Q124[twdata$Q124 >= 8] <- NA
twdata$Q125[twdata$Q125 >= 8] <- NA

# factors on policy outcome
twdata$nonmil.res <- factor(recode(as.numeric(twdata$Q124),
                            `1` = "Do nothing",
                            `2` = "Lodge a diplomatic protest",
                            `3` = "Impose economic sanctions on the attacker"),
                            ordered = TRUE,
                            levels = c("Do nothing", 
                                       "Lodge a diplomatic protest",
                                       "Impose economic sanctions on the attacker"
                            ))
table(twdata$nonmil.res)
# militart res 
twdata$mil.res <- factor(recode(as.numeric(twdata$Q125),
                                   `1` = "Launch airstrikes against Country A's forces",
                                   `2` = "Send ground troops to attack Country A's forces",
                                   `3` = "Use nuclear weapons against Country A's forces"),
                            ordered = TRUE,
                         )
table(twdata$mil.res)

# full response
twdata <- twdata %>% unite("full.res",
                            mil.res:nonmil.res, 
                            remove = FALSE,
                            na.rm = TRUE)
table(twdata$full.res)
twdata$full.res <- factor(twdata$full.res,
                    ordered = TRUE,
                    levels = c("Do nothing", 
                     "Lodge a diplomatic protest",
                     "Impose economic sanctions on the attacker",
                     "Launch airstrikes against Country A's forces",
                     "Send ground troops to attack Country A's forces",
                     "Use nuclear weapons against Country A's forces")
                          )


# full treatment content
# treatment content
twdata <- twdata %>%
  mutate(
    cas.text = paste("Casualties", casualties,  "|"),
    all.text = ifelse(alliance == 1, 
                      "Alliance |", ""),
    nuke.text = ifelse(nuclear == 1, 
                       "Nuclear |", ""),
    democ.text = ifelse(democracy == 1, 
                        "Democracy", ""),
    treat.text = paste(cas.text, all.text, 
                       nuke.text, democ.text,
                        sep = " ")
  ) 
table(twdata$treat.text)

treat.unique <- as.data.frame(unique(twdata$treat.text))
colnames(treat.unique) <- "treat.text"
treat.unique <- treat.unique %>%
                 mutate(
                   n = nchar(treat.text)
                 )
treat.unique <- treat.unique[with(treat.unique, 
                  order(n, treat.text)), ]

twdata$treat.text <- factor(twdata$treat.text,
                            ordered = TRUE,
                            levels = str_sort(unique(twdata$treat.text),
                            numeric = TRUE))
                            # levels = unique(twdata$treat.text)[
                            #      order(nchar(unique(twdata$treat.text)),
                            #         unique(twdata$treat.text))])



# manipulation checks
twdata <- twdata %>%
  mutate(
    cas.manip = ifelse((Q135 == 1 & casualties == 0) |
                       (Q135 == 2 & casualties == 9) |
                       (Q135 == 3 & casualties == 50) |
                       (Q135 == 4 & casualties == 250),
                       1, 0),
    cas.manip.dum = ifelse((Q135 > 1 & casualties > 0), 
                           1, 0), # any casualties
    nuke.manip = ifelse((Q136 == 1 & nuclear == 1) |
                          (Q136 == 2 & nuclear == 0),
                        1, 0),
    ally.manip = ifelse((Q137 == 1 & alliance == 1) |
                          (Q137 == 2 & alliance == 0),
                        1, 0),
    manip.cor = cas.manip + nuke.manip + ally.manip,
    # with cas dummy 
    manip.cor.dum = cas.manip.dum + nuke.manip + ally.manip
  )
table(twdata$manip.cor)
# % passed all MC
1967 / 3200



# filter to passed manipulations
twdata.pass <- filter(twdata, manip.cor == 3)


# clean july 2023 data 
july.23.data <- read.csv("data/Tripwires Experiment July 2023.csv")
glimpse(july.23.data)

# create policy choice labels
july.23.data <- july.23.data %>%
    mutate(
      response.txt = 
        factor(case_when(
          response == 1 ~ "Do\nNothing",
          response == 2 ~ "Diplomatic\nProtest",
          response == 3 ~ "Economic\nSanctions",
          response == 4 ~ "Cyber\nAttacks",
          response == 5 ~ "Airstrikes",
          response == 6 ~ "Ground\nTroops"),
          ordered = T,
          levels = c("Do\nNothing",
                     "Diplomatic\nProtest",
                     "Economic\nSanctions",
                     "Cyber\nAttacks",
                     "Airstrikes",
                     "Ground\nTroops")),
      country = case_when(
        china == 1 ~ "China & Taiwan",
        russia == 1 ~ "Russia & Estonia",
        .default = "Generic"
      )
    )




# load cleaned data
twdata <- read.csv("data/Tripwires Experiment - YouGov.csv") %>%
  mutate(
    # outcome dummy
    use.force = ifelse(q123 == "Yes", 1, 0),
    educ = as.numeric(factor(educ, ordered = TRUE,
                             levels = c(
                               "No HS", "High school graduate",
                               "Some college", "2-year", "4-year", 
                               "Post-grad"           
                             )))
  )

# add a couple other variables 
# factors on policy outcome
twdata$nonmil.res <- factor(twdata$q124,
                            ordered = TRUE,
                            levels = c("Do nothing", 
                                       "Lodge a diplomatic protest",
                                       "Impose economic sanctions on the attacker"
                            ))
table(twdata$nonmil.res)
# militart res 
twdata$mil.res <- factor(twdata$q125,
                         ordered = TRUE,
                         levels = c("Launch airstrikes against Country A's forces",
                                    "Send ground troops to attack Country A's forces",
                                    "Use nuclear weapons against Country A's forces")
)
table(twdata$mil.res)

# full response
twdata <- twdata %>% unite("full.res",
                           mil.res:nonmil.res, 
                           remove = FALSE,
                           na.rm = TRUE)
table(twdata$full.res)
twdata$full.res <- factor(twdata$full.res,
                          ordered = TRUE,
                          levels = c("Do nothing", 
                                     "Lodge a diplomatic protest",
                                     "Impose economic sanctions on the attacker",
                                     "Launch airstrikes against Country A's forces",
                                     "Send ground troops to attack Country A's forces",
                                     "Use nuclear weapons against Country A's forces")
)


# full treatment content
# treatment content
twdata <- twdata %>%
  mutate(
    cas.text = paste("Casualties", casualties,  "|"),
    all.text = ifelse(ally == 1, 
                      "Alliance |", ""),
    nuke.text = ifelse(nuclear == 1, 
                       "Nuclear |", ""),
    democ.text = ifelse(democracy == 1, 
                        "Democracy", ""),
    treat.text = paste(cas.text, all.text, 
                       nuke.text, democ.text,
                       sep = " "),
    treat.text.het = paste(all.text, 
                           nuke.text, democ.text,
                           sep = " "),
  ) 
table(twdata$treat.text)

twdata$treat.text <- factor(twdata$treat.text,
                            ordered = TRUE,
                            levels = str_sort(unique(twdata$treat.text),
                                              numeric = TRUE))

# passed all checks
twdata.pass <- filter(twdata, mc_all == 1)

