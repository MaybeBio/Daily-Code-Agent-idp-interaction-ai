library("dplyr")
library("ggplot2")
library("ggh4x")


#Load in data into data frames
ridaodat_whole <- read.table("/file/path/MDP_MODE/results.txt", header=TRUE) # proteins in entire secretome
ridaodat_blood <- read.table("/file/path/MDP_MODE/results.txt", header=TRUE) # only blood proteins in secretome
ridaodat_lim <- read.table("/file/path/MDP_MODE/results.txt", header=TRUE) # proteins of a limited secretome representing the kidney
ridaodat_hum <- read.table("/file/path/MDP_MODE/results.txt", header=TRUE) # entire human proteome
ridaodat_canon <- read.table("/file/path/MDP_MODE/results.txt", header=TRUE) # canonical kidney stone related proteins
ridaodat_restnet <- read.table("/file/path/MDP_MODE/results.txt", header=TRUE) # restricted expanded network of kidney stone related proteins
ridaodat_unrestnet <- read.table("/file/path/MDP_MODE/results.txt", header=TRUE) # unrestricted expanded network of kidney stone related proteins

#Clean out unused columns (leave CH CDF columns)
cleaned_whole <- ridaodat_whole %>%
  select(-SEQ_ID., -PER.VLXT.., -RNK.VLXT.., -PER.VSL2B.., -RNK.VSL2B.., -PER.VL3.., -RNK.VL3.., -PER.IUP_S.., -RNK.IUP_S.., -PER.IUP_L.., -RNK.IUP_L.., -PER.PFIT.., -RNK.PFIT.., -PER.MDP.., -RNK.MDP.., -len.seq.)
cleaned_lim <- ridaodat_lim %>%
  select(-SEQ_ID., -PER.VLXT.., -RNK.VLXT.., -PER.VSL2B.., -RNK.VSL2B.., -PER.VL3.., -RNK.VL3.., -PER.IUP_S.., -RNK.IUP_S.., -PER.IUP_L.., -RNK.IUP_L.., -PER.PFIT.., -RNK.PFIT.., -PER.MDP.., -RNK.MDP.., -len.seq.)
cleaned_blood <- ridaodat_blood %>%
  select(-SEQ_ID., -PER.VLXT.., -RNK.VLXT.., -PER.VSL2B.., -RNK.VSL2B.., -PER.VL3.., -RNK.VL3.., -PER.IUP_S.., -RNK.IUP_S.., -PER.IUP_L.., -RNK.IUP_L.., -PER.PFIT.., -RNK.PFIT.., -PER.MDP.., -RNK.MDP.., -len.seq.)
cleaned_hum <- ridaodat_hum %>%
  select(-SEQ_ID., -PER.VLXT.., -RNK.VLXT.., -PER.VSL2B.., -RNK.VSL2B.., -PER.VL3.., -RNK.VL3.., -PER.IUP_S.., -RNK.IUP_S.., -PER.IUP_L.., -RNK.IUP_L.., -PER.PFIT.., -RNK.PFIT.., -PER.MDP.., -RNK.MDP..)
cleaned_canon <- ridaodat_canon %>%
  select(-SEQ_ID., -PER.VLXT.., -RNK.VLXT.., -PER.VSL2B.., -RNK.VSL2B.., -PER.VL3.., -RNK.VL3.., -PER.IUP_S.., -RNK.IUP_S.., -PER.IUP_L.., -RNK.IUP_L.., -PER.PFIT.., -RNK.PFIT.., -PER.MDP.., -RNK.MDP.., -len.seq.)
cleaned_restnet <- ridaodat_restnet %>%
  select(-SEQ_ID., -PER.VLXT.., -RNK.VLXT.., -PER.VSL2B.., -RNK.VSL2B.., -PER.VL3.., -RNK.VL3.., -PER.IUP_S.., -RNK.IUP_S.., -PER.IUP_L.., -RNK.IUP_L.., -PER.PFIT.., -RNK.PFIT.., -PER.MDP.., -RNK.MDP.., -len.seq.)
cleaned_unrestnet <- ridaodat_unrestnet %>%
  select(-SEQ_ID., -PER.VLXT.., -RNK.VLXT.., -PER.VSL2B.., -RNK.VSL2B.., -PER.VL3.., -RNK.VL3.., -PER.IUP_S.., -RNK.IUP_S.., -PER.IUP_L.., -RNK.IUP_L.., -PER.PFIT.., -RNK.PFIT.., -PER.MDP.., -RNK.MDP.., -len.seq.)

#Remove commas at end of values
cleaned_whole <- cleaned_whole %>%
  mutate(across(c(dCDF., dCH.), ~ gsub(",", "", .)))
cleaned_lim <- cleaned_lim %>%
  mutate(across(c(dCDF., dCH.), ~ gsub(",", "", .)))
cleaned_blood <- cleaned_blood %>%
  mutate(across(c(dCDF., dCH.), ~ gsub(",", "", .)))
cleaned_hum <- cleaned_hum %>%
  mutate(across(c(dCDF., dCH.), ~ gsub(",", "", .)))
cleaned_canon <- cleaned_canon %>%
  mutate(across(c(dCDF., dCH.), ~ gsub(",", "", .)))
cleaned_restnet <- cleaned_restnet %>%
  mutate(across(c(dCDF., dCH.), ~ gsub(",", "", .)))
cleaned_unrestnet <- cleaned_unrestnet %>%
  mutate(across(c(dCDF., dCH.), ~ gsub(",", "", .)))

#Change chr to num for dCH and dCDF columns 
cleaned_whole$dCDF. <- as.numeric(cleaned_whole$dCDF.)
cleaned_whole$dCH. <- as.numeric(cleaned_whole$dCH.)
cleaned_lim$dCDF. <- as.numeric(cleaned_lim$dCDF.)
cleaned_lim$dCH. <- as.numeric(cleaned_lim$dCH.)
cleaned_blood$dCDF. <- as.numeric(cleaned_blood$dCDF.)
cleaned_blood$dCH. <- as.numeric(cleaned_blood$dCH.)
cleaned_hum$dCDF. <- as.numeric(cleaned_hum$dCDF.)
cleaned_hum$dCH. <- as.numeric(cleaned_hum$dCH.)
cleaned_canon$dCDF. <- as.numeric(cleaned_canon$dCDF.)
cleaned_canon$dCH. <- as.numeric(cleaned_canon$dCH.)
cleaned_restnet$dCDF. <- as.numeric(cleaned_restnet$dCDF.)
cleaned_restnet$dCH. <- as.numeric(cleaned_restnet$dCH.)
cleaned_unrestnet$dCDF. <- as.numeric(cleaned_unrestnet$dCDF.)
cleaned_unrestnet$dCH. <- as.numeric(cleaned_unrestnet$dCH.)


#Combine into one df
combined <- bind_rows(
  cleaned_whole %>% mutate(secretome = "Secretome: Whole"),
  cleaned_lim %>% mutate(secretome = "Secretome: Restricted"),
  cleaned_blood %>% mutate(secretome = "Secretome: Blood"),
  cleaned_hum %>% mutate(secretome = "Human Proteome"),
  cleaned_canon %>% mutate(secretome = "Stoneome: Canonical"),
  cleaned_restnet %>% mutate(secretome = "Stoneome: Restricted"),
  cleaned_unrestnet %>% mutate(secretome = "Stoneome: Unrestricted")
)

#Rearrange secretome plots to preferred order
combined$secretome <- factor(
  combined$secretome,
  levels = c(
    "Human Proteome",
    "Secretome: Whole",
    "Secretome: Restricted",
    "Secretome: Blood",
    "Stoneome: Canonical",
    "Stoneome: Restricted",
    "Stoneome: Unrestricted"
  )
)


#Create scatter plots
ggplot(combined, aes(dCDF., dCH., color = secretome)) +
  geom_point(size = 0.1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1)) +
  scale_color_manual(
    values = c(
      "Human Proteome"   = "#1b9e77",
      "Secretome: Whole" = "#7cb0f0", 
      "Secretome: Restricted"   = "#7570b3", 
      "Secretome: Blood" = "#d95f02",
      "Stoneome: Canonical" = "#1e3786",
      "Stoneome: Restricted" = "#ccccff", 
      "Stoneome: Unrestricted" = "#b4eeb4"
    )
  ) +
  facet_wrap(~ secretome, ncol = 4, axes = "all") +
  labs(title = "CH-CDF", x = expression(Delta*"CDF"), y = expression(Delta*"CH")) +
  theme_classic() +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "bold"), 
    panel.spacing = unit(1.2, "lines"), 
    plot.title = element_text(hjust = 0.5)
  )