library("dplyr")
library("ggplot2")

#Load in data into data frames (make sure to convert RIDAO results files to txt)
ridaodat_whole <- read.table("/file/path/MDP_MODE/results.txt", header=TRUE) # proteins in entire secretome
ridaodat_blood <- read.table("/file/path/MDP_MODE/results.txt", header=TRUE) # only blood proteins in secretome
ridaodat_lim <- read.table("/file/path/MDP_MODE/results.txt", header=TRUE) # proteins of a limited secretome representing the kidney
ridaodat_hum <- read.table("/file/path/MDP_MODE/results.txt", header=TRUE) # entire human proteome
ridaodat_canon <- read.table("/file/path/MDP_MODE/results.txt", header=TRUE) # canonical kidney stone related proteins
ridaodat_restnet <- read.table("/file/path/MDP_MODE/results.txt", header=TRUE) # restricted expanded network of kidney stone related proteins
ridaodat_unrestnet <- read.table("/file/path/MDP_MODE/results.txt", header=TRUE) # unrestricted expanded network of kidney stone related proteins

#Clean out unused columns (leave PER.VSL2B columns)
cleaned_whole <- ridaodat_whole %>%
  select(-SEQ_ID., -PER.VLXT.., -RNK.VLXT.., -RNK.VSL2B.., -PER.VL3.., -RNK.VL3.., -PER.IUP_S.., -RNK.IUP_S.., -PER.IUP_L.., -RNK.IUP_L.., -PER.PFIT.., -RNK.PFIT.., -PER.MDP.., -RNK.MDP.., -dCDF., -dCH., -len.seq.)
cleaned_lim <- ridaodat_lim %>%
  select(-SEQ_ID., -PER.VLXT.., -RNK.VLXT.., -RNK.VSL2B.., -PER.VL3.., -RNK.VL3.., -PER.IUP_S.., -RNK.IUP_S.., -PER.IUP_L.., -RNK.IUP_L.., -PER.PFIT.., -RNK.PFIT.., -PER.MDP.., -RNK.MDP.., -dCDF., -dCH., -len.seq.)
cleaned_blood <- ridaodat_blood %>%
  select(-SEQ_ID., -PER.VLXT.., -RNK.VLXT.., -RNK.VSL2B.., -PER.VL3.., -RNK.VL3.., -PER.IUP_S.., -RNK.IUP_S.., -PER.IUP_L.., -RNK.IUP_L.., -PER.PFIT.., -RNK.PFIT.., -PER.MDP.., -RNK.MDP.., -dCDF., -dCH., -len.seq.)
cleaned_hum <- ridaodat_hum %>%
  select(-SEQ_ID., -PER.VLXT.., -RNK.VLXT.., -RNK.VSL2B.., -PER.VL3.., -RNK.VL3.., -PER.IUP_S.., -RNK.IUP_S.., -PER.IUP_L.., -RNK.IUP_L.., -PER.PFIT.., -RNK.PFIT.., -PER.MDP.., -RNK.MDP.., -dCDF., -dCH.)
cleaned_canon <- ridaodat_canon %>%
  select(-SEQ_ID., -PER.VLXT.., -RNK.VLXT.., -RNK.VSL2B.., -PER.VL3.., -RNK.VL3.., -PER.IUP_S.., -RNK.IUP_S.., -PER.IUP_L.., -RNK.IUP_L.., -PER.PFIT.., -RNK.PFIT.., -PER.MDP.., -RNK.MDP.., -dCDF., -dCH., -len.seq.)
cleaned_restnet <- ridaodat_restnet %>%
  select(-SEQ_ID., -PER.VLXT.., -RNK.VLXT.., -RNK.VSL2B.., -PER.VL3.., -RNK.VL3.., -PER.IUP_S.., -RNK.IUP_S.., -PER.IUP_L.., -RNK.IUP_L.., -PER.PFIT.., -RNK.PFIT.., -PER.MDP.., -RNK.MDP.., -dCDF., -dCH., -len.seq.)
cleaned_unrestnet <- ridaodat_unrestnet %>%
  select(-SEQ_ID., -PER.VLXT.., -RNK.VLXT.., -RNK.VSL2B.., -PER.VL3.., -RNK.VL3.., -PER.IUP_S.., -RNK.IUP_S.., -PER.IUP_L.., -RNK.IUP_L.., -PER.PFIT.., -RNK.PFIT.., -PER.MDP.., -RNK.MDP.., -dCDF., -dCH., -len.seq.)

#Remove commas at end of values
cleaned_whole <- cleaned_whole %>%
  mutate(across(c(PER.VSL2B..), ~ gsub(",", "", .)))
cleaned_lim <- cleaned_lim %>%
  mutate(across(c(PER.VSL2B..), ~ gsub(",", "", .)))
cleaned_blood <- cleaned_blood %>%
  mutate(across(c(PER.VSL2B..), ~ gsub(",", "", .)))
cleaned_hum <- cleaned_hum %>%
  mutate(across(c(PER.VSL2B..), ~ gsub(",", "", .)))
cleaned_canon <- cleaned_canon %>%
  mutate(across(c(PER.VSL2B..), ~ gsub(",", "", .)))
cleaned_restnet <- cleaned_restnet %>%
  mutate(across(c(PER.VSL2B..), ~ gsub(",", "", .)))
cleaned_unrestnet <- cleaned_unrestnet %>%
  mutate(across(c(PER.VSL2B..), ~ gsub(",", "", .)))

#Change chr to num for VSL2 
cleaned_whole$PER.VSL2B.. <- as.numeric(cleaned_whole$PER.VSL2B..)
cleaned_lim$PER.VSL2B.. <- as.numeric(cleaned_lim$PER.VSL2B..)
cleaned_blood$PER.VSL2B.. <- as.numeric(cleaned_blood$PER.VSL2B..)
cleaned_hum$PER.VSL2B.. <- as.numeric(cleaned_hum$PER.VSL2B..)
cleaned_canon$PER.VSL2B.. <- as.numeric(cleaned_canon$PER.VSL2B..)
cleaned_restnet$PER.VSL2B.. <- as.numeric(cleaned_restnet$PER.VSL2B..)
cleaned_unrestnet$PER.VSL2B.. <- as.numeric(cleaned_unrestnet$PER.VSL2B..)


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


# Violin plot
proteome_median <- median(
  combined$PER.VSL2B..[combined$secretome == "Human Proteome"],
  na.rm = TRUE
)

ggplot(combined, aes(x = secretome, y = PER.VSL2B.., fill = secretome)) +
  geom_violin(alpha = 0.5, trim = FALSE, color = "grey30", linewidth = 0.3) +
  geom_boxplot(width = 0.1, outlier.shape = 16, alpha = 0.9, fill = "white") +
  geom_hline(
    yintercept = proteome_median,
    linetype = "dashed",
    color = "black",
    linewidth = 0.3
  ) +
  scale_fill_manual(
    values = c(
      "Human Proteome"   = "#1b9e77",
      "Secretome: Blood" = "#d95f02",
      "Secretome: Restricted"   = "#7570b3", 
      "Secretome: Whole" = "#7cb0f0", 
      "Stoneome: Canonical" = "#1e3786",
      "Stoneome: Restricted" = "#ccccff", 
      "Stoneome: Unrestricted" = "#b4eeb4"
    )
  ) +
  theme_classic() +
  labs(
    x = "",
    y = "Percent mean intrinsic disorder (VSL2B)"
  ) +
  coord_cartesian(ylim = c(0, 100)) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(face = "bold")
  )

#View statistics (mean, median, IQR)
stats_secretome <- combined %>%
  group_by(secretome) %>%
  summarise(
    n = n(),
    mean = mean(PER.VSL2B.., na.rm = TRUE),
    median = median(PER.VSL2B.., na.rm = TRUE),
    IQR = IQR(PER.VSL2B.., na.rm = TRUE)
  )

stats_secretome

# View range of values for different secretomes in combined df
combined %>%
  group_by(secretome) %>%
  summarise(
    min = min(PER.VSL2B..),
    max = max(PER.VSL2B..)
  )