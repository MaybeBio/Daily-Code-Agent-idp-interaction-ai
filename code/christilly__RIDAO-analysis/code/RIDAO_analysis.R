library("dplyr")

#Load in data into data frame
ridaodat <- read.table("/file/path/MDP_MODE/results.txt", header=TRUE) #make sure to convert RIDAO results file to txt


#Remove commas at end of values
ridaodat <- ridaodat %>%
  mutate(across(everything(), ~ gsub(",", "", .)))

#Change chr to num for VSL2, dCH, and dCDF columns
ridaodat$PER.VSL2B.. <- as.numeric(ridaodat$PER.VSL2B..)
ridaodat$dCDF. <- as.numeric(ridaodat$dCDF.)
ridaodat$dCH. <- as.numeric(ridaodat$dCH.)

#Find proteins which are in top left (disordered) quadrant of CH-CDF plot
chcdf <- ridaodat %>% filter(dCDF.<0, dCH.>0)
View(chcdf)  #

#Find which proteins have VSL2 percentage >50
VSL2 <- ridaodat %>% filter(PER.VSL2B.. > 50)
#Ordered from highest to lowest disorder 
ord_VSL2 <- VSL2 %>% 
  arrange(desc(PER.VSL2B..))

View(VSL2)  #
View(ord_VSL2)  #VSL2 proteins ranked least to greatest in disorder

#Export data frames to txt file
write.table(chcdf, file = "/file/path/chcdf.txt", row.names = TRUE,
            col.names = TRUE, quote = FALSE)
write.table(VSL2, file = "/file/path/VSL2.txt", row.names = TRUE,
            col.names = TRUE, quote = FALSE)

#Export ENSEMBL sequence IDs only
write.table(chcdf$SEQ_ID., file = "/file/path/chcdf_seq_ID.txt", row.names = TRUE,
            col.names = TRUE, quote = FALSE)
write.table(VSL2$SEQ_ID., file = "/file/path/VSL2_seq_ID.txt", row.names = TRUE,
            col.names = TRUE, quote = FALSE)
write.table(ord_VSL2$SEQ_ID., file = "/file/path/ordered_VSL2_seq_ID.txt", row.names = TRUE,
            col.names = TRUE, quote = FALSE)

#Excel-ready data files
write.csv(ridaodat, file = "/Users/tillytran/Desktop/ridaodata.csv", quote = FALSE) #CSV file of all RIDAO data
write.csv(ord_VSL2$SEQ_ID., file = "/Users/tillytran/Desktop/disordereddat.csv", quote = FALSE) #ENSEMBL/Uniprot ID of disordered proteins 
