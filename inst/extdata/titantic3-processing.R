# Get the variable levels from the Vanderbilt website
library(rvest)
parsed_level_table <- "http://biostat.mc.vanderbilt.edu/wiki/pub/Main/DataSets/Ctitanic3.html" %>%
  read_html %>% # Grab the html
  html_table %>% # parse the table
  `[[`(2) %>% # Grab the table with levels
  `names<-`(c("var", "levels")) # name the data frame columns
for(i in 2:nrow(parsed_level_table)){ # pass the var names down the rows
  if(parsed_level_table$var[i] == "")  parsed_level_table$var[i] <- parsed_level_table$var[i - 1] 
}
parsed_levels <- parsed_level_table %>%
  split(f = parsed_level_table$var) %>% #split by variable
  lapply(function(item) item$levels) %>% # Grab the levels
  lapply(function(item) setdiff(item, "")) # Remove ""
library(dplyr)
# Import the data and clean out the ""
fix_missing <- function(item){
  item[item == ""] <- NA
  item
}
titanic3 <- read.csv("data-raw/titanic3.csv", stringsAsFactors=F) %>%
  lapply(fix_missing) %>% as.data.frame %>% #Remove the ""
  mutate(boat = factor(boat, levels = parsed_levels$boat), # Convert to factors
         cabin = factor(cabin, levels = parsed_levels$cabin),
         embarked = factor(embarked),
         pclass = factor(ifelse(1, "1st", ifelse(2, "2nd", "3rd"))),
         sex = factor(sex))
devtools::use_data(titanic3)
