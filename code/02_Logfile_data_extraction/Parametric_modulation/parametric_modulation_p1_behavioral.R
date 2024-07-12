####Parametric modulation, working on the logfiles and making tsv file for each subject########
#Author Anna E. Fromm 
#Edited by Mohamed Abdelmotaleb
###Session Info
#R version 4.0.5 (2021-03-31) -- "Shake and Throw"
# only extraction of real task (not control condition)
######################################


rm(list = ls())

library(dplyr)
library(psycho)
library(tidyverse)
library(openxlsx)
library(data.table)
library(stringr)

# Define the root directory of your study
root_dir <- "/media/Data03/Studies/MeMoSLAP"

# List all subjects
subjects <- list.dirs(path = root_dir, recursive = FALSE, full.names = TRUE)

# Initialize an empty list to store log file paths
log_files <- list()


# Loop through each subject
for (subject in subjects) {
  # List all session directories for the current subject
  sessions <- list.dirs(path = subject, recursive = FALSE, full.names = TRUE)
  
  # Loop through each session
  for (session in sessions) {
    # Define the path to the behavioral (beh) folder
    beh_dir <- file.path(session, "beh")
    
    # Find the log file in the beh folder
    log_file <- list.files(path = beh_dir, pattern = "OLMM.*\\.log$", full.names = TRUE)
    
    # If a log file is found, store its path in the log_files list
    if (length(log_file) > 0) {
      log_files <- append(log_files, log_file)
    }
  }
}

# Print the list of log file paths
print(log_files)

# skip = 3 to exclude description 
df1 <- log_files %>%
  map_df(read.table,  sep="\t", header=TRUE, na.strings =  " ", fill=TRUE, skip=3)



############# done read all learning files of the long version #####


####working with the data set ##########

# t time is the difference between the response time - picture time: reaction time
# variable: "Code" contains information whether the response was correct/ incorrect/ too late
#h88_p16_i_down.jpg;block1;LS1 (block_stage) -> gives information of learning stage and block
#Picture	correct;h88_p35_k_down.jpg (feedback) -> gives information regarding the response of the participant 

#keep only relevant rows (Picture) ; because the Response is included within the feedback
df1 <- df1%>%
  filter(!(Code == "extra Pic Time"))

df1 <- df1 %>%
  group_by(Subject) %>%
  dplyr::mutate(RT = case_when(
    grepl("incorrect", Code) ~ shift(TTime, 1, type = "lag"),
    grepl("correct", Code) ~ shift(TTime, 1, type = "lag"),
    grepl("late", Code) ~ "NaN"
  ))
df1$RT <- as.numeric(df1$RT)

df2 <- df1 %>% filter(Event.Type %in% c("Picture"))

#drop fixation, drop instructions
df2 <- df2 %>%
  filter(!(Code == "fix")) %>%
  filter(!grepl("bubbles", Code)) %>%
  filter(!(Code == "Hello - Intro"))


#create new variable (block)
df3 <- df2 %>%
  group_by(Subject) %>%
  mutate(Block = case_when(
    grepl("block1", Code) ~ "1",
    grepl("block2", Code) ~ "2",
    grepl("block3", Code) ~ "3",
    grepl("block4", Code) ~ "4",
    grepl("block5", Code) ~ "5",
    grepl("block6", Code) ~ "6",
    grepl("crt1", Code) ~ "c1",
    grepl("crt2", Code) ~ "c2"))

df3 <- df3 %>%
  group_by(Subject) %>%
  fill(Block, .direction = "down")

#create new variable Stage
df3 <- df3 %>%
  group_by(Block) %>%
  mutate(Stage = case_when(
    grepl("LS1", Code) ~ "1",
    grepl("LS2", Code) ~ "2",
    grepl("LS3", Code) ~ "3",
    grepl("LS4", Code) ~ "4"))

df3$Task_Type = ifelse(grepl(c("c"), df3$Block), "control", "learning")


df3 <- df3 %>%
  group_by(Subject) %>%
  fill(Stage, .direction = "down")


#create new variable Block_Stage
df3$Block_Stage <- with(df3, paste(Block, Stage, sep = "_"))


# done forming the structure of for the analyses ##########


df4 <- df3 %>%
  filter(!grepl("block", Code)) %>%
  filter(!grepl("LS", Code))

#all these previous codes is to keep only the "feedback" rows to make the analysis of the correctness

#### calculate number of correct responses
#get correct, incorrect, to late
df4$Response <- gsub(";.*","",df4$Code)


df5 <- df4 %>%
  mutate(Response_num = case_when(Response == "correct" ~ '1',
                                  Response == "incorrect" ~ '0',
                                  Response == "to late" ~ '0'))

df5$Response_num <- as.numeric(df5$Response_num)

#extract the session(out of naming)
#sub001p1_s3_task1 (session 3, from subject 001)
df5 <- df5 %>%
  mutate(Session = case_when(
    grepl("s3", Subject, ignore.case=T) ~ "3",
    grepl("s4", Subject, ignore.case=T) ~ "4"))


## consider Presentation's way of formating Reaction Time (0.1 ms)
df5$RT <- (df5$RT/10)

# Calculate the mean of Response_num for each Stage of each learning block for each Subject

df5_learning <- df5 %>%
  filter(Block %in% 1:6) %>%
  group_by(Subject, Block, Stage) %>%
  summarize (mean_learning_response = mean(Response_num, na.rm = TRUE)*100) %>%
  ungroup()


# Calculate the mean of Response_num for each Stage of each control block for each Subject
df5_control <- df5 %>%
  filter(Block %in% c('c1', 'c2')) %>%
  group_by(Subject,Block, Stage) %>%
  summarize(mean_control_response = mean(Response_num, na.rm = TRUE)*100 ) %>%
  ungroup()

# Merge the mean values with the original data frame
df5 <- df5 %>%
  left_join(df5_learning, by = c("Subject", "Block", "Stage"), suffix = c("", "_learning")) %>%
  left_join(df5_control, by = c("Subject","Block", "Stage"), suffix = c("", "_control")) %>%
  mutate(mean_response_acc= case_when(
    Task_Type == "learning" ~ mean_learning_response,
    Task_Type == "control" ~ mean_control_response
  )) %>%
  select(-mean_learning_response, -mean_control_response)


#df5 <- df5 %>%
  #group_by(Subject,Block,Stage) %>%
  #mutate(mean_Response_block = round(mean(Response_num, na.rm = TRUE)*100,2 )) %>% 
  #ungroup()

#df5 <- df5 %>%
  #group_by(Subject,Block,Stage) %>%
  #mutate(mean_Response_block = round(mean(Response_num, na.rm = TRUE)*100,2 )) %>% 
  #ungroup()
#RT
df5_learning_RT <- df5 %>%
  filter(Block %in% 1:6) %>%
  group_by(Subject, Block, Stage) %>%
  summarize (mean_learning_RT = mean(RT, na.rm = TRUE)) %>%
  ungroup()


# Calculate the mean of Response_num for each Stage of each control block for each Subject
df5_control_RT <- df5 %>%
  filter(Block %in% c('c1', 'c2')) %>%
  group_by(Subject,Block, Stage) %>%
  summarize(mean_control_RT = mean(RT, na.rm = TRUE) ) %>%
  ungroup()

# Merge the mean values with the original data frame
df5 <- df5 %>%
  left_join(df5_learning_RT, by = c("Subject", "Block", "Stage"), suffix = c("", "_learning")) %>%
  left_join(df5_control_RT, by = c("Subject","Block", "Stage"), suffix = c("", "_control")) %>%
  mutate(mean_response_RT= case_when(
    Task_Type == "learning" ~ mean_learning_RT,
    Task_Type == "control" ~ mean_control_RT
  )) %>%
  select(-mean_learning_RT, -mean_control_RT)

#df5 <- df5 %>%
  #group_by(Stage) %>%
  #mutate(mean_RT = round(mean(RT, na.rm = TRUE), 3)) %>%
  #ungroup()

df5$ID <- gsub("p.*","",df5$Subject)
df5$ID <- gsub("[^0-9.]", "",df5$ID)

# split Subject column
df5 <- df5 %>% separate(Subject, into = c("Subject_", "Session", "Task"), sep = "_", remove = FALSE)

df5$Subject_ <- gsub("sub", "sub-", df5$Subject_)
df5$Subject_ <- gsub("p1", "", df5$Subject_)
df5$Session <- gsub("s", "ses-", df5$Session)
df5$Task <- gsub("task", "", df5$Task)

# Define the directory to save the TSV files
output_dir <- "/media/Data03/Studies/MeMoSLAP/derivatives/preprocessing/SPM-standard/"


# Get the unique subjects
unique_subjects <- unique(df5$Subject_)
unique_sessions <- unique(df5$Session)
unique_tasks <- unique(df5$Task)

# Loop through each unique subject and export corresponding rows to a TSV file
for (Subject_num in unique_subjects){
  for (Session_num in unique_sessions){
    # Filter the rows for the current subject
    Subject_df <- df5 %>%
      filter(Subject_ == Subject_num & Session == Session_num)
#for debugging 
    #Subject_df <- df5 %>%
    #  filter(Subject_ == unique_subjects[1] & Session == unique_sessions[1])

    Task <- unique(Subject_df$Task)[1]
    #Task <- Task_l[1]
    # Define the file name
    #exp /media/Data03/Studies/MeMoSLAP/derivatives/preprocessing/SPM-standard/sub-001/ses-3/beh
    file_name <- paste0(output_dir,
                                     Subject_num, "/", Session_num, "/beh/", Subject_num, "_",Session_num, "_", "task-OLMM","_","acq-",Task,"_", "desc-pmod.tsv")
    # debug
    #file_name <- paste0(output_dir,
    #                                 unique_subjects[1], "/", unique_sessions[1], "/beh/", unique_subjects[1], "_",unique_sessions[1], "_", "task-OLMM","_","acq-",Task,"_", "desc-pmod.tsv")
    
    # Write the data to a TSV file
    # View(file_name)
    write.table(Subject_df, file = file_name, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

    rm(Subject_df)
  }
}
# Print a message indicating completion
cat("TSV files have been successfully created in", output_dir, "\n")
