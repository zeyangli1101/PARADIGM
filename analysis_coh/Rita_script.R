library(dplyr)
library(stringr)
library(survival)
library(survminer)
library(broom)
library(ggplot2)



#1. filtering allo-HCT patients  
df_allo_clinical = read.csv('alloHCT_clinical_core.csv')
#select patients that are allo
df_allo_allo <- df_allo_clinical %>%
  filter(transplant_type == "allo",
         !transplant_classification %in% c("dli", "non-dli gvt",
                                           "other cellular therapy"))
#filter out the first transplant 
df_allo_first <- df_allo_allo %>%
  group_by(mrn) %>%
  summarise(
    transplant = if (all(is.na(transplant_date))) NA else min(transplant_date, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(transplant_type = "allo")
n_distinct(df_allo_first$mrn) #4596
#filter out allo-patients in medication exposure data
df_allo_cleaned = read.csv('alloHCT_medication_pcdm_transplant date included.csv')
df_allo_filtered <- semi_join(df_allo_cleaned, df_allo_first, by = "mrn")
n_distinct(df_allo_filtered$mrn, na.rm = TRUE) #--3950 patients 
#export file for patient bacteria response score calculation 
write.csv(df_allo_filtered, "alloHCT_medication_pcdm_allofiltered.csv", row.names = FALSE)





#2. calculate patient-specific bacteria reponse score 
df_medication_history = read.csv('alloHCT_medication_pcdm_transplant date included.csv')
df_drug_match = read.csv('summarized drug matching list_matched.csv')
#calculate order start date relative to transplant and order end date relative to transplant 
str(df_medication_history$order_start)
str(df_medication_history$transplant)
str(df_medication_history$order_end)
#convert everything in those columns to date 
df_medication_history$transplant <- as.Date(df_medication_history$transplant, format = "%Y-%m-%d")
df_medication_history$order_start <- as.Date(df_medication_history$order_start, format = "%Y-%m-%d")
df_medication_history$order_end <- as.Date(df_medication_history$order_end, format = "%Y-%m-%d")
#recalculate relative days to hct
df_medication_history$order_start_rel_hct <- as.numeric(df_medication_history$order_start - df_medication_history$transplant)
df_medication_history$order_end_rel_hct   <- as.numeric(df_medication_history$order_end - df_medication_history$transplant)
#select rows from drug matching list that are only"strong match"
df_strong_drug = df_drug_match[df_drug_match$match_category == "strong match",]
df_strong_drug_selected = df_strong_drug[,c("COH_medication_name","paradigm_name_matched")]
#merge df_strong_drug_selected with df_medication_history using medication name columns
df_merged = inner_join(
  df_strong_drug_selected,
  df_medication_history,
  by = c("COH_medication_name" = "order_description")
)
#check order_status values:
unique(df_merged$order_status) #----"Dispensed", "Discontinued", "Sent", NA, "Suspend", "Completed", "Verified" 
#for now just rule out "discontinued" and "suspend"??????????????
df_merged <- df_merged %>%
  filter(!order_status %in% c("Discontinued", "Suspend"))
df_merged_selected = df_merged[,c("COH_medication_name","paradigm_name_matched", "transplant","order_start","order_end",
                                  "order_start_rel_hct","order_end_rel_hct", "mrn")]
#select medication orders start and end date overlapping (day-14, day 14) time window relative to transplant 
#eg. day -20 to day 4, ---overlapping day -14 to day 4
#eg. day 5 to day 20, ---overlapping day 5 to day 14
#eg. day 2 to day 10, ---overlapping day 2 to day 10
#eg. day -20 to day -16, ---no oerlapping
#eg. day 15 to day 20, ---no overlapping 
#filter out order start date >=day 14--everything start before end of time window 
#filter out order end date <=day -14--everything end after beginning of time window 
df_merged_filtered <- df_merged_selected %>%
  filter(
    order_end_rel_hct   >= -14,
    order_start_rel_hct <=  14
  )
#calculate patient% for each exposure to limit 5-90% 
total_patients <- df_merged_filtered %>%
  distinct(mrn) %>%
  nrow()
# count unique patients per exposure
exposure_summary <- df_merged_filtered %>%
  distinct(mrn, paradigm_name_matched) %>%   # remove duplicates per patient-drug
  count(paradigm_name_matched, name = "n_exposed") %>%
  mutate(
    pct = n_exposed / total_patients * 100
  )
#keep 5-90% exposures
exposure_filtered <- exposure_summary %>%
  filter(pct >= 5 & pct <= 90)
df_merged_filtered <- df_merged_filtered %>%
  filter(paradigm_name_matched %in% exposure_filtered$paradigm_name_matched)
#calculate overlapping exposure days of (order start, oder end) to (day-14, day 14)
#eg. start date = -20, end date = 4
df_exposure <- df_merged_filtered %>%
  mutate(
    overlap_start = pmax(order_start_rel_hct, -14), #compare -20 and -14 -->-14
    overlap_end   = pmin(order_end_rel_hct,   14), #compare 4 and 14 --> 4
    exposure_days = overlap_end - overlap_start #use overlapend(4) - overlapstart(-14) to calculate exposure days
  )
#situations exist where exposure_day =0
#consider that exposure is only one time? 
#clean up the exposure days to have positive values 
df_exposure <- df_exposure %>%
  mutate(
    exposure_days = case_when(
      exposure_days == 0 ~ 1, #if exposure days =0, make it into 1
      TRUE ~ exposure_days
    )
  ) %>%
  filter(exposure_days >= 0) #keep only exposure days >=0
#sum exposure times per (drug,mrn) pair
df_exposure_total <- df_exposure %>%
  group_by(paradigm_name_matched, mrn) %>%
  summarise(
    total_exposure = sum(exposure_days, na.rm = TRUE),
    .groups = "drop"
  )
#incorporate PARADIGM drug specific response score 
df_bacteria_score <- read.csv("updated bacteria response score.csv")
#merge the two dt by exposure_name 
df_patient_score <- df_exposure_total %>%
  inner_join(df_bacteria_score, by = c("paradigm_name_matched" = "exposure_name"))
#change name for simpson_reciprical
df_patient_score <- df_patient_score %>%
  rename (alpha_diversity = simpson_reciprocal)
#calculate response score per patient per drug = exposure count * response score for that drug
#for every bacteria species, not just enterococcus
df_patient_score <- df_patient_score %>%
  mutate(
    across(
      c(Blautia, Enterococcus, Erysipelatoclostridium, Bifidobacterium, Collinsella,
      Eubacterium, Fusicatenibacter, Paraprevotella, Prevotella, Roseburia, 
      Stenotrophomonas, Sutterella, Veillonella, alpha_diversity),
      ~ .x * total_exposure,
      .names = "score_{.col}"
    )
  )
df_patient_score_final <- df_patient_score %>%
  select(
    mrn,
    starts_with("score_")
  )
# Sum across all scores per bacteria for each MRN
df_patient_score_final <- df_patient_score_final %>%
  group_by(mrn) %>%
  summarize(
    across(
      starts_with("score_"),
      ~ sum(.x, na.rm = TRUE),
      .names = "total_{.col}"
    ),
    .groups = "drop"
  )
#remove"total_"
df_patient_score_final <- df_patient_score_final %>%
  rename_with(~ sub("^total_score_", "", .))
#export file for cox plot 
write.csv(df_patient_score_final, "alloHCT_total patient response score.csv", row.names = FALSE)





#3. graph ocx plot 
df_medication = read.csv('alloHCT_total patient response score.csv')
df_survival = read.csv('alloHCT_clinical_post.csv')
#select allo-HCT patients included in cleaned medication list
df_survival_matched <- df_survival %>%
  semi_join(df_medication, by = "mrn")
n_distinct(df_survival_matched$mrn) #--3884
#select columns related to overall survival calculation
df_os <- df_survival_matched %>%
  select(
    mrn,
    date_of_transplant,
    vital_status,
    days_from_hct_to_last_contact_vital_status_or_death) %>%
  distinct()
#calculate survival days and assign survival/dead
df_os <- df_os %>%
  mutate(
  time_os  = as.numeric(days_from_hct_to_last_contact_vital_status_or_death),
  event_os = case_when(
    tolower(vital_status) %in% c("dead", "deceased", "d") ~ 1,
    tolower(vital_status) %in% c("alive", "living", "a") ~ 0,
    vital_status %in% c(1, "1") ~ 1,
    vital_status %in% c(0, "0") ~ 0,
    TRUE ~ NA_real_)
  )
#QC check
# Any weird times?
summary(df_os$time_os)
df_os %>% filter(is.na(time_os) | time_os < 0) %>% head()
#clean up patients alive but time_os = NA
df_os_clean <- df_os %>%
  filter(!is.na(time_os) & time_os >= 0)
#keeping mrns with max survival time and max event_os(death)
df_os_final <- df_os_clean %>%
  group_by(mrn) %>%
  summarise(
    date_of_transplant = min(date_of_transplant, na.rm = TRUE),
    time_os  = max(time_os, na.rm = TRUE),
    event_os = max(event_os, na.rm = TRUE),
    .groups = "drop"
  )
#merge covariables 
df_clinicalvariables = read.csv('alloHCT_clinical_core.csv')
df_os_full <- df_os_final %>%
  left_join(df_clinicalvariables, by = c("mrn", "date_of_transplant" = "transplant_date"))
df_os_full <- df_os_full %>%
  select(
    mrn,
    date_of_transplant,
    time_os,
    event_os,
    patient_age_at_hct,
    patient_gender,
    conditioning_regimens,
    conditioning_regimen_intensity,
    stem_cell_source,
    primary_diagnosis_at_transplant
  )
#group conditioning intensity
df_os_full <- df_os_full %>%
  mutate(
    conditioning_group = case_when(
      grepl("ablative", conditioning_regimen_intensity, ignore.case = TRUE) &
        !grepl("non", conditioning_regimen_intensity, ignore.case = TRUE) ~ "Ablative",
      grepl("non-myeloablative|reduced", conditioning_regimen_intensity, ignore.case = TRUE) ~ "Reduced/Nonmyeloablative",
      TRUE ~ NA_character_
    )
  )
#group graft source
df_os_full <- df_os_full %>%
  mutate(
    graft_type = case_when(
      grepl("cord blood", stem_cell_source, ignore.case = TRUE) ~ "cord",
      grepl("^bone marrow$", stem_cell_source, ignore.case = TRUE) ~ "BM unmodified",
      grepl("^peripheral stem cells$", stem_cell_source, ignore.case = TRUE) ~ "PBSC unmodified",
      TRUE ~ "Other"
    )
  )
#group disease type
df_os_full <- df_os_full %>%
  mutate(
    disease_group = case_when(
      grepl("\\(aml\\)|acute myeloid leukemia|acute promyelocytic leukemia|myeloid leukemia, nos", 
            primary_diagnosis_at_transplant, ignore.case = TRUE) ~ "AML",
      
      grepl("myelodysplastic syndrome|myeloproliferative neoplasms|chronic myelomonocytic leukemia", 
            primary_diagnosis_at_transplant, ignore.case = TRUE) ~ "MDS/MPN",
      
      grepl("non-hodgkin lymphoma|mycosis fungoides|lymphoblastic lymphoma", 
            primary_diagnosis_at_transplant, ignore.case = TRUE) ~ "NHL",
      
      grepl("acute lymphocytic leukemia|\\(all\\)|lymphoid leukemia, nos|acute biphenotypic leukemia|undifferentiated leukemia", 
            primary_diagnosis_at_transplant, ignore.case = TRUE) ~ "ALL",
      
      grepl("multiple myeloma", 
            primary_diagnosis_at_transplant, ignore.case = TRUE) ~ "Myeloma",
      
      grepl("chronic lymphocytic leukemia|\\(cll\\)|prolymphocytic leukemia", 
            primary_diagnosis_at_transplant, ignore.case = TRUE) ~ "CLL",
      
      grepl("hodgkin lymphoma", 
            primary_diagnosis_at_transplant, ignore.case = TRUE) ~ "Hodgkins",
      
      grepl("chronic myeloid leukemia|\\(cml\\)", 
            primary_diagnosis_at_transplant, ignore.case = TRUE) ~ "CML",
      
      grepl("aplastic anemia|fanconi anemia|pure red cell", 
            primary_diagnosis_at_transplant, ignore.case = TRUE) ~ "AA",
      
      TRUE ~ "other"
    )
  )
#categorize gender, conditioning intensity, graft source, and disease type
df_os_full <- df_os_full %>%
  mutate(
    conditioning_group = as.factor(conditioning_group),
    graft_type = as.factor(graft_type),
    disease_group = as.factor(disease_group),
    patient_gender = as.factor(patient_gender)
  )
df_os_selected <-df_os_full %>%
  select(
    mrn,
    date_of_transplant,
    time_os,
    event_os,
    patient_age_at_hct,
    patient_gender,
    conditioning_group,
    graft_type,
    disease_group
  )
#merge the overall survival to bacteria response score for patient 
df_score_full <- df_os_full %>%
  inner_join(df_medication, by = "mrn")
#list bacteria scores 
score_vars <- c(
  "alpha_diversity",
  "Blautia",
  "Enterococcus",
  "Erysipelatoclostridium",
  "Bifidobacterium",
  "Collinsella",
  "Eubacterium",
  "Fusicatenibacter",
  "Paraprevotella",
  "Prevotella",
  "Roseburia",
  "Stenotrophomonas",
  "Sutterella",
  "Veillonella"
)
#Fit a Cox model for each score 
covars <- "patient_age_at_hct + patient_gender + conditioning_group + graft_type + disease_group"
results_list <- lapply(score_vars, function(sc) {
  fit <- coxph(
    as.formula(
      paste0("Surv(time_os, event_os) ~ ", sc, " + ", covars)
    ),
    data = df_score_full
  )
  tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term == sc) %>%
    mutate(score = sc)
})
df_results <- bind_rows(results_list)
#Ensure scores are plotted in the order of score_vars
df_results$score <- factor(df_results$score, levels = rev(score_vars))
#Make a label for p-values (nice formatting)
df_results <- df_results %>%
  mutate(
    p_label = paste0("p=", format.pval(p.value, digits = 2, eps = 1e-3)),
    score   = factor(score, levels = rev(score_vars))
  )
#plot cox plot 
ggplot(df_results, aes(x = estimate, y = score, colour = score)) +
  geom_point(size = 4) +  # bigger dots
  geom_errorbarh(
    aes(xmin = conf.low, xmax = conf.high),
    height = 0.25,
    linewidth = 1.2       # thicker CI lines
  ) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 1.2       # thicker reference line
  ) +
  # p-values pushed further right
  geom_text(
    aes(x = Inf, label = p_label),
    colour = "black",
    hjust = -0.1,
    size = 5
  ) +
  scale_x_continuous(trans = "log10") +
  coord_cartesian(clip = "off") +
  theme_minimal() +
  theme(
    legend.position = "none",
    
    # axis text (Blautia / Enterococcus / etc.)
    axis.text.y = element_text(size = 20, face = "bold"),
    axis.text.x = element_text(size = 20, face = "bold"),
    
    # axis title
    axis.title.x = element_text(size = 20, face = "bold"),
    
    # plot title
    plot.title = element_text(
      size = 25,
      face = "bold",
      hjust = 0.5
    ),
    
    # extra space on right for p-values
    plot.margin = margin(5.5, 120, 5.5, 5.5)
  ) +
  labs(
    x = "Hazard ratio (95% CI)",
    y = NULL,
    title = "Overall survival vs. bacteria response scores"
  )
#limit to the most recent 500 patients?
df_recent <- df_score_full %>%
  arrange(desc(date_of_transplant)) %>%
  slice_head(n = 500)
recent_results_list <- lapply(score_vars, function(sc) {
  fit <- coxph(
    as.formula(
      paste0("Surv(time_os, event_os) ~ ", sc, " + ", covars)
    ),
    data = df_recent
  )
  tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term == sc) %>%
    mutate(score = sc)
})
df_recent_results <- bind_rows(recent_results_list)
#Ensure scores are plotted in the order of score_vars
df_recent_results$score <- factor(df_recent_results$score, levels = rev(score_vars))
#Make a label for p-values (nice formatting)
df_recent_results <- df_recent_results %>%
  mutate(
    p_label = paste0("p=", format.pval(p.value, digits = 2, eps = 1e-3)),
    score   = factor(score, levels = rev(score_vars))
  )
#plot cox plot 
ggplot(df_recent_results, aes(x = estimate, y = score, colour = score)) +
  geom_point(size = 4) +  # bigger dots
  geom_errorbarh(
    aes(xmin = conf.low, xmax = conf.high),
    height = 0.25,
    linewidth = 1.2       # thicker CI lines
  ) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 1.2       # thicker reference line
  ) +
  # p-values pushed further right
  geom_text(
    aes(x = Inf, label = p_label),
    colour = "black",
    hjust = -0.1,
    size = 5
  ) +
  scale_x_continuous(trans = "log10") +
  coord_cartesian(clip = "off") +
  theme_minimal() +
  theme(
    legend.position = "none",
    
    # axis text (Blautia / Enterococcus / etc.)
    axis.text.y = element_text(size = 20, face = "bold"),
    axis.text.x = element_text(size = 20, face = "bold"),
    
    # axis title
    axis.title.x = element_text(size = 20, face = "bold"),
    
    # plot title
    plot.title = element_text(
      size = 25,
      face = "bold",
      hjust = 0.5
    ),
    
    # extra space on right for p-values
    plot.margin = margin(5.5, 120, 5.5, 5.5)
  ) +
  labs(
    x = "Hazard ratio (95% CI)",
    y = NULL,
    title = "Overall survival vs. bacteria response scores"
  )

