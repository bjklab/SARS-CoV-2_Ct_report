#' ################################
#' load libraries and set seed
#' ################################
library(tidyverse)
library(readxl)
library(janitor)
library(stringr)
library(lubridate)
library(ggsci)
library(gt)
set.seed(16)




#' ################################
#' load and clean long-form data
#' ################################

#' filepath to long-format Ct data
long_data_filepath <- "../CT Values CSV 2021-03-16-10-51-18.csv"

#' possible instruments
instrument_search_string <- "BDMax|Cepheid|LIAT|m2000|Roche|Simplexa|Thermo"


vroom::vroom(long_data_filepath) %>%
  # clean variable names and select variables
  janitor::clean_names() %>%
  select(accession_number, discrete_task_display, contains("date"), contains("order_"), contains("result")) %>%
  select(-birth_date) %>%
  select(-perform_result_ascii_text) %>%
  # recognize NA values
  mutate(numeric_result_value = replace(numeric_result_value, numeric_result_value > 100 | numeric_result_value == 88, NA)) %>%
  # format dates to permit temporal analysis
  mutate_at(.vars = vars(contains("date")), .funs = ~ as.Date(.x, format = "%Y/%m/%d")) %>%
  mutate(discrete_task_display = gsub("SARS-CoV-2 \\(|\\)","",discrete_task_display)) %>%
  # ensure platform and amplicon are clearly documented
  rename(amplicon = discrete_task_display,
         platform = text_result_value) %>%
  group_by(accession_number) %>%
  mutate(platform = stringr::str_extract(string = paste0(platform, collapse = " "), pattern = instrument_search_string),
         platform = replace(platform, is.na(platform), "Platform Not Recorded")) %>%
  ungroup() %>%
  # add year and week values
  mutate(year = lubridate::year(container_in_lab_date_time),
         week = lubridate::week(container_in_lab_date_time),
         year_week = glue::glue("{year}-{stringr::str_pad(week, width = 2, side = 'left', pad = '0')}")) %>%
  # filter to include only recorded Ct values
  filter(!is.na(numeric_result_value)) %>%
  identity() -> ct

ct




#' ################################
#' summaries of Ct values from long-form data
#' ################################

# table: count of measured Ct values
ct %>%
  group_by(platform, amplicon) %>%
  summarise(n = n_distinct(accession_number)) %>%
  ungroup() %>%
  filter(grepl("Ct",amplicon)) %>%
  gt() %>%
  gt::cols_label(platform = "Platform", amplicon = "Amplicon", n = "# Specimens with Measured Ct Values") %>%
  gt::opt_table_lines(extent = "all")


# table: aggregate central tendency/spread of Ct values
ct %>%
  group_by(platform, amplicon) %>%
  summarise(median_ct = median(numeric_result_value, na.rm = TRUE),
            mean_ct = mean(numeric_result_value, na.rm = TRUE),
            sd_ct = sd(numeric_result_value, na.rm = TRUE),
            first_quartile_ct = quantile(numeric_result_value, probs = 0.25, na.rm = TRUE),
            third_quartile_ct = quantile(numeric_result_value, probs = 0.75, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(mean_sd = glue::glue("{round(mean_ct,1)} ({round(sd_ct,1)})"),
         median_iqr = glue::glue("{round(median_ct,1)} ({round(first_quartile_ct,1)} - {round(third_quartile_ct,1)})")) %>%
  select(platform, amplicon, median_iqr, mean_sd) %>%
  gt() %>%
  gt::cols_label(platform = "Platform",
                 amplicon = "Amplicon",
                 median_iqr = "Median (IQR) Ct",
                 mean_sd = "Mean (SD) Ct") %>%
  gt::opt_table_lines(extent = "all")


# figures
ct %>%
  ggplot(data = .) +
  geom_boxplot(aes(x = amplicon, y = numeric_result_value, fill = amplicon)) +
  facet_wrap(facets = ~ platform, scales = "free") +
  theme_bw() +
  theme(axis.text = element_text(color = "black"),
        strip.background = element_blank()) +
  ggsci::scale_fill_nejm() +
  labs(x = "", y = "Ct Value", fill = "RT-PCR\nTarget")


# weekly changes in central tendency/spread of Ct values
ct %>%
  # only include platforms with at least 4 weeks of data
  group_by(platform) %>%
  filter(n_distinct(year_week) >= 4) %>%
  ungroup() %>%
  ggplot(data = .) +
  geom_boxplot(aes(x = year_week, y = numeric_result_value, fill = amplicon)) +
  facet_grid(rows = vars(amplicon), cols = vars(platform), scales = "free") +
  theme_bw() +
  theme(axis.text = element_text(color = "black"),
        strip.background = element_blank()) +
  ggsci::scale_fill_nejm() +
  labs(x = "Year-Week", y = "Ct Value", fill = "RT-PCR\nTarget")




#' ################################
#' dropout analysis
#' ################################

# data
ct %>%
  mutate(amplicon_code = case_when(amplicon == "Ct ORF1ab" ~ "target_1",
                                   amplicon == "Ct S gene" ~ "target_2",
                                   amplicon == "Ct N gene" ~ "target_3",
                                   amplicon == "Ct N1" ~ "target_4",
                                   amplicon == "Ct N2" ~ "target_5",
                                   )) %>%
  select(accession_number, container_in_lab_date_time, year, week, year_week, platform, amplicon_code, numeric_result_value) %>%
  distinct() %>%
  filter(complete.cases(.)) %>%
  pivot_wider(id_cols = c(accession_number, container_in_lab_date_time, year, week, year_week, platform),
              names_from = amplicon_code,
              values_from = numeric_result_value) %>%
  group_by(accession_number) %>%
  mutate(dropout = (is.na(target_1) & !is.na(target_2)) | (is.na(target_2) & !is.na(target_1)) | (is.na(target_4) & !is.na(target_5)) | (is.na(target_5) & !is.na(target_4))) %>%
  mutate(abs_difference = abs(target_1 - target_2),
         abs_difference = replace(abs_difference, is.na(abs_difference), abs(target_4 - target_5))) %>%
  ungroup() %>%
  identity() -> ct_drop_dif

ct_drop_dif
  

# table  
ct_drop_dif %>%
  group_by(year_week) %>%
  summarise(n_two_target = n_distinct(accession_number),
            total_dropout = sum(dropout, na.rm = TRUE),
            s_gene_dropout = sum(dropout & is.na(target_2), na.rm = TRUE),
            tdo_prop = sum(dropout, na.rm = TRUE) / n_distinct(accession_number),
            sgd_prop = sum(dropout & is.na(target_2), na.rm = TRUE) / n_distinct(accession_number),
            total_report = glue::glue("{total_dropout} ({round(tdo_prop*100,1)}%)"),
            s_gene_report = glue::glue("{s_gene_dropout} ({round(sgd_prop*100,1)}%)")) %>%
  ungroup() %>%
  select(year_week, n_two_target, contains("report")) %>%
  gt() %>%
  gt::cols_label(total_report = "Any Dropout (%)",
                 s_gene_report = "S Gene Dropout (%)",
                 year_week = "Year - Week Number",
                 n_two_target = "Number of Two-Target Specimens") %>%
  gt::opt_table_lines(extent = "all")


# figure
ct_drop_dif %>%
  group_by(year_week) %>%
  summarise(n_two_target = n_distinct(accession_number),
            total_dropout = sum(dropout, na.rm = TRUE),
            s_gene_dropout = sum(dropout & is.na(target_2), na.rm = TRUE),
            tdo_prop = sum(dropout, na.rm = TRUE) / n_distinct(accession_number),
            sgd_prop = sum(dropout & is.na(target_2), na.rm = TRUE) / n_distinct(accession_number),
            total_report = glue::glue("{total_dropout} ({round(tdo_prop*100,1)}%)"),
            s_gene_report = glue::glue("{s_gene_dropout} ({round(sgd_prop*100,1)}%)")) %>%
  ungroup() %>%
  select(year_week, contains("prop")) %>%
  pivot_longer(cols = contains("prop"), names_to = "target_type", values_to = "proportion_dropout") %>%
  mutate(target_type = case_when(grepl("tdo",target_type) ~ "Any Dropout",
                                 grepl("sgd",target_type) ~ "S Gene Dropout")) %>%
  identity() %>%
  ggplot(data = .) +
  geom_line(aes(x = year_week, y = proportion_dropout, group = target_type, color = target_type), alpha = 0.8) +
  geom_point(aes(x = year_week, y = proportion_dropout, color = target_type), alpha = 0.8) +
  theme_bw() +
  theme(axis.text = element_text(color = "black"),
        strip.background = element_blank(),
        legend.position = "top") +
  ggsci::scale_color_nejm() +
  labs(x = "Year-Week", y = "Proportion Dropout", color = "")

  










