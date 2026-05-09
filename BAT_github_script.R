###R Script for "Blinding and Quality of Life Endpoints in Angina" study###----
#loading in packages and csv file----
#Load in packages
library(tidyverse)
library(lmtest)
library(sandwich)
library(scales)
library(flextable)
library(officer)
library(stringr)
#so far unused packages
library(forcats)
library(tidyr)
library(broom)
library(gt)
library(gtsummary)
#load in csv file
df <- read_csv("cleaned_data.csv", show_col_types = FALSE)
set.seed(67)
#conversion of types
analysis_df <- df %>%
  mutate(blinded_num = case_when(blinding_binary == "blinded" ~ 1,
      blinding_binary == "unblinded" ~ 0,
      TRUE ~ NA_real_),
    qol_endpoint = as.integer(qol_endpoint),
    year_of_publication = as.numeric(year_of_publication),
    log_n = as.numeric(log_n),
    treatment_type = factor(treatment_type),
    continent = factor(continent),
    decade = factor(decade))
#running poisson regression----
primary_df <- analysis_df %>%
  filter(!is.na(blinded_num),
    !is.na(qol_endpoint),
    !is.na(treatment_type),
    !is.na(year_of_publication),
    !is.na(log_n))


model_primary <- glm(blinded_num ~ qol_endpoint + treatment_type + year_of_publication + log_n,
  family = poisson(link = "log"),
  data = primary_df)

cov_primary <- vcovHC(model_primary, type = "HC0")
primary_coefs <- coeftest(model_primary, vcov = cov_primary)

print(primary_coefs)

primary_table <- as.data.frame(unclass(primary_coefs))
primary_table$term <- rownames(primary_table)
rownames(primary_table) <- NULL

names(primary_table)

names(primary_table)
primary_table <- primary_table %>%
  rename(estimate = Estimate,
    std_error = `Std. Error`,
    z_value = `z value`,
    p_value = `Pr(>|z|)`) %>%
  mutate(aPR = exp(estimate),
    CI_low = exp(estimate - 1.96 * std_error),
    CI_high = exp(estimate + 1.96 * std_error)) %>%
  select(term, estimate, std_error, z_value, p_value, aPR, CI_low, CI_high)

print(primary_table)

#Figure 1- temporal trends----
trend_table <- analysis_df %>%
  group_by(decade) %>%
  summarise(n_trials = n(),
    pct_qol = mean(qol_endpoint == 1, na.rm = TRUE),
    pct_blinded = mean(blinded_num == 1, na.rm = TRUE),
    .groups = "drop")

print(trend_table)

plot_trend_df <- trend_table %>%
  pivot_longer(cols = c(pct_qol, pct_blinded),
    names_to = "metric",
    values_to = "value") %>%
  mutate(metric = case_when(
      metric == "pct_qol" ~ "Quality of Life Endpoint use",
      metric == "pct_blinded" ~ "Blinded Trials"))

#Figure 2- temporal trend figure----
figure2 <- ggplot(
  plot_trend_df,
  aes(x = decade, y = value, group = metric, colour = metric)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3.2) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1),
    expand = expansion(mult = c(0.02, 0.02))) +
  scale_colour_manual(
    values = c(
      "Blinded Trials" = "#D55E00",
      "Quality of Life Endpoint use" = "#0072B2")) +
  labs(title = "Trends in Quality of Life Endpoint Use and Trial Blinding in Angina Trials",
    x = "Decade",
    y = "Percentage of Trials",
    colour = NULL) +
  theme_classic(base_size = 13) +
  theme(legend.position = "top",
    legend.text = element_text(size = 11),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11, colour = "black"),
    plot.title = element_text(face = "bold", size = 14, hjust = 0),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black")
  )

print(figure2)

ggsave(filename = "figure2_temporal_trends.png",plot = figure2,width = 9,height = 5,dpi = 300)

#Figure 3- Forest plot showing poisson regression results----
palette_main <- c("QoL endpoint" = "#D55E00", "Other variables" = "#0072B2")

forest_df <- primary_table %>%
  filter(term != "(Intercept)") %>%
  mutate(label = case_when(term == "qol_endpoint" ~ "QoL endpoint",
      term == "year_of_publication" ~ "Publication year",
      term == "log_n" ~ "Log sample size",
      term == "treatment_typeother" ~ "Treatment: other",
      term == "treatment_typepharmacological" ~ "Treatment: pharmacological",
      TRUE ~ term),
    colour_group = ifelse(term == "qol_endpoint", "QoL endpoint", "Other variables"))

forest_df <- forest_df |>
  mutate(label = recode(label,"Treatment: pharmacological" = "Pharmacological treatment",
      "Treatment: other" = "Other treatment",
      "Publication year" = "Publication year",
      "Log sample size" = "Sample size (log)",
      "QoL endpoint" = "Quality-of-life endpoint"),
    label = factor(label,levels = c("Quality-of-life endpoint",
        "Sample size (log)",
        "Publication year",
        "Other treatment",
        "Pharmacological treatment")),
    colour_group = recode(
      colour_group,
      "QoL endpoint" = "Quality-of-life endpoint",
      "Other variables" = "Other variables"))
figure3 <- ggplot(
  forest_df,
  aes(x = aPR, y = label, colour = colour_group)) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.7,
    colour = "grey40") +
  geom_errorbarh(
    aes(xmin = CI_low, xmax = CI_high),
    height = 0.18,
    linewidth = 0.9) +
  geom_point(size = 3.2) +
  scale_colour_manual(values = c("Other variables" = "#D55E00",
      "Quality-of-life endpoint" = "#0072B2")) +
  scale_x_continuous(
    limits = c(0, 5),
    breaks = seq(0, 5, by = 1)) +
  labs(title = "Factors Associated with Blinded Trial Design in Angina Trials",
    x = "Adjusted Prevalence Ratio",
    y = NULL,
    colour = NULL) +
  theme_classic(base_size = 13) +
  theme(legend.position = "none",
    axis.title.x = element_text(size = 12),
    axis.text.x = element_text(size = 11, colour = "black"),
    axis.text.y = element_text(size = 11, colour = "black", margin = margin(r = 6)),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    plot.title = element_text(face = "bold", size = 14, hjust = 0))
print(figure3)

ggsave("figure3_primary_forest_coloured.png",figure3,width = 8,height = 5,dpi = 300)
#Figure 4- stratified by endpoints----
stratified_table <- analysis_df %>%
  filter(!is.na(blinded_num)) %>%
  group_by(treatment_type, qol_endpoint) %>%
  summarise(n_blinded = sum(blinded_num == 1, na.rm = TRUE),
    n_unblinded = sum(blinded_num == 0, na.rm = TRUE),
    total = n(),
    pct_blinded = mean(blinded_num == 1, na.rm = TRUE),
    .groups = "drop")
plot_stratified_df <- stratified_table %>%
  mutate(qol_group = if_else(qol_endpoint == 1,"Quality of Life Endpoint","No Quality of Life Endpoint"),
    treatment_type = recode(treatment_type,
      "interventional" = "Interventional",
      "other" = "Other",
      "pharmacological" = "Pharmacological"))
treatment_sizes <- plot_stratified_df %>%
  group_by(treatment_type) %>%
  summarise(n = sum(total), .groups = "drop")

plot_stratified_df <- plot_stratified_df %>%
  left_join(treatment_sizes, by = "treatment_type")

plot_stratified_df <- plot_stratified_df %>%
  mutate(treatment_label = paste0(treatment_type, " (n = ", n, ")"),
    treatment_label = factor(treatment_label,levels = paste0(c("Interventional","Other","Pharmacological")," (n = ",treatment_sizes$n[match(c("Interventional","Other","Pharmacological"),treatment_sizes$treatment_type)],")")))

figure4 <- ggplot(plot_stratified_df,aes(x = treatment_label, y = pct_blinded, fill = qol_group)) +
  geom_col(position = position_dodge(width = 0.7),width = 0.65) +
  scale_fill_manual(values = c("No Quality of Life Endpoint" = "#D55E00","Quality of Life Endpoint" = "#0072B2")) +
  scale_y_continuous(labels = percent_format(accuracy = 1),limits = c(0, 1)) +
  labs(title = "Blinded Trial Prevalence by Treatment Type and Quality of Life Endpoint Use",
    x = "Treatment Type",y = "Percentage of Trials Blinded",fill = NULL) +
  theme_classic(base_size = 13) +
  theme(legend.position = "top",
    legend.text = element_text(size = 11),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11, colour = "black"),
    plot.title = element_text(face = "bold", size = 14, hjust = 0),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"))
print(figure4)

ggsave(
  filename = "figure4_stratified_bars.png",
  plot = figure4,
  width = 9,
  height = 5,
  dpi = 300
)
#prisma flowchart----
prisma_plot <- ggplot() +
  coord_cartesian(xlim = c(0, 14), ylim = c(0, 10), expand = FALSE) +
  theme_void() +
  annotate("rect",xmin = 2, xmax = 10,ymin = 7.8, ymax = 9.2,fill = "white", color = "black", linewidth = 0.7) +
  annotate("rect",xmin = 2, xmax = 10,ymin = 1.8, ymax = 3.2,fill = "white", color = "black", linewidth = 0.7) +
  annotate("rect",xmin = 10.8, xmax = 13.8,ymin = 3.0, ymax = 7.0,fill = "white", color = "black", linewidth = 0.7) +
  annotate("text",x = 6, y = 8.5,label = "Records identified through database searching\nn = 2410",size = 6, family = "sans", lineheight = 1.1) +
  annotate("text",x = 6, y = 2.5,label = "Studies included in analysis\nn = 381",size = 6, family = "sans", lineheight = 1.1) +
  annotate("text",x = 12.3, y = 6.35,label = "Excluded at screening (n = 2029)",size = 4.8, family = "sans", fontface = "bold") +
  annotate("text",x = 12.3, y = 5.5,label = paste("Not randomized controlled trial, \n or follow up of trial (n = 710)","Patient group not relevant (n = 1159)","Endpoint not angina-related (n = 146)","Article not available online (n = 14)",sep = "\n"),size = 4.5, family = "sans", lineheight = 1.08) +
  annotate("segment", x = 6, xend = 6,y = 7.8, yend = 3.2,linewidth = 0.7,arrow = arrow(length = unit(0.22, "inches"), type = "closed")) +
  annotate("segment",x = 6, xend = 10.8,y = 5.5, yend = 5.5,linewidth = 0.7,arrow = arrow(length = unit(0.22, "inches"), type = "closed"))
print(prisma_plot)
ggsave("figure1_prisma_diagram.png",prisma_plot,width = 15,height = 8,dpi = 600,bg = "white")


#table 1- grouping trials by QoL/non-QoL----
analysis_df <- df %>%
  mutate(sample_size = exp(log_n),qol_endpoint = factor(qol_endpoint,levels = c(0, 1),labels = c("No QoL endpoint", "QoL endpoint")),treatment_type = treatment_type %>% str_trim() %>% str_to_lower(),treatment_type = factor(treatment_type,levels = c("pharmacological", "interventional", "other"),labels = c("Pharmacological", "Interventional", "Other")),
    blinding_binary = blinding_binary %>% str_trim() %>% str_to_lower(),blinding_binary = factor(blinding_binary,levels = c("blinded", "unblinded"),labels = c("Blinded", "Unblinded")))

table1 <- analysis_df %>%select(qol_endpoint,sample_size,treatment_type,blinding_binary) %>%
  tbl_summary(by = qol_endpoint,statistic = list(sample_size ~ "{median} ({p25}, {p75})",all_categorical() ~ "{n} / {N} ({p}%)"),
    label = list(sample_size ~ "Sample size",treatment_type ~ "Treatment type",blinding_binary ~ "Blinding"),missing = "no") %>%add_overall(last = TRUE) %>% bold_labels() %>% modify_table_styling(columns = label,rows = TRUE,indent = 0)

table1_ft <- as_flex_table(table1)

table1_ft <- as_flex_table(table1) %>%
  add_header_lines(values = "Table 1. Characteristics of Included Trials according to Quality-of-Life Endpoint use") %>% bold(part = "header") %>% bold(i = ~ label %in% c("Sample size", "Treatment type", "Blinding"), bold = TRUE) %>%
  align(j = 1, align = "left", part = "all") %>% align(j = 2:4, align = "center", part = "all") %>%
  border_remove() %>% hline_top(part = "header", border = fp_border(width = 1.2)) %>% hline(part = "header", border = fp_border(width = 1.0)) %>% hline_bottom(part = "body", border = fp_border(width = 1.2)) %>% autofit()

table1_ft <- table1_ft %>% add_footer_lines(values = c("Blinding status was unavailable for 26 trials."))
save_as_docx("Table 1" = table1_ft,path = "table1.docx")

#table 2- poisson regression coded in same format as table 1----
analysis_df <- df %>%mutate(blinded_num = case_when(blinding_binary == "blinded" ~ 1,blinding_binary == "unblinded" ~ 0,TRUE ~ NA_real_),
    qol_endpoint = as.integer(qol_endpoint),treatment_type = treatment_type %>%stringr::str_trim() %>%stringr::str_to_lower(),
    treatment_type = factor(treatment_type,levels = c("interventional", "other", "pharmacological")),
    year_of_publication = as.numeric(year_of_publication),log_n = as.numeric(log_n))

primary_df <- analysis_df %>%filter(!is.na(blinded_num),!is.na(qol_endpoint),!is.na(treatment_type),!is.na(year_of_publication),!is.na(log_n))

model_primary <- glm(blinded_num ~ qol_endpoint + treatment_type + year_of_publication + log_n,family = poisson(link = "log"),data = primary_df)
cov_primary <- vcovHC(model_primary, type = "HC0")
primary_coefs <- coeftest(model_primary, vcov = cov_primary)

primary_table <- as.data.frame(unclass(primary_coefs))
primary_table$term <- rownames(primary_table)
rownames(primary_table) <- NULL

primary_table <- primary_table %>%rename(estimate = `Estimate`,std_error = `Std. Error`,z_value = `z value`,p_value = `Pr(>|z|)`) %>%
  mutate(aPR = exp(estimate),CI_low = exp(estimate - 1.96 * std_error),CI_high = exp(estimate + 1.96 * std_error)) %>%
  filter(term != "(Intercept)") %>%
  mutate(Variable = case_when(term == "qol_endpoint" ~ "QoL endpoint",term == "treatment_typeother" ~ "Treatment type: Other",term == "treatment_typepharmacological" ~ "Treatment type: Pharmacological",term == "year_of_publication" ~ "Publication year",term == "log_n" ~ "Log sample size",TRUE ~ term),
    `Adjusted prevalence ratio` = sprintf("%.2f", aPR),`95% CI` = paste0(sprintf("%.2f", CI_low), " to ", sprintf("%.2f", CI_high)),`P value` = ifelse(p_value < 0.001, "<0.001", sprintf("%.3f", p_value))) %>%
  select(Variable, `Adjusted prevalence ratio`, `95% CI`, `P value`)

primary_table <- primary_table %>%
  mutate(order = case_when(Variable == "QoL endpoint" ~ 1,Variable == "Treatment type: Pharmacological" ~ 2,Variable == "Treatment type: Other" ~ 3,Variable == "Publication year" ~ 4,Variable == "Log sample size" ~ 5,TRUE ~ 99)) %>%arrange(order) %>%select(-order)

forest_results_ft <- flextable(primary_table) %>%
  bold(part = "header") %>% bold(i = ~ Variable == "QoL endpoint", bold = TRUE) %>% align(j = 1, align = "left", part = "all") %>%
  align(j = 2:4, align = "center", part = "all") %>% border_remove() %>% hline_top(part = "header", border = fp_border(width = 1.2)) %>%
  hline(part = "header", border = fp_border(width = 1.0)) %>% hline_bottom(part = "body", border = fp_border(width = 1.2)) %>%
  autofit()

forest_results_ft <- add_header_lines(forest_results_ft,values = "Table 2. Adjusted Associations between Quality-of-Life Endpoint use and Trial Blinding")
doc <- read_docx()
doc <- body_add_flextable(doc, forest_results_ft)

print(doc, target = "Table2_regression_results.docx")

#sensitivity analysis ----
# Sensitivity analysis: regression restricted to trials with primary/PRO endpoint selected ----

analysis_df_sens <- df %>%
mutate(blinded_num = case_when(blinding_binary == "blinded" ~ 1,blinding_binary == "unblinded" ~ 0,TRUE ~ NA_real_),
    pro_endpoint = as.integer(pro_endpoint),qol_endpoint = as.integer(qol_endpoint),treatment_type = treatment_type %>%
      stringr::str_trim() %>%stringr::str_to_lower(),treatment_type = factor(treatment_type,levels = c("interventional", "other", "pharmacological")),year_of_publication = as.numeric(year_of_publication),log_n = as.numeric(log_n))

analysis_df_sens %>%
  count(stated_primary_endpoint)

sensitivity_df <- analysis_df_sens %>%filter(
    stated_primary_endpoint == 1, !is.na(blinded_num),!is.na(qol_endpoint),!is.na(treatment_type),!is.na(year_of_publication),!is.na(log_n)) %>%droplevels()

sensitivity_df %>%summarise(n_trials = n(), n_blinded = sum(blinded_num == 1), n_unblinded = sum(blinded_num == 0),pct_blinded = mean(blinded_num == 1))

sensitivity_df %>%count(treatment_type, qol_endpoint, blinded_num)

model_sensitivity <- glm(blinded_num ~ qol_endpoint + treatment_type + year_of_publication + log_n, family = poisson(link = "log"),data = sensitivity_df)
cov_sensitivity <- vcovHC(model_sensitivity, type = "HC0")
sensitivity_coefs <- coeftest(model_sensitivity, vcov = cov_sensitivity)
print(sensitivity_coefs)

sensitivity_table <- as.data.frame(unclass(sensitivity_coefs))
sensitivity_table$term <- rownames(sensitivity_table)
rownames(sensitivity_table) <- NULL

sensitivity_table <- sensitivity_table %>%
  rename(estimate = `Estimate`,std_error = `Std. Error`,z_value = `z value`,p_value = `Pr(>|z|)`) %>%
  mutate(aPR = exp(estimate),CI_low = exp(estimate - 1.96 * std_error),CI_high = exp(estimate + 1.96 * std_error)) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    Variable = case_when(term == "qol_endpoint" ~ "QoL endpoint",term == "treatment_typeother" ~ "Treatment type: Other",term == "treatment_typepharmacological" ~ "Treatment type: Pharmacological",term == "year_of_publication" ~ "Publication year",term == "log_n" ~ "Log sample size",TRUE ~ term),`Adjusted prevalence ratio` = sprintf("%.2f", aPR),`95% CI` = paste0(sprintf("%.2f", CI_low), " to ", sprintf("%.2f", CI_high)),`P value` = ifelse(p_value < 0.001, "<0.001", sprintf("%.3f", p_value))) %>%
  select(Variable, `Adjusted prevalence ratio`, `95% CI`, `P value`) %>%
  mutate(order = case_when(Variable == "QoL endpoint" ~ 1,Variable == "Treatment type: Pharmacological" ~ 2,Variable == "Treatment type: Other" ~ 3,Variable == "Publication year" ~ 4,Variable == "Log sample size" ~ 5,TRUE ~ 99)) %>% arrange(order) %>%select(-order)
print(sensitivity_table)


sensitivity_ft <- flextable(sensitivity_table) %>% bold(part = "header") %>% bold(i = ~ Variable == "QoL endpoint", bold = TRUE) %>% align(j = 1, align = "left", part = "all") %>% align(j = 2:4, align = "center", part = "all") %>% border_remove() %>% hline_top(part = "header", border = fp_border(width = 1.2)) %>% hline(part = "header", border = fp_border(width = 1.0)) %>%hline_bottom(part = "body", border = fp_border(width = 1.2)) %>% autofit()
 
sensitivity_ft <- add_header_lines( sensitivity_ft, values = "Table S1. Sensitivity Analysis Restricted to Trials with a Stated Primary Endpoint")

sensitivity_ft <- add_footer_lines(sensitivity_ft, values = c( "Analysis restricted to trials where stated_primary_endpoint = 1.", "Adjusted prevalence ratios were estimated using modified Poisson regression with robust standard errors.", "Model adjusted for quality-of-life endpoint use, treatment type, publication year, and log sample size."))

doc_sens <- read_docx()
doc_sens <- body_add_flextable(doc_sens, sensitivity_ft)

print(doc_sens, target = "TableS1_sensitivity_primary_endpoint.docx")

# Expanded Table 1 - grouping trials by QoL/non-QoL ----

analysis_df <- df %>%
  mutate(sample_size = exp(log_n),qol_endpoint = factor(qol_endpoint, levels = c(0, 1),labels = c("No QoL endpoint", "QoL endpoint")), treatment_type = treatment_type %>% str_trim() %>% str_to_lower(), treatment_type = factor( treatment_type, levels = c("pharmacological", "interventional", "other"), labels = c("Pharmacological", "Interventional", "Other")),blinding_binary = blinding_binary %>%str_trim() %>%str_to_lower(),blinding_binary = factor( blinding_binary,levels = c("blinded", "unblinded"), labels = c("Blinded", "Unblinded")),stated_primary_endpoint = factor(stated_primary_endpoint, levels = c(0, 1),labels = c("No", "Yes")),year_of_publication = as.numeric(year_of_publication),publication_period = case_when(year_of_publication < 1970 ~ "Before 1970",year_of_publication >= 1970 & year_of_publication < 1980 ~ "1970–1979",year_of_publication >= 1980 & year_of_publication < 1990 ~ "1980–1989",year_of_publication >= 1990 & year_of_publication < 2000 ~ "1990–1999",year_of_publication >= 2000 & year_of_publication < 2010 ~ "2000–2009",year_of_publication >= 2010 & year_of_publication < 2020 ~ "2010–2020",year_of_publication >= 2020 & year_of_publication <= 2026 ~ "2020–2026", TRUE ~ NA_character_),publication_period = factor( publication_period,levels = c("Before 1970", "1970–1979", "1980–1989", "1990–1999", "2000–2009", "2010–2020", "2020–2026" ) ) )

table1 <- analysis_df %>%select(qol_endpoint,publication_period,sample_size,stated_primary_endpoint,treatment_type,blinding_binary ) %>%tbl_summary(by = qol_endpoint,statistic = list(sample_size ~ "{mean} ({sd})",all_categorical() ~ "{n} / {N} ({p}%)" ),digits = list( sample_size ~ c(1, 1), all_categorical() ~ c(0, 0, 1)),label = list(publication_period ~ "Publication period",sample_size ~ "Sample size",stated_primary_endpoint ~ "Primary endpoint specified", treatment_type ~ "Treatment type", blinding_binary ~ "Blinding"),missing = "no") %>%add_overall(last = TRUE) %>%bold_labels() %>%modify_table_styling(columns = label,rows = TRUE,indent = 0)

table1_ft <- as_flex_table(table1) %>%add_header_lines(values = "Table 1. Characteristics of Included Trials according to Quality-of-Life Endpoint use") %>%bold(part = "header") %>%bold( i = ~ label %in% c("Publication period","Sample size","Primary endpoint specified", "Treatment type","Blinding"),bold = TRUE) %>%align(j = 1, align = "left", part = "all") %>%align(j = 2:4, align = "center", part = "all") %>%border_remove() %>%hline_top(part = "header", border = fp_border(width = 1.2)) %>%hline(part = "header", border = fp_border(width = 1.0)) %>%hline_bottom(part = "body", border = fp_border(width = 1.2)) %>%autofit()

table1_ft <- table1_ft %>%add_footer_lines( values = c("Values are n / N (%) unless otherwise stated.","Sample size is presented as mean (standard deviation).","Publication period 2010–2020 includes trials published from 2010 to 2019; 2020–2026 includes trials published from 2020 onwards.","Blinding status was unavailable for 26 trials."))

save_as_docx("Table 1" = table1_ft, path = "table1_expanded.docx")

#expanded sample size 

analysis_df <- df %>%
  mutate(sample_size = round(exp(log_n)),sample_size_category = case_when(sample_size < 50 ~ "<50",sample_size >= 50 & sample_size < 100 ~ "50–99",sample_size >= 100 & sample_size < 250 ~ "100–249",sample_size >= 250 & sample_size < 500 ~ "250–499",sample_size >= 500 ~ "≥500",TRUE ~ NA_character_),sample_size_category = factor( sample_size_category,levels = c("<50", "50–99", "100–249", "250–499", "≥500")),qol_endpoint = factor(qol_endpoint, levels = c(0, 1), labels = c("No QoL endpoint", "QoL endpoint")), treatment_type = treatment_type %>% str_trim() %>%str_to_lower(),treatment_type = factor(treatment_type, levels = c("pharmacological", "interventional", "other"),labels = c("Pharmacological", "Interventional", "Other")),blinding_binary = blinding_binary %>%str_trim() %>% str_to_lower(),blinding_binary = factor( blinding_binary, levels = c("blinded", "unblinded"), labels = c("Blinded", "Unblinded")),stated_primary_endpoint = factor( stated_primary_endpoint, levels = c(0, 1), labels = c("No", "Yes")),year_of_publication = as.numeric(year_of_publication),publication_period = case_when( year_of_publication < 1970 ~ "Before 1970", year_of_publication >= 1970 & year_of_publication < 1980 ~ "1970–1979", year_of_publication >= 1980 & year_of_publication < 1990 ~ "1980–1989", year_of_publication >= 1990 & year_of_publication < 2000 ~ "1990–1999", year_of_publication >= 2000 & year_of_publication < 2010 ~ "2000–2009", year_of_publication >= 2010 & year_of_publication < 2020 ~ "2010–2019", year_of_publication >= 2020 & year_of_publication <= 2026 ~ "2020–2026",TRUE ~ NA_character_),publication_period = factor( publication_period,levels = c("Before 1970","1970–1979","1980–1989","1990–1999","2000–2009","2010–2019", "2020–2026")))

table1 <- analysis_df %>%select(qol_endpoint,publication_period,sample_size_category,stated_primary_endpoint,treatment_type,blinding_binary) %>% tbl_summary(by = qol_endpoint,statistic = list(all_categorical() ~ "{n} / {N} ({p}%)"),digits = list(all_categorical() ~ c(0, 0, 1)),label = list(publication_period ~ "Publication period",sample_size_category ~ "Sample size",stated_primary_endpoint ~ "Primary endpoint specified",treatment_type ~ "Treatment type",blinding_binary ~ "Blinding"),missing = "no") %>% add_overall(last = TRUE) %>% bold_labels() %>% modify_table_styling(columns = label,rows = TRUE,indent = 0)

table1_ft <- as_flex_table(table1) %>%add_header_lines(values = "Table 1. Characteristics of Included Trials according to Quality-of-Life Endpoint use") %>% bold(part = "header") %>%bold(i = ~ label %in% c("Publication period", "Sample size","Primary endpoint specified","Treatment type", "Blinding" ),bold = TRUE) %>%align(j = 1, align = "left", part = "all") %>% align(j = 2:4, align = "center", part = "all") %>%border_remove() %>%hline_top(part = "header", border = fp_border(width = 1.2)) %>%hline(part = "header", border = fp_border(width = 1.0)) %>%hline_bottom(part = "body", border = fp_border(width = 1.2)) %>%autofit()

table1_ft <- table1_ft %>%add_footer_lines(values = c("Values are n / N (%).", "Sample size is presented as categories rather than median and interquartile range.","Publication period 2010–2019 includes trials published from 2010 to 2019; 2020–2026 includes trials published from 2020 onwards.", "Blinding status was unavailable for 26 trials." ))

save_as_docx("Table 1" = table1_ft,path = "table1_expanded.docx")

#updated figure 2----

analysis_df <- df %>%mutate(blinded_num = case_when(blinding_binary == "blinded" ~ 1,blinding_binary == "unblinded" ~ 0,TRUE ~ NA_real_),qol_endpoint = as.integer(qol_endpoint),year_of_publication = as.numeric(year_of_publication),decade_revised = case_when(year_of_publication < 1970 ~ "Before 1970",year_of_publication >= 1970 & year_of_publication < 1980 ~ "1970s",year_of_publication >= 1980 & year_of_publication < 1990 ~ "1980s",year_of_publication >= 1990 & year_of_publication < 2000 ~ "1990s",year_of_publication >= 2000 & year_of_publication < 2010 ~ "2000s",year_of_publication >= 2010 & year_of_publication < 2020 ~ "2010s",year_of_publication >= 2020 & year_of_publication <= 2026 ~ "2020s",TRUE ~ NA_character_),decade_revised = factor(decade_revised,levels = c("Before 1970","1970s","1980s","1990s","2000s","2010s","2020s")))

trend_table <- analysis_df %>%filter(!is.na(decade_revised)) %>%group_by(decade_revised) %>%summarise(n_trials = n(),pct_qol = mean(qol_endpoint == 1, na.rm = TRUE),pct_blinded = mean(blinded_num == 1, na.rm = TRUE),.groups = "drop")
plot_trend_df <- trend_table %>% pivot_longer(cols = c(pct_qol, pct_blinded),names_to = "metric",values_to = "value") %>% mutate(metric = case_when(metric == "pct_qol" ~ "Quality of Life Endpoint use",metric == "pct_blinded" ~ "Blinded Trials"))

figure2 <- ggplot(plot_trend_df,aes(x = decade_revised, y = value, group = metric, colour = metric)) + geom_line(linewidth = 1.2) + geom_point(size = 3.2) +scale_y_continuous(labels = percent_format(accuracy = 1),limits = c(0, 1),expand = expansion(mult = c(0.02, 0.02))) +scale_colour_manual(values = c("Blinded Trials" = "#D55E00","Quality of Life Endpoint use" = "#0072B2")) + labs(title = "Trends in Quality of Life Endpoint Use and Trial Blinding in Angina Trials",x = "Publication Decade",y = "Percentage of Trials",colour = NULL) +theme_classic(base_size = 13) +theme(legend.position = "top",legend.text = element_text(size = 11),axis.title = element_text(size = 12),axis.text.x = element_text(size = 10, colour = "black", angle = 0, hjust = 0.5),axis.text.y = element_text(size = 11, colour = "black"),plot.title = element_text(face = "bold", size = 14, hjust = 0),axis.line = element_line(colour = "black"),axis.ticks = element_line(colour = "black"))
ggsave(filename = "figure2_temporal_trends_updated.png",plot = figure2,width = 9,height = 5,dpi = 300)

appendix_table <- df
colnames(df)
#appendix table
yes_no <- function(x) {
  case_when(x == 1 ~ "Yes",x == 0 ~ "No",TRUE ~ "Missing")}
detailed_trial_table <- df %>%
  mutate(sample_size = round(exp(as.numeric(log_n))),title = str_squish(title),continent = str_to_title(str_squish(continent)),year_of_publication = as.numeric(year_of_publication),publication_decade = case_when(year_of_publication < 1970 ~ "Before 1970",year_of_publication >= 1970 & year_of_publication < 1980 ~ "1970s",year_of_publication >= 1980 & year_of_publication < 1990 ~ "1980s",year_of_publication >= 1990 & year_of_publication < 2000 ~ "1990s",year_of_publication >= 2000 & year_of_publication < 2010 ~ "2000s",year_of_publication >= 2010 & year_of_publication < 2020 ~ "2010s",year_of_publication >= 2020 ~ "2020s",TRUE ~ NA_character_),publication_decade = factor(publication_decade,levels = c("Before 1970","1970s","1980s","1990s","2000s", "2010s","2020s")),stated_primary_endpoint = yes_no(as.integer(stated_primary_endpoint)),pro_endpoint = yes_no(as.integer(pro_endpoint)),qol_endpoint = yes_no(as.integer(qol_endpoint)),physician_endpoint = yes_no(as.integer(physician_endpoint)),physiological_endpoint = yes_no(as.integer(physiological_endpoint)),holter_endpoint = yes_no(as.integer(holter_endpoint)),blinding_binary = blinding_binary %>%str_trim() %>%str_to_lower(),blinding_binary = case_when( blinding_binary == "blinded" ~ "Blinded", blinding_binary == "unblinded" ~ "Unblinded",is.na(blinding_binary) | blinding_binary == "" ~ "Missing",TRUE ~ "Missing"),
treatment_type = treatment_type %>%str_trim() %>%str_to_lower(),treatment_type = recode(treatment_type,"pharmacological" = "Pharmacological","interventional" = "Interventional","other" = "Other")) %>%select(`Trial title` = title,`Continent` = continent,`Decade` = publication_decade,`Sample size` = sample_size,`Treatment type` = treatment_type,`Blinding` = blinding_binary,`Primary endpoint stated` = stated_primary_endpoint,`Patient-reported endpoint` = pro_endpoint,`QoL endpoint` = qol_endpoint,`Physician endpoint` = physician_endpoint,`Exercise endpoint` = physiological_endpoint,`Holter endpoint` = holter_endpoint) %>% arrange( `Trial title`)

detailed_trial_ft <- flextable(detailed_trial_table) %>% add_header_lines(values = "Supplementary Table S1. Detailed Characteristics of Included Trials") %>% bold(part = "header") %>% fontsize(size = 7, part = "all") %>% fontsize(size = 8, part = "header") %>% align(j = "Trial title", align = "left", part = "all") %>% align(j = c("Continent","Decade","Sample size","Treatment type","Blinding","Primary endpoint stated","Patient-reported endpoint","QoL endpoint", "Physician endpoint", "Exercise endpoint","Holter endpoint"),align = "center",part = "all") %>% valign(valign = "top", part = "all") %>% border_remove() %>% hline_top(part = "header", border = fp_border(width = 1.2)) %>% hline(part = "header", border = fp_border(width = 1.0)) %>% hline_bottom(part = "body", border = fp_border(width = 1.2)) %>% set_table_properties(layout = "fixed", width = 1)
detailed_trial_ft <- detailed_trial_ft %>% width(j = "Trial title", width = 2.6) %>% width(j = "Continent", width = 0.65) %>% width(j = "Decade", width = 0.55) %>% width(j = "Sample size", width = 0.65) %>% width(j = "Treatment type", width = 0.85) %>% width(j = "Blinding", width = 0.65) %>% width(j = "Primary endpoint stated", width = 0.75) %>% width(j = "Patient-reported endpoint", width = 0.75) %>% width(j = "QoL endpoint", width = 0.55) %>% width(j = "Physician endpoint", width = 0.65) %>% width(j = "Exercise endpoint", width = 0.65) %>%width(j = "Holter endpoint", width = 0.55)
detailed_trial_ft <- detailed_trial_ft %>%add_footer_lines(values = c("Sample size was calculated by exponentiating log-transformed sample size and rounding to the nearest whole number.","The 2020s category includes trials published from 2020 onwards.","Missing blinding status is shown explicitly."))
doc <- read_docx()
landscape_section <- prop_section(page_size = page_size(orient = "landscape",width = 11.69,height = 8.27),page_margins = page_mar(top = 0.35,bottom = 0.35,left = 0.35,right = 0.35))
doc <- body_set_default_section(doc, value = landscape_section)
doc <- body_add_flextable(doc, detailed_trial_ft)
print( doc,target = "Supplementary_Table_S1_detailed_trial_characteristics.docx")
