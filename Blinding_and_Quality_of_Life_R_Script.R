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


forest_df <- forest_df %>% 
  mutate(label = recode(label,"Treatment: pharmacological" = "Pharmacological treatment","Treatment: other" = "Other treatment","Publication year" = "Publication year","Log sample size" = "Sample size (log)","QoL endpoint" = "Quality of Life endpoint"),
    label = factor(label,levels = c("Quality of Life endpoint","Sample size (log)","Publication year","Other treatment","Pharmacological treatment")),colour_group = recode(colour_group,"QoL endpoint" = "Quality of Life endpoint","Other variables" = "Other variables"))

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
      "Quality of Life endpoint" = "#0072B2")) +
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
  annotate("text",x = 12.3, y = 5.5,label = paste("Not randomised controlled trial \n or follow up of trial (n = 710)","Patient group not relevant (n = 1159)","Endpoint not angina related (n = 146)","Article not available online (n = 14)",sep = "\n"),size = 4.5, family = "sans", lineheight = 1.08) +
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

