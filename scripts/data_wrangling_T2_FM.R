### data wrangling for Marginal Zone compartmnet
rm(list = ls())
gc()

fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # remove + after exponent, if exists. E.g.: (e^+2 -> e^2)
  l <- gsub("e\\+","e",l)  
  # turn the 'e' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # convert 1x10^ or 1.000x10^ -> 10^
  l <- gsub("\\'1[\\.0]*\\'\\%\\*\\%", "", l)
  # return this as an expression
  parse(text=l)
}

log10minorbreaks=as.numeric(1:10 %o% 10^(4:8))


col_pal <- wesanderson::wes_palette(name = "Zissou1",  15, type = "continuous")
col_pal2 <- wesanderson::wes_palette(name = "Zissou1",  2, type = "continuous")


my_theme <- theme(axis.text = element_text(size = 11),
                  axis.title =  element_text(size = 11, face = "bold"),
                  plot.title = element_text(size = 11,  hjust = 0.5),
                  legend.text = element_text(size = 11),
                  legend.title = element_text(size = 11, face = "bold"))


library(tidyverse)

# total counts and donor fractions for the source poppulation
source_T2 <- read.csv("datafiles/T2_as_sourceSpline.csv") 

source_donor <- readxl::read_excel(path = "data/MZ_Master.xlsx", sheet = 2) %>%
  #filter(row_number() %% 2 == 1)%>% filter(is.na(notes)) %>%
  filter(is.na(Ki67))%>% filter(is.na(notes)) %>%
  select(contains("Lamis"), contains("days"), contains("age"), contains(parent_pop))%>% 
  mutate(Total_T1 = Transitional2) %>%  
  mutate_if(is.character, as.numeric) %>% unique()

## Ki67 data in the source compartmenmt
# total counts and donor fractions for the source poppulation
source_host_ki67 <- readxl::read_excel(path = "data/MZ_Master.xlsx", sheet = 1) %>%
  #filter(row_number() %% 2 == 1)%>% filter(is.na(notes)) %>%
  filter(!is.na(Ki67))%>% filter(is.na(notes)) %>%
  select(contains("Lamis"), contains("days"), contains("age"), contains(parent_pop))%>%
  mutate(ki67Pos_T1 = Transitional2) %>% 
  mutate_if(is.character, as.numeric) %>% unique()

source_donor_ki67 <- readxl::read_excel(path = "data/MZ_Master.xlsx", sheet = 2) %>%
  #filter(row_number() %% 2 == 1)%>% filter(is.na(notes)) %>%
  filter(!is.na(Ki67))%>% filter(is.na(notes)) %>%
  select(contains("Lamis"), contains("days"), contains("age"), contains(parent_pop))%>%
  mutate(ki67Pos_T1 = Transitional2) %>% 
  mutate_if(is.character, as.numeric) %>% unique()


# merging total counts for host and donor compartments
source_counts <- full_join(source_host, source_donor, by = c("Lamis.ID", "days.post.bmt", "age.at.S1K", "age.at.bmt"), suffix= c(".host", ".donor"))

# merging ki67Pos counts from host and donor
source_ki67 <- full_join(source_host_ki67, source_donor_ki67, by = c("Lamis.ID", "days.post.bmt", "age.at.S1K", "age.at.bmt"), suffix= c(".host", ".donor"))

# calculating total counts, donor fractions anf ki67Pos fractions
source_data <- full_join(source_ki67, source_counts, by = c("Lamis.ID", "days.post.bmt", "age.at.S1K", "age.at.bmt"), suffix= c(".host", ".donor"))%>%
  mutate(total_counts = Total_T1.host + Total_T1.donor,
         fd = (Total_T1.donor)/total_counts,
         # ki67 fractions without the immature FM compartment
         host_ki67_Tra = ki67Pos_T1.host/ Total_T1.host,
         donor_ki67_Tra = ki67Pos_T1.donor/ Total_T1.donor)%>%
  select(-contains(".host"), -contains(".donor"))

g1 <- ggplot(source_data) +
  geom_point(aes(x = age.at.S1K, y = total_counts), alpha = 0.8, size =3)+
  scale_y_log10(limits = c(1e5, 2e7), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  scale_x_continuous(limits = c(50, 800), trans = "log10", breaks = c(60, 200, 600)) + 
  labs(x = "Host age (days)", y = NULL, title = "T1 total counts") +
  theme_bw() + guides(colour = FALSE) + my_theme
  

g2 <- ggplot(source_data) +
  geom_point(aes(x = days.post.bmt, y = fd), alpha = 0.8, size =3)+
  geom_hline(yintercept = 1.00, linetype = 2)+
  ylim(0, 1.01) + scale_x_continuous(breaks=c(200, 400, 600))+
  labs(x = "Host age (days)", y = NULL, title = "Chimerism in T1") +
  theme_bw() + guides(colour = FALSE) +
  theme_bw() + guides(colour = FALSE) + my_theme


g3 <- ggplot(source_data)+
  geom_point(aes(x = age.at.S1K, y =  host_ki67_Tra), alpha = 0.8, size =2)+
  geom_point(aes(x = age.at.S1K, y = donor_ki67_Tra), alpha = 0.8, size =2)+
  scale_color_manual(values = col_pal2, name = "Subset", labels = c("Donor", "Host"))+
  ylim(0, 1) + scale_x_log10(limits = c(50, 750), breaks = c(60, 200, 600)) +  theme_bw()+
  labs(x = "Host Age (days)", y = NULL, title = "Proportions of Ki67Hi cells within T1") +
    theme_bw() + my_theme


# total counts, donor fractions and ki67Pos fractions  for the target poppulation (FM)
counts_host <- readxl::read_excel(path = "data/MZ_Master.xlsx", sheet = 1)%>%
  #mutate_at(c("immature_FM"), funs(lead), n = 1 )%>%
  filter(is.na(Ki67)) %>% filter(is.na(notes))%>%
  #filter(row_number() %% 2 == 1)%>% filter(is.na(notes)) %>%
  select(contains("Lamis"), contains("days"), contains("age"), contains(target_pop)) %>% 
  mutate(Total_MZ = MZ) %>% 
  mutate_if(is.character, as.numeric) %>% unique()

counts_donor <- readxl::read_excel(path = "data/MZ_Master.xlsx", sheet = 2)%>%
  #mutate_at(c("immature_FM"), funs(lead), n = 1 )%>%
  #filter(row_number() %% 2 == 1)%>% filter(is.na(notes)) %>%
  filter(is.na(Ki67)) %>% filter(is.na(notes))%>%
  select(contains("Lamis"), contains("days"), contains("age"), contains(target_pop))%>%
  mutate(Total_MZ = MZ) %>% 
  mutate_if(is.character, as.numeric) %>% unique()

# ki67 for FM
ki67_host <- readxl::read_excel(path = "data/MZ_Master.xlsx", sheet = 1)%>%
  filter(!is.na(Ki67)) %>% filter(is.na(notes))%>%
  select(contains("Lamis"), contains("days"), contains("age"), contains(target_pop))%>% unique()%>%
  mutate(ki67Pos_MZ = MZ) %>%
  mutate_if(is.character, as.numeric) %>% unique()

ki67_donor <- readxl::read_excel(path = "data/MZ_Master.xlsx", sheet = 2)%>%
  filter(!is.na(Ki67)) %>% filter(is.na(notes))%>%
  select(contains("Lamis"), contains("days"), contains("age"), contains(target_pop))%>% unique()%>%
  mutate(ki67Pos_MZ = MZ)%>%  select(-contains("SP"), -contains("LN"))%>%
  mutate_if(is.character, as.numeric) %>% unique()

# merging total counts for host and donor compartments
MZ_counts <- full_join(counts_host, counts_donor, by = c("Lamis.ID", "days.post.bmt", "age.at.S1K", "age.at.bmt"), suffix= c(".host", ".donor"))

# merging ki67Pos counts from host and donor
MZ_ki67 <- full_join(ki67_host, ki67_donor, by = c("Lamis.ID", "days.post.bmt", "age.at.S1K", "age.at.bmt"), suffix= c(".host", ".donor"))

# calculating total counts, donor fractions anf ki67Pos fractions
MZ_data <- full_join(MZ_ki67, MZ_counts, by = c("Lamis.ID", "days.post.bmt", "age.at.S1K", "age.at.bmt"), suffix= c(".host", ".donor"))%>%
  mutate(total_counts = Total_MZ.host + Total_MZ.donor,
         fd = Total_MZ.donor/ total_counts,
         # ki67 fractions without the immature FM compartment
         host_ki67_MZ = ki67Pos_MZ.host/ Total_MZ.host,
         donor_ki67_MZ = ki67Pos_MZ.donor/ Total_MZ.donor)%>%
  select(-contains(".host"), -contains(".donor"))

# normalising donor fraction in FM by dividing with the donor fractions in the source compartment
MZ_fd <- MZ_data%>%
  select(-contains("ki67")) %>%
  full_join(source_data, by = c("Lamis.ID", "days.post.bmt", "age.at.S1K", "age.at.bmt"), suffix= c(".MZ", ".source"))%>%
  mutate(Nfd = fd.MZ/ fd.source) %>% filter(Nfd < 1.5) %>% filter(Lamis.ID != 848700) %>%
  select(contains("Lamis"), contains("days"), contains("age"), contains("Nfd"))


counts_MZ <- MZ_data%>%
  right_join(MZ_fd, by = c("Lamis.ID", "days.post.bmt", "age.at.S1K", "age.at.bmt")) %>%
  select(contains("Lamis"), contains("days"), contains("age"), contains("total_counts"))

ki67_MZ <-  MZ_data %>%
  right_join(MZ_fd, by = c("Lamis.ID", "days.post.bmt", "age.at.S1K", "age.at.bmt")) %>%
  select(contains("Lamis"), contains("days"), contains("age"), contains("ki67"))%>%
  select(contains("Lamis"), contains("days"), contains("age"),  contains("MZ"))


ki67_plot_MZ <- ki67_MZ %>%
  gather(-c(Lamis.ID, days.post.bmt, age.at.S1K, age.at.bmt), key = subpopulation, value = ki67_fraction)


#plots

## open graphics device to save plots
pdf(file = file.path(getwd(), "figures", paste("MZ","Plots%03d.pdf", sep = "")),
    width = 8, height = 5.5, onefile = FALSE, useDingbats = FALSE)

g1 <- ggplot(counts_MZ) +
  geom_point(aes(x = age.at.S1K, y = total_counts, col = as.factor(age.at.bmt)), stroke = 0, alpha = 0.8, size =2)+
  scale_color_manual(values = col_pal)+
  scale_y_log10(limits = c(5e5, 1e7), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  scale_x_continuous(limits = c(50, 750), trans = "log10", breaks = c(60, 200, 600)) + 
  labs(x = "Host age (days)", y = NULL, title = "MZ total counts") +
  theme_bw() + guides(colour = FALSE) + my_theme


g2 <- ggplot(MZ_fd) +
  geom_point(aes(x = age.at.S1K, y = Nfd, color = as.factor(age.at.bmt)), stroke = 0, alpha = 0.8, size =2)+
  geom_hline(yintercept = 1.00, linetype = 2)+
  scale_color_manual(values = col_pal)+
  ylim(0, 1.25) + scale_x_continuous(breaks=c(200, 400, 600))+
  labs(x = "Host age (days)", y = NULL, title = "Donor fractions normalised to chimerism in T1") +
  theme_bw() + guides(colour = FALSE) + my_theme

g3 <- ggplot(ki67_plot_MZ)+
  geom_point(aes(x = age.at.S1K, y = ki67_fraction, color = subpopulation), stroke = 0, alpha = 0.8, size =2)+
  scale_color_manual(values = col_pal2, name = "Subset", labels = c("Donor", "Host"))+
  ylim(0, 1) + scale_x_log10(limits = c(50, 750), breaks = c(60, 200, 600)) +  theme_bw()+
  labs(x = "Host Age (days)", y = NULL, title = "Proportions of Ki67Hi cells within MZ") +
  theme_bw() + my_theme + theme(legend.position = c(0.8, 0.8), legend.background = element_blank())

lay1 <- rbind(c(1,1,2,2),
              c(NA,3,3,NA))


gridExtra::grid.arrange(g1, g2, g3, layout_matrix = lay1)

dev.off()

# writting files to the data directory 
write.csv(counts_MZ, "data/counts_MZ.csv", row.names = FALSE)
write.csv(MZ_fd, "data/Nfd_MZ.csv", row.names = FALSE)
write.csv(ki67_MZ, "data/Ki67_MZ.csv", row.names = FALSE)
write.csv(source_data, "data/T1_MZ.csv", row.names = FALSE)





