# genotyping full exp 2

# libraries
library(tidyverse)
library(forcats)
library(readxl)
library(plotly)
library(magrittr)
library(scales)

theme_set(theme_classic())
# data import -------------------------------------------------------------
# metadata
meta <-
  read_excel("data/Full exp 2/2023_Sept_20_ip_DFP_fullexp2_meta.xlsx") %>%
  mutate(
    fish_id = as.character(fish_id),
    HomeTank = as.factor(HomeTank),
    treatment = forcats::fct_inorder(treatment),
    sex = as.factor(sex),
    Trial = as.factor(Trial),
    ymazePosition = as.factor(ymazePosition),
    ymazeUnit = as.factor(ymazeUnit),
    day.inj = as.factor(day.inj),
    start.time = as.factor(start.time)
  )

# results file
data.results <-
  read_delim("data/Full exp 2/genotyping_data/KN_20230925_deferiproneExp2_nagluhomxhet8m_genotyp-data-for-KB.txt",
             delim = "\t",
             skip = 8) %>%
  dplyr::slice(1:96) %>%
  mutate(group = case_when(
    Well %in% c("A1", paste0("F", 9:12)) ~ "neg control",
    Well == "A2" ~ "wt control",
    Well == "A3" ~ "het control",
    Well == "A4" ~ "hom control",
    TRUE ~ "unknown"
  )) %>%
    dplyr::select(Well, group, Tm1, Tm2, Tm3)

#temperatures
temps <- read_delim("data/Full exp 2/genotyping_data/KN_20230925_deferiproneExp2_nagluhomxhet8m_genotyp-data-for-KB.txt",
           delim = "\t",
           skip = 114)  %>%
  dplyr::slice(1:96) %>%
  mutate(group = case_when(
    `Well Location` %in% c("A1", paste0("F", 9:12)) ~ "neg control",
    `Well Location` == "A2" ~ "wt control",
    `Well Location` == "A3" ~ "het control",
    `Well Location` == "A4" ~ "hom control",
    TRUE ~ "unknown"
  )) %>%
  dplyr::select(position.for.geno = `Well Location`,
                group,
                all_of(starts_with("Reading")))


# normailsed fluorescence
fluors <- read_delim("data/Full exp 2/genotyping_data/KN_20230925_deferiproneExp2_nagluhomxhet8m_genotyp-data-for-KB.txt",
           delim = "\t",
           skip = 214) %>%
  dplyr::slice(1:96) %>%
  mutate(group = case_when(
    `Well Location` %in% c("A1", paste0("F", 9:12)) ~ "neg control",
    `Well Location` == "A2" ~ "wt control",
    `Well Location` == "A3" ~ "het control",
    `Well Location` == "A4" ~ "hom control",
    TRUE ~ "unknown"
  )) %>%
  dplyr::select(position.for.geno = `Well Location`,
                group,
                all_of(starts_with("Reading")))

# derivities
dfs <- read_delim("data/Full exp 2/genotyping_data/KN_20230925_deferiproneExp2_nagluhomxhet8m_genotyp-data-for-KB.txt",
           delim = "\t",
           skip = 314) %>%
  dplyr::slice(1:96) %>%
  mutate(group = case_when(
    `Well Location` %in% c("A1", paste0("F", 9:12)) ~ "neg control",
    `Well Location` == "A2" ~ "wt control",
    `Well Location` == "A3" ~ "het control",
    `Well Location` == "A4" ~ "hom control",
    TRUE ~ "unknown"
  )) %>%
  dplyr::select(position.for.geno = `Well Location`,
                group,
                all_of(starts_with("Reading")))

# join together
data.meltcurve <- dfs %>%
  pivot_longer(names_to = "Reading", values_to = "df",
               starts_with("Reading")
  ) %>%
  left_join(temps %>%
              pivot_longer(names_to = "Reading", values_to = "temp",
            starts_with("Reading")
            )
  ) %>%
  left_join(fluors %>%
              pivot_longer(names_to = "Reading", values_to = "fluor",
                           starts_with("Reading") )
  ) %>%
  mutate(df = as.numeric(df),
         temp = as.numeric(temp),
         fluor = as.numeric(fluor),
         group = fct_inorder(group))


# fluorescence
data.meltcurve %>%
  left_join(meta %>%
              dplyr::filter(day.inj == 1)) %>%
  mutate(exp = case_when(
    group == "unknown" ~ "unknown",
    TRUE ~ "controls"
  )) %>%
  ggplot(aes(x = temp, y = fluor)) +
  geom_line(aes(group = position.for.geno,
                colour = group),
            alpha = 0.5) +
  scale_x_continuous(limits = c(78,86),
                     breaks = seq(70,90)) +
  scale_y_continuous(labels = comma,
                     limits = c(0, 200000)) +
  facet_wrap(~exp, ncol = 1, scales = "free_y") +
  ggtitle("Normalised fluorescence")


# plot derivitive ----------------------------------------------------------------


# plot out the genotyping results
data.meltcurve %>%
  left_join(meta %>%
              dplyr::filter(day.inj == 1)) %>%
  mutate(exp = case_when(
    group == "unknown" ~ "unknown",
    TRUE ~ "controls"
  )) %>%
  ggplot(aes(x = temp, y = df)) +
  geom_line(aes(group = position.for.geno,
                colour = group),
            alpha = 0.5) +
  scale_x_continuous(limits = c(78,86),
                     breaks = seq(70,90)) +
  scale_y_continuous(labels = comma,
                     limits = c(0, 200000)) +

  # mutant allele
  annotate("rect", fill = "red", alpha = 0.25,
           xmin = 79.3, xmax = 81.7,
           ymin  = 0, ymax = 200000
           ) +

  annotate("text", colour = "red",
           label = "mut allele",
           x = 80,
           y  = 190000
  ) +


  annotate("rect", fill = "blue", alpha = 0.25,
           xmin = 82, xmax = 85,
           ymin  = 0, ymax = 200000
  ) +

  annotate("text", colour = "blue",
           label = "wt allele",
           x = 82.6,
           y  = 190000
  ) +
  facet_wrap(~exp, ncol = 1, scales = "free_y") +
  ggtitle("derivitive fluorescence")


# define the hets and homs
data.results %<>%
  mutate(
    across(c(Tm1,Tm2,Tm3), as.numeric)
  ) %>%
  mutate(
    mut.allele = case_when(
      between(Tm1, left = 79.3, right = 81.7) ~ TRUE,
      between(Tm2, left = 79.3, right = 81.7) ~ TRUE,
      TRUE ~ FALSE
    ),
    wt.allele = case_when(
      between(Tm1, left = 82, right = 85) ~ TRUE,
      between(Tm2, left = 82, right = 85) ~ TRUE,
      TRUE ~ FALSE
    ),
    genotype = case_when(
      mut.allele == TRUE & wt.allele == TRUE ~ "het",
      mut.allele == FALSE & wt.allele == TRUE ~ "wt",
      mut.allele == TRUE & wt.allele == FALSE ~ "hom"
    )
  )



# repeat the plot
data.results %>%
  dplyr::rename(position.for.geno = Well) %>%
    full_join(data.meltcurve) %>%
    mutate(exp = case_when(
      group == "unknown" ~ "unknown",
      TRUE ~ "controls"
  )) %>%
  ggplot(aes(x = temp, y = df)) +
  geom_line(aes(group = position.for.geno,
                colour = genotype),
            alpha = 0.5) +
  scale_x_continuous(limits = c(78,86),
                     breaks = seq(70,90)) +
  scale_y_continuous(labels = comma,
                     limits = c(0, 200000)) +
  # mutant allele
  annotate("rect", fill = "red", alpha = 0.25,
           xmin = 79.3, xmax = 81.7,
           ymin  = 0, ymax = 200000
  ) +

  annotate("text", colour = "red",
           label = "mut allele",
           x = 80,
           y  = 190000
  ) +


  annotate("rect", fill = "blue", alpha = 0.25,
           xmin = 82, xmax = 85,
           ymin  = 0, ymax = 200000
  ) +

  annotate("text", colour = "blue",
           label = "wt allele",
           x = 82.6,
           y  = 190000
  ) +
 facet_wrap(~exp, ncol = 1, scales = "free_y") +
  ggtitle("derivitive fluorescence")



# # export data for the final document
 meta %>%
   dplyr::filter(day.inj == 1) %>%
   left_join(data.results %>%
               dplyr::select(position.for.geno = Well, genotype),
             by  = "position.for.geno")  %>%
   dplyr::select(1,2,3, genotype = genotype.y, everything(), -genotype.x) %>%
   saveRDS("data/Full exp 2/R_objects/metadata_withGenotype_dayinj1.rds")
#
 data.meltcurve %>%
   left_join(data.results %>%
               dplyr::rename("position.for.geno" = "Well")) %>%
    saveRDS("data/Full exp 2/R_objects/meltcurveWithGenotypes_dayinj1.rds")

# day 2 -------------------------------------------------------------------



 # data import -------------------------------------------------------------

 # results file
 data.results.2 <-
   read_delim("data/Full exp 2/genotyping_data/KB 20230925 DFP genotyping naglufish 57to113_data.txt",
              delim = "\t",
              skip = 8) %>%
   dplyr::slice(1:61) %>%
   mutate(group = case_when(
     Well %in% c("A1", paste0("F", 9:12)) ~ "neg control",
     Well == "A2" ~ "wt control",
     Well == "A3" ~ "het control",
     Well == "A4" ~ "hom control",
     TRUE ~ "unknown"
   )) %>%
   dplyr::select(Well, group, Tm1, Tm2, Tm3)

 #temperatures
 temps.2 <-
   read_delim("data/Full exp 2/genotyping_data/KB 20230925 DFP genotyping naglufish 57to113_data.txt",
                     delim = "\t",
                     skip = 84)  %>%
   dplyr::slice(1:61) %>%
   mutate(group = case_when(
     `Well Location` %in% c("A1", paste0("F", 9:12)) ~ "neg control",
     `Well Location` == "A2" ~ "wt control",
     `Well Location` == "A3" ~ "het control",
     `Well Location` == "A4" ~ "hom control",
     TRUE ~ "unknown"
   )) %>%
   dplyr::select(position.for.geno = `Well Location`,
                 group,
                 all_of(starts_with("Reading")))


 # normailsed fluorescence
 fluors.2 <-
   read_delim("data/Full exp 2/genotyping_data/KB 20230925 DFP genotyping naglufish 57to113_data.txt",
              delim = "\t",
              skip = 154) %>%
   dplyr::slice(1:61) %>%
   mutate(group = case_when(
     `Well Location` %in% c("A1", paste0("F", 9:12)) ~ "neg control",
     `Well Location` == "A2" ~ "wt control",
     `Well Location` == "A3" ~ "het control",
     `Well Location` == "A4" ~ "hom control",
     TRUE ~ "unknown"
   )) %>%
   dplyr::select(position.for.geno = `Well Location`,
                 group,
                 all_of(starts_with("Reading")))

 # derivities
 dfs.2 <-
   read_delim("data/Full exp 2/genotyping_data/KB 20230925 DFP genotyping naglufish 57to113_data.txt",
              delim = "\t",
              skip = 224) %>%
   dplyr::slice(1:61) %>%
   mutate(group = case_when(
     `Well Location` %in% c("A1", paste0("F", 9:12)) ~ "neg control",
     `Well Location` == "A2" ~ "wt control",
     `Well Location` == "A3" ~ "het control",
     `Well Location` == "A4" ~ "hom control",
     TRUE ~ "unknown"
   )) %>%
   dplyr::select(position.for.geno = `Well Location`,
                 group,
                 all_of(starts_with("Reading")))

 # join together
 data.meltcurve <- dfs %>%
   pivot_longer(names_to = "Reading", values_to = "df",
                starts_with("Reading")
   ) %>%
   left_join(temps %>%
               pivot_longer(names_to = "Reading", values_to = "temp",
                            starts_with("Reading")
               )
   ) %>%
   left_join(fluors %>%
               pivot_longer(names_to = "Reading", values_to = "fluor",
                            starts_with("Reading") )
   ) %>%
   mutate(df = as.numeric(df),
          temp = as.numeric(temp),
          fluor = as.numeric(fluor),
          group = fct_inorder(group))


 # fluorescence
 data.meltcurve %>%
   left_join(meta %>%
               dplyr::filter(day.inj == 1)) %>%
   mutate(exp = case_when(
     group == "unknown" ~ "unknown",
     TRUE ~ "controls"
   )) %>%
   ggplot(aes(x = temp, y = fluor)) +
   geom_line(aes(group = position.for.geno,
                 colour = group),
             alpha = 0.5) +
   scale_x_continuous(limits = c(78,86),
                      breaks = seq(70,90)) +
   scale_y_continuous(labels = comma,
                      limits = c(0, 200000)) +
   facet_wrap(~exp, ncol = 1, scales = "free_y") +
   ggtitle("Normalised fluorescence")


 # plot derivitive ----------------------------------------------------------------


 # plot out the genotyping results
 data.meltcurve %>%
   left_join(meta %>%
               dplyr::filter(day.inj == 1)) %>%
   mutate(exp = case_when(
     group == "unknown" ~ "unknown",
     TRUE ~ "controls"
   )) %>%
   ggplot(aes(x = temp, y = df)) +
   geom_line(aes(group = position.for.geno,
                 colour = group),
             alpha = 0.5) +
   scale_x_continuous(limits = c(78,86),
                      breaks = seq(70,90)) +
   scale_y_continuous(labels = comma,
                      limits = c(0, 200000)) +

   # mutant allele
   annotate("rect", fill = "red", alpha = 0.25,
            xmin = 79.3, xmax = 81.7,
            ymin  = 0, ymax = 200000
   ) +

   annotate("text", colour = "red",
            label = "mut allele",
            x = 80,
            y  = 190000
   ) +


   annotate("rect", fill = "blue", alpha = 0.25,
            xmin = 82, xmax = 85,
            ymin  = 0, ymax = 200000
   ) +

   annotate("text", colour = "blue",
            label = "wt allele",
            x = 82.6,
            y  = 190000
   ) +
   facet_wrap(~exp, ncol = 1, scales = "free_y") +
   ggtitle("derivitive fluorescence")


 # define the hets and homs
 data.results %<>%
   mutate(
     across(c(Tm1,Tm2,Tm3), as.numeric)
   ) %>%
   mutate(
     mut.allele = case_when(
       between(Tm1, left = 79.3, right = 81.7) ~ TRUE,
       between(Tm2, left = 79.3, right = 81.7) ~ TRUE,
       TRUE ~ FALSE
     ),
     wt.allele = case_when(
       between(Tm1, left = 82, right = 85) ~ TRUE,
       between(Tm2, left = 82, right = 85) ~ TRUE,
       TRUE ~ FALSE
     ),
     genotype = case_when(
       mut.allele == TRUE & wt.allele == TRUE ~ "het",
       mut.allele == FALSE & wt.allele == TRUE ~ "wt",
       mut.allele == TRUE & wt.allele == FALSE ~ "hom"
     )
   )



 # repeat the plot
 data.results %>%
   dplyr::rename(position.for.geno = Well) %>%
   full_join(data.meltcurve) %>%
   mutate(exp = case_when(
     group == "unknown" ~ "unknown",
     TRUE ~ "controls"
   )) %>%
   ggplot(aes(x = temp, y = df)) +
   geom_line(aes(group = position.for.geno,
                 colour = genotype),
             alpha = 0.5) +
   scale_x_continuous(limits = c(78,86),
                      breaks = seq(70,90)) +
   scale_y_continuous(labels = comma,
                      limits = c(0, 200000)) +
   # mutant allele
   annotate("rect", fill = "red", alpha = 0.25,
            xmin = 79.3, xmax = 81.7,
            ymin  = 0, ymax = 200000
   ) +

   annotate("text", colour = "red",
            label = "mut allele",
            x = 80,
            y  = 190000
   ) +


   annotate("rect", fill = "blue", alpha = 0.25,
            xmin = 82, xmax = 85,
            ymin  = 0, ymax = 200000
   ) +

   annotate("text", colour = "blue",
            label = "wt allele",
            x = 82.6,
            y  = 190000
   ) +
   facet_wrap(~exp, ncol = 1, scales = "free_y") +
   ggtitle("derivitive fluorescence")



 # # export data for the final document
 meta %>%
   dplyr::filter(day.inj == 1) %>%
   left_join(data.results %>%
               dplyr::select(position.for.geno = Well, genotype),
             by  = "position.for.geno")  %>%
   dplyr::select(1,2,3, genotype = genotype.y, everything(), -genotype.x) %>%
   saveRDS("data/Full exp 2/R_objects/metadata_withGenotype_dayinj1.rds")
 #
 data.meltcurve %>%
   left_join(data.results %>%
               dplyr::rename("position.for.geno" = "Well")) %>%
   saveRDS("data/Full exp 2/R_objects/meltcurveWithGenotypes_dayinj1.rds")
