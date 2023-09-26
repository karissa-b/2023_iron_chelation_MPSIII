glm_timebyturnsFeCit <- read_rds("~/Documents/2022_MPSII_qPCR_ironSupplementation/data/FeCit_fullexp_glmtimeBytotalturns.rds")

glm_altsByturnsFeCit <- read_rds("~/Documents/2022_MPSII_qPCR_ironSupplementation/output/R_objects/FeCit_fullexp.rds") %>%
  mutate(L_R_bias = factor(L_R_bias,
                           levels = c("Neither", "Left", "Right")),
         non_alts = total_turns - alts,
         bin = as.factor(bin)
) %>%
  glmmTMB(
    cbind(alts, non_alts) ~ (bin + genotype + treatment)^3 + L_R_bias + (1|start.time) + (1|fish_id),
    family = betabinomial(),
    data = .
   )

glm_alts_salbutamol <- read_rds("~/Documents/2023_MPSIII_salbutamol/output/R_objects/pilot/glmAlts7daysSalbutamol.rds")

glm_turnsSalbutamol <- read_rds("~/Documents/2023_MPSIII_salbutamol/output/R_objects/pilot/glmtotalTurnss7daysSalbutamol.rds")

theme_set(theme_classic())

ggarrange(

  print(emmeans(glm_turnsSalbutamol, ~ genotype * treatment), type ="response") %>%
    as_tibble() %>%
    dplyr::filter(treatment != "100 µM Salbutamol") %>%
    ggplot(aes(x = treatment, y = response, colour = genotype)) +
    geom_col(aes(fill = genotype),
             alpha = 0.75,
             width = 0.75,
             position = position_dodge()) +
    geom_errorbar(aes(ymin = asymp.LCL, ymax =asymp.UCL),
                  width = 0.125,
                  size = 1,
                  position = position_dodge(width = 0.75)) +
    scale_color_manual(values = c("#7F7F81", "#E6C6CD")) +
    scale_fill_manual(values = c("#7F7F81", "#E6C6CD")) +
    labs(y = "Average number of\narm entries per 10\nminutes in a Y-maze",
         x = " ")+
    theme(plot.background = element_blank(),
          panel.background =  element_blank(),
          legend.background =  element_blank(),
          panel.border =  element_blank(),
          text = element_text(colour = "white"),
          axis.text = element_text(colour = "white")
    ),

  print(emmeans(glm_alts_salbutamol, ~ genotype * treatment), type = "response") %>%
    as_tibble() %>%
    dplyr::filter(treatment != "100 µM Salbutamol") %>%
    ggplot(aes(x = treatment, y = prob, colour = genotype)) +
    geom_col(aes(fill = genotype),
             alpha = 0.75,
             width = 0.75,
             position = position_dodge()) +
    geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL ),
                  width = 0.125,
                  size = 1,
                  position = position_dodge(width = 0.75)) +
    scale_color_manual(values = c("#7F7F81", "#E6C6CD")) +
    scale_fill_manual(values = c("#7F7F81", "#E6C6CD")) +
    labs(y = "Estimated probability of\nalternation (LRLR + RLRL)",
         x = " " )+
    theme(plot.background = element_blank(),
          panel.background =  element_blank(),
          legend.background =  element_blank(),
          panel.border =  element_blank(),
          text = element_text(colour = "white"),
          axis.text = element_text(colour = "white")
    ),

  print(emmeans(glm_timebyturnsFeCit, ~ genotype * treatment), type ="response") %>%
    as_tibble() %>%
    ggplot(aes(x = treatment, y = response, colour = genotype)) +
    geom_col(aes(fill = genotype),
             alpha = 0.75,
             width = 0.75,
             position = position_dodge()) +
    geom_errorbar(aes(ymin = asymp.LCL, ymax =asymp.UCL),
                  width = 0.125,
                  size = 1,
                  position = position_dodge(width = 0.75)) +
    scale_color_manual(values = c("#7F7F81", "#E6C6CD")) +
    scale_fill_manual(values = c("#7F7F81", "#E6C6CD")) +
    labs(y = "Average number of\narm entries per 10\nminutes in a Y-maze",
         x = " ")+
    theme(plot.background = element_blank(),
          panel.background =  element_blank(),
          legend.background =  element_blank(),
          panel.border =  element_blank(),
          text = element_text(colour = "white"),
          axis.text = element_text(colour = "white")
    ),

  print(emmeans(glm_altsByturnsFeCit, ~ genotype * treatment), type = "response") %>%
    as_tibble() %>%
    ggplot(aes(x = treatment, y = prob, colour = genotype)) +
    geom_col(aes(fill = genotype),
             alpha = 0.75,
             width = 0.75,
             position = position_dodge()) +
    geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL ),
                  width = 0.125,
                  size = 1,
                  position = position_dodge(width = 0.75)) +
scale_color_manual(values = c("#7F7F81", "#E6C6CD")) +
  scale_fill_manual(values = c("#7F7F81", "#E6C6CD")) +
    labs(y = "Estimated probability of\nalternation (LRLR + RLRL)",
         x = " " )+
    theme(plot.background = element_blank(),
          panel.background =  element_blank(),
          legend.background =  element_blank(),
          panel.border =  element_blank(),
          text = element_text(colour = "white"),
          axis.text = element_text(colour = "white")
          ),



  print(emmeans(glm_timebyturns.males, ~ genotype * treatment), type ="response") %>%
    as_tibble() %>%
    mutate(treatment = case_when(
      grepl(treatment, pattern = "DFP") ~ "7 µg Deferiprone",
      TRUE ~ treatment)
    ) %>%
    ggplot(aes(x = treatment, y = response, colour = genotype)) +
    geom_col(aes(fill = genotype),
             alpha = 0.75,
             width = 0.75,
             position = position_dodge()) +
    geom_errorbar(aes(ymin = asymp.LCL, ymax =asymp.UCL),
                  width = 0.125,
                  size = 1,
                  position = position_dodge(width = 0.75)) +
scale_color_manual(values = c("#7F7F81", "#E6C6CD")) +
  scale_fill_manual(values = c("#7F7F81", "#E6C6CD")) +
    labs(y = "Average number of\narm entries per 10\nminutes in a Y-maze",
         x = " ") +
    theme(plot.background = element_blank(),
          panel.background =  element_blank(),
          legend.background =  element_blank(),
          panel.border =  element_blank(),
          text = element_text(colour = "white"),
          axis.text = element_text(colour = "white")
          ),

  print(emmeans(glm.alts.males.only, ~ genotype * treatment), type = "response") %>%
    as_tibble() %>%
    mutate(treatment = case_when(
      grepl(treatment, pattern = "DFP") ~ "7 µg Deferiprone",
      TRUE ~ treatment)
    ) %>%
    ggplot(aes(x = treatment, y = prob, colour = genotype)) +
    geom_col(aes(fill = genotype),
             alpha = 0.75,
             width = 0.75,
             position = position_dodge()) +
    geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL ),
                  width = 0.125,
                  size = 1,
                  position = position_dodge(width = 0.75)) +
  scale_color_manual(values = c("#7F7F81", "#E6C6CD")) +
scale_fill_manual(values = c("#7F7F81", "#E6C6CD")) +
    labs(y = "Estimated probability of\nalternation (LRLR + RLRL)",
         x = " ")+
    theme(plot.background = element_blank(),
          panel.background =  element_blank(),
          legend.background =  element_blank(),
          panel.border =  element_blank(),
          text = element_text(colour = "white"),
          axis.text = element_text(colour = "white")
    ),
  common.legend = T,
legend = "bottom",
nrow = 1)
ggsave("output/behavfirstpassfigrow.png",
       width = 35, height = 8,
       units = "cm",
       dpi = 400,
       scale = 1.3)
    )




~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



ggarrange(

    print(emmeans(glm_turnsSalbutamol, ~ genotype * treatment), type ="response") %>%
      as_tibble() %>%
      dplyr::filter(treatment != "100 µM Salbutamol") %>%
      ggplot(aes(x = treatment, y = response, colour = genotype)) +
      geom_col(aes(fill = genotype),
               alpha = 0.75,
               width = 0.75,
               position = position_dodge()) +
      geom_errorbar(aes(ymin = asymp.LCL, ymax =asymp.UCL),
                    width = 0.125,
                    size = 1,
                    position = position_dodge(width = 0.75)) +
      scale_color_manual(values = c("#7F7F81", "#E6C6CD")) +
      scale_fill_manual(values = c("#7F7F81", "#E6C6CD")) +
      labs(y = "Average number of\narm entries per 10\nminutes in a Y-maze",
           x = " ")+
      theme(plot.background = element_blank(),
            panel.background =  element_blank(),
            legend.background =  element_blank(),
            panel.border =  element_blank(),
            text = element_text(colour = "black", size = 15),
            axis.text = element_text(colour = "black", face = "bold")
      ) +
      ggtitle("p = 0.9"),

    print(emmeans(glm_alts_salbutamol, ~ genotype * treatment), type = "response") %>%
      as_tibble() %>%
      dplyr::filter(treatment != "100 µM Salbutamol") %>%
      ggplot(aes(x = treatment, y = prob, colour = genotype)) +
      geom_col(aes(fill = genotype),
               alpha = 0.75,
               width = 0.75,
               position = position_dodge()) +
      geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL ),
                    width = 0.125,
                    size = 1,
                    position = position_dodge(width = 0.75)) +
      scale_color_manual(values = c("#7F7F81", "#E6C6CD")) +
      scale_fill_manual(values = c("#7F7F81", "#E6C6CD")) +
      labs(y = "Estimated probability\nof alternation\n(LRLR + RLRL)",
           x = " " )+
      theme(plot.background = element_blank(),
            panel.background =  element_blank(),
            legend.background =  element_blank(),
            panel.border =  element_blank(),
            text = element_text(colour = "black", size = 15),
            axis.text = element_text(colour = "black", face = "bold")
      ) +
      ggtitle("p = 0.08"),

    print(emmeans(glm_timebyturnsFeCit, ~ genotype * treatment), type ="response") %>%
      as_tibble() %>%
      ggplot(aes(x = treatment, y = response, colour = genotype)) +
      geom_col(aes(fill = genotype),
               alpha = 0.75,
               width = 0.75,
               position = position_dodge()) +
      geom_errorbar(aes(ymin = asymp.LCL, ymax =asymp.UCL),
                    width = 0.125,
                    size = 1,
                    position = position_dodge(width = 0.75)) +
      scale_color_manual(values = c("#7F7F81", "#E6C6CD")) +
      scale_fill_manual(values = c("#7F7F81", "#E6C6CD")) +
      labs(y = "Average number of\narm entries per 10\nminutes in a Y-maze",
           x = " ")+
      theme(plot.background = element_blank(),
            panel.background =  element_blank(),
            legend.background =  element_blank(),
            panel.border =  element_blank(),
            text = element_text(colour = "black", size = 15),
            axis.text = element_text(colour = "black", face = "bold")
      ) +
      ggtitle("p = 0.03"),

    print(emmeans(glm_altsByturnsFeCit, ~ genotype * treatment), type = "response") %>%
      as_tibble() %>%
      ggplot(aes(x = treatment, y = prob, colour = genotype)) +
      geom_col(aes(fill = genotype),
               alpha = 0.75,
               width = 0.75,
               position = position_dodge()) +
      geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL ),
                    width = 0.125,
                    size = 1,
                    position = position_dodge(width = 0.75)) +
      scale_color_manual(values = c("#7F7F81", "#E6C6CD")) +
      scale_fill_manual(values = c("#7F7F81", "#E6C6CD")) +
      labs(y = "Estimated probability\nof alternation\n(LRLR + RLRL)",
           x = " " )+
      theme(plot.background = element_blank(),
            panel.background =  element_blank(),
            legend.background =  element_blank(),
            panel.border =  element_blank(),
            text = element_text(colour = "black", size = 15),
            axis.text = element_text(colour = "black", face = "bold")
      )+
      ggtitle("p = 0.1"),



    print(emmeans(glm_timebyturns.males, ~ genotype * treatment), type ="response") %>%
      as_tibble() %>%
      mutate(treatment = case_when(
        grepl(treatment, pattern = "DFP") ~ "7 µg Deferiprone",
        TRUE ~ treatment)
      ) %>%
      ggplot(aes(x = treatment, y = response, colour = genotype)) +
      geom_col(aes(fill = genotype),
               alpha = 0.75,
               width = 0.75,
               position = position_dodge()) +
      geom_errorbar(aes(ymin = asymp.LCL, ymax =asymp.UCL),
                    width = 0.125,
                    size = 1,
                    position = position_dodge(width = 0.75)) +
      scale_color_manual(values = c("#7F7F81", "#E6C6CD")) +
      scale_fill_manual(values = c("#7F7F81", "#E6C6CD")) +
      labs(y = "Average number of\narm entries per 10\nminutes in a Y-maze",
           x = " ") +
      theme(plot.background = element_blank(),
            panel.background =  element_blank(),
            legend.background =  element_blank(),
            panel.border =  element_blank(),
            text = element_text(colour = "black", size = 15),
            axis.text = element_text(colour = "black", face = "bold")
      )+
      ggtitle("p = 0.4"),

    print(emmeans(glm.alts.males.only, ~ genotype * treatment), type = "response") %>%
      as_tibble() %>%
      mutate(treatment = case_when(
        grepl(treatment, pattern = "DFP") ~ "7 µg Deferiprone",
        TRUE ~ treatment)
      ) %>%
      ggplot(aes(x = treatment, y = prob, colour = genotype)) +
      geom_col(aes(fill = genotype),
               alpha = 0.75,
               width = 0.75,
               position = position_dodge()) +
      geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL ),
                    width = 0.125,
                    size = 1,
                    position = position_dodge(width = 0.75)) +
      scale_color_manual(values = c("#7F7F81", "#E6C6CD")) +
      scale_fill_manual(values = c("#7F7F81", "#E6C6CD")) +
      labs(y = "Estimated probability\nof alternation\n(LRLR + RLRL)",
           x = " ")+
      theme(plot.background = element_blank(),
            panel.background =  element_blank(),
            legend.background =  element_blank(),
            panel.border =  element_blank(),
            text = element_text(colour = "black", size = 15),
            axis.text = element_text(colour = "black", face = "bold")
      )+
      ggtitle("p = 0.9"),
    common.legend = T,
    legend = "bottom",
    ncol = 2,
    nrow = 3) +
ggsave("output/behavfirstpassfigrowBlack.png",
       width = 20, height = 22,
       units = "cm",
       dpi = 400,
       scale = 1)
