# -------------------------------------------
# Fuzzy AHP (3 experts) - full reproducible script
# -------------------------------------------

# Criteria labels (11)
criteria <- c("Slope","Elevation","LandCover","WaterBodies","LandUse",
              "Roads","SoilType","Dumpsite","PowerLine","Population","Railway")

# -------------------------------
# Expert Saaty matrices (11x11)
# -------------------------------

exp1 <- matrix(c(
  1,1,1,1,1,2,2,5,7,5,7,
  1,1,1,1,1,1,3,7,7,3,5,
  1,1,1,1,1,1,2,4,5,4,7,
  1,1,1,1,1,1,3,7,5,3,5,
  1,1,1,1,1,1,2,5,7,3,5,
  0.742,1,1,1,1,1,1,5,5,5,5,
  0.5,0.33,0.5,0.33,0.5,1,1,3,4,4,3,
  0.2,0.14,0.25,0.14,0.2,0.2,0.33,1,1,0.5,0.5,
  0.14,0.14,0.2,0.2,0.14,0.2,0.25,1,1,1,1,
  0.2,0.33,0.25,0.33,0.33,0.2,0.25,2,1,1,3,
  0.14,0.2,0.14,0.2,0.2,0.2,0.33,2,1,0.33,1
), nrow=11, byrow=TRUE)

exp2 <- matrix(c(
  1,1,1,1,1,1,2,4,7,5,7,
  1,1,1,1,1,1,2,7,6,3,5,
  1,1,1,1,1,1,2,4,5,4,7,
  1,1,1,1,1,1,2,7,5,3,5,
  1,1,1,1,1,1,2,4,7,3,5,
  0.742,1,1,1,1,1,1,5,4,5,5,
  0.5,0.5,0.5,0.5,0.5,1,1,3,4,4,3,
  0.25,0.14,0.25,0.14,0.25,0.2,0.33,1,1,0.5,0.5,
  0.14,0.17,0.2,0.2,0.14,0.25,0.25,1,1,1,1,
  0.2,0.33,0.25,0.33,0.33,0.2,0.25,2,1,1,2,
  0.14,0.2,0.14,0.2,0.2,0.2,0.33,2,1,0.5,1
), nrow=11, byrow=TRUE)

exp3 <- matrix(c(
  1,1,1,1,1,2,2,5,7,5,7,
  1,1,1,1,1,1,3,7,6,3,4,
  1,1,1,1,1,1,2,4,5,4,6,
  1,1,1,1,1,1,3,7,5,3,5,
  1,1,1,1,1,1,2,5,7,3,5,
  0.742,1,1,1,1,1,1,5,5,4,5,
  0.5,0.33,0.5,0.33,0.5,1,1,3,4,3,4,
  0.2,0.14,0.25,0.14,0.2,0.2,0.33,1,2,0.5,0.5,
  0.14,0.17,0.2,0.2,0.14,0.2,0.25,0.5,1,1,1,
  0.2,0.33,0.25,0.33,0.33,0.25,0.33,2,1,1,2,
  0.14,0.25,0.17,0.2,0.2,0.2,0.25,2,1,0.5,1
), nrow=11, byrow=TRUE)

# -------------------------------
# Fuzzy conversion (triangular l,m,u)
# -------------------------------
toFuzzy <- function(mat) {
  lapply(1:nrow(mat), function(i) lapply(1:ncol(mat), function(j) {
    m <- mat[i,j]
    l <- min(m, 1/m)
    u <- max(m, 1/m)
    c(l, m, u)
  }))
}

fuzzy_exp1 <- toFuzzy(exp1)
fuzzy_exp2 <- toFuzzy(exp2)
fuzzy_exp3 <- toFuzzy(exp3)

# -------------------------------
# Fuzzy geometric mean aggregation
# -------------------------------
fuzzy_agg <- lapply(1:11, function(i) lapply(1:11, function(j) {
  l <- (fuzzy_exp1[[i]][[j]][1] * fuzzy_exp2[[i]][[j]][1] * fuzzy_exp3[[i]][[j]][1])^(1/3)
  m <- (fuzzy_exp1[[i]][[j]][2] * fuzzy_exp2[[i]][[j]][2] * fuzzy_exp3[[i]][[j]][2])^(1/3)
  u <- (fuzzy_exp1[[i]][[j]][3] * fuzzy_exp2[[i]][[j]][3] * fuzzy_exp3[[i]][[j]][3])^(1/3)
  c(l, m, u)
}))
print(fuzzy_agg)

# -------------------------------
# Defuzzification (Yager) to crisp matrix
# -------------------------------
crisp_matrix <- matrix(0, 11, 11)
for (i in 1:11) for (j in 1:11) {
  tri <- fuzzy_agg[[i]][[j]]
  crisp_matrix[i, j] <- (tri[1] + 2*tri[2] + tri[3]) / 4
}
rownames(crisp_matrix) <- criteria
colnames(crisp_matrix) <- criteria
print(crisp_matrix)
# Enforce reciprocity and diagonal ones
for (i in 1:11) for (j in 1:11) if (i != j) crisp_matrix[j, i] <- 1 / crisp_matrix[i, j]
diag(crisp_matrix) <- 1

# -------------------------------
# Weights via principal eigenvector (Saaty)
# -------------------------------
eig <- eigen(crisp_matrix)
lambda_max <- max(Re(eig$values))
principal_vec <- Re(eig$vectors[, which.max(Re(eig$values))])
weights <- abs(principal_vec)
weights <- weights / sum(weights)
names(weights) <- criteria

# Consistency ratio (CR)
n <- nrow(crisp_matrix)
CI <- (lambda_max - n) / (n - 1)
RI <- 1.51 # RI for n=11
CR <- CI / RI

cat("Aggregated, defuzzified pairwise matrix:\n")
print(round(crisp_matrix, 3))
cat("\nWeights (sum=1):\n")
print(round(weights, 4))
cat("\nConsistency Ratio (CR): ", round(CR, 4), "\n")

# -------------------------------
# Site performance values (A, C, D, E, F, G)
# Replace these with your normalized decision matrix if needed
# -------------------------------
criteria_short <- c("WB","EV","SP","ST","RD","RL","LU","LC","PD","ED","PL")
# If your values are ordered differently, map accordingly.
# Below uses the values you shared earlier (6 sites x 11 criteria)
values <- c(
  0.286,0.286,0.143,0.071,0.238,0.048,0.250,0.167,0.238,0.095,0.095,
  0.238,0.238,0.286,0.071,0.286,0.095,0.250,0.167,0.190,0.048,0.048,
  0.190,0.190,0.048,0.214,0.143,0.143,0.125,0.167,0.286,0.143,0.190,
  0.095,0.143,0.238,0.214,0.048,0.238,0.125,0.167,0.048,0.238,0.238,
  0.143,0.048,0.095,0.214,0.190,0.286,0.125,0.167,0.143,0.286,0.286,
  0.048,0.095,0.190,0.214,0.095,0.190,0.125,0.167,0.095,0.190,0.143
)
values <- matrix(values, nrow = 6, byrow = TRUE)
rownames(values) <- c("A","C","D","E","F","G")
colnames(values) <- criteria_short

# Ensure weight alignment: reorder weights to criteria_short order if needed
# Create a named mapping from long to short labels
map <- c("WaterBodies"="WB","Elevation"="EV","Slope"="SP","SoilType"="ST",
         "Roads"="RD","Railway"="RL","LandUse"="LU","LandCover"="LC",
         "Population"="PD","Dumpsite"="ED","PowerLine"="PL")

# Build weights_short in same order as values columns
weights_short <- sapply(criteria_short, function(k) {
  long_name <- names(map[map == k])
  weights[long_name]
})
names(weights_short) <- criteria_short

cat("\nWeights aligned to site values (order: WB, EV, SP, ST, RD, RL, LU, LC, PD, ED, PL):\n")
print(round(weights_short, 4))

# -------------------------------
# Site scoring and ranking
# -------------------------------
scores <- as.vector(values %*% weights_short)
rank <- rank(-scores, ties.method = "first")  # higher score = better

fullresults <- cbind(values, score = round(scores, 4), rank = rank)
cat("\nSite scores and ranking:\n")
print(fullresults)



#SENSITIVITY ANALYSIS
# Required libraries
library(dplyr)
library(ggplot2)

# 1. FUZZY AHP WEIGHTS (from manuscript Table 3)
original_weights <- c(0.1311, 0.1331, 0.1402, 0.0848, 0.1193, 0.0309, 0.1267, 0.1274, 0.0429, 0.0330,
                      0.0306 )
criteria <- c("WB", "EV", "SP", "ST", "RD", "RL", "LU", "LC", "PD", "ED", "PL")

# 2. NORMALIZED DECISION MATRIX (Table 5) - Sites A, C, D, E, F, G
decision_matrix <- matrix(c(
  0.286,0.286,0.143,0.071,0.238,0.048,0.250,0.167,0.238,0.095,0.095,  # A
  0.238,0.238,0.286,0.071,0.286,0.095,0.250,0.167,0.190,0.048,0.048,  # C
  0.190,0.190,0.048,0.214,0.143,0.143,0.125,0.167,0.286,0.143,0.190,  # D
  0.095,0.143,0.238,0.214,0.048,0.238,0.125,0.167,0.048,0.238,0.238,  # E
  0.143,0.048,0.095,0.214,0.190,0.286,0.125,0.167,0.143,0.286,0.286,  # F
  0.048,0.095,0.190,0.214,0.095,0.190,0.125,0.167,0.095,0.190,0.143   # G
), nrow = 6, byrow = TRUE)

rownames(decision_matrix) <- c("A","C","D","E","F","G")
colnames(decision_matrix) <- criteria

# 3. OAT SENSITIVITY FUNCTION
oat_sensitivity_analysis <- function(weights, decision_matrix, criteria, 
                                     change_levels = seq(-1, 2, by = 0.1)) {
  
  n_criteria <- length(weights)
  n_sites <- nrow(decision_matrix)
  results <- data.frame()
  
  # Baseline ranking
  baseline_scores <- as.vector(decision_matrix %*% weights)
  baseline_order <- order(baseline_scores, decreasing = TRUE)
  site_names <- rownames(decision_matrix)
  baseline_ranks <- setNames(match(seq_len(n_sites), baseline_order), site_names)
  
  c_index <- which(site_names == "C")
  c_baseline_rank <- as.integer(baseline_ranks["C"])
  
  cat("BASELINE RANKING:\n")
  print(data.frame(Site = site_names, Score = round(baseline_scores, 4), Rank = baseline_ranks))
  cat(sprintf("\nSite C baseline rank: %d\n\n", c_baseline_rank))
  
  for (i in seq_len(n_criteria)) {
    cat("=== ANALYZING", criteria[i], "===\n")
    
    for (change in change_levels) {
      new_weights <- weights
      new_weights[i] <- new_weights[i] * (1 + change)
      new_weights <- new_weights / sum(new_weights)
      
      new_scores <- as.vector(decision_matrix %*% new_weights)
      new_order <- order(new_scores, decreasing = TRUE)
      ranks <- match(seq_len(n_sites), new_order)
      c_rank <- ranks[c_index]
      
      results <- rbind(
        results,
        data.frame(
          Criterion = criteria[i],
          Change_pct = round(change * 100, 1),
          New_Weight = round(new_weights[i], 4),
          Site = site_names,
          Score = round(new_scores, 4),
          Rank = ranks,
          C_Rank_Stable = ifelse(site_names == "C" & c_rank == c_baseline_rank, "YES", "NO"),
          stringsAsFactors = FALSE
        )
      )
      
      if (c_rank != c_baseline_rank) {
        cat(sprintf("?????? CRITICAL: %s %+d%% ??? Site C rank %d (baseline %d)\n",
                    criteria[i], round(change * 100), c_rank, c_baseline_rank))
      }
    }
  }
  return(list(
    results = results,
    c_baseline_rank = c_baseline_rank
  ))
}

# 4. RUN OAT ANALYSIS
cat("Running OAT Sensitivity Analysis...\n")
oat <- oat_sensitivity_analysis(original_weights, decision_matrix, criteria)
sensitivity_results <- oat$results
c_baseline_rank <- oat$c_baseline_rank

# 5. CRITICALITY SUMMARY
criticality_summary <- sensitivity_results %>%
  filter(Site == "C") %>%
  group_by(Criterion) %>%
  summarise(
    Min_Change_To_Reverse_C = if (any(Rank > c_baseline_rank)) {
      min(Change_pct[Rank > c_baseline_rank])
    } else {
      NA_real_
    },
    Stability_pct = mean(C_Rank_Stable == "YES") * 100,
    Sensitivity_Index = ifelse(is.na(Min_Change_To_Reverse_C), NA_real_, 1 / abs(Min_Change_To_Reverse_C)),
    .groups = "drop"
  ) %>%
  arrange(desc(Sensitivity_Index))

cat("\n=== OAT SENSITIVITY SUMMARY ===\n")
print(criticality_summary)

# 6. VISUALIZATION
p1 <- ggplot(sensitivity_results %>% filter(Site == "C"), 
             aes(x = Change_pct, y = Rank, color = Criterion, group = Criterion)) +
  geom_line(size = 1.2) + 
  geom_point(size = 3) +
  labs(title = "OAT Sensitivity: Site C Rank Stability",
       subtitle = sprintf("Lower rank = better (baseline rank = %d)", c_baseline_rank),
       x = "Weight Change (%)", y = "Site C Rank") +
  scale_y_reverse(breaks = 1:nrow(decision_matrix)) +
  theme_minimal() + 
  theme(legend.position = "bottom")

p2 <- ggplot(criticality_summary, aes(x = reorder(Criterion, Sensitivity_Index, na.rm = TRUE), 
                                      y = Sensitivity_Index)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = "Sensitivity Ranking (higher = more sensitive)",
       subtitle = "Minimum change that worsens Site C's rank",
       x = "Criterion", y = "Sensitivity Index") +
  theme_minimal()

print(p1)
print(p2)

# 7. EXPORT RESULTS
write.csv(sensitivity_results, "oat_sensitivity_results.csv", row.names = FALSE)
write.csv(criticality_summary, "criticality_summary.csv", row.names = FALSE)

cat("\n??? Analysis complete! Files exported.\n")



# -------------------------------------------
# EXTENSION: Visualizations for Fuzzy-AHP
# -------------------------------------------

library(ggplot2)
library(fmsb)

# 1. Criterion Weights Bar Chart
weights_df <- data.frame(
  Criterion = names(weights_short),
  Weight = weights_short
)

ggplot(weights_df, aes(x = reorder(Criterion, Weight), y = Weight)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = "Fuzzy AHP Criterion Weights",
       x = "Criterion", y = "Weight (normalized)") +
  theme_minimal()

# 2. Site Scores Bar Chart
scores_df <- data.frame(
  Site = rownames(values),
  Score = scores,
  Rank = rank
)

ggplot(scores_df, aes(x = reorder(Site, Score), y = Score, fill = Site)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = paste0("Rank ", Rank)), vjust = -0.5) +
  labs(title = "Landfill Site Scores (Fuzzy AHP)",
       x = "Site", y = "Weighted Score") +
  theme_minimal()




##RADAR CHARTS


library(fmsb)

# Prepare data: decision_matrix already defined (sites A-G, criteria WB-PL)
max_vals <- apply(decision_matrix, 2, max)
min_vals <- apply(decision_matrix, 2, min)

# fmsb requires first two rows as max and min
radar_data <- rbind(max_vals, min_vals, decision_matrix)

# Convert to data frame
radar_df <- as.data.frame(radar_data)

# Colors for sites
colors <- c("red","blue","green","purple","orange","brown")
site_names <- rownames(decision_matrix)

# --- Individual radar charts ---
par(mfrow = c(2,3))  # arrange plots in grid
for(i in 1:nrow(decision_matrix)) {
  radarchart(radar_df[c(1,2,i+2), ],
             axistype = 1,
             pcol = colors[i],
             pfcol = adjustcolor(colors[i], alpha.f = 0.3),
             plwd = 2,
             title = paste("Site", site_names[i]))
}
par(mfrow = c(1,1))  # reset layout

# --- Combined radar chart ---
radarchart(radar_df,
           axistype = 1,
           pcol = colors,
           pfcol = scales::alpha(colors, 0.3),
           plwd = 2,
           plty = 1,
           title = "Combined Radar Chart: Site Profiles")

legend("topright", legend = site_names,
       col = colors, lty = 1, lwd = 2, bty = "n")









# -------------------------------------------
# EXTENSION: Sensitivity Plots for All Sites
# -------------------------------------------

library(ggplot2)
library(dplyr)

# 1. Line plots: rank stability for all sites
ggplot(sensitivity_results,
       aes(x = Change_pct, y = Rank, color = Site, group = Site)) +
  geom_line(size = 1.2) +
  facet_wrap(~Criterion, scales = "free_x") +
  scale_y_reverse(breaks = 1:nrow(decision_matrix)) +
  labs(title = "OAT Sensitivity: Rank Stability for All Sites",
       subtitle = "Lower rank = better",
       x = "Weight Change (%)", y = "Rank") +
  theme_minimal() +
  theme(legend.position = "bottom")

# 2. Sensitivity summary for all sites
criticality_summary_all <- sensitivity_results %>%
  group_by(Site, Criterion) %>%
  summarise(
    Min_Change_To_Worsen = if (any(Rank > min(Rank))) {
      min(Change_pct[Rank > min(Rank)])
    } else { NA_real_ },
    Stability_pct = mean(Rank == min(Rank)) * 100,
    Sensitivity_Index = ifelse(is.na(Min_Change_To_Worsen), NA_real_, 1 / abs(Min_Change_To_Worsen)),
    .groups = "drop"
  )

# 3. Bar chart: sensitivity indices by site
ggplot(criticality_summary_all,
       aes(x = reorder(Criterion, Sensitivity_Index, na.rm = TRUE),
           y = Sensitivity_Index, fill = Site)) +
  geom_col(position = "dodge") +
  coord_flip() +
  labs(title = "Sensitivity Indices by Site and Criterion",
       subtitle = "Higher = more sensitive to weight changes",
       x = "Criterion", y = "Sensitivity Index") +
  theme_minimal()
