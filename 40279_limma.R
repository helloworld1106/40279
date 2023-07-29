install.packages("remotes")
install.packages("ggplot2")
install.packages("patchwork")
install.packages("dplyr")
library(patchwork)
library(dplyr)
# Download modified GEOquery package from my github 
# by using function(install_github()) from 'remotes' package. 
library(remotes)
library(limma)
library(ggplot2)

install_github("curryhank08/GEOquery_with_modifiable_timeout_seconds", force = TRUE)
# Load modified GEOquery
library(GEOquery)
# Setting the max timeout_seconds
options(timeout=300000)
# Check the input timeout_seconds
getOption("timeout")
# Download GSE40279 by a fuction getGEO() from modified GEOquery package.
gse40279 <- getGEO("GSE40279", GSEMatrix = TRUE, AnnotGPL = TRUE)
gset <- gse40279


if (length(gset) > 1) idx <- grep("GPL13534", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
# Create age categories
age <- pData(gset)$characteristics_ch1
# Remove "age (y):" and convert to numeric
age <- sub("^\\s*age \\(y\\): ", "", age)
age <- as.numeric(age)
# Assign age values to a new column in pData of gset
pData(gset)$age <- age
# Define age categories based on specific age ranges
age_categories <- cut(age,
                      breaks = c(0, 30, 65, Inf),
                      labels = c("Young", "Middle", "Old"),
                      include.lowest = TRUE)
# Assign age categories to the pData of gset
pData(gset)$age_category <- age_categories
ex <- exprs(gset)
# log2 transform
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
ex <- log2(ex)}
conditions <- gset$age_category
f <- factor(conditions, levels = c("Young", "Middle", "Old"))
design <- model.matrix(~0+f)
colnames(design) <- c("Young", "Middle", "Old")
fit <- lmFit(gset, design)
contrast.matrix <- makeContrasts(Young-Middle, Young-Old, Middle-Old, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

result_limma_YM <- topTable(fit2, coef=1, number = Inf, adjust.method = "BH", sort.by = "logFC")
result_limma_YO <- topTable(fit2, coef=2, number = Inf, adjust.method = "BH", sort.by = "logFC")
result_limma_MO <- topTable(fit2, coef=3, number = Inf, adjust.method = "BH", sort.by = "logFC")

# Outcome of each hypothesis test
results <- decideTests(fit2,method="separate",adjust.method="BH",p.value=1e-5)

# Showing numbers of genes significant in each comparison
vennDiagram(results)

# Remove rows with missing values
result_limma_coef1_clean <- na.omit(result_limma_2)

# Get the row "cg00000029" from ex
cg00000029_row <- ex["cg00000029", ]

# Add the row as a new column to pData of gset
pData(gset)$cg00000029 <- cg00000029_row

# Create a boxplot for the expression levels of cg00000029 between "Young" and "Middle" age categories
box_plot_cg00000029 <- ggplot(pData(gset), aes(x = age_category, y = cg00000029)) +
  geom_boxplot() +
  labs(x = "Age Category", y = "Expression Level of cg00000029", title = "Expression of cg00000029 (Y-M)") +
  theme_minimal()

# Display the boxplot for cg00000029
print(box_plot_cg00000029)

# Load the ggplot2 library
library(ggplot2)

# Plot histogram of p-values
ggplot(result_limma_YO, aes(x = adj.P.Val)) +
  geom_histogram(binwidth = 0.01, color = "black", fill = "lightblue") +
  labs(x = "adj P-value", y = "Frequency", title = "Histogram of P-values") +
  theme_minimal()

# Create a volcano plot
plot_YM <- ggplot(result_limma_YM, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(aes(color = ifelse(-log10(P.Value) > -log10(0.05), "red", "black"))) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  labs(x = "Log Fold Change", y = "-log10(P-value)") +
  ggtitle("Volcano Plot: Young vs. Middle") +
  theme_minimal()
plot_YO <- ggplot(result_limma_YO, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(aes(color = ifelse(-log10(P.Value) > -log10(0.05), "red", "black"))) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  labs(x = "Log Fold Change", y = "-log10(P-value)") +
  ggtitle("Volcano Plot: Young vs. Old") +
  theme_minimal()
plot_MO <- ggplot(result_limma_MO, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(aes(color = ifelse(-log10(P.Value) > -log10(0.05), "red", "black"))) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  labs(x = "Log Fold Change", y = "-log10(P-value)") +
  ggtitle("Volcano Plot: Middle vs. Old") +
  theme_minimal()
# Display the plots using print()
print(plot_YM)
print(plot_YO)
print(plot_MO)

print("Summary of Differential Expression Tests:")
print(table(results))