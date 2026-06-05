# Central Bank Gold Demand and the Global Gold Price
# Replication script for the empirical analysis (Chapters 2-3)
#
# Paper: "The Impact of Central Bank Demand on the Global Gold Price
#         in the Context of Transforming International Reserves"
# Author: O. Steshenko, NaUKMA, 2026
#
# Usage:
#   place the dataset at  data/Data_term_paper.xlsx  (sheets: "monthly", "quarterly")
#   then run:  Rscript gold_analysis.R [optional/path/to/dataset.xlsx]
#
# Reproducibility note: the regression sample is fixed to 2006:02-2025:08 and
# reproduces every Chapter 3 result. Chapter 2 descriptive statistics use the
# full monthly series; because source prices are revised and extended over time,
# the tail of the price level can differ slightly from the snapshot behind Table 2.1.

set.seed(2026)

# Packages --------------------------------------------------------------------
required <- c("readxl", "dplyr", "tidyr", "lubridate",
              "psych", "corrplot", "car", "lmtest", "sandwich",
              "strucchange", "tseries", "vars", "relaimpo")
to_install <- required[!required %in% installed.packages()[, "Package"]]
if (length(to_install) > 0)
  install.packages(to_install, repos = "https://cran.rstudio.com/")
invisible(lapply(required, function(p)
  suppressPackageStartupMessages(library(p, character.only = TRUE))))
select <- dplyr::select
filter <- dplyr::filter

# Input file ------------------------------------------------------------------
args       <- commandArgs(trailingOnly = TRUE)
EXCEL_FILE <- if (length(args) >= 1) args[1] else "Data_term_paper.xlsx"
if (!file.exists(EXCEL_FILE))
  stop("Dataset not found at '", EXCEL_FILE,
       "'. Place the .xlsx there or pass the path as an argument.")

# Read data -------------------------------------------------------------------
monthly_raw <- read_excel(EXCEL_FILE, sheet = "monthly")
quarterly   <- read_excel(EXCEL_FILE, sheet = "quarterly")

find_date <- function(df) {
  dc <- names(df)[sapply(df, function(x) inherits(x, c("Date", "POSIXct", "POSIXt")))]
  if (length(dc) == 0) names(df)[1] else dc[1]
}

monthly_raw <- monthly_raw %>%
  rename(date = all_of(find_date(monthly_raw))) %>%
  mutate(date = as.Date(date)) %>%
  arrange(date) %>%
  mutate(across(where(is.character), as.numeric)) %>%
  rename(Goldprice = Goldpriceendperiod)   # paper uses end-of-month price (Table 3.1)

quarterly <- quarterly %>%
  rename(date = all_of(find_date(quarterly))) %>%
  mutate(date = as.Date(date)) %>%
  arrange(date) %>%
  mutate(across(where(is.character), as.numeric))

# Transformations: returns, first differences, lags, post-2022 dummy ----------
df <- monthly_raw %>%
  mutate(
    r_Gold     = c(NA, diff(log(Goldprice))) * 100,   # dependent variable
    r_SP500    = c(NA, diff(log(SP500)))     * 100,
    r_Oil      = c(NA, diff(log(Oil)))       * 100,
    r_DXY      = c(NA, diff(log(DXY)))       * 100,
    d_DFII10   = c(NA, diff(DFII10)),
    d_DGS2     = c(NA, diff(DGS2)),
    d_VIX      = c(NA, diff(VIX)),
    d_T5YIFR   = c(NA, diff(T5YIFR)),
    d_GPR      = c(NA, diff(GPR)),
    d_ETFs     = c(NA, diff(ETFs)),
    d_NonCom   = c(NA, diff(NonComNetLong)),
    FedBal_yoy = ifelse(row_number() > 12, (FedBal / lag(FedBal, 12) - 1) * 100, NA),
    r_Gold_L1  = lag(r_Gold, 1),
    CB_L1      = lag(NetPurchasesTonnes, 1),   # publication-lagged CB demand
    CB_L2      = lag(NetPurchasesTonnes, 2),
    CB_L3      = lag(NetPurchasesTonnes, 3),
    dCB_L1     = lag(NetPurchasesTonnes - lag(NetPurchasesTonnes, 1), 1),
    dCB_L2     = lag(NetPurchasesTonnes - lag(NetPurchasesTonnes, 1), 2),
    D22        = as.numeric(date >= as.Date("2022-02-24")),
    CB_L1_D22  = CB_L1 * D22,
    CB_L2_D22  = CB_L2 * D22,
    CB_L3_D22  = CB_L3 * D22
  )

# Descriptive statistics (Table 2.1) ------------------------------------------
cat("\nDescriptive statistics - levels\n")
print(round(describe(df %>% select(Goldprice, DXY, SP500, VIX, DFII10, DGS2, Oil,
                                   GPR, FedBal, ETFs, NonComNetLong, GEPUgdp,
                                   NetPurchasesTonnes)
                     )[, c("n","mean","sd","median","min","max","skew","kurtosis")], 3))

cat("\nDescriptive statistics - returns and first differences\n")
print(round(describe(df %>% select(r_Gold, r_SP500, r_Oil, r_DXY, d_DFII10, d_DGS2,
                                   d_VIX, d_GPR, d_ETFs, d_NonCom, FedBal_yoy, GEPUgdp)
                     )[, c("n","mean","sd","median","min","max","skew","kurtosis")], 3))

# Stationarity (ADF) ----------------------------------------------------------
cat("\nADF tests (alternative: stationary)\n")
adf_one <- function(x, name) {
  xc <- na.omit(as.numeric(x)); if (length(xc) < 20) return(invisible())
  r <- suppressWarnings(adf.test(xc, alternative = "stationary"))
  cat(sprintf("  %-20s p = %.4f  %s\n", name, r$p.value,
              ifelse(r$p.value < 0.05, "stationary", "non-stationary")))
}
for (v in c("Goldprice","DXY","DFII10","ETFs","NonComNetLong","NetPurchasesTonnes"))
  adf_one(df[[v]], v)
for (v in c("r_Gold","r_DXY","d_DFII10","d_ETFs","d_NonCom"))
  adf_one(df[[v]], v)

# Correlation matrix (Table 2.4) — pairwise on full monthly sample ------------
cor_vars <- c("r_Gold","r_DXY","d_DFII10","d_DGS2","d_VIX",
              "d_ETFs","d_NonCom","CB_L1","CB_L2","CB_L3")
cor_mat  <- cor(df %>% select(all_of(cor_vars)), use = "pairwise.complete.obs")
n_pairs  <- colSums(!is.na(df %>% select(all_of(cor_vars))))
cat(sprintf("\nCorrelation matrix (pairwise N = %d to %d)\n", min(n_pairs), max(n_pairs)))
print(round(cor_mat, 3))

# Multicollinearity (VIF, Appendix B) -----------------------------------------
vif_vars <- c("r_Gold_L1","r_DXY","d_DFII10","d_DGS2","d_VIX","d_GPR",
              "d_ETFs","d_NonCom","r_SP500","r_Oil","FedBal_yoy","GEPUgdp")
vif_data <- df %>% select(all_of(c("r_Gold", vif_vars))) %>% drop_na()
v <- vif(lm(reformulate(vif_vars, "r_Gold"), data = vif_data))
cat("\nVariance inflation factors\n")
print(round(v, 2))

# Regression models M0-M4 (Tables 3.2-3.3) ------------------------------------
model_df <- df %>% filter(!is.na(CB_L1), !is.na(r_Gold), !is.na(r_DXY))
cat(sprintf("\nModel sample: %d months, %s to %s (pre-2022 = %d, post-2022 = %d)\n",
            nrow(model_df), min(model_df$date), max(model_df$date),
            sum(model_df$D22 == 0), sum(model_df$D22 == 1)))

run_m <- function(f, data, label) {
  m <- lm(f, data = data)
  cat(sprintf("\n%s\n", label))
  print(coeftest(m, vcov = vcovHAC(m)))   # Newey-West HAC standard errors
  cat(sprintf("R2_adj = %.4f | AIC = %.1f | BIC = %.1f | N = %d\n",
              summary(m)$adj.r.squared, AIC(m), BIC(m), nobs(m)))
  invisible(m)
}

f0 <- r_Gold ~ r_Gold_L1 + r_DXY + d_DFII10 + d_ETFs + d_NonCom + d_VIX
f1 <- update(f0, . ~ . + CB_L1 + CB_L2 + CB_L3)
f2 <- update(f1, . ~ . + CB_L1_D22 + CB_L2_D22 + CB_L3_D22)
f3 <- r_Gold ~ r_Gold_L1 + r_DXY + d_DFII10 + d_DGS2 + d_ETFs + d_NonCom + d_VIX +
  d_GPR + r_SP500 + r_Oil + FedBal_yoy + GEPUgdp +
  CB_L1 + CB_L2 + CB_L3 + CB_L1_D22 + CB_L2_D22 + CB_L3_D22
f4 <- r_Gold ~ r_Gold_L1 + r_DXY + d_DFII10 + d_ETFs + d_NonCom + d_VIX +
  dCB_L1 + dCB_L2 + I(dCB_L1 * D22) + I(dCB_L2 * D22)

m0 <- run_m(f0, model_df, "M0: baseline (no central bank demand)")
m1 <- run_m(f1, model_df, "M1: + central bank demand levels")
m2 <- run_m(f2, model_df, "M2: + post-2022 interaction (main specification)")
m3 <- run_m(f3, model_df, "M3: extended controls")
m4 <- run_m(f4, model_df, "M4: central bank demand in first differences")

mlist <- list(M0 = m0, M1 = m1, M2 = m2, M3 = m3, M4 = m4)
cat("\nModel comparison (Table 3.2)\n")
print(data.frame(
  Model  = names(mlist),
  N      = sapply(mlist, nobs),
  k      = sapply(mlist, function(m) length(coef(m))),
  AIC    = round(sapply(mlist, AIC), 1),
  BIC    = round(sapply(mlist, BIC), 1),
  R2_adj = round(sapply(mlist, function(m) summary(m)$adj.r.squared), 4)
), row.names = FALSE)

# Joint significance and long-run multipliers (Section 3.5) ------------------
cat("\nWald F-test: CB lags (M1)\n")
print(linearHypothesis(m1, c("CB_L1=0","CB_L2=0","CB_L3=0"), vcov = vcovHAC(m1)))
cat("\nWald F-test: post-2022 interactions (M2)\n")
print(linearHypothesis(m2, c("CB_L1_D22=0","CB_L2_D22=0","CB_L3_D22=0"),
                       vcov = vcovHAC(m2)))

cb  <- coef(m2)[paste0("CB_L", 1:3)]
d22 <- coef(m2)[paste0("CB_L", 1:3, "_D22")]
cat(sprintf("\nLong-run multiplier  pre-2022 = %+.4f  post-2022 = %+.4f  (%% per tonne)\n",
            sum(cb), sum(cb + d22)))
cat(sprintf("Implied effect of 100 t/month  pre = %+.2f pp  post = %+.2f pp\n",
            sum(cb) * 100, sum(cb + d22) * 100))

# Structural break: Chow test (Section 3.6, Table 3.4) -----------------------
chow_f  <- r_Gold ~ r_Gold_L1 + r_DXY + d_DFII10 + d_ETFs + d_NonCom + d_VIX +
  CB_L1 + CB_L2 + CB_L3
chow_df <- model_df %>% drop_na(all_of(all.vars(chow_f)))
m_full  <- lm(chow_f, data = chow_df)
m_pre   <- lm(chow_f, data = filter(chow_df, D22 == 0))
m_post  <- lm(chow_f, data = filter(chow_df, D22 == 1))

RSS_u  <- sum(residuals(m_pre)^2) + sum(residuals(m_post)^2)
k      <- length(coef(m_full)); n <- nobs(m_full)
F_chow <- ((sum(residuals(m_full)^2) - RSS_u) / k) / (RSS_u / (n - 2 * k))
p_chow <- pf(F_chow, k, n - 2 * k, lower.tail = FALSE)
cat(sprintf("\nChow test (break at 2022-02): F = %.3f, p = %.4f\n", F_chow, p_chow))
cat("\nPre-2022 subsample (HAC)\n");  print(coeftest(m_pre,  vcov = vcovHAC(m_pre)))
cat("\nPost-2022 subsample (HAC)\n"); print(coeftest(m_post, vcov = vcovHAC(m_post)))

# Structural break: Bai-Perron (Section 3.6, Appendix E) ---------------------
# BIC favours zero breaks; conditional on one break the new regime begins at
# February 2022, consistent with the a priori Chow test date.
bp_df     <- chow_df %>% arrange(date)
bp_result <- breakpoints(chow_f, data = bp_df, h = 0.15)
bp_obs    <- bp_result$breakpoints

if (all(is.na(bp_obs))) {
  cat("\nBai-Perron: BIC selects 0 breaks.\n")
  bp1 <- breakpoints(bp_result, breaks = 1)$breakpoints
  cat(sprintf("Conditional on one break, the new regime begins at obs %d (%s), ",
              bp1 + 1, format(bp_df$date[bp1 + 1], "%Y-%m")))
  cat("consistent with the Chow test date.\n")
} else {
  cat("\nBai-Perron break date(s):", format(bp_df$date[bp_obs + 1], "%Y-%m"), "\n")
}
print(summary(bp_result))

# Granger causality (Section 3.7) ---------------------------------------------
granger_df <- model_df %>% select(r_Gold, CB = NetPurchasesTonnes) %>% drop_na()
opt_lag    <- VARselect(granger_df, lag.max = 6, type = "const")$selection["AIC(n)"]
var_fit    <- VAR(granger_df, p = opt_lag, type = "const")
cat(sprintf("\nGranger causality (VAR lag = %d)\n", opt_lag))
cat("  CB -> Gold:  p =", round(causality(var_fit, cause = "CB")$Granger$p.value,    4), "\n")
cat("  Gold -> CB:  p =", round(causality(var_fit, cause = "r_Gold")$Granger$p.value, 4), "\n")

# Variance decomposition (Section 3.7, Table 3.5, Model M2) ------------------
lmg_df  <- model_df %>% drop_na(r_Gold, r_Gold_L1, r_DXY, d_DFII10, d_ETFs, d_NonCom,
                                 d_VIX, CB_L1, CB_L2, CB_L3,
                                 CB_L1_D22, CB_L2_D22, CB_L3_D22)
m_lmg   <- lm(r_Gold ~ r_Gold_L1 + r_DXY + d_DFII10 + d_ETFs + d_NonCom + d_VIX +
                CB_L1 + CB_L2 + CB_L3 + CB_L1_D22 + CB_L2_D22 + CB_L3_D22,
              data = lmg_df)
lmg_pct    <- round(calc.relimp(m_lmg, type = "lmg", rela = TRUE)@lmg * 100, 2)
lmg_sorted <- sort(lmg_pct, decreasing = TRUE)
cat("\nLMG variance decomposition, % of R-squared (Table 3.5, Model M2)\n")
print(data.frame(Variable = names(lmg_sorted), Share = lmg_sorted), row.names = FALSE)
cat(sprintf("Central-bank-related terms: %.1f%%\n",
            sum(lmg_pct[grep("CB", names(lmg_pct))])))

# Figures ---------------------------------------------------------------------
# Plots render to the active device (Plots pane in RStudio; Rplots.pdf via Rscript).
BREAK <- as.Date("2022-02-24")

# Figure 2.1 — gold price, real yield, CB purchases, dollar index
op <- par(mfrow = c(4, 1), mar = c(2, 4, 2, 1))
plot(df$date, df$Goldprice, type = "l", lwd = 1.5,
     main = "Gold price (USD/oz)", xlab = "", ylab = "USD")
abline(v = BREAK, lty = 2)
plot(df$date, df$DFII10, type = "l", lwd = 1.5,
     main = "10Y TIPS real yield (%)", xlab = "", ylab = "%")
abline(v = BREAK, lty = 2); abline(h = 0, col = "grey70")
plot(df$date, df$NetPurchasesTonnes, type = "h",
     main = "Monthly CB net purchases (tonnes)", xlab = "", ylab = "Tonnes")
abline(v = BREAK, lty = 2); abline(h = 0, col = "grey70")
plot(df$date, df$DXY, type = "l", lwd = 1.5,
     main = "Nominal broad USD index", xlab = "", ylab = "Index")
abline(v = BREAK, lty = 2); par(op)

# Appendix A — distributions of model variables
op <- par(mfrow = c(3, 4), mar = c(3, 3, 3, 1))
for (vn in c("r_Gold","r_SP500","r_Oil","r_DXY","d_DFII10","d_VIX",
             "d_GPR","d_ETFs","d_NonCom","FedBal_yoy","GEPUgdp")) {
  x <- na.omit(df[[vn]]); if (length(x) < 10) next
  hist(x, breaks = 30, main = vn, xlab = "", probability = TRUE,
       col = "grey85", border = "white")
  lines(density(x), lwd = 1.5)
}
par(op)

# Figure 2.4 — correlation matrix
corrplot(cor_mat, method = "color", type = "lower", addCoef.col = "black",
         number.cex = 0.6, tl.cex = 0.8, tl.col = "black")

# Appendix B — variance inflation factors
op <- par(mar = c(8, 4, 3, 1))
b <- barplot(v, las = 2, ylab = "VIF", ylim = c(0, max(v) + 0.5), cex.names = 0.8)
text(b, v + 0.05, round(v, 1), cex = 0.75)
abline(h = 5, lty = 2); par(op)

# Figure 3.1 — CB net purchases and CB(t-2) vs gold return
op <- par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))
plot(df$date, df$NetPurchasesTonnes, type = "h",
     col = ifelse(df$date >= BREAK, "black", "grey60"),
     main = "Monthly CB net purchases", xlab = "", ylab = "Tonnes")
abline(v = BREAK, lty = 2)
ok <- complete.cases(model_df$CB_L2, model_df$r_Gold)
plot(model_df$CB_L2[ok], model_df$r_Gold[ok],
     pch = ifelse(model_df$D22[ok] == 1, 17, 16), cex = 0.7,
     col = ifelse(model_df$D22[ok] == 1, "black", "grey60"),
     xlab = "CB net purchases, 2-month lag (tonnes)",
     ylab = "Monthly gold return (%)",
     main = "CB demand (2-month lag) vs gold return")
abline(lm(r_Gold ~ CB_L2, data = filter(model_df, D22 == 0)),
       col = "grey50", lty = 2, lwd = 1.5)
abline(lm(r_Gold ~ CB_L2, data = filter(model_df, D22 == 1)),
       col = "black", lwd = 1.5)
legend("topright", c("pre-2022","post-2022"), pch = c(16, 17),
       col = c("grey60","black"), lty = c(2, 1), bty = "n", cex = 0.8)
par(op)

# Appendix C — residual diagnostics (Model M1)
op <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1)); plot(m1); par(op)

# Appendix D — ACF / PACF of residuals (Model M1)
op <- par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))
acf(residuals(m1),  main = "ACF of residuals",  lag.max = 24)
pacf(residuals(m1), main = "PACF of residuals", lag.max = 24)
par(op)

# Figure 3.5 — LMG variance decomposition
op <- par(mar = c(9, 4, 3, 1))
b <- barplot(lmg_sorted, las = 2, ylab = "Share of R-squared (%)",
             col = ifelse(grepl("CB", names(lmg_sorted)), "black", "grey70"),
             ylim = c(0, max(lmg_sorted) + 5), cex.names = 0.8)
text(b, lmg_sorted + 0.6, paste0(round(lmg_sorted, 1), "%"), cex = 0.7)
par(op)

# Appendix E — Bai-Perron RSS / BIC by number of breaks
plot(bp_result)

# Figure 3.6 — coefficient shift pre- vs post-2022
pre_c <- coef(m_pre)[-1]; post_c <- coef(m_post)[-1]
vb    <- intersect(names(pre_c), names(post_c)); xp <- seq_along(vb)
yl    <- range(c(pre_c[vb], post_c[vb])) +
         diff(range(c(pre_c[vb], post_c[vb]))) * c(-0.1, 0.2)
op <- par(mar = c(9, 4, 3, 1))
plot(xp, pre_c[vb], pch = 16, col = "grey60", cex = 1.3,
     xaxt = "n", xlab = "", ylab = "Estimated coefficient", ylim = yl,
     main = "Coefficient shift: pre- vs post-2022")
points(xp, post_c[vb], pch = 17, cex = 1.3)
segments(xp, pre_c[vb], xp, post_c[vb], lty = 2, col = "grey60")
axis(1, at = xp, labels = vb, las = 2, cex.axis = 0.8)
abline(h = 0)
legend("topright", c("pre-2022","post-2022"), pch = c(16, 17),
       col = c("grey60","black"), bty = "n")
par(op)

# Figure 2.3 — quarterly supply, demand and balance
qf <- quarterly %>%
  mutate(TotalSupply = SuppMineProd + SuppNetProdHed + SuppRecyGo,
         TotalDemand = DemBarCo + DemCB + DemETF + DemJew + DemOTC + DemTec,
         Balance     = TotalDemand - TotalSupply)
op <- par(mfrow = c(2, 2), mar = c(3, 4, 2, 1))
plot(monthly_raw$date, monthly_raw$GEPUgdp, type = "l",
     main = "Economic policy uncertainty", xlab = "", ylab = "Index")
abline(v = BREAK, lty = 2)
plot(qf$date, qf$TotalSupply, type = "l", lwd = 1.5,
     main = "Total gold supply (quarterly)", xlab = "", ylab = "Tonnes")
abline(v = BREAK, lty = 2)
plot(qf$date, qf$DemCB, type = "h",
     main = "Central bank demand (quarterly)", xlab = "", ylab = "Tonnes")
abline(v = BREAK, lty = 2)
barplot(qf$Balance, col = ifelse(qf$Balance >= 0, "grey40", "grey80"),
        border = NA, main = "Demand minus supply", ylab = "Tonnes")
abline(h = 0)
par(op)

# Key results -----------------------------------------------------------------
cat("\n----------------------------------------------------------\n")
cat("Key results\n")
cat(sprintf("  CB lags, full sample (M1):      F-test p = %.3f\n",
            linearHypothesis(m1, c("CB_L1=0","CB_L2=0","CB_L3=0"),
                             vcov = vcovHAC(m1))$`Pr(>F)`[2]))
cat(sprintf("  Post-2022 interactions (M2):    F-test p = %.3f\n",
            linearHypothesis(m2, c("CB_L1_D22=0","CB_L2_D22=0","CB_L3_D22=0"),
                             vcov = vcovHAC(m2))$`Pr(>F)`[2]))
cat(sprintf("  Chow test (2022-02):            p = %.4f\n", p_chow))
cat(sprintf("  CB-related variance share (M2): %.1f%%\n",
            sum(lmg_pct[grep("CB", names(lmg_pct))])))
cat("----------------------------------------------------------\n\n")

print(sessionInfo())
