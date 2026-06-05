# Central Bank Gold Demand and the Global Gold Price

Replication code for the paper:  
**"The Impact of Central Bank Demand on the Global Gold Price in the Context of Transforming International Reserves"**  
O. Steshenko, NaUKMA, 2026

## How to run

1. Clone or download this repository
2. Place the dataset in the same folder as the script
3. Open `gold_analysis.R` in RStudio and run, or from terminal:

```r
Rscript gold_analysis.R
```

## Repository structure

Gold_price_central_banks/
├── Data_term_paper.xlsx
├── gold_analysis.R
└── README.md

## Data sources

| Variable | Source |
|---|---|
| Gold price (end-of-month) | LBMA [27] |
| Dollar index, yields, VIX, oil, Fed balance | FRED [28–36] |
| Gold ETF holdings | WGC Goldhub [30] |
| CB net purchases | IMF IFS via DBnomics [37] |
| CFTC non-commercial net long | CFTC [31] |
| Geopolitical risk index | Caldara & Iacoviello (2022) [7] |
| Global EPU | Davis (2016) [26] |

## Requirements

R ≥ 4.6.0 Packages installed automatically on first run:  
`readxl`, `dplyr`, `tidyr`, `lubridate`, `psych`, `corrplot`,  
`car`, `lmtest`, `sandwich`, `strucchange`, `tseries`, `vars`, `relaimpo`
