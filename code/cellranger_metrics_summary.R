library(readr)
metrics_summary_SRR13081815 <- read_csv("/data/PRJNA679294/cellranger/SRR13081815/outs/metrics_summary.csv") %>% mutate(Run = "SRR13081815")
metrics_summary_SRR13081816 <- read_csv("/data/PRJNA679294/cellranger/SRR13081816/outs/metrics_summary.csv") %>% mutate(Run = "SRR13081816")
metrics_summary <-
  bind_rows(
    metrics_summary_SRR13081815,
    metrics_summary_SRR13081816)

metrics_summary |>
  select("Estimated Number of Cells", "Run")

write_tsv(metrics_summary, here("metrics_summary.tsv"))
