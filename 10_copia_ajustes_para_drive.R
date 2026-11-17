## Copia para para_drive/ajustes_v2/ so os ajustes (.rds) necessarios para
## reproduzir a analise sem reajustar nada, para compartilhar (ex.: Google Drive).
##   - modelos finais de grupos (g_sem e g_com, com as iteracoes finais)
##   - g_sem com inclinacao por UC (teste usado pelo 05)
##   - modelos de tamanho usados no 07 (_tam_dist e _tam_dist_sc)
##   - escolha_modelos.csv e comparacao_modelos.csv
##   - a pasta de resultados (resultados_v2/), pequena, lida pelo 09
## Ficam de fora os _tam.rds (nao usados) e ajustes substituidos por versoes
## com mais iteracoes.
##
## Uso (na pasta do projeto):  Rscript 10_copia_ajustes_para_drive.R
## Quem recebe: copiar as pastas ajustes_v2/ e resultados_v2/ para dentro da
## pasta do projeto (ver INSTRUCOES_WINDOWS.md).

library(readr)
source(Sys.getenv("CONFIG_ANALISE", "00_config.R"))

esc <- read_csv(file.path(dir_ajustes, "escolha_modelos.csv"), show_col_types = FALSE)
arqs <- file.path(dir_ajustes, c("escolha_modelos.csv", "comparacao_modelos.csv"))
for (i in seq_len(nrow(esc))) {
    sp <- esc$especie[i]
    est <- if (esc$inclinacao_uc[i]) "incl" else "noincl"
    pref <- file.path(dir_ajustes, gsub(" ", "_", sp))
    arqs <- c(arqs,
              arq_modelo(sp, "g_sem", est, iter_final(sp, "g_sem")),
              arq_modelo(sp, "g_com", est, iter_final(sp, "g_com")),
              arq_modelo(sp, "g_sem", "incl"),
              paste0(pref, c("_tam_dist.rds", "_tam_dist_sc.rds")))
}
arqs <- unique(arqs)
faltam <- arqs[!file.exists(arqs)]
if (length(faltam)) stop("Arquivos nao encontrados:\n", paste(faltam, collapse = "\n"))

destino <- file.path("para_drive", dir_ajustes)
dir.create(destino, recursive = TRUE, showWarnings = FALSE)
ok <- file.copy(arqs, destino, overwrite = TRUE)
stopifnot(all(ok))
cat(sprintf("%d arquivos copiados para %s (%.1f GB)\n",
            length(arqs), destino, sum(file.size(arqs)) / 1e9))

## Resultados combinados (07/08), lidos pelo 09; nao estao no git
res <- list.files(dir_resultados, full.names = TRUE)
dir.create(file.path("para_drive", dir_resultados), showWarnings = FALSE)
stopifnot(all(file.copy(res, file.path("para_drive", dir_resultados), overwrite = TRUE)))
cat(sprintf("%d arquivos copiados para %s\n", length(res), file.path("para_drive", dir_resultados)))
