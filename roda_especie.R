## Roda o fluxo completo de UMA especie: 05 (grupos), 06 (tamanho), 07 (combina)
## O 07 e rodado com o modelo principal e com o outro (sensibilidade).
## Uso no terminal:  Rscript roda_especie.R "Tayassu pecari"
## (normalmente chamado pelo 08_roda_todas.R)

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 1)
especie <- args[1]
source("00_config.R")   ## pastas e criterios

t0 <- Sys.time()
message("=== ", especie, " | inicio ", format(t0), " ===")

source("05_ajusta_grupos.R")
source("06_ajusta_tamanho.R")

## Modelo principal (escolhido no 05 pelo phi ajustado)
forcar_modelo <- NULL
source("07_combina_tendencias.R")
principal <- modelo_grupos

## Sensibilidade: o outro modelo de grupos
forcar_modelo <- setdiff(c("g_sem", "g_com"), principal)
source("07_combina_tendencias.R")

## Marca a especie como concluida
dir.create(dir_resultados, showWarnings = FALSE)
writeLines(c(especie, principal, format(Sys.time())),
           file.path(dir_resultados, paste0(gsub(" ", "_", especie), "_concluido.txt")))
message(sprintf("=== %s concluida em %.1f min ===", especie,
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
