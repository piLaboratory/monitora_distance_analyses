## Combina densidade de grupos (ubms) e tamanho medio de grupo (brms) em
## densidade e tendencia de INDIVIDUOS, propagando a incerteza pelos sorteios.
##
##   D_ind(t) = D_grp(t) * (1 + exp(eta_tam(t)))
##   tendencia (%/ano) = 100 * ((D_ind(fim) / D_ind(ini))^(1/R) - 1)
##
## O tamanho medio (1 + exp(eta)) nao e log-linear no ano, por isso a
## tendencia e calculada pela razao entre o ultimo e o primeiro ano, e nao
## somando inclinacoes.
## Os dois modelos sao ajustados separadamente; os sorteios sao pareados ao
## acaso (posteriores independentes).
##
## Uso:
##   especie <- "Tayassu pecari"
##   forcar_modelo <- "g_com"     ## opcional (sensibilidade); se ausente ou NULL,
##                                ## usa a escolha gravada pelo 05
##   crit_min_grupos <- 15; crit_min_anos <- 3   ## opcional (padrao)
##   source("07_combina_tendencias.R")
## Produz: tend_media (UC tipica), dens_ind e tend_uc (por UC), gravados em
## resultados/<Especie>_<modelo>_*.csv

library(ubms); library(brms); library(tidyverse)
## Criterio para uma UC ter tendencia propria classificada: numero minimo de
## grupos detectados e de anos com ao menos uma deteccao
if (!exists("crit_min_grupos")) crit_min_grupos <- 15
if (!exists("crit_min_anos"))   crit_min_anos   <- 3
## Modelo de grupos: o forcado pelo usuario ou o escolhido no 05 (phi)
if (exists("forcar_modelo") && !is.null(forcar_modelo)) {
    modelo_grupos <- forcar_modelo
} else {
    escolha <- read_csv("ajustes_grupos/escolha_modelos.csv", show_col_types = FALSE)
    modelo_grupos <- escolha$modelo_principal[escolha$especie == especie]
    stopifnot(length(modelo_grupos) == 1)
}
message(sprintf("07 (%s): modelo de grupos = %s", especie, modelo_grupos))
source("04_prepara_dados_uc_ano.R")

pref <- file.path("ajustes_grupos", gsub(" ", "_", especie))
m_grp <- readRDS(paste0(pref, "_", modelo_grupos, ".rds"))
## Modelo de tamanho com distancia (06): o tamanho medio e avaliado a
## distancia zero, corrigindo o vies de deteccao de grupos grandes.
m_tam <- readRDS(paste0(pref, "_tam_dist.rds"))
dist_c0 <- (0 - mean(dados_tam$dist.m)) / sd(dados_tam$dist.m)  ## mesma padronizacao do 04
R_anos <- ano_max - ano_min

################################################################################
## Sorteios do modelo de grupos
################################################################################
## Efeitos aleatorios completos (fixo + aleatorio) por nivel, em sorteios.
## ranef() do ubms devolve so a media quando summary = FALSE, entao os
## sorteios sao montados aqui seguindo a mesma logica interna do ubms.
re_sorteios <- function(fit, submodel = "state") {
    sm   <- fit[submodel]
    re   <- ubms:::get_reTrms(sm@formula, sm@data)
    b    <- extract(fit, paste0("b_", submodel))[[1]]
    beta <- extract(fit, paste0("beta_", submodel))[[1]]
    bn   <- ubms:::beta_names(sm)
    out <- list()
    for (i in seq_along(re$cnms)) {
        fac <- names(re$cnms)[i]; trm <- re$cnms[[i]]
        d <- b[, (re$Gp[i] + 1):re$Gp[i + 1], drop = FALSE]
        if (trm %in% bn) d <- d + beta[, match(trm, bn)]
        colnames(d) <- levels(re$flist[[fac]])
        out[[paste(fac, trm, sep = ":")]] <- d
    }
    out
}

beta_g <- extract(m_grp, "beta_state")[[1]]
bn_g   <- ubms:::beta_names(m_grp["state"])
g0 <- beta_g[, match("(Intercept)", bn_g)]   ## log densidade de grupos (por ha) em ano_min
g1 <- beta_g[, match("ano_c", bn_g)]         ## variacao no log ao longo do periodo
re_g  <- re_sorteios(m_grp, "state")
g0_uc <- re_g[["cnuc:(Intercept)"]]          ## sorteios x UCs
g1_uc <- re_g[["cnuc:ano_c"]]
ucs   <- colnames(g0_uc)

################################################################################
## Sorteios do modelo de tamanho
################################################################################
## Preditor linear em ano_min e ano_max, a distancia zero, para a UC tipica
## (re_formula = NA) e para cada UC do modelo de grupos. UCs sem dados de
## tamanho recebem um efeito aleatorio sorteado da distribuicao estimada.
## brms:: explicito porque o ubms tem uma funcao com o mesmo nome.
nd <- expand_grid(cnuc = ucs, ano_c = c(0, 1)) |> mutate(dist_c = dist_c0)
eta_uc  <- brms::posterior_linpred(m_tam, newdata = nd, allow_new_levels = TRUE,
                                   sample_new_levels = "gaussian")
eta_med <- brms::posterior_linpred(m_tam, newdata = tibble(ano_c = c(0, 1), dist_c = dist_c0),
                                   re_formula = NA)

## Pareia os sorteios ao acaso, com o mesmo numero nos dois modelos
set.seed(1)
nd_s <- min(nrow(beta_g), nrow(eta_uc))
ig <- sample(nrow(beta_g), nd_s); it <- sample(nrow(eta_uc), nd_s)

resume <- function(p) c(media = mean(p), ic_inf = unname(quantile(p, 0.025)),
                        ic_sup = unname(quantile(p, 0.975)), P_declinio = mean(p < 0))

################################################################################
## UC tipica (efeitos aleatorios = 0)
################################################################################
tam_ini <- 1 + exp(eta_med[it, 1]); tam_fim <- 1 + exp(eta_med[it, 2])
pct_grp <- 100 * (exp(g1[ig] / R_anos) - 1)
pct_tam <- 100 * ((tam_fim / tam_ini)^(1 / R_anos) - 1)
pct_ind <- 100 * ((exp(g1[ig]) * tam_fim / tam_ini)^(1 / R_anos) - 1)

tend_media <- bind_rows(
    c(componente = "grupos",     resume(pct_grp)),
    c(componente = "tamanho",    resume(pct_tam)),
    c(componente = "individuos", resume(pct_ind))) |>
    mutate(across(-componente, as.numeric), especie = especie, .before = 1)
print(as.data.frame(tend_media), digits = 3)

## Densidade de individuos (ind/ha) no primeiro e no ultimo ano, UC tipica
dens_ind <- tibble(
    especie = especie,
    ano = c(ano_min, ano_max),
    media  = c(mean(exp(g0[ig]) * tam_ini), mean(exp(g0[ig] + g1[ig]) * tam_fim)),
    ic_inf = c(quantile(exp(g0[ig]) * tam_ini, 0.025),
               quantile(exp(g0[ig] + g1[ig]) * tam_fim, 0.025)),
    ic_sup = c(quantile(exp(g0[ig]) * tam_ini, 0.975),
               quantile(exp(g0[ig] + g1[ig]) * tam_fim, 0.975)))
print(dens_ind)

################################################################################
## Por UC
################################################################################
col_ini <- match(paste(ucs, 0), paste(nd$cnuc, nd$ano_c))
col_fim <- match(paste(ucs, 1), paste(nd$cnuc, nd$ano_c))
tend_uc <- map_dfr(seq_along(ucs), function(k) {
    r_tam <- (1 + exp(eta_uc[it, col_fim[k]])) / (1 + exp(eta_uc[it, col_ini[k]]))
    p <- 100 * ((exp(g1_uc[ig, k]) * r_tam)^(1 / R_anos) - 1)
    as_tibble_row(c(cnuc = ucs[k], resume(p)))
}) |>
    mutate(across(-cnuc, as.numeric),
           evidencia = case_when(ic_sup < 0 ~ "declinio",
                                 ic_inf > 0 ~ "aumento",
                                 TRUE       ~ "inconclusivo"))

## Dados de cada UC: grupos detectados e anos com deteccao (apos truncagem)
dados_uc <-
    ua |>
    mutate(n_grp = rowSums(y_ua), cnuc = as.character(cnuc)) |>
    summarise(n_grupos = sum(n_grp), anos_amostrados = n(),
              anos_com_deteccao = sum(n_grp > 0), .by = cnuc)

## UCs que nao atingem o criterio: a tendencia estimada reflete sobretudo a
## media geral (encolhimento), e nao e classificada
tend_uc <-
    tend_uc |>
    left_join(dados_uc, by = "cnuc") |>
    mutate(classificavel = n_grupos >= crit_min_grupos &
                           anos_com_deteccao >= crit_min_anos,
           evidencia = if_else(classificavel, evidencia, "nao classificavel"),
           especie = especie, modelo = modelo_grupos) |>
    left_join(distinct(monitora.mam, cnuc, uc), by = "cnuc") |>
    relocate(especie, modelo, cnuc, uc) |>
    arrange(desc(classificavel), media)
print(as.data.frame(tend_uc), digits = 3)
message(sprintf("07 (%s): %d de %d UCs classificaveis (>= %d grupos em >= %d anos)",
                especie, sum(tend_uc$classificavel), nrow(tend_uc),
                crit_min_grupos, crit_min_anos))

################################################################################
## Grava os resultados
################################################################################
dir.create("resultados", showWarnings = FALSE)
pref_res <- file.path("resultados", paste0(gsub(" ", "_", especie), "_", modelo_grupos))
write_csv(mutate(tend_media, modelo = modelo_grupos), paste0(pref_res, "_tend_media.csv"))
write_csv(mutate(dens_ind,  modelo = modelo_grupos), paste0(pref_res, "_dens_ind.csv"))
write_csv(tend_uc, paste0(pref_res, "_tend_uc.csv"))
