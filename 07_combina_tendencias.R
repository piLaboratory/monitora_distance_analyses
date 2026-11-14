## Tendencias de GRUPOS (resultado principal) e de INDIVIDUOS (secundario),
## propagando a incerteza pelos sorteios dos modelos.
##
##   D_grp(t) = exp(g0 + g1 * t)                        (modelo de grupos, 05)
##   D_ind(t) = D_grp(t) * (1 + exp(eta_tam(t)))        (tamanho a distancia 0, 06)
##   tendencia (%/ano) = 100 * ((D(fim) / D(ini))^(1/R) - 1)
##
## Individuos sao calculados com dois modelos de tamanho:
##   "individuos"    : tamanho com contagens parciais censuradas
##   "individuos_sc" : tamanho so com contagens totais (sem censura)
## O tamanho medio (1 + exp(eta)) nao e log-linear no ano; por isso a tendencia
## e calculada pela razao entre o ultimo e o primeiro ano.
## Os modelos sao ajustados separadamente; os sorteios sao pareados ao acaso.
##
## Uso:
##   especie <- "Tayassu pecari"
##   forcar_modelo <- "g_com"     ## opcional (sensibilidade); se ausente ou NULL,
##                                ## usa a escolha gravada pelo 05
##   source("07_combina_tendencias.R")
## Produz tend_media, dens_ind e tend_uc, gravados em
## <dir_resultados>/<Especie>_<modelo>_*.csv

library(ubms); library(brms); library(tidyverse)
source("04_prepara_dados_uc_ano.R")   ## dados e configuracoes (00_config.R)

## Modelo de grupos: o forcado pelo usuario ou o escolhido no 05 (phi)
escolha <- read_csv(file.path(dir_ajustes, "escolha_modelos.csv"), show_col_types = FALSE) |>
    filter(especie == .env$especie)
stopifnot(nrow(escolha) == 1)
modelo_grupos <- if (exists("forcar_modelo") && !is.null(forcar_modelo)) forcar_modelo else
    escolha$modelo_principal
message(sprintf("07 (%s): modelo de grupos = %s", especie, modelo_grupos))

pref  <- file.path(dir_ajustes, gsub(" ", "_", especie))
estrutura <- if (escolha$inclinacao_uc) "incl" else "noincl"
m_grp <- readRDS(paste0(pref, "_", modelo_grupos, "_", estrutura, ".rds"))
m_tam <- list(individuos    = readRDS(paste0(pref, "_tam_dist.rds")),
              individuos_sc = readRDS(paste0(pref, "_tam_dist_sc.rds")))
dist_c0 <- (0 - mean(dados_tam$dist.m)) / sd(dados_tam$dist.m)  ## distancia 0, padronizacao do 04
R_anos  <- ano_max - ano_min

################################################################################
## Sorteios do modelo de grupos
################################################################################
## Efeitos aleatorios completos (fixo + aleatorio) por nivel, em sorteios.
## ranef() do ubms devolve so a media quando summary = FALSE, entao os
## sorteios sao montados aqui seguindo a mesma logica interna do ubms.
re_sorteios <- function(fit, submodel = "state") {
    sm   <- fit[submodel]
    re   <- ubms:::get_reTrms(sm@formula, sm@data)
    b    <- ubms::extract(fit, paste0("b_", submodel))[[1]]
    beta <- ubms::extract(fit, paste0("beta_", submodel))[[1]]
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

beta_g <- ubms::extract(m_grp, "beta_state")[[1]]
bn_g   <- ubms:::beta_names(m_grp["state"])
g0 <- beta_g[, match("(Intercept)", bn_g)]   ## log densidade de grupos (por ha) em ano_min
g1 <- beta_g[, match("ano_c", bn_g)]         ## variacao no log ao longo do periodo
re_g  <- re_sorteios(m_grp, "state")
ucs   <- colnames(re_g[["cnuc:(Intercept)"]])
## Inclinacao por UC: so existe nas especies ajustadas com (ano_c || cnuc)
inclinacao_uc <- "cnuc:ano_c" %in% names(re_g)
g1_uc <- if (inclinacao_uc) re_g[["cnuc:ano_c"]] else
    matrix(g1, nrow = length(g1), ncol = length(ucs), dimnames = list(NULL, ucs))

################################################################################
## Sorteios dos modelos de tamanho (a distancia zero)
################################################################################
## UC tipica (re_formula = NA) e cada UC do modelo de grupos. UCs sem dados de
## tamanho recebem um efeito sorteado da distribuicao estimada.
nd <- expand_grid(cnuc = ucs, ano_c = c(0, 1)) |> mutate(dist_c = dist_c0)
col_ini <- match(paste(ucs, 0), paste(nd$cnuc, nd$ano_c))
col_fim <- match(paste(ucs, 1), paste(nd$cnuc, nd$ano_c))
eta <- map(m_tam, function(m) list(
    uc  = brms::posterior_linpred(m, newdata = nd, allow_new_levels = TRUE,
                                  sample_new_levels = "gaussian"),
    med = brms::posterior_linpred(m, newdata = tibble(ano_c = c(0, 1), dist_c = dist_c0),
                                  re_formula = NA)))

## Pareia os sorteios ao acaso, com o mesmo numero em todos os modelos
set.seed(1)
nd_s <- min(nrow(beta_g), map_int(eta, ~ nrow(.x$med)))
ig <- sample(nrow(beta_g), nd_s)
it <- map(eta, ~ sample(nrow(.x$med), nd_s))

resume <- function(p) c(media = mean(p), ic_inf = unname(quantile(p, 0.025)),
                        ic_sup = unname(quantile(p, 0.975)), P_declinio = mean(p < 0))
pct <- function(razao) 100 * (razao^(1 / R_anos) - 1)

################################################################################
## UC tipica (efeitos aleatorios = 0)
################################################################################
tam <- imap(eta, function(e, k) list(ini = 1 + exp(e$med[it[[k]], 1]),
                                     fim = 1 + exp(e$med[it[[k]], 2])))
comp <- list(
    grupos        = pct(exp(g1[ig])),
    tamanho       = pct(tam$individuos$fim / tam$individuos$ini),
    individuos    = pct(exp(g1[ig]) * tam$individuos$fim / tam$individuos$ini),
    tamanho_sc    = pct(tam$individuos_sc$fim / tam$individuos_sc$ini),
    individuos_sc = pct(exp(g1[ig]) * tam$individuos_sc$fim / tam$individuos_sc$ini))
tend_media <- imap_dfr(comp, ~ c(componente = .y, resume(.x))) |>
    mutate(across(-componente, as.numeric),
           especie = especie, modelo = modelo_grupos,
           inclinacao_uc = inclinacao_uc, .before = 1)
print(as.data.frame(tend_media), digits = 3)

## Densidades (por ha) no primeiro e no ultimo ano, UC tipica
q <- function(x) c(media = mean(x), ic_inf = unname(quantile(x, 0.025)),
                   ic_sup = unname(quantile(x, 0.975)))
dens <- list(
    grupos        = list(exp(g0[ig]), exp(g0[ig] + g1[ig])),
    individuos    = list(exp(g0[ig]) * tam$individuos$ini,
                         exp(g0[ig] + g1[ig]) * tam$individuos$fim),
    individuos_sc = list(exp(g0[ig]) * tam$individuos_sc$ini,
                         exp(g0[ig] + g1[ig]) * tam$individuos_sc$fim))
dens_ind <- imap_dfr(dens, function(d, k) bind_rows(
    c(componente = k, ano = ano_min, q(d[[1]])),
    c(componente = k, ano = ano_max, q(d[[2]])))) |>
    mutate(across(-componente, as.numeric),
           especie = especie, modelo = modelo_grupos, .before = 1)
print(as.data.frame(dens_ind), digits = 3)

################################################################################
## Por UC (formato longo: uma linha por UC e componente)
################################################################################
tend_uc <- map_dfr(seq_along(ucs), function(k) {
    r_tam <- map2(eta, it, ~ (1 + exp(.x$uc[.y, col_fim[k]])) / (1 + exp(.x$uc[.y, col_ini[k]])))
    r_grp <- exp(g1_uc[ig, k])
    bind_rows(
        c(cnuc = ucs[k], componente = "grupos",        resume(pct(r_grp))),
        c(cnuc = ucs[k], componente = "individuos",    resume(pct(r_grp * r_tam$individuos))),
        c(cnuc = ucs[k], componente = "individuos_sc", resume(pct(r_grp * r_tam$individuos_sc))))
}) |>
    mutate(across(c(media, ic_inf, ic_sup, P_declinio), as.numeric),
           evidencia = case_when(ic_sup < 0 ~ "declinio",
                                 ic_inf > 0 ~ "aumento",
                                 TRUE       ~ "inconclusivo"))

## Dados de cada UC: grupos detectados e anos com deteccao (apos truncagem)
dados_uc <-
    ua |>
    mutate(n_grp = rowSums(y_ua), cnuc = as.character(cnuc)) |>
    summarise(n_grupos = sum(n_grp), anos_amostrados = n(),
              anos_com_deteccao = sum(n_grp > 0), .by = cnuc)

## Classificacao por UC so quando ha inclinacao por UC no modelo e a UC atinge
## o criterio de dados; caso contrario a estimativa reflete a media geral.
tend_uc <-
    tend_uc |>
    left_join(dados_uc, by = "cnuc") |>
    mutate(classificavel = inclinacao_uc & n_grupos >= crit_min_grupos &
                           anos_com_deteccao >= crit_min_anos,
           evidencia = case_when(!inclinacao_uc ~ "sem inclinacao por UC",
                                 !classificavel ~ "nao classificavel",
                                 TRUE           ~ evidencia),
           especie = especie, modelo = modelo_grupos) |>
    left_join(distinct(monitora.mam, cnuc, uc), by = "cnuc") |>
    relocate(especie, modelo, cnuc, uc, componente) |>
    arrange(componente, desc(classificavel), media)
print(as.data.frame(filter(tend_uc, componente == "grupos")), digits = 3)
message(sprintf("07 (%s): %d de %d UCs classificaveis (>= %d grupos em >= %d anos)%s",
                especie, sum(tend_uc$classificavel[tend_uc$componente == "grupos"]),
                length(ucs), crit_min_grupos, crit_min_anos,
                if (inclinacao_uc) "" else "; modelo sem inclinacao por UC"))

################################################################################
## Grava os resultados
################################################################################
dir.create(dir_resultados, showWarnings = FALSE)
pref_res <- file.path(dir_resultados, paste0(gsub(" ", "_", especie), "_", modelo_grupos))
write_csv(tend_media, paste0(pref_res, "_tend_media.csv"))
write_csv(dens_ind,   paste0(pref_res, "_dens_ind.csv"))
write_csv(tend_uc,    paste0(pref_res, "_tend_uc.csv"))
