## Modelo de distancia para a densidade de GRUPOS, dados agregados por UC-ano.
## Ajusta dois modelos e compara:
##   g_sem : abundancia ~ ano_c + (ano_c || cnuc)
##   g_com : abundancia ~ ano_c + (ano_c || cnuc) + (1 | cnuc_ano)
## Deteccao (ambos): half-normal, log(sigma) ~ ano_c
##   Sem efeito de UC na deteccao: com poucas deteccoes por UC, deteccao e
##   abundancia de cada UC nao sao separaveis (divergencias concentradas em
##   sigma_det alto no teste com Tayassu pecari).
##
## Uso:
##   especie <- "Tayassu pecari"
##   source("05_ajusta_grupos.R")
## Salva em ajustes_grupos/<Especie>_g_sem.rds e _g_com.rds

library(ubms)
source("04_prepara_dados_uc_ano.R")

## Configuracoes do amostrador
cfg <- list(chains = 6, iter = 4000, warmup = 2000,
            cores = min(6, parallel::detectCores()),
            adapt_delta = 0.99, max_treedepth = 12)

dir.create("ajustes_grupos", showWarnings = FALSE)
arq <- function(m) file.path("ajustes_grupos",
                             paste0(gsub(" ", "_", especie), "_", m, ".rds"))

ajusta <- function(formula) {
    stan_distsamp(formula, data = umf_ua,
                  keyfun = "halfnorm", output = "density", unitsOut = "ha",
                  chains = cfg$chains, iter = cfg$iter, warmup = cfg$warmup,
                  cores = cfg$cores, refresh = 0,
                  control = list(adapt_delta = cfg$adapt_delta,
                                 max_treedepth = cfg$max_treedepth))
}

## Formula dupla do ubms: ~ deteccao ~ abundancia
formulas <- list(
    g_sem = ~ ano_c ~ ano_c + (ano_c || cnuc),
    g_com = ~ ano_c ~ ano_c + (ano_c || cnuc) + (1 | cnuc_ano))

## Ajusta (ou carrega, se ja existir)
for (m in names(formulas)) {
    if (file.exists(arq(m))) {
        assign(m, readRDS(arq(m)))
    } else {
        t0 <- Sys.time()
        fit <- ajusta(formulas[[m]])
        message(sprintf("%s %s: %.1f min", especie, m,
                        as.numeric(difftime(Sys.time(), t0, units = "mins"))))
        saveRDS(fit, arq(m))
        assign(m, fit)
    }
}

################################################################################
## Comparacao
################################################################################
## 1. Convergencia. Criterios: Rhat < 1.01 e n_eff > 400 nos parametros
##    principais (beta_ e sigma_), nenhuma divergencia, BFMI > 0.3.
##    O n_eff minimo dos desvios (sigma_) e mostrado a parte: e o ponto
##    fraco esperado com dados esparsos (funil dos efeitos aleatorios).
convergencia <- function(fit) {
    s <- summary(fit@stanfit)$summary
    principais <- grepl("^beta_|^sigma_", rownames(s))
    sigmas     <- grepl("^sigma_", rownames(s))
    tibble(max_Rhat      = max(s[principais, "Rhat"], na.rm = TRUE),
           min_n_eff     = min(s[principais, "n_eff"], na.rm = TRUE),
           min_n_eff_sigma = if (any(sigmas)) min(s[sigmas, "n_eff"]) else NA_real_,
           divergencias  = rstan::get_num_divergent(fit@stanfit),
           treedepth_max = rstan::get_num_max_treedepth(fit@stanfit),
           min_BFMI      = min(rstan::get_bfmi(fit@stanfit)))
}

## 2. Tendencia media de grupos (%/ano), com IC 95%
R_anos <- ano_max - ano_min
tendencia <- function(fit) {
    ## coeficientes fixos de abundancia; a coluna de ano_c e localizada pelo nome
    b <- extract(fit, "beta_state")[[1]]
    b <- b[, match("ano_c", ubms:::beta_names(fit["state"]))]
    p <- 100 * (exp(b / R_anos) - 1)
    tibble(pct_ano = mean(p), ic_inf = quantile(p, 0.025),
           ic_sup = quantile(p, 0.975), P_declinio = mean(b < 0))
}

## 3. Dispersao de Pearson nos totais por UC-ano (esperado pela media
##    posterior preditiva). phi ~ 1 indica que a variancia esta bem descrita.
##    No g_com o esperado ja inclui o efeito de cada UC-ano, entao phi tende a
##    ficar <= 1 por construcao; o phi informativo e o do g_sem.
##    posterior_predict devolve colunas na ordem UC-ano 1 (faixas 1..J),
##    UC-ano 2 (faixas 1..J), ... (conferido no codigo do ubms).
dispersao <- function(fit) {
    yrep <- posterior_predict(fit, "y", draws = 1000)  ## draws x (UC-ano*faixa)
    e <- matrix(colMeans(yrep), nrow = nrow(y_ua), byrow = TRUE)
    obs <- rowSums(y_ua); esp <- rowSums(e)
    sum((obs - esp)^2 / esp) / (length(obs) - 1)
}

comparacao <-
    bind_rows(
        bind_cols(modelo = "sem cnuc_ano", convergencia(g_sem), tendencia(g_sem),
                  phi = dispersao(g_sem)),
        bind_cols(modelo = "com cnuc_ano", convergencia(g_com), tendencia(g_com),
                  phi = dispersao(g_com))) |>
    mutate(especie = especie, .before = 1)
print(as.data.frame(comparacao), digits = 3)

## 4. Ajuste (qui-quadrado posterior preditivo do ubms)
gof_sem <- gof(g_sem, draws = 500); gof_com <- gof(g_com, draws = 500)
print(gof_sem); print(gof_com)

## 5. LOO. Atencao: com efeito por observacao (cnuc_ano) o PSIS-LOO costuma
##    ter muitos k de Pareto altos; nesse caso a comparacao nao e confiavel e
##    vale mais olhar convergencia, phi e o IC da tendencia.
ms <- modSel(fitList(sem_cnuc_ano = g_sem, com_cnuc_ano = g_com))
print(ms)

## Desvio do efeito UC-ano (quanto maior, mais variacao extra-Poisson entre anos)
print(summary(g_com, "state"))
