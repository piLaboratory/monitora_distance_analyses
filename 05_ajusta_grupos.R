## Modelo de distancia para a densidade de GRUPOS, dados agregados por UC-ano.
## Ajusta dois modelos e compara:
##   m_sem : abundancia ~ ano_c + (ano_c || cnuc)
##   m_com : abundancia ~ ano_c + (ano_c || cnuc) + (1 | cnuc_ano)
## Deteccao (ambos): half-normal, log(sigma) ~ ano_c + (1 | cnuc)
##
## Uso:
##   especie <- "Tayassu pecari"
##   source("05_ajusta_grupos.R")
## Salva em ajustes_grupos/<Especie>_m_sem.rds e _m_com.rds

library(ubms)
source("04_prepara_dados_uc_ano.R")

## Configuracoes do amostrador (iguais as do modelo anterior)
cfg <- list(chains = 6, iter = 2000, warmup = 1000,
            cores = min(6, parallel::detectCores()),
            adapt_delta = 0.9, max_treedepth = 12)

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
f_sem <- ~ ano_c + (1 | cnuc) ~ ano_c + (ano_c || cnuc)
f_com <- ~ ano_c + (1 | cnuc) ~ ano_c + (ano_c || cnuc) + (1 | cnuc_ano)

## Ajusta (ou carrega, se ja existir)
for (m in c("m_sem", "m_com")) {
    if (file.exists(arq(m))) {
        assign(m, readRDS(arq(m)))
    } else {
        t0 <- Sys.time()
        fit <- ajusta(if (m == "m_sem") f_sem else f_com)
        message(sprintf("%s %s: %.1f min", especie, m,
                        as.numeric(difftime(Sys.time(), t0, units = "mins"))))
        saveRDS(fit, arq(m))
        assign(m, fit)
    }
}

################################################################################
## Comparacao
################################################################################
## 1. Convergencia: Rhat < 1.01 e n_eff > 400 nos parametros principais;
##    divergencias pos-aquecimento
convergencia <- function(fit) {
    s <- summary(fit@stanfit)$summary
    principais <- grepl("^beta_|^sigma_", rownames(s))
    np <- rstan::get_sampler_params(fit@stanfit, inc_warmup = FALSE)
    tibble(max_Rhat   = max(s[principais, "Rhat"], na.rm = TRUE),
           min_n_eff  = min(s[principais, "n_eff"], na.rm = TRUE),
           divergencias = sum(sapply(np, function(x) sum(x[, "divergent__"]))),
           treedepth_max = sum(sapply(np, function(x) sum(x[, "treedepth__"] >= cfg$max_treedepth))))
}

## 2. Tendencia media de grupos (%/ano), com IC 95%
R_anos <- ano_max - ano_min
tendencia <- function(fit) {
    ## coeficientes fixos de abundancia; a coluna de ano_c e localizada pelo nome
    b <- extract(fit, "beta_state")[[1]]
    b <- b[, match("ano_c", colnames(model.matrix(fit["state"])))]
    p <- 100 * (exp(b / R_anos) - 1)
    tibble(pct_ano = mean(p), ic_inf = quantile(p, 0.025),
           ic_sup = quantile(p, 0.975), P_declinio = mean(b < 0))
}

## 3. Dispersao de Pearson nos totais por UC-ano (esperado pela media
##    posterior preditiva). phi ~ 1 indica que a variancia esta bem descrita.
##    No m_com o esperado ja inclui o efeito de cada UC-ano, entao phi tende a
##    ficar <= 1 por construcao; o phi informativo e o do m_sem.
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
        bind_cols(modelo = "sem cnuc_ano", convergencia(m_sem), tendencia(m_sem),
                  phi = dispersao(m_sem)),
        bind_cols(modelo = "com cnuc_ano", convergencia(m_com), tendencia(m_com),
                  phi = dispersao(m_com))) |>
    mutate(especie = especie, .before = 1)
print(as.data.frame(comparacao), digits = 3)

## 4. Ajuste (qui-quadrado posterior preditivo do ubms)
gof_sem <- gof(m_sem, draws = 500); gof_com <- gof(m_com, draws = 500)
print(gof_sem); print(gof_com)

## 5. LOO. Atencao: com efeito por observacao (cnuc_ano) o PSIS-LOO costuma
##    ter muitos k de Pareto altos; nesse caso a comparacao nao e confiavel e
##    vale mais olhar convergencia, phi e o IC da tendencia.
ms <- modSel(fitList(sem_cnuc_ano = m_sem, com_cnuc_ano = m_com))
print(ms)

## Desvio do efeito UC-ano (quanto maior, mais variacao extra-Poisson entre anos)
print(summary(m_com, "state"))
