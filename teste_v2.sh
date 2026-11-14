#!/bin/bash
## Teste rapido da versao 2 antes da rodada completa (~30 min).
## Roda Cebus unicolor (especie rapida) pelos dois caminhos:
##   1. forcando "sem inclinacao por UC"
##   2. criterio normal (Cebus convergiu com inclinacao na v1)
## Resultados em ajustes_teste_*/ e resultados_teste_*/ (ignorados pelo git).
## Uso: bash teste_v2.sh
set -u
mkdir -p logs
ok=0
for cfg in teste_config_sem_incl.R teste_config_com_incl.R; do
    log="logs/${cfg%.R}.log"
    echo "== $cfg ($(date +%H:%M)) -> $log"
    if CONFIG_ANALISE=$cfg Rscript roda_especie.R "Cebus unicolor" > "$log" 2>&1; then
        echo "   OK"
        grep -E "inclinacao por UC|modelo de grupos =|classificaveis" "$log"
    else
        echo "   ERRO: veja $log"; tail -20 "$log"; ok=1
    fi
done
for d in resultados_teste_sem_incl resultados_teste_com_incl; do
    echo "== arquivos em $d:"; ls "$d" 2>/dev/null
done
exit $ok
