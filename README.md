# Pipeline-final-somaticas
## Esse pipeline tem por objetivo realizar a analise de variantes somaticas de 30 vcfs, primeiramente de forma bruta e pós aplicação de filtros de qualidade, no final realizando o comparativo entre às 2 coortes.
### 1) Clonar o repositório.
```bash
git clone https://github.com/Kaneka4850/Pipeline-final-somaticas.git
```
## Após a clonagem do repositório, é necessario a instalação do programa bcftools, para manipulação dos arquivos vcf.
```bash
sudo apt install bcftools
```
## Por fim, realizar a extração dos vcfs, utilizando o comando unzip
```bash
unzip -o /content/Pipeline-final-somaticas/liftOver-hg38-MF-annotVep.zip
```
## 2) Realizar a conversão dos arquivos vcf em tsv.
### Devido a dificuldade de realizar a manipulação de arquivos .vcf, iremos realizar a conversão dos arquivos para .tsv, garantindo uma padronização para utilizar o pandas
```bash
# Define o diretório onde os arquivos de entrada estão
DIR_INPUT="/content/liftOver-hg38-MF-annotVep"

# Define o diretório onde os arquivos de saída serão salvos
DIR_OUTPUT="outputs"

# Cria o diretório de saída (o -p garante que não dê erro se já existir)
mkdir -p "$DIR_OUTPUT"

# Inicia o loop
for ARQUIVO in ${DIR_INPUT}/liftOver_WP*_hg19ToHg38.vep.vcf; do

    # 1. Extrair o ID da amostra
    NOME_BASE=$(basename "$ARQUIVO")
    SAMPLE_ID=$(echo "$NOME_BASE" | cut -d'_' -f2)

    # Define o caminho completo do arquivo de saída (agora dentro da pasta outputs)
    OUTPUT_FILE="${DIR_OUTPUT}/liftOver_${SAMPLE_ID}_hg19ToHg38.vep.tsv"

    echo "Processando amostra: ${SAMPLE_ID} -> salvando em ${OUTPUT_FILE}"

    # 2. Criar o cabeçalho
    bcftools +split-vep -l "$ARQUIVO" | \
    cut -f2 | \
    tr '\n\r' '\t' | \
    awk '{print("CHROM\tPOS\tREF\tALT\t"$0"FILTER\tTumorID\tGT\tDP\tAD\tAF\tNormalID\tNGT\tNDP\tNAD\tNAF")}' \
    > "$OUTPUT_FILE"

    # 3. Adicionar as variantes
    bcftools +split-vep \
    -f '%CHROM\t%POS\t%REF\t%ALT\t%CSQ\t%FILTER\t[%SAMPLE\t%GT\t%DP\t%AD\t%AF\t]\n' \
    -i 'FMT/DP>=20 && FMT/AF>=0.05' -d -A tab "$ARQUIVO" \
    -p x >> "$OUTPUT_FILE"

done

echo "Processamento finalizado. Arquivos salvos em: ${DIR_OUTPUT}/"
```
## 3) Utilização da biblioteca pandas para gerar os arquivos de risco e de risco alto, lembrando que esse script tem um painel de genes fixos, caso for necessario alterar. Arrumar no código

```python
import pandas as pd
import glob
import os

# =========================
# CONFIGURAÇÃO (ACEITE)
# =========================
INPUT_PATTERN = "outputs/*.tsv"   # TSVs derivados do VCF anotado com VEP
OUTPUT_VARIANTS = "variants_high_risk.tsv"
OUTPUT_SAMPLES = "sample_risk.tsv"

# Painel fixo exigido
PANEL_GENES = {
    "TP53", "EZH2", "CBL", "U2AF1", "SRSF2",
    "IDH1", "IDH2", "NRAS", "KRAS"
}

# Consequences permitidas (exatas)
ALLOWED_CONSEQUENCES = [
    "missense_variant",
    "stop_gained",
    "frameshift_variant",
    "start_lost",
    "splice_"
]

# =========================
# FUNÇÕES AUXILIARES
# =========================
def smart_read_tsv(path):
    """Lê TSVs potencialmente malformados sem quebrar o pipeline."""
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        header = f.readline().rstrip("\n").split("\t")
        n_cols = len(header)
        rows = [line.rstrip("\n").split("\t")[:n_cols] for line in f]
    return pd.DataFrame(rows, columns=header)

def parse_float(val):
    try:
        return float(str(val).replace(",", "."))
    except:
        return 0.0

def consequence_ok(consequence):
    if pd.isna(consequence):
        return False
    return any(term in consequence for term in ALLOWED_CONSEQUENCES)

# =========================
# PROCESSAMENTO
# =========================
all_variants = []
summary = []

files = sorted(glob.glob(INPUT_PATTERN))
## Pega o nome da amostra, ex: WPXXX
for file in files:
    sample_id = os.path.basename(file).split("_")[1]

    df = smart_read_tsv(file)

    # Coluna de gene
    gene_col = "SYMBOL" if "SYMBOL" in df.columns else "Gene"
    df["GENE"] = df[gene_col].astype(str).str.upper().str.strip()

    # Conversões numéricas
    df["DP"] = df["DP"].apply(parse_float) if "DP" in df.columns else 0
    df["VAF"] = df["VAF"].apply(parse_float) if "VAF" in df.columns else 0

    # =========================
    # REGRAS DE FILTRO (EXATAS)
    # =========================
    mask_gene = df["GENE"].isin(PANEL_GENES)
    mask_filter = df["FILTER"] == "PASS"
    mask_impact = df["IMPACT"].isin(["MODERATE", "HIGH"])
    mask_consequence = df["Consequence"].apply(consequence_ok)
    mask_quality = (df["DP"] >= 20) | (df["VAF"] >= 0.05)

    high_risk = df[
        mask_gene &
        mask_filter &
        mask_quality &
        (mask_impact | mask_consequence)
    ].copy()

    high_risk["SAMPLEID"] = sample_id

    # =========================
    # VARIANTS TABLE
    # =========================
    if not high_risk.empty:
        all_variants.append(
            high_risk[
                [
                    "SAMPLEID",
                    "CHROM", "POS", "REF", "ALT",
                    "GENE",
                    "Consequence",
                    "IMPACT",
                    "FILTER",
                    "DP",
                    "VAF"
                ]
            ]
        )

    # =========================
    # SAMPLE SUMMARY
    # =========================
    summary.append({
        "SAMPLEID": sample_id,
        "MAIOR_RISCO": "SIM" if not high_risk.empty else "NÃO",
        "TP53_PRESENTE": "SIM" if any(high_risk["GENE"] == "TP53") else "NÃO",
        "GENES_ALTO_RISCO_ENCONTRADOS": ",".join(sorted(high_risk["GENE"].unique())),
        "N_VARIANTES_ALTO_RISCO": int(len(high_risk))
    })

# =========================
# SAÍDA FINAL
# =========================
df_variants = pd.concat(all_variants, ignore_index=True) if all_variants else pd.DataFrame()
df_samples = pd.DataFrame(summary)

df_variants.to_csv(OUTPUT_VARIANTS, sep="\t", index=False)
df_samples.to_csv(OUTPUT_SAMPLES, sep="\t", index=False)

# =========================
# RESUMO (<= 10 LINHAS)
# =========================
print("PROCESSAMENTO FINALIZADO")
print(f"Nº total de amostras processadas: {len(df_samples)}")
print(f"Nº de amostras com MAIOR_RISCO = SIM: {(df_samples['MAIOR_RISCO'] == 'SIM').sum()}")
print(f"Nº de amostras com TP53_PRESENTE = SIM: {(df_samples['TP53_PRESENTE'] == 'SIM').sum()}")
print(f"Arquivo gerado: {OUTPUT_VARIANTS}")
print(f"Arquivo gerado: {OUTPUT_SAMPLES}")
```

### Agora, para termos uma analise dos vcfs sem o pré processamento, iremos utilizar um script em python para gerar um dashboard interativo, permitindo assim a verificação das variantes somaticas identificadas. Lembrando que o painel foi submetido ao cgi previamente, devido ao número expressivo de variantes, a analise pode levar em média de 15-30 minutos. Por esse motivo, a analise foi feita previamente.

```python
import pandas as pd
import json
from google.colab import files

# 1. Carregamento
try:
    df_alt = pd.read_csv('/content/Pipeline-final-somaticas/alterations.tsv', sep='\t').fillna('')
    df_bio = pd.read_csv('/content/Pipeline-final-somaticas/biomarkers.tsv', sep='\t').fillna('')
    print("✅ Arquivos carregados!")
except:
    print("❌ Erro nos arquivos. Verifique os caminhos.")
    raise

# 2. Conversão Segura
json_alt = json.dumps(df_alt.to_dict(orient='records'))
json_bio = json.dumps(df_bio.to_dict(orient='records'))

# 3. Template HTML com Contagem Dinâmica de Biomarcadores
html_template = f"""
<!DOCTYPE html>
<html lang="pt-BR" class="light">
<head>
    <meta charset="UTF-8">
    <title>Dashboard Bioinformática | TCC</title>
    <script src="https://cdn.tailwindcss.com"></script>
    <script src="https://cdn.jsdelivr.net/npm/apexcharts"></script>
    <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.13.6/css/jquery.dataTables.min.css">
    <style>
        :root {{
            --bg: #f1f5f9; --card: #ffffff;
            --text-main: #0f172a; --border: #cbd5e1;
        }}
        .dark {{
            --bg: #020617; --card: #0f172a;
            --text-main: #f8fafc; --border: #1e293b;
        }}
        body {{ background: var(--bg); color: var(--text-main); transition: 0.3s; font-family: system-ui; }}
        .card {{ background: var(--card); border: 1px solid var(--border); border-radius: 1rem; padding: 1.5rem; }}

        /* Estilo DataTables */
        .dataTables_wrapper {{ color: var(--text-main) !important; }}
        table.dataTable thead th {{ background: var(--bg) !important; color: var(--text-main) !important; border-bottom: 2px solid var(--border) !important; font-weight: 800 !important; }}
        table.dataTable tbody td {{ border-bottom: 1px solid var(--border) !important; color: var(--text-main) !important; }}

        /* Badges de Oncogenicidade */
        .badge {{ padding: 4px 12px; border-radius: 99px; font-weight: bold; font-size: 0.72rem; }}
        .onco {{ background: #fee2e2; color: #991b1b; border: 1px solid #f87171; }}
        .pass {{ background: #f1f5f9; color: #475569; border: 1px solid #cbd5e1; }}
    </style>
</head>
<body class="p-4 md:p-10">

<div class="max-w-7xl mx-auto">
    <header class="flex flex-col md:flex-row justify-between items-center mb-10 gap-6">
        <div>
            <h1 class="text-4xl font-black italic text-blue-600">DASHBOARD FINAL</h1>
            <p class="text-slate-500 font-bold text-xs uppercase tracking-widest">Análise de Amostras e Biomarcadores</p>
        </div>
        <div class="flex flex-wrap gap-4 items-end">
            <div class="flex flex-col">
                <label class="text-[10px] font-bold uppercase text-slate-500 mb-1">Amostra</label>
                <select id="sampleSelect" class="card p-2 text-sm outline-none focus:ring-2 focus:ring-blue-500 bg-white dark:bg-slate-900"></select>
            </div>
            <div class="flex flex-col">
                <label class="text-[10px] font-bold uppercase text-slate-500 mb-1">Filtrar Gene</label>
                <input type="text" id="geneSearch" placeholder="Ex: TP53" class="card p-2 text-sm outline-none focus:ring-2 focus:ring-blue-500 bg-white dark:bg-slate-900">
            </div>
            <button onclick="toggleTheme()" class="card p-2.5 shadow-sm hover:bg-slate-50 transition">🌓</button>
        </div>
    </header>

    <div class="grid grid-cols-1 md:grid-cols-4 gap-6 mb-10 text-center">
        <div class="card border-l-4 border-l-blue-500 shadow-sm">
            <p class="text-[10px] font-bold text-slate-400 uppercase">Variantes</p>
            <h2 id="stat-total" class="text-3xl font-black">0</h2>
        </div>
        <div class="card border-l-4 border-l-red-500 shadow-sm">
            <p class="text-[10px] font-bold text-slate-400 uppercase">Oncogênicas</p>
            <h2 id="stat-onco" class="text-3xl font-black text-red-500">0</h2>
        </div>
        <div class="card border-l-4 border-l-emerald-500 shadow-sm">
            <p class="text-[10px] font-bold text-slate-400 uppercase">Biomarcadores</p>
            <h2 id="stat-bio" class="text-3xl font-black text-emerald-500">0</h2>
        </div>
        <div class="card flex items-center justify-center shadow-sm">
            <label class="flex items-center gap-2 cursor-pointer">
                <input type="checkbox" id="relevantOnly" class="w-4 h-4">
                <span class="text-[10px] font-bold text-slate-500 uppercase">Apenas Drivers</span>
            </label>
        </div>
    </div>

    <div class="grid grid-cols-1 lg:grid-cols-2 gap-8 mb-12">
        <div class="card shadow-md">
            <h3 class="font-bold mb-4 text-slate-400 uppercase text-xs">Distribuição</h3>
            <div id="oncoChart"></div>
        </div>
        <div class="card shadow-md">
            <h3 class="font-bold mb-4 text-slate-400 uppercase text-xs">Frequência Gênica (Top 10)</h3>
            <div id="geneChart"></div>
        </div>
    </div>

    <div class="space-y-10">
        <section class="card shadow-lg">
            <h2 class="text-xl font-black mb-6 text-blue-600">🔬 VARIANTES (ALTERATIONS)</h2>
            <table id="altTable" class="w-full text-sm">
                <thead><tr><th>Amostra</th><th>Gene</th><th>Proteína</th><th>Tipo</th><th>Sumário</th><th>Predição</th></tr></thead>
            </table>
        </section>

        <section class="card shadow-lg">
            <h2 class="text-xl font-black mb-6 text-emerald-600">💊 BIOMARCADORES E DROGAS</h2>
            <table id="bioTable" class="w-full text-sm">
                <thead><tr><th>Amostra</th><th>Alteração</th><th>Droga</th><th>Doença</th><th>Resposta</th><th>Evidência</th></tr></thead>
            </table>
        </section>
    </div>
</div>

<script src="https://code.jquery.com/jquery-3.7.0.min.js"></script>
<script src="https://cdn.datatables.net/1.13.6/js/jquery.dataTables.min.js"></script>

<script>
    const altData = {json_alt};
    const bioData = {json_bio};
    let altTable, bioTable, oncoChart, geneChart;

    function toggleTheme() {{
        const isDark = document.documentElement.classList.toggle('dark');
        const mode = isDark ? 'dark' : 'light';
        const color = isDark ? '#f8fafc' : '#0f172a';
        oncoChart.updateOptions({{ theme: {{ mode }}, chart: {{ foreColor: color }} }});
        geneChart.updateOptions({{ theme: {{ mode }}, chart: {{ foreColor: color }} }});
    }}

    function init() {{
        // Setup Amostras
        const samples = [...new Set(altData.map(d => d['SAMPLE']))].sort();
        const select = document.getElementById('sampleSelect');
        select.innerHTML = '<option value="">Todas</option>' + samples.map(s => `<option value="${{s}}">${{s}}</option>`).join('');

        // Tabela Alt
        altTable = $('#altTable').DataTable({{
            data: altData,
            columns: [
                {{ data: 'SAMPLE' }},
                {{ data: 'CGI-Gene', render: d => `<b class="text-blue-600 font-bold">${{d}}</b>` }},
                {{ data: 'CGI-Protein Change' }},
                {{ data: 'CGI-Type' }},
                {{ data: 'CGI-Oncogenic Summary', render: d => {{
                    const isO = (d||'').toLowerCase().includes('oncogenic');
                    return `<span class="badge ${{isO ? 'onco' : 'pass'}}">${{d}}</span>`;
                }} }},
                {{ data: 'CGI-Oncogenic Prediction' }}
            ]
        }});

        // Tabela Bio
        bioTable = $('#bioTable').DataTable({{
            data: bioData,
            columns: [
                {{ data: 'Sample ID' }}, {{ data: 'Alterations' }}, {{ data: 'Drugs' }}, {{ data: 'Diseases' }},
                {{ data: 'Response' }},
                {{ data: 'Evidence', render: d => `<span class="bg-blue-600 text-white w-6 h-6 flex items-center justify-center rounded-full text-[10px] font-bold mx-auto">${{d}}</span>` }}
            ]
        }});

        // Inicia Gráficos
        const isDark = document.documentElement.classList.contains('dark');
        const textColor = isDark ? '#f8fafc' : '#0f172a';

        oncoChart = new ApexCharts(document.querySelector("#oncoChart"), {{
            chart: {{ type: 'donut', height: 280, foreColor: textColor }},
            labels: ['Oncogênico', 'Outros'],
            series: [0, 0],
            colors: ['#ef4444', '#94a3b8']
        }});
        oncoChart.render();

        geneChart = new ApexCharts(document.querySelector("#geneChart"), {{
            chart: {{ type: 'bar', height: 280, toolbar: {{show:false}}, foreColor: textColor }},
            series: [{{ name: 'Mutas', data: [] }}],
            xaxis: {{ categories: [] }},
            colors: ['#3b82f6'],
            plotOptions: {{ bar: {{ borderRadius: 6, horizontal: true }} }}
        }});
        geneChart.render();

        // Listeners
        $('#sampleSelect, #relevantOnly, #geneSearch').on('change keyup', () => {{
            altTable.draw(); bioTable.draw(); updateStats();
        }});

        // Filtro DataTables
        $.fn.dataTable.ext.search.push((settings, data, idx, row) => {{
            const sample = $('#sampleSelect').val();
            const relevant = $('#relevantOnly').is(':checked');
            const geneSearch = $('#geneSearch').val().toUpperCase();

            const rSample = data[0];
            const rGene = data[1].toUpperCase();
            const rSum = data[4].toLowerCase();

            if (sample && rSample !== sample) return false;
            if (geneSearch && !rGene.includes(geneSearch)) return false;
            if (relevant && settings.nTable.id === 'altTable' && !rSum.includes('oncogenic')) return false;

            return true;
        }});

        updateStats();
    }}

    function updateStats() {{
        const sample = $('#sampleSelect').val();
        const gene = $('#geneSearch').val().toUpperCase();

        // FILTRO DINÂMICO DE VARIANTES
        let fAlt = altData;
        if (sample) fAlt = fAlt.filter(d => d['SAMPLE'] === sample);
        if (gene) fAlt = fAlt.filter(d => d['CGI-Gene'].toUpperCase().includes(gene));

        // FILTRO DINÂMICO DE BIOMARCADORES (Correção aqui!)
        let fBio = bioData;
        if (sample) fBio = fBio.filter(d => d['Sample ID'] === sample);
        if (gene) fBio = fBio.filter(d => d['Alterations'].toUpperCase().includes(gene));

        const oncoCount = fAlt.filter(d => (d['CGI-Oncogenic Summary']||'').toLowerCase().includes('oncogenic')).length;

        // Atualiza Cards
        document.getElementById('stat-total').innerText = fAlt.length;
        document.getElementById('stat-onco').innerText = oncoCount;
        document.getElementById('stat-bio').innerText = fBio.length;

        // Atualiza Gráficos
        oncoChart.updateSeries([oncoCount, fAlt.length - oncoCount]);
        const gFreq = {{}};
        fAlt.forEach(d => gFreq[d['CGI-Gene']] = (gFreq[d['CGI-Gene']] || 0) + 1);
        const top = Object.entries(gFreq).sort((a,b) => b[1]-a[1]).slice(0, 10);
        geneChart.updateOptions({{ xaxis: {{ categories: top.map(x => x[0]) }} }});
        geneChart.updateSeries([{{ data: top.map(x => x[1]) }}]);
    }}

    $(document).ready(init);
</script>
</body>
</html>
"""

# 4. Injeção e Download
final_html = html_template.replace('{json_alt}', json_alt).replace('{json_bio}', json_bio)
with open('dashboard_tcc_final_sem_filtro.html', 'w', encoding='utf-8') as f:
    f.write(final_html)
files.download('dashboard_tcc_final_sem_filtro.html')
```

# Injeção de arquivos no CGI. 
## Após realizar a filtragem dos vcfs em arquivos .tsv, iremos fazer a verificação das variantes relevantes no contexto clínico de mielofibrose. Para injeção no GCI, lembrando que o CGI é uma API, sendo necessario realizar as requisições necessarias
### Para obter a chave da API, entre por esse link, realize seu cadastro e gere seu token, lembrando que o Token é pessoal:
### https://www.cancergenomeinterpreter.org/rest_api

```bash
OUTPUT="df_final-cgi.txt"
INPUT="/content/variants_high_risk.tsv"

# 1. CRIAR CABEÇALHO PADRÃO CGI
echo -e "CHR\tPOS\tREF\tALT\tSAMPLE" > "$OUTPUT"

# 2. EXTRAIR E REORDENAR AS COLUNAS DO ARQUIVO JÁ FILTRADO
# Colunas no high_risk: 1:SAMPLEID, 2:CHROM, 3:POS, 4:REF, 5:ALT
# Queremos para o CGI: CHR(2), POS(3), REF(4), ALT(5), SAMPLE(1)
tail -n +2 "$INPUT" | awk -F'\t' -v OFS="\t" '{print $2, $3, $4, $5, $1}' >> "$OUTPUT"

# 3. VERIFICAÇÃO
echo "✅ Arquivo gerado com sucesso: $OUTPUT"
echo "--------------------------------------"
echo "Primeiras 5 linhas do arquivo final:"
head -n 5 "$OUTPUT"
echo "--------------------------------------"
echo "Total de variantes para o CGI: $(tail -n +2 $OUTPUT | wc -l)"
```

## 2 - Geração do Job_id do CGI

```python
import requests
# cabeçalho
headers = {'Authorization': '# Insira seu email aqui # insira seu token aqui'}
payload = {'cancer_type': 'HEMATO', 'title': 'Somatic filtered', 'reference': 'hg38'}
# requisição
r = requests.post('https://www.cancergenomeinterpreter.org/api/v1',
                headers=headers,
                files={
                        'mutations': open('/content/df_final-cgi.txt', 'rb')
                        },
                data=payload)
r.json() # Formato json para verificação
```
3- Verificação das logs do job do CGI
```python
import requests

job_id = input("Digite seu job_id (Obtido na 3º célula)")
headers = {'Authorization': '# Insira seu email aqui # insira seu token aqui'}
r = requests.get('https://www.cancergenomeinterpreter.org/api/v1/%s' % job_id, headers=headers)
r.json()
```

## 4 - Enriquecimento dos dados via CGI
```python 
import requests
job_id = input("Digite seu job_id (Obtido na 3º célula)")

headers = {'Authorization': '# Insira seu email aqui # insira seu token aqui'}
payload={'action':'logs'}
r = requests.get('https://www.cancergenomeinterpreter.org/api/v1/%s' % job_id, headers=headers, params=payload)
r.json()
```
## 5 - Download do arquivo zipados. (Esperar uns 5 minutos entre o passo 4 e o passo 5)
```python
import requests
job_id = input("Digite seu job_id (Obtido na 3º célula)")

headers = {'Authorization': '# Insira seu email aqui # insira seu token aqui'} # permissões do CGI
payload={'action':'download'} # passando o que é pra ele fazer
r = requests.get('https://www.cancergenomeinterpreter.org/api/v1/%s' % job_id, headers=headers, params=payload) # requisições
with open('sample01.zip', 'wb') as fd:
    fd.write(r._content)
```
## 6 - Unzip do request gerado
```bash
unzip -o sample01.zip
```
## 7 - Verificação do arquivo alterations.tsv
```python
import pandas as pd
pd.read_csv('/content/alterations.tsv',sep='\t',index_col=False, engine= 'python')
```
## 8 - Verificação do arquivo biomarkers.tsv
```python
import pandas as pd
pd.read_csv('/content/biomarkers.tsv',sep='\t',index_col=False, engine= 'python')
```
## 9 - Alterar o nome do arquivo alterations.tsv, para não misturar com o gerado sem filtro
```python
mv alterations.tsv alterations_filtered.tsv #troca de nome do arquivo final, para não confundir
```

## 10 - Comparativo entre a amostra sem filtro vs amostra filtrada, gerando um gráfico de interceção de tabelas
```python
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# 1. DEFINIÇÃO DOS CAMINHOS (Conforme você passou)
path_geral = '/content/Pipeline-final-somaticas/alterations.tsv'
path_filtrado = '/content/alterations_filtered.tsv'

# Verificação de segurança
if not os.path.exists(path_geral) or not os.path.exists(path_filtrado):
    print("⚠️ Atenção: Verifique se os nomes dos arquivos e pastas estão corretos no seu Drive/Colab.")
else:
    # 2. CARREGAMENTO DOS DADOS
    # O CGI usa tabulação (\t) como separador
    df_geral = pd.read_csv(path_geral, sep='\t')
    df_filtrado = pd.read_csv(path_filtrado, sep='\t')

    # 3. CRIAÇÃO DE UMA CHAVE ÚNICA PARA COMPARAÇÃO
    # Usamos Amostra + CHR + POS + REF + ALT para garantir que o "match" seja perfeito
    def create_key(df, sample_col, chr_col, pos_col, ref_col, alt_col):
        return (df[sample_col].astype(str) + "_" +
                df[chr_col].astype(str).str.replace('chr', '') + "_" +
                df[pos_col].astype(str) + "_" +
                df[ref_col].astype(str) + "_" +
                df[alt_col].astype(str))

    # Adaptando nomes de colunas (o CGI costuma variar entre CHROMOSOME e CHR)
    chr_col_geral = 'CHR' if 'CHR' in df_geral.columns else 'CHROMOSOME'
    chr_col_filt = 'CHR' if 'CHR' in df_filtrado.columns else 'CHROMOSOME'

    df_geral['key'] = create_key(df_geral, 'SAMPLE', chr_col_geral, 'POS', 'REF', 'ALT')
    df_filtrado['key'] = create_key(df_filtrado, 'SAMPLE', chr_col_filt, 'POS', 'REF', 'ALT')

    # 4. ANÁLISE DE INTERSECÇÃO
    intersecao = df_geral[df_geral['key'].isin(df_filtrado['key'])]

    # 5. CONTAGEM POR AMOSTRA
    count_geral = df_geral.groupby('SAMPLE').size().reset_index(name='Total_Geral')
    count_filt = df_filtrado.groupby('SAMPLE').size().reset_index(name='Total_Filtrado')
    count_inter = intersecao.groupby('SAMPLE').size().reset_index(name='Interseccao')

    # Unindo os resultados em uma única tabela comparativa
    resumo = pd.merge(count_geral, count_filt, on='SAMPLE', how='outer')
    resumo = pd.merge(resumo, count_inter, on='SAMPLE', how='outer').fillna(0)

    # Converter para inteiros
    for col in ['Total_Geral', 'Total_Filtrado', 'Interseccao']:
        resumo[col] = resumo[col].astype(int)

    # 6. EXIBIÇÃO DOS RESULTADOS
    print("📊 RESUMO GLOBAL DA ANÁLISE")
    print(f"Total de variantes no arquivo GERAL: {len(df_geral)}")
    print(f"Total de variantes no arquivo FILTRADO: {len(df_filtrado)}")
    print(f"Total de variantes na INTERSECÇÃO: {len(intersecao)}")
    print("-" * 50)
    print("\n📋 TABELA COMPARATIVA POR AMOSTRA (Top 10):")
    print(resumo.sort_values(by='Total_Filtrado', ascending=False).head(10))

    # 7. VISUALIZAÇÃO (Para o seu Portfólio)
    plt.figure(figsize=(12, 6))
    resumo_melt = resumo.melt(id_vars='SAMPLE', value_vars=['Total_Filtrado', 'Interseccao'],
                              var_name='Tipo', value_name='Quantidade')

    sns.barplot(data=resumo_melt, x='SAMPLE', y='Quantidade', hue='Tipo')
    plt.xticks(rotation=45)
    plt.title('Comparação: Variantes Filtradas vs Intersecção com Geral')
    plt.tight_layout()
    plt.show()

    # Salvar o resultado para baixar
    resumo.to_csv('comparativo_detalhado_cgi.csv', index=False)
    print("\n✅ Arquivo 'comparativo_detalhado_cgi.csv' gerado com sucesso!")
```
Geração do dashboard com as variantes filtradas.

```python
import pandas as pd
import json
from google.colab import files

# 1. Carregamento
try:
    df_alt = pd.read_csv('/content/alterations_filtered.tsv', sep='\t').fillna('')
    df_bio = pd.read_csv('/content/biomarkers.tsv', sep='\t').fillna('')
    print("✅ Arquivos carregados!")
except:
    print("❌ Erro nos arquivos. Verifique os caminhos.")
    raise

# 2. Conversão Segura
json_alt = json.dumps(df_alt.to_dict(orient='records'))
json_bio = json.dumps(df_bio.to_dict(orient='records'))

# 3. Template HTML com Contagem Dinâmica de Biomarcadores
html_template = f"""
<!DOCTYPE html>
<html lang="pt-BR" class="light">
<head>
    <meta charset="UTF-8">
    <title>Dashboard Bioinformática | TCC</title>
    <script src="https://cdn.tailwindcss.com"></script>
    <script src="https://cdn.jsdelivr.net/npm/apexcharts"></script>
    <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.13.6/css/jquery.dataTables.min.css">
    <style>
        :root {{
            --bg: #f1f5f9; --card: #ffffff;
            --text-main: #0f172a; --border: #cbd5e1;
        }}
        .dark {{
            --bg: #020617; --card: #0f172a;
            --text-main: #f8fafc; --border: #1e293b;
        }}
        body {{ background: var(--bg); color: var(--text-main); transition: 0.3s; font-family: system-ui; }}
        .card {{ background: var(--card); border: 1px solid var(--border); border-radius: 1rem; padding: 1.5rem; }}

        /* Estilo DataTables */
        .dataTables_wrapper {{ color: var(--text-main) !important; }}
        table.dataTable thead th {{ background: var(--bg) !important; color: var(--text-main) !important; border-bottom: 2px solid var(--border) !important; font-weight: 800 !important; }}
        table.dataTable tbody td {{ border-bottom: 1px solid var(--border) !important; color: var(--text-main) !important; }}

        /* Badges de Oncogenicidade */
        .badge {{ padding: 4px 12px; border-radius: 99px; font-weight: bold; font-size: 0.72rem; }}
        .onco {{ background: #fee2e2; color: #991b1b; border: 1px solid #f87171; }}
        .pass {{ background: #f1f5f9; color: #475569; border: 1px solid #cbd5e1; }}
    </style>
</head>
<body class="p-4 md:p-10">

<div class="max-w-7xl mx-auto">
    <header class="flex flex-col md:flex-row justify-between items-center mb-10 gap-6">
        <div>
            <h1 class="text-4xl font-black italic text-blue-600">DASHBOARD FINAL</h1>
            <p class="text-slate-500 font-bold text-xs uppercase tracking-widest">Análise de Amostras e Biomarcadores</p>
        </div>
        <div class="flex flex-wrap gap-4 items-end">
            <div class="flex flex-col">
                <label class="text-[10px] font-bold uppercase text-slate-500 mb-1">Amostra</label>
                <select id="sampleSelect" class="card p-2 text-sm outline-none focus:ring-2 focus:ring-blue-500 bg-white dark:bg-slate-900"></select>
            </div>
            <div class="flex flex-col">
                <label class="text-[10px] font-bold uppercase text-slate-500 mb-1">Filtrar Gene</label>
                <input type="text" id="geneSearch" placeholder="Ex: TP53" class="card p-2 text-sm outline-none focus:ring-2 focus:ring-blue-500 bg-white dark:bg-slate-900">
            </div>
            <button onclick="toggleTheme()" class="card p-2.5 shadow-sm hover:bg-slate-50 transition">🌓</button>
        </div>
    </header>

    <div class="grid grid-cols-1 md:grid-cols-4 gap-6 mb-10 text-center">
        <div class="card border-l-4 border-l-blue-500 shadow-sm">
            <p class="text-[10px] font-bold text-slate-400 uppercase">Variantes</p>
            <h2 id="stat-total" class="text-3xl font-black">0</h2>
        </div>
        <div class="card border-l-4 border-l-red-500 shadow-sm">
            <p class="text-[10px] font-bold text-slate-400 uppercase">Oncogênicas</p>
            <h2 id="stat-onco" class="text-3xl font-black text-red-500">0</h2>
        </div>
        <div class="card border-l-4 border-l-emerald-500 shadow-sm">
            <p class="text-[10px] font-bold text-slate-400 uppercase">Biomarcadores</p>
            <h2 id="stat-bio" class="text-3xl font-black text-emerald-500">0</h2>
        </div>
        <div class="card flex items-center justify-center shadow-sm">
            <label class="flex items-center gap-2 cursor-pointer">
                <input type="checkbox" id="relevantOnly" class="w-4 h-4">
                <span class="text-[10px] font-bold text-slate-500 uppercase">Apenas Drivers</span>
            </label>
        </div>
    </div>

    <div class="grid grid-cols-1 lg:grid-cols-2 gap-8 mb-12">
        <div class="card shadow-md">
            <h3 class="font-bold mb-4 text-slate-400 uppercase text-xs">Distribuição</h3>
            <div id="oncoChart"></div>
        </div>
        <div class="card shadow-md">
            <h3 class="font-bold mb-4 text-slate-400 uppercase text-xs">Frequência Gênica (Top 10)</h3>
            <div id="geneChart"></div>
        </div>
    </div>

    <div class="space-y-10">
        <section class="card shadow-lg">
            <h2 class="text-xl font-black mb-6 text-blue-600">🔬 VARIANTES (ALTERATIONS)</h2>
            <table id="altTable" class="w-full text-sm">
                <thead><tr><th>Amostra</th><th>Gene</th><th>Proteína</th><th>Tipo</th><th>Sumário</th><th>Predição</th></tr></thead>
            </table>
        </section>

        <section class="card shadow-lg">
            <h2 class="text-xl font-black mb-6 text-emerald-600">💊 BIOMARCADORES E DROGAS</h2>
            <table id="bioTable" class="w-full text-sm">
                <thead><tr><th>Amostra</th><th>Alteração</th><th>Droga</th><th>Doença</th><th>Resposta</th><th>Evidência</th></tr></thead>
            </table>
        </section>
    </div>
</div>

<script src="https://code.jquery.com/jquery-3.7.0.min.js"></script>
<script src="https://cdn.datatables.net/1.13.6/js/jquery.dataTables.min.js"></script>

<script>
    const altData = {json_alt};
    const bioData = {json_bio};
    let altTable, bioTable, oncoChart, geneChart;

    function toggleTheme() {{
        const isDark = document.documentElement.classList.toggle('dark');
        const mode = isDark ? 'dark' : 'light';
        const color = isDark ? '#f8fafc' : '#0f172a';
        oncoChart.updateOptions({{ theme: {{ mode }}, chart: {{ foreColor: color }} }});
        geneChart.updateOptions({{ theme: {{ mode }}, chart: {{ foreColor: color }} }});
    }}

    function init() {{
        // Setup Amostras
        const samples = [...new Set(altData.map(d => d['SAMPLE']))].sort();
        const select = document.getElementById('sampleSelect');
        select.innerHTML = '<option value="">Todas</option>' + samples.map(s => `<option value="${{s}}">${{s}}</option>`).join('');

        // Tabela Alt
        altTable = $('#altTable').DataTable({{
            data: altData,
            columns: [
                {{ data: 'SAMPLE' }},
                {{ data: 'CGI-Gene', render: d => `<b class="text-blue-600 font-bold">${{d}}</b>` }},
                {{ data: 'CGI-Protein Change' }},
                {{ data: 'CGI-Type' }},
                {{ data: 'CGI-Oncogenic Summary', render: d => {{
                    const isO = (d||'').toLowerCase().includes('oncogenic');
                    return `<span class="badge ${{isO ? 'onco' : 'pass'}}">${{d}}</span>`;
                }} }},
                {{ data: 'CGI-Oncogenic Prediction' }}
            ]
        }});

        // Tabela Bio
        bioTable = $('#bioTable').DataTable({{
            data: bioData,
            columns: [
                {{ data: 'Sample ID' }}, {{ data: 'Alterations' }}, {{ data: 'Drugs' }}, {{ data: 'Diseases' }},
                {{ data: 'Response' }},
                {{ data: 'Evidence', render: d => `<span class="bg-blue-600 text-white w-6 h-6 flex items-center justify-center rounded-full text-[10px] font-bold mx-auto">${{d}}</span>` }}
            ]
        }});

        // Inicia Gráficos
        const isDark = document.documentElement.classList.contains('dark');
        const textColor = isDark ? '#f8fafc' : '#0f172a';

        oncoChart = new ApexCharts(document.querySelector("#oncoChart"), {{
            chart: {{ type: 'donut', height: 280, foreColor: textColor }},
            labels: ['Oncogênico', 'Outros'],
            series: [0, 0],
            colors: ['#ef4444', '#94a3b8']
        }});
        oncoChart.render();

        geneChart = new ApexCharts(document.querySelector("#geneChart"), {{
            chart: {{ type: 'bar', height: 280, toolbar: {{show:false}}, foreColor: textColor }},
            series: [{{ name: 'Mutas', data: [] }}],
            xaxis: {{ categories: [] }},
            colors: ['#3b82f6'],
            plotOptions: {{ bar: {{ borderRadius: 6, horizontal: true }} }}
        }});
        geneChart.render();

        // Listeners
        $('#sampleSelect, #relevantOnly, #geneSearch').on('change keyup', () => {{
            altTable.draw(); bioTable.draw(); updateStats();
        }});

        // Filtro DataTables
        $.fn.dataTable.ext.search.push((settings, data, idx, row) => {{
            const sample = $('#sampleSelect').val();
            const relevant = $('#relevantOnly').is(':checked');
            const geneSearch = $('#geneSearch').val().toUpperCase();

            const rSample = data[0];
            const rGene = data[1].toUpperCase();
            const rSum = data[4].toLowerCase();

            if (sample && rSample !== sample) return false;
            if (geneSearch && !rGene.includes(geneSearch)) return false;
            if (relevant && settings.nTable.id === 'altTable' && !rSum.includes('oncogenic')) return false;

            return true;
        }});

        updateStats();
    }}

    function updateStats() {{
        const sample = $('#sampleSelect').val();
        const gene = $('#geneSearch').val().toUpperCase();

        // FILTRO DINÂMICO DE VARIANTES
        let fAlt = altData;
        if (sample) fAlt = fAlt.filter(d => d['SAMPLE'] === sample);
        if (gene) fAlt = fAlt.filter(d => d['CGI-Gene'].toUpperCase().includes(gene));

        // FILTRO DINÂMICO DE BIOMARCADORES (Correção aqui!)
        let fBio = bioData;
        if (sample) fBio = fBio.filter(d => d['Sample ID'] === sample);
        if (gene) fBio = fBio.filter(d => d['Alterations'].toUpperCase().includes(gene));

        const oncoCount = fAlt.filter(d => (d['CGI-Oncogenic Summary']||'').toLowerCase().includes('oncogenic')).length;

        // Atualiza Cards
        document.getElementById('stat-total').innerText = fAlt.length;
        document.getElementById('stat-onco').innerText = oncoCount;
        document.getElementById('stat-bio').innerText = fBio.length;

        // Atualiza Gráficos
        oncoChart.updateSeries([oncoCount, fAlt.length - oncoCount]);
        const gFreq = {{}};
        fAlt.forEach(d => gFreq[d['CGI-Gene']] = (gFreq[d['CGI-Gene']] || 0) + 1);
        const top = Object.entries(gFreq).sort((a,b) => b[1]-a[1]).slice(0, 10);
        geneChart.updateOptions({{ xaxis: {{ categories: top.map(x => x[0]) }} }});
        geneChart.updateSeries([{{ data: top.map(x => x[1]) }}]);
    }}

    $(document).ready(init);
</script>
</body>
</html>
"""

# 4. Injeção e Download
final_html = html_template.replace('{json_alt}', json_alt).replace('{json_bio}', json_bio)
with open('dashboard_tcc_filtrado.html', 'w', encoding='utf-8') as f:
    f.write(final_html)
files.download('dashboard_tcc_filtrado.html')
```
