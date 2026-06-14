"""
pubmed_monitor_app.py — PubMed Monitor v4
Seleção visual de periódicos e termos MeSH por especialidade.
Janela padrão: 90–120 dias antes de hoje.
Botões selecionar/desmarcar todos via session_state.
"""

import time
import xml.etree.ElementTree as ET
from datetime import date, timedelta

import pandas as pd
import streamlit as st
from Bio import Entrez

# ─────────────────────────────────────────────────────────────────────────────
# DADOS EMBUTIDOS
# ─────────────────────────────────────────────────────────────────────────────

PERIODICOS = {
    "📰 Radiologia Geral": {
        "Radiology":                    "Radiology",
        "RadioGraphics":                "Radiographics",
        "AJR (Am J Roentgenol)":        "AJR Am J Roentgenol",
        "European Radiology":           "Eur Radiol",
        "British Journal of Radiology": "Br J Radiol",
        "Radiologia Brasileira":        "Radiol Bras",
        "Insights into Imaging":        "Insights Imaging",
        "European Journal of Radiology":"Eur J Radiol",
        "Radiologic Clinics N Am":      "Radiol Clin North Am",
        "Korean Journal of Radiology":  "Korean J Radiol",
    },
    "🧠 Neuroradiologia": {
        "AJNR":                         "AJNR Am J Neuroradiol",
        "Neuroradiology":               "Neuroradiology",
        "Neuroimaging Clinics N Am":    "Neuroimaging Clin N Am",
        "Brain":                        "Brain",
        "Neurology":                    "Neurology",
        "Lancet Neurology":             "Lancet Neurol",
        "Neurosurgery":                 "Neurosurgery",
    },
    "🦴 MSK": {
        "Skeletal Radiology":           "Skeletal Radiol",
        "Semin Musculoskeletal Radiol": "Semin Musculoskelet Radiol",
    },
    "🧲 Ressonância Magnética": {
        "JMRI":                         "J Magn Reson Imaging",
        "Magnetic Resonance in Medicine":"Magn Reson Med",
        "MRI Clinics N Am":             "Magn Reson Imaging Clin N Am",
    },
    "🔊 Ultrassom": {
        "Ultrasound in Medicine & Biology": "Ultrasound Med Biol",
    },
    "🤖 IA / Tecnologia": {
        "Radiology: Artificial Intelligence": "Radiol Artif Intell",
        "European Radiology Experimental":    "Eur Radiol Exp",
        "Artificial Intelligence in Medicine":"Artif Intell Med",
    },
    "👶 Pediátrico": {
        "Pediatric Radiology":          "Pediatr Radiol",
    },
    "🫀 Vascular / Neurointervenção": {
        "Stroke (AHA)":                       "Stroke",
        "Journal of Stroke (Korean)":         "J Stroke",
        "Stroke and Vascular Neurology":      "Stroke Vasc Neurol",
        "JVIR":                               "J Vasc Interv Radiol",
        "J NeuroInterventional Surgery":      "J Neurointerv Surg",
        "Circulation":                        "Circulation",
        "European Heart Journal":             "Eur Heart J",
        "Seminars in Neurology (Thieme)":     "Semin Neurol",
    },
    "🌍 Alto Impacto Geral": {
        "Lancet":                          "Lancet",
        "New England Journal of Medicine": "N Engl J Med",
    },
}

MESH = {
    "🧠 Neuroradiologia": [
        "Magnetic Resonance Imaging", "Diffusion Magnetic Resonance Imaging",
        "Diffusion Tensor Imaging", "Perfusion Imaging",
        "Spectroscopy, Magnetic Resonance", "Neuroimaging",
        "Tomography, X-Ray Computed", "Positron-Emission Tomography",
        "Cerebral Angiography",
        "Stroke", "Cerebral Infarction", "Intracranial Hemorrhages",
        "Subarachnoid Hemorrhage", "Intracranial Aneurysm",
        "Cerebrovascular Disorders", "Cerebral Small Vessel Diseases",
        "Arteriovenous Malformations",
        "Brain Neoplasms", "Glioma", "Glioblastoma", "Meningioma", "Brain Metastasis",
        "Multiple Sclerosis", "Leukoencephalopathies", "White Matter",
        "Demyelinating Diseases",
        "Dementia", "Alzheimer Disease", "Parkinson Disease",
        "Spinal Cord Diseases", "Spinal Cord Neoplasms",
        "Intervertebral Disc Degeneration",
        "Encephalitis", "Meningitis", "Brain Abscess",
        "Epilepsy", "Brain Malformations",
    ],
    "🦴 MSK": [
        "Magnetic Resonance Imaging", "Ultrasonography, Musculoskeletal",
        "Tomography, X-Ray Computed", "Radiography",
        "Knee Joint", "Shoulder Joint", "Hip Joint", "Ankle Joint",
        "Wrist Joint", "Temporomandibular Joint",
        "Tendons", "Tendon Injuries", "Ligaments, Articular",
        "Rotator Cuff", "Rotator Cuff Injuries",
        "Anterior Cruciate Ligament", "Anterior Cruciate Ligament Injuries",
        "Cartilage, Articular", "Bone Marrow", "Osteoarthritis",
        "Bone Neoplasms", "Fractures, Bone", "Stress Fractures",
        "Osteonecrosis", "Osteoporosis",
        "Intervertebral Disc Degeneration", "Spinal Stenosis",
        "Lumbar Vertebrae", "Cervical Vertebrae", "Spondylolisthesis",
        "Soft Tissue Neoplasms", "Muscle, Skeletal", "Bursa, Synovial", "Synovitis",
        "Arthritis, Rheumatoid", "Spondylarthritis", "Gout",
    ],
    "🤖 IA em Radiologia": [
        "Artificial Intelligence", "Deep Learning", "Machine Learning",
        "Neural Networks, Computer", "Diagnosis, Computer-Assisted",
        "Radiography", "Tomography, X-Ray Computed", "Magnetic Resonance Imaging",
        "Positron-Emission Tomography", "Ultrasonography",
        "Image Interpretation, Computer-Assisted",
        "Image Processing, Computer-Assisted",
        "Radiomics",
        "Sensitivity and Specificity", "ROC Curve",
        "Predictive Value of Tests", "Observer Variation",
        "Reproducibility of Results",
        "Workflow", "Clinical Decision Support Systems",
        "Natural Language Processing", "Radiology Information Systems",
        "Early Detection of Cancer", "Lung Neoplasms", "Breast Neoplasms",
        "Pulmonary Nodules", "Fractures, Bone",
    ],
}

# ─────────────────────────────────────────────────────────────────────────────
# FUNÇÕES
# ─────────────────────────────────────────────────────────────────────────────

def janela_padrao():
    hoje = date.today()
    return hoje - timedelta(days=120), hoje - timedelta(days=90)


def montar_query(termos, journals, d1, d2):
    q_mesh = '(' + ' OR '.join(f'"{t}"[MeSH Terms]' for t in termos) + ')'
    q_jour = '(' + ' OR '.join(f'"{j}"[Journal]' for j in journals) + ')'
    q_data = (f'("{d1.strftime("%Y/%m/%d")}"[Date - Publication] : '
              f'"{d2.strftime("%Y/%m/%d")}"[Date - Publication])')
    return f'{q_mesh} AND {q_jour} AND {q_data}'


def extrair_texto(elem, xpath):
    node = elem.find(xpath)
    return node.text.strip() if node is not None and node.text else ''


def parse_article(art):
    r = {}
    pmid = art.find('.//PMID')
    r['PMID'] = pmid.text if pmid is not None else ''
    title = art.find('.//ArticleTitle')
    r['Título'] = ''.join(title.itertext()) if title is not None else ''
    abs_parts = art.findall('.//AbstractText')
    r['Abstract'] = ' '.join(''.join(a.itertext()) for a in abs_parts)
    r['Periódico'] = extrair_texto(art, './/Journal/Title')
    r['Abreviação'] = extrair_texto(art, './/Journal/ISOAbbreviation')
    pub = art.find('.//PubDate')
    if pub is not None:
        y = extrair_texto(pub, 'Year')
        m = extrair_texto(pub, 'Month')
        r['Data'] = f'{y}/{m}' if m else y
    else:
        r['Data'] = ''
    autores = []
    for au in art.findall('.//Author'):
        last = extrair_texto(au, 'LastName')
        fore = extrair_texto(au, 'ForeName')
        if last:
            autores.append(f'{last} {fore}'.strip())
    r['Autores'] = '; '.join(autores[:6]) + (' et al.' if len(autores) > 6 else '')
    pts = [pt.text for pt in art.findall('.//PublicationType') if pt.text]
    r['Tipo'] = '; '.join(pts)
    r['_is_review'] = any('review' in p.lower() or 'systematic' in p.lower() for p in pts)
    r['DOI'] = ''
    for aid in art.findall('.//ArticleId'):
        if aid.get('IdType') == 'doi':
            r['DOI'] = aid.text
            break
    r['URL'] = f'https://pubmed.ncbi.nlm.nih.gov/{r["PMID"]}/'
    return r


@st.cache_data(show_spinner=False)
def buscar_pubmed(email, query, max_results):
    Entrez.email = email
    handle = Entrez.esearch(db='pubmed', term=query, retmax=max_results, usehistory='y')
    rec = Entrez.read(handle)
    handle.close()
    ids   = rec['IdList']
    total = int(rec['Count'])
    if not ids:
        return pd.DataFrame(), total
    artigos = []
    for inicio in range(0, len(ids), 50):
        lote = ids[inicio:inicio + 50]
        h = Entrez.efetch(db='pubmed', id=','.join(lote), rettype='xml', retmode='xml')
        xml = h.read()
        h.close()
        time.sleep(0.4)
        root = ET.fromstring(xml)
        for art in root.findall('PubmedArticle'):
            artigos.append(parse_article(art))
    df = pd.DataFrame(artigos)
    df = df.sort_values(['_is_review', 'Data'], ascending=[False, False]).reset_index(drop=True)
    return df, total


# ─────────────────────────────────────────────────────────────────────────────
# INICIALIZAÇÃO DO SESSION STATE
# ─────────────────────────────────────────────────────────────────────────────

# Periódicos — padrão: todos marcados
for _grupo, _revistas in PERIODICOS.items():
    for _nome, _nlm in _revistas.items():
        _key = f"j_{_nlm}"
        if _key not in st.session_state:
            st.session_state[_key] = True

# MeSH — padrão: todos marcados
for _esp, _termos in MESH.items():
    for _termo in _termos:
        _key = f"m_{_esp}_{_termo}"
        if _key not in st.session_state:
            st.session_state[_key] = True

# ─────────────────────────────────────────────────────────────────────────────
# CALLBACKS — alteram session_state ANTES da renderização
# ─────────────────────────────────────────────────────────────────────────────

def sel_todos_j():
    for g, r in PERIODICOS.items():
        for n, nlm in r.items():
            st.session_state[f"j_{nlm}"] = True

def des_todos_j():
    for g, r in PERIODICOS.items():
        for n, nlm in r.items():
            st.session_state[f"j_{nlm}"] = False

def sel_todos_m():
    for esp, termos in MESH.items():
        for termo in termos:
            st.session_state[f"m_{esp}_{termo}"] = True

def des_todos_m():
    for esp, termos in MESH.items():
        for termo in termos:
            st.session_state[f"m_{esp}_{termo}"] = False

def sel_esp_m(esp):
    for termo in MESH[esp]:
        st.session_state[f"m_{esp}_{termo}"] = True

def des_esp_m(esp):
    for termo in MESH[esp]:
        st.session_state[f"m_{esp}_{termo}"] = False

# ─────────────────────────────────────────────────────────────────────────────
# UI
# ─────────────────────────────────────────────────────────────────────────────

st.set_page_config(page_title="PubMed Monitor", page_icon="🔬", layout="wide")

st.title("🔬 PubMed Monitor")
st.caption("Busca periódica de artigos — janela padrão: 90 a 120 dias antes de hoje")

# ── PASSO 1: E-mail ───────────────────────────────────────────────────────────
with st.container(border=True):
    st.subheader("① E-mail cadastrado no NCBI")
    st.caption("Obrigatório pela política de uso da API Entrez. "
               "Cadastre-se em https://www.ncbi.nlm.nih.gov/account/")
    email_ncbi = st.text_input(
        "E-mail",
        placeholder="seu.email@exemplo.com",
        label_visibility="collapsed",
    )

st.divider()

# ── PASSO 2: Período ──────────────────────────────────────────────────────────
with st.container(border=True):
    st.subheader("② Período de publicação")
    d_ini_pad, d_fim_pad = janela_padrao()
    modo = st.radio(
        "Modo",
        ["🗓 Padrão — 90 a 120 dias atrás", "✏️ Personalizado"],
        horizontal=True,
    )
    if modo.startswith("🗓"):
        data_inicio, data_fim = d_ini_pad, d_fim_pad
        st.info(
            f"**De:** {data_inicio.strftime('%d/%m/%Y')}  \u00a0\u00a0"
            f"**Até:** {data_fim.strftime('%d/%m/%Y')}"
        )
    else:
        c1, c2 = st.columns(2)
        data_inicio = c1.date_input("De", value=d_ini_pad)
        data_fim    = c2.date_input("Até", value=d_fim_pad)

st.divider()

# ── PASSO 3: Periódicos ───────────────────────────────────────────────────────
with st.container(border=True):
    st.subheader("③ Periódicos")
    st.caption("Sobreposições entre grupos são resolvidas automaticamente.")

    # Botões globais
    bj1, bj2, bj3 = st.columns([1, 1, 4])
    bj1.button("☑️ Selecionar todos", key="btn_sel_j",
               on_click=sel_todos_j, use_container_width=True)
    bj2.button("🔲 Desmarcar todos", key="btn_des_j",
               on_click=des_todos_j, use_container_width=True)

    journals_selecionados = {}
    for grupo, revistas in PERIODICOS.items():
        with st.expander(grupo, expanded=True):
            cols = st.columns(3)
            for i, (nome_exib, nlm) in enumerate(revistas.items()):
                chave = f"j_{nlm}"
                if cols[i % 3].checkbox(nome_exib, key=chave):
                    journals_selecionados[nlm] = True

    journals_lista = list(journals_selecionados.keys())
    st.caption(f"**{len(journals_lista)}** periódico(s) selecionado(s)")

st.divider()

# ── PASSO 4: MeSH ────────────────────────────────────────────────────────────
with st.container(border=True):
    st.subheader("④ Termos MeSH por especialidade")
    st.caption("Termos duplicados entre especialidades são removidos automaticamente.")

    # Botões globais MeSH
    bm1, bm2, bm3 = st.columns([1, 1, 4])
    bm1.button("☑️ Selecionar todos", key="btn_sel_m",
               on_click=sel_todos_m, use_container_width=True)
    bm2.button("🔲 Desmarcar todos", key="btn_des_m",
               on_click=des_todos_m, use_container_width=True)

    mesh_selecionados = {}
    for esp, termos_lista in MESH.items():
        with st.expander(esp, expanded=True):
            # Botões por especialidade
            ce1, ce2, ce3 = st.columns([1, 1, 4])
            ce1.button(f"☑️ Todos",
                       key=f"btn_sel_{esp}",
                       on_click=sel_esp_m, args=(esp,),
                       use_container_width=True)
            ce2.button(f"🔲 Nenhum",
                       key=f"btn_des_{esp}",
                       on_click=des_esp_m, args=(esp,),
                       use_container_width=True)

            cols = st.columns(2)
            for i, termo in enumerate(termos_lista):
                chave = f"m_{esp}_{termo}"
                if cols[i % 2].checkbox(termo, key=chave):
                    mesh_selecionados[termo] = True

    mesh_lista = list(mesh_selecionados.keys())
    st.caption(f"**{len(mesh_lista)}** termo(s) MeSH selecionado(s)")

    # Termos adicionais livres
    with st.expander("➕ Termos adicionais (opcional)", expanded=False):
        extra = st.text_area(
            "Um termo MeSH por linha",
            height=100,
            placeholder="Ex: Lung Neoplasms\nPulmonary Embolism",
        )
        if extra.strip():
            for linha in extra.splitlines():
                t = linha.strip()
                if t and not t.startswith('#'):
                    mesh_selecionados[t] = True
            mesh_lista = list(mesh_selecionados.keys())

st.divider()

# ── PASSO 5: Parâmetros e busca ───────────────────────────────────────────────
with st.container(border=True):
    st.subheader("⑤ Parâmetros e busca")
    max_results = st.slider("Máx. resultados", 50, 500, 200, step=50)
    buscar = st.button("🔍 Buscar no PubMed", type="primary", use_container_width=True)

# ── RESULTADOS ────────────────────────────────────────────────────────────────
if buscar:
    erros = []
    if not email_ncbi or '@' not in email_ncbi:
        erros.append("Informe um e-mail válido no Passo 1.")
    if not journals_lista:
        erros.append("Selecione ao menos um periódico no Passo 3.")
    if not mesh_lista:
        erros.append("Selecione ao menos um termo MeSH no Passo 4.")
    if data_fim < data_inicio:
        erros.append("Data final deve ser posterior à data inicial.")

    for e in erros:
        st.error(e)
    if erros:
        st.stop()

    query = montar_query(mesh_lista, journals_lista, data_inicio, data_fim)

    with st.expander("🔎 Query PubMed gerada", expanded=False):
        st.code(query, language='text')

    with st.spinner("Consultando PubMed…"):
        try:
            df, total = buscar_pubmed(email_ncbi, query, max_results)
        except Exception as ex:
            st.error(f"Erro PubMed: {ex}")
            st.stop()

    if df.empty:
        st.warning("Nenhum artigo encontrado para este período e seleção.")
        st.stop()

    n_rev = int(df['_is_review'].sum())
    st.success(
        f"**{total}** artigos no PubMed | **{len(df)}** retornados | "
        f"**{n_rev}** Reviews ⭐"
    )

    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Total retornado", len(df))
    m2.metric("Reviews ⭐", n_rev)
    m3.metric("Outros", len(df) - n_rev)
    m4.metric("Periódicos distintos", df['Abreviação'].nunique())

    st.divider()
    st.subheader("📋 Resultados")
    f1, f2, f3 = st.columns(3)
    filtro_tipo = f1.selectbox("Tipo", ["Todos", "Apenas Reviews ⭐", "Excluir Reviews"])
    filtro_per  = f2.multiselect("Periódico", sorted(df['Abreviação'].unique().tolist()))
    busca_txt   = f3.text_input("Palavra no título")

    df_f = df.copy()
    if filtro_tipo == "Apenas Reviews ⭐":
        df_f = df_f[df_f['_is_review']]
    elif filtro_tipo == "Excluir Reviews":
        df_f = df_f[~df_f['_is_review']]
    if filtro_per:
        df_f = df_f[df_f['Abreviação'].isin(filtro_per)]
    if busca_txt.strip():
        df_f = df_f[df_f['Título'].str.contains(busca_txt.strip(), case=False, na=False)]

    st.caption(f"Exibindo **{len(df_f)}** artigos")

    df_exib = df_f[['PMID','Título','Autores','Abreviação','Data','Tipo','URL']].copy()
    df_exib['Tipo'] = df_f.apply(
        lambda r: ('⭐ ' if r['_is_review'] else '') + r['Tipo'], axis=1
    )
    st.dataframe(
        df_exib,
        use_container_width=True,
        column_config={
            "URL":    st.column_config.LinkColumn("Link", display_text="Abrir"),
            "Título": st.column_config.TextColumn("Título", width="large"),
            "Tipo":   st.column_config.TextColumn("Tipo", width="medium"),
        },
        hide_index=True,
    )

    st.divider()
    if st.toggle("📄 Mostrar abstracts", value=False):
        for _, row in df_f.iterrows():
            label = ('⭐ ' if row['_is_review'] else '') + row['Título'][:110] + '…'
            with st.expander(label):
                st.markdown(f"**Autores:** {row['Autores']}")
                st.markdown(f"**{row['Periódico']}** | {row['Data']} | {row['Tipo']}")
                if row['DOI']:
                    st.markdown(f"[https://doi.org/{row['DOI']}](https://doi.org/{row['DOI']})")
                st.markdown("---")
                st.write(row['Abstract'] or "_Abstract não disponível._")

    st.divider()
    cols_dl = ['PMID','Título','Autores','Periódico','Abreviação','Data','Tipo','DOI','URL','Abstract']
    csv = df_f[cols_dl].to_csv(index=False, encoding='utf-8-sig').encode('utf-8-sig')
    st.download_button(
        "⬇️ Baixar resultados (.csv)",
        data=csv,
        file_name=f"pubmed_{data_inicio.strftime('%Y%m%d')}_{data_fim.strftime('%Y%m%d')}.csv",
        mime='text/csv',
        use_container_width=True,
    )
