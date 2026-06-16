"""
pubmed_monitor_app.py — PubMed Monitor v7
Seleção visual de periódicos e termos MeSH por especialidade.
Janela padrão: 90–120 dias antes de hoje.

Especialidades:
  - Radiologia Geral, Neuroradiologia, MSK, Ressonância Magnética
  - Ultrassom, IA / Tecnologia, Pediátrico, Vascular / Neurointervenção
  - Alto Impacto Geral, Tórax / Pulmão, Cardiológica, Oncologia
  - Medicina Nuclear / PET, Radiologia de Emergência, Genômica / Radiomômica
  - Odontologia / Radiologia Oral e Maxilofacial
  - Subespecialidades Cirúrgicas, Epidemiologia / Saúde Pública

Enriquecimento:
  - Citações via Semantic Scholar (batch) → OpenAlex (batch) como fallback
  - Open Access + links PMC
  - Painel de métricas enriquecido + ordenação por citações

[NOVO v7]
  - CORREÇÃO: janela_padrao() estava com o corpo quebrado (return inacessível). Reparada.
  - PERFORMANCE: OpenAlex agora é consultado em LOTE (batch de 50 DOIs por chamada),
    eliminando o loop de requisições unitárias com time.sleep por artigo.
  - PERSISTÊNCIA: resultados da busca e do enriquecimento vivem em st.session_state;
    filtros e ordenação não disparam nova busca no PubMed nem "piscam" a tela.
"""

import time
import xml.etree.ElementTree as ET
from datetime import date, timedelta

import pandas as pd
import requests
import streamlit as st
from Bio import Entrez



PERIODICOS = {
    # ── ORIGINAIS ─────────────────────────────────────────────────────────────
    "📰 Radiologia Geral": {
        "Radiology":                     "Radiology",
        "RadioGraphics":                 "Radiographics",
        "AJR (Am J Roentgenol)":         "AJR Am J Roentgenol",
        "European Radiology":            "Eur Radiol",
        "British Journal of Radiology":  "Br J Radiol",
        "Radiologia Brasileira":         "Radiol Bras",
        "Insights into Imaging":         "Insights Imaging",
        "European Journal of Radiology": "Eur J Radiol",
        "Radiologic Clinics N Am":       "Radiol Clin North Am",
        "Korean Journal of Radiology":   "Korean J Radiol",
    },
    "🧠 Neuroradiologia": {
        "AJNR":                          "AJNR Am J Neuroradiol",
        "Neuroradiology":                "Neuroradiology",
        "Neuroimaging Clinics N Am":     "Neuroimaging Clin N Am",
        "Brain":                         "Brain",
        "Neurology":                     "Neurology",
        "Lancet Neurology":              "Lancet Neurol",
        "Neurosurgery":                  "Neurosurgery",
    },
    "🦴 MSK": {
        "Skeletal Radiology":            "Skeletal Radiol",
        "Semin Musculoskeletal Radiol":  "Semin Musculoskelet Radiol",
    },
    "🧲 Ressonância Magnética": {
        "JMRI":                          "J Magn Reson Imaging",
        "Magnetic Resonance in Medicine":"Magn Reson Med",
        "MRI Clinics N Am":              "Magn Reson Imaging Clin N Am",
    },
    "🔊 Ultrassom": {
        "Ultrasound in Medicine & Biology": "Ultrasound Med Biol",
        "Journal of Ultrasound in Medicine": "J Ultrasound Med",
        "Ultrasound in Obstetrics & Gynecology": "Ultrasound Obstet Gynecol",
    },
    "🤖 IA / Tecnologia": {
        "Radiology: Artificial Intelligence":  "Radiol Artif Intell",
        "European Radiology Experimental":     "Eur Radiol Exp",
        "Artificial Intelligence in Medicine": "Artif Intell Med",
        "npj Digital Medicine":                "NPJ Digit Med",
        "Journal of Digital Imaging":          "J Digit Imaging",
    },
    "👶 Pediátrico": {
        "Pediatric Radiology":           "Pediatr Radiol",
        "American Journal of Perinatology": "Am J Perinatol",
    },
    "🫀 Vascular / Neurointervenção": {
        "Stroke (AHA)":                        "Stroke",
        "Journal of Stroke (Korean)":          "J Stroke",
        "Stroke and Vascular Neurology":       "Stroke Vasc Neurol",
        "JVIR":                                "J Vasc Interv Radiol",
        "J NeuroInterventional Surgery":       "J Neurointerv Surg",
        "Circulation":                         "Circulation",
        "European Heart Journal":              "Eur Heart J",
        "Seminars in Neurology (Thieme)":      "Semin Neurol",
    },
    "🌍 Alto Impacto Geral": {
        "Lancet":                          "Lancet",
        "New England Journal of Medicine": "N Engl J Med",
        "JAMA":                            "JAMA",
        "BMJ":                             "BMJ",
        "Nature Medicine":                 "Nat Med",
    },

    # ── NOVOS v5 ──────────────────────────────────────────────────────────────
    "🫁 Tórax / Pulmão": {
        "Journal of Thoracic Imaging":        "J Thorac Imaging",
        "Chest":                              "Chest",
        "Thorax":                             "Thorax",
        "American Journal of Respiratory and Critical Care Medicine": "Am J Respir Crit Care Med",
        "European Respiratory Journal":       "Eur Respir J",
        "Respiratory Medicine":               "Respir Med",
        "Lung Cancer":                        "Lung Cancer",
        "Journal of Thoracic and Cardiovascular Surgery": "J Thorac Cardiovasc Surg",
    },
    "❤️ Cardiológica": {
        "JACC (J Am Coll Cardiol)":           "J Am Coll Cardiol",
        "JACC Cardiovascular Imaging":        "JACC Cardiovasc Imaging",
        "Circulation: Cardiovascular Imaging":"Circ Cardiovasc Imaging",
        "European Heart Journal – Cardiovascular Imaging": "Eur Heart J Cardiovasc Imaging",
        "International Journal of Cardiovascular Imaging": "Int J Cardiovasc Imaging",
        "Journal of Cardiovascular Computed Tomography": "J Cardiovasc Comput Tomogr",
        "Journal of Cardiovascular Magnetic Resonance": "J Cardiovasc Magn Reson",
        "Heart":                              "Heart",
    },
    "🔬 Medicina Nuclear / PET": {
        "Journal of Nuclear Medicine":        "J Nucl Med",
        "European Journal of Nuclear Medicine and Molecular Imaging": "Eur J Nucl Med Mol Imaging",
        "EJNMMI Research":                    "EJNMMI Res",
        "Clinical Nuclear Medicine":          "Clin Nucl Med",
        "Nuclear Medicine and Biology":       "Nucl Med Biol",
        "Seminars in Nuclear Medicine":       "Semin Nucl Med",
    },
    "🏥 Oncologia": {
        "Radiology: Imaging Cancer":          "Radiol Imaging Cancer",
        "Cancer Imaging":                     "Cancer Imaging",
        "Journal of Clinical Oncology":       "J Clin Oncol",
        "Annals of Oncology":                 "Ann Oncol",
        "CA: A Cancer Journal for Clinicians":"CA Cancer J Clin",
        "International Journal of Radiation Oncology Biology Physics": "Int J Radiat Oncol Biol Phys",
        "Radiotherapy and Oncology":          "Radiother Oncol",
        "European Journal of Cancer":         "Eur J Cancer",
    },
    "🩺 Radiologia de Emergência": {
        "Emergency Radiology":                "Emerg Radiol",
        "Annals of Emergency Medicine":       "Ann Emerg Med",
        "Academic Emergency Medicine":        "Acad Emerg Med",
        "Journal of Emergency Medicine":      "J Emerg Med",
        "Injury":                             "Injury",
        "Journal of Trauma and Acute Care Surgery": "J Trauma Acute Care Surg",
    },
    "🧬 Genômica / Radiomômica": {
        "Nature Genetics":                    "Nat Genet",
        "Genome Medicine":                    "Genome Med",
        "European Radiology (Radiomics)":     "Eur Radiol",
        "Cancers (MDPI)":                     "Cancers (Basel)",
        "Frontiers in Oncology":              "Front Oncol",
        "Diagnostics (MDPI)":                 "Diagnostics (Basel)",
        "Quantitative Imaging in Medicine and Surgery": "Quant Imaging Med Surg",
    },
    "🦷 Odontologia / Radiologia Oral e Maxilofacial": {
        "Oral Surgery, Oral Medicine, Oral Pathology and Oral Radiology": "Oral Surg Oral Med Oral Pathol Oral Radiol",
        "Dentomaxillofacial Radiology":       "Dentomaxillofac Radiol",
        "Journal of Oral and Maxillofacial Surgery": "J Oral Maxillofac Surg",
        "Clinical Oral Investigations":       "Clin Oral Investig",
        "Oral Radiology":                     "Oral Radiol",
    },
    "🔪 Subespecialidades Cirúrgicas": {
        "Annals of Surgery":                  "Ann Surg",
        "Surgery":                            "Surgery",
        "Journal of the American College of Surgeons": "J Am Coll Surg",
        "Surgical Endoscopy":                 "Surg Endosc",
        "Journal of Surgical Oncology":       "J Surg Oncol",
        "British Journal of Surgery":         "Br J Surg",
        "World Journal of Surgery":           "World J Surg",
        "Journal of Hepato-Biliary-Pancreatic Sciences": "J Hepatobiliary Pancreat Sci",
    },
    "📊 Epidemiologia / Saúde Pública": {
        "American Journal of Epidemiology":   "Am J Epidemiol",
        "Epidemiology":                       "Epidemiology",
        "International Journal of Epidemiology": "Int J Epidemiol",
        "PLOS Medicine":                      "PLoS Med",
        "Bulletin of the World Health Organization": "Bull World Health Organ",
        "Cadernos de Saúde Pública":          "Cad Saude Publica",
        "Revista de Saúde Pública":           "Rev Saude Publica",
        "Journal of Public Health":           "J Public Health (Oxf)",
    },
}

MESH = {
    # ── ORIGINAIS ─────────────────────────────────────────────────────────────
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

    # ── NOVOS v5 ──────────────────────────────────────────────────────────────
    "🫁 Tórax / Pulmão": [
        "Tomography, X-Ray Computed", "Radiography, Thoracic",
        "Magnetic Resonance Imaging", "Positron-Emission Tomography",
        "Lung Neoplasms", "Pulmonary Nodules", "Lung Diseases",
        "Pulmonary Embolism", "Pneumonia", "COVID-19",
        "Pulmonary Fibrosis", "Interstitial Lung Diseases",
        "Bronchitis, Chronic", "Asthma", "Emphysema",
        "Pleural Effusion", "Pleural Neoplasms", "Pneumothorax",
        "Mediastinal Neoplasms", "Thymus Neoplasms",
        "Bronchiectasis", "Sarcoidosis",
        "Pulmonary Hypertension", "Respiratory Distress Syndrome",
        "Tuberculosis, Pulmonary", "Lung Transplantation",
        "Adenocarcinoma of Lung", "Carcinoma, Non-Small-Cell Lung",
    ],
    "❤️ Cardiológica": [
        "Magnetic Resonance Imaging", "Tomography, X-Ray Computed",
        "Echocardiography", "Positron-Emission Tomography",
        "Coronary Artery Disease", "Myocardial Infarction",
        "Heart Failure", "Cardiomyopathies",
        "Atrial Fibrillation", "Arrhythmias, Cardiac",
        "Aortic Diseases", "Aortic Aneurysm", "Aortic Dissection",
        "Heart Valve Diseases", "Aortic Valve Stenosis",
        "Coronary Angiography", "Cardiac Catheterization",
        "Pericardial Effusion", "Pericarditis",
        "Congenital Heart Defects", "Ventricular Dysfunction",
        "Atherosclerosis", "Plaque, Atherosclerotic",
        "Cardiac MRI", "Coronary CT Angiography",
    ],
    "🔬 Medicina Nuclear / PET": [
        "Positron-Emission Tomography", "Tomography, Emission-Computed, Single-Photon",
        "Radiopharmaceuticals", "Fluorodeoxyglucose F18",
        "Radionuclide Imaging", "Molecular Imaging",
        "Bone Density", "Bone Scintigraphy",
        "Thyroid Neoplasms", "Thyroid Diseases",
        "Neuroendocrine Tumors", "Carcinoid Tumor",
        "PSMA", "Prostate Neoplasms",
        "Lymphoma", "Multiple Myeloma",
        "Amyloid", "Tau Proteins", "Dementia",
        "Myocardial Perfusion Imaging", "Cardiac Scintigraphy",
        "Lutetium", "Yttrium Radioisotopes",
        "Theranostics", "Targeted Radionuclide Therapy",
    ],
    "🏥 Oncologia": [
        "Neoplasms", "Tumor Microenvironment",
        "Lung Neoplasms", "Breast Neoplasms", "Colorectal Neoplasms",
        "Liver Neoplasms", "Pancreatic Neoplasms", "Kidney Neoplasms",
        "Brain Neoplasms", "Prostate Neoplasms", "Uterine Neoplasms",
        "Ovarian Neoplasms", "Lymphoma", "Melanoma",
        "Neoplasm Staging", "Neoplasm Metastasis", "Neoplasm Recurrence, Local",
        "Radiotherapy", "Chemotherapy", "Immunotherapy",
        "Antineoplastic Agents", "Biomarkers, Tumor",
        "Survival Analysis", "Disease-Free Survival",
        "Treatment Outcome", "Neoadjuvant Therapy",
        "Radiation Oncology", "Stereotactic Body Radiotherapy",
    ],
    "🩺 Radiologia de Emergência": [
        "Tomography, X-Ray Computed", "Ultrasonography",
        "Radiography", "Magnetic Resonance Imaging",
        "Wounds and Injuries", "Fractures, Bone",
        "Traumatic Brain Injuries", "Spinal Cord Injuries",
        "Abdominal Injuries", "Thoracic Injuries",
        "Aortic Rupture", "Pneumothorax", "Hemothorax",
        "Appendicitis", "Bowel Obstruction", "Intestinal Perforation",
        "Pulmonary Embolism", "Stroke", "Intracranial Hemorrhages",
        "Foreign Bodies", "Burns",
        "Multiple Trauma", "Polytrauma",
        "Point-of-Care Systems", "Emergency Medical Services",
        "Triage", "Critical Care",
    ],
    "🧬 Genômica / Radiomômica": [
        "Radiomics", "Machine Learning", "Deep Learning",
        "Artificial Intelligence", "Image Processing, Computer-Assisted",
        "Genomics", "Proteomics", "Metabolomics",
        "Biomarkers, Tumor", "Genetic Markers",
        "Tumor Microenvironment", "Tumor Heterogeneity",
        "EGFR", "KRAS", "TP53", "IDH1",
        "DNA Methylation", "Gene Expression Profiling",
        "High-Throughput Nucleotide Sequencing",
        "Prognosis", "Disease-Free Survival", "Overall Survival",
        "Neoplasm Grading", "Pathology, Molecular",
        "Texture Analysis", "Feature Extraction",
    ],
    "🦷 Odontologia / Radiologia Oral e Maxilofacial": [
        "Cone-Beam Computed Tomography", "Radiography, Dental",
        "Radiography, Panoramic", "Radiography, Dental, Digital",
        "Dental Implants", "Alveolar Bone Loss",
        "Periodontal Diseases", "Periodontitis",
        "Dental Caries", "Tooth Diseases",
        "Mouth Neoplasms", "Mandibular Neoplasms", "Maxillary Neoplasms",
        "Jaw Cysts", "Ameloblastoma",
        "Temporomandibular Joint Disorders",
        "Jaw, Edentulous", "Alveolar Ridge Augmentation",
        "Orthognathic Surgery", "Cleft Palate", "Cleft Lip",
        "Salivary Gland Diseases", "Parotid Neoplasms",
        "Facial Bones", "Skull Base",
        "Obstructive Sleep Apnea", "Airway",
    ],
    "🔪 Subespecialidades Cirúrgicas": [
        "General Surgery", "Minimally Invasive Surgical Procedures",
        "Laparoscopy", "Robotic Surgical Procedures",
        "Liver Diseases", "Hepatectomy", "Liver Transplantation",
        "Pancreatic Diseases", "Pancreatectomy", "Pancreatic Neoplasms",
        "Colorectal Surgery", "Colectomy", "Inflammatory Bowel Diseases",
        "Cholecystitis", "Cholelithiasis", "Cholangiocarcinoma",
        "Hernia", "Hernia Repair",
        "Kidney Neoplasms", "Nephrectomy",
        "Thyroid Neoplasms", "Thyroidectomy",
        "Breast Neoplasms", "Mastectomy",
        "Postoperative Complications", "Surgical Wound Infection",
        "Intraoperative Complications",
        "Endoscopy", "Endoscopic Ultrasound",
    ],
    "📊 Epidemiologia / Saúde Pública": [
        "Epidemiology", "Public Health", "Prevalence", "Incidence",
        "Risk Factors", "Cohort Studies", "Case-Control Studies",
        "Cross-Sectional Studies", "Randomized Controlled Trials as Topic",
        "Meta-Analysis as Topic", "Systematic Reviews as Topic",
        "Health Surveys", "Population Surveillance",
        "Mortality", "Morbidity", "Disease Burden",
        "Socioeconomic Factors", "Health Disparities",
        "Health Services Research", "Healthcare Quality",
        "Global Health", "Communicable Diseases",
        "Non-Communicable Diseases", "Chronic Disease",
        "Vaccination", "Infectious Disease Transmission",
        "Environmental Health", "Occupational Diseases",
    ],
}


# ─────────────────────────────────────────────────────────────────────────────
# FUNÇÕES — datas, query, parsing e busca PubMed
# ─────────────────────────────────────────────────────────────────────────────

def janela_padrao():
    """Retorna o intervalo padrão de 90 a 120 dias antes de hoje."""
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
# ENRIQUECIMENTO — Citações (Semantic Scholar batch → OpenAlex batch) + OA + PMC
# ─────────────────────────────────────────────────────────────────────────────

HEADERS_S2  = {"User-Agent": "PubMedMonitor/7.0"}
# Para entrar no "Polite Pool" do OpenAlex (mais rápido/estável), troque o e-mail abaixo:
HEADERS_OA  = {"User-Agent": "PubMedMonitor/7.0 (mailto:seu.email@exemplo.com)"}


# ── Semantic Scholar (batch por DOI) ──────────────────────────────────────────

def _s2_batch(dois: list[str]) -> dict[str, dict]:
    """Busca citações e OA em lote no Semantic Scholar (máx 500 por chamada)."""
    resultado = {}
    for i in range(0, len(dois), 500):
        lote = dois[i : i + 500]
        try:
            r = requests.post(
                "https://api.semanticscholar.org/graph/v1/paper/batch",
                params={"fields": "externalIds,citationCount,isOpenAccess,openAccessPdf"},
                json={"ids": [f"DOI:{d}" for d in lote]},
                headers=HEADERS_S2,
                timeout=20,
            )
            if r.status_code == 200:
                for item in r.json():
                    if not item:
                        continue
                    doi = (item.get("externalIds") or {}).get("DOI", "").lower()
                    if doi:
                        resultado[doi] = {
                            "citacoes":   item.get("citationCount"),
                            "is_oa":      item.get("isOpenAccess", False),
                            "oa_pdf_url": (item.get("openAccessPdf") or {}).get("url", ""),
                            "fonte_cit":  "Semantic Scholar",
                        }
        except Exception:
            pass
        time.sleep(0.3)
    return resultado


# ── OpenAlex (batch por DOI) ──────────────────────────────────────────────────

def _openalex_batch(dois: list[str]) -> dict[str, dict]:
    """
    Busca citações e OA em lote no OpenAlex.
    Batch fixado em 50 DOIs por chamada para evitar HTTP 414 (URI Too Long).
    """
    resultado = {}
    for i in range(0, len(dois), 50):
        lote = dois[i : i + 50]
        filtro_dois = "|".join(lote)
        try:
            r = requests.get(
                "https://api.openalex.org/works",
                params={
                    "filter": f"doi:{filtro_dois}",
                    "select": "doi,cited_by_count,open_access",
                    "per-page": 50,
                },
                headers=HEADERS_OA,
                timeout=20,
            )
            if r.status_code == 200:
                for work in r.json().get("results", []):
                    # OpenAlex devolve o DOI como URL completa (https://doi.org/10.xxx)
                    doi_url = (work.get("doi") or "").lower()
                    if doi_url:
                        doi_clean = doi_url.replace("https://doi.org/", "")
                        oa = work.get("open_access", {})
                        resultado[doi_clean] = {
                            "citacoes":   work.get("cited_by_count"),
                            "is_oa":      oa.get("is_oa", False),
                            "oa_pdf_url": oa.get("oa_url", "") or "",
                            "fonte_cit":  "OpenAlex",
                        }
        except Exception:
            pass
        time.sleep(0.2)
    return resultado


# ── PMC link (batch por PMID) ─────────────────────────────────────────────────

def _pmc_links(pmids: list[str]) -> dict[str, str]:
    """Retorna dict pmid → URL PMC para artigos disponíveis no PMC."""
    links = {}
    for i in range(0, len(pmids), 200):
        lote = pmids[i : i + 200]
        try:
            r = requests.get(
                "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/",
                params={"ids": ",".join(lote), "format": "json", "tool": "PubMedMonitor"},
                timeout=15,
            )
            if r.status_code == 200:
                for rec in r.json().get("records", []):
                    pmcid = rec.get("pmcid", "")
                    pmid  = rec.get("pmid", "")
                    if pmcid and pmid:
                        links[pmid] = f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/"
        except Exception:
            pass
        time.sleep(0.3)
    return links


# ── Função principal ──────────────────────────────────────────────────────────

def enriquecer_df(df: "pd.DataFrame") -> "pd.DataFrame":
    """
    Adiciona ao DataFrame as colunas:
      Citações, Fonte_Cit, Is_OA, URL_OA, URL_PMC, URL_Acesso

    Estratégia de citações:
      1. Semantic Scholar em lote (todos os DOIs)
      2. OpenAlex em lote, apenas para os DOIs que o S2 não encontrou
      3. Mapeamento em memória sobre o DataFrame (sem requisições no laço)
    """
    df = df.copy()

    # 1. Colunas-padrão
    for col, val in [
        ("Citações",   None),
        ("Fonte_Cit",  ""),
        ("Is_OA",      False),
        ("URL_OA",     ""),
        ("URL_PMC",    ""),
        ("URL_Acesso", ""),
    ]:
        if col not in df.columns:
            df[col] = val

    # 2. PMC links (usa PMIDs)
    pmids = df["PMID"].dropna().astype(str).tolist()
    pmc_map = _pmc_links(pmids)
    df["URL_PMC"] = df["PMID"].astype(str).map(pmc_map).fillna("")

    # 3. DOIs válidos
    dois_validos = df["DOI"].dropna().str.lower().str.strip()
    dois_validos = dois_validos[dois_validos != ""].tolist()

    # 4. Semantic Scholar em lote
    s2_map = _s2_batch(dois_validos) if dois_validos else {}

    # 5. OpenAlex em lote — só para os faltantes
    dois_faltantes = [d for d in dois_validos if d not in s2_map]
    oa_map = _openalex_batch(dois_faltantes) if dois_faltantes else {}

    # 6. Unificação (S2 tem prioridade sobre OpenAlex)
    info_map = {**oa_map, **s2_map}

    # 7. Mapeamento direto no DataFrame (sem rede, sem sleeps)
    for idx, row in df.iterrows():
        doi = str(row.get("DOI", "")).lower().strip()
        info = info_map.get(doi)
        if info:
            df.at[idx, "Citações"]  = info["citacoes"]
            df.at[idx, "Fonte_Cit"] = info["fonte_cit"]
            df.at[idx, "Is_OA"]     = info["is_oa"]
            df.at[idx, "URL_OA"]    = info["oa_pdf_url"]

    # 8. URL_Acesso — melhor link disponível: OA PDF > PMC > DOI > PubMed
    def melhor_link(row):
        if row.get("URL_OA"):
            return row["URL_OA"]
        if row.get("URL_PMC"):
            return row["URL_PMC"]
        if row.get("DOI"):
            return f"https://doi.org/{row['DOI']}"
        return row.get("URL", "")

    df["URL_Acesso"] = df.apply(melhor_link, axis=1)

    return df


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

# Persistência de resultados (Fase 1) — sobrevive a interações com filtros
if "df_base" not in st.session_state:
    st.session_state["df_base"] = None
if "df_enrich" not in st.session_state:
    st.session_state["df_enrich"] = None
if "total_resultados" not in st.session_state:
    st.session_state["total_resultados"] = 0
if "query_executada" not in st.session_state:
    st.session_state["query_executada"] = ""

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

def sel_grupo_j(grupo):
    for n, nlm in PERIODICOS[grupo].items():
        st.session_state[f"j_{nlm}"] = True

def des_grupo_j(grupo):
    for n, nlm in PERIODICOS[grupo].items():
        st.session_state[f"j_{nlm}"] = False
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
        with st.expander(grupo, expanded=False):
            # Botões por grupo
            cg1, cg2, cg3 = st.columns([1, 1, 4])
            cg1.button("☑️ Todos", key=f"btn_sel_j_{grupo}",
                       on_click=sel_grupo_j, args=(grupo,), use_container_width=True)
            cg2.button("🔲 Nenhum", key=f"btn_des_j_{grupo}",
                       on_click=des_grupo_j, args=(grupo,), use_container_width=True)

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
        with st.expander(esp, expanded=False):
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


# ─────────────────────────────────────────────────────────────────────────────
# PROCESSAMENTO DO GATILHO DE BUSCA (apenas executa e salva no estado)
# ─────────────────────────────────────────────────────────────────────────────

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
    st.session_state["query_executada"] = query

    with st.spinner("Consultando PubMed…"):
        try:
            df_resultado, total = buscar_pubmed(email_ncbi, query, max_results)
            st.session_state["df_base"]          = df_resultado
            st.session_state["df_enrich"]        = None   # reseta enriquecimento anterior
            st.session_state["total_resultados"] = total
        except Exception as ex:
            st.error(f"Erro PubMed: {ex}")
            st.stop()


# ─────────────────────────────────────────────────────────────────────────────
# RENDERIZAÇÃO PERSISTENTE (observa o estado, não o evento do botão)
# ─────────────────────────────────────────────────────────────────────────────

if st.session_state["df_base"] is not None:

    with st.expander("🔎 Query PubMed gerada", expanded=False):
        st.code(st.session_state["query_executada"], language='text')

    if st.session_state["df_base"].empty:
        st.warning("Nenhum artigo encontrado para este período e seleção.")
        st.stop()

    # DataFrame ativo: enriquecido se existir, senão o base
    enriquecido = st.session_state["df_enrich"] is not None
    df_work = st.session_state["df_enrich"] if enriquecido else st.session_state["df_base"]

    # ── Métricas gerais ───────────────────────────────────────────────────────
    n_rev = int(df_work['_is_review'].sum())
    n_oa  = int(df_work.get("Is_OA", pd.Series(False)).sum()) if enriquecido else None

    st.success(
        f"**{st.session_state['total_resultados']}** artigos no PubMed | "
        f"**{len(df_work)}** retornados | **{n_rev}** Reviews ⭐"
        + (f" | **{n_oa}** Open Access 🔓" if n_oa is not None else "")
    )

    cols_m = st.columns(5 if enriquecido else 4)
    cols_m[0].metric("Total retornado",      len(df_work))
    cols_m[1].metric("Reviews ⭐",           n_rev)
    cols_m[2].metric("Outros",               len(df_work) - n_rev)
    cols_m[3].metric("Periódicos distintos", df_work['Abreviação'].nunique())
    if enriquecido:
        cols_m[4].metric("Open Access 🔓", n_oa)

    if enriquecido and "Citações" in df_work.columns:
        validas = df_work["Citações"].dropna()
        if not validas.empty:
            c1, c2, c3 = st.columns(3)
            c1.metric("Citações — mediana", int(validas.median()))
            c2.metric("Citações — máx.",    int(validas.max()))
            c3.metric("Artigos com dados",  f"{len(validas)}/{len(df_work)}")

    st.divider()

    # ── Controle de enriquecimento ────────────────────────────────────────────
    if not enriquecido:
        st.info(
            "💡 Clique em **Enriquecer** para buscar citações (Semantic Scholar / OpenAlex), "
            "status Open Access e links PMC para os artigos listados."
        )
        if st.button("🔗 Enriquecer selecionados — citações + Open Access",
                     type="secondary", use_container_width=True):
            with st.spinner("Buscando dados nas APIs científicas…"):
                st.session_state["df_enrich"] = enriquecer_df(st.session_state["df_base"])
            st.rerun()
    else:
        st.success("✅ Dados de citação e acesso carregados.")
        if st.button("🔄 Limpar enriquecimento", use_container_width=True):
            st.session_state["df_enrich"] = None
            st.rerun()

    st.divider()

    # ── Filtros (persistentes — não disparam nova busca) ──────────────────────
    st.subheader("📋 Resultados")

    f1, f2, f3, f4 = st.columns([1.2, 1.5, 1.5, 1.5])
    filtro_tipo = f1.selectbox("Tipo", ["Todos", "Apenas Reviews ⭐", "Excluir Reviews"])
    filtro_per  = f2.multiselect("Periódico", sorted(df_work['Abreviação'].unique().tolist()))
    busca_txt   = f3.text_input("Palavra no título")

    opcoes_oa = ["Todos"]
    if enriquecido:
        opcoes_oa += ["Apenas Open Access 🔓", "Apenas PMC disponível"]
    filtro_oa = f4.selectbox("Acesso", opcoes_oa)

    opcoes_ord = ["Padrão (Review > Data)"]
    if enriquecido and "Citações" in df_work.columns:
        opcoes_ord += ["↓ Mais citados", "↑ Menos citados"]
    ordenacao = st.selectbox("Ordenar por", opcoes_ord)

    df_f = df_work.copy()

    if filtro_tipo == "Apenas Reviews ⭐":
        df_f = df_f[df_f['_is_review']]
    elif filtro_tipo == "Excluir Reviews":
        df_f = df_f[~df_f['_is_review']]
    if filtro_per:
        df_f = df_f[df_f['Abreviação'].isin(filtro_per)]
    if busca_txt.strip():
        df_f = df_f[df_f['Título'].str.contains(busca_txt.strip(), case=False, na=False)]
    if enriquecido:
        if filtro_oa == "Apenas Open Access 🔓":
            df_f = df_f[df_f["Is_OA"] == True]
        elif filtro_oa == "Apenas PMC disponível":
            df_f = df_f[df_f["URL_PMC"] != ""]

    if ordenacao == "↓ Mais citados" and enriquecido:
        df_f = df_f.sort_values("Citações", ascending=False, na_position="last")
    elif ordenacao == "↑ Menos citados" and enriquecido:
        df_f = df_f.sort_values("Citações", ascending=True, na_position="last")

    st.caption(f"Exibindo **{len(df_f)}** artigos")

    # ── Tabela ────────────────────────────────────────────────────────────────
    colunas_base   = ['PMID', 'Título', 'Autores', 'Abreviação', 'Data', 'Tipo']
    colunas_enrich = ['Citações', 'Fonte_Cit'] if enriquecido else []
    colunas_link   = ['URL_Acesso'] if enriquecido else ['URL']

    df_exib = df_f[colunas_base + colunas_enrich + colunas_link].copy()
    df_exib['Tipo'] = df_f.apply(
        lambda r: ('⭐ ' if r['_is_review'] else '') + r['Tipo'], axis=1
    )

    if enriquecido:
        df_exib['Título'] = df_f.apply(
            lambda r: ('🔓 ' if r.get('Is_OA') else '') + r['Título'], axis=1
        )

    col_cfg = {
        "Título": st.column_config.TextColumn("Título", width="large"),
        "Tipo":   st.column_config.TextColumn("Tipo", width="medium"),
    }
    if enriquecido:
        col_cfg["Citações"]   = st.column_config.NumberColumn("Cit. 📊", format="%d")
        col_cfg["Fonte_Cit"]  = st.column_config.TextColumn("Fonte", width="small")
        col_cfg["URL_Acesso"] = st.column_config.LinkColumn("Acesso 🔗", display_text="Abrir")
    else:
        col_cfg["URL"] = st.column_config.LinkColumn("Link", display_text="Abrir")

    st.dataframe(
        df_exib,
        use_container_width=True,
        column_config=col_cfg,
        hide_index=True,
    )

    # ── Top 10 mais citados ───────────────────────────────────────────────────
    if enriquecido and "Citações" in df_f.columns:
        top = df_f.dropna(subset=["Citações"]).nlargest(10, "Citações")
        if not top.empty:
            with st.expander("🏆 Top 10 mais citados", expanded=False):
                for _, row in top.iterrows():
                    oa_badge  = "🔓 " if row.get("Is_OA") else ""
                    rev_badge = "⭐ " if row["_is_review"] else ""
                    link = row.get("URL_Acesso") or row.get("URL", "")
                    cit  = int(row["Citações"]) if pd.notna(row["Citações"]) else "–"
                    st.markdown(
                        f"{rev_badge}{oa_badge}**{row['Título'][:120]}**  \n"
                        f"_{row['Autores']}_  |  {row['Abreviação']}  |  "
                        f"{row['Data']}  |  📊 {cit} cit.  "
                        + (f"  |  [Abrir]({link})" if link else "")
                    )

    # ── Abstracts ─────────────────────────────────────────────────────────────
    st.divider()
    if st.toggle("📄 Mostrar abstracts", value=False):
        for _, row in df_f.iterrows():
            oa_badge  = "🔓 " if enriquecido and row.get("Is_OA") else ""
            rev_badge = "⭐ " if row["_is_review"] else ""
            label = rev_badge + oa_badge + row['Título'][:110] + '…'
            with st.expander(label):
                st.markdown(f"**Autores:** {row['Autores']}")
                st.markdown(f"**{row['Periódico']}** | {row['Data']} | {row['Tipo']}")

                links_col = st.columns(3)
                if row.get("DOI"):
                    links_col[0].markdown(f"[🔗 DOI](https://doi.org/{row['DOI']})")
                if enriquecido and row.get("URL_PMC"):
                    links_col[1].markdown(f"[📄 PMC Full Text]({row['URL_PMC']})")
                if enriquecido and row.get("URL_OA"):
                    links_col[2].markdown(f"[🔓 PDF Open Access]({row['URL_OA']})")

                if enriquecido and pd.notna(row.get("Citações")):
                    st.caption(
                        f"📊 {int(row['Citações'])} citações ({row.get('Fonte_Cit', '')})"
                    )
                st.markdown("---")
                st.write(row['Abstract'] or "_Abstract não disponível._")

    # ── Download CSV ──────────────────────────────────────────────────────────
    st.divider()
    cols_dl = ['PMID', 'Título', 'Autores', 'Periódico', 'Abreviação',
               'Data', 'Tipo', 'DOI', 'URL', 'Abstract']
    if enriquecido:
        cols_dl += ['Citações', 'Fonte_Cit', 'Is_OA', 'URL_OA', 'URL_PMC', 'URL_Acesso']

    cols_dl_ok = [c for c in cols_dl if c in df_f.columns]
    csv = df_f[cols_dl_ok].to_csv(index=False, encoding='utf-8-sig').encode('utf-8-sig')
    st.download_button(
        "⬇️ Baixar resultados (.csv)",
        data=csv,
        file_name=f"pubmed_{data_inicio.strftime('%Y%m%d')}_{data_fim.strftime('%Y%m%d')}.csv",
        mime='text/csv',
        use_container_width=True,
    )
