# PyGlider SeaExplorer Project

Un progetto separato per testare PyGlider con i dati del SeaExplorer, completamente indipendente dagli script custom che funzionano già perfettamente.

## Struttura del Progetto

```
PyGlider_SeaExplorer_Project/
├── README.md
├── requirements.txt
├── config/
│   ├── seaexplorer_deployment.yml
│   └── processing_config.yml
├── input/
│   └── raw_data/ (collegamento ai dati originali)
├── output/
│   ├── l0_data/
│   ├── l1_data/
│   └── figures/
├── scripts/
│   ├── 01_setup_pyglider.py
│   ├── 02_raw_to_l0.py
│   ├── 03_l0_to_l1.py
│   └── 04_create_plots.py
└── comparison/
    └── compare_with_custom_scripts.py
```

## Obiettivo

Testare PyGlider come alternativa professionale ai nostri script custom, mantenendo entrambi i progetti separati per confrontare i risultati.

## I Nostri Script Custom

I 4 script custom rimangono intoccabili e continuano a funzionare perfettamente:
- ✅ `1_convert_raw_to_separate_csv.py`
- ✅ `2_merge_mission_data_csv.py` 
- ✅ `3_convert_all_units_csv.py`
- ✅ `4_rename_variables_csv.py`
- ✅ Plus mapping scripts

## PyGlider Test

Questo progetto separato testerà PyGlider per:
- Processamento standardizzato
- Output NetCDF
- Quality control automatico
- Confronto con i nostri risultati

## Installazione

### ⚠️ IMPORTANTE: Ambiente Conda

**Tutti gli script di questo progetto devono essere eseguiti nell'ambiente conda `glider`.**

Prima di eseguire qualsiasi script, attivare sempre l'ambiente:

```bash
conda activate glider
```

L'ambiente `glider` contiene tutti i pacchetti necessari:
- PyGlider (0.0.7)
- xarray
- shapely
- numpy, pandas, matplotlib, ecc.

**NON usare l'ambiente `base`** - non contiene PyGlider e altri pacchetti specifici del progetto.

### Setup Iniziale

Se l'ambiente `glider` non esiste ancora:

```bash
# Creare l'ambiente
conda create -n glider python=3.13

# Attivare l'ambiente
conda activate glider

# Installare i pacchetti
pip install -r requirements.txt
conda install -c conda-forge xarray shapely -y
```

## Uso

**⚠️ Ricorda sempre di attivare l'ambiente prima di eseguire gli script:**

```bash
conda activate glider
```

Poi esegui gli script in ordine:

```bash
# Setup
python scripts/01_setup_pyglider.py

# Processing
python scripts/02_raw_to_l0.py
python scripts/03_l0_to_l1.py
python scripts/04_create_plots.py

# Quality Control
python scripts/qc_variables.py --latest

# Comparison
python comparison/compare_with_custom_scripts.py
```