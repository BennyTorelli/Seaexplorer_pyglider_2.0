#!/usr/bin/env python3
"""
PYGLIDER SEAEXPLORER - PIPELINE COMPLETA
Script unico che processa tutti i dati SeaExplorer con PyGlider
"""

import os
import sys
import gzip
import shutil
import re
import xarray as xr
import numpy as np
import pandas as pd
import pyglider.seaexplorer as seaexplorer
import pyglider.ncprocess as ncprocess
from pathlib import Path

class SeaExplorerProcessor:
    """Processore completo per dati SeaExplorer con PyGlider"""
    
    def __init__(self):
        """Inizializza percorsi e configurazione"""
        self.rawdir = 'input/raw/'
        self.sanitized_dir = 'input/sanitized/'
        self.deploymentyaml = 'config/seaexplorer_0067.yml'
        
        # Output directories
        self.rawncdir = 'output/l0_data/rawnc/'
        self.l0_tsdir = 'output/l0_data/timeseries/'
        self.profiles_dir = 'output/l1_data/profiles/'
        self.grid_dir = 'output/l1_data/gridded/'
        self.analysis_dir = 'output/analysis/'
        
        # Crea tutte le directory
        for directory in [self.rawncdir, self.l0_tsdir, 
                         self.profiles_dir, self.grid_dir, self.analysis_dir, self.sanitized_dir]:
            os.makedirs(directory, exist_ok=True)
    
    def run_complete_pipeline(self):
        """Esegue l'intera pipeline PyGlider"""
        
        print("PYGLIDER SEAEXPLORER - PIPELINE COMPLETA")
        print("=" * 47)
        print("Elaborazione completa dei dati SeaExplorer con PyGlider")
        print()
        
        # STEP 0: Decompressione file .gz (se necessaria)
        print("STEP 0: Verifica e decompressione file .gz")
        print("-" * 42)
        decompressed = self._decompress_gz_files()
        if decompressed > 0:
            print(f" {decompressed} file decompressi")
        else:
            print("  Nessun file .gz trovato o già decompresso")
        
        print()
        
        # STEP 1: Raw to parquet
        if not self._check_parquet_exists():
            print("STEP 1: Conversione raw → parquet")
            print("-" * 35)
            success1 = self._convert_raw_to_parquet()
            if not success1:
                print("Pipeline fermata - errore Step 1")
                return False
        else:
            print("STEP 1: File parquet già esistenti - saltato")
        
        print()
        
        # STEP 2: Merge parquet
        print(" STEP 2: Merge file parquet")
        print("-" * 27)
        success2 = self._merge_parquet_files()
        print()
        
        # STEP 3: L0 Timeseries (siempre regenera)
        print("STEP 3: Creazione L0 timeseries")
        print("-" * 33)
        
        # Elimina file L0 esistente se presente
        existing_l0 = self._get_l0_file()
        if existing_l0 and os.path.exists(existing_l0):
            os.remove(existing_l0)
            print(f"File L0 esistente rimosso: {os.path.basename(existing_l0)}")
        
        l0_file = self._create_l0_timeseries()
        if not l0_file:
            print(" Pipeline fermata - errore Step 3")
            return False
        
        # STEP 3b: Applica conversioni personalizzate ai dati
        print("\nSTEP 3b: Applicazione conversioni personalizzate")
        print("-" * 47)
        self._apply_custom_conversions(l0_file)

        # Non-invasive QC: produce a separate CSV for temperature QC
        # DISABILITATO - usa lo script qc_variables.py separatamente quando necessario
        # try:
        #     from scripts.qc_variables import range_qc_temperature, export_qc_to_csv
        #     print()
        #     print("STEP 3b: Generazione CSV QC temperatura (separato)")
        #     print("-" * 45)
        #     try:
        #         ds_l0 = xr.open_dataset(l0_file)
        #         ds_l0 = range_qc_temperature(ds_l0)
        #         qc_csv = f"{self.analysis_dir}seaexplorer_qc_temperature.csv"
        #         export_qc_to_csv(ds_l0, varname='temperature', csv_path=qc_csv)
        #         print(f"CSV QC temperatura creato: {qc_csv}")
        #         ds_l0.close()
        #     except Exception as e:
        #         print(f"Impossibile generare CSV QC temperatura: {e}")
        # except Exception:
        #     # If the qc module isn't available, skip quietly
        #     print("Modulo QC non disponibile; salto generazione CSV QC temperatura")
        
        print()
        
        # STEP 4: Conversione CSV (siempre regenera)
        csv_file = f"{self.analysis_dir}seaexplorer_data_complete.csv"
        print("STEP 4: Conversione NetCDF → CSV")
        print("-" * 33)
        
        # Elimina CSV esistente se presente
        if os.path.exists(csv_file):
            os.remove(csv_file)
            print(f"File CSV esistente rimosso: {os.path.basename(csv_file)}")
        
        self._convert_to_csv(l0_file, csv_file)
        
        print()
        
        # STEP 5: Analisi qualità
        print("STEP 5: Analisi qualità dati")
        print("-" * 30)
        self._analyze_data_quality(l0_file)
        
        print()
        
        # STEP 6: Rinomina variabili con nomi standard
        print("STEP 6: Rinomina variabili con nomi standard")
        print("-" * 45)
        standardized_files = self._standardize_variable_names(l0_file, csv_file)
        
        print()
        print("PIPELINE COMPLETATA CON SUCCESSO!")
        print("RISULTATI DISPONIBILI:")
        print(f"L0 timeseries: {l0_file}")
        print(f"CSV completo: {csv_file}")
        if standardized_files:
            print("File con nomi standard:")
            for std_file in standardized_files:
                print(f"      - {std_file}")
        
        return True
    
    def _check_parquet_exists(self):
        """Controlla se i file parquet esistono già"""
        return len([f for f in os.listdir(self.rawncdir) if f.endswith('.parquet')]) > 0
    
    def _decompress_gz_files(self):
        """Decomprime tutti i file .gz nella directory raw"""
        decompressed_count = 0
        
        # Controlla la directory raw principale
        if not os.path.exists(self.rawdir):
            return decompressed_count
                
        # Trova tutti i file .gz nella directory raw
        gz_files = [f for f in os.listdir(self.rawdir) if f.endswith('.gz')]
        
        for gz_file in gz_files:
            gz_path = os.path.join(self.rawdir, gz_file)
            
            # Nome del file decompresso (rimuovi .gz)
            decompressed_name = gz_file[:-3]  # Rimuovi .gz
            decompressed_path = os.path.join(self.rawdir, decompressed_name)
            
            # Decomprime solo se il file non esiste già
            if not os.path.exists(decompressed_path):
                try:
                    with gzip.open(gz_path, 'rb') as f_in:
                        with open(decompressed_path, 'wb') as f_out:
                            shutil.copyfileobj(f_in, f_out)
                    decompressed_count += 1
                    print(f"       {gz_file} → {decompressed_name}")
                except Exception as e:
                    print(f"       Errore decomprimendo {gz_file}: {e}")
        
        return decompressed_count
    
    def _sanitize_raw_logs(self, input_dir, output_dir, timestamp_col='PLD_REALTIMECLOCK'):
        """Copia i file da input_dir a output_dir, sostituendo valori vuoti/invalidi
        della colonna timestamp specificata con un placeholder sicuro così il parsing
        (es. di Polars) non fallisce.
        """
        Path(output_dir).mkdir(parents=True, exist_ok=True)

        # regex per corrispondere dd/mm/YYYY HH:MM:SS(.ms)
        ts_re = re.compile(r"^\d{1,2}/\d{1,2}/\d{4}\s+\d{1,2}:\d{2}:\d{2}(?:\.\d{1,6})?$")
        fixed_count = 0

        for fname in os.listdir(input_dir):
            src = os.path.join(input_dir, fname)
            dest = os.path.join(output_dir, fname)

            # salta directory e file compressi
            if os.path.isdir(src) or fname.endswith('.gz'):
                continue

            try:
                with open(src, 'r', encoding='utf-8', errors='replace') as fi:
                    lines = fi.readlines()
            except Exception:
                # fallback: copia file binari o non leggibili
                try:
                    shutil.copy2(src, dest)
                except Exception:
                    pass
                continue

            if not lines:
                open(dest, 'w').close()
                continue

            header = lines[0].rstrip('\n\r')
            cols = header.split(';')
            if timestamp_col not in cols:
                # copia così com'è
                with open(dest, 'w', encoding='utf-8') as fo:
                    fo.writelines(lines)
                continue

            idx = cols.index(timestamp_col)
            fixed = [header + '\n']
            file_modified = False

            for line in lines[1:]:
                row = line.rstrip('\n\r')
                parts = row.split(';')
                if len(parts) <= idx:
                    parts += [''] * (idx + 1 - len(parts))
                rawval = parts[idx].strip()
                if not rawval or rawval.upper() in ('NULL', 'NONE', '-'):
                    parts[idx] = '01/01/1971 00:00:00.000'
                    file_modified = True
                else:
                    if ts_re.match(rawval):
                        parts[idx] = rawval
                    else:
                        parts[idx] = '01/01/1971 00:00:00.000'
                        file_modified = True

                fixed.append(';'.join(parts) + '\n')

            with open(dest, 'w', encoding='utf-8') as fo:
                fo.writelines(fixed)
            
            if file_modified:
                fixed_count += 1

        return fixed_count
    
    def _convert_raw_to_parquet(self):
        """Converte file raw in parquet"""
        try:
            success = False
            
            # Processa tutti i file dalla directory raw
            if os.path.exists(self.rawdir):
                # Conta tutti i file non compressi (sia payload che navigazione)
                all_files = [f for f in os.listdir(self.rawdir) if not f.endswith('.gz')]
                pld_files = [f for f in all_files if not f.endswith('.gli.sub.')]
                nav_files = [f for f in all_files if f.endswith('.gli.sub.')]
                
                if all_files:
                    print(f"   Processando {len(all_files)} file totali ({len(pld_files)} payload, {len(nav_files)} navigazione)...")
                    
                    # Sanitize raw files to avoid datetime conversion failures
                    print(f"   Sanitizing files per timestamp vuoti...")
                    fixed_count = self._sanitize_raw_logs(self.rawdir, self.sanitized_dir, timestamp_col='PLD_REALTIMECLOCK')
                    if fixed_count > 0:
                        print(f"       {fixed_count} file corretti per timestamp vuoti")
                    else:
                        print(f"        Nessun file necessitava correzioni timestamp")
                    
                    # Ora chiama pyglider sulla directory sanitizzata invece che su raw directory
                    success = seaexplorer.raw_to_rawnc(
                        self.sanitized_dir, self.rawncdir, self.deploymentyaml
                    )
                else:
                    print(f"     Nessun file decompresso trovato in {self.rawdir}")
            else:
                print(f"     Directory {self.rawdir} non trovata")
            
            parquet_count = len([f for f in os.listdir(self.rawncdir) if f.endswith('.parquet')])
            print(f"Step 1 completato: {parquet_count} file parquet creati totali")
            
            return success
            
        except Exception as e:
            print(f"Step 1 errore: {e}")
            return False
    
    def _merge_parquet_files(self):
        """Merge dei file parquet"""
        try:
            success = seaexplorer.merge_parquet(
                self.rawncdir, self.rawncdir, self.deploymentyaml
            )
            print(f"Step 2 completato: file parquet uniti")
            return success
        except Exception as e:
            print(f"Step 2 warning: {e}")
            return False
    
    def _create_l0_timeseries(self):
        """Crea L0 timeseries"""
        try:
            l0_file = seaexplorer.raw_to_timeseries(
                self.rawncdir, self.l0_tsdir, self.deploymentyaml, kind='raw'
            )
            
            if l0_file and os.path.exists(l0_file):
                # No additional lat/lon conversions: rely on pyglider/netcdf data
                pass
                size = os.path.getsize(l0_file) / (1024*1024)
                print(f"Step 3 completato: {os.path.basename(l0_file)} ({size:.1f} MB)")
                return l0_file
            return None
            
        except Exception as e:
            print(f"Step 3 errore: {e}")
            return None

    def _apply_custom_conversions(self, nc_file):
        """Applica conversioni personalizzate ai dati nel NetCDF"""
        try:
            import xarray as xr
            import gsw
            import numpy as np
            
            ds = xr.open_dataset(nc_file)
            modified = False
            
            # Conversione TURB: backscatter_700 → NTU
            # Formula: NTU = beta(m-1sr-1) / 0.002727
            if 'backscatter_700' in ds:
                conversion_factor = 1.0 / 0.002727  # = 366.8329
                ds['backscatter_700'] = ds['backscatter_700'] * conversion_factor
                ds['backscatter_700'].attrs['units'] = 'NTU'
                ds['backscatter_700'].attrs['long_name'] = 'turbidity from 700 nm backscatter'
                ds['backscatter_700'].attrs['conversion_applied'] = f'NTU = beta / 0.002727'
                print(f"   Conversione TURB: backscatter → NTU (fattore {conversion_factor:.4f})")
                modified = True
            
            # Conversione CNDC: mS/cm → S/m
            # Formula: S/m = mS/cm / 10
            if 'conductivity' in ds:
                ds['conductivity'] = ds['conductivity'] / 10.0
                ds['conductivity'].attrs['units'] = 'S m-1'
                ds['conductivity'].attrs['long_name'] = 'water conductivity'
                ds['conductivity'].attrs['conversion_applied'] = 'S/m = mS/cm / 10'
                print(f"   Conversione CNDC: mS/cm → S/m (divisione per 10)")
                modified = True
            
            # REMOVED: Erroneous mmol/L → µmol/L conversion
            # PyGlider already provides DOXY in µmol/L (confirmed by RAW data analysis)
            # The low values (0.3-1 µmol/L) are data errors/interpolation artifacts, not mmol/L
            # These should be flagged by QC, not "fixed" with x1000 multiplication
            
            # Ricalcolo SALINITÀ con conducibilità corretta
            # Usa gsw (TEOS-10) per calcolare salinità pratica
            if 'salinity' in ds and 'conductivity' in ds and 'temperature' in ds and 'pressure' in ds:
                # Estrai variabili necessarie
                cndc = ds['conductivity'].values  # Ora in S/m
                temp = ds['temperature'].values   # °C
                pres = ds['pressure'].values      # dbar
                
                # Converti conducibilità S/m → mS/cm per gsw (gsw vuole mS/cm!)
                cndc_mscm = cndc * 10.0
                
                # Calcola salinità pratica con gsw
                # SP = gsw.SP_from_C(C, t, p) dove C in mS/cm
                salinity_new = gsw.SP_from_C(cndc_mscm, temp, pres)
                
                # Sostituisci nel dataset
                ds['salinity'] = (ds['salinity'].dims, salinity_new)
                ds['salinity'].attrs['units'] = '1'
                ds['salinity'].attrs['long_name'] = 'practical salinity'
                ds['salinity'].attrs['standard_name'] = 'sea_water_practical_salinity'
                ds['salinity'].attrs['conversion_applied'] = 'Recalculated with gsw.SP_from_C using corrected conductivity'
                
                # Statistiche per verifica
                sal_valid = salinity_new[~np.isnan(salinity_new)]
                print(f"   Ricalcolo SALINITÀ: min={np.min(sal_valid):.2f}, median={np.median(sal_valid):.2f}, max={np.max(sal_valid):.2f} PSU")
                modified = True
            
            # Ricalcolo DENSITÀ con salinità corretta
            # Usa gsw per calcolare densità in-situ
            if 'density' in ds and 'salinity' in ds and 'temperature' in ds and 'pressure' in ds:
                sal = ds['salinity'].values  # PSU (practical salinity)
                temp = ds['temperature'].values  # °C
                pres = ds['pressure'].values  # dbar
                
                # Calcola salinità assoluta (SA) da salinità pratica (SP)
                # Serve lat/lon per correzione geografica
                lat = ds['latitude'].values if 'latitude' in ds else np.full_like(sal, 28.6)  # Default La Palma
                lon = ds['longitude'].values if 'longitude' in ds else np.full_like(sal, -17.9)
                
                SA = gsw.SA_from_SP(sal, pres, lon, lat)
                
                # Calcola temperatura conservativa
                CT = gsw.CT_from_t(SA, temp, pres)
                
                # Calcola densità in-situ (kg/m³)
                density_new = gsw.rho(SA, CT, pres)
                
                # Sostituisci nel dataset
                ds['density'] = (ds['density'].dims, density_new)
                ds['density'].attrs['units'] = 'kg m-3'
                ds['density'].attrs['long_name'] = 'in-situ density'
                ds['density'].attrs['standard_name'] = 'sea_water_density'
                ds['density'].attrs['conversion_applied'] = 'Recalculated with gsw.rho using corrected salinity'
                
                dens_valid = density_new[~np.isnan(density_new)]
                print(f"   Ricalcolo DENSITÀ: min={np.min(dens_valid):.2f}, median={np.median(dens_valid):.2f}, max={np.max(dens_valid):.2f} kg/m³")
                modified = True
                
                # Conversione OSSIGENO: µmol/L → µmol/kg usando densità corretta
                if 'oxygen_concentration' in ds:
                    doxy_umol_L = ds['oxygen_concentration'].values  # µmol/L from PyGlider
                    
                    # Conversione: µmol/kg = µmol/L / (densità_kg/m³ / 1000)
                    # Equivale a: µmol/kg = µmol/L * 1000 / densità_kg/m³
                    doxy_umol_kg = doxy_umol_L * 1000.0 / density_new
                    
                    # Sostituisci nel dataset
                    ds['oxygen_concentration'] = (ds['oxygen_concentration'].dims, doxy_umol_kg)
                    ds['oxygen_concentration'].attrs['units'] = 'umol kg-1'
                    ds['oxygen_concentration'].attrs['long_name'] = 'oxygen concentration'
                    ds['oxygen_concentration'].attrs['conversion_applied'] = 'Converted from µmol/L to µmol/kg using TEOS-10 in-situ density'
                    
                    doxy_valid = doxy_umol_kg[~np.isnan(doxy_umol_kg)]
                    print(f"   Conversione DOXY: µmol/L → µmol/kg, range {np.min(doxy_valid):.1f}-{np.max(doxy_valid):.1f} µmol/kg")
                    modified = True
            
            # Chiudi prima di riaprire per scrittura
            if modified:
                # Salva in un file temporaneo
                temp_file = nc_file + '.tmp'
                ds.to_netcdf(temp_file)
                ds.close()
                
                # Sostituisci il file originale
                import shutil
                shutil.move(temp_file, nc_file)
                print(f"   NetCDF aggiornato con conversioni")
            else:
                ds.close()
            
        except Exception as e:
            print(f"   Avviso conversioni: {e}")

    def _convert_to_csv(self, nc_file, csv_file):
        """Converte NetCDF in CSV con conversione coordinate automatica"""
        try:
            ds = xr.open_dataset(nc_file)
            df = ds.to_dataframe().reset_index()
            
            # CONVERSIONE AUTOMATICA COORDINATE
            # Ensure we keep original coordinate columns from dataset if present
            # Some datasets use coords rather than data_vars; promote them if needed
            if 'latitude' not in df.columns and 'latitude' in ds.coords:
                df['latitude'] = ds.latitude.values
            if 'longitude' not in df.columns and 'longitude' in ds.coords:
                df['longitude'] = ds.longitude.values

            # Salva CSV con coordinate disponibili e NaN espliciti
            df.to_csv(csv_file, index=False, na_rep='NaN')
            
            size = os.path.getsize(csv_file) / (1024*1024)
            coord_info = ""
            if 'latitude' in df.columns and 'longitude' in df.columns:
                coord_info = " (con coordinate dal NetCDF)"
            print(f"Step 4 completato: CSV creato ({size:.1f} MB){coord_info}")
            
            ds.close()
            
        except Exception as e:
            print(f"Step 4 error: {e}")
    
    def _standardize_variable_names(self, l0_file, csv_file):
        """Rinomina variabili con nomi standard specificati"""
        standardized_files = []
        
        # Mappatura nomi variabili originali → standard
        variable_mapping = {
            'temperature': 'TEMP',
            'time': 'TIME', 
            'oxygen_concentration': 'DOXY',
            'backscatter_700': 'TURB',
            'chlorophyll': 'CHLA',
            'cdom': 'CDOM',
            'conductivity': 'CNDC',
            'pressure': 'PRES',
            'latitude': 'LATITUDE',
            'longitude': 'LONGITUDE',
            'salinity': 'PSAL',
            'depth': 'DEPTH'  # Aggiungiamo depth se viene calcolata
        }
        
        # Processa file NetCDF (solo L0)
        nc_files = []
        if l0_file and os.path.exists(l0_file):
            nc_files.append(('L0', l0_file))
            
        for level, nc_file in nc_files:
            try:
                # Prima aggiungi depth calcolata
                self._add_depth_variable(nc_file)
                # Poi rinomina variabili
                std_file = self._rename_variables_in_netcdf(nc_file, variable_mapping, level)
                if std_file:
                    standardized_files.append(std_file)
            except Exception as e:
                print(f"   Errore rinominando {level}: {e}")
        
        # Processa file CSV
        if csv_file and os.path.exists(csv_file):
            try:
                std_csv = self._rename_variables_in_csv(csv_file, variable_mapping)
                if std_csv:
                    standardized_files.append(std_csv)
            except Exception as e:
                print(f"   Errore rinominando CSV: {e}")
                
        return standardized_files
    
    def _add_depth_variable(self, nc_file):
        """Aggiunge variabile depth calcolata da pressure"""
        try:
            # Carica dataset in memoria prima di modificare per evitare file lock
            ds = xr.open_dataset(nc_file)
            try:
                ds.load()
            finally:
                try:
                    ds_mem = xr.Dataset.from_dict(ds.to_dict())
                finally:
                    try:
                        ds.close()
                    except Exception:
                        pass

            ds = ds_mem

            # Controlla se depth esiste già
            if 'depth' in ds.data_vars:
                return

            # Controlla se pressure esiste
            if 'pressure' not in ds.data_vars:
                print(f"   Pressure non trovata, impossibile calcolare depth")
                return

            # Calcolo depth da pressure
            pressure = ds.pressure
            depth = pressure.copy()
            depth.values = pressure.values.copy()

            # Imposta attributi depth
            depth.attrs = {
                'units': 'm',
                'long_name': 'depth below sea surface',
                'standard_name': 'depth',
                'positive': 'down',
                'comment': 'Calculated from pressure using approximation: depth(m) ≈ pressure(dbar)',
                'accuracy': 'Approximation valid for ocean waters, ±1% accuracy'
            }

            # Aggiungi depth al dataset
            ds = ds.assign(depth=depth)

            # Sovrascrivi il file con depth aggiunta (scriviamo il dataset in memoria)
            try:
                ds.to_netcdf(nc_file, mode='w')
            except PermissionError as e:
                outdir = os.path.dirname(nc_file) or '.'
                alt = os.path.join(outdir, Path(nc_file).stem + '.depth.nc')
                try:
                    ds.to_netcdf(alt, mode='w')
                    print(f"   Wrote depth to alternate file due to permission: {alt}")
                except Exception as e2:
                    print(f"   Failed writing depth to alternate file: {e2}")

            print(f"   DEPTH aggiunta: range {float(depth.min().values):.1f}-{float(depth.max().values):.1f} m")
            
        except Exception as e:
            print(f"   Errore calcolo depth: {e}")
    
    def _rename_variables_in_netcdf(self, input_file, mapping, level):
        """Rinomina variabili in un file NetCDF"""
        try:
            # Nome file standardizzato
            base_name = os.path.splitext(os.path.basename(input_file))[0]
            output_file = f"{self.analysis_dir}{base_name}_standard_names.nc"
            
            # Elimina file esistente se presente (per forzare rigenerazione)
            if os.path.exists(output_file):
                os.remove(output_file)
                print(f"   Eliminato file esistente: {os.path.basename(output_file)}")
            
            # Apri dataset
            ds = xr.open_dataset(input_file)
            
            # Rinomina variabili presenti
            renamed_vars = {}
            available_vars = list(ds.data_vars) + list(ds.coords)
            
            for old_name, new_name in mapping.items():
                if old_name in available_vars:
                    renamed_vars[old_name] = new_name
            
            if renamed_vars:
                # Applica rinominazioni
                ds_renamed = ds.rename(renamed_vars)
                
                # Salva file con nomi standard
                ds_renamed.to_netcdf(output_file)
                
                size = os.path.getsize(output_file) / (1024*1024)
                print(f"   {level} rinominato: {os.path.basename(output_file)} ({size:.1f} MB)")
                print(f"      Variabili rinominate: {len(renamed_vars)}")
                for old, new in renamed_vars.items():
                    print(f"         {old} → {new}")
                
                ds_renamed.close()
            else:
                print(f"   {level}: Nessuna variabile da rinominare trovata")
                output_file = None
            
            ds.close()
            return output_file
            
        except Exception as e:
            print(f"   Errore NetCDF {level}: {e}")
            return None
    
    def _rename_variables_in_csv(self, input_csv, mapping):
        """Rinomina colonne in un file CSV"""
        try:
            # Nome file standardizzato
            base_name = os.path.splitext(os.path.basename(input_csv))[0]
            output_csv = f"{self.analysis_dir}{base_name}_standard_names.csv"
            
            # Elimina file esistente se presente (per forzare rigenerazione)
            if os.path.exists(output_csv):
                os.remove(output_csv)
                print(f"   Eliminato file esistente: {os.path.basename(output_csv)}")
            
            # Leggi CSV
            df = pd.read_csv(input_csv)

            # Applica la mappatura completa (incluse coordinate)
            rename_dict = {}
            for old_name, new_name in mapping.items():
                if old_name in df.columns:
                    rename_dict[old_name] = new_name

            if rename_dict:
                # Applica rinominazioni a tutte le colonne (comprese coordinate)
                df_renamed = df.rename(columns=rename_dict)

                # Salva CSV con nomi standard e NaN espliciti
                df_renamed.to_csv(output_csv, index=False, na_rep='NaN')

                size = os.path.getsize(output_csv) / (1024*1024)
                print(f"   CSV rinominato: {os.path.basename(output_csv)} ({size:.1f} MB)")
                print(f"      Colonne rinominate: {len(rename_dict)}")
                for old, new in rename_dict.items():
                    print(f"         {old} → {new}")

                return output_csv
            else:
                # Anche se non ci sono colonne da rinominare, assicurati che le colonne coordinate
                # originali e le versioni decimali/uppercase siano presenti nel CSV di output
                if 'LATITUDE' in df.columns and 'latitude' not in df.columns:
                    df['latitude'] = df['LATITUDE']
                if 'LONGITUDE' in df.columns and 'longitude' not in df.columns:
                    df['longitude'] = df['LONGITUDE']

                # Se non esiste una colonna LATITUDE/LONGITUDE ma esistono latitude/longitude, promuovile
                if 'latitude' in df.columns and 'LATITUDE' not in df.columns:
                    df['LATITUDE'] = df['latitude']
                if 'longitude' in df.columns and 'LONGITUDE' not in df.columns:
                    df['LONGITUDE'] = df['longitude']

                # Salva il CSV (anche se nulla da rinominare) per coerenza con NaN espliciti
                df.to_csv(output_csv, index=False, na_rep='NaN')
                print(f"   CSV: Nessuna colonna da rinominare trovata; creato output con colonne coordinate preservate")
                return output_csv
                
        except Exception as e:
            print(f"   Errore CSV: {e}")
            return None
    
    def _analyze_data_quality(self, l0_file):
        """Analizza la qualità dei dati"""
        try:
            print("Analisi L0:")
            self._analyze_single_file(l0_file, "L0")
                
        except Exception as e:
            print(f"Analisi errore: {e}")
    
    def _analyze_single_file(self, nc_file, level):
        """Analizza un singolo file NetCDF"""
        ds = xr.open_dataset(nc_file)
        
        print(f"   Campioni: {len(ds.time):,}")
        
        if 'salinity' in ds:
            sal = ds.salinity.values
            sal_valid = sal[~np.isnan(sal)]
            outliers = np.sum(sal_valid > 50)
            print(f"   Salinita: {sal_valid.min():.2f}-{sal_valid.max():.2f} PSU ({outliers:,} outlier)")
            
        if 'temperature' in ds:
            temp = ds.temperature.values
            temp_valid = temp[~np.isnan(temp)]
            print(f"   Temperatura: {temp_valid.min():.2f}-{temp_valid.max():.2f} degC")
        
        ds.close()
    
    def _get_l0_file(self):
        """Trova il file L0 esistente"""
        files = list(Path(self.l0_tsdir).glob("*.nc"))
        return str(files[0]) if files else None


def main():
    """Funzione principale"""
    print("Inizializzazione processore SeaExplorer...")
    processor = SeaExplorerProcessor()
    
    success = processor.run_complete_pipeline()
    
    if success:
        print("\nMISSIONE COMPLETATA!")
        print("Tutti i dati SeaExplorer sono stati processati con PyGlider")
    else:
        print("\nERRORE NELLA PIPELINE")
        sys.exit(1)


if __name__ == '__main__':
    main()