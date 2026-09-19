import zipfile
import json
import re
import pandas as pd

def get_sample_ids_from_qza(qza_path):
    """Extract sample IDs directly from a QIIME 2 feature table .qza file."""
    with zipfile.ZipFile(qza_path, 'r') as z:
        for name in z.namelist():
            if name.endswith('sample-frequency-detail.json') or name.endswith('sample-frequencies.json'):
                with z.open(name) as f:
                    data = json.load(f)
                    return list(data.keys())
    # Fallback using qiime2 artifact view
    import qiime2
    tbl = qiime2.Artifact.load(qza_path).view(pd.DataFrame)
    return list(tbl.index)

def parse_sample_metadata(sample_id):
    """Parse sample naming convention, isolating x/y extra samples."""
    domain = "Bacteria" if sample_id.startswith("B_") else "Fungi"
    code = sample_id[2:]  # Strip 'B_' or 'F_' prefix
    
    # Check if the sample code ends with 'x' or 'y' (e.g., J2x, J2y, J13x)
    is_extra = bool(re.search(r'[xy]$', code))
    sample_type = "Extra_Sample" if is_extra else "Standard"

    # Base Treatment Parsing
    if code.startswith("jM"):
        group = "GJ_Control_Extra" if is_extra else "GJ_Control"
        soil = "Base_Soil"
        inoculum = "MES_Buffer_Alone"
    elif code.startswith("M"):
        group = "GB_Control_Extra" if is_extra else "GB_Control"
        soil = "Base_Soil"
        inoculum = "MES_Buffer_Alone"
    elif code.startswith("B"):
        group = "GB_Treatment_Extra" if is_extra else "GB_Treatment"
        soil = "Gijang_B"
        inoculum = "GB_Soil_Suspension"
    elif code.startswith("J"):
        group = "GJ_Treatment_Extra" if is_extra else "GJ_Treatment"
        soil = "Gyeongju"
        inoculum = "GJ_Soil_Suspension"
    else:
        group = "Unknown"
        soil = "Unknown"
        inoculum = "Unknown"

    # Batch extraction
    if is_extra:
        batch = "GJ3_Extra"
    else:
        match = re.search(r'^(jM\d+|M\d+|B\d+|J\d+)', code)
        batch = match.group(1) if match else "Other"

    return {
        "sample-id": sample_id,
        "Domain": domain,
        "Sample_Type": sample_type,         # 'Standard' vs 'Extra_Sample'
        "Treatment_Group": group,          # Keeps main groups clean
        "Target_Soil": soil,
        "Inoculum": inoculum,
        "Sequencing_Batch": batch,
        "Generation": "P"
    }

# Read sample IDs from both cleaned tables
samples_16S = get_sample_ids_from_qza("05_taxonomy_16S/table_16S_clean.qza")
samples_ITS = get_sample_ids_from_qza("05_taxonomy_ITS/table_ITS_clean.qza")

all_samples = sorted(list(set(samples_16S + samples_ITS)))

# Build DataFrame
rows = [parse_sample_metadata(sid) for sid in all_samples]
df = pd.DataFrame(rows)

# Add QIIME 2 type declaration row (#q2:types)
q2_types = {
    "sample-id": "#q2:types",
    "Domain": "categorical",
    "Sample_Type": "categorical",
    "Treatment_Group": "categorical",
    "Target_Soil": "categorical",
    "Inoculum": "categorical",
    "Sequencing_Batch": "categorical",
    "Generation": "categorical"
}

df_q2 = pd.concat([pd.DataFrame([q2_types]), df], ignore_index=True)

# Export to TSV
output_file = "metadata_P_generation.tsv"
df_q2.to_csv(output_file, sep='\t', index=False)
print(f"Successfully generated '{output_file}' with {len(all_samples)} samples.")