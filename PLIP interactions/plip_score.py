from pathlib import Path
from lxml import etree
import csv

# Configuration
xml_dir = Path('plip_results')

# reference: plant homolog with PCO
ref_xml = Path('1MVN_report.xml')
ref_ligand_id = "PCO:A:1001"

# reference: human PPCDC with PPC
#ref_xml = Path('1QZU_pred_report.xml')
#ref_ligand_id = "l01:B:1"

plip_ligand_id = "UNL:Z:1"  # ligand ID to extract from PLIP results
output_csv = Path("1MVN_PCO_plip_scores.csv")

def extract_pattern(xmlfile, only_member=None):
    tree = etree.parse(str(xmlfile))
    pattern = set()
    for bs in tree.xpath('//bindingsite'):
        if only_member:
            members = [m.text.strip() for m in bs.xpath('identifiers/members/member')]
            if only_member not in members:
                continue
        for itype in [
            'hydrophobic_interactions', 'hydrogen_bonds', 'water_bridges',
            'salt_bridges', 'pi_stacks', 'pi_cation_interactions',
            'halogen_bonds', 'metal_complexes'
        ]:
            for ia in bs.xpath(f'interactions/{itype}/*'):
                resnr = ia.findtext('resnr')
                if resnr is not None:
                    pattern.add((itype, int(resnr)))
    return pattern

def plip_score(ligand_set, ref_set):
    if not ref_set:
        return 0
    return len(ligand_set & ref_set) / len(ref_set) # features in both reference and docked ligand / features in reference

# === Extract reference pattern ===
ref_pattern = extract_pattern(ref_xml, only_member=ref_ligand_id)

# === Analyze all ligands ===
similar_ligands = []
for xml_file in xml_dir.glob('*.xml'):
    ligand_pattern = extract_pattern(xml_file, only_member=plip_ligand_id)
    score = plip_score(ligand_pattern, ref_pattern) # calculate score
    similar_ligands.append((xml_file.stem, score))

# Sort by plip score descending
similar_ligands.sort(key=lambda x: x[1], reverse=True)


# === Save to CSV ===
with open(output_csv, mode='w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["ID", "PLIP_score"])
    for ligand_name, score in similar_ligands:
        writer.writerow([ligand_name, f"{score:.4f}"])

print(f"\n[INFO] Results saved to {output_csv}")