from pathlib import Path
from collections import Counter, defaultdict
from lxml import etree
import matplotlib.pyplot as plt

class XMLStorage:
    @staticmethod
    def getdata(tree, xpath, force_string=False):
        found = tree.xpath(f"{xpath}/text()")
        if not found:
            return None
        val = found[0]
        if force_string:
            return val
        if val in ('True','False'):
            return val == 'True'
        try:
            return int(val)
        except ValueError:
            try:
                return float(val)
            except ValueError:
                return val

def parse_interactions_by_type(xmlfile, only_member=None):
    tree = etree.parse(str(xmlfile))
    interaction_data = defaultdict(list)
    for bs in tree.xpath('//bindingsite'):
        if only_member:
            members = [m.text.strip() for m in bs.xpath('identifiers/members/member')]
            if only_member not in members:
                continue

        for interaction_type in [
            'hydrophobic_interactions',
            'hydrogen_bonds',
            'water_bridges',
            'salt_bridges',
            'pi_stacks',
            'pi_cation_interactions',
            'halogen_bonds',
            'metal_complexes'
        ]:
            for ia in bs.xpath(f'interactions/{interaction_type}/*'):
                resnr_nodes = ia.xpath('resnr/text()')
                if resnr_nodes:
                    resnr = int(resnr_nodes[0])
                    interaction_data[interaction_type].append(resnr)
    return interaction_data


def parse_residue_names(xmlfile, only_member=None):
    tree = etree.parse(str(xmlfile))
    resname_map = defaultdict(dict)
    for bs in tree.xpath('//bindingsite'):
        if only_member:
            members = [m.text.strip() for m in bs.xpath('identifiers/members/member')]
            if only_member not in members:
                continue

        for interaction_type in [
            'hydrophobic_interactions', 'hydrogen_bonds', 'water_bridges',
            'salt_bridges', 'pi_stacks', 'pi_cation_interactions',
            'halogen_bonds', 'metal_complexes'
        ]:
            for ia in bs.xpath(f'interactions/{interaction_type}/*'):
                resnr = XMLStorage.getdata(ia, 'resnr')
                restype = XMLStorage.getdata(ia, 'restype', force_string=True)
                if resnr is not None and restype:
                    resname_map[interaction_type][resnr] = restype # maps each residue to interaction types
    return resname_map

# === Configuration ===
# plant homolog PCO reference
xml_dir = Path('plip_results') # folder with plip results
ref_xml = Path('1MVN_report.xml') # plip xml file of reference
ref_ligand_id = "PCO:A:1001" # reference ID
plip_ligand_id = "UNL:Z:1" # unnamed ID within docking

# human PPCDC PPC reference
# ref_xml = Path('1QZU_pred_report.xml')
# ref_ligand_id = "l01:B:1"

# 6AIM PPC reference
#ref_xml = Path('6AIM_report.xml')
#ref_ligand_id = "9Z3:A:401"


output_dir = Path('sorted_plots_with_reference')
output_dir.mkdir(exist_ok=True)

# === Parse all non-reference xmls ===
interaction_counters = defaultdict(Counter)
for xml_file in xml_dir.glob('*.xml'):
    interactions = parse_interactions_by_type(xml_file, only_member=plip_ligand_id) # parses docked ligands interactions
    for itype, residues in interactions.items():
        interaction_counters[itype].update(residues) # counts residues of specific interactions

# === Parse reference xml ===
reference_residues = defaultdict(set)
ref_interactions = parse_interactions_by_type(ref_xml, only_member=ref_ligand_id) # parses reference interactions
for itype, residues in ref_interactions.items():
    reference_residues[itype].update(residues) # is not a counter -> only registers the residues (not how often)

print(f"\n[INFO] Reference residues for ligand {ref_ligand_id}:")
for itype, resset in reference_residues.items():
    print(f"  {itype}: {sorted(resset)}")

"""
output substrate interaction:    
[INFO] Reference residues for ligand PCO:A:1001:
hydrophobic_interactions: [31, 34]
hydrogen_bonds: [174, 183]
water_bridges: [34]
"""

# === residue names ===
residue_names_all = defaultdict(dict)
for xml_file in xml_dir.glob('*.xml'):
    names = parse_residue_names(xml_file, only_member=plip_ligand_id) # get residue names for all ligands
    for itype in names:
        residue_names_all[itype].update(names[itype]) # safe residues for each interaction type

residue_names_ref = parse_residue_names(ref_xml, only_member=ref_ligand_id) # get residue names for reference
residue_labels = defaultdict(dict)
# gets overlap between total residues and reference residues
for itype in set(residue_names_all) | set(residue_names_ref):
    residue_labels[itype] = {**residue_names_all.get(itype, {}), **residue_names_ref.get(itype, {})}

# === plots ===
for itype in sorted(set(interaction_counters.keys()) | set(reference_residues.keys())):
    counter_all = interaction_counters.get(itype, Counter())
    top_all = set([res for res, _ in counter_all.most_common(20)]) # only the top 20 residues 
    top_ref = set(reference_residues.get(itype, [])) 
    combined_residues = sorted(top_all | top_ref) # bars that are both in docked ligands and reference

    if not combined_residues:
        continue

    freqs_all = [counter_all.get(res, 0) for res in combined_residues]
    labels = [f"{residue_labels[itype].get(r, 'UNK')} {r}" for r in combined_residues] # residue labels on x-axis, UNK if unknown

    # Color bars: red if in reference, else light blue
    colors = ["red" if res in top_ref else "lightblue" for res in combined_residues]

    plt.figure(figsize=(17, 4))
    plt.bar(combined_residues, freqs_all, color=colors)
    plt.xticks(combined_residues, labels, rotation=90, ha='center', fontsize=6)
    plt.xlabel('Residue (Amino Acid and Position)')
    plt.ylabel('Number of Interactions')
    plt.title(f'Interaction Type: {itype.replace("_", " ").title()}')
    plt.tight_layout()
    plt.savefig(output_dir / f'{itype}_ref_colored.png')
    plt.show()
