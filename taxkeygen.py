
from typing import Optional


class Taxon:
    def __init__(self):
        self.name = ""
        self.characters = []


class TraitVariant:
    def __init__(self):
        self.trait = None
        self.id = 0
        self.description = ""
        self.figures = []


class Trait:
    def __init__(self):
        self.title = ""
        self.variants = {}
        self.priority = 1
        self.id = 0

    def __len__(self):
        return len(self.variants)

    def add_variant(self, variant: TraitVariant) -> None:
        self.variants[variant.id] = variant
        variant.trait = self


def read_taxa_data(inname: str) -> dict:
    taxa_dict = {}
    with open(inname, "r") as infile:
        full_file = infile.readlines()
        for line in full_file[1:]:
            data = line.strip().split("\t")
            taxon = Taxon()
            taxa_dict[data[0]] = taxon
            taxon.name = data[0]
            taxon.characters = data[1:]
    return taxa_dict


def read_trait_data(inname: str) -> dict:
    trait_dict = {}
    with open(inname, "r") as infile:
        full_file = infile.readlines()
        for line in full_file[1:]:
            data = line.strip().split("\t")
            if data[0] in trait_dict:
                trait = trait_dict[data[0]]
            else:
                trait = Trait()
                trait_dict[data[0]] = trait
                trait.id = data[0]
                trait.title = data[1]
                trait.priority = int(data[2])
            tv = TraitVariant()
            tv.id = data[3]
            tv.description = data[4]
            tv.figures = data[5]
            trait.add_variant(tv)
    return trait_dict


def match_traits_to_taxa(trait_data: dict, taxa_data: dict) -> None:
    """
    replace trait keys from taxon data file with references to matching trait objects
    """
    for t in taxa_data:
        taxon = taxa_data[t]
        for i, c in enumerate(taxon.characters):
            trait_id, _ = c.split(".")
            trait = trait_data[trait_id]
            taxon.characters[i] = trait.variants[c]


def generate_taxonomic_key(trait_name: str, taxa_name: str, out_name: Optional[str], verbose: bool = True) -> list:
    trait_data = read_trait_data(trait_name)
    if verbose:
        print()
        print("Read {} traits from {}".format(len(trait_data), trait_name))
    taxa_data = read_taxa_data(taxa_name)
    if verbose:
        print("Read {} taxa from {}".format(len(taxa_data), taxa_name))
    match_traits_to_taxa(trait_data, taxa_data)
    output = []
    if out_name is not None:
        with open(out_name, "w") as outfile:
            outfile.writelines(output)
    return output


def input_query(query: str, default: str) -> str:
    x = input("{} [default: {}]: ".format(query, default))
    if x == "":
        x = default
    return x


def main():
    trait_file = input_query("Trait data file", "trait_data.txt")
    taxa_file = input_query("Taxa data file", "taxa_data.txt")
    output_file = input_query("Output HTML file", "output.html")
    generate_taxonomic_key(trait_file, taxa_file, output_file)


if __name__ == "__main__":
    main()
