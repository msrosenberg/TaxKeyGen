
from typing import Optional


class Taxon:
    def __init__(self):
        self.name = ""
        self.characters = []

    def find_variant(self, trait):
        t = None
        for c in self.characters:
            if trait == c.trait:
                t = c
        return t


class TraitVariant:
    def __init__(self):
        self.trait = None
        self.id = 0
        self.description = ""
        self.figures = []

    def __lt__(self, other):
        if self.id < other.id:
            return True
        else:
            return False


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


class VariantDist:
    def __init__(self):
        self.trait = None
        self.variants = set()
        self.vfreq = {}
        self.pattern = ""

    def __len__(self):
        return len(self.variants)

    def __repr__(self):
        freqs = []
        for v in self.vfreq:
            freqs.append("{}: {}".format(v.description, self.vfreq[v]))
        return "Variant Distribution of Trait: {}\n\tFreqs: {}\n\rPattern: {}".format(self.trait.title,
                                                                                      ", ".join(freqs),
                                                                                      self.pattern)

    def add_variant(self, variant):
        self.variants.add(variant)
        if variant in self.vfreq:
            self.vfreq[variant] += 1
        else:
            self.vfreq[variant] = 1


def sorted_taxa_keys(tax_dict: dict) -> list:
    return sorted(tax_dict.keys())


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


def determine_variant_freqs(taxa_data: dict, trait_data: dict) -> list:
    """
    determine which variants are present in the taxa for each trait and determine the frequency
    of each as well
    """
    var_freqs = []
    for t_key in trait_data:
        trait = trait_data[t_key]
        var_dist = VariantDist()
        var_dist.trait = trait
        for tax_key in sorted_taxa_keys(taxa_data):
            taxon = taxa_data[tax_key]
            tax_var = taxon.find_variant(trait)
            var_dist.add_variant(tax_var)
        var_freqs.append(var_dist)
    return var_freqs


def filter_var_freqs(var_freqs: list) -> list:
    """
    filter out traits without exactly two variants
    """
    new_freqs = []
    for vf in var_freqs:
        if len(vf) == 2:
            new_freqs.append(vf)
    return new_freqs


def determine_var_pattern(var_freqs: list, taxa_data: dict) -> None:
    """
    determine the pattern of variants across taxa

    at this stage only binary traits are examined, so the patten is stored as a string of 0's and 1's
    """
    for vf in var_freqs:
        variants = sorted(vf.variants)
        for tk in sorted_taxa_keys(taxa_data):
            taxon = taxa_data[tk]
            tv = taxon.find_variant(vf.trait)
            vf.pattern += str(variants.index(tv))


def reverse_pattern(x: str) -> str:
    """
    returns the inverse of a binary string

    thus 00011100 -> 11100011
    """
    x.replace("0", "a")
    x.replace("1", "0")
    return x.replace("a", "1")


def cluster_traits(var_freqs: list) -> list:
    clusters = [[var_freqs[0]]]  # seed with first variant
    for vf in var_freqs[1:]:
        match = False
        for c in clusters:
            cpattern = c[0].pattern
            if (cpattern == vf.pattern) or (cpattern == reverse_pattern(vf.pattern)):
                c.append(vf)
                match = True
        if not match:
            clusters.append([vf])
    return clusters


def create_key_tree(taxa_data: dict, trait_data: dict):
    var_freqs = determine_variant_freqs(taxa_data, trait_data)
    var_freqs = filter_var_freqs(var_freqs)
    determine_var_pattern(var_freqs, taxa_data)
    trait_clusters = cluster_traits(var_freqs)
    # for i, c in enumerate(trait_clusters):
    #     print("Cluster", i)
    #     for cc in c:
    #         print(cc)
    #     print()


def generate_taxonomic_key(trait_name: str, taxa_name: str, out_name: Optional[str], verbose: bool = True) -> list:
    trait_data = read_trait_data(trait_name)
    if verbose:
        print()
        print("Read {} traits from {}".format(len(trait_data), trait_name))
    taxa_data = read_taxa_data(taxa_name)
    if verbose:
        print("Read {} taxa from {}".format(len(taxa_data), taxa_name))
    match_traits_to_taxa(trait_data, taxa_data)
    key_tree = create_key_tree(taxa_data, trait_data)
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
