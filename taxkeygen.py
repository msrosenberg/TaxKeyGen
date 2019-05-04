
from typing import Optional, Tuple, Union


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

    def __str__(self):
        return self.description


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
        self.cluster = None

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

    def add_to_cluster(self, cluster):
        cluster.append(self)
        self.cluster = cluster

    def freq_ratio(self) -> float:
        c0 = self.pattern.count("0")
        c1 = self.pattern.count("1")
        n = c0 + c1
        return min(c0/n, c1/n)

    def __lt__(self, other):
        if self.trait.priority < other.trait.priority:
            return True
        elif self.trait.priority > other.trait.priority:
            return False
        elif self.freq_ratio() < other.freq_ratio():
            return True
        elif self.freq_ratio() > other.freq_ratio():
            return False
        elif len(self.cluster) < len(other.cluster):
            return True
        else:
            return False


class KeyNode:
    def __init__(self):
        self.parent = None
        # self.children = []
        self.child0 = None
        self.child1 = None
        self.child0variants = []
        self.child1variants = []
        self.traits = []
        self.number = 0

    def new_child_node(self, c):
        child = KeyNode()
        if c == 0:
            self.child0 = child
        else:
            self.child1 = child
        # self.children.append(child)
        child.parent = self
        return child


class TraitFigure:
    def __init__(self, image: str, caption: str):
        self.image = image
        self.caption = caption


def extract_figures(x: str) -> list:
    figs = []
    if x != ".":
        fs = x.split("||")
        for f in fs:
            fig, cap = f.split("|")
            figs.append(TraitFigure(fig, cap))
    return figs


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
            tv.figures = extract_figures(data[5])
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
    x = x.replace("0", "a")
    x = x.replace("1", "0")
    return x.replace("a", "1")


def cluster_traits(var_freqs: list) -> list:
    # seed with first variant
    vf = var_freqs[0]
    c = []
    vf.add_to_cluster(c)
    clusters = [c]
    for vf in var_freqs[1:]:
        match = False
        for c in clusters:
            cpattern = c[0].pattern
            if (cpattern == vf.pattern) or (cpattern == reverse_pattern(vf.pattern)):
                vf.add_to_cluster(c)
                match = True
        if not match:
            c = []
            clusters.append(c)
            vf.add_to_cluster(c)
    return clusters


def split_taxa(taxa_data: dict, key_vf: VariantDist) -> Tuple[list, list]:
    """
    split taxa into two groups based on the key variant
    """
    taxa0 = []
    taxa1 = []
    variants = sorted(key_vf.variants)
    for tk in sorted_taxa_keys(taxa_data):
        taxon = taxa_data[tk]
        tv = taxon.find_variant(key_vf.trait)
        i = variants.index(tv)
        if i == 0:
            taxa0.append(taxon)
        else:
            taxa1.append(taxon)
    return taxa0, taxa1


def trait_list(cluster: list) -> list:
    tlist = []
    for vf in cluster:
        tlist.append(vf.trait)
    return tlist


def match_variants(t_list: list, taxon: Taxon):
    return [taxon.find_variant(t) for t in t_list]


def create_key_tree(taxa_data: dict, trait_data: dict, node: KeyNode):
    """
    primary key tree creation function
    """
    var_freqs = determine_variant_freqs(taxa_data, trait_data)
    var_freqs = filter_var_freqs(var_freqs)
    if len(var_freqs) > 0:
        determine_var_pattern(var_freqs, taxa_data)
        cluster_traits(var_freqs)
        var_freqs.sort(reverse=True)
        key_vf = var_freqs[0]
        taxa0, taxa1 = split_taxa(taxa_data, key_vf)
        node.traits = trait_list(key_vf.cluster)  # assign all traits which align with this split to node
        node.child0variants = match_variants(node.traits, taxa0[0])
        node.child1variants = match_variants(node.traits, taxa1[0])
        if len(taxa0) > 1:
            new_node = node.new_child_node(0)
            create_key_tree({t.name: t for t in taxa0}, trait_data, new_node)
        else:
            node.child0 = taxa0[0]
        if len(taxa1) > 1:
            new_node = node.new_child_node(1)
            create_key_tree({t.name: t for t in taxa1}, trait_data, new_node)
        else:
            node.child1 = taxa1[0]
    else:  # taxa are undivisible
        if node == node.parent.child0:
            node.parent.child0 = taxa_data
        else:
            node.parent.child1 = taxa_data


def number_nodes(tree: KeyNode, node_number: int) -> int:
    """
    number the tree nodes

    this is used for the numbering rules in the key, 'if xxxxx, go to NUMBER'
    """
    tree.number = node_number
    if isinstance(tree.child0, KeyNode):
        node_number = number_nodes(tree.child0, node_number+1)
    if isinstance(tree.child1, KeyNode):
        node_number = number_nodes(tree.child1, node_number+1)
    return node_number


def start_output(output: list, nnodes: int) -> None:
    output.append("<html>\n")
    output.append("  <head>\n")
    output.append("    <style>\n")
    output.append("      .tree-grid {\n")
    output.append("                    display: grid;\n")
    output.append("                    grid-template-columns: 4ch 1fr;\n")
    output.append("                    grid-template-areas: \"fork-number fork-option\";\n")
    output.append("                    grid-row-gap: 10px;\n")
    output.append("                 }\n")
    for n in range(nnodes):
        output.append("      #key-fork-n-{} {{grid-area: {} / fork-number }}\n".format(n+1, n*2 + 1))
        output.append("      #key-fork-a-{} {{grid-area: {} / fork-option }}\n".format(n+1, n*2 + 1))
        output.append("      #key-fork-b-{} {{grid-area: {} / fork-option }}\n".format(n+1, (n+1)*2))
    output.append("      .key-fork-n { font-weight: bold }\n")
    output.append("      .key-fork-a { padding-left: 3ch; text-indent: -3ch }\n")
    output.append("      .key-fork-b { padding-left: 3ch; text-indent: -3ch; padding-bottom: 1em }\n")
    output.append("      .variant-fig { margin: 0.25em; display: inline-block; text-indent: 0 }\n")
    output.append("      .variant-fig img { height: 200px }\n")
    output.append("      .variant-fig figcaption { font-style: italic; font-size: 0.75em; text-align: center; "
                  "margin-top: 0.25em }\n")
    output.append("    </style>\n")
    output.append("  </head>\n")
    output.append("  <body>\n")
    output.append("    <div class=\"tree-grid\">\n")


def end_output(output: list) -> None:
    output.append("    </div>\n")
    output.append("  </body>\n")
    output.append("</html>\n")


def get_var_figs(variants: list) -> str:
    fig_str = ""
    for v in variants:
        for f in v.figures:
            fig_str += "<figure class=\"variant-fig\"><img src=\"images/{}\" />".format(f.image)
            if f.caption != ".":
                fig_str += "<figcaption>{}</figcaption>".format(f.caption)
            fig_str += "</figure> "
    return fig_str

def write_key(tree: KeyNode, output: list) -> None:
    # # --- original output code ---
    # def fork_str(letter: str, variants: list, tip: Union[KeyNode, Taxon]):
    #     var_strs = [str(v) for v in variants]
    #     outstr = "    <p>{}. ".format(letter) + "; ".join(var_strs) + ". &mdash; "
    #     if isinstance(tip, KeyNode):
    #         outstr += "<a href=\"#key-node-{0}\">Go to {0}</a>".format(tip.number)
    #     elif isinstance(tip, Taxon):
    #         outstr += tip.name
    #     elif isinstance(tip, dict):
    #         pass
    #         taxa_names = sorted_taxa_keys(tip)
    #         if len(taxa_names) == 2:
    #             outstr += taxa_names[0] + " or " + taxa_names[1]
    #         else:
    #             outstr += ", ".join(taxa_names[:-1]) + ", or " + taxa_names[len(taxa_names)]
    #     else:
    #         print("ERROR: Child node of invalid type:", tip)
    #     return outstr + "</p>\n"
    #
    # output.append("    <p><a name=\"key-node-{0}\">{0}.</a></p>\n".format(tree.number))
    # output.append(fork_str("a", tree.child0variants, tree.child0))
    # output.append(fork_str("b", tree.child1variants, tree.child1))
    # output.append("    <p>&nbsp;</p>\n")
    # if isinstance(tree.child0, KeyNode):
    #     write_key(tree.child0, output)
    # if isinstance(tree.child1, KeyNode):
    #     write_key(tree.child1, output)

    def fork_str(letter: str, variants: list, tip: Union[KeyNode, Taxon], n: int):
        var_strs = [str(v) for v in variants]
        var_figs = get_var_figs(variants)
        outstr = "    <div id=\"key-fork-{0}-{1}\" class=\"key-fork-{0}\">{0}. ".format(letter, n) + \
                 "; ".join(var_strs) + ". &mdash; "
        if isinstance(tip, KeyNode):
            outstr += "<a href=\"#key-node-{0}\">Go to {0}</a>".format(tip.number)
        elif isinstance(tip, Taxon):
            outstr += tip.name
        elif isinstance(tip, dict):
            pass
            taxa_names = sorted_taxa_keys(tip)
            if len(taxa_names) == 2:
                outstr += taxa_names[0] + " or " + taxa_names[1]
            else:
                outstr += ", ".join(taxa_names[:-1]) + ", or " + taxa_names[len(taxa_names)]
        else:
            print("ERROR: Child node of invalid type:", tip)
        if var_figs != "":
            outstr += "<br/>" + var_figs
        return outstr + "</div>\n"

    output.append("    <div id=\"key-fork-n-{0}\" class=\"key-fork-n\">"
                  "<a name=\"key-node-{0}\">{0}.</a></div>\n".format(tree.number))
    output.append(fork_str("a", tree.child0variants, tree.child0, tree.number))
    output.append(fork_str("b", tree.child1variants, tree.child1, tree.number))
    if isinstance(tree.child0, KeyNode):
        write_key(tree.child0, output)
    if isinstance(tree.child1, KeyNode):
        write_key(tree.child1, output)


def generate_taxonomic_key(trait_name: str, taxa_name: str, out_name: Optional[str], verbose: bool = True) -> list:
    trait_data = read_trait_data(trait_name)
    if verbose:
        print()
        print("Read {} traits from {}".format(len(trait_data), trait_name))
    taxa_data = read_taxa_data(taxa_name)
    if verbose:
        print("Read {} taxa from {}".format(len(taxa_data), taxa_name))
    match_traits_to_taxa(trait_data, taxa_data)
    key_tree = KeyNode()
    create_key_tree(taxa_data, trait_data, key_tree)
    total_nodes = number_nodes(key_tree, 1)
    output = []
    start_output(output, total_nodes)
    write_key(key_tree, output)
    end_output(output)
    if out_name is not None:
        with open(out_name, "w") as outfile:
            outfile.writelines(output)
        if verbose:
            print("Key written to {}".format(out_name))
    if verbose:
        print("Finished")
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
