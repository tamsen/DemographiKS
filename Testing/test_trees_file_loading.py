import os
import string
import unittest
import msprime
import numpy as np
import tskit
from matplotlib import pyplot as plt


class TestTreesFile(unittest.TestCase):


    def test_getting_started(self):

        test_output_folder="../test_out"
        if not os.path.exists(test_output_folder):
            os.makedirs(test_output_folder)

        test_output_folder
        pop_size=10_000
        seq_length=10_000_000

        sweep_model = msprime.SweepGenicSelection(
            position=seq_length/2, start_frequency=0.0001, end_frequency=0.9999, s=0.25, dt=1e-6)

        ts = msprime.sim_ancestry(
            20,
            model=[sweep_model, msprime.StandardCoalescent()],
            population_size=pop_size,
            sequence_length=seq_length,
            recombination_rate=1e-8,
            #recombination_rate=0, <- the trees go away if there is no recombination
            random_seed=1234,  # only needed for repeatabilty
            )
        # Optionally add finite-site mutations to the ts using the Jukes & Cantor model, creating SNPs
        ts = msprime.sim_mutations(ts, rate=1e-8, random_seed=4321)

        #Since we simulated the ancestry of 20 diploid individuals, 
        #our tree sequence contains 40 sample nodes, one for each genome.
        
        #each tiny part of the genome will be represented by a tree,
        #ie, the topology for how the first 10 bases of genomes 1-40 are connected.
        
        print("Full treeseq")
        for tree in ts.trees():
            print(f"Tree {tree.index} covers {tree.interval}")
            if tree.index >= 4:
                print("...")
                break
        print(f"Tree {ts.last().index} covers {ts.last().interval}")

        #if we only care about genomes 4 and 13, we could do:
        reduced_ts = ts.simplify([4,13])  # simplify to the first 10 samples
        print("Reduced treeseq: 4 and 13 only")
        for tree in reduced_ts.trees():
            print(f"Tree {tree.index} covers {tree.interval}")
            tree_txt = tree.draw_text()
            mrca = ts.node(tree.root)
            print("mrca time: {0}".format(mrca.time))
            print(tree_txt)
            if tree.index >= 4:
                print("...")
                break
            #if tree.index == 1: #first locus?
            #    tree_txt=tree.draw_text()
            #    print(tree_txt)

        print(f"Tree {ts.last().index} covers {ts.last().interval}")
        print("Note the intervals got larger..")

        #reduced_ts.first().draw_text()

        swept_tree = ts.at(5_000_000)  # or you can get e.g. the nth tree using ts.at_index(n)
        intvl = swept_tree.interval
        print(f"Tree number {swept_tree.index}, which runs from position {intvl.left} to {intvl.right}:")
        # Draw it at a wide size, to make room for all 40 tips
        tree_pic= swept_tree.draw_svg(size=(1000, 200))
        svg_out=os.path.join(test_output_folder,"tree.svg")
        with open(svg_out, "w") as f:
            f.write(tree_pic)
        tree_txt = swept_tree.draw_text()
        print(tree_txt)
        
        #get tree for every gene
        num_genes=int(5_000_000/1000)
        gene_starts=[i*1000 for i in range(0,num_genes)]
        gene_centers=[g+500 for g in gene_starts]
        mrcas=[]
        for gene_center in gene_centers:
            tree_for_gene=reduced_ts.at(gene_center)
            mrca = reduced_ts.node(tree_for_gene.root)
            mrcas.append(mrca.time)
            print("gene {0}, mrca time: {1}".format(gene_center,mrca.time))

        
        fig = plt.figure(figsize=(10, 10), dpi=350)
        x = mrcas
        png_out=os.path.join(test_output_folder,"mrca_hist.png")
        label="hist to test coalescence plotting"
        print(label)
        max_mrca= max(mrcas)
        bin_size=1000
        bins = np.arange(0, max_mrca, bin_size)
        n, bins, patches = plt.hist(x, bins=bins, facecolor='b', alpha=0.25,
                                        label='label', density=False)

        print('foo')
        plt.xlim([0, max_mrca * (1.1)])

        plt.title(label)
        #plt.xlim([0, max_mrca])
        #plt.ylim([0, 300])
        plt.legend()
        plt.xlabel("MRCA time")
        plt.ylabel("# genes in bin")
    
        plt.savefig(png_out)
        plt.clf()
        plt.close()
        
    def test_load_one_tree_file(self):

        out_dir="/home/tamsen/Data/DemographiKS_output/DemographiKS_m12d12y2024_h16m41s25/SLiM_output"
        trees_file = os.path.join("test_trees_file.txt")
        trees_file = os.path.join(out_dir,"allotetraploid_bottleneck_trees_at_div.txt")

        ts = tskit.load(trees_file)
        reduced_ts = ts.simplify([0,1])  # simplify to the first 10 samples


        for tree in reduced_ts.trees():
            print(tree.index)
            if tree.has_multiple_roots:
                #if False:
                print("not coalesced")
            else:
                mrca = ts.node(tree.root)
                print("mrca time: {0}".format(mrca.time))
                tree.draw_text()


        #reduced_ts.draw_svg()
        #reduced_ts.draw_text()
        print("foo")

    def test_loading_tress_file(self):

        trees_file = os.path.join("test_trees_file.txt")
        ts = tskit.load(trees_file)
        metadata = ts.metadata["SLiM"]
        # log.write_to_log("SLiM metadata dict:\t" + str(metadata))
        print("size SLiM population:\t" + str((ts.individuals_population.size)))
        print("size SLiM samples:\t" + str((ts.num_samples)))

        #genomes
        samp_ids = ts.samples()
        print("  ID of diploid individual: ", " ".join([f"{ts.node(s).individual:3}" for s in samp_ids]))
        print("       ID of (sample) node: ", " ".join([f"{s:3}" for s in samp_ids]))
        for v in ts.variants():
            site = v.site
            alleles = np.array(v.alleles)
            print(f"Site {site.id} (ancestral state '{site.ancestral_state}')",  alleles[v.genotypes])
            if site.id >= 4:  # only print up to site ID 4
                print("...")
            break

        for sample in samp_ids:
            individual= ts.node(sample).individual
            print(individual)
            
        for tree in ts.trees():

            if tree.has_multiple_roots:
                print("Tree {0} has not coalesced".format(tree.index))
            else:
                print("Tree {0} has coalesced".format(tree.index))
                genome_interval=tree.interval
                individual= ts.node(tree.root).individual
                mrca = ts.node(tree.root)
                print("mrca time: {0}".format(mrca.time))
                print("interval:\t" + str(tree.interval))
                print("individual:\t" + str(individual))
        #else:
        #elapsed = time.time() - elapsed
        #print(f"All {ts.num_trees} trees coalesced")
        #print(f"Checked in {elapsed:.6g} secs")
        
        #labels = {i: string.ascii_lowercase[i] for i in range(ts.num_nodes)}
        labels = {i: str(i) for i in range(ts.num_nodes)}

        genome_order = [n for n in ts.first().nodes(order="minlex_postorder") if ts.node(n).is_sample()]
        labels.update({n: labels[i] for i, n in enumerate(genome_order)})
        style1 = (
        ".node:not(.sample) > .sym, .node:not(.sample) > .lab {visibility: hidden;}"
        ".mut {font-size: 12px} .y-axis .tick .lab {font-size: 85%}")
        #sz = (800, 250)  # size of the plot, slightly larger than the default
        ticks = [0, 5]
        #ts.draw_svg()
        ts.draw_text()
        #ts.draw_svg(
        #size=sz, node_labels=labels, style=style1, y_label="Time ago",
        #y_axis=True, y_ticks=ticks)
        
        self.assertEqual(True, False)  # add assertion here

def test_TS_kit_example(self):

        mutated_ts = tskit.load("data/whatis_example.trees")
        ts = mutated_ts.delete_sites(list(range(mutated_ts.num_sites)))
        # Extra code to label and order the tips alphabetically rather than numerically
        labels = {i: string.ascii_lowercase[i] for i in range(ts.num_nodes)}
        genome_order = [n for n in ts.first().nodes(order="minlex_postorder") if ts.node(n).is_sample()]
        labels.update({n: labels[i] for i, n in enumerate(genome_order)})
        style1 = (
            ".node:not(.sample) > .sym, .node:not(.sample) > .lab {visibility: hidden;}"
            ".mut {font-size: 12px} .y-axis .tick .lab {font-size: 85%}")
        sz = (800, 250)  # size of the plot, slightly larger than the default
        ticks = [0, 5000, 10000, 15000, 20000]
        ts.draw_svg(
            size=sz, node_labels=labels, style=style1, y_label="Time ago",
            y_axis=True, y_ticks=ticks)

if __name__ == '__main__':
    unittest.main()
