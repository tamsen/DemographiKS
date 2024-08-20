import os
import unittest

from matplotlib import pyplot as plt


class TestSteadyState(unittest.TestCase):
    def test_steady_state(self):

        out_folder="../SLiM_output"
        in_file_name = os.path.join(out_folder, "diploid_snm.txt")
        data=read_data(in_file_name)

        #cycle,Coalesced,numMutations,pi_site,theoretical_theta
        Ne=[int(d) for d in data["Ne"]]
        cycles_x=[int(d) for d in data["cycle"]]
        numMut_y=[int(d) for d in data["numMutations"]]
        piSite_y=[float(d) for d in data["pi_site"]]
        theta_y=[float(d) for d in data["theoretical_theta"]]
        coalesced=[d for d in data["Coalesced"]]
        initial_N=Ne[0]

        cycle_coalescence_achieved=False
        for i in range(0,len(coalesced)):
            if coalesced[i]=="T":
                cycle_coalescence_achieved=cycles_x[i]
                break


        fig, ax = plt.subplots(2, 1, figsize=(5, 5))
        fig.suptitle("Simple Diploid Neutral Model\n Ne=" + str(initial_N))
        coalescence_label="coalescence at cycle " + str(cycle_coalescence_achieved)
        ax[0].set_title("Total num mutations", fontsize=10)
        ax[0].set(ylabel="num mutations")
        ax[0].set(xlabel="cycle (num generations)")
        ax[0].plot(cycles_x,numMut_y, label="num mutations")

        if cycle_coalescence_achieved:
            ax[0].axvline(cycle_coalescence_achieved, label=coalescence_label, color="k")
        ax[0].legend()

        ax[1].set_title("Empirical Pi and Theoretical Theta", fontsize=10)
        ax[1].set(ylabel="diversity")
        ax[1].set(xlabel="cycle (num generations)")
        ax[1].plot(cycles_x,piSite_y, label="Pi (empirical)")
        ax[1].plot(cycles_x,theta_y, label="Theta (theoretical)")
        if cycle_coalescence_achieved:
            ax[1].axvline(cycle_coalescence_achieved, label=coalescence_label, color="k")

        ax[1].set(ylim=(0, max(theta_y)*1.25))
        ax[1].legend()

        plt.tight_layout()
        out_file_name = os.path.join(out_folder, "steady_state_N"+str(initial_N)+".png")
        plt.savefig(out_file_name)
        plt.close()

if __name__ == '__main__':
    unittest.main()

def read_data(csv_file):

    data_dict={}
    header_to_index={}
    index_to_header={}
    with open(csv_file, "r") as f:

        while True:

            line = f.readline()
            if not line:
                break

            line=line.replace("\n","")
            data_splat=line.split(',')
            if "cycle" in line:
                for i in range(0,len(data_splat)):
                    col_header_i=data_splat[i]
                    data_dict[col_header_i]=[]
                    header_to_index[col_header_i]=i
                    index_to_header[i]=col_header_i
                continue

            for i in range(0, len(data_splat)):
                col_header_i = index_to_header[i]
                data_dict[col_header_i].append(data_splat[i])

    return data_dict
